"""Purged jobs need positive submission/query evidence, not an empty lookup."""
import json
import time
from datetime import datetime, timezone, timedelta
from types import SimpleNamespace

import pytest

from kinbot import live_control as control, run_policy
from kinbot.anl import dispatch, live_control as owner, lost_job
from kinbot.execution import atomic_json, file_sha256
from tests.test_live_control import graph, request_path
from tests.test_anl_dispatch import _mock_execution


def setup(tmp_path, monkeypatch):
    control.configure(tmp_path)
    g = graph(tmp_path)
    _mock_execution(g, 'done')
    dispatch.advance(g)
    state = dispatch._load(g)[2]
    state['tasks']['change'].update(status='failed', job_id='42', failure_phase='scheduler',
        error='Slurm job disappeared without execution.json.', missing_since=time.time()-3600)
    state['tasks']['live'].update(status='submitted', job_id='43')
    atomic_json(g/'state.json', state)
    d = g/'tasks/change'
    admission = tmp_path/'.kinbot_policy_submissions/42.json'
    atomic_json(admission, dict(job_id='42', status='submitted', directory=str(d),
        command=['sbatch','--parsable','job.slurm'],
        created_at=(datetime.now(timezone.utc)-timedelta(hours=2)).isoformat()))
    # Use the normal test policy so all dispatcher policy reads remain valid.
    run_policy.update(tmp_path, max_jobs=23, exclude=['bad'])
    monkeypatch.setattr(owner, 'slurm_job_record', lambda job: {})
    monkeypatch.setattr(dispatch, '_job_active', lambda job: job == '43')
    def query(cmd, **kwargs):
        if cmd[0] == 'scontrol':
            return SimpleNamespace(returncode=1, stdout='', stderr='slurm_load_jobs error: Invalid job id specified\n')
        return SimpleNamespace(returncode=0, stdout='43\n', stderr='')
    monkeypatch.setattr(lost_job.subprocess, 'run', query)
    return g, d, admission


def test_explicit_lost_job_retry_archives_and_preserves_siblings(tmp_path, monkeypatch):
    g, d, _ = setup(tmp_path, monkeypatch)
    siblings = {p: file_sha256(p) for n in ('live','done')
                for p in (g/'tasks'/n).iterdir() if p.is_file()}
    r = control.request(tmp_path, d, reason='Reviewed lost allocation; retry unchanged task')
    owner.process(g)
    request = json.loads(request_path(tmp_path,r).read_text())
    assert request['status'] == 'applied'
    assert request['lost_job_evidence']['original_exit'] == 'unknown'
    state = dispatch.advance(g)
    assert state['tasks']['change']['status'] == 'staged'
    assert state['tasks']['change']['attempt'] == 2
    assert siblings == {p: file_sha256(p) for p in siblings}
    owner.process(g)
    assert dispatch._load(g)[2]['tasks']['change']['attempt'] == 2


@pytest.mark.parametrize('bad', ['no_receipt','wrong_directory','wrong_job','future_receipt',
    'output','empty_stderr','guard','changed_input','live','queue_down','controller_down','empty_query','grace'])
def test_lost_retry_refuses_uncertainty_or_worker_evidence(tmp_path, monkeypatch, bad):
    g, d, admission = setup(tmp_path, monkeypatch)
    if bad == 'no_receipt':
        admission.unlink()
    if bad in ('wrong_directory','wrong_job','future_receipt'):
        a=json.loads(admission.read_text())
        a.update({'wrong_directory': {'directory':str(g)}, 'wrong_job': {'job_id':'99'},
                  'future_receipt': {'created_at':datetime.now(timezone.utc).isoformat()}}[bad])
        atomic_json(admission,a)
    if bad == 'output': (d/'change.out').write_text('native result')
    if bad == 'empty_stderr': (d/'slurm.stderr').touch()
    if bad == 'guard': (d/'.kinbot_slurm_io_retries').mkdir()
    if bad == 'changed_input': (d/'change.inp').write_text('edited')
    if bad == 'grace':
        s=dispatch._load(g)[2]; s['tasks']['change']['missing_since']=time.time()
        atomic_json(g/'state.json',s)
    original=lost_job.subprocess.run
    def query(cmd,**kwargs):
        if cmd[0]=='scontrol' and bad in ('controller_down','empty_query'):
            return SimpleNamespace(returncode=1,stdout='',stderr='Controller unavailable' if bad=='controller_down' else '')
        if cmd[0]=='squeue' and bad in ('live','queue_down'):
            return SimpleNamespace(returncode=int(bad=='queue_down'),stdout='42\n',stderr='')
        return original(cmd,**kwargs)
    monkeypatch.setattr(lost_job.subprocess,'run',query)
    r=control.request(tmp_path,d,reason='Explicit retry')
    owner.process(g)
    assert json.loads(request_path(tmp_path,r).read_text())['status'] != 'applied'
    assert (d/'task.json').is_file()
    assert not (g/'.live_control_transaction.json').exists()
