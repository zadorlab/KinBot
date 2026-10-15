"""Recover Slurm's pre-shell I/O failure without hiding chemistry failures."""
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import molpro_reconnaissance as recon
from kinbot.execution import atomic_json, file_sha256, slurm_job_record


def _manifest(root, identifier='71'):
    directory = root / 'probe'
    directory.mkdir(parents=True)
    inp = directory / 'probe.inp'
    inp.write_text('immutable chemistry\n')
    (directory / 'probe.sbatch').write_text('#!/bin/bash\n#SBATCH -n 8\ntrue\n')
    job = {'name': 'probe', 'directory': str(directory), 'status': 'incomplete',
           'input_sha256': file_sha256(inp), 'submission_status': 'submitted',
           'missing_since': 0., 'job_id': identifier}
    manifest = {'jobs': [job], 'status': 'incomplete'}
    path = atomic_json(root / 'reconnaissance.json', manifest)
    record = {'JobId': identifier, 'WorkDir': str(directory), 'JobState': 'FAILED',
              'ExitCode': '0:53', 'Reason': 'RaisedSignal:53(Real-time_signal_19)',
              'BatchHost': 'faulty-host', 'RunTime': '00:00:01'}
    return path, manifest, job, record


@pytest.mark.parametrize('reason', ['unknown', 'cancelled', 'signal9', 'wrong_job',
    'wrong_directory', 'no_node', 'native', 'output', 'execution', 'stderr', 'guard',
    'guard_started', 'live'])
def test_no_automatic_retry_without_proof_of_prestart_failure(tmp_path, monkeypatch, reason):
    path, manifest, job, record = _manifest(tmp_path)
    if reason == 'unknown':
        record = {}
    elif reason == 'cancelled':
        record['JobState'] = 'CANCELLED'
    elif reason == 'signal9':
        record['ExitCode'] = '0:9'
    elif reason == 'wrong_job':
        record['JobId'] = '999'
    elif reason == 'wrong_directory':
        record['WorkDir'] = '/another/campaign'
    elif reason == 'no_node':
        record['BatchHost'] = '(null)'
    elif reason == 'guard_started':
        (Path(job['directory']) / '.kinbot_slurm_io_retries/token/71').mkdir(parents=True)
    elif reason in ('native', 'output', 'execution', 'stderr', 'guard'):
        filename = {'native': 'probe.native.inp', 'output': 'probe.out',
                    'execution': 'execution.json', 'stderr': 'probe.stderr',
                    'guard': '.kinbot_slurm_io_retries/token/71/execution.json'}[reason]
        artifact = Path(job['directory']) / filename
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text('native or shell has started')
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: record)
    live = {'71'} if reason == 'live' else set()
    recon._recover_slurm_startup_failures(path, manifest, live)
    assert job['job_id'] == '71'
    assert not job.get('startup_retry_history')
    assert not manifest.get('slurm_excluded_nodes')


def test_startup_retry_journals_evidence_and_survives_controller_restart(tmp_path, monkeypatch):
    path, manifest, job, record = _manifest(tmp_path)
    job.pop('missing_since')
    now = [100.]
    monkeypatch.setattr(recon.time, 'time', lambda: now[0])
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: record.copy())
    recon._recover_slurm_startup_failures(path, manifest, set())
    assert job['job_id'] == '71'  # Publication grace; save scheduler evidence now.
    assert job['scheduler_record'] == record
    # A new controller can recover even after scontrol purges the old record.
    manifest = json.loads(path.read_text())
    job = manifest['jobs'][0]
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: {})
    now[0] = 161.
    recon._recover_slurm_startup_failures(path, manifest, set())
    assert 'job_id' not in job
    assert job['submission_status'] == 'retry_ready'
    assert manifest['slurm_excluded_nodes'] == ['faulty-host']
    assert job['startup_retry_history'][0]['scheduler_record'] == record
    assert json.loads(path.read_text()) == manifest
    assert (Path(job['directory']) / 'probe.inp').read_text() == 'immutable chemistry\n'


def test_startup_retries_are_bounded_and_do_not_lose_final_failure(tmp_path, monkeypatch):
    path, manifest, job, record = _manifest(tmp_path)
    monkeypatch.setattr(recon, 'slurm_job_record', lambda identifier: dict(record, JobId=identifier))
    for attempt in range(3):
        job.update(job_id=str(71 + attempt), missing_since=0.)
        recon._recover_slurm_startup_failures(path, manifest, set())
    assert len(job['startup_retry_history']) == 2
    assert job['job_id'] == '73'
    with pytest.raises(recon.ProbeFailure, match=r'retry limit \(2\) exhausted'):
        recon._finished(path, manifest, set())


def test_requeued_job_is_reserved_even_before_squeue_sees_it(tmp_path, monkeypatch):
    path, manifest, job, record = _manifest(tmp_path)
    record.update(JobState='PENDING', Reason='BeginTime', ExitCode='0:0')
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: record)
    live = set()
    recon._recover_slurm_startup_failures(path, manifest, live)
    assert live == {'71'} and 'missing_since' not in job
    assert not recon._finished(path, manifest, live)
    assert job['job_id'] == '71'


def test_startup_retry_refuses_changed_input(tmp_path, monkeypatch):
    path, manifest, job, record = _manifest(tmp_path)
    (Path(job['directory']) / 'probe.inp').write_text('changed chemistry')
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: record)
    with pytest.raises(ValueError, match='input changed'):
        recon._recover_slurm_startup_failures(path, manifest, set())
    assert job['job_id'] == '71'


def test_cached_active_record_cannot_keep_a_missing_job_live_forever(tmp_path, monkeypatch):
    path, manifest, job, record = _manifest(tmp_path)
    job['scheduler_record'] = dict(record, JobState='RUNNING', ExitCode='0:0')
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: {})
    live = set()
    recon._recover_slurm_startup_failures(path, manifest, live)
    assert not live and job['job_id'] == '71'
    with pytest.raises(recon.ProbeFailure):
        recon._finished(path, manifest, live)


def test_exclusions_preserve_site_hostlists_without_changing_parent_environment(monkeypatch):
    monkeypatch.setenv('SBATCH_EXCLUDE', 'site[01-03]')
    env = recon._slurm_submission_environment(
        '#!/bin/bash\n#SBATCH -x original\n#SBATCH --exclude=other[04-06]\ntrue\n',
        ['faulty-host'])
    assert env['SBATCH_EXCLUDE'] == 'site[01-03],faulty-host,original,other[04-06]'
    assert recon.os.environ['SBATCH_EXCLUDE'] == 'site[01-03]'


@pytest.mark.parametrize('limit', [2, 20])
def test_shared_controller_retries_prestart_failure_adopts_live_jobs_and_honors_cap(
        tmp_path, monkeypatch, limit):
    paths, records = [], {}
    for index in range(limit + 2):
        path, manifest, job, record = _manifest(tmp_path / str(index), str(71 + index))
        paths.append(path)
        if index == 0:
            records[job['job_id']] = record
        elif index > 1:
            job.pop('job_id')
            atomic_json(path, manifest)
    live, submitted, tick = {'72'}, [], [0]
    protected = {p: p.read_bytes() for p in tmp_path.rglob('*.inp')}
    monkeypatch.setattr(recon, 'slurm_job_record', lambda identifier: records.get(identifier, {}))
    monkeypatch.setattr(recon, '_audit', lambda path, **kw: json.loads(Path(path).read_text()))

    def scheduler(command, **kwargs):
        if command[0] == 'squeue':
            return SimpleNamespace(stdout='\n'.join(live), returncode=0)
        assert command[0] == 'sbatch'
        assert 'faulty-host' in kwargs['env']['SBATCH_EXCLUDE']
        identifier = str(100 + len(submitted))
        submitted.append(kwargs['cwd'])
        live.add(identifier)
        assert len(live) <= limit
        return SimpleNamespace(stdout=identifier, stderr='', returncode=0)

    def complete_tick(_):
        tick[0] += 1
        assert tick[0] < 5
        for path in paths:
            manifest = json.loads(path.read_text())
            job = manifest['jobs'][0]
            if job.get('job_id') in live:
                job['status'] = 'complete'
                manifest['status'] = 'complete_review_required'
                atomic_json(path, manifest)
        live.clear()

    monkeypatch.setattr(recon.subprocess, 'run', scheduler)
    monkeypatch.setattr(recon.time, 'sleep', complete_tick)
    reports = recon.run_many(paths, max_jobs=limit, startup_checked=True)
    assert reports[1]['jobs'][0]['job_id'] == '72'  # Healthy job was adopted.
    assert len(submitted) == len(paths) - 1 and len(set(submitted)) == len(submitted)
    assert len(reports[0]['jobs'][0]['startup_retry_history']) == 1
    assert all(m['slurm_excluded_nodes'] == ['faulty-host'] for m in reports)
    assert all(p.read_bytes() == original for p, original in protected.items())


def test_scontrol_evidence_needs_no_accounting(monkeypatch):
    calls = []
    def scheduler(command, **kwargs):
        calls.append(command)
        assert kwargs['timeout'] == 15
        return SimpleNamespace(returncode=0, stdout=(
            'JobId=71 JobState=FAILED Reason=RaisedSignal:53(Real-time_signal_19) '
            'ExitCode=0:53 WorkDir=/campaign/job BatchHost=faulty-host'))
    monkeypatch.setattr(recon.subprocess, 'run', scheduler)
    assert slurm_job_record('71')['ExitCode'] == '0:53'
    assert calls == [['scontrol', 'show', 'job', '-o', '71']]


def test_documented_incident_repair_preserves_live_and_started_jobs(tmp_path, monkeypatch):
    import sys
    from kinbot import execution
    records, systems, paths = {}, [], []
    for index, identifier in enumerate(range(52907719, 52907726)):
        path, manifest, job, record = _manifest(tmp_path / str(index), str(identifier))
        paths.append(path)
        systems.append({'pilot_manifest': str(path), 'manifest': str(path)})
        if index == 0:
            record = {}  # Retained scheduler record has expired: explicit incident retry.
        elif index == 3:
            (Path(job['directory']) / 'probe.native.inp').write_text('worker started')
        elif index == 4:
            record['JobState'] = 'CANCELLED'
        elif index == 5:
            record['WorkDir'] = '/another/campaign'
        elif index == 6:
            job['status'] = 'complete'
            atomic_json(path, manifest)
        records[str(identifier)] = record
    atomic_json(tmp_path / 'campaign.json', {'systems': systems})
    monkeypatch.setattr(sys, 'argv', ['-', str(tmp_path)])
    monkeypatch.setattr(execution, 'slurm_job_record', lambda jid: records[jid])
    monkeypatch.setattr(recon.subprocess, 'run', lambda *a, **kw:
                        SimpleNamespace(stdout='52907721\n', returncode=0))
    retried = []
    def retry(path, *, include_interrupted):
        assert include_interrupted
        manifest = json.loads(path.read_text())
        retried.append(manifest['jobs'][0].pop('job_id'))
        atomic_json(path, manifest)
        return {}
    monkeypatch.setattr(recon, 'retry_failed', retry)
    document = (Path(__file__).parents[1] /
                'docs/SULFUR_REMAINING_CHECKS_2026-10-15.md').read_text()
    section = document.split('## Resume Batch A after the node21 startup failures')[1]
    snippet = section.split('<<\'PY\'\n', 1)[1].split('\nPY\n', 1)[0]
    exec(compile(snippet, '<documented incident repair>', 'exec'), {})
    exec(compile(snippet, '<repeated incident repair>', 'exec'), {})
    assert retried == ['52907719', '52907720']
    assert all(json.loads(p.read_text())['slurm_excluded_nodes'] == ['node21'] for p in paths)
    assert json.loads(paths[2].read_text())['jobs'][0]['job_id'] == '52907721'
