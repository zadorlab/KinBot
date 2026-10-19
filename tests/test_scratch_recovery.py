"""Scratch faults retry on another allocation without losing native work."""
import errno
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import qc_scratch, run_recover
from kinbot.anl import dispatch, scratch_recovery
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task
from tests.test_anl_dispatch import _mock_execution


ERROR = ('[Errno 28] No verified local QC scratch is usable: '
         "account_scratch=/scratch/tester: [Errno 13] Permission denied: '/scratch/tester'; "
         'TMPDIR=/tmp: requires 63000 MB; 37946 MB free')


def failed(graph, *, job='42', node='badnode', legacy=True, error=ERROR):
    state = dispatch._load(graph)[2]
    entry = state['tasks']['bad']
    entry.update(status='failed', job_id=job)
    atomic_json(graph/'state.json', state)
    directory = graph/'tasks/bad'
    record = {'schema': 1, 'task_id': 'bad', 'geometry_sha256': entry['geometry_sha256'],
              'status': 'failed', 'failure_phase': 'execution', 'error': error,
              'traceback': 'File "/kinbot/qc_scratch.py", in allocate\nraise OSError'}
    if legacy:
        atomic_json(directory/f'.kinbot_slurm_io_retries/hash/{job}/execution.json', {
            'identity': 'bad', 'slurm_job_id': job, 'status': 'failed', 'host': node})
    else:
        record.update(failure_phase='scratch_setup', retryable_scratch_node=True,
                      slurm_job_id=job, execution_node=node)
    atomic_json(directory/'execution.json', record)


@pytest.fixture
def graph(tmp_path, monkeypatch):
    monkeypatch.setenv('SBATCH_EXCLUDE', 'sitebad')
    monkeypatch.setenv('KINBOT_QC_SCRATCH_RETRIES', '2')
    monkeypatch.setattr(scratch_recovery, 'slurm_job_record', lambda job: {})
    spec = dispatch_spec()
    spec['tasks'] = [molpro_task(n, 'rhf\n', geometry_from='initial')
                     for n in ('bad', 'live', 'done', 'ready', 'child')]
    spec['tasks'][-1]['depends_on'] = ['bad']
    root = dispatch.prepare(spec, tmp_path/'graph')
    failed(root)
    _mock_execution(root, 'done')
    state = dispatch._load(root)[2]
    state['tasks']['live'].update(status='submitted', job_id='43')
    state['tasks']['child'] = {'status': 'blocked', 'blocked_by': ['bad']}
    atomic_json(root/'state.json', state)
    monkeypatch.setattr(dispatch, '_job_active', lambda job: str(job) == '43')
    return root


@pytest.mark.parametrize('legacy', [True, False])
def test_retry_archives_only_failed_attempt_and_preserves_science_and_jobs(graph, legacy):
    failed(graph, legacy=legacy)
    specification = (graph/'workflow.json').read_bytes()
    unchanged = {p: file_sha256(p) for name in ('live', 'done')
                 for p in (graph/'tasks'/name).iterdir() if p.is_file()}
    state = dispatch.advance(graph)
    entry = state['tasks']['bad']
    assert entry['status'] == 'staged' and entry['attempt'] == 2
    assert state['tasks']['live']['job_id'] == '43'
    assert state['tasks']['live']['status'] == 'submitted'
    assert state['tasks']['done']['status'] == 'complete'
    assert state['tasks'].get('child', {}).get('status') != 'blocked'
    assert json.loads((Path(entry['previous_attempt'])/'execution.json').read_text())['error'] == ERROR
    assert (graph/'workflow.json').read_bytes() == specification
    assert unchanged == {p: file_sha256(p) for p in unchanged}
    script = (graph/'tasks/bad/job.slurm').read_text()
    assert '#SBATCH --exclude=badnode,sitebad\n' in script
    assert '--exclusive' in script
    assert state['scratch_node_retries']['bad'][0]['node'] == 'badnode'
    assert dispatch.advance(graph) == state


def test_cap_is_two_replacements_and_exhaustion_is_hard_failure(graph, monkeypatch):
    dispatch.advance(graph)
    failed(graph, job='44', node='otherbad')
    state = dispatch.advance(graph)
    assert state['tasks']['bad']['attempt'] == 3
    assert state['scratch_excluded_nodes'] == ['badnode', 'otherbad']
    failed(graph, job='45', node='thirdbad')
    state = dispatch.advance(graph)
    assert state['tasks']['bad']['status'] == 'failed'
    assert state['tasks']['bad']['failure_phase'] == 'scratch_setup'
    assert 'limit reached (2)' in state['tasks']['bad']['scratch_retry_blocked']
    assert len(state['scratch_node_retries']['bad']) == 2
    from kinbot.anl import automatic
    # The independent L3 retry mechanism must not extend this budget.
    monkeypatch.setattr(automatic, '_required_geometries', lambda spec: {'bad'})
    monkeypatch.setattr(dispatch, 'retry_failed', lambda *a: pytest.fail('extra L3 retry'))
    assert automatic._retry_l3_execution_failures(graph, {'intent': {'l3_max_retries': 10}}, state) == state


def test_disabled_retry_stays_hard_failure(graph, monkeypatch):
    monkeypatch.setenv('KINBOT_QC_SCRATCH_RETRIES', '0')
    state = dispatch.advance(graph)
    assert state['tasks']['bad']['failure_phase'] == 'scratch_setup'
    assert state['tasks']['bad']['status'] == 'failed'
    assert not (graph/'scratch_recovery.json').exists()


@pytest.mark.parametrize('change', ['other_permission', 'missing_trace', 'wrong_provenance', 'missing_node'])
def test_unverified_or_unrelated_errors_do_not_retry(graph, change):
    path = graph/'tasks/bad/execution.json'
    record = json.loads(path.read_text())
    if change == 'other_permission':
        record['error'] = '[Errno 13] Permission denied: /opt/qc/program'
    elif change == 'missing_trace':
        record['traceback'] = ''
    elif change == 'wrong_provenance':
        record['geometry_sha256'] = 'wrong'
    else:
        (graph/'tasks/bad/.kinbot_slurm_io_retries/hash/42/execution.json').unlink()
    atomic_json(path, record)
    if change == 'wrong_provenance':
        with pytest.raises(RuntimeError, match='execution record does not match'):
            dispatch.advance(graph)
    else:
        state = dispatch.advance(graph)
        assert state['tasks']['bad']['status'] == 'failed'
        if change == 'missing_node':
            assert 'cannot be verified' in state['tasks']['bad']['scratch_retry_blocked']
    assert not (graph/'scratch_recovery.json').exists()


def test_active_intermediate_failure_waits_for_queue_exit(graph, monkeypatch):
    state = dispatch._load(graph)[2]
    state['tasks']['bad']['status'] = 'submitted'
    atomic_json(graph/'state.json', state)
    monkeypatch.setattr(dispatch, '_job_active', lambda job: True)
    assert dispatch.advance(graph)['tasks']['bad']['status'] == 'submitted'
    assert not (graph/'scratch_recovery.json').exists()
    monkeypatch.setattr(dispatch, '_job_active', lambda job: str(job) == '43')
    assert dispatch.advance(graph)['tasks']['bad']['status'] == 'staged'


def test_archive_interruption_resumes_without_second_retry(graph, monkeypatch):
    original = run_recover.atomic_json
    def interrupted(path, data):
        if Path(path) == graph/'state.json':
            raise OSError('interrupted after moving attempt')
        return original(path, data)
    monkeypatch.setattr(run_recover, 'atomic_json', interrupted)
    with pytest.raises(OSError, match='interrupted'):
        dispatch.advance(graph)
    assert not (graph/'tasks/bad').exists()
    monkeypatch.setattr(run_recover, 'atomic_json', original)
    state = dispatch.advance(graph)
    assert state['tasks']['bad']['attempt'] == 2
    assert state['tasks']['live']['job_id'] == '43'
    assert len(list((graph/'attempts').glob('scratch_node_*'))) == 1


def test_previously_staged_script_gets_new_exclusions_at_submission(graph, monkeypatch):
    science = {p: file_sha256(p) for p in (graph/'tasks/ready').iterdir()
               if p.name != 'job.slurm'}
    calls = []
    def submit(command, **kwargs):
        calls.append(command)
        assert command == ['sbatch', '--parsable', 'job.slurm']
        script = (Path(kwargs['cwd'])/'job.slurm').read_text()
        assert '#SBATCH --exclude=badnode,sitebad\n' in script
        assert '--exclusive' in script
        return SimpleNamespace(returncode=0, stdout='44\n', stderr='')
    monkeypatch.setattr(dispatch.subprocess, 'run', submit)
    state = dispatch.advance(graph, submit=True, submit_only={'ready'})
    assert len(calls) == 1
    assert state['tasks']['ready']['job_id'] == '44'
    assert state['tasks']['live']['job_id'] == '43'
    assert science == {p: file_sha256(p) for p in science}
    dispatch._verify_stage_files(graph, next(t for t in dispatch._load(graph)[1]['tasks']
                                            if t['id'] == 'ready'), state['tasks']['ready'])


def test_exclusion_merging_is_idempotent_and_preserves_hostlists(monkeypatch):
    monkeypatch.setenv('SBATCH_EXCLUDE', 'node4')
    original = '#!/bin/bash\n#SBATCH --exclusive\n#SBATCH --exclude=node[1,2],node3\ncommand\n'
    state = {'scratch_excluded_nodes': ['node3', 'node5']}
    updated = scratch_recovery.apply_exclusions(original, state)
    assert 'node[1,2]' in updated and 'node4' in updated and 'node5' in updated
    assert updated.count('--exclude') == 1
    assert scratch_recovery.apply_exclusions(updated, state) == updated


@pytest.mark.parametrize('kind', ['external', 'ase_optimize'])
def test_worker_records_node_and_structured_scratch_failure(tmp_path, monkeypatch, kind):
    spec = dispatch_spec()
    spec['tasks'] = ([molpro_task('bad', 'rhf\n', geometry_from='initial')]
                     if kind == 'external' else [spec['tasks'][0]])
    graph = dispatch.prepare(spec, tmp_path/'graph')
    ident = spec['tasks'][0]['id']
    monkeypatch.setenv('SLURM_JOB_ID', '42')
    monkeypatch.setenv('SLURMD_NODENAME', 'badnode')
    def failure(*args):
        try:
            raise qc_scratch.ScratchUnavailableError(errno.EACCES, ERROR, retryable_node=True)
        except OSError as error:
            raise RuntimeError('wrapped scratch failure') from error
    monkeypatch.setattr(dispatch, '_run_external' if kind == 'external' else '_run_ase_optimize', failure)
    with pytest.raises(RuntimeError, match='wrapped scratch failure'):
        dispatch.run_task(graph/'tasks'/ident/'task.json')
    result = json.loads((graph/'tasks'/ident/'execution.json').read_text())
    assert result['slurm_job_id'] == '42' and result['execution_node'] == 'badnode'
    assert result['failure_phase'] == 'scratch_setup' and result['retryable_scratch_node']
    assert scratch_recovery.eligible(result)
