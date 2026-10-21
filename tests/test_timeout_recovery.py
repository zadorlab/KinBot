"""Confirmed walltime failures retry once without duplicating or mixing work."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

from ase import Atoms
from ase.db import connect
import pytest

from kinbot import qc_timeout, run_policy, run_recover, timeout_retry
from kinbot.anl import dispatch
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task
from tests.test_anl_dispatch import _mock_execution


def allocation(directory, job='42', status='TIMEOUT', walltime='1-00:00:00'):
    return {'JobId': job, 'JobState': status, 'TimeLimit': walltime,
            'Partition': 'long', 'WorkDir': str(directory), 'StartTime': '2026-10-20T01:00:00'}


@pytest.fixture
def pool(tmp_path, monkeypatch):
    run_policy.update(tmp_path, max_jobs=23, exclude=['bad'], timeout_retries=1,
                      minimum_walltime='7-00:00:00')
    monkeypatch.setattr(timeout_retry, 'partition_limits', lambda **kw: {'long': 31*86400})
    return tmp_path


@pytest.mark.parametrize('backend', ['molpro', 'mrcc', 'cfour', 'gaussian', 'geometry'])
def test_dispatch_all_backends_preserve_inputs_geometry_siblings_and_archive(pool, monkeypatch, backend):
    spec = dispatch_spec()
    bad = molpro_task('bad', 'rhf\n', geometry_from='initial')
    if backend == 'geometry':
        bad = deepcopy(spec['tasks'][0]); bad['id'] = 'bad'
    else:
        bad.update(backend=backend, command=[backend, '{input}'])
        if backend == 'mrcc': bad['input_name'] = 'MINP'
        if backend == 'cfour': bad['input_name'] = 'ZMAT'
        if backend == 'gaussian': bad.update(input_name='bad.com', stdin='bad.com')
    spec['tasks'] = [bad, molpro_task('live', 'rhf\n', geometry_from='initial'),
                     molpro_task('done', 'rhf\n', geometry_from='initial')]
    graph = dispatch.prepare(spec, pool/'graph')
    _mock_execution(graph, 'done')
    state = dispatch._load(graph)[2]
    state['tasks']['bad'].update(status='submitted', job_id='42', missing_since=0)
    state['tasks']['live'].update(status='submitted', job_id='43')
    atomic_json(graph/'state.json', state)
    directory = graph/'tasks/bad'
    (directory/'partial.chk').write_text('unfinished')
    original_geom = (directory/'geometry.xyz').read_bytes()
    original_input = {p.name: p.read_bytes() for p in directory.glob('*.inp')}
    sibling = {p: file_sha256(p) for name in ('live', 'done')
               for p in (graph/'tasks'/name).iterdir() if p.is_file()}
    monkeypatch.setattr(dispatch, '_job_active', lambda job: job == '43')
    monkeypatch.setattr(timeout_retry, 'slurm_job_record', lambda job: allocation(directory, job))
    result = dispatch.advance(graph)
    entry = result['tasks']['bad']
    assert entry['status'] == 'staged' and entry['attempt'] == 2
    assert result['tasks']['done']['status'] == 'complete'
    assert result['tasks']['live']['job_id'] == '43'
    assert sibling == {p: file_sha256(p) for p in sibling}
    assert (Path(entry['previous_attempt'])/'partial.chk').read_text() == 'unfinished'
    assert not (directory/'partial.chk').exists()
    assert (directory/'geometry.xyz').read_bytes() == original_geom
    assert original_input == {p.name: p.read_bytes() for p in directory.glob('*.inp')}
    assert '#SBATCH --time=7-00:00:00' in (directory/'job.slurm').read_text()
    assert dispatch.advance(graph)['timeout_retries'] == result['timeout_retries']
    # A second timeout is a hard failure; optional fallback and geometry
    # retries cannot silently extend the independent timeout budget.
    state = dispatch._load(graph)[2]
    state['tasks']['bad'].update(status='submitted', job_id='44', missing_since=0)
    atomic_json(graph/'state.json', state)
    result = dispatch.advance(graph)
    assert result['tasks']['bad']['failure_phase'] == 'walltime'
    assert 'limit reached' in result['tasks']['bad']['timeout_retry_blocked']


def small_graph(pool, monkeypatch):
    spec = dispatch_spec()
    spec['tasks'] = [molpro_task('bad', 'rhf\n', geometry_from='initial')]
    graph = dispatch.prepare(spec, pool/'graph')
    state = dispatch._load(graph)[2]
    state['tasks']['bad'].update(status='submitted', job_id='42', missing_since=0)
    atomic_json(graph/'state.json', state)
    monkeypatch.setattr(dispatch, '_job_active', lambda job: False)
    monkeypatch.setattr(timeout_retry, 'slurm_job_record', lambda job: allocation(graph/'tasks/bad', job))
    return graph


def test_transaction_resumes_after_archive_before_state_write(pool, monkeypatch):
    graph = small_graph(pool, monkeypatch)
    original = run_recover.atomic_json
    def interrupt(path, data):
        if Path(path) == graph/'state.json':
            raise OSError('interrupted')
        return original(path, data)
    monkeypatch.setattr(run_recover, 'atomic_json', interrupt)
    with pytest.raises(OSError, match='interrupted'):
        dispatch.advance(graph)
    monkeypatch.setattr(run_recover, 'atomic_json', original)
    state = dispatch.advance(graph)
    assert state['tasks']['bad']['attempt'] == 2
    assert len(state['timeout_retries']['bad']) == 1


@pytest.mark.parametrize('kind', ['active', 'unreachable', 'completed', 'cancelled', 'wrong_directory', 'unknown'])
def test_no_retry_without_inactive_failed_owned_timeout(pool, monkeypatch, kind):
    graph = small_graph(pool, monkeypatch)
    if kind == 'active':
        monkeypatch.setattr(dispatch, '_job_active', lambda job: True)
    elif kind == 'unreachable':
        def outage(job): raise dispatch.SchedulerUnavailableError('outage')
        monkeypatch.setattr(dispatch, '_job_active', outage)
        with pytest.raises(dispatch.SchedulerUnavailableError): dispatch.advance(graph)
    elif kind == 'completed':
        _mock_execution(graph, 'bad')
    else:
        record = allocation(graph/'tasks/bad')
        if kind == 'cancelled': record['JobState'] = 'CANCELLED'
        elif kind == 'wrong_directory': record['WorkDir'] = str(pool/'another')
        else: record = {}
        monkeypatch.setattr(timeout_retry, 'slurm_job_record', lambda job: record)
    if kind != 'unreachable': dispatch.advance(graph)
    assert not (graph/'timeout_recovery.json').exists()


def test_partition_cap_and_exhausted_budget_remain_hard_failure(pool, monkeypatch):
    graph = small_graph(pool, monkeypatch)
    monkeypatch.setattr(timeout_retry, 'partition_limits', lambda **kw: {'long': 86400})
    state = dispatch.advance(graph)
    assert state['tasks']['bad']['failure_phase'] == 'walltime'
    assert 'Partition maximum' in state['tasks']['bad']['timeout_retry_blocked']
    assert not (graph/'timeout_recovery.json').exists()


def test_scheduler_diagnostic_requires_exact_id_and_observed_workdir(pool, monkeypatch):
    record = allocation(pool, status='RUNNING')
    monkeypatch.setattr(timeout_retry, 'slurm_job_record', lambda job: record)
    timeout_retry.observe(pool, '42')
    path = pool/'slurm.stderr'
    monkeypatch.setattr(timeout_retry, 'slurm_job_record', lambda job: {})
    path.write_text('slurmstepd: error: *** JOB 420 ON n1 CANCELLED AT 2026-10-21T01:00:00 DUE TO TIME LIMIT ***\n')
    assert timeout_retry.evidence(pool, '42', [path]) is None
    path.write_text('slurmstepd: error: *** JOB 42 ON n1 CANCELLED AT 2026-10-27T01:00:00 DUE TO TIME LIMIT ***\n')
    assert timeout_retry.evidence(pool, '42', [path])['source'] == str(path)
    assert timeout_retry.next_walltime(run_policy.read(pool), timeout_retry.evidence(pool, '42', [path])) == '14-00:00:00'
    monkeypatch.setattr(timeout_retry, 'slurm_job_record', lambda job: allocation(pool, status='CANCELLED'))
    assert timeout_retry.evidence(pool, '42', [path]) is None


def test_growth_uses_actual_old_limit_and_never_shortens_or_changes_partition(pool, monkeypatch):
    settings = run_policy.read(pool)
    assert timeout_retry.next_walltime(settings, allocation_proof('7-00:00:00')) == '14-00:00:00'
    assert timeout_retry.next_walltime(settings, allocation_proof('20-00:00:00')) == '31-00:00:00'
    assert timeout_retry.next_walltime(settings, allocation_proof('31-00:00:00')) is None


def test_pending_allocation_does_not_freeze_an_unknown_start_time(pool, monkeypatch):
    record = allocation(pool, status='PENDING'); record['StartTime'] = 'Unknown'
    monkeypatch.setattr(timeout_retry, 'slurm_job_record', lambda job: record)
    timeout_retry.observe(pool, '42')
    assert not (pool/'.kinbot_timeout_allocations/42.json').exists()
    timeout_retry._OBSERVED.clear()
    record.update(JobState='RUNNING', StartTime='2026-10-20T01:00:00')
    timeout_retry.observe(pool, '42')
    assert json.loads((pool/'.kinbot_timeout_allocations/42.json').read_text())['StartTime'] == record['StartTime']


def allocation_proof(time):
    return {'previous_walltime': time, 'partition': 'long'}


@pytest.fixture
def ordinary(pool, monkeypatch):
    monkeypatch.chdir(pool)
    monkeypatch.setattr(qc_timeout, '_GRACE', 0)
    from kinbot import qc
    monkeypatch.setattr(qc, '_slurm_job_active', lambda job: False)
    monkeypatch.setattr(timeout_retry, 'slurm_job_record', lambda job: allocation(pool, job))
    Path('well.py').write_text('print("worker")\n')
    Path('well.sbatch').write_text('#!/bin/bash\n#SBATCH --time=1-00:00:00\n#SBATCH --exclusive\npython well.py\n')
    Path('well.chk').write_bytes(b'pre-allocation checkpoint')
    token = hashlib.sha256(b'well').hexdigest()[:20]
    baseline = pool/'.kinbot_slurm_io_retries'/token/'42'
    snapshots = {}
    for name in ('well.py', 'well.sbatch', 'well.chk'):
        saved = baseline/'initial'/name; saved.parent.mkdir(parents=True, exist_ok=True)
        saved.write_bytes(Path(name).read_bytes()); snapshots[name] = file_sha256(saved)
    atomic_json(baseline/'execution.json', {'identity': 'well', 'slurm_job_id': '42',
                'random_seed': '123456', 'initial_artifacts': snapshots, 'status': 'running'})
    Path('well.chk').write_bytes(b'interrupted checkpoint')
    Path('well.log').write_text('partial native result')
    Path('other.chk').write_text('sibling')
    db = connect(str(pool/'kinbot.db'))
    db.write(Atoms('H'), name='well', data={'status': 'normal', 'energy': -1.})
    worker = SimpleNamespace(job_ids={'well': '42'}, job_generation_tagged={'well': True}, db=db)
    calls = []
    def submit(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(returncode=0, stdout='Submitted batch job 45\n', stderr='')
    monkeypatch.setattr(run_policy, 'submit', submit)
    return worker, calls


def test_gaussian_retry_restores_only_verified_seeds_and_adopts_new_generation(ordinary):
    qc, calls = ordinary
    assert qc_timeout.recover(qc, 'well') == 'running'
    assert qc.job_ids['well'] == '45' and len(calls) == 1
    assert calls[0][1]['env']['SBATCH_TIMELIMIT'] == '7-00:00:00'
    assert Path('well.chk').read_bytes() == b'pre-allocation checkpoint'
    assert Path('other.chk').read_text() == 'sibling'
    assert not Path('well.log').exists()
    assert list(qc.db.select(name='well'))[-1].data['status'] == 0
    record = json.loads(qc_timeout._marker('well').read_text())
    assert (Path(record['attempts'][0]['archive'])/'files/well.chk').read_bytes() == b'interrupted checkpoint'
    subprocess.run(['bash', '-n', record['script']], check=True)
    qc.job_ids['well'] = '42'  # Simulate a controller with an older route ledger.
    assert qc_timeout.resume(qc, 'well') is None
    assert qc.job_ids['well'] == '45' and len(calls) == 1
    assert qc_timeout.recover(qc, 'well') is None
    assert 'limit reached' in qc._timeout_retry_blocked
    qc.job_ids['well'] = '99'  # A later constrained-path step is not this retry.
    Path('well.py').write_text('print("next step")')
    assert qc_timeout.resume(qc, 'well') is None
    assert qc.job_ids['well'] == '99'


def test_pending_admission_survives_controller_restart(ordinary, monkeypatch):
    qc, calls = ordinary
    original = run_policy.submit
    monkeypatch.setattr(run_policy, 'submit', lambda *a, **kw: None)
    assert qc_timeout.recover(qc, 'well') == 'running'
    assert qc.job_ids['well'] == '42'
    assert json.loads(qc_timeout._marker('well').read_text())['phase'] == 'ready'
    monkeypatch.setattr(run_policy, 'submit', original)
    qc_timeout.resume(qc, 'well')
    assert qc.job_ids['well'] == '45' and len(calls) == 1


def test_uncertain_submission_never_automatically_repeats(ordinary, monkeypatch):
    qc, calls = ordinary
    def lost_response(*args, **kwargs):
        raise subprocess.TimeoutExpired('sbatch', 30)
    monkeypatch.setattr(run_policy, 'submit', lost_response)
    with pytest.raises(subprocess.TimeoutExpired): qc_timeout.recover(qc, 'well')
    with pytest.raises(RuntimeError, match='uncertain'): qc_timeout.resume(qc, 'well')
    assert qc.job_ids['well'] == '42'
    record = json.loads(qc_timeout._marker('well').read_text())
    atomic_json(Path.cwd()/'.kinbot_policy_submissions/45.json', {
        'job_id': '45', 'directory': str(Path.cwd()), 'command': ['sbatch', record['script']]})
    qc_timeout.resume(qc, 'well')
    assert qc.job_ids['well'] == '45'


def test_corrupt_seed_does_not_submit(ordinary):
    qc, calls = ordinary
    seed = next(Path('.kinbot_slurm_io_retries').glob('*/42/initial/well.chk'))
    seed.write_bytes(b'changed')
    with pytest.raises(RuntimeError, match='seed is invalid'): qc_timeout.recover(qc, 'well')
    assert not calls and Path('well.chk').read_bytes() == b'interrupted checkpoint'


def test_ordinary_publication_grace_then_check_qc_adopts_retry(ordinary, monkeypatch):
    from kinbot.qc import QuantumChemistry
    qc, calls = ordinary
    qc.queuing = 'slurm'
    monkeypatch.setattr(QuantumChemistry, '_check_qc_with_grace', lambda *a: 'error')
    monkeypatch.setattr(qc_timeout, '_GRACE', 60.)
    assert QuantumChemistry.check_qc(qc, 'well') == 'running'
    assert not calls
    qc._timeout_publication_waits['42'] -= 61
    assert QuantumChemistry.check_qc(qc, 'well') == 'running'
    assert qc.job_ids['well'] == '45'
    assert len(calls) == 1
    # Finished results take precedence over even positive timeout evidence.
    monkeypatch.setattr(QuantumChemistry, '_check_qc_with_grace', lambda *a: 'normal')
    assert QuantumChemistry.check_qc(qc, 'well') == 'normal'
    assert len(calls) == 1


def test_ordinary_archive_transaction_is_resumable(ordinary, monkeypatch):
    qc, calls = ordinary
    original = qc_timeout.atomic_json
    def interrupt(path, value):
        if value.get('phase') == 'restoring':
            raise OSError('interrupt after archive')
        return original(path, value)
    monkeypatch.setattr(qc_timeout, 'atomic_json', interrupt)
    with pytest.raises(OSError, match='interrupt'): qc_timeout.recover(qc, 'well')
    assert not calls and not Path('well.chk').exists()
    monkeypatch.setattr(qc_timeout, 'atomic_json', original)
    qc_timeout.resume(qc, 'well')
    assert len(calls) == 1 and qc.job_ids['well'] == '45'
    assert Path('well.chk').read_bytes() == b'pre-allocation checkpoint'


def test_ordinary_post_submission_io_error_is_uncertain(ordinary, monkeypatch):
    qc, calls = ordinary
    def io_failure(*a, **kw): raise OSError('receipt write failed')
    monkeypatch.setattr(run_policy, 'submit', io_failure)
    with pytest.raises(OSError): qc_timeout.recover(qc, 'well')
    assert json.loads(qc_timeout._marker('well').read_text())['phase'] == 'submitting'
    with pytest.raises(RuntimeError, match='uncertain'): qc_timeout.resume(qc, 'well')


def test_tampered_dispatch_input_blocks_retry(pool, monkeypatch):
    graph = small_graph(pool, monkeypatch)
    (graph/'tasks/bad/bad.inp').write_text('changed chemistry')
    with pytest.raises(RuntimeError): dispatch.advance(graph)
    assert not (graph/'timeout_recovery.json').exists()


def test_dispatch_replacement_uses_pool_and_forces_increased_budget(pool, monkeypatch):
    graph = small_graph(pool, monkeypatch)
    commands = []
    monkeypatch.setenv('SBATCH_TIMELIMIT', '01:00:00')
    def submit(command, **kwargs):
        commands.append(command)
        assert kwargs['env']['SBATCH_TIMELIMIT'] == '7-00:00:00'
        assert kwargs['env']['SBATCH_PARTITION'] == 'long'
        assert command == ['sbatch', '--parsable', 'job.slurm']
        return SimpleNamespace(returncode=0, stdout='45\n', stderr='')
    monkeypatch.setattr(run_policy, 'submit', submit)
    state = dispatch.advance(graph, submit=True)
    assert state['tasks']['bad']['status'] == 'submitted'
    assert state['tasks']['bad']['job_id'] == '45'
    assert len(commands) == 1


def test_timeout_never_falls_back_to_optional_component_or_extra_l3_retry(pool, monkeypatch):
    graph = small_graph(pool, monkeypatch)
    spec = json.loads((graph/'workflow.json').read_text())
    spec['tasks'][0].update(failure_policy='skip_component', component_families=['test'])
    atomic_json(graph/'workflow.json', spec)
    state = json.loads((graph/'state.json').read_text())
    state['spec_sha256'] = file_sha256(graph/'workflow.json')
    # Restage the corresponding task.json so its policy matches the spec.
    task_file = graph/'tasks/bad/task.json'
    task = json.loads(task_file.read_text()); task['task'] = spec['tasks'][0]
    atomic_json(task_file, task)
    state['tasks']['bad']['task_sha256'] = file_sha256(task_file)
    atomic_json(graph/'state.json', state)
    monkeypatch.setattr(timeout_retry, 'partition_limits', lambda **kw: {'long': 86400})
    result = dispatch.advance(graph)
    assert result['tasks']['bad']['status'] == 'failed'
    from kinbot.anl import automatic
    monkeypatch.setattr(automatic, '_required_geometries', lambda spec: {'bad'})
    monkeypatch.setattr(dispatch, 'retry_failed', lambda *a: pytest.fail('extra L3 retry'))
    assert automatic._retry_l3_execution_failures(graph, {'intent': {'l3_max_retries': 10}}, result) == result


def test_old_unavailable_timeout_recovers_only_unambiguous_guard_generation(pool, monkeypatch):
    graph = small_graph(pool, monkeypatch)
    state = dispatch._load(graph)[2]
    state['tasks']['bad'].update(status='unavailable')
    state['tasks']['bad'].pop('job_id')
    atomic_json(graph/'state.json', state)
    path = graph/'tasks/bad/.kinbot_slurm_io_retries/token/42/execution.json'
    atomic_json(path, {'identity': 'bad', 'slurm_job_id': '42'})
    other = path.parent.parent/'43/execution.json'
    atomic_json(other, {'identity': 'bad', 'slurm_job_id': '43'})
    assert dispatch.advance(graph)['tasks']['bad']['status'] == 'unavailable'
    other.unlink()
    result = dispatch.advance(graph)
    assert result['tasks']['bad']['status'] == 'staged'
    assert result['timeout_retries']['bad'][0]['job_id'] == '42'


def test_dispatch_io_error_after_sbatch_keeps_uncertain_state(pool, monkeypatch):
    graph = small_graph(pool, monkeypatch)
    monkeypatch.setattr(run_policy, '_queue', lambda root: {})
    calls = []
    def sbatch(command, **kw):
        calls.append(command)
        return SimpleNamespace(returncode=0, stdout='45\n', stderr='')
    monkeypatch.setattr(run_policy.subprocess, 'run', sbatch)
    original = run_policy.atomic_json
    def fail_receipt(path, data):
        if Path(path).parent.name == '.kinbot_policy_submissions':
            raise OSError('filesystem failed after sbatch accepted')
        return original(path, data)
    monkeypatch.setattr(run_policy, 'atomic_json', fail_receipt)
    with pytest.raises(RuntimeError, match='submission outcome uncertain'):
        dispatch.advance(graph, submit=True)
    assert dispatch._load(graph)[2]['tasks']['bad']['status'] == 'submitting'
    with pytest.raises(RuntimeError, match='submission outcome uncertain'):
        dispatch.advance(graph, submit=True)
    assert len(calls) == 1


@pytest.mark.parametrize('key,value', [('timeout_retries', True), ('timeout_retries', -1),
    ('timeout_retries', 11), ('timeout_walltime_factor', 1), ('timeout_walltime_factor', float('nan'))])
def test_input_and_live_policy_reject_invalid_retry_controls(pool, key, value):
    with pytest.raises(ValueError): run_policy.update(pool, **{key: value})
    with pytest.raises(ValueError): run_policy.validate_options({'max_jobs': 23, key: value})
