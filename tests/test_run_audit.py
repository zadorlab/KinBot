"""An audit must observe live work without advancing or repairing it."""
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import run_audit
from kinbot.anl import dispatch
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task


def snapshot(root):
    return {str(p.relative_to(root)): file_sha256(p) for p in root.rglob('*') if p.is_file()}


@pytest.fixture
def run(tmp_path, monkeypatch):
    root = tmp_path/'run'; case = root/'example'; worker = case/'123'
    worker.mkdir(parents=True)
    params = dict(queue_exclusive=1, rotdpy_exclusive=1, queue_job_limit=10,
                  slurm_exclude_nodes=['node1', 'node9', 'node10'],
                  qc_scratch_fallbacks=[], qc_scratch_min_mb=4096)
    atomic_json(case/'input.json', params)
    atomic_json(worker/'123.json', params)
    atomic_json(root/'prepared.json', {'inputs': [str(case/'input.json')]})
    spec = dispatch_spec()
    spec['tasks'] = [molpro_task('seed', 'rhf\n', geometry_from='initial'),
                     molpro_task('opt', 'rhf\n', geometry_from='seed', geometry_output='final.xyz'),
                     molpro_task('energy', 'rhf\n', geometry_from='opt')]
    spec['tasks'][0]['geometry_output'] = 'final.xyz'
    spec['tasks'][0]['required_outputs'].append('final.xyz')
    graph = worker/'anl_automatic/states/species/composite'
    monkeypatch.setenv('SBATCH_EXCLUDE', 'node1,node9,node10')
    dispatch.prepare(spec, graph)
    state = dispatch._load(graph)[2]
    state['tasks']['seed']['status'] = 'complete'
    atomic_json(graph/'state.json', state)
    atomic_json(root/'outage_recovery.json', dict(phase='launched', time='then',
        graphs=[dict(graph=str(graph), complete=['seed'], active={}, reset=['opt', 'energy'])]))
    monkeypatch.setattr('kinbot.run_cleanup._processes', lambda root: [])
    for name in ('advance', 'prepare', 'preflight_until_scheduler_ready'):
        monkeypatch.setattr(dispatch, name, lambda *a, **kw: pytest.fail('audit must not mutate'))
    return SimpleNamespace(root=root, graph=graph, worker=worker, state=state, params=params)


def test_full_audit_is_read_only_and_explains_dependencies(run):
    before = snapshot(run.root)
    report = run_audit.audit(run.root, scheduler=([], None), details={})
    assert before == snapshot(run.root)
    rows = {r['task']: r for r in report['cases'][0]['tasks']}
    assert rows['seed']['salvaged_at_recovery'] and rows['seed']['status'] == 'complete'
    assert rows['opt']['status'] == 'ready_to_stage'
    assert rows['energy']['status'] == 'dependency_wait'
    assert rows['energy']['waiting_on'] == {'opt': 'unstarted'}
    assert rows['energy']['dependency_level'] == 2
    assert rows['opt']['scratch_policy']['status'] == 'unobserved'
    assert 'not a global queue' in run_audit.render(report, full=True)
    assert not report['cases'][0]['errors']


def test_scheduler_ownership_and_policies_are_observed_separately(run):
    directory = run.graph/'tasks/seed'
    job = dict(job_id='77', directory=str(directory), state='RUNNING', nodes='node2')
    detail = dict(JobId='77', OverSubscribe='NO', ExcNodeList='node[1,9-10]', NodeList='node2')
    before = snapshot(run.root)
    report = run_audit.audit(run.root, scheduler=([job], None), details={'77': detail})
    row = report['cases'][0]['tasks'][0]
    assert row['status'] == 'allocated'
    assert row['script_policy']['status'] == row['scheduler_policy']['status'] == 'pass'
    assert row['scratch_policy']['status'] == 'unobserved'
    assert before == snapshot(run.root)
    detail['ExcNodeList'] = 'node1'
    assert run_audit.scheduler_policy(job, detail, run.params['slurm_exclude_nodes'])['status'] == 'fail'
    detail.update(ExcNodeList='node[1,9-10]', OverSubscribe='YES')
    assert run_audit.scheduler_policy(job, detail, run.params['slurm_exclude_nodes'])['status'] == 'fail'


def test_missing_scheduler_does_not_mean_finished(run):
    run.state['tasks']['seed'].update(status='submitted', job_id='77')
    atomic_json(run.graph/'state.json', run.state)
    report = run_audit.audit(run.root, scheduler=([], 'timed out'), details={})
    assert report['cases'][0]['tasks'][0]['status'] == 'scheduler_unknown'
    report = run_audit.audit(run.root, scheduler=([], None), details={})
    assert report['cases'][0]['tasks'][0]['status'] == 'awaiting_reconciliation'


@pytest.mark.parametrize('value,expected', [
    ('node[1,9-10]', {'node1', 'node9', 'node10'}),
    ('rack[01-02]n[1,3]', {'rack01n1', 'rack01n3', 'rack02n1', 'rack02n3'}),
    ('node2,node[4-5]', {'node2', 'node4', 'node5'}), ('(null)', set()),
])
def test_hostlist(value, expected):
    assert run_audit.hostnames(value) == expected


def test_script_policy_ignores_late_directives_and_detects_changed_hash(tmp_path):
    path = tmp_path/'job.slurm'
    path.write_text('#!/bin/bash\n#SBATCH --exclude=node1\necho start\n#SBATCH --exclusive\n')
    assert run_audit.script_policy(path, ['node1'], {})['status'] == 'fail'
    path.write_text('#!/bin/bash\n#SBATCH --exclude=node1\n#SBATCH --exclusive\n')
    record = {'job_sha256': file_sha256(path)}
    assert run_audit.script_policy(path, ['node1'], record)['status'] == 'pass'
    path.write_text(path.read_text()+'echo changed\n')
    assert 'script differs' in run_audit.script_policy(path, ['node1'], record)['problems'][0]


def test_scratch_is_remote_evidence_not_login_probe(tmp_path, monkeypatch):
    from kinbot import qc_scratch
    monkeypatch.setattr(qc_scratch, 'allocate', lambda *a, **k: pytest.fail('no scratch allocation'))
    monkeypatch.setattr(qc_scratch, 'filesystem', lambda *a, **k: pytest.fail('no login mount probe'))
    runtime = dict(scratch_root='/scratch/alice', scratch_directory='/scratch/alice/job',
                   native_work_directory='/scratch/alice/job/work', scratch_host='node2',
                   scratch_io_checked=True, scratch_mount={'filesystem': 'ext4'})
    path = tmp_path/'energy.runtime.json'; atomic_json(path, runtime)
    policy = run_audit.scratch_policy(tmp_path, home='/home/alice', account='alice', nodes={'node2'})
    assert policy['status'] == 'preferred_local'
    assert run_audit.scratch_policy(tmp_path, home='/home/alice', account='alice', nodes={'node3'})['status'] == 'unknown'
    runtime.update(scratch_root='/home/alice', scratch_directory='/home/alice/job',
                   native_work_directory='/home/alice/job/work')
    atomic_json(path, runtime)
    assert run_audit.scratch_policy(tmp_path, home='/home/alice', account='alice')['status'] == 'fail'
    runtime.update(scratch_root='/tmp', scratch_directory='/tmp/job', native_work_directory='/tmp/job/work')
    atomic_json(path, runtime)
    assert run_audit.scratch_policy(tmp_path, home='/home/alice', account='alice')['status'] == 'alternate_local'


def test_conditional_fallback_is_not_announced_as_next(run):
    spec = json.loads((run.graph/'workflow.json').read_text())
    task = molpro_task('fallback', 'rhf\n', geometry_from='seed')
    task['fallback_on'] = ['opt']
    spec['tasks'].append(task)
    atomic_json(run.graph/'workflow.json', spec)
    run.state['spec_sha256'] = file_sha256(run.graph/'workflow.json')
    atomic_json(run.graph/'state.json', run.state)
    report = run_audit.audit(run.root, scheduler=([], None), details={})
    row = report['cases'][0]['tasks'][-1]
    assert row['status'] == 'conditional_fallback' and row['fallback_on'] == {'opt': 'unstarted'}


def test_output_cannot_change_run_and_complete_inventory_is_saved(run, tmp_path, monkeypatch):
    before = snapshot(run.root)
    with pytest.raises(SystemExit):
        run_audit.main([str(run.root), '--output', str(run.root)])
    assert snapshot(run.root) == before
    output = tmp_path/'report'; output.mkdir()
    monkeypatch.setattr('kinbot.run_monitor.queue', lambda: ([], None))
    run_audit.main([str(run.root), '--output', str(output), '--limit', '0'])
    assert snapshot(run.root) == before
    assert 'seed' in (output/'audit.txt').read_text()
    saved = json.loads((output/'audit.json').read_text())
    assert saved['recovery_inventory'][0]['completed'] == ['seed']
    with pytest.raises(SystemExit):
        run_audit.main([str(run.root), '--output', str(output)])


def test_other_live_allocations_not_hidden(run):
    job = dict(job_id='88', directory=str(run.worker/'rotdPy'), state='PENDING', nodes='(null)')
    report = run_audit.audit(run.root, scheduler=([job], None), details={})
    assert report['cases'][0]['other_allocations'][0]['job_id'] == '88'
    assert report['cases'][0]['other_allocations'][0]['scheduler_policy']['status'] == 'unknown'
    assert 'Other live allocations=1' in run_audit.render(report)


def test_missing_legacy_runtime_evidence_is_unknown(tmp_path):
    atomic_json(tmp_path/'old.runtime.json', {'scratch_directory': '/scratch/alice/old'})
    assert run_audit.scratch_policy(tmp_path, home='/home/alice', account='alice')['status'] == 'unknown'


def test_wrong_owner_and_disabled_worker_policy_are_visible(run):
    run.state['tasks']['seed']['job_id'] = 'old'
    atomic_json(run.graph/'state.json', run.state)
    atomic_json(run.worker/'123.json', {**run.params, 'queue_exclusive': 0})
    job = dict(job_id='99', directory=str(run.graph/'tasks/seed'), state='RUNNING', nodes='node2')
    report = run_audit.audit(run.root, scheduler=([job], None), details={})
    assert report['cases'][0]['tasks'][0]['status'] == 'ownership_conflict'
    assert report['cases'][0]['policy_problems']


def test_pinned_policy_compares_recorded_source_without_importing(run):
    import shutil
    source = Path(run_audit.__file__).parent
    pin = run.root/'.salvage_code.test'
    names = ['qc_scratch.py', 'qc.py', 'anl/runtime.py', 'anl/dispatch.py',
             'ase_modules/calculators/molpro.py', 'ase_modules/calculators/gaussian.py']
    hashes = {}
    for name in names:
        dest = pin/'kinbot'/name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source/name, dest)
        hashes['kinbot/' + name] = file_sha256(dest)
    atomic_json(pin/'manifest.json', hashes)
    atomic_json(run.worker.parent/'launch.json', {'code': str(pin)})
    assert run_audit.pinned_policy(run.root, run.worker.parent)['status'] == 'matches_audit_code'
    (pin/'kinbot/qc_scratch.py').write_text('raise RuntimeError("do not execute")')
    assert run_audit.pinned_policy(run.root, run.worker.parent)['status'] == 'differs'
