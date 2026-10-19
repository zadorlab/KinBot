"""Outage recovery must preserve owned work and never fabricate completion."""
from copy import deepcopy
import fcntl
import json
from pathlib import Path
import socket
from types import SimpleNamespace

import pytest

from kinbot import run_cleanup, run_recover, run_refresh, run_salvage
from kinbot.anl import dispatch
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task
from tests.test_anl_dispatch import _mock_execution


def ledger(graph, ident, **values):
    state = json.loads((graph/'state.json').read_text())
    state['tasks'][ident].update(values)
    atomic_json(graph/'state.json', state)
    return state['tasks'][ident]


def success(graph, ident, *, geometry=None):
    if geometry is None:
        (graph/'tasks'/ident/f'{ident}.out').write_text('Molpro calculation terminated\n')
    _mock_execution(graph, ident, geometry=geometry)


def file_snapshot(root):
    return {str(p.relative_to(root)): file_sha256(p) for p in root.rglob('*') if p.is_file()}


@pytest.fixture
def outage(tmp_path, monkeypatch):
    root = tmp_path/'run'; root.mkdir()
    pin = root/'.salvage_code.new'; pin.mkdir()
    inputs, graphs, jobs, processes, launches = [], [], [], [], []
    chemistry = {'ipea': 0.25, 'shift': 0.1, 'active_orbitals': [17, 18, 19, 20, 21],
                 'queue_exclusive': 0, 'rotdpy_exclusive': 0}
    for case in ('alpha', 'beta'):
        directory = root/case; directory.mkdir()
        path = directory/'input.json'; atomic_json(path, chemistry)
        inputs.append(str(path))
        worker = directory/'123'; worker.mkdir()
        (worker/'.kinbot_run.json').write_text('{}')
        (worker/'kinbot.db').write_bytes(b'accepted L1/L2')
        atomic_json(worker/'kinbot_status.json', {'status': 'failed', 'host': socket.gethostname()})
        atomic_json(worker/'.kinbot_post_reaction_branches.json', {'branches': {'cbh_anl': {'status': 'failed'}}})
        atomic_json(directory/'launch.json', {'pid': 123, 'host': socket.gethostname()})
        (directory/'driver.log').write_text('old failed controller')
        spec = dispatch_spec()
        spec['tasks'] = [molpro_task('energy', 'rhf\n', geometry_from='initial')]
        graph = worker/'anl_automatic/states/test/composite'
        dispatch.prepare(spec, graph)
        graphs.append(graph)
    success(graphs[0], 'energy')
    ledger(graphs[0], 'energy', status='submitting')  # Lost sbatch response.
    ledger(graphs[1], 'energy', status='submitted', job_id='555')
    jobs.append(dict(job_id='555', directory=str(graphs[1]/'tasks/energy'), state='RUNNING'))
    atomic_json(root/'prepared.json', dict(inputs=inputs, sha256={
        str(Path(p).relative_to(root)): file_sha256(p) for p in inputs}))
    atomic_json(root/'salvage.json', dict(status='launched', inputs=inputs, code='old-code'))
    monkeypatch.setattr(run_cleanup, '_jobs', lambda root: deepcopy(jobs))
    monkeypatch.setattr(run_cleanup, '_processes', lambda root: deepcopy(processes))
    monkeypatch.setattr(run_refresh, '_preflight', lambda root: None)
    monkeypatch.setattr(run_salvage, '_pin_code', lambda *a, **kw: pin)
    monkeypatch.setattr(run_salvage, '_preflight_worker_routes', lambda root: None)
    monkeypatch.setattr('kinbot.parameters.Parameters', lambda path: SimpleNamespace(par={}))
    monkeypatch.setattr('kinbot.anl.automatic.validate_automatic_environment', lambda par: None)
    def popen(command, **kwargs):
        launches.append((command, kwargs))
        return SimpleNamespace(pid=900+len(launches))
    monkeypatch.setattr(run_salvage.subprocess, 'Popen', popen)
    return SimpleNamespace(root=root, graphs=graphs, jobs=jobs, processes=processes,
                           inputs=inputs, launches=launches, chemistry=chemistry)


def test_read_only_plan_recovers_lost_response_and_preserves_live_job(outage):
    t = outage; before = file_snapshot(t.root)
    result = run_recover.recover(t.root)
    assert result['phase'] == 'reviewed'
    assert result['graphs'][0]['complete'] == ['energy']
    assert result['graphs'][1]['active'] == {'energy': '555'}
    assert not any(g['reset'] for g in result['graphs'])
    assert file_snapshot(t.root) == before


def test_apply_keeps_native_work_and_chemistry_then_launches_once(outage):
    t = outage
    native = [file_snapshot(g/'tasks') for g in t.graphs]
    result = run_recover.recover(t.root, apply=True, launch=True)
    assert result['phase'] == 'launched' and len(t.launches) == 2
    assert [file_snapshot(g/'tasks') for g in t.graphs] == native
    assert dispatch._load(t.graphs[0])[2]['tasks']['energy']['status'] == 'complete'
    assert dispatch._load(t.graphs[1])[2]['tasks']['energy']['job_id'] == '555'
    for filename in t.inputs:
        values = json.loads(Path(filename).read_text())
        assert values['queue_exclusive'] == values['rotdpy_exclusive'] == 1
        assert {k: values[k] for k in ('ipea', 'shift', 'active_orbitals')} == {
            k: t.chemistry[k] for k in ('ipea', 'shift', 'active_orbitals')}
        assert (Path(filename).parent/'123/kinbot.db').read_bytes() == b'accepted L1/L2'
        assert not (Path(filename).parent/'123/.kinbot_post_reaction_branches.json').exists()
    assert all(kw['env']['PYTHONPATH'].startswith(result['code']) for _, kw in t.launches)
    assert (Path(result['archive'])/'alpha/driver.log').read_text() == 'old failed controller'
    assert run_recover.recover(t.root, apply=True, launch=True) == result
    assert len(t.launches) == 2 and len(t.jobs) == 1


@pytest.mark.parametrize('problem', ['controller', 'scheduler', 'unowned', 'duplicate',
                                     'wrong_owner', 'changed_live_input', 'preflight', 'new_job'])
def test_unsafe_inventory_fails_before_archive(outage, monkeypatch, problem):
    t = outage
    if problem == 'controller':
        t.processes.append({'pid': 100})
    elif problem == 'scheduler':
        monkeypatch.setattr(run_cleanup, '_jobs', lambda root: (_ for _ in ()).throw(RuntimeError('scheduler unavailable')))
    elif problem == 'unowned':
        t.jobs.append(dict(job_id='556', directory=str(t.root/'unknown'), state='RUNNING'))
    elif problem == 'duplicate':
        t.jobs.append({**t.jobs[0], 'job_id': '556'})
    elif problem == 'wrong_owner':
        ledger(t.graphs[1], 'energy', job_id='556')
    elif problem == 'changed_live_input':
        (t.graphs[1]/'tasks/energy/energy.inp').write_text('changed')
    elif problem == 'preflight':
        monkeypatch.setattr(run_refresh, '_preflight', lambda root: (_ for _ in ()).throw(RuntimeError('bad environment')))
    elif problem == 'new_job':
        count = 0
        def inventory(root):
            nonlocal count
            count += 1
            return deepcopy(t.jobs) + ([dict(job_id='556', directory=str(t.root/'new'))] if count > 1 else [])
        monkeypatch.setattr(run_cleanup, '_jobs', inventory)
    before = file_snapshot(t.root)
    with pytest.raises((RuntimeError, ValueError)):
        run_recover.recover(t.root, apply=True)
    assert file_snapshot(t.root) == before
    assert not t.launches and not list(t.root.parent.glob('*.outage_archive.*'))


def test_corrupt_success_is_archived_and_never_retained(outage):
    t = outage; graph = t.graphs[0]
    (graph/'tasks/energy/energy.out').write_text('truncated')
    result = run_recover.recover(t.root, apply=True)
    assert result['graphs'][0]['complete'] == []
    assert result['graphs'][0]['reset'] == ['energy']
    assert not (graph/'tasks/energy').exists()
    archived = Path(result['archive'])/'graphs'/graph.relative_to(t.root)/'tasks/energy/energy.out'
    assert archived.read_text() == 'truncated'


def test_live_vpt2_recovery_input_is_not_rewritten(tmp_path):
    spec = dispatch_spec()
    task = spec['tasks'][-1]
    task.pop('depends_on')
    task['input_template'] = task['input_template'].replace('NoSymm', 'Symmetry=(PG=C1)')
    task['recovery'] = {'framework_symmetry': 'C1'}
    spec['tasks'] = [spec['tasks'][0], task]
    graph = tmp_path/'graph'; dispatch.prepare(spec, graph)
    success(graph, 'l2_geometry', geometry='final.xyz'); dispatch.advance(graph)
    ledger(graph, 'gaussian_vpt2', status='submitted', job_id='555', attempt=2)
    before = file_snapshot(graph/'tasks')
    item = run_recover.graph_plan(graph, [dict(job_id='555', directory=str(graph/'tasks/gaussian_vpt2'))])
    run_recover._apply_graph(item, tmp_path, tmp_path/'archive')
    assert item['active'] == {'gaussian_vpt2': '555'}
    assert item['spec']['tasks'] == spec['tasks']
    assert file_snapshot(graph/'tasks') == before


def test_live_submission_with_lost_response_is_adopted_without_resubmission(outage):
    graph = outage.graphs[1]
    entry = ledger(graph, 'energy', status='submitting', job_id=None)
    item = run_recover.graph_plan(graph, outage.jobs)
    assert item['active'] == {'energy': '555'}
    assert item['state']['tasks']['energy']['status'] == 'submitted'
    assert item['state']['tasks']['energy']['task_sha256'] == entry['task_sha256']
    assert not item['reset']


def test_missing_outcome_is_retry_not_success(outage):
    graph = outage.graphs[0]
    (graph/'tasks/energy/execution.json').unlink()
    item = run_recover.graph_plan(graph, [])
    assert item['reset'] == ['energy'] and not item['complete']


@pytest.mark.parametrize('folder', ['molpro_pivots/fragment', 'molpro/dispatch', 'molpro/dispatch_123'])
def test_inventory_includes_later_vrc_graphs(outage, folder):
    t = outage
    graph = t.root/'alpha/123/vrctst'/folder
    spec = dispatch_spec()
    spec['tasks'] = [molpro_task('pivot', 'rhf\n', geometry_from='initial')]
    dispatch.prepare(spec, graph)
    ledger(graph, 'pivot', status='submitted', job_id='556')
    t.jobs.append(dict(job_id='556', directory=str(graph/'tasks/pivot'), state='RUNNING'))
    result = run_recover.plan(t.root)
    assert {j for g in result['graphs'] for j in g['active'].values()} == {'555', '556'}


def test_vrc_scan_after_archive_reaches_dispatcher_for_restaging(tmp_path, monkeypatch):
    import numpy as np
    from kinbot.vrc_tst_scan import VTS
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr('kinbot.anl.site._partitions', lambda: [
        dict(name='short', cores=8, memory_mb=32000, seconds=3600, default=True)])
    geom = np.asarray([[0., 0., 0.], [0., 0., 2.5]])
    reaction = SimpleNamespace(scan_coo=[0, 1], equiv=[[0], [1]], products=[],
                 species=SimpleNamespace(atom=['H', 'H'], charge=0, mult=1))
    par = dict(vrc_tst_scan_method='caspt2(2,2)', vrc_tst_scan_basis='vdz',
               vrc_tst_walltime='01:00:00', single_point_ppn=4,
               vrc_tst_min_stack_mw=50, queue_name='short')
    original = dispatch.prepare
    def archived(spec, graph):
        original(spec, graph)
        ledger(graph, 'relaxed_scan', status='failed', failure_phase='execution')
        item = run_recover.graph_plan(graph, [])
        run_recover._apply_graph(item, tmp_path, tmp_path/'archive')
    class Restaged(Exception):
        pass
    def drive(graph, tasks, **kwargs):
        state = dispatch.advance(graph)
        assert tasks == ('relaxed_scan',)
        assert state['tasks']['relaxed_scan']['status'] == 'staged'
        assert state['tasks']['relaxed_scan']['attempt'] == 2
        raise Restaged
    monkeypatch.setattr(dispatch, 'prepare', archived)
    monkeypatch.setattr(dispatch, 'preflight_until_scheduler_ready', lambda *a, **k: None)
    monkeypatch.setattr('kinbot.vrc_tst_scan._drive_vrc_dispatch', drive)
    with pytest.raises(Restaged):
        VTS(SimpleNamespace(), par, None)._molpro_multireference_scan_task(
            reaction, geom, 2.5, 'test_point')


def test_native_failure_marker_even_with_matching_artifacts_is_not_success(tmp_path):
    spec = dispatch_spec()
    task = molpro_task('energy', 'rhf\n', geometry_from='initial')
    task['failure_markers'] = [{'file': 'energy.out', 'contains': 'NO CONVERGENCE'}]
    spec['tasks'] = [task]
    graph = tmp_path/'graph'; dispatch.prepare(spec, graph)
    (graph/'tasks/energy/energy.out').write_text('NO CONVERGENCE\nMolpro calculation terminated\n')
    _mock_execution(graph, 'energy')
    item = run_recover.graph_plan(graph, [])
    assert item['reset'] == ['energy'] and not item['complete']


def test_wrong_geometry_single_point_is_not_salvaged(tmp_path):
    spec = dispatch_spec()
    spec['tasks'] = spec['tasks'][:2] + [molpro_task('sp', 'rhf\n')]
    graph = tmp_path/'graph'; dispatch.prepare(spec, graph)
    success(graph, 'l2_geometry', geometry='final.xyz'); dispatch.advance(graph)
    success(graph, 'l3_geometry', geometry='final.xyz'); dispatch.advance(graph)
    success(graph, 'sp'); dispatch.advance(graph)
    ledger(graph, 'sp', geometry_from='l2_geometry')
    item = run_recover.graph_plan(graph, [])
    assert item['reset'] == ['sp']
    assert item['complete'] == ['l2_geometry', 'l3_geometry']


def test_failed_geometry_defers_descendants_but_keeps_accepted_seed(tmp_path):
    spec = dispatch_spec()
    spec['tasks'] = spec['tasks'][:2] + [molpro_task('sp', 'rhf\n')]
    graph = tmp_path/'composite'; dispatch.prepare(spec, graph)
    success(graph, 'l2_geometry', geometry='final.xyz')
    dispatch.advance(graph)
    success(graph, 'l3_geometry', geometry='final.xyz')
    dispatch.advance(graph)
    success(graph, 'sp'); dispatch.advance(graph)
    entry = ledger(graph, 'l3_geometry', status='submitted', job_id='old', attempt=2)
    atomic_json(graph/'tasks/l3_geometry/execution.json', dict(schema=1, task_id='l3_geometry',
                geometry_sha256=entry['geometry_sha256'], status='failed',
                failure_phase='execution', error='copy-back failed'))
    item = run_recover.graph_plan(graph, [])
    assert item['complete'] == ['l2_geometry']
    assert set(item['reset']) == {'l3_geometry', 'sp'}
    seed = file_snapshot(graph/'tasks/l2_geometry')
    archive = tmp_path/'archive'
    run_recover._apply_graph(item, tmp_path, archive)
    run_recover._apply_graph(item, tmp_path, archive)  # Resumable graph commit.
    assert file_snapshot(graph/'tasks/l2_geometry') == seed
    state = dispatch.advance(graph)
    assert state['tasks']['l3_geometry']['attempt'] == 3
    assert state['tasks']['l3_geometry']['status'] == 'staged'
    assert 'sp' not in state['tasks'] and not (graph/'tasks/sp').exists()
    success(graph, 'l3_geometry', geometry='final.xyz')
    state = dispatch.advance(graph)
    assert state['tasks']['sp']['status'] == 'staged'
    assert state['tasks']['sp']['attempt'] == 2
    assert Path(state['tasks']['sp']['previous_attempt']).is_dir()
    assert 'sp' not in state['pending_retries']


@pytest.mark.parametrize('interrupt_at', ['graph_write', 'input_hash_update', 'archive_logs'])
def test_interrupted_migration_resumes_without_overwriting_evidence(outage, monkeypatch, interrupt_at):
    t = outage
    # Unsubmitted work needs restaging, while the other graph remains live.
    (t.graphs[0]/'tasks/energy/execution.json').unlink()
    ledger(t.graphs[0], 'energy', status='staged')
    native = file_snapshot(t.graphs[1]/'tasks')
    original_atomic = run_recover.atomic_json
    original_archive = run_refresh._archive_controllers
    raised = False
    def write(path, value):
        nonlocal raised
        hit = ((interrupt_at == 'graph_write' and Path(path) == t.graphs[0]/'state.json') or
               (interrupt_at == 'input_hash_update' and Path(path) == t.root/'prepared.json'))
        if hit and not raised:
            raised = True
            raise OSError('simulated outage')
        return original_atomic(path, value)
    def archive(root, saved):
        nonlocal raised
        original_archive(root, saved)
        if interrupt_at == 'archive_logs' and not raised:
            raised = True
            raise OSError('simulated outage')
    monkeypatch.setattr(run_recover, 'atomic_json', write)
    monkeypatch.setattr(run_refresh, '_archive_controllers', archive)
    with pytest.raises(OSError, match='simulated outage'):
        run_recover.recover(t.root, apply=True)
    result = run_recover.recover(t.root, apply=True, launch=True)
    assert result['phase'] == 'launched' and len(t.launches) == 2
    assert file_snapshot(t.graphs[1]/'tasks') == native
    assert len(list(t.root.parent.glob('*.outage_archive.*'))) == 1
    archived = Path(result['archive'])/'graphs'/t.graphs[0].relative_to(t.root)/'tasks/energy'
    assert archived.is_dir() and not (t.graphs[0]/'tasks/energy').exists()


def test_cli_releases_worker_locks_and_partial_launch_can_resume(outage, monkeypatch):
    t = outage
    first = run_salvage.subprocess.Popen
    calls = 0
    def launch(command, **kw):
        nonlocal calls
        # The future worker must be able to acquire its controller lock.
        with (Path(kw['cwd'])/'123/.kinbot_run.json').open() as handle:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        calls += 1
        if calls == 2:
            raise OSError('second controller launch interrupted')
        return first(command, **kw)
    monkeypatch.setattr(run_salvage.subprocess, 'Popen', launch)
    with pytest.raises(OSError, match='launch interrupted'):
        run_recover.main([str(t.root), '--apply', '--launch'])
    assert len(t.launches) == 1
    # The first controller now owns its lock while we resume the second.
    with (t.root/'alpha/123/.kinbot_run.json').open() as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        run_recover.main([str(t.root), '--apply', '--launch'])
    assert len(t.launches) == 2
    run_recover.main([str(t.root), '--apply', '--launch'])
    assert len(t.launches) == 2
