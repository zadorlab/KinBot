"""Saved-results migration preserves accepted QC and scopes every mutation."""
import importlib.util
import json
from pathlib import Path

import pytest

from kinbot.anl import dispatch
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec
from tests.test_anl_dispatch import _mock_execution
from tests.test_anl_results import _GAUSSIAN_VPT2


@pytest.fixture
def repair(monkeypatch):
    package = Path(__file__).parents[1] / 'examples/anl/tools'
    monkeypatch.syspath_prepend(str(package))
    spec = importlib.util.spec_from_file_location('resume_test', package / 'resume_downstream.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def graph_fixture(root):
    spec = dispatch_spec()
    spec['tasks'] = [t for t in spec['tasks'] if t['id'] in
                     ('l2_geometry', 'l3_geometry', 'gaussian_vpt2')]
    spec['tasks'][-1]['input_template'] = spec['tasks'][-1]['input_template'].replace(
        'NoSymm', 'Symmetry=(PG=C1)')
    graph = dispatch.prepare(spec, root / 'composite')
    _mock_execution(graph, 'l2_geometry', geometry='final.xyz')
    dispatch.advance(graph)
    _mock_execution(graph, 'l3_geometry', geometry='final.xyz')
    dispatch.advance(graph)
    return graph


def fail_vpt2(graph, native):
    directory = graph / 'tasks/gaussian_vpt2'
    (directory / 'vpt2.log').write_text(native)
    (directory / 'vpt2.err').write_text('')
    task = json.loads((directory / 'task.json').read_text())
    atomic_json(directory / 'execution.json', {
        'schema': 1, 'task_id': 'gaussian_vpt2', 'status': 'failed',
        'geometry_sha256': task['geometry_sha256'], 'error': 'Old parser failure',
        'traceback': 'parse_result(native)',
    })
    state = json.loads((graph / 'state.json').read_text())
    state['tasks']['gaussian_vpt2']['status'] = 'failed'
    atomic_json(graph / 'state.json', state)


def test_native_success_reparsed_without_qc_and_history_external(repair, tmp_path):
    root = tmp_path / 'run'
    graph = graph_fixture(root)
    native = _GAUSSIAN_VPT2.replace('B3LYP', 'B2PLYP').replace(
        'NoSymm', 'Symmetry=(PG=C1) EmpiricalDispersion=GD3BJ')
    native = native.replace(' Fundamental Bands', ' Full point group C*V\n Fundamental Bands')
    fail_vpt2(graph, native)
    accepted = (graph / 'tasks/l3_geometry/execution.json').read_bytes()
    archive = repair.Archive(root)
    result = repair.repair_graph(graph, archive)
    assert result['complete_tasks_preserved'] == 3
    assert (graph / 'tasks/l3_geometry/execution.json').read_bytes() == accepted
    history = graph / 'tasks/gaussian_vpt2/execution.failed.json'
    assert history.is_symlink() and history.resolve().is_relative_to(archive.path)
    assert dispatch.advance(graph)['tasks']['gaussian_vpt2']['status'] == 'complete'


def test_framework_retry_preserves_geometry_and_finished_siblings(repair, tmp_path):
    root = tmp_path / 'run'
    graph = graph_fixture(root)
    fail_vpt2(graph, 'ERROR: Inconsistency found in framework group definition:\n'
              'New: CS - Old: C1\nError termination via Lnk1e in l717.exe\n')
    before = (graph / 'tasks/gaussian_vpt2/geometry.xyz').read_bytes()
    archive = repair.Archive(root)
    result = repair.repair_graph(graph, archive)
    assert result['recovered_vpt2'] == ['gaussian_vpt2']
    assert (graph / 'tasks/gaussian_vpt2/geometry.xyz').read_bytes() == before
    assert 'Symmetry=None' in (graph / 'tasks/gaussian_vpt2/vpt2.com').read_text()
    assert not (graph / 'tasks/gaussian_vpt2/execution.json').exists()
    assert not (graph / 'attempts').exists()
    assert list(archive.path.rglob('vpt2.log'))


def test_unknown_failure_stops_instead_of_erasing_evidence(repair, tmp_path):
    root = tmp_path / 'run'
    graph = graph_fixture(root)
    fail_vpt2(graph, 'Error termination via Lnk1e: unrelated failure\n')
    with pytest.raises(RuntimeError, match='not the framework-group'):
        repair.repair_graph(graph, repair.Archive(root))
    assert (graph / 'tasks/gaussian_vpt2/vpt2.log').is_file()


def test_unstarted_template_migrates_without_rerunning_geometries(repair, tmp_path):
    root = tmp_path / 'run'
    graph = graph_fixture(root)
    result = repair.repair_graph(graph, repair.Archive(root))
    assert result['updated_unstarted_vpt2'] == ['gaussian_vpt2']
    assert result['complete_tasks_preserved'] == 2
    assert dispatch.advance(graph)['tasks']['gaussian_vpt2']['status'] == 'staged'


def test_inactive_migration_retains_db_and_rejected_marker_and_updates_only_levels(repair, tmp_path, monkeypatch):
    import run_examples
    from kinbot import run_cleanup
    root = run_examples.prepare(tmp_path / 'run')
    # Simulate the old prepared input, before measured molecular levels.
    for case in ('ch3s', 'ch3ss'):
        for name in ('input.json', 'resolved_parameters.json'):
            path = root / case / name
            data = json.loads(path.read_text())
            data['electronic_state_models'] = [m for m in data['electronic_state_models']
                                               if m['smiles'] != '[SH]']
            for m in data['electronic_state_models']:
                if m['smiles'] == '[CH3][S]':
                    m.pop('levels', None)
            atomic_json(path, data)
    record = json.loads((root / 'prepared.json').read_text())
    record['sha256'] = {n: file_sha256(root/n) for n in record['sha256']}
    atomic_json(root / 'prepared.json', record)
    worker = root / 'ch3s/well'
    worker.mkdir()
    (worker / 'kinbot.db').write_text('accepted and rejected result cache')
    (worker / 'failed_ts.exitcode').write_text('100 1\n')
    (worker / 'failed_ts.log').write_text('native failure')
    (worker / 'kinbot.log').write_text('old controller failure')
    (worker / 'accepted.log').write_text('Normal termination')
    (root / 'ch3s/launch.json').write_text('{}')
    vrc = worker / 'vrctst/molpro_scan/point'
    vrc.mkdir(parents=True)
    atomic_json(vrc / 'state.json', {'tasks': {'relaxed_scan': {'status': 'failed'}}})
    monkeypatch.setattr(run_cleanup, 'inspect', lambda p: {'jobs': [], 'processes': []})
    original = json.loads((root / 'ch3s/input.json').read_text())
    result = repair.repair_before_resume(root, run_examples.verify(root))
    assert Path(result['archive']).is_dir()
    assert not (worker / 'kinbot.log').exists()
    assert not (worker / 'failed_ts.log').exists()
    assert not (worker / 'vrctst').exists()
    assert (worker / 'failed_ts.exitcode').read_text() == '100 1\n'
    assert (worker / 'accepted.log').read_text() == 'Normal termination'
    assert (worker / 'kinbot.db').read_text() == 'accepted and rejected result cache'
    assert not (root / 'ch3s/launch.json').exists()
    new = json.loads((root / 'ch3s/input.json').read_text())
    assert new.pop('electronic_state_models') != original.pop('electronic_state_models')
    assert new == original
    run_examples.verify(root)


def test_active_run_refused_before_archiving(repair, tmp_path, monkeypatch):
    from kinbot import run_cleanup
    root = tmp_path / 'run'
    case = root / 'ch3s'
    case.mkdir(parents=True)
    monkeypatch.setattr(run_cleanup, 'inspect', lambda p: {'jobs': [{'id': '1'}], 'processes': []})
    with pytest.raises(RuntimeError, match='inactive'):
        repair.repair_before_resume(root, {'inputs': [str(case / 'input.json')]})
    assert not list(tmp_path.glob('*.repair_archive.*'))
