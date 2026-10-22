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


def graph_fixture(root, *, optional=False):
    spec = dispatch_spec()
    spec['tasks'] = [t for t in spec['tasks'] if t['id'] in
                     ('l2_geometry', 'l3_geometry', 'gaussian_vpt2')]
    spec['tasks'][-1]['input_template'] = spec['tasks'][-1]['input_template'].replace(
        'NoSymm', 'Symmetry=(PG=C1)')
    if optional:
        spec['tasks'][-1].update(failure_policy='skip_component',
                                component_families=['vpt2_correction'])
    graph = dispatch.prepare(spec, root / 'composite')
    _mock_execution(graph, 'l2_geometry', geometry='final.xyz')
    dispatch.advance(graph)
    _mock_execution(graph, 'l3_geometry', geometry='final.xyz')
    dispatch.advance(graph)
    return graph


def fail_vpt2(graph, native, *, failure_phase=None):
    directory = graph / 'tasks/gaussian_vpt2'
    (directory / 'vpt2.log').write_text(native)
    (directory / 'vpt2.err').write_text('')
    task = json.loads((directory / 'task.json').read_text())
    atomic_json(directory / 'execution.json', {
        'schema': 1, 'task_id': 'gaussian_vpt2', 'status': 'failed',
        'geometry_sha256': task['geometry_sha256'], 'error': 'Old parser failure',
        'traceback': 'parse_result(native)',
        **({'failure_phase': failure_phase} if failure_phase is not None else {}),
    })
    state = json.loads((graph / 'state.json').read_text())
    state['tasks']['gaussian_vpt2']['status'] = 'failed'
    atomic_json(graph / 'state.json', state)


@pytest.mark.parametrize('unavailable', [False, True])
def test_native_success_reparsed_without_qc_and_history_external(repair, tmp_path, unavailable):
    root = tmp_path / 'run'
    graph = graph_fixture(root)
    native = _GAUSSIAN_VPT2.replace('B3LYP', 'B2PLYP').replace(
        'NoSymm', 'Symmetry=(PG=C1) EmpiricalDispersion=GD3BJ')
    native = native.replace(' Fundamental Bands', ' Full point group C*V\n Fundamental Bands')
    fail_vpt2(graph, native)
    if unavailable:
        state = json.loads((graph / 'state.json').read_text())
        state['tasks']['gaussian_vpt2'].update(
            status='unavailable', failure_phase='execution',
            unavailable_reason='old launcher failure after native success')
        atomic_json(graph / 'state.json', state)
    accepted = (graph / 'tasks/l3_geometry/execution.json').read_bytes()
    archive = repair.Archive(root)
    result = repair.repair_graph(graph, archive)
    assert result['complete_tasks_preserved'] == 3
    assert (graph / 'tasks/l3_geometry/execution.json').read_bytes() == accepted
    history = graph / 'tasks/gaussian_vpt2/execution.failed.json'
    assert history.is_symlink() and history.resolve().is_relative_to(archive.path)
    assert dispatch.advance(graph)['tasks']['gaussian_vpt2']['status'] == 'complete'
    entry = dispatch.advance(graph)['tasks']['gaussian_vpt2']
    assert 'unavailable_reason' not in entry and 'failure_phase' not in entry


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
    previous = dispatch._load(graph)[2]['tasks']['gaussian_vpt2']['previous_attempt']
    assert (Path(previous)/'vpt2.log').is_file()
    assert (graph/'.live_control_transaction.json').is_file()


def test_native_reparse_releases_derived_skips_but_not_unrelated_work(tmp_path):
    from copy import deepcopy
    spec = dispatch_spec()
    base = next(t for t in spec['tasks'] if t['id'] == 'gaussian_vpt2')
    base = deepcopy(base)
    base.pop('depends_on', None)
    base['geometry_from'] = 'l2_geometry'
    base.update(id='source', failure_policy='skip_component',
                component_families=['vpt2_correction'],
                failure_skips=['larger'])
    larger = deepcopy(base)
    larger.update(id='larger', failure_skips=[])
    unrelated = deepcopy(base)
    unrelated.update(id='unrelated', failure_skips=['other_larger'])
    other = deepcopy(larger)
    other['id'] = 'other_larger'
    child = deepcopy(larger)
    child.update(id='child', depends_on=['source'])
    geometry = next(t for t in spec['tasks'] if t['id'] == 'l2_geometry')
    spec['tasks'] = [geometry, base, larger, unrelated, other, child]
    graph = dispatch.prepare(spec, tmp_path / 'graph')
    _mock_execution(graph, 'l2_geometry', geometry='final.xyz')
    dispatch.advance(graph)
    for ident in ('source', 'unrelated'):
        task_dir = graph / 'tasks' / ident
        entry = json.loads((graph / 'state.json').read_text())['tasks'][ident]
        atomic_json(task_dir / 'execution.json', {
            'schema': 1, 'task_id': ident, 'geometry_sha256': entry['geometry_sha256'],
            'status': 'failed', 'failure_phase': 'execution', 'error': 'native failure'})
    state = dispatch.advance(graph)
    assert state['tasks']['larger']['resume_status'] == 'staged'
    assert state['tasks']['child']['status'] == 'skipped'
    # Simulate the native file becoming visible after an old launcher failure;
    # the real parser and provenance checks must accept it before skips clear.
    native = _GAUSSIAN_VPT2.replace('B3LYP', 'B2PLYP').replace(
        'NoSymm', 'NoSymm EmpiricalDispersion=GD3BJ')
    (graph / 'tasks/source/vpt2.log').write_text(native)
    (graph / 'tasks/source/vpt2.err').write_text('')
    state = dispatch.advance(graph)
    assert state['tasks']['source']['status'] == 'complete'
    assert state['tasks']['source']['automatically_reparsed']
    assert state['tasks']['larger']['status'] == 'staged'
    assert state['tasks']['child']['status'] == 'staged'
    assert state['tasks']['unrelated']['status'] == 'unavailable'
    assert state['tasks']['other_larger']['status'] == 'skipped'


def test_optional_native_framework_failure_retries_from_unavailable(repair, tmp_path):
    root = tmp_path / 'run'
    graph = graph_fixture(root, optional=True)
    fail_vpt2(graph, 'ERROR: Inconsistency found in framework group definition:\n'
              'New: CS - Old: C1\nError termination via Lnk1e in l717.exe\n',
              failure_phase='execution')
    # An older controller may already have made the component unavailable.
    # The new common dispatcher recovers it before the legacy helper runs.
    state = dispatch._load(graph)[2]
    state['tasks']['gaussian_vpt2']['status'] = 'unavailable'
    atomic_json(graph/'state.json', state)
    accepted = (graph / 'tasks/l3_geometry/execution.json').read_bytes()
    geometry = (graph / 'tasks/gaussian_vpt2/geometry.xyz').read_bytes()
    result = repair.repair_graph(graph, repair.Archive(root))
    assert result['recovered_vpt2'] == ['gaussian_vpt2']
    assert dispatch.advance(graph)['tasks']['gaussian_vpt2']['status'] == 'staged'
    assert (graph / 'tasks/gaussian_vpt2/geometry.xyz').read_bytes() == geometry
    assert (graph / 'tasks/l3_geometry/execution.json').read_bytes() == accepted
    assert 'Symmetry=None' in (graph / 'tasks/gaussian_vpt2/vpt2.com').read_text()
    # Rerunning an interrupted migration does not create another QC attempt.
    again = repair.repair_graph(graph, repair.Archive(root))
    assert again['recovered_vpt2'] == []
    assert dispatch.advance(graph)['tasks']['gaussian_vpt2']['attempt'] == 2


def test_preflight_unavailable_component_is_not_reenabled(repair, tmp_path):
    root = tmp_path / 'run'
    graph = graph_fixture(root, optional=True)
    spec = json.loads((graph / 'workflow.json').read_text())
    state = json.loads((graph / 'state.json').read_text())
    dispatch._mark_preflight_unavailable(spec, state, 'gaussian_vpt2', 'test unavailable backend')
    atomic_json(graph / 'state.json', state)
    before = (graph / 'tasks/gaussian_vpt2/vpt2.com').read_bytes()
    repair.repair_graph(graph, repair.Archive(root))
    state = dispatch.advance(graph)['tasks']['gaussian_vpt2']
    assert state['status'] == 'unavailable'
    assert state['unavailable_reason'] == 'test unavailable backend'
    assert (graph / 'tasks/gaussian_vpt2/vpt2.com').read_bytes() == before


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
                    m.pop('spin_free_reference', None)
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
    composite = worker / 'anl_automatic'
    (composite / 'profiles').mkdir(parents=True)
    (composite / 'profiles/old.json').write_text('old derived energy')
    (composite / 'network_request.json').write_text('frozen QC request')
    (worker / 'me').mkdir()
    (worker / 'me/mess.inp').write_text('old populations and energies')
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
    assert (composite / 'network_request.json').read_text() == 'frozen QC request'
    assert not (composite / 'profiles').exists()
    assert not (worker / 'me').exists()
    assert (Path(result['archive']) / 'old_postprocessing/ch3s/well/anl_automatic/profiles/old.json').read_text() == 'old derived energy'
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


def test_vrc_fragments_survive_scan_archive_and_old_archive_is_recoverable(repair, tmp_path, monkeypatch):
    from ase import Atoms
    from ase.db import connect
    from kinbot.qc import QuantumChemistry
    root = tmp_path / 'run'
    case = root / 'ch3s'
    worker = case / 'well'
    vrc = worker / 'vrctst'
    (vrc / 'molpro_scan/point').mkdir(parents=True)
    atomic_json(vrc / 'molpro_scan/point/state.json',
                {'tasks': {'relaxed_scan': {'status': 'failed'}}})
    db = connect(str(worker / 'kinbot.db'))
    db.write(Atoms('HH', positions=[[0, 0, 0], [0, 0, .74]]),
             name='vrctst/fragment_vts', data={'energy': -1., 'status': 'normal'})
    expected = {'.log': 'Normal termination of Gaussian 16\ndone\n',
                '.chk': 'accepted checkpoint', '.exitcode': '142 0\n',
                '.com': 'accepted input', '.py': 'original worker'}
    for suffix, text in expected.items():
        (vrc / ('fragment_vts' + suffix)).write_text(text)
    # Reproduce the previous helper's broad archive exactly, then a resume
    # that rewrites only the worker script and sees an old database job ID.
    old = repair.Archive(root)
    old.save(vrc, move=True, group='invalid_vrc')
    vrc.mkdir()
    (vrc / 'fragment_vts.py').write_text('regenerated worker')
    current = repair.Archive(root)
    restored = repair.restore_archived_vrc_fragments(root, case, current)
    assert len(restored) == 4
    assert (vrc / 'fragment_vts.py').read_text() == 'regenerated worker'
    for suffix in ('.log', '.chk', '.exitcode', '.com'):
        assert (vrc / ('fragment_vts' + suffix)).read_text() == expected[suffix]
    assert not (vrc / 'molpro_scan').exists()  # wrong-geometry scans not restored
    assert repair.restore_archived_vrc_fragments(root, case, current) == []
    (vrc / 'molpro_scan/point').mkdir(parents=True)
    atomic_json(vrc / 'molpro_scan/point/state.json',
                {'tasks': {'relaxed_scan': {'status': 'failed'}}})
    moved = repair.archive_invalid_vrc(vrc, current)
    assert len(moved) == 1 and not (vrc / 'molpro_scan').exists()
    assert (vrc / 'fragment_vts.chk').read_text() == expected['.chk']
    # The ordinary QC cache reader must accept the restored generation
    # without querying Slurm, polling 60 seconds, or submitting another job.
    monkeypatch.chdir(worker)
    qc = QuantumChemistry.__new__(QuantumChemistry)
    qc.db, qc.qc, qc.queuing, qc.use_sella = db, 'gauss', 'slurm', True
    qc.job_ids = {'vrctst/fragment_vts': '142'}
    qc.job_generation_tagged = {'vrctst/fragment_vts': True}
    monkeypatch.setattr('kinbot.qc._slurm_job_active',
                        lambda *args: pytest.fail('completed fragment queried in Slurm'))
    assert qc.check_qc('vrctst/fragment_vts') == 'normal'


def test_vrc_restore_refuses_to_overwrite_new_native_output(repair, tmp_path):
    from ase import Atoms
    from ase.db import connect
    root = tmp_path / 'run'
    case = root / 'ch3ss'
    vrc = case / 'well/vrctst'
    vrc.mkdir(parents=True)
    connect(str(vrc.parent / 'kinbot.db')).write(
        Atoms('H'), name='vrctst/fragment_vts', data={'status': 'normal'})
    (vrc / 'fragment_vts.exitcode').write_text('142 0\n')
    (vrc / 'fragment_vts.log').write_text('Normal termination of Gaussian\n')
    old = repair.Archive(root)
    old.save(vrc, move=True, group='invalid_vrc')
    vrc.mkdir()
    (vrc / 'fragment_vts.log').write_text('new calculation evidence')
    with pytest.raises(ValueError, match='differs from archived'):
        repair.restore_archived_vrc_fragments(root, case, repair.Archive(root))
    assert (vrc / 'fragment_vts.log').read_text() == 'new calculation evidence'


def test_interrupted_two_case_repair_finishes_after_inputs_already_updated(repair, tmp_path, monkeypatch):
    import run_examples
    from kinbot import run_cleanup
    root = run_examples.prepare(tmp_path / 'run')
    for case in ('ch3s', 'ch3ss'):
        path = root / case / 'input.json'
        data = json.loads(path.read_text())
        data['electronic_state_models'] = [m for m in data['electronic_state_models']
                                           if m['smiles'] != '[SH]']
        atomic_json(path, data)
        base = root / case / 'well/anl_automatic'
        (base / 'profiles').mkdir(parents=True)
        (base / 'profiles/old.json').write_text('old derived energy')
        (base / 'network_request.json').write_text('frozen QC request')
    graph = graph_fixture(root / 'ch3s/well/anl_automatic/states/state_h2s', optional=True)
    fail_vpt2(graph, 'ERROR: Inconsistency found in framework group definition:\n'
              'New: CS - Old: C1\nError termination via Lnk1e in l717.exe\n',
              failure_phase='execution')
    accepted = (graph / 'tasks/l3_geometry/execution.json').read_bytes()
    record = json.loads((root / 'prepared.json').read_text())
    record['sha256'] = {n: file_sha256(root/n) for n in record['sha256']}
    atomic_json(root / 'prepared.json', record)
    monkeypatch.setattr(run_cleanup, 'inspect', lambda p: {'jobs': [], 'processes': []})
    real_repair = repair.repair_graph

    def interrupted(*args):
        # Fail at the same point as the reported helper, after both input
        # files changed but before CH3SS postprocessing was archived.
        raise RuntimeError('simulated unavailable migration interruption')

    monkeypatch.setattr(repair, 'repair_graph', interrupted)
    with pytest.raises(RuntimeError, match='simulated unavailable'):
        repair.repair_before_resume(root, run_examples.verify(root))
    assert not (root / 'ch3s/well/anl_automatic/profiles').exists()
    assert (root / 'ch3ss/well/anl_automatic/profiles').exists()
    monkeypatch.setattr(repair, 'repair_graph', real_repair)
    result = repair.repair_before_resume(root, run_examples.verify(root))
    assert result['electronic_inputs_updated'] == []
    assert not (root / 'ch3ss/well/anl_automatic/profiles').exists()
    assert (graph / 'tasks/l3_geometry/execution.json').read_bytes() == accepted
    assert dispatch.advance(graph)['tasks']['gaussian_vpt2']['attempt'] == 2
    assert (root / 'ch3ss/well/anl_automatic/network_request.json').read_text() == 'frozen QC request'
    assert len(list(tmp_path.glob('run.repair_archive.*'))) == 2
    run_examples.verify(root)
