from copy import deepcopy
import json
from pathlib import Path

import pytest

from kinbot import vrc_graph, vrc_deduplicate
from kinbot.anl import dispatch
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task
from tests.test_anl_dispatch import _mock_execution


def spec(name='asymptote', walltime='24:00:00'):
    value = dispatch_spec()
    value['name'] = name
    value['tasks'] = [molpro_task('point', 'rhf\n', geometry_from='initial', walltime=walltime)]
    value['intent'] = {'vrc_scan_fingerprint': walltime,
                       'vrc_scan_distance_angstrom': 20.,
                       'cap_walltime_to_explicit_partition': True}
    return value


def graph(path, value, status='complete', job='10'):
    dispatch.prepare(value, path)
    state = dispatch._load(path)[2]
    state['tasks']['point'].update(status=status, job_id=job)
    if status == 'complete':
        _mock_execution(path, 'point')
    atomic_json(path/'state.json', state)
    return path


@pytest.mark.parametrize('status', ['complete', 'submitted', 'staged', 'failed'])
def test_scheduling_change_adopts_legacy_path_and_preserves_job_files(tmp_path, status):
    original = graph(tmp_path/'asymptote_legacy', spec(), status)
    before = {p: file_sha256(p) for p in original.rglob('*') if p.is_file()}
    wanted = spec(walltime='7-00:00:00')
    wanted['tasks'][0]['resources'].update(cores=8, memory_mb=32000, partition='long')
    wanted['limits']['max_nodes'] = 23
    assert vrc_graph.select_graph(tmp_path, 'asymptote', wanted) == original
    assert before == {p: file_sha256(p) for p in before}
    assert len(list(tmp_path.glob('*/workflow.json'))) == 1
    assert vrc_graph.select_graph(tmp_path, 'asymptote', wanted) == original


@pytest.mark.parametrize('change', ['geometry', 'method', 'input', 'distance', 'checkpoint'])
def test_scientific_changes_do_not_reuse_results(tmp_path, change):
    old = graph(tmp_path/'asymptote_old', spec())
    wanted = spec(walltime='7-00:00:00')
    if change == 'geometry': wanted['molecule']['positions'][0][0] += 1.e-8
    if change == 'method': wanted['tasks'][0]['input_template'] += 'uccsd(t)\n'
    if change == 'input': wanted['tasks'][0]['input_template'] += 'basis=vdz\n'
    if change == 'distance': wanted['intent']['vrc_scan_distance_angstrom'] = 21.
    if change == 'checkpoint': wanted['tasks'][0]['input_template'] += 'file,2,changed.wfu\n'
    assert vrc_graph.select_graph(tmp_path, 'asymptote', wanted) != old


def test_original_complete_preferred_and_choice_does_not_flip(tmp_path):
    old = graph(tmp_path/'asymptote_old', spec(), job='10')
    new = graph(tmp_path/'asymptote_new', spec(walltime='7-00:00:00'), 'submitted', '20')
    wanted = spec(walltime='14-00:00:00')
    assert vrc_graph.select_graph(tmp_path, 'asymptote', wanted) == old
    state = dispatch._load(new)[2]
    state['tasks']['point']['status'] = 'complete'
    atomic_json(new/'state.json', state)
    assert vrc_graph.select_graph(tmp_path, 'asymptote', wanted) == old


def test_changed_manifest_is_not_silently_accepted(tmp_path):
    old = graph(tmp_path/'asymptote_old', spec())
    vrc_graph.select_graph(tmp_path, 'asymptote', spec())
    saved = json.loads((old/'workflow.json').read_text())
    saved['tasks'][0]['resources']['walltime'] = '7-00:00:00'
    atomic_json(old/'workflow.json', saved)
    with pytest.raises(RuntimeError, match='Workflow changed'):
        vrc_graph.select_graph(tmp_path, 'asymptote', spec())


def test_mr_json_identity_preserves_root_weights_and_orbital_provenance():
    from kinbot.molpro import validate_vrc_state_model
    model = validate_vrc_state_model({
        'state_average': [{'multiplicity': 2, 'roots': 2, 'weights': [.5, .5]}],
        'selected_roots': {'2': 1}, 'energy_selection': {'2': 'lowest'}})
    requested = {'tasks': [{'profile': {'calculator_kwargs': {'state_model': model}}}]}
    saved = json.loads(json.dumps(requested))
    assert model['selected_roots'] == {2: 1}
    assert saved['tasks'][0]['profile']['calculator_kwargs']['state_model']['selected_roots'] == {'2': 1}
    assert vrc_graph.scientific_spec(requested) == vrc_graph.scientific_spec(saved)
    assert vrc_graph.identity(requested) == vrc_graph.identity(saved)
    # Existing binding hashes of persisted specifications remain valid.
    import hashlib
    expected = hashlib.sha256(json.dumps(saved, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    assert vrc_graph.identity(requested) == expected
    for field, value in [('selected_roots', {'2': 2}),
                         ('energy_selection', {}),
                         ('state_average', [{'multiplicity': 2, 'roots': 2, 'weights': [.4, .6]}]),
                         ('orbital_selection', {'mode': 'molpro_record', 'path': '/different.wfu', 'sha256': 'b'*64})]:
        changed = deepcopy(saved)
        changed['tasks'][0]['profile']['calculator_kwargs']['state_model'][field] = value
        assert vrc_graph.identity(changed) != vrc_graph.identity(requested)
    assert model['selected_roots'] == {2: 1}  # Caller remains unmodified.


def test_ambiguous_or_nonfinite_identity_is_rejected():
    for model in ({'selected_roots': {2: 1, '2': 2}}, {'weight': float('nan')}):
        with pytest.raises(ValueError):
            vrc_graph.identity({'tasks': [{'model': model}]})


@pytest.fixture
def duplicates(tmp_path):
    root = tmp_path/'run'
    worker = root/'case'/'123'
    scans = worker/'vrctst/molpro_scan'
    old = graph(scans/'asymptote_old', spec(), job='10')
    new = graph(scans/'asymptote_new', spec(walltime='7-00:00:00'), 'submitted', '20')
    old_point = spec('point09'); old_point['tasks'][0]['input_template'] += f'file,2,{old}/seed.wfu\n'
    point = graph(scans/'point09_old', old_point, 'submitted', '21')
    new_point = spec('point00'); new_point['tasks'][0]['input_template'] += f'file,2,{new}/seed.wfu\n'
    derived = graph(scans/'point00_new', new_point, 'submitted', '22')
    pivot = graph(worker/'vrctst/molpro_pivots/fragment_0_new', spec('pivot'), job='23')
    evidence = worker/'vrctst/pivot_cubes/leaf/fragment_0.orbitals.json'
    atomic_json(evidence, {'run_dir': str(pivot), 'geometry_source': {'dispatch': str(new)}})
    evidence.with_suffix('').with_suffix('.cube').write_text('cube from repeated asymptote')
    jobs = [{'job_id': job, 'directory': str(p/'tasks/point'), 'state': 'RUNNING'}
            for p, job in [(new, '20'), (point, '21'), (derived, '22')]]
    return root, old, new, point, derived, pivot, evidence, jobs


def test_plan_traces_changed_geometry_pivots_and_wfu_chain_only(duplicates):
    root, old, new, point, derived, pivot, evidence, jobs = duplicates
    plan = vrc_deduplicate.plan(root, jobs)
    assert set(plan['keep']) == {str(old)}
    assert set(plan['remove']) == {str(new), str(derived), str(pivot)}
    assert {j['job_id'] for j in plan['selected_jobs']} == {'20', '22'}
    assert str(evidence) in plan['artifacts']
    assert point.exists() and new.exists()  # Preview is read-only.


def test_original_corruption_blocks_retirement(duplicates):
    root, old, *rest = duplicates
    (old/'tasks/point/point.out').write_text('corrupted')
    with pytest.raises(RuntimeError, match='artifact.*changed'):
        vrc_deduplicate.plan(root, rest[-1])


def test_unfinished_duplicate_group_requires_review(duplicates):
    root, old, *rest = duplicates
    state = dispatch._load(old)[2]
    state['tasks']['point']['status'] = 'submitted'
    atomic_json(old/'state.json', state)
    with pytest.raises(RuntimeError, match='no completed original'):
        vrc_deduplicate.plan(root, rest[-1])


def test_repair_archive_is_resumable_and_retains_original_chain(duplicates, monkeypatch, tmp_path):
    root, old, new, point, derived, pivot, evidence, jobs = duplicates
    vrc_graph.bind(new.parent, dispatch._load(new)[1], new)
    stale = vrc_graph.binding(new.parent, dispatch._load(new)[1])
    plan = vrc_deduplicate.plan(root, jobs)
    assert str(stale) in plan['artifacts']
    archive = tmp_path/'archive'
    journal = {'vrc_deduplication': plan}; marker = root/'refresh.json'
    calls = []
    def cancel(record, journal, marker):
        calls.append([j['job_id'] for j in record['selected_jobs']])
    monkeypatch.setattr('kinbot.vrc_restart.cancel_scans', cancel)
    monkeypatch.setattr('kinbot.run_refresh._owned_jobs', lambda root: [jobs[1]])
    before = {p: file_sha256(p) for g in (old, point) for p in g.rglob('*') if p.is_file()}
    rename = Path.rename
    interrupted = False
    def stop(path, destination):
        nonlocal interrupted
        if path == derived and not interrupted:
            interrupted = True
            raise OSError('interrupted archive')
        return rename(path, destination)
    monkeypatch.setattr(Path, 'rename', stop)
    with pytest.raises(OSError, match='interrupted archive'):
        vrc_deduplicate.apply(root, plan, archive, journal, marker)
    vrc_deduplicate.apply(root, plan, archive, journal, marker)
    vrc_deduplicate.apply(root, plan, archive, journal, marker)
    assert all((archive/p.relative_to(root)/'workflow.json').exists() for p in (new, derived, pivot))
    assert not new.exists() and not derived.exists() and not pivot.exists()
    assert before == {p: file_sha256(p) for p in before}
    assert vrc_graph.select_graph(old.parent, 'asymptote', spec(walltime='7-00:00:00')) == old
    assert vrc_graph.select_graph(point.parent, 'point09', dispatch._load(point)[1]) == point


def test_unknown_or_changed_ownership_never_cancelled(duplicates):
    root, old, new, point, derived, pivot, evidence, jobs = duplicates
    jobs[0]['job_id'] = '999'
    with pytest.raises(RuntimeError, match='ownership changed'):
        vrc_deduplicate.plan(root, jobs)


def test_interrupted_binding_after_archive_resumes_from_committed_journal(duplicates, monkeypatch, tmp_path):
    root, old, new, point, derived, pivot, evidence, jobs = duplicates
    vrc_graph.bind(new.parent, dispatch._load(new)[1], new)
    plan = vrc_deduplicate.plan(root, jobs)
    journal = {'vrc_deduplication': plan}; marker = root/'refresh.json'
    monkeypatch.setattr('kinbot.vrc_restart.cancel_scans', lambda *a: None)
    monkeypatch.setattr('kinbot.run_refresh._owned_jobs', lambda root: [jobs[1]])
    original_bind = vrc_deduplicate.bind
    def interrupted(*a):
        original_bind(*a)
        raise OSError('interrupted binding')
    monkeypatch.setattr(vrc_deduplicate, 'bind', interrupted)
    with pytest.raises(OSError, match='interrupted binding'):
        vrc_deduplicate.apply(root, plan, tmp_path/'archive', journal, marker)
    restored = json.loads(marker.read_text())
    assert restored['vrc_deduplication']['phase'] == 'archived'
    monkeypatch.setattr(vrc_deduplicate, 'bind', original_bind)
    vrc_deduplicate.apply(root, restored['vrc_deduplication'], tmp_path/'archive', restored, marker)
    assert vrc_graph.select_graph(old.parent, 'asymptote', spec(walltime='7-00:00:00')) == old
    assert point.exists() and not new.exists()


def test_completed_corrections_require_review(duplicates):
    root, old, new, point, derived, pivot, evidence, jobs = duplicates
    (old.parent.parent/'corr_leaf.json').write_text('{}')
    with pytest.raises(RuntimeError, match='correction potentials already exist'):
        vrc_deduplicate.plan(root, jobs)


@pytest.mark.parametrize('prebound', [False, True])
@pytest.mark.parametrize('initial_status', ['complete', 'submitted'])
def test_actual_mr_scan_resumes_legacy_checkpoint_chain_after_walltime_change(
        tmp_path, monkeypatch, initial_status, prebound):
    from types import SimpleNamespace
    import numpy as np
    from kinbot.vrc_tst_scan import VTS
    from kinbot.molpro import automatic_vrc_state_model
    monkeypatch.chdir(tmp_path)
    geom = np.asarray([[0., 0., 0.], [0., 0., 2.5]])
    reaction = SimpleNamespace(scan_coo=[0, 1], equiv=[[0], [1]], products=[],
        maps=[[0], [1]], irc_prod=SimpleNamespace(bondlist=[]),
        species=SimpleNamespace(atom=['H', 'H'], charge=0, mult=1))
    par = dict(vrc_tst_scan_method='caspt2(2,2)', vrc_tst_scan_basis='vdz',
               vrc_tst_walltime='24:00:00', single_point_ppn=4,
               vrc_tst_min_stack_mw=50, queue_name='short')
    par['vrc_tst_state_model'] = automatic_vrc_state_model(
        ['H', 'H'], 0, [1, 3], roots={1: 2, 3: 1},
        selected_roots={1: 1, 3: 1},
        active_space={'active_electrons': 2, 'active_orbitals': 2})
    captured = {}
    class Prepared(Exception): pass
    prepare = dispatch.prepare
    def capture(value, path):
        captured.update(value)
        raise Prepared
    monkeypatch.setattr(dispatch, 'prepare', capture)
    with pytest.raises(Prepared):
        VTS(SimpleNamespace(), par, None)._molpro_multireference_scan_task(reaction, geom, 2.5, 'point')
    monkeypatch.setattr(dispatch, 'prepare', prepare)
    captured['tasks'][0]['resources'] = {'cores': 4, 'memory_mb': 16000, 'walltime': '24:00:00'}
    old = prepare(captured, tmp_path/'vrctst/molpro_scan/point_legacy')
    _mock_execution(old, 'relaxed_scan', geometry='final.xyz')
    work = old/'tasks/relaxed_scan'
    wfu = work/'accepted.wfu'; wfu.write_text('accepted native orbitals')
    execution = json.loads((work/'execution.json').read_text())
    model = {'orbital_selection': {'mode': 'molpro_record', 'path': str(wfu), 'sha256': file_sha256(wfu)}}
    execution['details'].update(orbital_checkpoint_model=model, energy_ev=-1.)
    execution['artifacts'][wfu.name] = file_sha256(wfu)
    atomic_json(work/'execution.json', execution)
    state = dispatch._load(old)[2]
    state['tasks']['relaxed_scan'].update(status=initial_status, job_id='42',
        final_geometry_sha256=execution['details']['final_geometry_sha256'])
    atomic_json(old/'state.json', state)
    if prebound:
        # The repair binds the JSON-loaded original; the restarted controller
        # constructs an in-memory MR model with integer multiplicity keys.
        vrc_graph.bind(old.parent, dispatch._load(old)[1], old)
    monkeypatch.setattr(dispatch, 'preflight_until_scheduler_ready', lambda *a, **kw: None)
    adopted = []
    def drive(path, *args, **kwargs):
        assert path == old and dispatch._load(path)[2]['tasks']['relaxed_scan']['job_id'] == '42'
        adopted.append('42')
        state['tasks']['relaxed_scan']['status'] = 'complete'
        atomic_json(old/'state.json', state)
        return 0, {}
    monkeypatch.setattr('kinbot.vrc_tst_scan._drive_vrc_dispatch', drive)
    par['vrc_tst_walltime'] = '7-00:00:00'
    result = VTS(SimpleNamespace(), par, None)._molpro_multireference_scan_task(reaction, geom, 2.5, 'point')
    assert result[2] == old and result[4] == model
    assert adopted == (['42'] if initial_status == 'submitted' else [])
    assert len(list(old.parent.glob('*/workflow.json'))) == 1


def test_actual_fixed_pivot_reuses_original_cube_after_resource_change(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from ase import Atoms
    from kinbot.parameters import Parameters
    from kinbot.vrc_orbitals import fixed_pivot_spec, prepare_fixed_pivot
    from tests.test_orbital_cube import write_cube
    monkeypatch.chdir(tmp_path)
    inp = tmp_path/'input.json'
    inp.write_text('{"barrier_threshold": 100, "theory_preset": "b3lyp-b2plyp-anl"}')
    par = Parameters(str(inp)).par
    par['vrc_tst_walltime'] = '24:00:00'
    atoms = Atoms('CO', positions=[[0, 0, 0], [1.13, 0, 0]])
    product = SimpleNamespace(atom=atoms.get_chemical_symbols(), geom=atoms.positions, mult=1, charge=0)
    original = fixed_pivot_spec(product, par)
    original['tasks'][0]['resources'] = {'cores': 8, 'memory_mb': 4096, 'walltime': '24:00:00'}
    old = dispatch.prepare(original, tmp_path/'vrctst/molpro_pivots/fragment_0_fixed_old')
    cube = write_cube(old/'tasks/pivot/homo.cube', atoms=[(6, atoms.positions[0]), (8, atoms.positions[1])])
    _mock_execution(old, 'pivot')
    state = dispatch._load(old)[2]
    state['tasks']['pivot'].update(status='complete', job_id='10')
    atomic_json(old/'state.json', state)
    monkeypatch.setattr(dispatch, 'preflight_until_scheduler_ready', lambda *a, **kw: None)
    selected = []
    def drive(path, *a, **kw):
        assert path == old
        vrc_deduplicate.verify_complete(path)
        selected.append(path)
        return 0, {}
    monkeypatch.setattr('kinbot.vrc_tst_scan._drive_vrc_dispatch', drive)
    par['vrc_tst_walltime'] = '7-00:00:00'
    destination = tmp_path/'vrctst/pivot_cubes/channel/fragment_0'
    result = prepare_fixed_pivot(product, par, destination, {'dispatch': 'original_asymptote'})
    assert selected == [old] and Path(result).read_bytes() == cube.read_bytes()
    assert json.loads(destination.with_suffix('.orbitals.json').read_text())['run_dir'] == str(old)
    assert len(list(old.parent.glob('*/workflow.json'))) == 1
