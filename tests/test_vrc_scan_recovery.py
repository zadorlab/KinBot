"""No native QC: exercise interrupted scans, bounded retry and review gates."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.db import connect
from ase.io import write

from kinbot import vrc_scan_recovery as recovery
from kinbot.execution import atomic_json, file_sha256
from kinbot.vrc_tst_scan import VTS


@pytest.fixture
def scan(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    geometry = np.array([[0., 0., 0.], [0., 0., .74], [0., 0., 3.24], [0., 0., 3.98]])
    products = [SimpleNamespace(atom=['H', 'H'], natom=2, mult=1,
        geom=geometry[part].copy(), bond01=np.array([[0, 1], [1, 0]]))
        for part in (slice(0, 2), slice(2, 4))]
    reaction = SimpleNamespace(scan_coo=[1, 2], maps=[[0, 1], [2, 3]], equiv=[[1], [2]],
        products=products, do_vdW=False,
        species=SimpleNamespace(atom=['H']*4, geom=geometry.copy(), charge=0, mult=1),
        irc_prod=SimpleNamespace(geom=geometry.copy(), bondlist=[[0, 1], [2, 3]]))
    parameters = dict(vrc_tst_fragment_geometry='asymptotic_scan', pp_orient='geometric',
        vrc_tst_scan_method='caspt2(2,2)', vrc_tst_scan_basis='avdz-f12',
        vrc_tst_scan_points=[2.5, 2.7, 2.9], vrc_tst_asymptote_distance=30.)
    fixture = SimpleNamespace(calls=[], retry_fails=False, retry_native_fails=False,
                              repair_primary=False, interrupt=None,
                              extra_failures=set(),
                              results={}, original=deepcopy(reaction))
    def optimize(reaction, geom, distance, logical_name, **kwargs):
        fixture.calls.append((logical_name, np.array(geom), deepcopy(kwargs)))
        if fixture.interrupt and logical_name.endswith('pt02'):
            raise fixture.interrupt
        graph = tmp_path/'vrctst/molpro_scan'/logical_name
        task = graph/'tasks/relaxed_scan'; task.mkdir(parents=True, exist_ok=True)
        atomic_json(task/'task.json', {'task': {'id': 'relaxed_scan'}})
        if logical_name.endswith('pt01') and not fixture.repair_primary:
            raise recovery.ScanPointFailure('H transfer', graph, connectivity=True)
        if logical_name in fixture.extra_failures:
            raise recovery.ScanPointFailure('H transfer', graph, connectivity=True)
        if logical_name.endswith('fragment_retry') and fixture.retry_fails:
            raise recovery.ScanPointFailure('still H transfer', graph, connectivity=True)
        if logical_name.endswith('fragment_retry') and fixture.retry_native_fails:
            raise recovery.ScanPointFailure('Optimization budget exhausted', graph)
        a = Atoms('H4', positions=geom)
        a.positions[0] = a.positions[1] - [0, 0, .8 if distance == 30. else .81]
        a.positions[3] = a.positions[2] + [0, 0, .82 if distance == 30. else .83]
        model = {'seed': logical_name}
        write(task/'final.xyz', a)
        atomic_json(task/'execution.json', {'status': 'executed'})
        atomic_json(graph/'workflow.json', {'name': logical_name})
        result = a, -2., graph, {'details': {'final_geometry_sha256': file_sha256(task/'final.xyz')}}, model
        fixture.results[logical_name] = result
        return result
    def make():
        vts = VTS(SimpleNamespace(), deepcopy(parameters), None)
        vts.scan_reac['r'] = deepcopy(fixture.original)
        monkeypatch.setattr(vts, '_molpro_multireference_scan_task', optimize)
        return vts
    fixture.make = make
    fixture.review = tmp_path/'vrctst/scan_review/r/scan_review.json'
    return fixture


def test_retry_uses_asymptotic_coordinates_and_orbitals_then_continues(scan):
    vts = scan.make(); vts.do_scan(['r'])
    assert [c[0] for c in scan.calls] == ['r_vts_fragments_asymptote', 'r_vts_pt00',
        'r_vts_pt01', 'r_vts_pt01_fragment_retry', 'r_vts_pt02']
    retry = scan.calls[3]
    assert np.linalg.norm(retry[1][1]-retry[1][0]) == pytest.approx(.8)
    assert retry[2] == dict(fragment_retry=True, state_model_override={'seed': 'r_vts_fragments_asymptote'})
    assert scan.calls[4][2]['state_model_override'] == {'seed': 'r_vts_pt01_fragment_retry'}
    assert vts._scan_review_requests  # Whole path held even with a complete grid.
    record = json.loads(scan.review.read_text())
    assert record['complete_grid'] and record['status'] == 'needs_review'
    assert record['points']['r_vts_pt01']['attempts'][0]['connectivity']
    recovery.approve(scan.review, 'Reviewed every geometry and orbital record.')
    assert recovery.review_changed(vts._scan_review_requests)
    next_vts = scan.make(); next_vts.do_scan(['r'])
    assert not next_vts._scan_review_requests
    assert json.loads(scan.review.read_text())['status'] == 'resolved'


@pytest.mark.parametrize('native_failure', [False, True])
def test_failed_retry_skips_gap_and_next_point_restarts_from_fragments(scan, native_failure):
    scan.retry_fails = not native_failure
    scan.retry_native_fails = native_failure
    with pytest.raises(recovery.ScanReviewRequired):
        scan.make().do_scan(['r'])
    next_point = scan.calls[-1]
    assert next_point[0] == 'r_vts_pt02'
    assert next_point[2]['state_model_override'] == {'seed': 'r_vts_fragments_asymptote'}
    assert np.linalg.norm(next_point[1][1]-next_point[1][0]) == pytest.approx(.8)
    db = connect('kinbot.db')
    assert not list(db.select(name='vrctst/r_vts_pt01'))
    assert list(db.select(name='vrctst/r_vts_pt02'))
    assert list(db.select(name='vrctst/r_vts_pt_asymptote'))
    record = json.loads(scan.review.read_text())
    assert record['points']['r_vts_pt01']['status'] == 'needs_review'
    assert record['points']['r_vts_pt01']['attempts'][-1]['connectivity'] is not native_failure
    with pytest.raises(ValueError, match='complete accepted grid'):
        recovery.approve(scan.review, 'Cannot waive a gap.')


def test_later_successful_seed_is_retained_when_earlier_gap_is_repaired(scan):
    scan.retry_fails = True
    with pytest.raises(recovery.ScanReviewRequired): scan.make().do_scan(['r'])
    before = scan.calls[-1]
    scan.repair_primary = True
    vts = scan.make(); vts.do_scan(['r'])
    after = scan.calls[-1]
    assert after[1] == pytest.approx(before[1])
    assert after[2] == before[2]
    assert vts._scan_review_requests  # Prior connectivity event cannot disappear.


def test_after_gap_fragment_primary_does_not_repeat_identical_retry(scan):
    scan.retry_fails = True
    scan.extra_failures.add('r_vts_pt02')
    with pytest.raises(recovery.ScanReviewRequired):
        scan.make().do_scan(['r'])
    names = [c[0] for c in scan.calls]
    assert names[-1] == 'r_vts_pt02'
    assert 'r_vts_pt01_fragment_retry' in names  # This was a different seed.
    assert 'r_vts_pt02_fragment_retry' not in names
    record = json.loads(scan.review.read_text())
    point = record['points']['r_vts_pt02']
    assert point['status'] == 'needs_review'
    assert point['fragment_retry_skipped'] == 'identical_geometry_and_state_model'
    assert len(point['attempts']) == 1
    assert not record['complete_grid']
    with pytest.raises(ValueError, match='complete accepted grid'):
        recovery.approve(scan.review, 'A skipped duplicate cannot waive the gap.')
    scan.calls.clear()
    with pytest.raises(recovery.ScanReviewRequired):
        scan.make().do_scan(['r'])
    assert 'r_vts_pt02_fragment_retry' not in [c[0] for c in scan.calls]


def test_existing_same_seed_retry_is_reattached_even_without_progress_receipt(scan):
    scan.retry_fails = True
    scan.extra_failures.add('r_vts_pt02')
    directory = Path('vrctst/molpro_scan/r_vts_pt02_fragment_retry_legacy')
    directory.mkdir(parents=True)
    state = directory/'state.json'
    state.write_text('{"tasks":{"relaxed_scan":{"status":"submitted","job_id":"42"}}}')
    original = state.read_bytes()
    with pytest.raises(recovery.ScanReviewRequired):
        scan.make().do_scan(['r'])
    assert 'r_vts_pt02_fragment_retry' in [c[0] for c in scan.calls]
    assert state.read_bytes() == original
    record = json.loads(scan.review.read_text())
    assert len(record['points']['r_vts_pt02']['attempts']) == 2


def test_retry_seed_comparison_preserves_changed_geometry_or_orbitals(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    a = np.array([[0., 0., 0.], [0., 0., 1.]])
    model = {'selected_roots': {'2': 1}, 'orbital_selection': {'sha256': 'original'}}
    assert recovery.redundant_fragment_retry('r', a, model, a.copy(), deepcopy(model))
    b = a.copy(); b[1, 2] += 1.e-9
    assert not recovery.redundant_fragment_retry('r', a, model, b, model)
    changed = deepcopy(model); changed['orbital_selection']['sha256'] = 'different'
    assert not recovery.redundant_fragment_retry('r', a, model, a, changed)
    changed = deepcopy(model); changed['selected_roots']['2'] = 2
    assert not recovery.redundant_fragment_retry('r', a, model, a, changed)
    b[0, 0] = np.nan
    assert not recovery.redundant_fragment_retry('r', b, model, b, model)


def test_explicit_primary_fallback_is_not_deduplicated_from_logical_seed(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    graph = tmp_path/'primary'
    a = [[0., 0., 0.]]
    model = {'seed': 'fragment'}
    # Missing native task provenance is also insufficient for suppression.
    assert not recovery.redundant_fragment_retry('r', a, model, a, model, primary_graph=graph)
    task = graph/'tasks/relaxed_scan'; task.mkdir(parents=True)
    atomic_json(task/'task.json', {'task': {'id': 'relaxed_scan', 'manual_fallback': {
        'options': {'calculator_kwargs': {'correlation_convergence': 1.e-9}}}}})
    assert not recovery.redundant_fragment_retry('r', a, model, a, model, primary_graph=graph)


def test_refresh_interrupt_is_not_a_point_failure_and_seed_plan_survives(scan):
    from kinbot.pes_phase import ExplorationReopened
    scan.interrupt = ExplorationReopened('refresh')
    with pytest.raises(ExplorationReopened): scan.make().do_scan(['r'])
    before = scan.calls[-1]
    scan.interrupt = None
    vts = scan.make(); vts.do_scan(['r'])
    assert scan.calls[-1][1] == pytest.approx(before[1])
    assert scan.calls[-1][2] == before[2]
    assert vts._scan_review_requests


def test_changed_or_incomplete_evidence_cannot_be_approved(scan):
    scan.make().do_scan(['r'])
    graph = scan.results['r_vts_pt02'][2]
    (graph/'tasks/relaxed_scan/final.xyz').write_text('changed')
    with pytest.raises(ValueError, match='artifacts changed'):
        recovery.approve(scan.review, 'Should not pass')


def test_seed_tampering_fails_without_dispatch(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    progress = recovery.ScanProgress('r', {'method': 'caspt2'})
    progress.seed('pt', [[0, 0, 0]], None)
    path = progress.directory/'pt.json'
    record = json.loads(path.read_text()); record['seed']['geometry'][0][0] = 1.
    atomic_json(path, record)
    with pytest.raises(RuntimeError, match='seed changed'):
        progress.seed('pt', [[0, 0, 0]], None)


def test_old_correction_cannot_bypass_review_at_rotd_generation_or_result_read(scan):
    from kinbot import pes, rotdpy
    scan.make().do_scan(['r'])
    correction = Path('vrctst/corr_r.json')
    atomic_json(correction, {})
    par = dict(rotdpy_run=0, vrc_tst_scan={'parent': ['r']}, vrc_tst_noscan={})
    with pytest.raises(recovery.ScanReviewRequired):
        pes.create_rotdpy_inputs(par, [['parent', 'r', ['a', 'b'], 0.]], [], correction_root='.')
    atomic_json(Path('rotdPy/r.rotdpy.json'), dict(status='complete',
        validation={'correction': {'source': str(correction.resolve())}}))
    with pytest.raises(recovery.ScanReviewRequired):
        rotdpy.read_result('rotdPy/r.py')
    recovery.approve(scan.review, 'All geometries and orbital evidence reviewed.')
    scan.make().do_scan(['r'])  # Only the owner resolves the approved evidence.
    with pytest.raises(ValueError, match='not bound to the current'):
        recovery.require_ready_for_rates(correction)
    old_bytes = correction.read_bytes()
    recovery.bind_correction_review(correction)
    assert correction.read_bytes() == old_bytes
    recovery.require_ready_for_rates(correction)
    with pytest.raises(ValueError, match='sampled correction differs'):
        recovery.require_ready_for_rates(correction, expected_sha256='old')


def test_scan_hold_not_hidden_by_an_accepted_orbital_model(scan):
    from kinbot import run_monitor
    scan.make().do_scan(['r'])
    path = Path('vrc_model_review.json')
    data = dict(status='needs_review', channels=[dict(kind='scan_connectivity', reaction='r')])
    atomic_json(path, data)
    atomic_json(Path('vrctst/orbital_review/r/request.json'), dict(reaction='r', status='model_accepted'))
    assert run_monitor.vrc_review_status(path, data, []) == {}


def test_cli_retry_queues_only_gap_and_wakes_held_owner(scan, monkeypatch):
    from kinbot import live_control
    scan.retry_fails = True
    with pytest.raises(recovery.ScanReviewRequired): scan.make().do_scan(['r'])
    record = json.loads(scan.review.read_text())
    graph = Path(record['points']['r_vts_pt01']['attempts'][-1]['graph'])
    graph.mkdir(parents=True, exist_ok=True)
    live_control.configure(Path.cwd())
    monkeypatch.setattr(live_control, 'generation', lambda *args: {'attempt': 1})
    fallback = dict(schema=1, optimizer={'steps': 200})
    request = recovery.retry(scan.review, 2.7, 'Retry with reviewed optimizer budget', fallback)
    assert request['target'] == str(graph/'tasks/relaxed_scan')
    assert request['fallback'] == fallback
    assert recovery.review_changed([dict(kind='scan_connectivity', request=str(scan.review))])
    with pytest.raises(ValueError, match='accepted points are retained'):
        recovery.retry(scan.review, 2.5, 'Must not replace an accepted point')
    with pytest.raises(ValueError, match='already pending'):
        recovery.retry(scan.review, 2.7, 'Must not queue duplicate retries')


def test_live_vrc_owner_reenters_after_review_without_starting_a_driver(tmp_path, monkeypatch):
    from kinbot import kb, live_control, pes_phase
    monkeypatch.chdir(tmp_path)
    live_control.configure(tmp_path)
    monkeypatch.setattr(pes_phase, 'active_worker', lambda path: True)
    monkeypatch.setattr(pes_phase, 'checkpoint', lambda: None)
    monkeypatch.setattr(recovery, 'review_changed', lambda channels: True)
    calls = []
    class Surface:
        def __init__(self, *args):
            self.review_requests = [dict(reaction='r', kind='scan_connectivity')]
        def calculate_correction_potentials(self):
            calls.append(1)
            if len(calls) == 1:
                raise recovery.ScanReviewRequired(self.review_requests)
    monkeypatch.setattr(kb, 'VTS', Surface)
    qc = SimpleNamespace(queue_job_limit=3, queue_job_root='root', queue_job_root_recursive=True)
    assert kb._execute_vrc_branch(SimpleNamespace(reac_obj=[]), {}, qc) == {'status': 'complete'}
    assert len(calls) == 2
    assert json.loads(Path('vrc_model_review.json').read_text())['status'] == 'resolved'
    assert qc.queue_job_limit == 3 and qc.queue_job_root == 'root' and qc.queue_job_root_recursive


def test_bounded_graph_does_not_get_another_automatic_allocation(tmp_path, monkeypatch):
    from kinbot.anl import dispatch
    from tests.test_optimization_budget import failed_graph
    graph, name, _ = failed_graph(tmp_path, monkeypatch)
    _, spec, state = dispatch._load(graph)
    spec['intent'] = {'vrc_fragment_retry': True}
    atomic_json(graph/'workflow.json', spec)
    state['spec_sha256'] = file_sha256(graph/'workflow.json'); atomic_json(graph/'state.json', state)
    before = (graph/'tasks'/name/'execution.json').read_bytes()
    result = dispatch.advance(graph)
    assert result['tasks'][name]['status'] == state['tasks'][name]['status']
    assert not result.get('optimizer_continuations')
    assert (graph/'tasks'/name/'execution.json').read_bytes() == before


def test_explicit_live_retry_can_archive_scientifically_rejected_native_completion(tmp_path, monkeypatch):
    from kinbot import live_control
    from kinbot.anl import dispatch, live_control as owner
    from tests.anl_fixture import dispatch_spec
    from tests.test_anl_dispatch import _mock_execution
    live_control.configure(tmp_path)
    spec = dispatch_spec(); spec['tasks'] = [spec['tasks'][0]]
    spec['tasks'][0]['id'] = 'relaxed_scan'
    spec['intent'] = {'vrc_scan_distance_angstrom': 2.7}
    graph = dispatch.prepare(spec, tmp_path/'vrctst/molpro_scan/g')
    _mock_execution(graph, 'relaxed_scan', geometry='final.xyz')
    state = dispatch.advance(graph)
    state['tasks']['relaxed_scan']['job_id'] = '42'
    atomic_json(graph/'state.json', state)
    recovery.reject_point(graph, 'fragment_connectivity')
    assert recovery.rejected_point(graph, 'relaxed_scan')
    target = graph/'tasks/relaxed_scan'
    old_execution = (target/'execution.json').read_bytes()
    queued = live_control.request(tmp_path, target, reason='Review found wrong fragment minimum',
        fallback={'schema': 1, 'optimizer': {'steps': 120}})
    monkeypatch.setattr(owner, 'slurm_job_record', lambda *a:
        dict(JobId='42', JobState='COMPLETED', WorkDir=str(target)))
    owner.process(graph)
    state = dispatch.advance(graph)
    assert state['tasks']['relaxed_scan']['status'] == 'staged'
    assert state['tasks']['relaxed_scan']['attempt'] == 2
    old = Path(state['tasks']['relaxed_scan']['previous_attempt'])
    assert (old/'execution.json').read_bytes() == old_execution
    assert (old/'scan_rejection.json').is_file()
    assert json.loads((tmp_path/'.kinbot_controls/requests'/f'{queued["id"]}.json').read_text())['status'] == 'applied'
    owner.process(graph)
    assert dispatch._load(graph)[2]['tasks']['relaxed_scan']['attempt'] == 2


def test_changed_rejected_result_cannot_be_used_to_retry_complete_work(tmp_path):
    graph = tmp_path/'g'; task = graph/'tasks/relaxed_scan'; task.mkdir(parents=True)
    for p in (graph/'workflow.json', task/'execution.json', task/'final.xyz'):
        p.write_text('original')
    recovery.reject_point(graph, 'fragment_connectivity')
    (task/'final.xyz').write_text('different')
    assert not recovery.rejected_point(graph, 'relaxed_scan')


def test_scan_failure_does_not_stop_other_leaves_or_publish_bad_channel(tmp_path, monkeypatch):
    from kinbot import vrc_multistate
    from tests.test_vrc_parallel import surfaces
    monkeypatch.chdir(tmp_path)
    vts, choices = surfaces(tmp_path)
    def scan(child, names):
        if names[0].startswith('cs'):
            raise recovery.ScanReviewRequired([dict(kind='scan_connectivity', reaction=names[0],
                request='review.json', error='gap')])
    def energies(child, names):
        Path(f'vrctst/corr_{names[0]}.json').write_text('{}')
    monkeypatch.setattr(VTS, 'do_scan', scan)
    monkeypatch.setattr(VTS, 'energies', energies)
    with pytest.raises(recovery.ScanReviewRequired) as raised:
        vrc_multistate.calculate_channels(vts, choices, limit=4)
    assert {c['channel'] for c in raised.value.channels} == {'cs'}
    assert Path('vrctst/corr_ss.json').is_file()
    assert not Path('vrctst/corr_cs.json').exists()
    states = json.loads(Path('vrctst/parallel_surfaces.json').read_text())['surfaces']
    assert {v['status'] for v in states.values()} == {'needs_review', 'complete'}


def test_complete_but_flagged_scan_computes_corrections_and_holds_rates(tmp_path, monkeypatch):
    from kinbot import vrc_multistate
    from tests.test_vrc_parallel import surfaces
    monkeypatch.chdir(tmp_path)
    vts, choices = surfaces(tmp_path)
    computed = []
    def scan(child, names):
        child._scan_review_requests = [dict(kind='scan_connectivity', reaction=names[0],
                request='review.json', error='review whole path')]
    def energies(child, names):
        computed.append(names[0]); Path(f'vrctst/corr_{names[0]}.json').write_text('{}')
    monkeypatch.setattr(VTS, 'do_scan', scan)
    monkeypatch.setattr(VTS, 'energies', energies)
    with pytest.raises(recovery.ScanReviewRequired):
        vrc_multistate.calculate_channels(vts, choices, limit=4)
    assert len(computed) == 4
    assert not Path('vrctst/corr_cs.json').exists()
    assert not Path('vrctst/corr_ss.json').exists()
