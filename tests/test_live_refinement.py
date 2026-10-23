"""Live refinement changes one scientific dependency chain, not its siblings."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from kinbot import live_control as control, run_control
from kinbot.anl import dispatch, live_control as owner, refinement
from kinbot.anl.validation import configure_task_convergence
from kinbot.execution import atomic_json, file_sha256
from kinbot.manual_fallback import logical_task
from tests.anl_fixture import dispatch_spec, molpro_task
from tests.test_anl_dispatch import _mock_execution

OPTIONS = dict(fmax=.0001, scf=1.e-11, cc=1.e-11)


@pytest.fixture
def graph(tmp_path, monkeypatch):
    control.configure(tmp_path)
    spec = dispatch_spec()
    geom = deepcopy(spec['tasks'][1]); geom['geometry_from'] = 'initial'
    geom['optimizer']['fmax'] = .0005
    harmonic = deepcopy(spec['tasks'][2])
    harmonic['input_template'] = harmonic['input_template'].replace('rhf\nuccsd(t)',
        'gthresh,energy=1.d-8,coeff=1.d-7,civec=1.d-7\nrhf\n{uccsd(t),thrden=1.d-7,thrvar=1.d-10}')
    harmonic['result_parser']['stationary_point'] = 'transition_state'
    configure_task_convergence(geom, scf_convergence=1.e-8, cc_convergence=1.e-7)
    qz = deepcopy(geom); qz.update(id='qz_geometry', geometry_from='l3_geometry')
    spec['tasks'] = [geom, harmonic, qz,
        molpro_task('downstream', 'rhf\n', geometry_from='qz_geometry'),
        molpro_task('independent', 'rhf\n', geometry_from='initial'),
        molpro_task('done', 'rhf\n', geometry_from='initial')]
    g = dispatch.prepare(spec, tmp_path/'graph')
    _mock_execution(g, 'l3_geometry', geometry='final.xyz')
    _mock_execution(g, 'done')
    dispatch.advance(g)
    state = dispatch._load(g)[2]
    state['tasks']['harmonic'].update(status='failed', job_id='42')
    state['tasks']['qz_geometry'].update(status='submitted', job_id='43')
    state['tasks']['independent'].update(status='submitted', job_id='99')
    e = state['tasks']['harmonic']
    atomic_json(g/'tasks/harmonic/execution.json', dict(schema=1, task_id='harmonic',
        status='failed', failure_phase='validation', error='wrong stationary-point order',
        geometry_sha256=e['geometry_sha256'], slurm_job_id='42'))
    atomic_json(g/'state.json', state)
    monkeypatch.setattr(dispatch, '_job_active', lambda j: j in ('43', '99'))
    return g


def queue(g, cancel=False):
    p = refinement.plan(g, 'harmonic', OPTIONS)
    p.pop('replacement')
    r = control.request(g.parent, g/'tasks/harmonic', kind='refinement', reason='tighten geometry and Hessian',
        cancel_running=cancel, options={'refinement_plan': p})
    return g.parent/'.kinbot_controls/requests'/(r['id']+'.json')


def test_preview_is_read_only_and_covers_all_dependent_results(graph, capsys):
    before = {p: file_sha256(p) for p in graph.rglob('*') if p.is_file()}
    run_control.main(['refine', str(graph.parent), '--target', 'graph/tasks/harmonic',
        '--fmax', '.0001', '--scf', '1e-11', '--cc', '1e-11', '--reason', 'review'])
    plan = json.loads(capsys.readouterr().out)
    assert plan['reset'] == ['downstream', 'harmonic', 'l3_geometry', 'qz_geometry']
    assert plan['changed'] == ['l3_geometry', 'harmonic', 'qz_geometry']
    assert before == {p: file_sha256(p) for p in before}
    assert not list((graph.parent/'.kinbot_controls/requests').glob('*.json'))


def test_default_waits_without_cancellation_and_freezes_new_descendants(graph, monkeypatch):
    p = queue(graph)
    monkeypatch.setattr(owner, 'stop_attempt', lambda *a, **k: pytest.fail('No cancellation authorized'))
    owner.process(graph)
    assert json.loads(p.read_text())['status'] == 'applying'
    assert dispatch._load(graph)[2]['tasks']['qz_geometry']['job_id'] == '43'
    # Even if the source finishes during the wait, no old-source descendant
    # may newly stage/submit while the refinement transaction is pending.
    _mock_execution(graph, 'qz_geometry', geometry='final.xyz')
    monkeypatch.setattr(dispatch, '_job_active', lambda j: j == '99')
    dispatch.advance(graph)
    assert not (graph/'tasks/downstream').exists()
    owner.process(graph)
    assert json.loads(p.read_text())['status'] == 'applied'
    state = dispatch.advance(graph)
    assert state['tasks']['l3_geometry']['status'] == 'staged'


def test_explicit_replacement_keeps_unrelated_jobs_and_verified_seed(graph, monkeypatch):
    original = dispatch._load(graph)[1]
    hashes = {p: file_sha256(p) for n in ('independent', 'done')
              for p in (graph/'tasks'/n).iterdir() if p.is_file()}
    seed = (graph/'tasks/l3_geometry/final.xyz').read_bytes()
    p = queue(graph, cancel=True)
    calls = []
    def stop(directory, job, path, record, **kw):
        calls.append((directory.name, job, kw))
        return True
    monkeypatch.setattr(owner, 'stop_attempt', stop)
    owner.process(graph)
    assert calls == [('qz_geometry', '43', {'allow_completed': True})]
    assert json.loads(p.read_text())['status'] == 'applied'
    state = dispatch.advance(graph)
    assert state['tasks']['independent']['job_id'] == '99'
    assert state['tasks']['done']['status'] == 'complete'
    assert hashes == {p: file_sha256(p) for p in hashes}
    assert (graph/'tasks/l3_geometry/continuation.xyz').read_bytes() == seed
    spec = dispatch._load(graph)[1]
    for old, new in zip(original['tasks'], spec['tasks']):
        assert logical_task(new) == old
    geom = spec['tasks'][0]
    assert geom['optimizer']['fmax'] == .0001
    assert geom['profile']['calculator_kwargs']['scf_convergence'] == 1.e-11
    assert 'thrden=1.d-11,thrvar=1.d-11' in spec['tasks'][1]['input_template']
    assert spec['tasks'][1]['result_parser']['stationary_point'] == 'transition_state'
    # Idempotent owner polling cannot perform the request a second time.
    owner.process(graph)
    assert len(calls) == 1 and len(list((graph/'attempts').glob('live_*'))) == 1


@pytest.mark.parametrize('change', ['generation', 'geometry', 'settings'])
def test_changed_generation_or_provenance_never_cancels(graph, monkeypatch, change):
    p = queue(graph, cancel=True)
    if change == 'generation':
        state = dispatch._load(graph)[2]
        state['tasks']['qz_geometry']['job_id'] = 'new'
        atomic_json(graph/'state.json', state)
    elif change == 'geometry':
        with (graph/'tasks/l3_geometry/final.xyz').open('a') as f: f.write('changed')
    else:
        spec = dispatch._load(graph)[1]
        spec['tasks'][0]['optimizer']['fmax'] = .00001
        atomic_json(graph/'workflow.json', spec)
        state = json.loads((graph/'state.json').read_text())
        state['spec_sha256'] = file_sha256(graph/'workflow.json')
        atomic_json(graph/'state.json', state)
    monkeypatch.setattr(owner, 'stop_attempt', lambda *a, **k: pytest.fail('Must validate first'))
    owner.process(graph)
    assert json.loads(p.read_text()).get('error')
    assert not list((graph/'attempts').glob('live_*'))


def test_interrupted_commit_resumes_once(graph, monkeypatch):
    from kinbot import run_recover
    p = queue(graph, cancel=True)
    monkeypatch.setattr(owner, 'stop_attempt', lambda *a, **k: True)
    save = run_recover.atomic_json
    def interrupt(path, data):
        if path == graph/'state.json': raise OSError('interrupted')
        save(path, data)
    monkeypatch.setattr(run_recover, 'atomic_json', interrupt)
    with pytest.raises(OSError): owner.process(graph)
    monkeypatch.setattr(run_recover, 'atomic_json', save)
    owner.process(graph)
    assert json.loads(p.read_text())['status'] == 'applied'
    assert len(list((graph/'attempts').glob('live_*'))) == 1


def test_stricter_existing_thresholds_are_preserved(graph):
    spec = dispatch._load(graph)[1]
    task = deepcopy(spec['tasks'][1])
    task['input_template'] = task['input_template'].replace('thrvar=1.d-10', 'thrvar=1.d-13')
    assert 'thrvar=1.d-13' in refinement._frequency_task(task, OPTIONS)['input_template']
    with pytest.raises(ValueError, match='loosen'):
        refinement.plan(graph, 'harmonic', dict(OPTIONS, fmax=.1))


def test_request_is_invisible_to_old_single_task_owners(graph):
    p = queue(graph)
    assert not control.requests(graph/'tasks', kinds={'task'})
    assert [path for path, _ in control.requests(graph/'tasks', kinds={'refinement'})] == [p]


@pytest.mark.parametrize('order', [0, 1])
def test_qz_refinement_archives_completed_consumers_without_resetting_tz(graph, monkeypatch, order):
    spec = dispatch._load(graph)[1]
    spec['tasks'][2]['optimizer']['order'] = order
    qh = deepcopy(spec['tasks'][1])
    qh.update(id='qz_harmonic', geometry_from='qz_geometry')
    qh['result_parser']['stationary_point'] = 'transition_state' if order else 'minimum'
    spec['tasks'].append(qh)
    # Build a distinct graph whose QZ geometry and energy have finished, but
    # whose QZ harmonic result fails. Both minima and saddles use the same DAG.
    g = dispatch.prepare(spec, graph.parent/('qz'+str(order)))
    _mock_execution(g, 'l3_geometry', geometry='final.xyz')
    dispatch.advance(g)
    _mock_execution(g, 'qz_geometry', geometry='final.xyz')
    dispatch.advance(g)
    _mock_execution(g, 'downstream')
    state = dispatch.advance(g)
    state['tasks']['qz_harmonic']['status'] = 'failed'
    atomic_json(g/'state.json', state)
    old_tz = {p: file_sha256(p) for p in (g/'tasks/l3_geometry').iterdir() if p.is_file()}
    p = refinement.plan(g, 'qz_harmonic', OPTIONS)
    assert p['reset'] == ['downstream', 'qz_geometry', 'qz_harmonic']
    assert refinement.plan(g, 'downstream', OPTIONS)['geometry'] == 'qz_geometry'
    p.pop('replacement')
    r = control.request(g.parent, g/'tasks/qz_harmonic', kind='refinement',
                        reason='QZ curvature validation', options={'refinement_plan': p})
    monkeypatch.setattr(dispatch, '_job_active', lambda job: False)
    owner.process(g)
    record = json.loads((g.parent/'.kinbot_controls/requests'/(r['id']+'.json')).read_text())
    assert record['status'] == 'applied'
    archived = Path(record['archive'])/'graphs/tasks/downstream'
    assert (archived/'execution.json').is_file()
    assert old_tz == {p: file_sha256(p) for p in old_tz}
    state = dispatch.advance(g)
    assert state['tasks']['qz_geometry']['status'] == 'staged'
    assert not (g/'tasks/downstream').exists()
    # The energy becomes eligible only after the replacement QZ geometry.
    _mock_execution(g, 'qz_geometry', geometry='final.xyz')
    state = dispatch.advance(g)
    assert state['tasks']['downstream']['status'] == 'staged'
    assert state['tasks']['downstream']['geometry_sha256'] == state['tasks']['qz_geometry']['final_geometry_sha256']
    assert not (g/'tasks/downstream/execution.json').exists()


def test_final_handoff_is_fenced_by_pending_refinement(graph):
    refinement.window(graph)
    p = queue(graph)
    assert not refinement.window(graph, close=True)
    request = json.loads(p.read_text())
    control.acknowledge(p, request, 'blocked', error='test review cancelled')
    assert refinement.window(graph, close=True)
    with pytest.raises(ValueError, match='handed downstream'):
        refinement.plan(graph, 'harmonic', OPTIONS)


def test_multiple_affected_allocations_get_separate_cancel_receipts(graph, monkeypatch):
    state = dispatch._load(graph)[2]
    state['tasks']['harmonic']['status'] = 'submitted'
    # The geometry failure supplies the diagnostic while two dependent jobs
    # from a previous generation are still live.
    state['tasks']['l3_geometry']['status'] = 'failed'
    atomic_json(graph/'state.json', state)
    monkeypatch.setattr(dispatch, '_job_active', lambda j: j in ('42', '43', '99'))
    p = queue(graph, cancel=True)
    seen = []
    def cancel(directory, job, path, record, **kw):
        assert 'cancel_requested' not in record
        record['cancel_requested'] = True
        seen.append((job, str(path)))
        return True
    monkeypatch.setattr(owner, 'stop_attempt', cancel)
    owner.process(graph)
    assert [j for j, _ in seen] == ['42', '43']
    assert len({path for _, path in seen}) == 2
    assert json.loads(p.read_text())['status'] == 'applied'
