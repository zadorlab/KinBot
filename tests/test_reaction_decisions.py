from copy import deepcopy
import json

import pytest

from kinbot import reaction_decisions as decision, live_control as control
from kinbot.anl import dispatch, live_control
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task
from tests.test_pes_selection import evidence


def setup(tmp_path, *, shared=False):
    worker = tmp_path/'case'/'well'
    worker.mkdir(parents=True)
    atomic_json(worker/'accepted_l2_network.json', evidence('well', [('well_hom_sci_1_2',['p','q'],0.,None)]))
    spec = dispatch_spec()
    spec['tasks'] = [molpro_task('well_hom_sci_1_2_vts_pt00', 'rhf\n', geometry_from='initial'),
                     molpro_task('unrelated', 'rhf\n', geometry_from='initial')]
    folder = worker/'vrctst'/'molpro'
    graph = dispatch.prepare(spec, folder/('dispatch' if shared else 'well_hom_sci_1_2_vts_pt00_hash'))
    control.configure(tmp_path)
    request = dict(id='test', kind='network', target=str(worker), status='queued', created_at='now',
        reason='Reviewed fragment transfer', cancel_running=True, pending_only=False)
    path = tmp_path/'.kinbot_controls/requests/test.json'
    atomic_json(path, request)
    return worker, graph, path, request


def test_durable_exclusion_preserves_raw_chemistry_and_native_files(tmp_path, monkeypatch):
    worker, graph, path, request = setup(tmp_path)
    monkeypatch.setattr(dispatch, '_job_active', lambda _:False)
    original = (worker/'accepted_l2_network.json').read_bytes()
    native = {p:file_sha256(p) for p in (graph/'tasks').glob('*/*') if p.is_file()}
    old = dispatch._load(graph)[2]
    assert decision.apply(worker, {'well_hom_sci_1_2':'reject'}, path, request)
    assert decision.rejected(worker, 'well_hom_sci_1_2')
    assert decision.graph_rejected(graph)
    assert {e['status'] for e in dispatch.advance(graph)['tasks'].values()} == {'skipped'}
    assert original == (worker/'accepted_l2_network.json').read_bytes()
    assert native == {p:file_sha256(p) for p in native}
    # Re-entry and future controller polling cannot reopen it.
    assert decision.apply(worker, {'well_hom_sci_1_2':'reject'}, path, request)
    assert decision.graph_rejected(graph)
    include = dict(request, id='include', reason='Reconsidered', status='queued')
    assert decision.apply(worker, {'well_hom_sci_1_2':'include'}, path.with_name('include.json'), include)
    assert not decision.graph_rejected(graph)
    assert dispatch._load(graph)[2] == old
    assert not decision.rejected(worker, 'well_hom_sci_1_2')


def test_shared_correction_batch_only_excludes_selected_channel(tmp_path, monkeypatch):
    worker, graph, path, request = setup(tmp_path, shared=True)
    monkeypatch.setattr(dispatch, '_job_active', lambda _:False)
    before = deepcopy(dispatch._load(graph)[2]['tasks']['unrelated'])
    assert decision.apply(worker, {'well_hom_sci_1_2':'reject'}, path, request)
    assert not decision.graph_rejected(graph)
    state = dispatch.advance(graph)
    assert state['tasks']['unrelated'] == before
    assert state['tasks']['well_hom_sci_1_2_vts_pt00']['status'] == 'skipped'


def test_cancellation_is_owned_and_waited_for_before_publication(tmp_path, monkeypatch):
    worker, graph, path, request = setup(tmp_path)
    state = dispatch._load(graph)[2]
    state['tasks']['well_hom_sci_1_2_vts_pt00'].update(status='submitted', job_id='42')
    atomic_json(graph/'state.json', state)
    monkeypatch.setattr(dispatch, '_job_active', lambda _:True)
    owned = []
    monkeypatch.setattr(live_control, '_ownership', lambda d,j: owned.append((d,j)))
    cancelled = []
    monkeypatch.setattr(live_control, 'stop_attempt', lambda d,j,*args,**kw: cancelled.append(j) or False)
    assert not decision.apply(worker, {'well_hom_sci_1_2':'reject'}, path, request)
    assert cancelled == ['42'] and owned
    assert not decision.rejected(worker, 'well_hom_sci_1_2')
    monkeypatch.setattr(dispatch, '_job_active', lambda _:False)
    assert decision.apply(worker, {'well_hom_sci_1_2':'reject'}, path, request)
    assert decision.graph_rejected(graph)


@pytest.mark.parametrize('problem', ['no_cancel_authorization', 'foreign_allocation'])
def test_live_rejection_preflight_preserves_everything_on_failure(tmp_path, monkeypatch, problem):
    worker, graph, path, request = setup(tmp_path)
    state = dispatch._load(graph)[2]
    state['tasks']['well_hom_sci_1_2_vts_pt00'].update(status='submitted', job_id='42')
    atomic_json(graph/'state.json', state)
    monkeypatch.setattr(dispatch, '_job_active', lambda _: True)
    def ownership(*args):
        if problem == 'foreign_allocation':
            raise ValueError('Allocation belongs to another directory')
    monkeypatch.setattr(live_control, '_ownership', ownership)
    if problem == 'no_cancel_authorization':
        request['cancel_running'] = False
    def stop(*args, **kwargs):
        pytest.fail('Failed preflight must not cancel anything')
    monkeypatch.setattr(live_control, 'stop_attempt', stop)
    with pytest.raises(ValueError):
        decision.apply(worker, {'well_hom_sci_1_2':'reject'}, path, request)
    assert dispatch._load(graph)[2] == state
    assert not (worker/decision.NAME).exists()
    assert not (graph/decision.MARKER).exists()


def test_interrupted_marker_publication_finishes_without_losing_original_state(tmp_path, monkeypatch):
    worker, graph, path, request = setup(tmp_path)
    monkeypatch.setattr(dispatch, '_job_active', lambda _:False)
    original = dispatch._load(graph)[2]
    write = decision.atomic_json
    def interrupt(p, value):
        if p == graph/'state.json': raise OSError('interrupted')
        write(p, value)
    monkeypatch.setattr(decision, 'atomic_json', interrupt)
    with pytest.raises(OSError, match='interrupted'):
        decision.apply(worker, {'well_hom_sci_1_2':'reject'}, path, request)
    monkeypatch.setattr(decision, 'atomic_json', write)
    assert decision.apply(worker, {'well_hom_sci_1_2':'reject'}, path, request)
    marker = json.loads((graph/decision.MARKER).read_text())['reactions']['well_hom_sci_1_2']
    from pathlib import Path
    assert json.loads((Path(marker['archive'])/'state.json').read_text()) == original
    assert all(e['status'] == 'skipped' for e in dispatch._load(graph)[2]['tasks'].values())


def test_failed_chemistry_cannot_be_force_included(tmp_path):
    worker, *_ = setup(tmp_path)
    with pytest.raises(ValueError, match='previously accepted'):
        decision.validate(worker, {'failed_or_unknown':'include'})


def test_selection_excludes_manually_without_erasing_snapshot():
    from kinbot.pes_selection import choose
    snapshot = evidence('well', [('well_hom_sci_1_2',['p','q'],0.,None)])
    original = deepcopy(snapshot)
    selected = choose({'well':snapshot}, {'well':{'well_hom_sci_1_2':{'decision':'reject','reason':'Transfer'}}})
    assert selected['workers']['well']['selected'] == []
    assert selected['workers']['well']['omitted']['well_hom_sci_1_2']['reason'] == 'operator_rejected'
    assert snapshot == original
