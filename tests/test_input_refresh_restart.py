"""Accepted model inheritance across PES generations and live reconciliation."""
from copy import deepcopy
import json
import os
from types import SimpleNamespace

import pytest

from kinbot import input_refresh as refresh, pes_phase as phase, live_control as control
from kinbot.execution import atomic_json
from tests.test_input_refresh import run, change, surface


@pytest.fixture(autouse=True)
def isolated_phase(monkeypatch):
    # Explicit simulated Linux identity keeps owner tests portable to macOS.
    monkeypatch.setattr(phase, 'process_identity', lambda pid: dict(
        host='test-host', pid=pid, start='test-start', boot='test-boot') if pid == os.getpid() else None)
    for name in ('_ROOT', '_WORKER', '_EPOCH'):
        monkeypatch.setattr(phase, name, None)


def published_then_restarted(run):
    root, inp, original = run
    m = json.loads((root/phase.NAME).read_text())
    m['phase'] = 'exploring'
    for worker in m['workers'].values():
        worker['status'] = 'exploring'
    atomic_json(root/phase.NAME, m)
    updated = change(inp, original)
    first = refresh.submit(inp, apply=True)
    assert first['workers'] == 0
    phase.tick(root, ['123', '456'], max_searches=2)
    assert refresh.receipt(root)['status'] == 'applied'
    old = (root/refresh.JOURNAL).read_bytes()
    phase.initialize(root, inp)
    for name in ('123', '456'):
        phase.register(root, name, os.getpid())
    m = json.loads((root/phase.NAME).read_text())
    for worker in m['workers'].values():
        worker['status'] = 'vrc'
    atomic_json(root/phase.NAME, m)
    return root, inp, original, updated, first, old


def generator(root, name, original):
    species = SimpleNamespace(name=name, chemid=int(name), smiles='fixture',
        optical_reference={'status': 'assigned', 'canonical_graphs': []},
        reac_obj=[object()], reac_ts_done=[-1], reac_step=[5])
    return SimpleNamespace(inp=str(root/name/(name+'.json')),
                           par=deepcopy(original), species=species)


def test_accepted_models_survive_controller_generation_change(run):
    root, inp, original, updated, first, old = published_then_restarted(run)
    g = generator(root, '456', original)
    refresh.inherit_models(root, g)
    assert g.par['vrc_tst_surface_models'] == updated['vrc_tst_surface_models']
    assert json.loads((root/'456/.kinbot_search_options.json').read_text())[
        'vrc_tst_surface_models'] == updated['vrc_tst_surface_models']
    assert (root/refresh.JOURNAL).read_bytes() == old
    before = (root/'456/456.json').read_bytes()
    refresh.inherit_models(root, g)
    assert (root/'456/456.json').read_bytes() == before


def test_inheritance_uses_accepted_input_without_a_refresh_journal(run):
    root, inp, original = run
    updated = change(inp, original)
    phase.initialize(root, inp)
    assert not (root/refresh.JOURNAL).exists()
    # An edited but unsubmitted proposal must not enter a ready worker.
    proposed = deepcopy(updated)
    proposed['vrc_tst_surface_models']['789'] = surface()
    atomic_json(inp, proposed)
    g = generator(root, '456', original)
    refresh.inherit_models(root, g)
    assert g.par['vrc_tst_surface_models'] == updated['vrc_tst_surface_models']


def test_external_accepted_model_retains_relative_orbital_origin(run, tmp_path):
    root, _, original = run
    external = tmp_path/'inputs'; external.mkdir()
    inp = external/'case.json'
    accepted = deepcopy(original)
    stage = accepted['vrc_tst_surface_models']['123']['channels'][0]['bottlenecks'][0]
    stage['parameters']['vrc_tst_state_model_file'] = 'review/model.json'
    atomic_json(inp, accepted)
    phase.initialize(root, inp)
    g = generator(root, '456', dict(original, vrc_tst_surface_models={}))
    refresh.inherit_models(root, g)
    stage = g.par['vrc_tst_surface_models']['123']['channels'][0]['bottlenecks'][0]
    assert stage['parameters']['vrc_tst_state_model_file'] == str(external/'review/model.json')


def test_reconcile_applies_through_existing_owner_preserving_qc_and_policy(run, monkeypatch):
    from kinbot.search_control import revise_network
    root, inp, original, updated, first, old = published_then_restarted(run)
    native = root/'123/native.out'; native.write_bytes(b'accepted expensive result')
    live = root/'123/state.json'
    atomic_json(live, {'tasks': {'relaxed_scan': {'status': 'submitted', 'job_id': '908'}}})
    policy = root.parent/'.kinbot_queue_policy.json'
    atomic_json(policy, {'max_jobs': 23, 'excluded_nodes': ['unusable'], 'revision': 14})
    retained = {p: p.read_bytes() for p in (native, live, policy, inp)}
    # No process creation or scheduler command is part of reconciliation.
    monkeypatch.setattr(phase, 'owner_alive', lambda m: True)
    monkeypatch.setattr(phase, 'active_worker', lambda p: True)
    import subprocess
    monkeypatch.setattr(subprocess, 'Popen', lambda *a, **k: pytest.fail('started process'))
    preview = refresh.submit(inp, reconcile_models=True)
    assert preview['status'] == 'reviewed' and preview['workers'] == 2
    assert (root/refresh.JOURNAL).read_bytes() == old
    request = refresh.submit(inp, apply=True, reconcile_models=True)
    assert request['id'] != first['id']
    assert (root/'.kinbot_input_history'/(first['id']+'.json')).read_bytes() == old
    assert refresh.submit(inp, apply=True, reconcile_models=True)['id'] == request['id']
    phase.tick(root, ['123', '456'], max_searches=2)
    requests = control.requests(root, kinds={'network'})
    assert len(requests) == 2
    assert all(not r['cancel_running'] for _, r in requests)
    for name in ('123', '456'):
        g = generator(root, name, original)
        reaction = g.species.reac_obj[0]
        with monkeypatch.context() as context:
            context.chdir(root/name)
            revise_network(g)
        assert g.par['vrc_tst_surface_models'] == updated['vrc_tst_surface_models']
        assert g.species.reac_obj == [reaction] and g.species.reac_ts_done == [-1]
    assert refresh.receipt(root)['status'] == 'applied'
    assert refresh.submit(inp, apply=True, reconcile_models=True)['id'] == request['id']
    for path, content in retained.items():
        assert path.read_bytes() == content


def test_reconcile_refuses_unaccepted_edits_and_conflicting_models(run):
    root, inp, original, updated, first, old = published_then_restarted(run)
    proposal = deepcopy(updated)
    proposal['vrc_tst_surface_models']['789'] = surface()
    atomic_json(inp, proposal)
    with pytest.raises(ValueError, match='accepted input'):
        refresh.submit(inp, apply=True, reconcile_models=True)
    assert (root/refresh.JOURNAL).read_bytes() == old
    atomic_json(inp, updated)
    incompatible = deepcopy(original)
    incompatible['vrc_tst_surface_models']['123']['channels'][0]['bottlenecks'][0][
        'parameters']['vrc_tst_scan_basis'] = 'avtz-f12'
    atomic_json(root/'123/123.json', incompatible)
    with pytest.raises(ValueError, match='scientific migration'):
        refresh.submit(inp, apply=True, reconcile_models=True)
    assert (root/refresh.JOURNAL).read_bytes() == old


def test_inheritance_during_partial_publication_is_deferred(run, monkeypatch):
    root, inp, original = run
    change(inp, original)
    refresh.submit(inp, apply=True)
    real = refresh.atomic_json
    def interrupted(path, value):
        if str(path).endswith(refresh.JOURNAL) and value.get('status') == 'published':
            raise OSError('interrupted publication')
        real(path, value)
    monkeypatch.setattr(refresh, 'atomic_json', interrupted)
    with pytest.raises(OSError, match='interrupted publication'):
        phase.tick(root, ['123', '456'], max_searches=2)
    # A new worker has no pending network request, but must still wait until
    # the entire accepted-input publication transaction has finished.
    w = root/'789'; w.mkdir(); atomic_json(w/'789.json', original)
    g = generator(root, '789', original)
    refresh.inherit_models(root, g)
    assert g.par == original


def test_reconcile_waits_for_search_boundary_without_republishing(run):
    root, inp, original, updated, first, old = published_then_restarted(run)
    m = json.loads((root/phase.NAME).read_text())
    for worker in m['workers'].values():
        worker['status'] = 'exploring'
    atomic_json(root/phase.NAME, m)
    request = refresh.submit(inp, apply=True, reconcile_models=True)
    assert request['workers'] == 0
    phase.tick(root, ['123', '456'], max_searches=2)
    before = (root/refresh.JOURNAL).read_bytes()
    assert refresh.submit(inp, apply=True, reconcile_models=True)['id'] == request['id']
    assert (root/refresh.JOURNAL).read_bytes() == before
    g = generator(root, '456', original)
    refresh.inherit_models(root, g)
    assert g.par['vrc_tst_surface_models'] == updated['vrc_tst_surface_models']


def test_reconcile_includes_worker_becoming_ready_during_publication(run):
    root, inp, original, updated, first, old = published_then_restarted(run)
    m = json.loads((root/phase.NAME).read_text())
    m['workers']['456']['status'] = 'exploring'
    atomic_json(root/phase.NAME, m)
    request = refresh.submit(inp, apply=True, reconcile_models=True)
    assert request['workers'] == 1
    m['workers']['456']['status'] = 'vrc'
    atomic_json(root/phase.NAME, m)
    phase.tick(root, ['123', '456'], max_searches=2)
    assert len(control.requests(root, kinds={'network'})) == 2
