"""A controller replacement must not create another reaction-search attempt."""
from datetime import datetime, timezone
import json
import os
from pathlib import Path
from types import SimpleNamespace as NS

import pytest

from kinbot import search_state
from kinbot.execution import atomic_json
from kinbot.qc import QuantumChemistry


@pytest.fixture
def tree(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    inp = tmp_path/'input.json'
    atomic_json(inp, {'smiles': 'CC'})
    def make(status=0, step=0, job='41', family='test', scan=False):
        obj = NS(instance_name='reaction', family_name=family, instance=[0, 1, 2],
                 max_step=12, scan=scan, skip=False, products=[])
        species = NS(reac_obj=[obj], reac_ts_done=[status], reac_step=[step],
                     reac_scan_energy=[[]], reac_type=[family], reac_inst=[[0, 1, 2]])
        qc = NS(job_ids={'reaction': job} if job else {})
        generator = NS(species=species, qc=qc, inp=str(inp))
        return generator
    return NS(root=tmp_path, inp=inp, make=make)


def checkpoint(g):
    g.search_checkpoint = search_state.SearchState(g)
    return g.search_checkpoint


def test_controller_rejection_survives_successful_intermediate_qc(tree):
    g = tree.make()
    cp = checkpoint(g)
    g.species.reac_step[0] = -1
    g.species.reac_ts_done[0] = -999
    cp.save(0)
    for _ in range(3):
        restored = tree.make()
        restored.qc.check_qc = lambda *a: pytest.fail('QC success must not reopen a rejected search')
        checkpoint(restored).save(0)
        assert restored.species.reac_ts_done == [-999]
        assert restored.species.reac_step == [-1]


def test_generate_keeps_imported_rejection_without_calling_qc(tree, monkeypatch):
    from kinbot.parameters import Parameters
    from kinbot.reaction_generator import ReactionGenerator
    (tree.root/'kinbot_monitor.out').write_text('-999 -1 reaction\n')
    g = tree.make()
    g.qc.check_qc = lambda *a: pytest.fail('Rejected search was reopened')
    atomic_json(tree.inp, {'smiles': 'CC', 'barrier_threshold': 100})
    par = Parameters(str(tree.inp)).par
    par['delete_intermediate_files'] = 0
    generator = ReactionGenerator(g.species, par, g.qc, str(tree.inp))
    monkeypatch.setattr('kinbot.reaction_generator.time.sleep', lambda *a: None)
    generator.generate()
    assert g.species.reac_ts_done == [-999]
    assert (tree.root/'kinbot_monitor.out').read_text().split() == ['-999', '-1', 'reaction']


@pytest.mark.parametrize('stage', [1, 2, 3, 4, 5, -1])
def test_completed_ts_search_is_not_repeated_while_rebuilding_objects(tree, stage):
    g = tree.make()
    cp = checkpoint(g)
    g.species.reac_ts_done[0] = stage
    g.species.reac_step[0] = 13
    cp.save(0)
    restored = tree.make()
    cp = checkpoint(restored)
    assert restored.species.reac_ts_done == [1]
    assert restored.species.reac_step == [13]
    if stage == -1:
        cp.enter(0)
        assert restored.qc._search_replay == 'reaction'
        cp.save(0)  # another refresh in the middle of object reconstruction
        assert json.loads(cp.path.read_text())['reactions']['reaction']['status'] == -1
        again = tree.make()
        assert checkpoint(again).replaying == {'reaction'}
        cp.leave()
        assert restored.qc._search_replay is None


def test_scission_rebuilds_fragments_without_saddle_search(tree):
    g = tree.make(family='hom_sci')
    g.species.reac_obj[0].family_name = 'hs'
    cp = checkpoint(g)
    g.species.reac_ts_done[0] = -1
    cp.save(0)
    g = tree.make(family='hom_sci')
    g.species.reac_obj[0].family_name = 'hs'
    cp = checkpoint(g)
    assert g.species.reac_ts_done == [0]
    assert cp.replaying == {'reaction'}


@pytest.mark.parametrize('acknowledged', [False, True])
def test_path_submission_acknowledgement_distinguishes_point_from_allocation_retry(tree, acknowledged):
    g = tree.make(step=5)
    cp = checkpoint(g)
    cp.prepare_submission(g.species.reac_obj[0], 5)
    g.qc.job_ids['reaction'] = '42'
    if acknowledged:
        search_state.submitted(g.qc, 'reaction', '42')
        cp.save(0)  # exception/finally before carry_out_reaction returned
    restored = tree.make(job='42')
    checkpoint(restored)
    assert restored.species.reac_step == [6 if acknowledged else 5]
    again = tree.make(job='42')
    checkpoint(again)
    assert again.species.reac_step == [6 if acknowledged else 5]  # advance once only


def test_waiting_for_pool_does_not_advance_path(tree):
    g = tree.make(step=5)
    cp = checkpoint(g)
    cp.prepare_submission(g.species.reac_obj[0], 5)
    restored = tree.make()
    checkpoint(restored)
    assert restored.species.reac_step == [5]


def test_scan_history_survives_restart_and_does_not_duplicate_current_point(tree):
    g = tree.make(step=2, scan=True)
    cp = checkpoint(g)
    cp.scan_energy(0, -1.0)
    cp.prepare_submission(g.species.reac_obj[0], 2)
    restored = tree.make(scan=True)
    cp = checkpoint(restored)
    cp.scan_energy(0, -1.0)
    assert restored.species.reac_scan_energy == [[-1.0]]
    restored.qc.job_ids['reaction'] = '42'
    restored.species.reac_step[0] = 3
    cp.scan_energy(0, -0.9)
    assert restored.species.reac_scan_energy == [[-1.0, -0.9]]


def test_explicit_retry_archives_search_decision_and_resets_only_target(tree):
    g = tree.make()
    cp = checkpoint(g)
    g.species.reac_ts_done[0] = -999
    cp.save(0)
    cp.entries['unrelated'] = {'status': -999, 'step': 4}
    g.species.reac_ts_done[0] = 0
    g.species.reac_step[0] = 0
    cp.reset(0, 'operator-request')
    entry = cp.entries['reaction']
    assert entry['status'] == 0 and entry['step'] == 0
    assert entry['history'][0]['previous']['status'] == -999
    assert cp.entries['unrelated'] == {'status': -999, 'step': 4}
    restored = tree.make()
    checkpoint(restored)
    assert restored.species.reac_ts_done == [0]


@pytest.mark.parametrize('change', ['input', 'profile', 'reaction'])
def test_incompatible_restart_does_not_silently_reopen_or_reuse_search(tree, change):
    checkpoint(tree.make())
    g = tree.make()
    if change == 'input':
        atomic_json(tree.inp, {'smiles': 'CCC'})
    elif change == 'profile':
        atomic_json(tree.root/'.kinbot_theory_jobs.json', {'profile_sha256': 'different'})
    else:
        g.species.reac_obj[0].max_step += 1
    with pytest.raises(RuntimeError, match='changed'):
        checkpoint(g)


def test_explicit_policy_revision_revalidates_without_implicitly_retrying_failures(tree):
    g = tree.make()
    cp = checkpoint(g)
    g.species.reac_ts_done[0] = -1
    cp.save(0)
    cp.entries['rejected'] = {'status': -999, 'step': 4}
    cp.write()
    atomic_json(tree.inp, {'smiles': 'CC', 'vdW_detection': 2.0})
    search_state.rebase(tree.inp, 'policy-change', reassess=True)
    g = tree.make()
    cp = checkpoint(g)
    assert not cp.replaying  # explicitly changed policy can need additional QC
    assert cp.entries['rejected']['status'] == -999


@pytest.mark.parametrize('new_submission', [False, True])
def test_legacy_monitor_uses_admission_not_worker_render_time(tree, new_submission):
    from kinbot import run_policy
    monitor = tree.root/'kinbot_monitor.out'
    monitor.write_text('0 5 reaction\n-999 -1 rejected\n')
    os.utime(monitor, (1000, 1000))
    (tree.root/'reaction.py').write_text('newly rendered input while waiting\n')
    atomic_json(tree.root/'.kinbot_theory_jobs.json', {
        'jobs': {'reaction': {'job_id': '41'}}})
    atomic_json(tree.root/run_policy.NAME, {})
    atomic_json(tree.root/'.kinbot_policy_submissions/41.json', dict(
        directory=str(tree.root), job_id='41',
        command=['sbatch', 'reaction.sbatch'],
        created_at=datetime.fromtimestamp(1001 if new_submission else 999, timezone.utc).isoformat()))
    search_state.capture_legacy(tree.inp)
    g = tree.make()
    checkpoint(g)
    assert g.species.reac_step == [6 if new_submission else 5]
    assert json.loads((tree.root/search_state.NAME).read_text())['reactions']['rejected']['status'] == -999


def test_accepted_reconstruction_cannot_submit_missing_or_replacement_qc():
    qc = NS(_search_replay='accepted', check_qc=lambda job: 0)
    with pytest.raises(RuntimeError, match='no calculation was submitted'):
        QuantumChemistry._submit_qc(qc, 'product_well', 1)
    qc.check_qc = lambda job: 'normal'
    assert QuantumChemistry._submit_qc(qc, 'product_well', 1) == 0
    with pytest.raises(RuntimeError, match='no calculation was submitted'):
        QuantumChemistry._submit_qc(qc, 'reaction', 1, singlejob=0)


def test_accepted_reconstruction_guards_both_theory_backends(tree):
    g = tree.make()
    cp = checkpoint(g)
    g.species.reac_ts_done[0] = -1
    cp.save(0)
    g = tree.make()
    g.qc.backends = {'l1': NS(), 'l2': NS()}
    cp = checkpoint(g)
    cp.enter(0)
    assert all(b._search_replay == 'reaction' for b in g.qc.backends.values())
    cp.leave()
    assert all(b._search_replay is None for b in g.qc.backends.values())


def test_resuming_step_zero_never_rewrites_a_running_worker(tree):
    from kinbot.reac_family import carry_out_reaction
    worker = tree.root/'reaction.py'
    worker.write_text('immutable queued worker input\n')
    reaction = NS(instance_name='reaction', qc=NS(check_qc=lambda name: 'running'))
    assert carry_out_reaction(reaction, 0, 'backend') == 0
    assert worker.read_text() == 'immutable queued worker input\n'


@pytest.mark.parametrize('source', ['snapshot', 'legacy_log'])
def test_restore_embedding_requires_exact_full_input_hash_and_preserves_search(tree, source):
    original = {'smiles': '', 'structure': ['H', 0., 0., 0.], 'title': '123', 'barrier_threshold': 100.}
    atomic_json(tree.inp, original)
    cp = checkpoint(tree.make())
    cp.entries['reaction']['status'] = -999
    if source == 'legacy_log':
        cp.record.pop('input_snapshot')
    cp.write()
    saved = cp.path.read_bytes()
    changed = dict(original, structure=['H', 1., 2., 3.])
    atomic_json(tree.inp, changed)
    before = tree.inp.read_bytes()
    log = tree.root/'archived.log'
    log.write_text('INFO: Input parameters\n\tstructure: '+repr(original['structure'])+
                   '\n\ttitle: 123\nINFO: Starting KinBot\n')
    plan = search_state.input_restore_plan(tree.inp, [log])
    assert plan['input'] == original
    assert tree.inp.read_bytes() == before  # planning does not edit files
    archive = tree.root/'archive/input.json'
    search_state.restore_input(tree.inp, plan, archive)
    assert archive.read_bytes() == before
    assert cp.path.read_bytes() == saved
    assert json.loads(tree.inp.read_text()) == original
    search_state.restore_input(tree.inp, plan, archive)  # interrupted handoff
    assert checkpoint(tree.make()).generator.species.reac_ts_done == [-999]


@pytest.mark.parametrize('change', ['theory', 'policy', 'missing_evidence', 'changed_after_plan'])
def test_embedding_recovery_never_rebases_unproven_changes(tree, change):
    original = {'structure': ['H', 0., 0., 0.], 'barrier_threshold': 100.}
    atomic_json(tree.inp, original)
    cp = checkpoint(tree.make())
    atomic_json(tree.inp, dict(original, structure=['H', 1., 2., 3.]))
    if change == 'theory':
        atomic_json(tree.root/'.kinbot_theory_jobs.json', {'profile_sha256': 'different'})
    elif change == 'policy':
        atomic_json(tree.inp, {'structure': ['H', 1., 2., 3.], 'barrier_threshold': 50.})
    elif change == 'missing_evidence':
        cp.record.pop('input_snapshot'); cp.write()
    if change == 'changed_after_plan':
        plan = search_state.input_restore_plan(tree.inp)
        atomic_json(tree.inp, dict(original, structure=['H', 4., 5., 6.]))
        with pytest.raises(RuntimeError, match='changed since recovery planning'):
            search_state.restore_input(tree.inp, plan, tree.root/'archive.json')
        assert not (tree.root/'archive.json').exists()
    else:
        before = tree.inp.read_bytes()
        with pytest.raises(RuntimeError, match='inputs/theory changed'):
            search_state.input_restore_plan(tree.inp)
        assert tree.inp.read_bytes() == before


def test_pes_restart_retains_input_despite_new_smiles_embedding(tmp_path, monkeypatch):
    from kinbot import pes
    from kinbot.stationary_pt import StationaryPoint
    from kinbot.species_routing import routing_name
    import numpy as np
    monkeypatch.chdir(tmp_path)
    source = tmp_path/'input.json'
    atomic_json(source, {'smiles': 'C', 'barrier_threshold': 100.})
    species = StationaryPoint('initial', 0, 1, smiles='C')
    species.characterize()
    pes.write_input(str(source), species, 100., None, str(tmp_path), 0)
    inp = tmp_path/routing_name(species)/(routing_name(species)+'.json')
    (inp.parent/'kinbot_monitor.out').write_text('-999 -1 rejected\n')
    search_state.capture_legacy(inp)
    original = inp.read_bytes()
    # A rigid reorientation is the same species, but a different raw JSON hash.
    species.geom = np.asarray(species.geom) @ np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    for _ in range(2):
        pes.write_input(str(source), species, 100., None, str(tmp_path), 0)
        assert inp.read_bytes() == original
        search_state.validate_input(inp)
        pes.write_input_keep(str(source), routing_name(species), str(tmp_path))
        assert inp.read_bytes() == original
