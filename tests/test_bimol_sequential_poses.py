"""Stop pose exploration at intended IRC evidence and resume failed pathways."""
import copy
import json
from pathlib import Path
from types import SimpleNamespace as NS
from unittest.mock import Mock

import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule

from kinbot import bimol, constants
from kinbot.reactions.reac_abstraction import Abstraction
from kinbot.stationary_pt import StationaryPoint


def reactants():
    points = []
    for formula, mult in (('CH4', 1), ('OH', 2)):
        atoms = molecule(formula)
        point = StationaryPoint(formula, 0, mult, atom=atoms.get_chemical_symbols(), geom=atoms.positions)
        point.characterize()
        points.append(point)
    species = StationaryPoint('pair', 0, 2, atom=list(points[0].atom) + list(points[1].atom),
        geom=np.vstack((points[0].geom, points[1].geom + [6., 0., 0.])),
        fragA=points[0], fragB=points[1])
    species.characterize()
    return species


def validation(species, expected=True):
    return dict(valid=True, product_matches_expected=expected,
        reactant_configuration_matches=True, endpoints=[
            dict(state='reactant', direction='forward', converged=True,
                 geometry_angstrom=species.geom.tolist()),
            dict(state='product', direction='reverse', converged=True,
                 geometry_angstrom=species.geom.tolist())])


def search_fixture(tmp_path, monkeypatch, validations):
    species = reactants()
    atoms = Atoms(species.atom, positions=species.geom)
    poses = lambda *args: ((atoms.copy(), {'pose': i}) for i in range(4))
    seeds = [dict(atoms=atoms, mode=np.ones(3 * species.natom), source='fixture')]
    force = Mock(return_value=(seeds, {'status': 'ok'}))
    monkeypatch.setattr(bimol, '_collision_free', lambda *args: True)
    monkeypatch.setattr(bimol, '_force_path', force)
    def refine(candidate, charge, mult, calculator, settings, directory):
        directory.mkdir(parents=True, exist_ok=True)
        refined = atoms.copy()
        pose = int(directory.parent.name[1:])
        seed = int(directory.name.split('_')[1])
        refined.positions[species.fragA.natom:] += [.5 * (pose + seed), 0., 0.]
        return dict(saddle=True, energy=-10., converged=True), refined
    monkeypatch.setattr(bimol, '_refine_candidate', refine)
    def validate(*args):
        Path(args[5]).mkdir(parents=True, exist_ok=True)
        return validations.pop(0)
    monkeypatch.setattr(bimol, 'validate_candidate', validate)
    def run(**kwargs):
        return bimol.find_ts_candidates(species, [species.fragA.geom], [species.fragB.geom],
            [5, 1, 0], [[5, 1]], [[1, 0]], poses, None, tmp_path, **kwargs)
    return species, run, force


def test_wrong_product_does_not_stop_search_and_resume_keeps_remaining_pose(tmp_path, monkeypatch):
    species = reactants()
    _, run, force = search_fixture(tmp_path, monkeypatch,
        [validation(species, False), validation(species), validation(species)])
    status, geometries, saved = run()
    assert status == 0 and len(geometries) == 1
    assert saved['next_cursor'] == [0, 0, 1, 0]
    assert len(saved['observations']) == 2
    assert force.call_count == 2
    _, _, resumed = run(start_after=saved['next_cursor'])
    assert resumed['next_cursor'] == [0, 0, 2, 0]
    # Saved seeds avoid repeating the force search at the previous pose.
    assert force.call_count == 3


def test_resume_uses_the_remaining_candidate_in_the_same_pose(tmp_path, monkeypatch):
    species = reactants()
    _, run, force = search_fixture(tmp_path, monkeypatch,
                                   [validation(species), validation(species)])
    seeds, path = force.return_value
    seeds.append(copy.deepcopy(seeds[0]))
    _, _, first = run()
    _, _, resumed = run(start_after=first['next_cursor'])
    assert first['next_cursor'] == [0, 0, 0, 0]
    assert resumed['next_cursor'] == [0, 0, 0, 1]
    assert force.call_count == 1


def test_wrong_configuration_does_not_stop_search(tmp_path, monkeypatch):
    species = reactants()
    _, run, _ = search_fixture(tmp_path, monkeypatch, [validation(species), validation(species)])
    monkeypatch.setattr(bimol, '_requested_reactant_configuration', Mock(side_effect=[False, True]))
    _, _, saved = run()
    assert saved['next_cursor'] == [0, 0, 1, 0]
    assert saved['observations'][0]['validation']['reactant_configuration_matches'] is False


def test_configured_reactant_gate_respects_specified_and_racemic_population():
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles('CC[C@H](O)C'))
    AllChem.EmbedMolecule(mol, randomSeed=13)
    fragment = StationaryPoint('butanol', 0, 1, atom=[a.GetSymbol() for a in mol.GetAtoms()],
                               geom=mol.GetConformer().GetPositions())
    fragment.characterize()
    species = reactants()
    species.fragA = fragment
    endpoint = np.vstack((fragment.geom * [-1., 1., 1.], species.fragB.geom + [20., 0., 0.]))
    evidence = dict(endpoints=[dict(state='reactant', geometry_angstrom=endpoint.tolist())])
    assert not bimol._requested_reactant_configuration(species, evidence, 'specified')
    assert bimol._requested_reactant_configuration(species, evidence, 'racemic')


def test_racemic_reactants_allow_global_reflection_not_one_fragment_only():
    from rdkit import Chem
    from rdkit.Chem import AllChem
    fragments = []
    for smiles in ('CC[C@H](O)C', 'CC[C@H](F)C'):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol, randomSeed=13)
        point = StationaryPoint('configured', 0, 1, atom=[a.GetSymbol() for a in mol.GetAtoms()],
                                geom=mol.GetConformer().GetPositions())
        point.characterize()
        fragments.append(point)
    species = NS(fragA=fragments[0], fragB=fragments[1])
    for reflect_second, expected in ((False, False), (True, True)):
        second = fragments[1].geom * ([-1., 1., 1.] if reflect_second else [1., 1., 1.])
        endpoint = np.vstack((fragments[0].geom * [-1., 1., 1.], second + [20., 0., 0.]))
        evidence = dict(endpoints=[dict(state='reactant', geometry_angstrom=endpoint.tolist())])
        assert bimol._requested_reactant_configuration(species, evidence, 'racemic') == expected


def controller(tmp_path, monkeypatch, count=2):
    monkeypatch.chdir(tmp_path)
    species = reactants()
    par = dict(scan_step=20, bimol_search_kwargs={}, conformer_search=1,
               bimol_ts_rmsd=.1, bimol_energy_tolerance=.25, multi_conf_tst=0)
    rows, statuses = {}, {}
    qc = NS(db=NS(select=lambda name: [NS(data=rows[name])] if name in rows else []),
            check_qc=lambda name: statuses.get(name, 0), qc_bimol_search=Mock())
    for field in bimol._REACTION_FIELDS:
        setattr(species, field, [])
    for index in range(count):
        reaction = Abstraction(species, qc, par, [5, 1, 0], f'path_{index}')
        reaction.stereopath_id = f'stereo_{index}'
        values = dict(reac_type='abstraction', reac_inst=[5, 1, 0], reac_step=0,
                      reac_scan_energy=[], reac_ts_done=0, reac_ts_geom=0, reac_ts_freq=0,
                      reac_name=reaction.instance_name, reac_obj=reaction)
        for field in bimol._REACTION_FIELDS:
            getattr(species, field).append(values[field])
    bimol.prepare_bimolecular_reactions(species, par, qc)
    def finish(path, attempt=0, pose=0, geom=None, energy=-10.):
        job = bimol.search_job(f'path_{path}', attempt, sequential=True)
        candidate = dict(geometry=(species.geom if geom is None else geom).tolist(),
            directory=f'source/path_{path}/pose_{pose}',
            result={'energy': energy}, cursor=[0, 0, pose, 0], validation=validation(species),
            validation_source=f'source/path_{path}/pose_{pose}/validation')
        rows[job] = {'search': dict(candidates=[candidate])}
        statuses[job] = 'normal'
    return species, par, qc, finish, statuses


@pytest.mark.parametrize('failure', ['refinement', 'hir_refinement', 'final_irc',
                                     'final_barrier_rejected', 'ts_frequency_failed'])
def test_failed_path_resumes_without_repeating_successful_path(tmp_path, monkeypatch, failure):
    species, par, qc, finish, statuses = controller(tmp_path, monkeypatch)
    finish(0); finish(1)
    assert not bimol.advance_pose_searches(species, par, qc)
    assert len(species.reac_obj) == 4
    assert species.reac_obj[2].stereopath_id != species.reac_obj[3].stereopath_id
    failed = species.reac_obj[2]
    if failure in {'final_barrier_rejected', 'ts_frequency_failed'}:
        failed.final_irc_processed = True
        failed.bimol_final_status = failure
    if failure == 'final_irc':
        failed.final_irc_attempted = True
        failed.final_ts_cluster_id = 'old_cluster'
        failed.final_ts_rank = 0
    species.reac_ts_done[2:] = [-999, -1]
    assert bimol.advance_pose_searches(species, par, qc)
    assert qc.qc_bimol_search.call_args.kwargs == dict(attempt=1, start_after=[0, 0, 0, 0])
    # A genuinely different saddle can still replace the failed attempt.
    geom = species.geom.copy()
    geom[5:] += [1., 0., 0.]
    finish(0, attempt=1, pose=1, geom=geom)
    assert not bimol.advance_pose_searches(species, par, qc)
    assert species.reac_ts_done == [-998, -998, -999, -1, 0]
    assert species.reac_obj[4].instance_name != failed.instance_name
    saved = json.loads(Path('bimol/validated_l0/index.json').read_text())['validated_l0']
    assert len(saved) == 3
    # A new ready TS must reach IRC despite an older failed cluster.
    if failure == 'final_irc':
        point = copy.deepcopy(species)
        point.energy, point.zpe = -10., 0.
        species.reac_obj[3].final_ts_cluster_id = 'successful_cluster'
        species.reac_obj[3].final_ts_rank = 0
        species.reac_obj[3].ts = point
        species.reac_obj[4].ts = point
        species.reac_ts_done[4] = 6
        bimol.release_final_ts_clusters(species, par)
        assert species.reac_ts_done[4] == 7


def test_pending_search_keeps_controller_alive_and_exhaustion_finishes(tmp_path, monkeypatch):
    species, par, qc, finish, statuses = controller(tmp_path, monkeypatch, count=1)
    job = bimol.search_job('path_0', sequential=True)
    statuses[job] = 'running'
    assert bimol.advance_pose_searches(species, par, qc)
    statuses[job] = 'error'
    assert not bimol.advance_pose_searches(species, par, qc)
    assert species.bimol_pose_searches[0]['done']


def test_post_irc_product_failure_does_not_repeat_ts_search(tmp_path, monkeypatch):
    species, par, qc, finish, statuses = controller(tmp_path, monkeypatch, count=1)
    finish(0)
    bimol.advance_pose_searches(species, par, qc)
    species.reac_ts_done[1] = -999
    species.reac_obj[1].final_irc_processed = True
    qc.qc_bimol_search.reset_mock()
    assert not bimol.advance_pose_searches(species, par, qc)
    qc.qc_bimol_search.assert_not_called()


def test_final_wrong_product_is_saved_and_reopens_the_pose_search(tmp_path, monkeypatch):
    from kinbot import reaction_generator
    species, par, qc, finish, statuses = controller(tmp_path, monkeypatch, count=1)
    par['bimol'] = 1
    finish(0)
    bimol.advance_pose_searches(species, par, qc)
    reaction = species.reac_obj[1]
    species.reac_ts_done[1] = 7
    result = NS(product_matches_expected=False, product=NS(geom=species.geom.copy()))
    monkeypatch.setattr(reaction_generator, 'IRC', Mock(return_value=NS(
        check_irc=lambda: ['normal', 'normal'], irc2stationary_pt=lambda: result)))
    generator = reaction_generator.ReactionGenerator(species, par, qc, 'input.json')
    generator._advance_irc(1, reaction, 10.)
    assert species.reac_ts_done[1] == -999
    assert reaction.bimol_final_status == 'final_irc_wrong_product'
    saved = json.loads(Path('bimol/validated_l0/index.json').read_text())['validated_l0'][0]
    assert saved['final_status'] == 'final_irc_wrong_product'
    np.testing.assert_allclose(saved['observed_product_geometry'], result.product.geom)
    assert bimol.advance_pose_searches(species, par, qc)
    assert qc.qc_bimol_search.call_args.kwargs == dict(attempt=1, start_after=[0, 0, 0, 0])


def test_failed_final_irc_ignores_old_endpoint_and_reopens_pose_search(tmp_path, monkeypatch):
    from kinbot import reaction_generator
    species, par, qc, finish, statuses = controller(tmp_path, monkeypatch, count=1)
    par['bimol'] = 1
    finish(0)
    bimol.advance_pose_searches(species, par, qc)
    reaction = species.reac_obj[1]
    species.reac_ts_done[1] = 7
    old_endpoint = Mock(side_effect=AssertionError('Failed IRC must not reuse an old endpoint'))
    monkeypatch.setattr(reaction_generator, 'IRC', Mock(return_value=NS(
        check_irc=lambda: ['normal', 'error'], irc2stationary_pt=old_endpoint)))
    generator = reaction_generator.ReactionGenerator(species, par, qc, 'input.json')
    generator._advance_irc(1, reaction, 10.)
    old_endpoint.assert_not_called()
    assert species.reac_ts_done[1] == -999
    assert reaction.bimol_final_status == 'final_irc_failed'
    assert bimol.advance_pose_searches(species, par, qc)
    assert qc.qc_bimol_search.call_args.kwargs == dict(attempt=1, start_after=[0, 0, 0, 0])


@pytest.mark.parametrize('change', ['same', 'rotate_translate', 'equivalent_hydrogens'])
def test_equivalent_retry_skips_l1_and_preserves_failed_evidence(tmp_path, monkeypatch,
                                                              caplog, change):
    species, par, qc, finish, statuses = controller(tmp_path, monkeypatch, count=1)
    finish(0)
    bimol.advance_pose_searches(species, par, qc)
    first = species.reac_obj[1]
    species.reac_ts_done[1] = -999
    bimol.update_bimol_status(first, l1_valid=False, final_status='l1_refinement_rejected')
    index = Path('bimol/validated_l0/index.json')
    before = index.read_bytes()
    bimol.advance_pose_searches(species, par, qc)
    geom = species.geom.copy()
    if change == 'rotate_translate':
        geom = geom[:, [1, 2, 0]] + [4., -2., 1.]
    elif change == 'equivalent_hydrogens':
        geom[[2, 3]] = geom[[3, 2]]
    finish(0, attempt=1, pose=1, geom=geom)
    with caplog.at_level('INFO', logger='KinBot'):
        assert bimol.advance_pose_searches(species, par, qc)
    assert 'already promoted to L1' in caplog.text
    assert species.reac_ts_done == [-998, -999]
    assert index.read_bytes() == before
    assert len(list(Path('bimol/unique_candidates').glob('candidate_*'))) == 1
    state = species.bimol_pose_searches[0]
    assert state['active'] is None and not state['done']
    assert state['attempt'] == 2 and state['cursor'] == [0, 0, 1, 0]
    # The saved L0 observation remains available even though it was not promoted.
    job = bimol.search_job('path_0', attempt=1, sequential=True)
    assert bimol._job_data(qc, job)['search']['candidates'][0]['cursor'] == [0, 0, 1, 0]
    assert bimol.advance_pose_searches(species, par, qc)
    assert qc.qc_bimol_search.call_args.kwargs == dict(attempt=2, start_after=[0, 0, 1, 0])
    statuses[bimol.search_job('path_0', attempt=2, sequential=True)] = 'error'
    assert not bimol.advance_pose_searches(species, par, qc)
    assert state['done']


@pytest.mark.parametrize('difference, duplicate', [(-.1, True), (.1, True), (.5, False)])
def test_retry_uses_existing_l0_energy_tolerance(tmp_path, monkeypatch, difference, duplicate):
    species, par, qc, finish, statuses = controller(tmp_path, monkeypatch, count=1)
    finish(0)
    bimol.advance_pose_searches(species, par, qc)
    species.reac_ts_done[1] = -999
    bimol.advance_pose_searches(species, par, qc)
    energy = -10. + difference * constants.KCALtoHARTREE / constants.EVtoHARTREE
    finish(0, attempt=1, pose=1, energy=energy)
    assert bimol.advance_pose_searches(species, par, qc) == duplicate
    assert len(species.reac_obj) == (2 if duplicate else 3)


@pytest.mark.parametrize('allow_reflection', [False, True])
def test_retry_respects_existing_reflection_policy(tmp_path, monkeypatch, allow_reflection):
    species, par, qc, finish, statuses = controller(tmp_path, monkeypatch, count=1)
    par['bimol_ts_rmsd'] = .01
    species.reac_obj[0].allow_reflection = allow_reflection
    geom = species.geom + np.random.default_rng(7).normal(0., .13, species.geom.shape)
    finish(0, geom=geom)
    bimol.advance_pose_searches(species, par, qc)
    species.reac_ts_done[1] = -999
    bimol.advance_pose_searches(species, par, qc)
    finish(0, attempt=1, pose=1, geom=geom * [-1., 1., 1.])
    assert bimol.advance_pose_searches(species, par, qc) == allow_reflection
    assert len(species.reac_obj) == (2 if allow_reflection else 3)


def test_retry_checks_all_previous_l0_saddles_and_replays_saved_jobs(tmp_path, monkeypatch):
    species, par, qc, finish, statuses = controller(tmp_path, monkeypatch, count=1)
    finish(0)
    different = species.geom.copy()
    different[5:] += [1., 0., 0.]
    finish(0, attempt=1, pose=1, geom=different)
    finish(0, attempt=2, pose=2)  # returns to the first saddle, not the second
    finish(0, attempt=3, pose=3)  # another duplicate must also advance
    finish(0, attempt=4, pose=4, geom=different)
    statuses[bimol.search_job('path_0', attempt=5, sequential=True)] = 'error'
    original = {field: list(getattr(species, field)) for field in bimol._REACTION_FIELDS}
    for restart in (False, True):
        if restart:
            # On restart the ordinary controller reconstructs promotions from
            # completed QC jobs, without needing another saved saddle registry.
            for field, values in original.items():
                setattr(species, field, list(values))
            bimol._initialize_pose_searches(species, par)
        assert not bimol.advance_pose_searches(species, par, qc)
        species.reac_ts_done[1] = -999  # completed first L1 failure
        assert not bimol.advance_pose_searches(species, par, qc)
        species.reac_ts_done[2] = -999  # completed distinct L1 failure
        for attempt in (3, 4, 5):
            assert bimol.advance_pose_searches(species, par, qc)
            assert species.bimol_pose_searches[0]['attempt'] == attempt
            assert len(species.reac_obj) == 3
        assert not bimol.advance_pose_searches(species, par, qc)
        assert species.bimol_pose_searches[0]['done']
    qc.qc_bimol_search.assert_called_once()  # initial submission from the fixture
