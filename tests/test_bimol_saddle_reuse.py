"""Repeated L0 saddles reuse completed checks without losing reaction paths."""
import json
from pathlib import Path
from unittest.mock import Mock

import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule

from kinbot import bimol, constants
from kinbot.bimol_network import BimolecularNetwork
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.reactions.reac_abstraction import Abstraction
from kinbot.stationary_pt import StationaryPoint
from tests.test_bimol_sequential_poses import (controller, reactants, search_fixture,
                                              validation)


def repeated_search(tmp_path, monkeypatch, checks, geometries=None, energies=None):
    species, run, force = search_fixture(tmp_path, monkeypatch, checks)
    geometries = geometries or [species.geom] * 4
    energies = energies or [-10.] * 4

    def refine(candidate, charge, mult, calculator, settings, directory):
        directory.mkdir(parents=True, exist_ok=True)
        pose = int(directory.parent.name[1:])
        return dict(saddle=True, energy=energies[pose], converged=True), Atoms(
            species.atom, positions=geometries[pose])

    monkeypatch.setattr(bimol, '_refine_candidate', refine)
    check = Mock(wraps=bimol.validate_candidate)
    monkeypatch.setattr(bimol, 'validate_candidate', check)
    return species, run, check


def rejected_validation(species):
    result = validation(species, False)
    result.update(valid=False, reactant_configuration_matches=False)
    return result


def test_repeated_rejected_saddle_skips_irc_but_distinct_saddle_succeeds(tmp_path, monkeypatch):
    species = reactants()
    shifted = species.geom.copy()
    shifted[5:] += [1., 0., 0.]
    rotated = species.geom[:, [1, 2, 0]] + [3., -2., 1.]
    _, run, check = repeated_search(tmp_path, monkeypatch,
        [rejected_validation(species), validation(species)],
        [species.geom, rotated, shifted, shifted])
    status, _, saved = run()
    assert status == 0 and saved['next_cursor'] == [0, 0, 2, 0]
    assert check.call_count == 2
    first, duplicate, accepted = saved['observations']
    assert duplicate['irc_reused_from'] == first['directory']
    assert 'validation' not in duplicate  # Do not relabel permuted endpoint atoms.
    assert accepted['validation']['valid']


@pytest.mark.parametrize('problem', ['endpoint', 'irc_exception', 'traceback'])
def test_incomplete_check_does_not_exclude_a_later_saddle(tmp_path, monkeypatch, problem):
    species = reactants()
    failed = rejected_validation(species)
    if problem == 'endpoint':
        failed['endpoints'][0]['converged'] = False
    elif problem == 'irc_exception':
        failed['endpoints'][0]['irc_error'] = 'calculator failed'
    else:
        failed['traceback'] = 'validation interrupted'
    _, run, check = repeated_search(tmp_path, monkeypatch, [failed, validation(species)])
    assert run()[2]['next_cursor'] == [0, 0, 1, 0]
    assert check.call_count == 2


def test_completed_checks_survive_retry_and_skip_previously_accepted_saddle(tmp_path, monkeypatch):
    species = reactants()
    first = validation(species)
    for endpoint in first['endpoints']:
        endpoint['irc_converged'] = False  # A point limit is not an exception.
    _, run, check = repeated_search(tmp_path, monkeypatch, [first])
    accepted = run()[2]
    status, _, resumed = run(start_after=accepted['next_cursor'])
    assert status == -1 and resumed['exhausted']
    assert check.call_count == 1
    assert len(resumed['observations']) == 4
    assert all(item['irc_reused_from'] == accepted['observations'][0]['directory']
               for item in resumed['observations'][1:])


def test_retry_ignores_checks_beyond_requested_cursor(tmp_path, monkeypatch):
    species = reactants()
    _, run, check = repeated_search(tmp_path, monkeypatch, [validation(species)])
    Path(tmp_path, 'search.json').write_text(json.dumps({'observations': [dict(
        geometry=species.geom.tolist(), result={'energy': -10.}, cursor=[0, 0, 3, 0],
        directory='future', validation=validation(species))]}))
    assert run(start_after=[0, 0, 0, 0])[2]['next_cursor'] == [0, 0, 1, 0]
    assert check.call_count == 1


def test_skipped_saddles_do_not_extend_the_energy_tolerance(tmp_path, monkeypatch):
    species = reactants()
    ev = constants.KCALtoHARTREE / constants.EVtoHARTREE
    _, run, check = repeated_search(tmp_path, monkeypatch,
        [rejected_validation(species), validation(species)],
        energies=[-10., -10. + .2 * ev, -10. + .4 * ev, -10.])
    assert run()[2]['next_cursor'] == [0, 0, 2, 0]
    assert check.call_count == 2


@pytest.mark.parametrize('allow_reflection', [False, True])
def test_worker_uses_reflection_policy(tmp_path, monkeypatch, allow_reflection):
    species = reactants()
    geom = species.geom + np.random.default_rng(7).normal(0., .13, species.geom.shape)
    _, run, check = repeated_search(tmp_path, monkeypatch,
        [rejected_validation(species), validation(species)],
        geometries=[geom, geom * [-1., 1., 1.], geom, geom])
    status, _, _ = run(search_kwargs=dict(allow_reflection=allow_reflection, bimol_ts_rmsd=.01))
    assert check.call_count == (1 if allow_reflection else 2)
    assert status == (-1 if allow_reflection else 0)


def equivalent_sites(tmp_path, monkeypatch):
    species, par, qc, finish, statuses = controller(tmp_path, monkeypatch)
    for reaction in species.reac_obj:
        reaction.stereopath_id = 'same_hydrogen_class'
        reaction.site_degeneracy = 4
    original_finish = finish

    def finish(path, **kwargs):
        original_finish(path, **kwargs)
        job = bimol.search_job(f'path_{path}', kwargs.get('attempt', 0), sequential=True)
        saved = bimol._job_data(qc, job)['search']['candidates'][0]
        for endpoint in saved['validation']['endpoints']:
            endpoint['configured_identity'] = dict(status='assigned', id=endpoint['state'])
        return saved

    return species, par, qc, finish, statuses


@pytest.mark.parametrize('state', ['active', 'successful', 'failed'])
def test_equivalent_site_uses_existing_l1_and_keeps_both_retry_cursors(
        tmp_path, monkeypatch, state):
    species, par, qc, finish, statuses = equivalent_sites(tmp_path, monkeypatch)
    finish(0, pose=1)
    bimol.advance_pose_searches(species, par, qc)
    representative = species.reac_obj[2]
    species.reac_ts_done[2] = {'active': 0, 'successful': -1, 'failed': -999}[state]
    finish(1, pose=2)
    assert bimol.advance_pose_searches(species, par, qc)
    assert len(species.reac_obj) == 3
    assert {member['reaction'] for member in representative.search_members} == {'path_0', 'path_1'}
    assert representative.site_degeneracy == 4
    assert BimolecularNetwork._site_degeneracy(None, [representative]) == 4
    record = json.loads(Path('bimol/validated_l0/index.json').read_text())['validated_l0'][0]
    assert record['members'] == representative.search_members
    if state == 'active':
        assert [item['active'] for item in species.bimol_pose_searches] == [2, 2]
        species.reac_ts_done[2] = -999
    elif state == 'successful':
        bimol.advance_pose_searches(species, par, qc)
        assert all(item['done'] for item in species.bimol_pose_searches)
        return
    qc.qc_bimol_search.reset_mock()
    bimol.advance_pose_searches(species, par, qc)
    # A rejected common TS leaves each site's next pose independent.
    assert [item['cursor'] for item in species.bimol_pose_searches] == [[0, 0, 1, 0], [0, 0, 2, 0]]
    assert all(item['attempt'] == 1 and item['active'] is None and not item['done']
               for item in species.bimol_pose_searches)
    submitted = {call.args[0].instance_name: call.kwargs['start_after']
                 for call in qc.qc_bimol_search.call_args_list}
    assert submitted == {'path_0': [0, 0, 1, 0], 'path_1': [0, 0, 2, 0]}


@pytest.mark.parametrize('saved_status', ['normal', 'running', 'error'])
def test_completed_or_started_duplicate_l1_is_restored(tmp_path, monkeypatch, saved_status):
    species, par, qc, finish, statuses = equivalent_sites(tmp_path, monkeypatch)
    finish(0)
    bimol.advance_pose_searches(species, par, qc)
    candidate = finish(1)
    entry = dict(reaction=species.reac_obj[1], candidate=candidate)
    digest = bimol._candidate_digest(species, [entry], entry, append=True)
    statuses[f'path_1_l0_ts_{digest}'] = saved_status
    assert not bimol.advance_pose_searches(species, par, qc)
    assert len(species.reac_obj) == 4
    assert [item['active'] for item in species.bimol_pose_searches] == [2, 3]


@pytest.mark.parametrize('saved_status', ['normal', 'running', 'error'])
@pytest.mark.parametrize('attempt', [0, 1])
def test_restart_restores_later_listed_owner_before_equivalent_site(
        tmp_path, monkeypatch, saved_status, attempt):
    species, par, qc, finish, statuses = equivalent_sites(tmp_path, monkeypatch)
    original = {field: list(getattr(species, field)) for field in bimol._REACTION_FIELDS}
    finish(1)
    bimol.advance_pose_searches(species, par, qc)
    first = species.reac_obj[2]
    finish(0)
    bimol.advance_pose_searches(species, par, qc)
    assert len(species.reac_obj) == 3 and first.search_parent == 'path_1'
    statuses[first.instance_name] = saved_status if attempt == 0 else 'error'
    if attempt:
        distinct = species.geom.copy()
        distinct[5:] += [1., 0., 0.]
        candidate = finish(0, attempt=1, pose=1, geom=distinct)
        finish(1, attempt=1, pose=1, geom=distinct)
        # Ownership changes on retry: site 0 has the only new L1 job.
        entry = dict(reaction=species.bimol_pose_searches[0]['reaction'], candidate=candidate)
        digest = bimol._candidate_digest(species, [entry], entry, append=True)
        statuses[f'path_0_l0_ts_{digest}'] = saved_status
    for field, values in original.items():
        setattr(species, field, list(values))
    bimol._initialize_pose_searches(species, par)
    bimol.advance_pose_searches(species, par, qc)
    assert len(species.reac_obj) == 3
    assert species.reac_obj[2].instance_name == first.instance_name
    assert [state['active'] for state in species.bimol_pose_searches] == [2, 2]
    if attempt:
        species.reac_ts_done[2] = -999
        bimol.advance_pose_searches(species, par, qc)
        assert len(species.reac_obj) == 4
        assert species.reac_obj[3].search_parent == 'path_0'
        assert [state['active'] for state in species.bimol_pose_searches] == [3, 3]


def test_different_configured_products_keep_independent_l1(tmp_path, monkeypatch):
    species, par, qc, finish, statuses = equivalent_sites(tmp_path, monkeypatch)
    finish(0)
    candidate = finish(1)
    candidate['validation']['endpoints'][1]['configured_identity']['id'] = 'different_product'
    assert not bimol.advance_pose_searches(species, par, qc)
    assert len(species.reac_obj) == 4


def test_exchanged_ethyl_fragments_keep_both_sites_and_one_degeneracy(tmp_path, monkeypatch):
    atoms = molecule('C2H5')
    fragments = []
    for name in ('ethyl_a', 'ethyl_b'):
        point = StationaryPoint(name, 0, 2, atom=atoms.get_chemical_symbols(), geom=atoms.positions)
        point.characterize()
        fragments.append(point)
    pair = StationaryPoint('ethyl_pair', 0, 3,
        atom=list(fragments[0].atom) + list(fragments[1].atom),
        geom=np.vstack((atoms.positions, atoms.positions + [6., 0., 0.])),
        fragA=fragments[0], fragB=fragments[1])
    pair.characterize()
    monkeypatch.setattr('tests.test_bimol_sequential_poses.reactants', lambda: pair)
    species, par, qc, finish, statuses = equivalent_sites(tmp_path, monkeypatch)
    instances = [[1, 9, 7], [8, 2, 0]]
    for index, instance in enumerate(instances):
        species.reac_obj[index].instance = instance
        species.reac_obj[index].site_degeneracy = 3
        species.reac_inst[index] = instance
        species.bimol_pose_searches[index]['values']['reac_inst'] = instance
    finish(0)
    finish(1, geom=np.vstack((pair.geom[7:], pair.geom[:7])))
    assert bimol.advance_pose_searches(species, par, qc)
    assert len(species.reac_obj) == 3
    representative = species.reac_obj[2]
    assert representative.search_instances == instances
    assert BimolecularNetwork._site_degeneracy(None, [representative]) == 3
    assert [item['active'] for item in species.bimol_pose_searches] == [2, 2]


def test_generated_worker_uses_separated_reactants_and_requested_comparison(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('input.json').write_text(json.dumps(dict(barrier_threshold=100., queuing='local',
        qc='fc', fc_model_path='test-model', bimol_ts_rmsd=.03, bimol_energy_tolerance=.07)))
    par = Parameters('input.json', show_warnings=False).par
    par['bimol'] = 1
    qc = QuantumChemistry(par).l0
    monkeypatch.setattr(qc, 'submit_qc', Mock())
    species = reactants()
    reaction = Abstraction(species, qc, par, [5, 1, 0], 'methane_oh')
    reaction.allow_reflection = False
    qc.qc_bimol_search(reaction)
    monkeypatch.setattr(bimol, 'create_calculator', lambda *args, **kwargs: None)

    def search(worker_species, *args, **kwargs):
        np.testing.assert_array_equal(worker_species.bond01, species.bond01)
        np.testing.assert_array_equal(worker_species.atomid, species.atomid)
        assert not worker_species.bond01[:5, 5:].any()
        settings = args[-1]
        assert settings['bimol_ts_rmsd'] == .03
        assert settings['bimol_energy_tolerance'] == .07
        assert settings['allow_reflection'] is False
        return -1, [], dict(status='error', exhausted=True)

    check = Mock(side_effect=search)
    monkeypatch.setattr(bimol, 'find_ts_candidates', check)
    job = bimol.search_job(reaction.instance_name, sequential=True)
    exec(compile(Path(job + '.py').read_text(), job + '.py', 'exec'), {})
    check.assert_called_once()
    saved = json.loads(Path(job + '_search.json').read_text())
    assert saved == dict(status='error', exhausted=True)
    assert Path(job + '.done').is_file()
