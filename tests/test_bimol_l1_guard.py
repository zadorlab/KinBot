"""Reject changed L1 saddles before starting their conformer calculations."""
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from kinbot import bimol, reaction_generator
from kinbot.optimize import same_ts_geometry, valid_ts_frequencies
from kinbot.stereo_identity import configured_geometry_allowed
from tests.test_bimol_stereo_integration import donor_and_hydrogen, abstraction
from tests.test_bimol_sequential_poses import controller


def promoted_reaction(smiles='CO'):
    reaction = abstraction(donor_and_hydrogen(smiles))
    reaction.search_candidate = {'geometry': reaction.species.geom.tolist()}
    bimol.build_bimolecular_ts(reaction, reaction.qc)
    return reaction


def test_unchanged_l1_ts_passes_without_a_saved_l0_mode():
    reaction = promoted_reaction()
    bimol.check_l1_ts_refinement(reaction, {'imagfreq_threshold': 50.})


def test_forming_bond_is_checked_against_l0_instead_of_l1():
    reaction = promoted_reaction()
    acceptor, hydrogen, _ = reaction.instance
    # The incoming H is unbonded in the reactant. Its forming H-H bond must
    # nevertheless be part of the saddle comparison.
    assert reaction.species.bond[acceptor, hydrogen] == 0
    assert reaction.ts.bond[acceptor, hydrogen] > 0
    reaction.ts.geom[acceptor] += [5., 0., 0.]
    with pytest.raises(ValueError, match='bond lengths'):
        bimol.check_l1_ts_refinement(reaction, {'imagfreq_threshold': 50.})


@pytest.mark.parametrize('freq, accepted', [
    ([-1000., 100., 200.], True), ([-1000., -20., 200.], True),
    ([-1000., -50., 200.], True), ([-1000., -51., 200.], False),
    ([-1000., -20., -10.], False), ([100., 200., 300.], False)])
def test_l1_frequencies_keep_the_ordinary_ts_tolerance(freq, accepted):
    reaction = promoted_reaction()
    reaction.ts.freq = freq
    assert valid_ts_frequencies(freq, 50.) == accepted
    if accepted:
        bimol.check_l1_ts_refinement(reaction, {'imagfreq_threshold': 50.})
    else:
        with pytest.raises(ValueError, match='imaginary frequencies'):
            bimol.check_l1_ts_refinement(reaction, {'imagfreq_threshold': 50.})


def test_l1_stereochemistry_is_checked_against_l0_not_itself():
    reaction = promoted_reaction('CC[C@H](O)C')
    l0_geometry = reaction.species.geom.copy()
    reflected = l0_geometry * [-1., 1., 1.]
    reaction.qc.get_qc_geom = lambda *args, **kwargs: (0, reflected)
    row = next(reaction.qc.db.select())
    row.positions = reflected
    # Build exactly as the production L1 handoff does, with L1 stereo metadata.
    bimol.build_bimolecular_ts(reaction, reaction.qc)
    assert configured_geometry_allowed(reaction.ts, reaction.ts.geom)
    with pytest.raises(ValueError, match='configured stereochemistry'):
        bimol.check_l1_ts_refinement(reaction, {'imagfreq_threshold': 50.})
    # The ordinary racemic setting permits the global mirror.
    bimol.check_l1_ts_refinement(reaction,
        {'imagfreq_threshold': 50., 'optical_population': 'racemic'})
    np.testing.assert_array_equal(reaction.ts.geom, reflected)


def test_ordinary_mode_assisted_geometry_arm_is_preserved():
    reference = SimpleNamespace(natom=3, wellorts=1,
        geom=np.array([[-1., 0., 0.], [0., 0., 0.], [1., 0., 0.]]),
        bond=np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]))
    refined = SimpleNamespace(geom=reference.geom * 1.25, bond=reference.bond)
    mode = np.array([[0., 1., 0.], [1., 0., 0.], [0., 0., 1.]])
    assert same_ts_geometry(reference, refined, mode, mode)
    assert not same_ts_geometry(reference, refined)
    refined.geom = reference.geom * 1.15
    assert same_ts_geometry(reference, refined)
    refined.bond = np.ones((3, 3)) - np.eye(3)
    assert not same_ts_geometry(reference, refined, mode, mode)


def test_rejected_l1_ts_skips_conformers_and_reopens_remaining_poses(tmp_path, monkeypatch):
    species, par, qc, finish, statuses = controller(tmp_path, monkeypatch, count=1)
    par.update(bimol=1, imagfreq_threshold=50.)
    finish(0)
    bimol.advance_pose_searches(species, par, qc)
    reaction = species.reac_obj[1]
    invalid = promoted_reaction()
    invalid.ts.freq = [100.] * len(invalid.ts.freq)
    reaction.search_candidate['geometry'] = invalid.search_candidate['geometry']
    monkeypatch.setattr(bimol, 'build_bimolecular_ts',
                        lambda target, qc: setattr(target, 'ts', invalid.ts))
    optimize = Mock(side_effect=AssertionError('Rejected TS must not start conformer work'))
    monkeypatch.setattr(reaction_generator, 'Optimize', optimize)
    generator = reaction_generator.ReactionGenerator(species, par, qc, 'input.json')
    generator._advance_bimol_ts_selection(1, reaction)
    assert species.reac_ts_done[1] == -999
    assert reaction.bimol_final_status == 'l1_refinement_rejected'
    optimize.assert_not_called()
    assert bimol.advance_pose_searches(species, par, qc)
    assert qc.qc_bimol_search.call_args.kwargs == {
        'attempt': 1, 'start_after': [0, 0, 0, 0]}


@pytest.mark.parametrize('level', ['l1', 'l2'])
def test_final_irc_diagnostic_uses_its_assigned_profile(level, monkeypatch, caplog):
    irc = SimpleNamespace(qc=SimpleNamespace(theory_level=level),
        check_irc=lambda: [0, 0], do_irc_calculations=Mock())
    monkeypatch.setattr(reaction_generator, 'IRC', lambda *args: irc)
    reaction = SimpleNamespace(instance_name='hydrogen_abstraction')
    generator = reaction_generator.ReactionGenerator(None, {'bimol': 1}, None, None)
    with caplog.at_level('INFO', logger='KinBot'):
        generator._advance_irc(0, reaction, 12.3)
    assert f'at {level.upper()} for hydrogen_abstraction' in caplog.text
