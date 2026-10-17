import warnings

import numpy as np
import pytest

from kinbot import geometry
from kinbot.modify_geom import modify_coordinates
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.reactions.reac_h2_elim import H2Elim
from kinbot.reactions.reac_r12_insertion_R import R12Insertion
from kinbot.stationary_pt import StationaryPoint


@pytest.fixture(scope='module')
def acetaldehyde_path_context(tmp_path_factory):
    import json
    path = tmp_path_factory.mktemp('path_input')/'input.json'
    path.write_text(json.dumps({'smiles': 'CC=O', 'charge': 0, 'mult': 1,
        'barrier_threshold': 100., 'theory_preset': 'uma-b2plyp-anl',
        'composite_method': 'ANL0-F12', 'queue_job_limit': 8}))
    parameters = Parameters(str(path))
    qc = QuantumChemistry(parameters.par)
    species = StationaryPoint('acetaldehyde', 0, 1, smiles='CC=O')
    species.characterize()
    return parameters.par, qc, species


def _zero_based(changes):
    converted = []
    for change in changes:
        if change[0] == 'L':
            item = ['L'] + [index - 1 for index in change[1:-1]]
        else:
            item = [index - 1 for index in change[:-1]]
        item.append(change[-1])
        converted.append(item)
    return converted


@pytest.mark.parametrize(
    'reaction_type,instance,name',
    [(H2Elim, [3, 0, 1, 6], 'acetaldehyde_h2_elim_4_7'),
     (R12Insertion, [6, 1, 0],
      'acetaldehyde_r12_insertion_R_7_2_1')])
def test_selected_acetaldehyde_path_has_finite_verified_coordinate_guesses(
        acetaldehyde_path_context, reaction_type, instance, name):
    parameters, qc, species = acetaldehyde_path_context
    reaction = reaction_type(species, qc, parameters, instance, name)
    current = species.geom.copy()

    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        for requested_step in range(reaction.max_step):
            _, _, changes, _ = reaction.get_constraints(
                requested_step, current)
            success, modified = modify_coordinates(
                species, name, current, _zero_based(changes), species.bond)
            assert success == 1, requested_step
            assert modified.shape == current.shape
            assert np.isfinite(modified).all()
            assert not np.allclose(modified, current)
            for change in _zero_based(changes):
                if len(change) == 3:
                    actual = np.linalg.norm(
                        modified[change[0]] - modified[change[1]])
                    assert actual == pytest.approx(change[2], abs=0.05)
                elif len(change) == 4:
                    actual = geometry.calc_angle(
                        modified[change[0]], modified[change[1]],
                        modified[change[2]])
                    assert actual == pytest.approx(
                        np.radians(change[3]), abs=0.05)
                elif len(change) == 5:
                    actual = geometry.calc_dihedral(
                        modified[change[0]], modified[change[1]],
                        modified[change[2]], modified[change[3]])[0]
                    difference = (actual - change[4] + 180.) % 360. - 180.
                    assert abs(difference) < 0.05 * 180. / np.pi
            current = modified
