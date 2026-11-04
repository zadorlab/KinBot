import copy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from kinbot.motion_character import characterize_modes, tangent_subspaces
from kinbot.stationary_pt import StationaryPoint


DATA = json.loads((Path(__file__).parent/'reference/soft_motion_native.json').read_text())


def native(label):
    record = DATA[label]
    point = StationaryPoint(label, record['charge'], record['multiplicity'],
        atom=np.asarray(record['atoms']), geom=np.asarray(record['geometry']))
    point.characterize()
    return point, np.asarray(record['modes'])


@pytest.mark.parametrize('label,lo,hi', [('rxn22', .98, .995), ('rxn29', .96, .98)])
def test_native_saddle_torsional_character_is_not_a_global_path_verdict(label, lo, hi):
    point, modes = native(label)
    result = characterize_modes(point, modes)
    assert lo < result['modes'][0]['torsional_fraction'] < hi
    assert not result['lowest_path_established']
    assert not result['alternative_paths_excluded']
    assert not result['anharmonic_partition_function_established']
    assert not result['population_merge_authorized']


@pytest.mark.parametrize('label', ['w7', 'w14', 'w15'])
def test_soft_non_torsional_native_mode_is_not_lost_when_axes_exist(label):
    point, modes = native(label)
    result = characterize_modes(point, modes)
    assert result['torsional_rank'] == 2
    assert result['modes'][0]['torsional_fraction'] > .97
    assert result['modes'][1]['residual_fraction'] > .97
    assert DATA[label]['frequencies_cm1'][1] < 300.


def test_w3_comparison_case_has_mixed_soft_motion_without_implied_extra_well():
    point, modes = native('w3')
    result = characterize_modes(point, modes)
    assert result['torsional_rank'] == 2
    assert .90 < result['modes'][0]['torsional_fraction'] < .93
    assert .23 < result['modes'][1]['torsional_fraction'] < .27
    assert not result['alternative_paths_excluded']
    assert not result['population_merge_authorized']


def test_projection_invariant_under_frame_atom_order_sign_and_scale():
    point, modes = native('rxn22')
    result = characterize_modes(point, modes)
    rng = np.random.default_rng(20)
    for _ in range(6):
        rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        order = rng.permutation(len(point.atom))
        inverse = np.argsort(order)
        view = copy.copy(point)
        view.atom = point.atom[order]
        view.geom = (point.geom @ rotation + [3., -4., .5])[order]
        view.bond = point.bond[np.ix_(order, order)]
        view.dihed = [inverse[rotor].tolist() for rotor in point.dihed[::-1]]
        observed = characterize_modes(view, (modes @ rotation)[:, order] * -7.)
        np.testing.assert_allclose([v['torsional_fraction'] for v in observed['modes']],
                                   [v['torsional_fraction'] for v in result['modes']], atol=1.e-12)


def test_pure_bend_with_no_axis_stays_non_torsional():
    point = SimpleNamespace(atom=['H','O','H'], geom=np.array([[1.,0.,0.],[0.,0.,0.],[0.,1.,0.]]),
        bond=np.array([[0,1,0],[1,0,1],[0,1,0]]), dihed=[])
    result = characterize_modes(point, [[[0.,1.,0.],[0.,0.,0.],[1.,0.,0.]]])
    assert result['torsional_rank'] == 0
    assert result['modes'][0]['residual_fraction'] == 1.


def test_bad_axes_and_unbound_masses_are_not_silently_interpreted():
    point, modes = native('w7')
    with pytest.raises(ValueError, match='Duplicate'):
        tangent_subspaces(point.geom, [12,32,32,1,1,1], point.bond, [[0,1],[1,0]])
    with pytest.raises(ValueError, match='splits'):
        tangent_subspaces(point.geom, [12,32,32,1,1,1], point.bond, [[1,2]])
    point.isotopes = [13,0,0,0,0,0]
    with pytest.raises(ValueError, match='isotope'):
        characterize_modes(point, modes)
    assert characterize_modes(point, modes, masses=[13,32,32,1,1,1])['masses_amu'][0] == 13.


def test_axis_without_off_axis_top_is_not_an_extra_mode():
    ext, tor, individual, axes = tangent_subspaces([[0,0,0],[1,0,0]], [1,1], [[0,1],[1,0]], [[0,1]])
    assert len(ext) == 5 and len(tor) == 0
    assert not axes[0]['nonzero']
