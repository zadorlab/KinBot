"""Counterexamples to independently fitting mode vectors; no QC required."""
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np
import pytest

from kinbot import constants, geometry
from kinbot.spatial import proper_alignment, mode_overlap, aligned_mode_overlap
from kinbot.qc import QuantumChemistry
from kinbot.optimize import Optimize
from kinbot.stationary_pt import StationaryPoint


GEOM = np.array([[0., 0., 0.], [.9, 0., 0.], [-.2, .9, 0.]])
MODE = np.array([[.2, -.1, 0.], [-.3, .6, 0.], [.1, -.5, 0.]])
ROT = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
MASS = [constants.exact_mass[s] for s in ['O', 'H', 'H']]


def test_geometry_frame_rotation_and_vector_sign_do_not_change_mode_identity():
    for sign in (-1, 1):
        assert aligned_mode_overlap(GEOM, GEOM @ ROT + 30., MODE,
            sign * MODE @ ROT, masses=MASS) == pytest.approx(1.)
    order = [2, 0, 1]
    assert aligned_mode_overlap(GEOM[order], (GEOM @ ROT + 30.)[order],
        MODE[order], (MODE @ ROT)[order], masses=np.array(MASS)[order]) == pytest.approx(1.)


def test_different_motion_cannot_be_fitted_independently_of_its_geometry():
    # Old routine can rotate this 90-degree-different displacement back onto
    # the first. The molecular frame is unchanged, so its true dot product is 0.
    other = MODE @ ROT
    assert geometry.matrix_corr(MODE, other) == pytest.approx(1.)
    assert aligned_mode_overlap(GEOM, GEOM, MODE, other, masses=MASS) == pytest.approx(0., abs=1.e-14)


def test_vector_translation_is_not_a_coordinate_translation():
    assert mode_overlap(MODE, MODE + 2., masses=MASS) < .9
    assert mode_overlap(MODE, -MODE, masses=MASS) == pytest.approx(1.)
    assert mode_overlap(MODE, -MODE, masses=MASS, signed=True) == pytest.approx(-1.)
    for bad in (np.zeros_like(MODE), np.full_like(MODE, np.nan)):
        with pytest.raises(ValueError):
            mode_overlap(MODE, bad, masses=MASS)


def test_proper_fit_cannot_reflect_a_chiral_frame():
    tetra = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
    fit = proper_alignment(tetra, tetra * [-1, 1, 1])
    assert np.linalg.det(fit.rotation) == pytest.approx(1.)
    assert fit.rmsd > .1


def test_ase_native_modes_are_rotated_with_their_recorded_frame(tmp_path):
    job = tmp_path/'native'
    Path(str(job)+'_vib.xyz').write_text('3\nnative vibration\n'+''.join(
        f'{s} '+ ' '.join(map(str, [*g, *v]))+'\n'
        for s, g, v in zip(['O', 'H', 'H'], GEOM, MODE)))
    qc = SimpleNamespace(qc='fc', use_sella=True, frequency_mode='finite_difference')
    result = QuantumChemistry.read_imag_mode(qc, job, 3,
        geometry=GEOM @ ROT + 4., symbols=['O', 'H', 'H'])
    np.testing.assert_allclose(result, MODE @ ROT, atol=1.e-12)
    with pytest.raises(ValueError, match='differs'):
        QuantumChemistry.read_imag_mode(qc, job, 3,
            geometry=GEOM * 1.1, symbols=['O', 'H', 'H'])


@pytest.mark.parametrize('motion,close,accepted', [
    ('same', False, True), ('different', False, False),
    ('missing', False, False), ('missing', True, True)])
def test_real_l1_l2_acceptance_uses_bound_overlap_and_keeps_geometry_fallback(motion, close, accepted):
    p = StationaryPoint('ts', 0, 1, atom=['O', 'H', 'H'], geom=GEOM.copy(), wellorts=1)
    p.characterize()
    opt = object.__new__(Optimize)
    opt.species, opt.name, opt.wait = p, p.name, 0
    opt.par = dict(conformer_search=0, multi_conf_tst=0, imagfreq_threshold=50.)
    opt.log_name = lambda high, **kw: 'higher' if high else 'lower'
    other = MODE @ ROT if motion == 'same' else MODE @ ROT @ ROT
    modes = Mock(side_effect=[MODE, other] if motion != 'missing' else OSError('no artifact'))
    opt.qc = SimpleNamespace(read_imag_mode=modes,
        get_qc_geom=lambda job, *a, **kw: (0, GEOM.copy() if job == 'lower' else GEOM @ ROT + 2.),
        get_qc_freq=lambda *a: (0, [-1000., 400., 800.]))
    with patch('kinbot.optimize.geometry.equal_geom', side_effect=lambda p, q, t: close or t == .3), \
         patch('kinbot.optimize.configured_geometry_allowed', return_value=True), \
         patch('kinbot.optimize.load_calculation_record') as load:
        opt.compare_structures()
    assert opt.shigh == (1 if accepted else -999)
    assert load.called is accepted
    assert opt.mode_comparison['status'] == ('unavailable' if motion == 'missing' else 'calculated')
    # No QC submission/retry exists on this verification path.
    np.testing.assert_allclose(modes.call_args_list[0].kwargs['geometry'], GEOM)
