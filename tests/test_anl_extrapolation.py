"""Original ANL two-point coefficients and correction arithmetic."""

import math

import pytest

from kinbot.anl.extrapolation import (
    ORIGINAL_ANL_POWER, core_valence_correction,
    extrapolation_coefficient, two_point_cbs,
)
from kinbot.anl.recipes import recipe


@pytest.mark.parametrize(('upper', 'reported'), [
    (3, 0.29), (4, 0.53), (5, 0.78), (6, 1.04),
])
def test_original_paper_coefficients_are_derived_not_rounded(upper, reported):
    coefficient = extrapolation_coefficient(upper)
    assert coefficient == pytest.approx(
        (upper - 1) ** ORIGINAL_ANL_POWER /
        (upper ** ORIGINAL_ANL_POWER - (upper - 1) ** ORIGINAL_ANL_POWER))
    assert round(coefficient, 2) == reported


def test_ch4_native_f12_totals_give_reproducible_tq_cbs_value():
    # Exact F12b totals from the licensed first CH4 dispatch, not F12a.
    result = two_point_cbs(-40.454906199189, -40.456608306474,
                           upper_cardinal=4)
    assert result == pytest.approx(-40.457504545051066, abs=1e-12)


def test_core_valence_is_all_electron_cbs_minus_frozen_core_cbs():
    correction = core_valence_correction(
        all_electron_lower=-100.20, all_electron_upper=-100.25,
        frozen_core_lower=-100.17, frozen_core_upper=-100.21)
    alpha = extrapolation_coefficient(4)
    assert correction == pytest.approx(-0.04 - alpha * 0.01)


def test_recipe_pins_each_reference_pair_and_the_tq_corrections():
    for name, cardinal in (('ANL0', 5), ('ANL0-F12', 4), ('ANL1', 6)):
        requirements = {item.key: item for item in recipe(name).requirements}
        assert requirements['reference_cbs'].settings['upper_cardinal'] == cardinal
        assert requirements['reference_cbs'].settings['extrapolation_power'] == 3.7
        assert requirements['core_valence_cbs'].settings['upper_cardinal'] == 4
        assert requirements['core_valence_cbs'].method == \
            'CCSD(T,full)-CCSD(T,frozen-core)'
        assert requirements['scalar_relativistic'].basis == 'aug-cc-pcVTZ-DK'
    assert {item.key: item for item in recipe('ANL1').requirements}[
        'harmonic_zpe'].settings['upper_cardinal'] == 4


@pytest.mark.parametrize('bad', [math.nan, math.inf, True, '1.0'])
def test_extrapolation_rejects_nonphysical_inputs(bad):
    with pytest.raises(ValueError):
        two_point_cbs(bad, -1.0, upper_cardinal=4)
    with pytest.raises(ValueError):
        two_point_cbs(-1.0, bad, upper_cardinal=4)
