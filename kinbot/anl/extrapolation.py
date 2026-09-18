"""Two-point CBS arithmetic from Klippenstein, Harding, and Ruscic (2017).

The original ANL paper extrapolates *total* CCSD(T) energies with an
inverse-angular-momentum power of 3.7. Later ANL-family variants may use
different fitted coefficients and must request those explicitly.
"""

from __future__ import annotations

import math


ORIGINAL_ANL_POWER = 3.7


def extrapolation_coefficient(upper_cardinal, *, power=ORIGINAL_ANL_POWER):
    """Return alpha in E_CBS = E_upper + alpha*(E_upper - E_lower)."""
    if isinstance(upper_cardinal, bool) or not isinstance(upper_cardinal, int) \
            or upper_cardinal < 3:
        raise ValueError('Upper cardinal must be an integer of at least 3.')
    if isinstance(power, bool) or not isinstance(power, (int, float)) \
            or not math.isfinite(power) or power <= 0:
        raise ValueError('CBS power must be positive and finite.')
    lower = upper_cardinal - 1
    return lower ** power / (upper_cardinal ** power - lower ** power)


def two_point_cbs(lower_hartree, upper_hartree, *, upper_cardinal,
                  power=ORIGINAL_ANL_POWER):
    """Extrapolate adjacent cardinal-basis values using the paper's formula."""
    for value in (lower_hartree, upper_hartree):
        if isinstance(value, bool) or not isinstance(value, (int, float)) \
                or not math.isfinite(value):
            raise ValueError('CBS inputs must be finite Hartree values.')
    alpha = extrapolation_coefficient(upper_cardinal, power=power)
    result = math.fsum((upper_hartree,
                        alpha * (upper_hartree - lower_hartree)))
    if not math.isfinite(result):
        raise ValueError('CBS result is nonfinite.')
    return result


def core_valence_correction(*, all_electron_lower, all_electron_upper,
                            frozen_core_lower, frozen_core_upper,
                            upper_cardinal=4, power=ORIGINAL_ANL_POWER):
    """Subtract frozen-core from all-electron CCSD(T) CBS energies."""
    full = two_point_cbs(all_electron_lower, all_electron_upper,
                         upper_cardinal=upper_cardinal, power=power)
    frozen = two_point_cbs(frozen_core_lower, frozen_core_upper,
                           upper_cardinal=upper_cardinal, power=power)
    return math.fsum((full, -frozen))
