"""Counting oracles, not calculated H + methanol rate coefficients.

Fernandez-Ramos et al., Theor Chem Acc 118, 813 (2007),
doi:10.1007/s00214-007-0328-0, Eq. 10 and section 3.5.
"""
import itertools
import unittest

import numpy as np
from scipy.integrate import quad


class TestMethanolCounting(unittest.TestCase):
    def test_labelled_and_representative_harmonic_counts_agree(self):
        # With one reactant RRHO configuration (sigma_ext = 1), the two
        # distinct mirror saddles give 2*q_TS. Three methyl labels do not
        # turn that into 6*q_TS: the labelled convention normalizes by 3.
        labelled = list(itertools.product(range(3), ('mirror_a', 'mirror_b')))
        for temperature in (200., 500., 1500.):
            one_saddle = np.exp(-7. / (.001987204258 * temperature))
            explicit = sum(one_saddle / 3 for _ in labelled)
            representative = 3 * (2 * one_saddle / 3)
            conventional = 2 * one_saddle
            self.assertAlmostEqual(explicit, representative)
            self.assertAlmostEqual(representative, conventional)

    def test_full_mirror_covering_rotor_needs_no_extra_optical_multiplier(self):
        # Analytic model with two energetically equivalent *distinct* TS
        # wells. Potential periodicity two does not imply sigma_int = 2.
        for temperature in (200., 500., 1500.):
            beta = 1 / (.001987204258 * temperature)
            density = lambda angle: np.exp(-beta * 2 * (1 - np.cos(2 * angle)))
            full = quad(density, 0., 2 * np.pi)[0]  # TS sigma_int = 1
            sectors = sum(quad(density, i * np.pi, (i + 1) * np.pi)[0]
                          for i in range(2))
            representative = 2 * quad(density, 0., np.pi)[0]
            self.assertAlmostEqual(full, sectors)
            self.assertAlmostEqual(full, representative)
            self.assertNotAlmostEqual(full, 2 * sectors)

    def test_reactant_internal_symmetry_removes_labelled_methyl_repetitions(self):
        for temperature in (200., 500., 1500.):
            beta = 1 / (.001987204258 * temperature)
            density = lambda angle: np.exp(-beta * .5 * (1 - np.cos(3 * angle)))
            labelled = quad(density, 0., 2 * np.pi)[0]
            physical = labelled / 3  # methanol C-O sigma_int = 3
            fundamental = quad(density, 0., 2 * np.pi / 3)[0]
            self.assertAlmostEqual(physical, fundamental)
            # Symmetry-only thought experiment, not a physical isotope rate:
            # distinguishing the methyl labels triples the torsional Q_R.
            self.assertAlmostEqual(labelled / physical, 3.)


if __name__ == '__main__':
    unittest.main()
