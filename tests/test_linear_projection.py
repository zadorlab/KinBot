"""Near-linear molecules must keep 3N-5 vibrations after external-mode projection."""

import unittest

from ase.build import molecule
import numpy as np

from kinbot import frequencies
from kinbot.stationary_pt import StationaryPoint


def projected_count(atoms, hessian=None):
    """Number of frequencies get_frequencies keeps; independent of the Hessian values."""
    natom = len(atoms)
    if hessian is None:
        hessian = np.eye(3 * natom)
    point = StationaryPoint('probe', 0, 1, atom=atoms.get_chemical_symbols(),
                            geom=atoms.positions)
    point.characterize()
    freqs, _ = frequencies.get_frequencies(point, hessian, atoms.positions)
    return len(freqs), freqs


def bent_co2(degrees):
    """CO2 with both oxygens displaced sideways, as an optimiser might leave it."""
    atoms = molecule('CO2')
    shift = np.array([0., 1.16 * np.sin(np.radians(degrees)), 0.])
    atoms.positions[1] += shift
    atoms.positions[2] += shift
    return atoms


class TestLinearProjection(unittest.TestCase):
    def test_exactly_linear_molecules_keep_3n_minus_5(self):
        for name in ('CO2', 'HCN', 'C2H2', 'N2'):
            with self.subTest(molecule=name):
                atoms = molecule(name)
                count, _ = projected_count(atoms)
                self.assertEqual(count, 3 * len(atoms) - 5)

    def test_optimiser_residual_bend_still_counts_as_linear(self):
        # 0.01 degrees is far below any optimiser's convergence; 2 degrees is a
        # generous residual. Both must still be treated as linear.
        for degrees in (0.01, 0.1, 0.5, 2.0):
            with self.subTest(bend=degrees):
                count, _ = projected_count(bent_co2(degrees))
                self.assertEqual(count, 4)

    def test_genuinely_bent_molecules_keep_3n_minus_6(self):
        for name in ('H2O', 'NH3', 'CH3OH', 'C2H6'):
            with self.subTest(molecule=name):
                atoms = molecule(name)
                count, _ = projected_count(atoms)
                self.assertEqual(count, 3 * len(atoms) - 6)

    def test_frequencies_are_real_for_degenerate_spectra(self):
        # A Hessian with exactly degenerate eigenvalues is where a general
        # eigensolver may return complex pairs; the symmetric solver must not.
        atoms = molecule('CO2')
        rng = np.random.default_rng(7)
        block = rng.normal(size=(9, 9))
        hessian = block @ block.T          # symmetric positive semidefinite
        count, freqs = projected_count(atoms, hessian)
        self.assertEqual(count, 4)
        self.assertTrue(all(isinstance(f, float) for f in freqs), freqs)


if __name__ == '__main__':
    unittest.main()
