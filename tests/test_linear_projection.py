"""Near-linear molecules must keep 3N-5 vibrations after external-mode projection,
and elongated but non-linear molecules must keep 3N-6."""

import unittest

from ase import Atoms
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


def co2_with_angle(degrees):
    """CO2 with the given O-C-O angle, carbon at the origin."""
    atoms = molecule('CO2')
    half = np.radians((180. - degrees) / 2.)
    atoms.positions[1] = [1.16 * np.cos(half), 1.16 * np.sin(half), 0.]
    atoms.positions[2] = [-1.16 * np.cos(half), 1.16 * np.sin(half), 0.]
    return atoms


def methyl_cyanodiyne():
    """CH3-C#C-C#C-C#N: a long, prolate but non-linear molecule."""
    x = [0.0, 1.46, 2.67, 4.05, 5.26, 6.64, 7.80]
    symbols = ['C', 'C', 'C', 'C', 'C', 'C', 'N']
    positions = [[xi, 0., 0.] for xi in x]
    for k in range(3):
        angle = 2 * np.pi * k / 3
        positions.append([-0.36, 1.03 * np.cos(angle), 1.03 * np.sin(angle)])
        symbols.append('H')
    return Atoms(symbols, positions=positions)


class TestLinearProjection(unittest.TestCase):
    def test_exactly_linear_molecules_keep_3n_minus_5(self):
        for name in ('CO2', 'HCN', 'C2H2', 'N2'):
            with self.subTest(molecule=name):
                atoms = molecule(name)
                count, _ = projected_count(atoms)
                self.assertEqual(count, 3 * len(atoms) - 5)

    def test_optimiser_residual_bend_still_counts_as_linear(self):
        for degrees in (179.99, 179.5, 178., 175.):
            with self.subTest(oco_angle=degrees):
                count, _ = projected_count(co2_with_angle(degrees))
                self.assertEqual(count, 4)

    def test_clearly_bent_co2_is_not_linear(self):
        for degrees in (160., 150., 120.):
            with self.subTest(oco_angle=degrees):
                count, _ = projected_count(co2_with_angle(degrees))
                self.assertEqual(count, 3)

    def test_prolate_molecule_with_off_axis_hydrogens_is_not_linear(self):
        # Its smallest moment of inertia is under 1% of the largest, so a
        # moment-ratio test calls it linear; the off-axis methyl hydrogens say no.
        atoms = methyl_cyanodiyne()
        count, _ = projected_count(atoms)
        self.assertEqual(count, 3 * len(atoms) - 6)

    def test_genuinely_bent_molecules_keep_3n_minus_6(self):
        for name in ('H2O', 'NH3', 'CH3OH', 'C2H6'):
            with self.subTest(molecule=name):
                atoms = molecule(name)
                count, _ = projected_count(atoms)
                self.assertEqual(count, 3 * len(atoms) - 6)

    def test_frequencies_are_real_for_degenerate_spectra(self):
        atoms = molecule('CO2')
        rng = np.random.default_rng(7)
        block = rng.normal(size=(9, 9))
        hessian = block @ block.T
        count, freqs = projected_count(atoms, hessian)
        self.assertEqual(count, 4)
        self.assertTrue(all(isinstance(f, float) for f in freqs), freqs)


if __name__ == '__main__':
    unittest.main()
