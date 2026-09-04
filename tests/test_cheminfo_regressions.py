"""Regression tests for optional chemistry helpers and reaction images."""

import importlib.util
import unittest

import numpy as np

from kinbot import cheminfo


@unittest.skipUnless(importlib.util.find_spec('rdkit'), 'RDKit is optional')
class TestRDKitHelpers(unittest.TestCase):
    def test_molecular_formula_loads_rdkit_locally(self):
        self.assertEqual(cheminfo.get_molecular_formula('C=O'), 'CH2O')

    def test_rdkit_structure_contains_coordinates_and_bonds(self):
        molecule, structure, bond = cheminfo.generate_3d_structure('CO', obabel=0)
        self.assertEqual(molecule.GetNumAtoms(), 6)
        self.assertEqual(len(structure), 24)
        self.assertEqual(sorted(structure[::4]), ['C', 'H', 'H', 'H', 'H', 'O'])
        coordinates = np.asarray([structure[i + 1:i + 4] for i in range(0, 24, 4)])
        self.assertTrue(np.isfinite(coordinates).all())
        self.assertGreater(np.linalg.norm(coordinates[0] - coordinates[1]), 0.5)
        np.testing.assert_array_equal(bond, bond.T)
        self.assertEqual(np.count_nonzero(bond), 10)


if __name__ == '__main__':
    unittest.main()
