"""Regression tests for optional chemistry helpers and reaction images."""

import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from PIL import Image

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


class TestReactionDepiction(unittest.TestCase):
    def test_reactant_arrow_and_product_are_pasted_in_order(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'tpl').mkdir()
            with Image.new('RGB', (3, 5), 'green') as arrow:
                arrow.save(root / 'tpl/arrow.png')

            def readstring(kind, smiles):
                size, color = ((2, 3), 'red') if smiles == 'C' else ((4, 7), 'blue')

                def draw(show, filename):
                    with Image.new('RGB', size, color) as panel:
                        panel.save(filename)

                return SimpleNamespace(draw=draw)

            with patch.object(cheminfo, 'kb_path', directory), \
                    patch.object(cheminfo, 'pybel', SimpleNamespace(readstring=readstring), create=True):
                cheminfo.create_rxn_depiction('C', 'CO', directory, 'reaction')

            with Image.open(root / 'reaction.png') as result:
                self.assertEqual(result.size, (9, 7))
                self.assertEqual(result.getpixel((0, 2)), (255, 0, 0))
                self.assertEqual(result.getpixel((2, 1)), (0, 128, 0))
                self.assertEqual(result.getpixel((5, 0)), (0, 0, 255))
                self.assertEqual(result.getpixel((0, 0)), (255, 255, 255))


if __name__ == '__main__':
    unittest.main()
