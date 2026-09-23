"""Endpoint interpretation distinguishes stereo without running extra IRCs."""
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from kinbot.irc import IRC
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity


def point(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=312) == 0
    species = StationaryPoint('endpoint', 0, 1,
        atom=[a.GetSymbol() for a in mol.GetAtoms()],
        geom=mol.GetConformer().GetPositions())
    species.characterize()
    return species


class TestIRCStereoIdentity(unittest.TestCase):
    def interpret(self, reactant, other, population='specified'):
        qc = SimpleNamespace(get_qc_geom=lambda name, *a, **kw:
            (0, reactant.geom if '_F_' in name else other.geom))
        reaction = SimpleNamespace(species=reactant, instance_name='test', qc=qc)
        return IRC(reaction, {'bimol': 0, 'optical_population': population}).irc2stationary_pt()

    def test_ez_reaction_is_not_discarded_as_two_initial_endpoints(self):
        trans, cis = point('C/C=C/C'), point('C/C=C\\C')
        self.assertEqual(trans.chemid, cis.chemid)
        np.testing.assert_array_equal(trans.chiral, cis.chiral)
        product = self.interpret(trans, cis)
        self.assertIsInstance(product, StationaryPoint)
        self.assertEqual(canonical_identity(product)['id'], canonical_identity(cis)['id'])

    def test_mirror_scope_does_not_include_diastereomers(self):
        first = point('C[C@H](F)[C@H](Cl)C')
        mirror = point('C[C@@H](F)[C@@H](Cl)C')
        other = point('C[C@H](F)[C@@H](Cl)C')
        self.assertIsInstance(self.interpret(first, mirror), StationaryPoint)
        self.assertEqual(self.interpret(first, mirror, 'racemic'), 0)
        self.assertIsInstance(self.interpret(first, other, 'racemic'), StationaryPoint)

    def test_ordinary_identical_endpoints_and_unknown_identity_keep_old_behavior(self):
        ethanol = point('CCO')
        self.assertEqual(self.interpret(ethanol, ethanol), 0)
        first, second = point('C[C@H](F)Cl'), point('C[C@@H](F)Cl')
        with patch('kinbot.irc.canonical_identity', return_value={'status': 'unavailable'}):
            self.assertEqual(self.interpret(ethanol, ethanol), 0)
            self.assertIsInstance(self.interpret(first, second), StationaryPoint)


if __name__ == '__main__':
    unittest.main()
