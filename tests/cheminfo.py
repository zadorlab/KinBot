###################################################
##                                               ##
## This file is part of the KinBot code v2.0     ##
##                                               ##
## The contents are covered by the terms of the  ##
## BSD 3-clause license included in the LICENSE  ##
## file, found at the root.                      ##
##                                               ##
## Copyright 2018 National Technology &          ##
## Engineering Solutions of Sandia, LLC (NTESS). ##
## Under the terms of Contract DE-NA0003525 with ##
## NTESS, the U.S. Government retains certain    ##
## rights to this software.                      ##
##                                               ##
## Authors:                                      ##
##   Judit Zador                                 ##
##   Ruben Van de Vijver                         ##
##                                               ##
###################################################
"""
This class tests the resonance structure generation algorithm

The data is a dictionary of which the keys are the smiles
and the values are the expected number of resonance isomers
"""
import os
import unittest
import numpy as np

import kinbot.cheminfo as cheminfo


class TestCheminfo(unittest.TestCase):
    def setUp(self):
        pass

    def testMolFrom(self):
        """
        Test the molecular formula creator, which is
        basically a RDKit feature, but used in KinBot
        for postprocessing
        """
        data = {'C1=CC=CC=C1': 'C6H6',
                'C1CCCCC1': 'C6H12',
                'CCCO[O]': 'C3H7O2',
                'C=C': 'C2H4',
                'S=S': 'S2',
                }

        for smi in data:
            cal = cheminfo.get_molecular_formula(smi)
            exp = data[smi]
            warn = smi + ': expected: {}, calculated: {}'.format(exp, cal)
            self.assertEqual(exp, cal, warn)

    def testReactionDepiction(self):
        """
        Test the depiction of a reaction
        """
        reacts = 'CCC[CH2]'
        prods = 'CCC=C.[H]'
        dir = 'output'
        name = 'cheminfo_rxn_depiction'
        path = '{}/{}.png'.format(dir, name)

        # make the directory if not yet present
        if not os.path.exists(dir):
            os.mkdir(dir)
        # remove potential old png files
        if os.path.exists(path):
            os.remove(path)

        # make a new png file
        cheminfo.create_rxn_depiction(reacts, prods, dir, name)

        # check for the existance of the new png file
        self.assertTrue(os.path.exists(path), 'Reaction depiction failed')

    def test3DStructureGeneration(self):
        """
        Test the generation of a 3D structure
        """
        data = {'CCC': 11, '[CH3]': 4, 'C1CCC1': 12}
        for smi in data:
            # do a generation with open babel
            obmol, struct, bond = cheminfo.generate_3d_structure(smi)
            warn = 'Wrong 3D generation with OpenBabel for {}'.format(smi)
            self.assertEqual(len(struct), 4 * data[smi], warn)
            self.assertEqual(len(bond), data[smi], warn)

            # do a generation with RDKit
            rdmol, struct, bond = cheminfo.generate_3d_structure(smi, obabel=0)
            warn = 'Wrong 3D generation with RDKit for {}'.format(smi)
            self.assertEqual(len(struct), 4 * data[smi], warn)
            self.assertEqual(len(bond), data[smi], warn)

    def testCreateOBMol(self):
        """
        Test the creation of an OpenBabel molecule from its smiles
        """
        smis = ['C', 'CCOO', 'CCCO[O]', 'C1CCCCC1']
        for smi in smis:
            obmol = cheminfo.create_ob_mol(smi)
            warn = 'Unable to make OBabel molecule from smiles {}'.format(smi)
            self.assertIsNotNone(obmol, warn)

    def testCreateRDKitMol(self):
        """
        Test the creation of an RDKit molecule from its
        bond matrix and atom vector
        """
        bond = [[0, 1, 1, 1, 1],
                [1, 0, 0, 0, 0],
                [1, 0, 0, 0, 0],
                [1, 0, 0, 0, 0],
                [1, 0, 0, 0, 0]]
        atom = ['C', 'H', 'H', 'H', 'H']
        rdmol, smi, structure = cheminfo.create_rdkit_mol(bond, atom)
        warn = 'Unable to create RDKit mol from bond matrix'
        self.assertIsNotNone(rdmol, warn)
        self.assertEqual(len(atom) * 4, len(structure))

    def testInchiFromGeom(self):
        """
        Test the inchi creation from a geometry
        """
        smi = 'CCCC'
        inchi_expected = cheminfo.create_inchi_from_smi(smi)
        obmol, structure, bond = cheminfo.generate_3d_structure(smi)
        natom = len(structure) // 4
        structure = np.reshape(structure, (natom, 4))
        atom = structure[:, 0]
        geom = structure[:, 1:4].astype(float)
        for i, at in enumerate(atom):
            x, y, z = geom[i]
        inchi = cheminfo.create_inchi_from_geom(atom, geom)
        warn = 'Inchi generation from geometry failed.'
        self.assertEqual(inchi, inchi_expected, warn)

    def testInchiFromFile(self):
        """
        Test the inchi creation from a geometry
        """
        inchi_expected = 'InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3'
        file = 'input/cheminfo_xyz.xyz'
        inchi = cheminfo.create_inchi('', '', file)
        warn = 'Inchi generation from geometry failed.'
        self.assertEqual(inchi, inchi_expected, warn)

    def testInchiFromSmi(self):
        """
        Test the inchi creation from a smiles
        """
        smis = ['C', 'CCOO', 'CCCO[O]', 'C1CCCCC1']
        for smi in smis:
            inchi = cheminfo.create_inchi_from_smi(smi)
            warn = 'Unable to make inchi from smiles {}'.format(smi)
            self.assertIsNotNone(inchi, warn)

    def testSmiFromInchi(self):
        """
        Test the smiles creation from an inchi
        """
        inchis = ['InChI=1S/C6H12/c1-2-4-6-5-3-1/h1-6H2',
                  'InChI=1S/C3H8O/c1-3-4-2/h3H2,1-2H3',
                  'InChI=1S/C4H6/c1-3-4-2/h3-4H,1-2H2'
                  ]
        for inchi in inchis:
            smi = cheminfo.create_smiles(inchi)
            warn = 'Unable to make smiles from inchi {}'.format(inchi)
            self.assertIsNotNone(smi, warn)


if __name__ == "__main__":
    unittest.main()
