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
import unittest
import numpy as np

import kinbot.geometry as geometry
from kinbot.stationary_pt import StationaryPoint


class TestGeometry(unittest.TestCase):
    def setUp(self):
        pass

    def testCalcAngle(self):
        """
        Test the angle calculation
        """
        # calculate a random anle
        np.random.seed(1)
        a = np.random.uniform(size=3)
        b = np.random.uniform(size=3)
        c = np.random.uniform(size=3)
        cal = geometry.calc_angle(a, b, c)
        exp = 1.24358021797
        warn = 'Angle calc faild: exp: {}, calc: {}'.format(exp, cal)
        self.assertAlmostEqual(exp, cal, places=10, msg=warn)

        # calculate an angle between to collinear vectors
        np.random.seed(1)
        a = np.random.uniform(size=3)
        b = np.array([0., 0., 0.])
        cal = geometry.calc_angle(a, b, a)
        exp = 0.
        warn = 'Angle calc faild: exp: {}, calc: {}'.format(exp, cal)
        self.assertAlmostEqual(exp, cal, places=10, msg=warn)

        # calculate an angle between to inverse vectors
        np.random.seed(1)
        a = np.random.uniform(size=3)
        b = np.array([0., 0., 0.])
        c = -a
        cal = geometry.calc_angle(a, b, c)
        exp = np.pi
        warn = 'Angle calc faild: exp: {}, calc: {}'.format(exp, cal)
        self.assertAlmostEqual(exp, cal, places=10, msg=warn)

        # calculate an angle between to orthogonal vectors
        a = np.array([1., 0., 0.])
        b = np.array([0., 0., 0.])
        c = np.array([0., 1., 0.])
        cal = geometry.calc_angle(a, b, c)
        exp = np.pi / 2
        warn = 'Angle calc faild: exp: {}, calc: {}'.format(exp, cal)
        self.assertAlmostEqual(exp, cal, places=10, msg=warn)

    def testLinear(self):
        """
        Test the method to see if a geometry is linear
        """
        np.random.seed(1)
        a = np.random.uniform(size=3)
        b = np.random.uniform(size=3)
        d1 = np.random.rand()
        d2 = np.random.rand()
        c = a + b / np.linalg.norm(b) * d1
        d = a + b / np.linalg.norm(b) * d2
        geom = [a, c, d]
        bond = [[0, 1, 0], [1, 0, 1], [0, 0, 1]]
        dummy = geometry.is_linear(geom, bond)
        exp = 1
        cal = len(dummy)
        warn = 'Linear molecule not well perceived.'
        self.assertEqual(exp, cal, msg=warn)

    def testOutOfPlaneAngle(self):
        """
        Test the out-of-plane angle calculation
        """
        np.random.seed(1)
        a = np.random.uniform(size=3)
        b = np.random.uniform(size=3)
        c = np.random.uniform(size=3)
        d = np.random.uniform(size=3)
        cal, collinear = geometry.calc_out_of_plane_angle(a, b, c, d)
        # verify the angle
        exp = -0.602981717275
        warn = 'Out-of-plane angle not correct: '
        warn += 'expected {}, calculated {}'.format(exp, cal)
        self.assertAlmostEqual(exp, cal, places=10, msg=warn)
        # verify the collinear boolean
        exp = 0
        cal = collinear
        warn = 'Out-of-plane collinear boolean is not correct: '
        warn += 'expected {}, calculated {}'.format(exp, cal)
        self.assertEqual(exp, cal, msg=warn)

    def testCalcDihedral(self):
        """
        Test the dihedral angle calculation
        """
        np.random.seed(1)
        a = np.random.uniform(size=3)
        b = np.random.uniform(size=3)
        c = np.random.uniform(size=3)
        d = np.random.uniform(size=3)
        cal, collinear = geometry.calc_dihedral(a, b, c, d)
        # verify the angle
        exp = -83.9898211154
        warn = 'Dihedral angle not correct: '
        warn += 'expected {}, calculated {}'.format(exp, cal)
        self.assertAlmostEqual(exp, cal, places=10, msg=warn)
        # verify the collinear boolean
        exp = 0
        cal = collinear
        warn = 'Dihedral collinear boolean is not correct: '
        warn += 'expected {}, calculated {}'.format(exp, cal)
        self.assertEqual(exp, cal, msg=warn)

    def testNewRingDihedrals(self):
        """
        Test the calculation of new dihedrals necessary for the update
        """
        smi = 'CCCC'
        mol = StationaryPoint(smi, 0, 1, smiles=smi)
        mol.characterize()
        ins = [0, 1, 2, 3]  # change the C-C-C-C dihedral
        step_nr = 0  # This corresponds to the first dihedral update
        # there are 12 steps done in total, this means that the dihedral
        # angle should be changed by 1/10 of the total change
        total_nr_of_steps = 10
        # final dihedral value we are shooting for after 10 updates
        # this is a default value for a instance of 5 or less atoms
        final_val = 15.
        # initial dihedral
        ini = geometry.calc_dihedral(mol.geom[ins[0]],
                                     mol.geom[ins[1]],
                                     mol.geom[ins[2]],
                                     mol.geom[ins[3]])[0]
        # this is the new dihedral angle after one step
        update = geometry.new_ring_dihedrals(mol, ins, step_nr,
                                             total_nr_of_steps)
        exp = ini - (ini - final_val) / 10
        cal = update[0]
        warn = 'Dihedral update is not correct: '
        warn += 'expected {}, calculated {}'.format(exp, cal)
        self.assertEqual(exp, cal, msg=warn)

    def testNewBondLength(self):
        """
        Test the calculation of new bond length necessary for the update
        """
        smi = 'CCCC'
        mol = StationaryPoint(smi, 0, 1, smiles=smi)
        mol.characterize()
        ati = 0
        atj = 1
        step_nr = 1  # This corresponds to the first bond length update
        # there are 12 steps done in total, this means that the bond
        # length should be changed by 1/10 of the total change
        total_nr_of_steps = 10
        # final bond length we are shooting for after 10 updates
        final_val = 2.0
        # initial bond length
        ini = np.linalg.norm(mol.geom[ati] - mol.geom[atj])
        # this is the new bond length one step
        cal = geometry.new_bond_length(mol, ati, atj, step_nr,
                                       total_nr_of_steps, final_val)
        exp = ini + (final_val - ini) / 10
        warn = 'Dihedral update is not correct: '
        warn += 'expected {}, calculated {}'.format(exp, cal)
        self.assertEqual(exp, cal, msg=warn)

    def testTranslateAndRotate(self):
        """
        Test the translation and rotation of a molecule
        """
        cart = np.array([[2., 2., 2.],
                         [3., 2., 2.],
                         [2., 3., 2.],
                         [2., 2., 3.],
                         [1., 1., 1.]])
        atom = []  # this argument is deprecated
        i = 0
        j = 1
        cal = geometry.translate_and_rotate(cart, atom, i, j)
        exp = np.array([[0., 0., 0.],
                        [0., 0., 1.],
                        [0., 1., 0.],
                        [-1., 0., 0.],
                        [1., -1., -1.]])
        for i, ci in enumerate(cal):
            for j, cij in enumerate(ci):
                self.assertAlmostEqual(exp[i][j], cij, places=10)

    def testCenterOfMass(self):
        """
        Test the center of mass calculation
        """
        atom = ['C', 'C', 'C', 'C']
        geom = np.array([[0., 0., 0.],
                         [1., 0., 0.],
                         [-1., 0., 0.],
                         [0., 2., 0.]])
        cal = geometry.get_center_of_mass(geom, atom)
        exp = [0., 0.5, 0.]
        for i, ci in enumerate(cal):
            self.assertAlmostEqual(exp[i], ci, places=10)

    def testRotateAtom(self):
        """
        Test the rotation of an atom
        """
        # vector to rotate
        v = np.array([1., 2., 2.])
        # normal vector to rotation
        n = np.array([1., 0., 0.])
        # angle of rotation
        th = np.radians(90.)
        cal = geometry.rotate_atom(v, n, th)
        exp = [1., -2., 2.]
        for i, ci in enumerate(cal):
            self.assertAlmostEqual(exp[i], ci, places=10)

    def testGetMomentsOfInertia(self):
        """
        Test the moments of inertia calculator
        """
        atom = ['C', 'C', 'C', 'C']
        geom = np.array([[0., 0., 0.],
                         [1., 0., 0.],
                         [-1., 0., 0.],
                         [0., 2., 0.]])
        eigvals, eigvecs = geometry.get_moments_of_inertia(geom, atom)
        exp_eigvals = [24, 36, 60]
        exp_eigvecs = np.array([[0., 1., 0.],
                                [1., 0., 0.],
                                [0., 0., 1.]])
        for i, ci in enumerate(eigvecs):
            for j, cij in enumerate(ci):
                self.assertAlmostEqual(exp_eigvecs[i][j], cij, places=10)
        for i, ci in enumerate(eigvals):
            self.assertAlmostEqual(exp_eigvals[i], ci, places=10)

    def testEqualGeom(self):
        """
        Test the algorithm to verify if two geometries are identical
        considering a bond-length cutoff of 10%
        """
        bond = [[0, 1, 1, 1],
                [1, 0, 1, 0],
                [1, 1, 0, 0],
                [1, 0, 0, 0]]
        geom1 = np.array([[0., 0., 0.],
                          [1., 0., 0.],
                          [-1., 0., 0.],
                          [0., 2., 0.]])
        geom2 = np.array([[0., 0., 0.],
                          [1.09, 0., 0.],
                          [-1., 0., 0.],
                          [0., 2., 0.]])
        geom3 = np.array([[0., 0., 0.],
                          [1.11, 0., 0.],
                          [-1., 0., 0.],
                          [0., 2., 0.]])
        warn = 'Verifying if two geometries are equal failed'
        equal = geometry.equal_geom(bond, geom1, geom2, 0.1)
        self.assertTrue(equal, warn)
        equal = geometry.equal_geom(bond, geom1, geom3, 0.1)
        self.assertFalse(equal, warn)


if __name__ == "__main__":
    unittest.main()
