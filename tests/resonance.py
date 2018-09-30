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

import sys
import os
import unittest
import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0,os.path.expanduser('~/ml-kinbot/code/kinbot'))

from stationary_pt import *
from cheminfo import *
from par import *
from symmetry import *


class TestResonance(unittest.TestCase):
    def setUp(self):
        pass
    
    def testAll(self):
        f = open('resonance_data.inp')
        par = imp.load_source('par', '', f)
        data = par.data
        
        for name in data:
            mol = stationary_pt(name)
            obmol, structure = generate_3d_structure(name)
            
            natom = len(structure)/4
            structure = np.reshape(structure, ( natom,4))
            atom = structure[:,0]
            mol.geom = structure[:,1:4].astype(float)
            
            mol.natom = natom
            mol.atom = atom
            mol.charge = 0
            mol.mult = mol.calc_multiplicity(atom)

            mol.bond_mx(mol.natom,mol.atom)
            cal = len(mol.bonds)
            exp = data[name]
            self.assertEqual(exp ,cal ,name + ': expected: {}, calculated: {}'.format(exp,cal))

if __name__ == "__main__":
    unittest.main()
    