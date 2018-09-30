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
This class tests the multimolecular search functionality of KinBot.
The input is a geometry from either a molecule or a set of molecules.
The return value is True if it found the expected number of molecules per geometry. 
"""

import os, sys, imp
import unittest

sys.dont_write_bytecode = True
sys.path.insert(0,os.path.expanduser('~/ml-kinbot/code/kinbot'))

from stationary_pt import *
from par import *

class TestMultimolecular(unittest.TestCase):
    def setUp(self):
        pass
        
    def testAll(self):
        f = open('multimolecular_data.inp')
        par = imp.load_source('par', '', f)
        data = par.data
        for name in data:
            mol = stationary_pt(name)
            natom = len(data[name]['structure'])/4
            structure = np.reshape(data[name]['structure'], ( natom,4))
            atom = structure[:,0]
            
            mol.geom = structure[:,1:4].astype(float)
            mol.natom = natom
            mol.atom = atom
            mol.charge = 0
            mol.mult = mol.calc_multiplicity(atom)
            
            par.mult = mol.mult
            par.charge = mol.charge
            par.natom = natom
            par.atom = atom

            mols = mol.start_multi_molecular(natom,atom)
            calculated = len(mols)
            expected = data[name]['expected_value']
            self.assertEqual(calculated,expected, name + ': expected: {}, calculated: {}'.format(expected,calculated))

if __name__ == '__main__':
    unittest.main()
