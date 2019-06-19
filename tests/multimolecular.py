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
import json
import unittest
import numpy as np

from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint

class TestMultimolecular(unittest.TestCase):
    def setUp(self):
        pass
        
    def testAll(self):
        with open('multimolecular_data.json') as f:
            data = json.load(f)
        for name in data:
            print(name)
            par = Parameters()
            qc = QuantumChemistry(par)
            structure = data[name]['structure']
            mol = StationaryPoint(name,0,1,structure = structure)
            mol.characterize()

            mols = mol.start_multi_molecular()
            calculated = len(mols)
            expected = data[name]['expected_value']
            self.assertEqual(calculated,expected, name + ': expected: {}, calculated: {}'.format(expected,calculated))

if __name__ == '__main__':
    unittest.main()
