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
This class tests the symmetry calculation algorithm

The data is a dictionary of which the keys are the smiles
and the values are a dictionary with the structure and the 
expected symmetry numbers: the external rotational symmetry 
number, the internal rotational symmtry number and the 
number of optical isomers.
"""

import sys
import os
import unittest
import json
import numpy as np

from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint
import kinbot.symmetry


class TestSymmetry(unittest.TestCase):
    def setUp(self):
        pass
    
    def testAll(self):
        with open('symmetry_data.json') as f:
            data = json.load(f)
        messages = [
        'Expected external symmetry: {}, calculated: {}',
        'Expected internal symmetry: {}, calculated: {}',
        'Expected number of single events symmetry: {}, calculated: {}',]
        for name in data:
            par = Parameters()
            qc = QuantumChemistry(par)
            mol = StationaryPoint(name,0,1,smiles = name)
            mol.characterize()
            kinbot.symmetry.calculate_symmetry(mol)

            sigma_int = 1
            for row in mol.sigma_int:
                for at in row:
                    sigma_int *= at
            calc = [mol.sigma_ext, sigma_int, mol.nopt]
            
            for i in range(3):
                cal = calc[i]
                exp = data[name]['expected_values'][i]
                self.assertEqual(exp ,cal ,name + ': ' + messages[i].format(exp,cal))

if __name__ == "__main__":
    unittest.main()
    