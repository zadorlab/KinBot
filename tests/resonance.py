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

from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint


class TestResonance(unittest.TestCase):
    def setUp(self):
        pass
    
    def testAll(self):
        data = {"C1=CC=CC=C1":2,
            "C1=CC=C(C)C=C1":2,
            "C=C[CH2]":2,
            "C=C=C":1,
            "C#C[CH2]":2,
            "S=S":1,
            "O=S=C":1,
            "O=S(C)[CH2]":3,
            "C1CC=CC=C1":1
            }

        for name in data:
            par = Parameters()
            qc = QuantumChemistry(par)
            mol = StationaryPoint(name,0,1,smiles = name)
            mol.characterize()

            cal = len(mol.bonds)
            exp = data[name]
            self.assertEqual(exp ,cal ,name + ': expected: {}, calculated: {}'.format(exp,cal))

if __name__ == "__main__":
    unittest.main()
    