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
This class tests the conformational an hindered rotor dihedral search

The data is a dictionary of which the keys are the smiles
and the values are the expected number of resonance isomers
"""

import sys
import os
import unittest
import imp
import numpy as np

from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint

class TestDihedrals(unittest.TestCase):
    def setUp(self):
        pass
    
    def testAll(self):
        f = open('dihedral_data.inp')
        par = imp.load_source('par', '', f)
        data = par.data
        
        for name in data:
            par = Parameters()
            qc = QuantumChemistry(par)
            mol = StationaryPoint(name,0,1,smiles = name)
            mol.characterize()
            
            hir_exp = data[name][0]
            conf_exp = data[name][1]
            hir_calc = len(mol.dihed)
            conf_calc = len(mol.conf_dihed)
            self.assertEqual(hir_exp ,hir_calc ,name + ': HIR, expected: {}, calculated: {}'.format(hir_exp,hir_calc))
            self.assertEqual(conf_exp ,conf_calc ,name + ': CONF, expected: {}, calculated: {}'.format(conf_exp,conf_calc))
if __name__ == "__main__":
    unittest.main()
    