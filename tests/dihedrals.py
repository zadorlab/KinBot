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
import unittest

from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint

class TestDihedrals(unittest.TestCase):
    def setUp(self):
        pass
    
    def testAll(self):
        data = {
            'CC':[1,0],
            'CCC':[2,0],
            'CCCC':[3,1],
            'C=C':[0,0],
            'C=CC':[1,0],
            'C=C[CH2]':[0,0],
            'CC=C[CH2]':[1,0],
            'C1CCCC1':[0,0],
            'CO':[1,0],
            'C=CO':[1,1],
            }
        
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
    