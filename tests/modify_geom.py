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
This class tests the geometry modification algorithm 

"""

import sys
import os
import unittest
import numpy as np

import kinbot.modify_geom
from kinbot.stationary_pt import StationaryPoint
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry

class TestGeometryModification(unittest.TestCase):
    def setUp(self):
        self.write_files = 0

    def testDihedralChangeHexane(self):
        """
        The generation of ring conformers requires qc calculations and is therefore slow!
        """
        if not os.path.exists('conf/'):
            os.mkdir('conf/')
        par = Parameters()
        qc = QuantumChemistry(par)
        smi = 'CCCCCC'
        mol = StationaryPoint(smi,0,1,smiles = smi)
        mol.characterize()
        changes = [
        [0,1,2,3,15.],
        [1,2,3,4,15.],
        [2,3,4,5,15.],
        ]
        name = 'hexane_dihedral'
        success, new_geom = kinbot.modify_geom.modify_coordinates(mol,name,mol.geom,changes,mol.bond,write_files = self.write_files)

    def testBondChangeEthane(self):
        """
        The generation of a longer C-C bond in ethane
        """
        if not os.path.exists('conf/'):
            os.mkdir('conf/')
        par = Parameters()
        qc = QuantumChemistry(par)
        smi = 'CC'
        mol = StationaryPoint(smi,0,1,smiles = smi)
        mol.characterize()
        changes = [
        [0,1,1.8],
        ]
        name = 'ethane_bond_length_test'
        success, new_geom = kinbot.modify_geom.modify_coordinates(mol,name,mol.geom,changes,mol.bond,write_files = self.write_files)

    def testDihedralChangeHeptyl(self):
        """
        The generation of ring conformers requires qc calculations and is therefore slow!
        """
        if not os.path.exists('conf/'):
            os.mkdir('conf/')
        par = Parameters()
        qc = QuantumChemistry(par)
        smi = '[CH2]CCCCC'
        mol = StationaryPoint(smi,0,2,smiles = smi)
        mol.characterize()
        changes = [
        [0,3,4,5,25.],
        [3,4,5,6,25.],
        [4,5,6,7,25.],
        [5,6,7,16,25.],
        ]
        name = 'hexyl_dihedral'
        success, new_geom = kinbot.modify_geom.modify_coordinates(mol,name,mol.geom,changes,mol.bond,write_files = self.write_files)

    def testDihedralAndBondChangeHeptyl(self):
        """
        The generation of ring conformers requires qc calculations and is therefore slow!
        """
        if not os.path.exists('conf/'):
            os.mkdir('conf/')
        par = Parameters()
        qc = QuantumChemistry(par)
        smi = '[CH2]CCCCC'
        mol = StationaryPoint(smi,0,2,smiles = smi)
        mol.characterize()
        changes = [
        [0,3,4,5,20.],
        [3,4,5,6,20.],
        [4,5,6,7,20.],
        [5,6,7,16,20.],
        ]
        name = 'hexyl_dihedral_and_bond'
        success, new_geom = kinbot.modify_geom.modify_coordinates(mol,name,mol.geom,changes,mol.bond,write_files = self.write_files)
        changes = [
        [0,16,1.35],
        [7,16,1.35],
        ]
        success, new_geom2 = kinbot.modify_geom.modify_coordinates(mol,name,new_geom,changes,mol.bond,write_files = self.write_files)

    def testAngleDihedralAndBondChangeHeptyl(self):
        """
        The generation of ring conformers requires qc calculations and is therefore slow!
        """
        if not os.path.exists('conf/'):
            os.mkdir('conf/')
        par = Parameters()
        qc = QuantumChemistry(par)
        smi = '[CH2]CCCCC'
        mol = StationaryPoint(smi,0,2,smiles = smi)
        mol.characterize()
        changes = [
        [0,3,4,128.],
        [3,4,5,128.],
        [4,5,6,128.],
        [5,6,7,128.],
        [6,7,16,128.],
        ]
        name = 'hexyl_angle_dihedral_and_bond'
        success, new_geom = kinbot.modify_geom.modify_coordinates(mol,name,mol.geom,changes,mol.bond,write_files = self.write_files)
        changes = [
        [0,3,4,5,15.],
        [3,4,5,6,15.],
        [4,5,6,7,15.],
        [5,6,7,16,15.],
        ]
        success, new_geom2 = kinbot.modify_geom.modify_coordinates(mol,name,new_geom,changes,mol.bond,write_files = self.write_files)
        changes = [
        [0,16,1.35],
        [7,16,1.35],
        ]
        success, new_geom3 = kinbot.modify_geom.modify_coordinates(mol,name,new_geom2,changes,mol.bond,write_files = self.write_files)

    def testAngleChangeHeptyl(self):
        """
        The generation of ring conformers requires qc calculations and is therefore slow!
        """
        if not os.path.exists('conf/'):
            os.mkdir('conf/')
        par = Parameters()
        qc = QuantumChemistry(par)
        smi = '[CH2]CCCCC'
        mol = StationaryPoint(smi,0,2,smiles = smi)
        mol.characterize()
        changes = [
        [0,3,4,128.],
        [3,4,5,128.],
        [4,5,6,128.],
        [5,6,7,128.],
        [6,7,16,128.],
        ]
        name = 'hexyl_angle'
        success, new_geom = kinbot.modify_geom.modify_coordinates(mol,name,mol.geom,changes,mol.bond,write_files = self.write_files)


if __name__ == "__main__":
    unittest.main()
    