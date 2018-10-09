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
This class tests the conformer generation algorithm

The data is a dictionary of which the keys are the smiles
and the values are a dictionary with the structure and the 
expected symmetry numbers: the external rotational symmetry 
number, the internal rotational symmtry number and the 
number of optical isomers.
"""

import sys
import os
import unittest
import copy
import time
import numpy as np

from kinbot.stationary_pt import StationaryPoint
from kinbot.conformers import Conformers
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry

class TestConformers(unittest.TestCase):
    def setUp(self):
        pass
    
    def testRingConformers(self):
        """
        The generation of ring conformers requires qc calculations and is therefore slow!
        """
        if not os.path.exists('conf/'):
            os.mkdir('conf/')
        par = Parameters()
        qc = QuantumChemistry(par)
        
        smis = ['C1C(C)CCCC1']
        
        for smi in smis:
            mol = StationaryPoint(smi,0,1,smiles = smi)
            mol.confs = Conformers(mol,par,qc)
            mol.characterize()
            
            scycconf = -1
            while 1:
                if scycconf == -1: 
                    #start the ring conf search
                    if len(mol.cycle_chain) > 0: 
                        #there are rings in the molecule, do a search
                        mol.confs.generate_ring_conformers(copy.deepcopy(mol.geom))
                        #set the cyclic conf status to running
                        scycconf = 0 
                    else: 
                        #there are no rings in the molecule, continue from the current geometry
                        mol.confs.cyc_conf_geoms = [copy.deepcopy(mol.geom)]
                        #no ring conf search has to be done, set status to finished
                        scycconf = 1 
                if scycconf == 0: 
                    #ring conf search is running, check if finished
                    #TODO: change the cyc conf architecture
                    #ring conf search is finished
                    status,geoms = mol.confs.check_ring_conformers()
                    if status == 1:
                        scycconf = 1 
                if scycconf == 1:
                    print 'breaking'
                    break
                time.sleep(1)
        self.assertEqual(1 ,1)
        
    """

    def testAll(self):
        f = open('symmetry_data.inp')
        par = imp.load_source('par', '', f)
        data = par.data
        messages = [
        'Expected external symmetry: {}, calculated: {}',
        'Expected internal symmetry: {}, calculated: {}',
        'Expected number of single events symmetry: {}, calculated: {}',]
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

            mol.characterize(mol.natom,mol.atom,mol.mult,mol.charge)
            calculate_symmetry(mol, natom, atom)

            sigma_int = 1
            for row in mol.sigma_int:
                for at in row:
                    sigma_int *= at
            calc = [mol.sigma_ext, sigma_int, mol.nopt]
            
            for i in range(3):
                exp = data[name]['expected_values'][i]
                cal = calc[i]
                self.assertEqual(exp ,cal ,name + ': ' + messages[i].format(exp,cal))
    """
if __name__ == "__main__":
    unittest.main()
    