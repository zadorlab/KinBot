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
This is the main class to run KinBot.
It includes 

1. Reading the options defined by the user
2. Optimizing the reactant, including frequency calculations,
hindered rotor calculations, high-level calculations, according
to the user input file
3. call the find_reactions method to find all reactions 
KinBot will search for
4. call the reac_generator method to search for the reactions
on the PES

"""
from __future__ import print_function
import sys
import os
import numpy as np
import re
import logging 
import datetime
import copy
import time

import license_message
import postprocess
from conformers import Conformers
from hindered_rotors import HIR
from homolytic_scissions import HomolyticScissions
from parameters import Parameters
from mesmer import MESMER
from mess import MESS
from optimize import Optimize
from reaction_finder import ReactionFinder
from reaction_generator import ReactionGenerator
from stationary_pt import StationaryPoint
from qc import QuantumChemistry

def main(input_file):
    #print the license message to the console
    print(license_message.message)
    
    #initialize the parameters for this run
    par = Parameters(input_file)
    
    # set up the logging environment 
    if par.par['verbose']:
        logging.basicConfig(filename='kinbot.log',level=logging.DEBUG)
    else:
        logging.basicConfig(filename='kinbot.log',level=logging.INFO)
    
    #write the license message to the log file
    logging.info(license_message.message)
    #time stamp of the KinBot start
    logging.info('Starting KinBot at %s'%(datetime.datetime.now()))

    #Make the necessary directories
    if not os.path.exists('perm'):
        os.makedirs('perm')
    if not os.path.exists('scratch'):
        os.makedirs('scratch')
    if not os.path.exists('molpro'):
        os.mkdir('molpro')
    if par.par['rotor_scan'] == 1:
        if not os.path.exists('hir'):
            os.mkdir('hir')
        if not os.path.exists('hir_profiles'):
            os.mkdir('hir_profiles')
        if not os.path.exists('perm/hir/'):
            os.makedirs('perm/hir/')
    if par.par['conformer_search'] == 1:
        if not os.path.exists('conf'):
            os.mkdir('conf')
        if not os.path.exists('perm/conf'):
            os.makedirs('perm/conf')
    if par.par['me'] == 1:
        if not os.path.exists('me'):
            os.mkdir('me')

    #initialize the reactant
    well0 = StationaryPoint('well0', par.par['charge'], par.par['mult'], smiles = par.par['smiles'], structure = par.par['structure'])
    well0.short_name = 'w1'
    
    #wrtie the initial reactant geometry to a file for visualization
    geom_out = open('geometry.xyz','w')
    geom_out.write('{}\n\n'.format(well0.natom))
    for i,at in enumerate(well0.atom):
        geom_out.write('{} {:.6f} {:.6f} {:.6f}\n'.format(at,well0.geom[i][0],well0.geom[i][1],well0.geom[i][2]))
    geom_out.write('\n\n')
    geom_out.close()

    #characterize the initial reactant
    well0.characterize()

    #initialize the qc instance
    qc = QuantumChemistry(par)

    #start the initial optimization of the reactant
    logging.info('Starting optimization of intial well')
    qc.qc_opt(well0, well0.geom)
    err, well0.geom = qc.get_qc_geom(str(well0.chemid) + '_well', well0.natom, wait=1)
    err, well0.freq = qc.get_qc_freq(str(well0.chemid) + '_well', well0.natom, wait=1)
    if err < 0: 
        logging.error('Error with initial structure optimization.')
        return 
    if any(well0.freq[i] <= 0 for i in range(len(well0.freq))):
        logging.error('Found imaginary frequency for initial structure.')
        return

    #do an MP2 optimization of the reactant, to compare Beta scission barrier heigths to
    logging.info('Starting MP2 optimization of intial well')
    qc.qc_opt(well0, well0.geom, mp2 = 1)
    err, geom = qc.get_qc_geom(str(well0.chemid) + '_well_mp2', well0.natom, 1)
    
    #characterize again and look for differences
    well0.characterize()
    well0.name = str(well0.chemid)
    
    #read the energy and the zpe corrected energy
    err, well0.energy = qc.get_qc_energy(str(well0.chemid) + '_well', 1)
    err, well0.zpe = qc.get_qc_zpe(str(well0.chemid) + '_well', 1)
    
    well_opt = Optimize(well0, par, qc, wait = 1)
    well_opt.do_optimization()
    
    #do the reaction search using heuristics
    if par.par['reaction_search'] == 1:
        logging.info('Starting reaction searches of intial well')
        rf = ReactionFinder(well0,par,qc)
        rf.find_reactions()
        rg = ReactionGenerator(well0,par,qc)
        rg.generate()
    #do the homolytic scission products search
    if par.par['homolytic_scissions'] == 1:
        logging.info('Starting the search for homolytic scission products')
        well0.homolytic_scissions = HomolyticScissions(well0,par,qc)
        well0.homolytic_scissions.find_homolytic_scissions()
    if par.par['me'] == 1:
        logging.info('Starting Master Equation calculations')
        #initialize the master equation instance
        if par.par['me_code'] == 'mess':
            me = MESS(par,well0)
        elif par.par['me_code'] == 'mesmer':
            me = MESMER(par,well0)
        else:
            logging.error('Cannot recognize me code {}'.format(par.par['me_code']))
        me.write_input()
        me.run()
        
    
    #postprocess the calculations
    postprocess.createSummaryFile(well0,qc,par)
    postprocess.createPESViewerInput(well0,qc,par)


    logging.info('Finished KinBot at %s'%datetime.datetime.now())
    print("Done!")


if __name__ == "__main__":
    input_file = sys.argv[1]
    main(input_file)
