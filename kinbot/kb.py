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
import subprocess
import logging 
import datetime
import copy
import time

import license_message
from conformers import Conformers
from hindered_rotors import HIR
from parameters import Parameters
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
    par = Parameters('input.json')
    
    # set up the logging environment 
    if par.par['verbose']:
        logging.basicConfig(filename='kinbot.log',level=logging.DEBUG)
    else:
        logging.basicConfig(filename='kinbot.log',level=logging.INFO)
    
    #write the license message to the log file
    logging.info(license_message.message)
    #time stamp of the KinBot start
    logging.info('Starting KinBot at %s'%(datetime.datetime.now()))
    
    #initialize the parameters for this run
    par = Parameters('input.json')

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
    
    #initialize the master equation instance
    if par.par['me_code'] == 'mess':
        me = MESS(par)
    elif par.par['me_code'] == 'mesmer':
        me = MESMER(par)
    else:
        logging.error('Cannot recognize me code {}'.format(par.par['me_code']))

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
    
    """
    #do the conformational search
    if par.par['conformer_search'] == 1:
        logging.info('Starting conformational search of intial well')
        well0.confs = Conformers(well0,par,qc)
        #generate cyclic conformers, if any
        if len(well0.cycle_chain) > 0:
            geoms = well0.confs.generate_ring_conformers(copy.deepcopy(well0.geom))
            #wait for the cyclic conformational serach to finish
            #status, geoms = well0.confs.check_conformers(conf, wait = 1, ring = 1)
        else:
            geoms = [copy.deepcopy(well0.geom)]
        
        for geom in geoms:
            #for each cyclic conformer, generate and submit all non-cyclic conformers
            well0.confs.generate_conformers(0, geom)
        #wait for the conformational search to finish and select the lowest energy conformer
        status, well0.geom = well0.confs.check_conformers(conf, wait = 1, ring = 0)
    
    for it in range(3): #restart a high level calculation 3 times tops
        #high level calculations
        if par.par['high_level'] == 1:
            #re-optimize at high level of theory
            logging.info('Starting high level optimization of intial well')
            qc.qc_opt(well0, well0.geom, high_level = 1)
            err, well0.geom = qc.get_qc_geom(str(well0.chemid) + '_well_high', well0.natom, wait=1)
            err, well0.energy = qc.get_qc_energy(str(well0.chemid) + '_well_high', 1)
            #re-calculate the frequencies at high level of theory
            logging.info('Starting high level frequency calculation of intial well')
            qc.qc_freq(well0, well0.geom, high_level = 1)
            err, well0.freq = qc.get_qc_freq(str(well0.chemid) + '_fr_high', well0.natom, wait=1)
            err, well0.zpe = qc.get_qc_zpe(str(well0.chemid) + '_fr_high')
        
        # re-calculate the frequencies without internal rotations:
        if par.par['high_level'] == 1:
            fr_file = str(well0.chemid) + '_fr_high'
        else:
            fr_file = str(well0.chemid) + '_fr'
        hess = qc.read_qc_hess(fr_file,well0.natom)
        well0.kinbot_freqs, well0.reduced_freqs = frequencies.get_frequencies(well0, hess, well0.geom)
        

        
        
        #1D hindered rotor scans
        if par.par['rotor_scan'] == 1:
            logging.info('Starting hindered rotor calculations of intial well')
            well0.hir = HIR(well0,qc,par)
            well0.hir.generate_hir_geoms(copy.deepcopy(well0.geom))
            well0.hir.check_hir(wait = 1)
            if len(well0.hir.hir_energies) > 0:
                #check if along the hir potential a structure was found with a lower energy
                min = well0.hir.hir_energies[0][0]
                min_rotor = -1
                min_ai = -1
                for rotor in range(len(well0.dihed)):
                    for ai in range(well0.hir.nrotation):
                        if well0.hir.hir_energies[rotor][ai] < min - 1.6E-4: #use a 0.1kcal/mol cutoff for numerical noise 
                            min = well0.hir.hir_energies[rotor][ai]
                            min_rotor = rotor
                            min_ai = ai
                if min_rotor > -1:
                    #lower energy structure found
                    logging.info("Lower energy found during hindered rotor scan for initial well")
                    logging.info("Rotor: " + str(min_rotor))
                    logging.info("Scan point: " + str(min_ai))
                    job = 'hir/' + str(well0.chemid) + '_hir_' + str(min_rotor) + '_' + str(min_ai).zfill(2)
                    err,well0.geom = qc.get_qc_geom(job, well0.natom)
                    #delete the high_level log file and the hir log files
                    if os.path.exists(str(well0.chemid) + '_well_high.log'):
                        logging.info("Removing file " + str(well0.chemid) + '_well_high.log')
                        os.remove(str(well0.chemid) + '_well_high.log')
                    if os.path.exists(str(well0.chemid) + '_fr_high.log'):
                        logging.info("Removing file " + str(well0.chemid) + '_fr_high.log')
                        os.remove(str(well0.chemid) + '_fr_high.log')
                    for rotor in range(len(well0.dihed)):
                        for ai in range(well0.hir.nrotation):
                            if os.path.exists('hir/' + str(well0.chemid) + '_hir_' + str(rotor) + '_' + str(ai).zfill(2) + '.log'):
                                logging.info("Removing file " + 'hir/' + str(well0.chemid) + '_hir_' + str(rotor) + '_' + str(ai).zfill(2) + '.log')
                                os.remove('hir/' + str(well0.chemid) + '_hir_' + str(rotor) + '_' + str(ai).zfill(2) + '.log')
                else:
                    break
            else:
                break
        else:
            break
    """
    
        
    #write a mess block for the well
    me.write_well(well0, well0, {str(well0.chemid):well0.short_name})
    

    #do the reaction search using heuristics
    if par.par['reaction_search'] == 1:
        logging.info('Starting reaction searches of intial well')
        rf = ReactionFinder(well0,par,qc)
        rf.find_reactions()
        rg = ReactionGenerator(well0,par,qc)
        rg.generate()
        me.write(well0)


    logging.info('Finished KinBot at %s'%datetime.datetime.now())
    print("Done!")


if __name__ == "__main__":
    input_file = sys.argv[1]
    main(input_file)
