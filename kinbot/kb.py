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
import sys
import os
import numpy as np
import re
import subprocess
import logging 
import datetime
import copy
import time


from parameters import Parameters
import license_message

"""
from vector import *
from motif import *
from reaction import *
from stationary_pt import *
from conformers import *
from zmat import *
import constants
from write_mess import *
from high_level import *
from geom import *
from hir import *
from symmetry import *
from frequencies import *
from molpro import *
"""


def main():

    print license_message.message

    #sys.path.append(os.path.expanduser('~/KinBot/code'))
    
    
    #os.chdir(os.path.expanduser(workdir))
    
    # set up the logging environment 
    logging.basicConfig(filename='kinbot.log',level=logging.INFO)
    
    logging.info(license_message.message)
    logging.info('Starting KinBot at %s'%(datetime.datetime.now()))
    
    par = Parameters('input.json')
    print 'title:'
    print par.par['title']
    print ''
    
    print par.print_parameters()
    sys.exit(-1)
    if not os.path.exists('perm'):
        os.makedirs('perm')

    if not os.path.exists('scratch'):
        os.makedirs('scratch')
        
    if not os.path.exists('molpro'):
        os.mkdir('molpro')


    well0 = stationary_pt('well0')
    well0.short_name = 'w1'
    
    well0.geom = read_input('input.inp')
    
    if par.rotor_scan == 1:
        if not os.path.exists('hir'):
            os.mkdir('hir')
        if not os.path.exists('hir_profiles'):
            os.mkdir('hir_profiles')
        if not os.path.exists('perm/hir/'):
            os.makedirs('perm/hir/')
    
    if par.conformer_search == 1:
        if not os.path.exists('conf'):
            os.mkdir('conf')
        if not os.path.exists('perm/conf'):
            os.makedirs('perm/conf')
        
    geom_out = open('geometry.xyz','w')
    
    geom_out.write('{}\n\n'.format(par.natom))
    for i,at in enumerate(par.atom):
        geom_out.write('{} {:.6f} {:.6f} {:.6f}\n'.format(at,well0.geom[i][0],well0.geom[i][1],well0.geom[i][2]))
    geom_out.write('\n\n')
    
    geom_out.close()

    well0.characterize(par.natom, par.atom, par.mult, par.charge)

    logging.info('Starting optimization of intial well')
    qc_opt(well0, well0.geom, 0, par.natom, par.atom, par.mult, par.charge)
    
    
    err, well0.geom = get_qc_geom(str(well0.chemid) + '_well', par.natom, 1)
    if err < 0: 
        logging.error('Error with initial structure optimization.')
        return 
    

    logging.info('Starting frequency calculation of intial well')
    qc_freq(well0, well0.geom, par.natom, par.atom, par.mult, par.charge)
    err, well0.freq = get_qc_freq(str(well0.chemid) + '_fr', par.natom, 1)
    if err < 0:
        logging.error('Error with initial structure freqency calculations.')
        return
    if any(well0.freq[i] <= 0 for i in range(len(well0.freq))):
        logging.error('Found imaginary frequency for initial structure.')
        return
    
    err, well0.zpe = get_qc_zpe(str(well0.chemid) + '_fr')

    well0.characterize(par.natom, par.atom, par.mult, par.charge)

    err, well0.energy = get_qc_energy(str(well0.chemid) + '_well', 1)
    err, well0.zpe = get_qc_zpe(str(well0.chemid) + '_fr', 1)
    #err, well0.nmode = get_qc_nmode(str(well0.chemid) + '_fr', par.natom)
    #well0.scaled_nmode = nmode_scale(well0.freq, well0.nmode)

    #conformational search
    if par.conformer_search == 1:
        logging.info('Starting conformational search of intial well')
        #generate cyclic conformers, if any
        if len(well0.cycle_chain) > 0:
            geoms = generate_ring_conformers(well0, par.natom, par.atom, par.mult, par.charge, copy.deepcopy(well0.geom), 0)
            #wait for the cyclic conformational serach to finish
            #status, geoms = check_conformers(well0, conf, par.natom, par.atom, par.mult, par.charge,wait = 1, ring = 1)
        else:
            geoms = [copy.deepcopy(well0.geom)]
        
        conf = 0
        for geom in geoms:
            #for each cyclic conformer, generate and submit all non-cyclic conformers
            conf = generate_conformers(well0, par.natom, par.atom, par.mult, par.charge, 0, geom, conf, 0)
        #wait for the conformational search to finish and select the lowest energy conformer
        status, well0.geom = check_conformers(well0, conf, par.natom, par.atom, par.mult, par.charge,wait = 1)
    
    #do an MP2 optimization of the reactant, to compare Beta scission barrier heigths to
    logging.info('Starting MP2 optimization of intial well')
    qc_opt(well0, well0.geom, 0, par.natom, par.atom, par.mult, par.charge, mp2 = 1)
    err, geom = get_qc_geom(str(well0.chemid) + '_well_mp2', par.natom, 1)
    
    for it in range(3): #restart a high level calculation 3 times tops
        #high level calculations
        if par.high_level == 1:
            #re-optimize at high level of theory
            logging.info('Starting high level optimization of intial well')
            qc_opt(well0, well0.geom, 0, par.natom, par.atom, par.mult, par.charge, high_level = 1)
            err, well0.geom = get_qc_geom(str(well0.chemid) + '_well_high', par.natom, 1)
            err, well0.energy = get_qc_energy(str(well0.chemid) + '_well_high', 1)
            #re-calculate the frequencies at high level of theory
            logging.info('Starting high level frequency calculation of intial well')
            qc_freq(well0, well0.geom, par.natom, par.atom, par.mult, par.charge, high_level = 1)
            err, well0.freq = get_qc_freq(str(well0.chemid) + '_fr_high', par.natom, 1)
            err, well0.zpe = get_qc_zpe(str(well0.chemid) + '_fr_high')
        
        # re-calculate the frequencies without internal rotations:
        if par.high_level == 1:
            fr_file = str(well0.chemid) + '_fr_high'
        else:
            fr_file = str(well0.chemid) + '_fr'
        hess = read_qc_hess(fr_file,par.natom)
        well0.kinbot_freqs, well0.reduced_freqs = get_frequencies(well0, hess, par.natom, par.atom, well0.geom)
        

        #calculate the symmetry numbers
        calculate_symmetry(well0, par.natom, par.atom)
        
        #1D hindered rotor scans
        if par.rotor_scan == 1:
            logging.info('Starting hindered rotor calculations of intial well')
            generate_hir_geoms(well0, par.natom, par.atom, par.mult, par.charge, copy.deepcopy(well0.geom), 0)
            check_hir(well0,par.natom, par.atom, par.mult, par.charge, 0,wait = 1)
            if len(well0.hir_energies) > 0:
                #check if along the hir potential a structure was found with a lower energy
                min = well0.hir_energies[0][0]
                min_rotor = -1
                min_ai = -1
                for rotor in range(len(well0.dihed)):
                    for ai in range(par.nrotation):
                        if well0.hir_energies[rotor][ai] < min - 1.6E-4: #use a 0.1kcal/mol cutoff for numerical noise 
                            min = well0.hir_energies[rotor][ai]
                            min_rotor = rotor
                            min_ai = ai
                if min_rotor > -1:
                    #lower energy structure found
                    logging.info("Lower energy found during hindered rotor scan for initial well")
                    logging.info("Rotor: " + str(min_rotor))
                    logging.info("Scan point: " + str(min_ai))
                    job = 'hir/' + str(well0.chemid) + '_hir_' + str(min_rotor) + '_' + str(min_ai).zfill(2)
                    err,well0.geom = get_qc_geom(job, par.natom)
                    #delete the high_level log file and the hir log files
                    if os.path.exists(str(well0.chemid) + '_well_high.log'):
                        logging.info("Removing file " + str(well0.chemid) + '_well_high.log')
                        os.remove(str(well0.chemid) + '_well_high.log')
                    if os.path.exists(str(well0.chemid) + '_fr_high.log'):
                        logging.info("Removing file " + str(well0.chemid) + '_fr_high.log')
                        os.remove(str(well0.chemid) + '_fr_high.log')
                    for rotor in range(len(well0.dihed)):
                        for ai in range(par.nrotation):
                            if os.path.exists('hir/' + str(well0.chemid) + '_hir_' + str(rotor) + '_' + str(ai).zfill(2) + '.log'):
                                logging.info("Removing file " + 'hir/' + str(well0.chemid) + '_hir_' + str(rotor) + '_' + str(ai).zfill(2) + '.log')
                                os.remove('hir/' + str(well0.chemid) + '_hir_' + str(rotor) + '_' + str(ai).zfill(2) + '.log')
                else:
                    break
            else:
                break
        else:
            break

    #write a mess block for the well
    write_mess_well(well0, par.natom, par.atom, par.mult, par.charge, well0, {str(well0.chemid):well0.short_name})
    
    #write a zmat for EStoKTp
    estoktp_zmat(well0,str(well0.chemid) + '_well',par.natom,par.atom)
    
    #write the molpro input and read the molpro energy, if available
    create_molpro_input(well0, par.natom, par.atom, par.mult, par.charge, 0)
    status, molpro_energy = get_molpro_energy(well0, 0)
    if status:
        well0.energy = molpro_energy

    #do the reaction search using heuristics
    if par.reaction_search == 1:
        logging.info('Starting reaction searches of intial well')
        well0.find_reactions(par.natom, par.atom)
        reac_generator(well0)


    logging.info('Finished KinBot at %s'%datetime.datetime.now())
    print "Done!"


if __name__ == "__main__":
    #workdir = '~/' + sys.argv[1]
    #main(workdir)
    main()
