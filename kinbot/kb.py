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
import logging
import datetime

from kinbot import license_message
from kinbot import postprocess
from kinbot.homolytic_scissions import HomolyticScissions
from kinbot.parameters import Parameters
from kinbot.mesmer import MESMER
from kinbot.mess import MESS
from kinbot.optimize import Optimize
from kinbot.reaction_finder import ReactionFinder
from kinbot.reaction_generator import ReactionGenerator
from kinbot.stationary_pt import StationaryPoint
from kinbot.qc import QuantumChemistry


def main():
    try:
        input_file = sys.argv[1]
    except IndexError:
        print('To use KinBot, supply one argument being the input file!')
        sys.exit(-1)

    # print the license message to the console
    print(license_message.message)

    # initialize the parameters for this run
    par = Parameters(input_file)

    # set up the logging environment
    if par.par['verbose']:
        logging.basicConfig(filename='kinbot.log', level=logging.DEBUG)
    else:
        logging.basicConfig(filename='kinbot.log', level=logging.INFO)

    # write the license message to the log file
    logging.info(license_message.message)
    # time stamp of the KinBot start
    logging.info('Starting KinBot at {}'.format(datetime.datetime.now()))

    # Make the necessary directories
    if not os.path.exists('perm'):
        os.makedirs('perm')
    if not os.path.exists('scratch'):
        os.makedirs('scratch')
    if not os.path.exists(par.par['single_point_qc']):
        os.mkdir(par.par['single_point_qc'])
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
    if not os.path.exists('me'):
        os.mkdir('me')


    if par.par['pes'] and par.par['specific_reaction']:
        logging.error('Specific reaction cannot be searched in PES mode.')
        return

    if !par.par['break_bonds'] and !par.par['form.bonds'] and par.par['specific_reaction']:
        logging.error('For pecific reaction search you need to define at least one breaking or forming bond.')
        return

    # initialize the reactant
    well0 = StationaryPoint('well0',
                            par.par['charge'],
                            par.par['mult'],
                            smiles=par.par['smiles'],
                            structure=par.par['structure'])
    well0.short_name = 'w1'

    # wrtie the initial reactant geometry to a file for visualization
    geom_out = open('geometry.xyz', 'w')
    geom_out.write('{}\n\n'.format(well0.natom))
    for i, at in enumerate(well0.atom):
        x, y, z = well0.geom[i]
        geom_out.write('{} {:.6f} {:.6f} {:.6f}\n'.format(at, x, y, z))
    geom_out.write('\n\n')
    geom_out.close()

    # characterize the initial reactant
    well0.characterize(par.par['dimer'])
    well0.name = str(well0.chemid)
    start_name = well0.name

    # initialize the qc instance
    qc = QuantumChemistry(par)

    # start the initial optimization of the reactant
    logging.info('Starting optimization of intial well')
    qc.qc_opt(well0, well0.geom)
    err, well0.geom = qc.get_qc_geom(str(well0.chemid) + '_well',
                                     well0.natom, wait=1)
    err, well0.freq = qc.get_qc_freq(str(well0.chemid) + '_well',
                                     well0.natom, wait=1)
    if err < 0:
        logging.error('Error with initial structure optimization.')
        return
    if any(well0.freq[i] <= 0 for i in range(len(well0.freq))):
        logging.error('Found imaginary frequency for initial structure.')
        return

    # characterize again and look for differences
    well0.characterize(par.par['dimer'])
    well0.name = str(well0.chemid)
    if well0.name != start_name:
        logging.error('The first well optimized to a structure different from the input.')
        return

    # do an MP2 optimization of the reactant,
    # to compare Beta scission barrier heigths to
    logging.info('Starting MP2 optimization of intial well')
    qc.qc_opt(well0, well0.geom, mp2=1)
    err, geom = qc.get_qc_geom(str(well0.chemid) + '_well_mp2', well0.natom, 1)

    # characterize again and look for differences
    well0.characterize(par.par['dimer'])
    well0.name = str(well0.chemid)

    # read the energy and the zpe corrected energy
    err, well0.energy = qc.get_qc_energy(str(well0.chemid) + '_well', 1)
    err, well0.zpe = qc.get_qc_zpe(str(well0.chemid) + '_well', 1)

    well_opt = Optimize(well0, par, qc, wait=1)
    well_opt.do_optimization()
    if well_opt.shigh == -999:
        logging.error('Error with high level optimization of initial structure.')
        return

    # do the reaction search using heuristics
    if par.par['reaction_search'] == 1:
        logging.info('Starting reaction searches of intial well')
        rf = ReactionFinder(well0, par, qc)
        rf.find_reactions()
        rg = ReactionGenerator(well0, par, qc)
        rg.generate()
    # do the homolytic scission products search
    if par.par['homolytic_scissions'] == 1:
        logging.info('Starting the search for homolytic scission products')
        well0.homolytic_scissions = HomolyticScissions(well0, par, qc)
        well0.homolytic_scissions.find_homolytic_scissions()
    # initialize the master equation instance
    mess = MESS(par, well0)
    mess.write_input()
    mesmer = MESMER(par, well0)
    mesmer.write_input()
    if par.par['me'] == 1:
        logging.info('Starting Master Equation calculations')
        if par.par['me_code'] == 'mess':
            mess.run()
        elif par.par['me_code'] == 'mesmer':
            mesmer.run()
        else:
            logging.error('Cannot recognize me code {}'.format(par.par['me_code']))

    # postprocess the calculations
    postprocess.createSummaryFile(well0, qc, par)
    postprocess.createPESViewerInput(well0, qc, par)
    postprocess.creatMLInput(well0, qc, par)

    logging.info('Finished KinBot at {}'.format(datetime.datetime.now()))
    print("Done!")


if __name__ == "__main__":
    main()
