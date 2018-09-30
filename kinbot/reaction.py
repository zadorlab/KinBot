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
# -*- coding: utf-8 -*-
import numpy as np
import sys
import copy
import time
import logging

from geom import *
from vector import *
from qc import *
from constants import *
from kinbot import *
from stationary_pt import *
from reac_intra_H_migration import *
from reac_intra_R_migration import *
from reac_cpd_H_migration import *
from reac_intra_OH_migration import *
from reac_Intra_RH_Add_Endocyclic_F import *
from reac_Intra_RH_Add_Endocyclic_R import *
from reac_Cyclic_Ether_Formation import *
from reac_Intra_RH_Add_Exocyclic_F import *
from reac_Intra_RH_Add_Exocyclic_R import *
from reac_Intra_R_Add_Endocyclic_F import *
from reac_Intra_R_Add_ExoTetCyclic_F import *
from reac_Intra_R_Add_Exocyclic_F import *
from reac_Retro_Ene import *
from reac_Korcek_step2 import *
from reac_r22_cycloaddition import *
from reac_r12_cycloaddition import *
from reac_r12_insertion_R import *
from reac_r13_insertion_CO2 import *
from reac_r13_insertion_ROR import *
from reac_r14_birad_scission import *
from reac_r14_cyclic_birad_scission_R import *
from reac_birad_recombination_F import *
from reac_birad_recombination_R import *
from reac_Intra_disproportionation_F import *
from reac_Intra_disproportionation_R import *
from reac_Diels_alder_addition import *
from reac_Intra_Diels_alder_R import *
from reac_ketoenol import *
from reac_HO2_Elimination_from_PeroxyRadical import *
from reac_R_Addition_COm3_R import *
from reac_R_Addition_MultipleBond import *
from reac_12_shift_S_F import *
from reac_12_shift_S_R import *
from reac_r13_insertion_RSR import *
from reac_R_Addition_CSm_R import *
from reac_family import *
from irc import *
from postprocess import *
from zmat import *
from pes import *
from hir import *
import par



def reac_generator(species):
    """ 
    Creates the input for each reaction, runs them, and tests for success.
    If successful, it creates the barrier and product objects.
    It also then does the conformational search, and finally, the hindered rotor scans.
    To make the code the most efficient, all of these happen in parallel, in a sense that
    the jobs are not waiting for each other. E.g., one reaction can still be in the stage
    of TS search, while the other can be already at the hindered rotor scan. This way, 
    all cores are occupied efficiently.

    The switching between the various stages are done via the reac_ts_done variable.
    0: initiates the TS search
    1: checks barrier height and errors in TS, and initiates normal mode displacement test, start the irc calculations 
    2: submits product optimization
    3: submit the frequency calculation 
    4: reads frequencies, and if there is not error, writes MESS block
    5: do a conformational analysis of the TS
    6: do a high level calculation of the TS
    7: do the hindered rotor calculations of the TS
    8: do a conformational analysis of the products
    9: do a high level calculation of the products
    10: do the hindered rotor calculations of the products
    11: final calculations:
        * zmats of ts and products
        * symmetry numbers of ts and products
        * frequencies without internal rotations of ts and products

    """

    barriers = [stationary_pt(str(i)) for i in range(len(species.reac_inst))]
    products = [stationary_pt(str(i)) for i in range(len(species.reac_inst))]

    #short names for the MESS input files
    well_short_names = {}
    well_short_names[str(species.chemid)] = species.short_name
    
    bimol_short_names = {}
    fragment_short_names = {}
    ts_short_names = {}
    
    #product_bonds keeps track of the bond matrix of bimolecular products before
    #they get divided in separate products. This is used to build the bond matrix
    #of the transition state
    product_bonds = [[] for i in range(len(species.reac_inst))]
    
    #keep track of the conformational search and hindered rotor potentials of the
    #transition states and products
    confs = [-1 for i in range(len(species.reac_inst))]
    hirs = [-1 for i in range(len(species.reac_inst))]
    hir_restart = [0 for i in range(len(species.reac_inst))]
    prod_confs = [[] for i in range(len(species.reac_inst))]
    prod_hir = [[] for i in range(len(species.reac_inst))]
    prod_hir_restart = [0 for i in range(len(species.reac_inst))]
    
    if len(species.reac_inst) > 0:
        alldone = 1
    else: 
        alldone = 0

    while alldone:
        for index, instance in enumerate(species.reac_inst):
            # START REATION SEARCH
            if species.reac_ts_done[index] == 0 and species.reac_step[index] == 0:
                #verify after restart if search has failed in previous kinbot run
                if check_qc(species.reac_name[index]) == 'error' or check_qc(species.reac_name[index]) == 'killed':
                    logging.info('\tRxn search failed (error or killed) for %s'%species.reac_name[index])
                    species.reac_ts_done[index] = -999

            
            if species.reac_ts_done[index] == 0: # ts search is ongoing
                
                if species.reac_type[index] == 'combinatorial': 
                    if species.reac_step[index] == 47:
                        if read_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif read_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_combinatorial(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'cpd_H_migration': 
                    if species.reac_step[index] == 3:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_cpd_H_migration(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'intra_H_migration':
                    if species.reac_step[index] == 15:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_intra_H_migration(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'intra_R_migration':
                    if species.reac_step[index] == 15:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_intra_R_migration(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'intra_OH_migration':       
                    if species.reac_step[index] == 15:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_intra_OH_migration(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'Intra_RH_Add_Endocyclic_F':
                    if species.reac_step[index] == 23:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_Intra_RH_Add_Endocyclic_F(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'Intra_RH_Add_Endocyclic_R':
                    if species.reac_step[index] == 13:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_Intra_RH_Add_Endocyclic_R(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'Intra_R_Add_Endocyclic_F':
                    if species.reac_step[index] == 15:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_Intra_R_Add_Endocyclic_F(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'Intra_R_Add_Endocyclic_R':
                    if species.reac_step[index] == 15:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_Intra_R_Add_Endocyclic_R(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'Intra_R_Add_ExoTetCyclic_F':
                    if species.reac_step[index] == 23:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_Intra_R_Add_ExoTetCyclic_F(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'Intra_R_Add_Exocyclic_F':
                    if species.reac_step[index] == 15:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_Intra_R_Add_Exocyclic_F(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'Cyclic_Ether_Formation':       
                    if species.reac_step[index] == 15:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_Cyclic_Ether_Formation(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'Intra_RH_Add_Exocyclic_F':       
                    if species.reac_step[index] == 23:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_Intra_RH_Add_Exocyclic_F(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'Intra_RH_Add_Exocyclic_R':
                    if species.reac_step[index] == 13:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_Intra_RH_Add_Exocyclic_R(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'Retro_Ene':       
                    if species.reac_step[index] == 23:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_Retro_Ene(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'Korcek_step2':       
                    if species.reac_step[index] == 13:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_Korcek_step2(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'r22_cycloaddition':       
                    if species.reac_step[index] == 2:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_r22_cycloaddition(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'r12_cycloaddition':       
                    if species.reac_step[index] == par.scan_step + 1:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else:        
                        if species.reac_step[index] == 0:
                            species.reac_step[index] = do_r12_cycloaddition(species, instance, species.reac_step[index], species.reac_name[index])
                        elif species.reac_step[index] > 0:
                            err, energy = read_qc_energy(species.reac_name[index])
                            if err == 0:
                                species.reac_scan_energy[index].append(energy)
                                if len(species.reac_scan_energy[index]) > 1:
                                    if species.reac_scan_energy[index][-1] < species.reac_scan_energy[index][-2]:
                                        species.reac_step[index] = par.scan_step 
                                species.reac_step[index] = do_r12_cycloaddition(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'r12_insertion_R':       
                    if species.reac_step[index] == 13:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_r12_insertion_R(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'r13_insertion_CO2':       
                    if species.reac_step[index] == 23:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_r13_insertion_CO2(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'r13_insertion_ROR':       
                    if species.reac_step[index] == 23:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_r13_insertion_ROR(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'r14_birad_scission':       
                    if species.reac_step[index] == par.scan_step + 1:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else:        
                        if species.reac_step[index] == 0:
                            species.reac_step[index] = do_r14_birad_scission(species, instance, species.reac_step[index], species.reac_name[index])
                        elif species.reac_step[index] > 0:
                            err, energy = read_qc_energy(species.reac_name[index])
                            if err == 0:
                                species.reac_scan_energy[index].append(energy)
                                if len(species.reac_scan_energy[index]) > 1:
                                    if species.reac_scan_energy[index][-1] < species.reac_scan_energy[index][-2]:
                                        species.reac_step[index] = par.scan_step 
                                species.reac_step[index] = do_r14_birad_scission(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'r14_cyclic_birad_scission_R':       
                    if species.reac_step[index] == 23:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_r14_cyclic_birad_scission_R(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'birad_recombination_F':       
                    if species.reac_step[index] == 23:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_birad_recombination_F(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'birad_recombination_R':       
                    if species.reac_step[index] == par.scan_step + 1:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else:        
                        if species.reac_step[index] == 0:
                            species.reac_step[index] = do_birad_recombination_R(species, instance, species.reac_step[index], species.reac_name[index])
                        elif species.reac_step[index] > 0:
                            err, energy = read_qc_energy(species.reac_name[index])
                            if err == 0:
                                species.reac_scan_energy[index].append(energy)
                                if len(species.reac_scan_energy[index]) > 1:
                                    if species.reac_scan_energy[index][-1] < species.reac_scan_energy[index][-2]:
                                        species.reac_step[index] = par.scan_step 
                                species.reac_step[index] = do_birad_recombination_R(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'Intra_disproportionation_F':       
                    if species.reac_step[index] == 23:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_Intra_disproportionation_F(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'Intra_disproportionation_R':       
                    if species.reac_step[index] == 23:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_Intra_disproportionation_R(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'Diels_alder_addition':       
                    if species.reac_step[index] == 2:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_Diels_alder_addition(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'Intra_Diels_alder_R':       
                    if species.reac_step[index] == 23:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_Intra_Diels_alder_R(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'ketoenol':                     
                    if species.reac_step[index] == 15:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_ketoenol(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'H2_Elimination':
                    if species.reac_step[index] == 15:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_H2_Elimination(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'HO2_Elimination_from_PeroxyRadical':                    
                    if species.reac_step[index] == 15:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_HO2_Elimination_from_PeroxyRadical(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'R_Addition_COm3_R':       
                    if species.reac_step[index] == 2:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_R_Addition_COm3_R(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'R_Addition_MultipleBond':       
                    if species.reac_step[index] == 13:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_R_Addition_MultipleBond(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == '12_shift_S_F':       
                    if species.reac_step[index] == 13:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_12_shift_S_F(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == '12_shift_S_R':       
                    if species.reac_step[index] == 13:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_12_shift_S_R(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'r13_insertion_RSR':       
                    if species.reac_step[index] == 23:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_r13_insertion_RSR(species, instance, species.reac_step[index], species.reac_name[index])

                if species.reac_type[index] == 'R_Addition_CSm_R':       
                    if species.reac_step[index] == 2:
                        if get_qc_freq(species.reac_name[index], par.natom)[0] == 0:
                            species.reac_ts_done[index] = 1
                        elif get_qc_freq(species.reac_name[index], par.natom)[0] == -1:
                            logging.info('\tRxn search failed for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else: 
                        species.reac_step[index] = do_R_Addition_CSm_R(species, instance, species.reac_step[index], species.reac_name[index])

            elif species.reac_ts_done[index] == 1:
                if check_qc(species.reac_name[index]) == 'running': continue
                elif check_qc(species.reac_name[index]) == 'error': 
                    logging.info('\tRxn search failed (gaussian error) for %s'%species.reac_name[index])
                    species.reac_ts_done[index] = -999
                else: 
                    #check the barrier height:
                    if species.reac_type[index] == 'R_Addition_MultipleBond':
                        sp_energy = get_qc_energy(str(species.chemid) + '_well_mp2')[1]
                        barrier = (get_qc_energy(species.reac_name[index])[1] - sp_energy) * AUtoKCAL
                    else:
                        barrier = (get_qc_energy(species.reac_name[index])[1] - species.energy) * AUtoKCAL
                    if barrier > par.barrier_threshold:
                        logging.info('\tRxn barrier too high for %s'%species.reac_name[index])
                        species.reac_ts_done[index] = -999
                    else:
                        #test_normal_mode_displacement(species, index, barriers, products)
                        irc_status = check_irc(species,index)
                        if 0 in irc_status:
                            # No IRC started yet, start the IRC now
                            logging.info('\tStarting IRC calculations for %s'%species.reac_name[index])
                            do_irc_calculations(species,index)
                        elif irc_status[0] == 'running' or irc_status[1] == 'running':
                            continue
                        else: 
                            #IRC's have succesfully finished, have an error or were killed, in any case
                            #read the geometries and try to make products out of them
                            #verify which of the ircs leads back to the reactant, if any
                            prod = irc2stationary_pt(species,index)
                            if prod == 0:
                                logging.info('\t\tNo product found for %s'%species.reac_name[index])
                                species.reac_ts_done[index] = -999
                            else:
                                #IRC's are done, continue with the optimization of the products
                                products[index] = prod
                                product_bonds[index] = prod.bond
                                #qc_opt(products[index], products[index].geom, 0, par.natom, par.atom, par.mult, par.charge)
                                species.reac_ts_done[index] = 2
            elif species.reac_ts_done[index] == 2:
                #identify bimolecular products and wells
                fragments =products[index].start_multi_molecular(par.natom,par.atom)
                products[index] = []
                for i,frag in enumerate(fragments):
                    products[index].append(frag)
                    qc_opt(frag, frag.geom, 0, frag.natom, frag.atom, frag.mult, frag.charge)
                if par.pes:
                    #verify if product is monomolecular, and if it is new
                    if len(fragments) ==1:
                        chemid = products[index][0].chemid
                        dir = os.path.dirname(os.getcwd()) 
                        jobs = open(dir+'/chemids','r').read().split('\n')
                        
                        if not str(chemid) in jobs:
                            write_input(par,products[index][0],dir)
                            f = open(dir+'/chemids','a')
                            f.write(str(chemid)+'\n')
                            f.close()
                species.reac_ts_done[index] = 3
            elif species.reac_ts_done[index] == 3:
                #wait for the optimization to finish and optionally start freq calculations
                err = 0
                for st_pt in products[index]:
                    chemid = st_pt.chemid
                    orig_geom = copy.deepcopy(st_pt.geom)
                    e, st_pt.geom = get_qc_geom(str(st_pt.chemid) + '_well', st_pt.natom)
                    e2, st_pt.energy = get_qc_energy(str(st_pt.chemid) + '_well')
                    if e < 0:
                        logging.info('\tProduct optimization failed for %s, product %s'%(species.reac_name[index],st_pt.chemid))
                        species.reac_ts_done[index] = -999
                        err = -1
                    elif e != 0:
                        err = -1
                    else:
                        st_pt.bond_mx(st_pt.natom, st_pt.atom)
                        st_pt.characterize(st_pt.natom, st_pt.atom, st_pt.mult, st_pt.charge)
                        st_pt.calc_chemid(st_pt.natom,st_pt.atom,st_pt.mult)
                        if chemid != st_pt.chemid:
                            #product was optimized to another structure, give warning and remove this reaction
                            logging.info('\tProduct optimizatied to other structure for %s, product %s to %s'%(species.reac_name[index],chemid,st_pt.chemid))
                            species.reac_ts_done[index] = -999
                            err = -1
                if err == 0:
                    for p,st_pt in enumerate(products[index]):
                        if par.high_level == 0:
                            # only do this frequency calculation is no high-level calculations are done
                            qc_freq(st_pt, st_pt.geom, st_pt.natom, st_pt.atom, st_pt.mult, st_pt.charge)
                    species.reac_ts_done[index] = 4
            elif species.reac_ts_done[index] == 4:
                err = 0
                for st_pt in products[index]:
                    if par.high_level == 0:
                        e, st_pt.freq = get_qc_freq(str(st_pt.chemid) + '_fr', st_pt.natom)
                        e2, st_pt.zpe = get_qc_zpe(str(st_pt.chemid) + '_fr')
                        if any(st_pt.freq[i] < 0 for i in range(len(st_pt.freq))) or e < 0:
                            logging.info('\tProduct frequency calculations failed for %s, product %s'%(species.reac_name[index],st_pt.chemid))
                            species.reac_ts_done[index] = -999
                            err = -1
                        elif e != 0:
                            e, st_pt.zpe = get_qc_zpe(str(st_pt.chemid) + '_fr')
                            err = -1
                if err == 0:
                    logging.info('\tSuccessfully finished reaction %s'%species.reac_name[index])
                    
                    #make a stationary point object of the ts
                    bond_mx = np.zeros((par.natom, par.natom), dtype=int)
                    for i in range(par.natom):
                        for j in range(par.natom):
                            bond_mx[i][j] = max(species.bond[i][j],product_bonds[index][i][j])
                    ts = stationary_pt(species.reac_name[index])
                    err, ts.geom = get_qc_geom(species.reac_name[index], par.natom)
                    err, ts.energy = get_qc_energy(species.reac_name[index])
                    ts.bond = bond_mx
                    ts.find_cycle(par.natom,par.atom)
                    ts.find_conf_dihedral(par.natom,par.atom,par.mult)
                    barriers[index] = ts
                
                    species.reac_ts_done[index] = 5 
                    for st_pt in products[index]:
                        st_pt.find_cycle(st_pt.natom,st_pt.atom) #needed for the conformational search later on
                        st_pt.find_conf_dihedral(st_pt.natom,st_pt.atom,st_pt.mult) #needed for the conformational search later on
            elif species.reac_ts_done[index] == 5:
                # do a conformational search of the ts
                ts = barriers[index]
                if par.conformer_search == 1:
                    if confs[index] == -1:
                        #generate cyclic conformers, if any
                        if len(ts.cycle_chain) > 0:
                            geoms = generate_ring_conformers(ts, par.natom, par.atom, par.mult, par.charge, copy.deepcopy(ts.geom), 1)
                        else:
                            geoms = [copy.deepcopy(ts.geom)]
                        conf = 0
                        for geom in geoms:
                            conf = generate_conformers(ts, par.natom, par.atom, par.mult, par.charge, 0, copy.deepcopy(ts.geom), conf, 1)
                        confs[index] = conf
                    else:
                        status, geom = check_conformers(ts, confs[index], par.natom, par.atom, par.mult, par.charge,wait = 0, ts = 1)
                        if status == 1:
                            ts.geom = geom
                            species.reac_ts_done[index] = 6
                else:
                    species.reac_ts_done[index] = 6
                    
            elif species.reac_ts_done[index] == 6:
                #do a high level calculation of the ts
                if par.high_level == 1:
                    job = species.reac_name[index] + '_high'
                    status = check_qc(job)
                    
                    if status == 0:
                        kwargs = get_qc_arguments(job,par.mult,ts = 1,step = 1,max_step=1,high_level = 1)
                        carry_out_reaction(species, instance, 1, job, 1, barriers[index].geom, kwargs, [], [], [])
                    elif status == 'running':
                        continue
                    elif status == 'error':
                        logging.info('\tHigh level optimization failed for %s'%species.reac_name[index])
                        species.reac_ts_done[index] = -999
                    elif status == 'normal':
                        ts = barriers[index]
                        err, new_geom = get_qc_geom(job, par.natom)
                        if equal_geom(ts.bond,ts.geom,new_geom,0.1):
                            err, ts.geom = get_qc_geom(job, par.natom)
                            err, ts.freq = get_qc_freq(job, par.natom)
                            err, ts.zpe = get_qc_zpe(job)
                            err, ts.energy = get_qc_energy(job)
                            species.reac_ts_done[index] = 7
                        else:
                            logging.info('\tHigh level ts optimization converged to different structure for %s'%species.reac_name[index])
                            species.reac_ts_done[index] = -999
                    else:
                        logging.info('\tDont know what happend with %s'%species.reac_name[index])
                        species.reac_ts_done[index] = -999
                else:
                    species.reac_ts_done[index] = 7
            elif species.reac_ts_done[index] == 7:
                # do the hindered rotor calculations of the TS
                if par.rotor_scan == 1:
                    ts = barriers[index]
                    if hirs[index] == -1:
                        generate_hir_geoms(ts, par.natom, par.atom, par.mult, par.charge, copy.deepcopy(ts.geom), 1)
                        hirs[index] = 1
                    else:
                        status = check_hir(ts,par.natom, par.atom, par.mult, par.charge, 1,wait = 0)
                        if status == 1:
                            if len(ts.hir_energies) > 0:
                                #check if along the hir potential a structure was found with a lower energy
                                min = ts.hir_energies[0][0]
                                min_rotor = -1
                                min_ai = -1
                                for rotor in range(len(ts.dihed)):
                                    for ai in range(par.nrotation):
                                        if ts.hir_energies[rotor][ai] < min - 1.6E-4: #use a 0.1kcal/mol cutoff for numerical noise 
                                            min = ts.hir_energies[rotor][ai]
                                            min_rotor = rotor
                                            min_ai = ai
                                if min_rotor > -1 and hir_restart[index] < 3:
                                    #lower energy structure found
                                    logging.info("Lower energy found during hindered rotor scan for " + species.reac_name[index])
                                    logging.info("Rotor: " + str(min_rotor))
                                    logging.info("Scan point: " + str(min_ai))
                                    job = 'hir/' + species.reac_name[index] + '_hir_' + str(min_rotor) + '_' + str(min_ai).zfill(2)
                                    err,ts.geom = get_qc_geom(job, par.natom)
                                    #delete the high_level log file and the hir log files
                                    if os.path.exists(species.reac_name[index] + '_high.log'):
                                        os.remove(species.reac_name[index] + '_high.log') 
                                    for rotor in range(len(ts.dihed)):
                                        for ai in range(par.nrotation):
                                            if os.path.exists('hir/' + species.reac_name[index] + '_hir_' + str(rotor) + '_' + str(ai).zfill(2) + '.log'):
                                                os.remove('hir/' + species.reac_name[index] + '_hir_' + str(rotor) + '_' + str(ai).zfill(2) + '.log') 
                                    hirs[index] = -1
                                    hir_restart[index] += 1
                                    species.reac_ts_done[index] = 6
                                else:
                                    species.reac_ts_done[index] = 8
                            else:
                                species.reac_ts_done[index] = 8
                else:
                    species.reac_ts_done[index] = 8
            elif species.reac_ts_done[index] == 8:
                # do a conformational search of the products
                if par.conformer_search == 1:
                    tot_status = [-1 for i in range(len(products[index]))]
                    for p,st_pt in enumerate(products[index]):
                        if len(prod_confs[index]) <= p:
                            #generate cyclic conformers, if any
                            if len(st_pt.cycle_chain) > 0:
                                geoms = generate_ring_conformers(st_pt, st_pt.natom, st_pt.atom, st_pt.mult, st_pt.charge, copy.deepcopy(st_pt.geom), 0)
                            else:
                                geoms = [copy.deepcopy(st_pt.geom)]
                            conf = 0
                            for geom in geoms:
                                conf = generate_conformers(st_pt, st_pt.natom, st_pt.atom, st_pt.mult, st_pt.charge, 0, copy.deepcopy(st_pt.geom), 0, 0)
                            confs[index] = conf
                            prod_confs[index].append(conf)
                        else:
                            status, geom = check_conformers(st_pt, prod_confs[index][p], st_pt.natom, st_pt.atom, st_pt.mult, st_pt.charge,wait = 0)
                            tot_status[p] = status
                            if status == 1:
                                st_pt.geom = geom
                    
                    if all([si == 1 for si in tot_status]):
                        species.reac_ts_done[index] = 9
                else:
                    species.reac_ts_done[index] = 9
            elif species.reac_ts_done[index] == 9:
                #do a high level calculation of the products
                if par.high_level == 1:
                    opt_stats = [check_qc(str(st_pt.chemid) + '_well_high') for st_pt in products[index]]
                    fr_stats = [check_qc(str(st_pt.chemid) + '_fr_high') for st_pt in products[index]]
                    for i, st_pt in enumerate(products[index]):
                        opt_status = opt_stats[i]
                        fr_status = fr_stats[i]
                        if opt_status == 0:
                            qc_opt(st_pt, st_pt.geom, 0, st_pt.natom, st_pt.atom, st_pt.mult, st_pt.charge, high_level = 1)
                        elif opt_status == 'running':
                            continue
                        elif opt_status == 'error':
                            logging.info('\tHigh level optimization failed for %s, product %s'%(species.reac_name[index],st_pt.chemid))
                            species.reac_ts_done[index] = -999
                        elif opt_status == 'normal':
                            err, st_pt.geom = get_qc_geom(str(st_pt.chemid) + '_well_high', par.natom)
                            err, st_pt.energy = get_qc_energy(str(st_pt.chemid) + '_well_high')
                            if fr_status == 0:
                                qc_freq(st_pt, st_pt.geom, st_pt.natom, st_pt.atom, st_pt.mult, st_pt.charge, high_level = 1)
                            elif fr_status == 'running':
                                continue
                            elif fr_status == 'error':
                                logging.info('\tHigh level frequency calculation failed for %s, product %s'%(species.reac_name[index],st_pt.chemid))
                                species.reac_ts_done[index] = -999
                            elif fr_status == 'normal':
                                #read the frequency output
                                err, st_pt.freq = get_qc_freq(str(st_pt.chemid) + '_fr_high', st_pt.natom)
                                err, st_pt.zpe = get_qc_zpe(str(st_pt.chemid) + '_fr_high')
                            else:
                                logging.info('\tDont know what happend with %s, product %s'%(species.reac_name[index],st_pt.chemid))
                                species.reac_ts_done[index] = -999
                        else:
                            logging.info('\tDont know what happend with %s, product %s'%(species.reac_name[index],st_pt.chemid))
                            species.reac_ts_done[index] = -999
                    if all([fi == 'normal' for fi in fr_stats]):
                        species.reac_ts_done[index] = 10
                else:
                    species.reac_ts_done[index] = 10
            elif species.reac_ts_done[index] == 10:
                # do the hindered rotor calculations of the products
                
                if par.rotor_scan == 1:
                    tot_status = [-1 for i in range(len(products[index]))]
                    for p, st_pt in enumerate(products[index]):
                        if len(prod_hir[index]) <= p:
                            generate_hir_geoms(st_pt, st_pt.natom, st_pt.atom, st_pt.mult, st_pt.charge, copy.deepcopy(st_pt.geom), 0)
                            prod_hir[index].append(1)
                        else:
                            status = check_hir(st_pt,st_pt.natom, st_pt.atom, st_pt.mult, st_pt.charge, 0,wait = 0)
                            tot_status[p] = status
                    if all([si == 1 for si in tot_status]):
                        lower = 0
                        #check for lower energy structures along the hir scans
                        for p, st_pt in enumerate(products[index]):
                            if len(st_pt.hir_energies) > 0:
                                #check if along the hir potential a structure was found with a lower energy
                                min = st_pt.hir_energies[0][0]
                                min_rotor = -1
                                min_ai = -1
                                for rotor in range(len(st_pt.dihed)):
                                    for ai in range(par.nrotation):
                                        if st_pt.hir_energies[rotor][ai] < min - 1.6E-4: #use a 0.1kcal/mol cutoff for numerical noise 
                                            min = st_pt.hir_energies[rotor][ai]
                                            min_rotor = rotor
                                            min_ai = ai
                                
                                if min_rotor > -1 and prod_hir_restart[index] < 3:
                                    #lower energy structure found
                                    lower = 1
                                    logging.info("Lower energy found during hindered rotor scan for " + str(st_pt.chemid))
                                    job = 'hir/' + str(st_pt.chemid) + '_hir_' + str(min_rotor) + '_' + str(min_ai).zfill(2)
                                    err,st_pt.geom = get_qc_geom(job, st_pt.natom)
                                    #delete the high_level log file and the hir log files
                                    if os.path.exists(str(st_pt.chemid) + '_well_high.log'):
                                        os.remove(str(st_pt.chemid) + '_well_high.log') 
                                    if os.path.exists(str(st_pt.chemid) + '_fr_high.log'):
                                        os.remove(str(st_pt.chemid) + '_fr_high.log') 
                                    for rotor in range(len(st_pt.dihed)):
                                        for ai in range(par.nrotation):
                                            if os.path.exists('hir/' + str(st_pt.chemid) + '_hir_' + str(rotor) + '_' + str(ai).zfill(2) + '.log'):
                                                os.remove('hir/' + str(st_pt.chemid) + '_hir_' + str(rotor) + '_' + str(ai).zfill(2) + '.log')
                        if lower:
                            prod_hir[index] = []
                            prod_hir_restart[index] += 1
                            species.reac_ts_done[index] = 9
                        else:
                            species.reac_ts_done[index] = 11
                else:
                    species.reac_ts_done[index] = 11
            elif species.reac_ts_done[index] == 11:
                # write the zmatrices
                ts = barriers[index]
                if len(ts.dihed) >= 0:
                    estoktp_zmat(ts, species.reac_name[index] , par.natom, par.atom)
                
                #calculate the symmetry numbers
                calculate_symmetry(ts, par.natom, par.atom)
                
                #calculate the new frequencies by projecting out the internal rotational ones
                if par.high_level == 1:
                    fr_file = species.reac_name[index] + '_high'
                else:
                    fr_file = species.reac_name[index]
                hess = read_qc_hess(fr_file,par.natom)
                ts.kinbot_freqs, ts.reduced_freqs = get_frequencies(ts, hess, par.natom, par.atom, ts.geom)
                
                #write the molpro input and read the molpro energy, if available
                create_molpro_input(ts, par.natom, par.atom, par.mult, par.charge, 1)
                status, molpro_energy = get_molpro_energy(ts, 1)
                if status:
                    ts.energy = molpro_energy
                
                for st_pt in products[index]:
                    # write the zmatrices
                    if len(st_pt.dihed) >= 0:
                        estoktp_zmat(st_pt, str(st_pt.chemid) + '_well', st_pt.natom, st_pt.atom)
                    
                    #calculate the symmetry numbers
                    calculate_symmetry(st_pt, st_pt.natom, st_pt.atom)
                    
                    #calculate the new frequencies by projecting out the internal rotational ones
                    if par.high_level == 1:
                        fr_file = str(st_pt.chemid) + '_fr_high'
                    else:
                        fr_file = str(st_pt.chemid) + '_fr'
                    hess = read_qc_hess(fr_file,st_pt.natom)
                    st_pt.kinbot_freqs, st_pt.reduced_freqs = get_frequencies(st_pt, hess, st_pt.natom, st_pt.atom, st_pt.geom)
                    
                    #write the molpro input and read the molpro energy, if available
                    create_molpro_input(st_pt, st_pt.natom, st_pt.atom, st_pt.mult, st_pt.charge, 0)
                    status, molpro_energy = get_molpro_energy(st_pt, 0)
                    if status:
                        st_pt.energy = molpro_energy
                
                neg_freq = 0
                for st_pt in products[index]:
                    if any([fi < 0. for fi in st_pt.reduced_freqs]):
                        neg_freq = 1
                if any([fi < 0. for fi in ts.reduced_freqs[1:]]): #check for more than one neg freq
                    neg_freq = 1
                
                if neg_freq:
                    logging.info('\tFound negative frequency for ' + species.reac_name[index])
                    species.reac_ts_done[index] = -999
                else:
                    #generate names for the mess calculations
                    if len(products[index]) == 1:
                        if not str(products[index][0].chemid) in well_short_names:
                            short_name = 'w'+str(len(well_short_names)+1)
                            well_short_names[str(products[index][0].chemid)] = short_name
                        else:
                            short_name = well_short_names[str(products[index][0].chemid)]
                        products[index][0].short_name = short_name
                    else:
                        name = '_'.join(sorted([str(prod.chemid) for prod in products[index]]))
                        if not name in bimol_short_names:
                            bimol_short_names[name] = 'b'+str(len(bimol_short_names)+1)
                        for prod in products[index]:
                            if not str(prod.chemid) in fragment_short_names:
                                name = 'f' + str(len(fragment_short_names)+1)
                                fragment_short_names[str(prod.chemid)] = name
                            else:
                                name = fragment_short_names[str(prod.chemid)]
                            prod.short_name = name
                    ts_short_names[species.reac_name[index]] = 'ts' + str(len(ts_short_names)+1)
                    
                    #write the mess block for the ts
                    write_mess_barrier(species, index, ts, products[index], par.natom, par.atom, par.mult, par.charge, well_short_names, bimol_short_names, ts_short_names)
                    
                    #write the mess block for the products
                    if len(products[index]) > 1:
                        write_mess_bimol(products[index],species, bimol_short_names, fragment_short_names)
                    else:
                        st_pt = products[index][0]
                        write_mess_well(st_pt, st_pt.natom, st_pt.atom, st_pt.mult, st_pt.charge, species, well_short_names)

                    #the reaction search is finished
                    species.reac_ts_done[index] = -1 # this is the success code

        
        alldone = 1
        for index, instance in enumerate(species.reac_inst):
            if any(species.reac_ts_done[i] >= 0 for i in range(len(species.reac_inst))):
                alldone = 1
                break 
            else: 
                alldone = 0
        
        
        # write a small summary while running
        wr = 1
        if wr:
            f_out = open('kinbot_monitor.out','w')
            for index, instance in enumerate(species.reac_inst):
                f_out.write('{}\t{}\t{}\n'.format(species.reac_ts_done[index],species.reac_step[index],species.reac_name[index]))
            f_out.close()
        time.sleep(1)
    
    if par.me:
        #write a complete mess file
        write_mess_all(species, barriers, products, well_short_names, bimol_short_names, fragment_short_names, ts_short_names)
        if not os.path.exists('mess/all.out'):
            run_mess()
        all_rates = read_mess()
        write_high_p_rates(all_rates, well_short_names,bimol_short_names)
    
    createSummaryFile(species,barriers,products)
    createPESViewerInput(species,barriers,products)
    logging.info("Done!")

def test_normal_mode_displacement(species, index, barriers, products):
    """
    Test for each reaction type whether the normal mode displacement 
    suggests that it is the right chemical reaction.
    """
    
    err, tsgeom = get_qc_geom(species.reac_name[index], par.natom)
    err, imode = get_qc_geom(species.reac_name[index], par.natom) 
    
    if err != 0:
        species.reac_ts_done[index] = -999
        return -1

    if species.reac_type[index] == 'intra_H_migration':
        rad = species.reac_inst[index][0]
        h = species.reac_inst[index][-1]
        site = species.reac_inst[index][-2]
        
        test_new_well = stationary_pt('test_new_well')
        test_new_well.geom = imode_shift(rad, h, par.atom[rad], par.atom[h], tsgeom, imode)
        if np.array_equal(test_new_well.geom, tsgeom): # nothing happened, because the equations didn have a solution
            species.reac_ts_done[index] = -999
            return -1

        test_new_well.bond_mx(par.natom, par.atom)
        bond_diff = test_new_well.bond - species.bond

        if bond_diff[h][rad] == 1 and bond_diff[h][site] == -1 and np.count_nonzero(bond_diff) == 4:
            test_old_well = stationary_pt('test_old_well')
            test_old_well.geom = imode_shift(site, h, par.atom[site], par.atom[h], tsgeom, imode)
            if np.array_equal(test_old_well.geom, tsgeom):
                species.reac_ts_done[index] = -999
                return -1
        
            test_old_well.bond_mx(par.natom, par.atom)
            bond_diff = test_old_well.bond - species.bond
            if np.count_nonzero(bond_diff) == 0:
                
                barriers[index] = stationary_pt('barrier' + str(index))
                barriers[index].geom = copy.copy(tsgeom)
                
                products[index] = stationary_pt('product' + str(index))
                products[index].geom = copy.copy(test_new_well.geom)
                
                species.reac_ts_done[index] = 2
                return 0
                
    species.reac_ts_done[index] = -999
    return -1
        


def imode_shift(atom1, atom2, element1, element2, tsgeom, imode):
    """
    Calculate the geometry newgeom relative to tsgeom that we get if
    the atoms are shifted by imode normal mode in a way that atom1 and atom2
    are at their canonical distance.   
    """            


    A = 0.
    B = 0.
    C = 0.
    for i in range(3):
        A += (imode[atom1][i] - imode[atom2][i])**2.
        B += (tsgeom[atom1][i] - tsgeom[atom2][i]) * (imode[atom1][i] - imode[atom2][i])
        C += (tsgeom[atom1][i] - tsgeom[atom2][i])**2.

    B = 2. * B
    
    desired_dist = st_bond[''.join(sorted(element1 + element2))] / 1.2

    min_dist_geom = tsgeom - B / (2. * A) * imode
    min_dist = np.linalg.norm(min_dist_geom[atom1] - min_dist_geom[atom2])

    if min_dist > desired_dist:
        return tsgeom

    C = C - desired_dist**2
                    
    lamb1 = (-B + np.sqrt(B**2. - 4. * A * C)) / (2. * A)
    lamb2 = (-B - np.sqrt(B**2. - 4. * A * C)) / (2. * A)
    lamb = min(abs(lamb1), abs(lamb2)) * np.sign(lamb1) 
    
    return tsgeom + imode * lamb



if __name__ == "__main__":
    main()
