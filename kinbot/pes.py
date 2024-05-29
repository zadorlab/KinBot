"""
This is the main class to run KinBot to explore
a full PES instead of only the reactions of one well
"""
import sys
import os
import re
import stat
import shutil
import datetime
import time
import subprocess
import json
import networkx as nx
import numpy as np
import getpass

from copy import deepcopy
from ase.db import connect

from kinbot import kb_path
from kinbot import constants
from kinbot import license_message
from kinbot import pp_settings
from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint
from kinbot.fragments import Fragment
from kinbot.mess import MESS
from kinbot.uncertaintyAnalysis import UQ
from kinbot.config_log import config_log


def main():
    if sys.version_info.major < 3:
        print(f'KinBot only runs with python 3.10 or higher. You have python {sys.version_info.major}.{sys.version_info.minor}. Bye!')
        sys.exit(-1)
    elif sys.version_info.minor < 10:
        print(f'KinBot only runs with python 3.10 or higher. You have python {sys.version_info.major}.{sys.version_info.minor}. Bye!')
        sys.exit(-1)

    try:
        input_file = sys.argv[1]
    except IndexError:
        print('To use the pes script, supply one argument being the input file!')
        sys.exit(-1)

    # TODO: write information about the arguments
    # change this to nice argument parsers with
    # dashes etc.
    no_kinbot = 0
    task = 'all'
    names = []
    if len(sys.argv) > 2:
        if sys.argv[2] == 'no-kinbot':
            no_kinbot = 1
    if len(sys.argv) > 3:
        # possible tasks are:
        # 1. all: This is the default showing all pathways
        # 2. lowestpath: show the lowest path between the species
        # corresponding to the names
        # 3. allpaths: show all paths between the species
        # corresponding to the names
        # 4. wells: show all reactions of one wells
        # corresponding to the names
        task = sys.argv[3]
        names = sys.argv[4:]

    # print the license message to the console
    print(license_message.message)
    global logger
    logger = config_log('KinBot', mode='pes')

    # initialize the parameters
    par = Parameters(input_file, show_warnings=True).par

    # set up the logging environment
    if par['verbose']:
        logger = config_log('KinBot', mode='pes', level='debug')

    msg = 'Starting the PES search at {}'.format(datetime.datetime.now())
    logger.info(msg)

    well0 = StationaryPoint('well0',
                            par['charge'],
                            par['mult'],
                            smiles=par['smiles'],
                            structure=par['structure'])
    well0.characterize()
    write_input(input_file, well0, par['barrier_threshold'], par['barrier_threshold_L2'], os.getcwd(), par['me'])

    # add the initial well to the chemids
    with open('chemids', 'w') as f:
        f.write(str(well0.chemid) + '\n')

    # create a directory for the L3 single point calculations
    # directory has the name of the code, e.g., molpro
    try:
        os.mkdir(par['single_point_qc'])
    except OSError:
        pass

    # jobs that are running
    running = []
    # jobs that are finished
    finished = []
    # jobs that are waiting 
    waiting = []
    # list of all jobs
    jobs = []
    # dict of the pid's for all jobs
    pids = {}
    a = 0
    b = 0
    c = 0

    if 'none' not in par['keep_chemids']:
        with open('chemids', 'w') as f:
            for keep in par['keep_chemids']:
                f.write(keep + '\n')
                write_input_keep(input_file, keep, os.getcwd())

    if 'none' in par['skip_chemids']:
        logger.info('No KinBot runs to be skipped.')
    if 'none' in par['keep_chemids']:
        logger.info('All valid explored KinBot runs are kept.')

    while 1:
        j = len(jobs)
        if j != a:
            logger.info('{0} {1} {2}'.format("len(jobs): ", j, "\n"))
        a = j
        with open('chemids', 'r') as f:
            jobs = f.read().split('\n')
            jobs = [ji for ji in jobs if ji != '']

        if len(jobs) > j:
            logger.info('Picked up new jobs: ' + ' '.join(jobs[j:]))

        k = len(running)
        l = len(finished)
        if b != k:
            logger.info('{0} {1} {2}'.format("len(running): ", len(running), "\n"))
        b = k
        if c != l:
            logger.info('{0} {1} {2}'.format("len(finished): ", len(finished), "\n"))
        c = l
        if len(finished) == len(jobs):
            time.sleep(2)
            if len(finished) == len(jobs):
                break

        while (len(running) < par['simultaneous_kinbot'] and
               len(running) + len(finished) < len(jobs)):
            # start a new job
            job = jobs[len(running) + len(finished)]
            kb = 1
            logger.info('Job: {}'.format(job))
            if job in par['skip_chemids'] and 'none' not in par['skip_chemids']:
                kb = 0
            if job not in par['keep_chemids'] and 'none' not in par['keep_chemids']:
                kb = 0
            logger.info(f'kb: {kb} for {job}')
            if kb == 1:
                pid = 0
                if not no_kinbot:
                    pid = submit_job(job, par)  # kinbot is submitted here
                else:
                    get_wells(job)
                pids[job] = pid
                t = datetime.datetime.now()
                logger.info('Started job {} at {}'.format(job, t))
                running.append(job)
            elif kb == 0:
                logger.info('Skipping Kinbot for {}'.format(job))
                finished.append(job)
            else:
                logger.info('kb value not 0 or 1')

        # check if a thread is done
        for job in running:
            if not check_status(job, pids[job]):
                t = datetime.datetime.now()
                logger.info('Finished job {} at {}'.format(job, t))
                finished.append(job)
                if not no_kinbot:
                    # write a temporary pes file
                    # remove old xval and im_extent files
                    try:
                        os.remove('{}_xval.txt'.format(par['title']))
                    except OSError:
                        pass
                    try:
                        os.remove('{}_im_extent.txt'.format(par['title']))
                    except OSError:
                        pass
        # remove the finished threads
        for job in finished:
            if job in running:
                running.remove(job)
        if not no_kinbot:
            # write a summary of what is running and finished
            summary_lines = []
            summary_lines.append('Total\t\t{}'.format(len(jobs)))
            summary_lines.append('Running\t\t{}'.format(len(running)))
            summary_lines.append('Finished\t{}'.format(len(finished)))
            summary_lines.append('')
            summary_lines.append('Running:')
            for job in running:
                summary_lines.append('\t{}'.format(job))
            summary_lines.append('')
            summary_lines.append('Finished:')
            for job in finished:
                summary_lines.append('\t{}'.format(job))
            waiting = jobs[len(running) + len(finished):]
            summary_lines.append('')
            summary_lines.append('Waiting:')
            for job in waiting:
                summary_lines.append('\t{}'.format(job))
            with open('pes_summary.txt', 'w') as f:
                f.write('\n'.join(summary_lines) + '\n')
            time.sleep(1)

    # delete skipped jobs from the jobs before sending to postprocess
    for skip in par['skip_chemids']:
        try:
            jobs.pop(jobs.index(skip))
        except ValueError:
            pass

    # only keep the jobs we wanted
    if 'none' not in par['keep_chemids']:
        jobs = par['keep_chemids']

    postprocess(par, jobs, task, names, well0.mass)
    # make molpro inputs for all keys above
    # place submission script in the directory for offline submission
    # read in the molpro energies for the keys in the above three dicts
    # for key in newdict.keys():
    #      print(key)
    # if all energies are there
    # do something like postprocess, but with new energies
    # postprocess_L3(saddle_zpe, well_zpe, prod_zpe, saddle_energy, well_energy, prod_energy, conn)

    # Notify user the search is done
    logger.info('PES search done!')
    try:
        print('PES search done!')
    except OSError:
        pass


def get_wells(job):
    """
    Read the summary file and add the wells to the chemid list
    """
    try:
        summary = open(job + '/summary_' + job + '.out', 'r').readlines()
    except:
        return 0
    with open('chemids', 'r') as f:
        jobs = f.read().split('\n')
    jobs = [ji for ji in jobs if ji != '']

    new_wells = []
    for line in summary:
        if (line.startswith('SUCCESS') or line.startswith("HOMOLYTIC_SCISSION")):
            #Unpack the succesfull lines
            if 'vdW' not in line :
                success, ts_energy, reaction_name, *prod = line.split()
            elif 'vdW' in line : #Unpack differently when a vdW well is in line
                success, ts_energy, reaction_name, *prod, vdW_energy, vdW_direction = line.split()
            
            if (len(prod) == 1 and
                    prod[0] not in jobs and
                    prod[0] not in new_wells):
                new_wells.append(prod[0])
    if len(new_wells) > 0:
        with open('chemids', 'a') as f:
            f.write('\n'.join(new_wells) + '\n')
    

def postprocess(par, jobs, task, names, mass):
    """
    postprocess a pes search
    par: parameters of the search
    jobs: all of the kinbot jobs that were run
    temp: this is a temporary output file writing
    """

    l3done = 1  # flag for L3 calculations to be complete

    # base of the energy is the first well, these are L2 energies
    base_energy, base_zpe = get_energy(jobs, jobs[0], 0, par['high_level'],
                                       conf=par['conformer_search'])
    # L3 energies
    status, base_l3energy = get_l3energy(jobs[0], par)
    if not status:
        l3done = 0
    # list of lists with four elements
    # 0. reactant chemid
    # 1. reaction name
    # 2. products chemid list
    # 3. reaction barrier height
    reactions = []

    # list of the parents for each calculation
    # the key is the name of the calculation
    # the value is the parent directory,
    # i.e. the well kinbot started from to find
    # this calculation
    parent = {}

    wells = []
    failedwells = []
    bimol_products = []


    # list of reactions for which mp2 energies should be used at L1
    mp2_list = ['R_Addition_MultipleBond', 'reac_birad_recombination_R', 
                'reac_r12_cycloaddition', 'reac_r14_birad_scission']
    
    #list of booleans, length is number of wells
    #True if the well is a vdW well
    do_vdW = []

    # read all the jobs
    for ji in jobs:
        try:
            summary = open(f"{ji}/summary_{ji}.out", "r", newline='\n').readlines()
        except:
            failedwells.append(ji)
            continue
        # read the summary file from after corporate message
        for line in summary[5:]:
            if line.startswith("SUCCESS"):
                #Unpack the succesfull lines
                if 'vdW' not in line :
                    ts_energy, reaction_name, *products = line.split()[1:]
                elif 'vdW' in line : #Unpack differently when a vdW well is in line
                    ts_energy, reaction_name, *products, vdW_energy, vdW_direction = line.split()[1:]
                    vdW_well = f"{reaction_name}{vdW_direction.split('vdW')[1]}"
                if reaction_name in par['skip_reactions']:
                    continue

                reactant = ji
                #products this is the chemid of the product
                if 'none' not in par['keep_chemids']:
                    if len(products) == 1 and products[0] not in par['keep_chemids']:
                        continue
                # calculate the barrier based on the new energy base
                barrier = 0. - base_energy - base_zpe

                # overwrite energies with mp2 energy if needed
                mp2_list = ['R_Addition_MultipleBond', 'reac_birad_recombination_R', 
                        'reac_r12_cycloaddition', 'reac_r14_birad_scission']
                if any([mm in reaction_name for mm in mp2_list]) \
                       and not par['high_level'] \
                       and par['qc'] != 'nn_pes':
                    mp2_energies = get_energy(jobs, jobs[0], 0, par['high_level'], 
                                              mp2=1, conf=par['conformer_search'])
                    base_energy_mp2, base_zpe_mp2 = mp2_energies
                    barrier = 0. - base_energy_mp2 - base_zpe_mp2

                # overwrite energies with bls energy if needed
                if 'barrierless_saddle' in reaction_name and not par['high_level']:
                    bls_energies = get_energy(jobs, jobs[0], 0, par['high_level'],
                                              bls=1, conf=par['conformer_search'])
                    base_energy_bls, base_zpe_bls = bls_energies
                    barrier = 0. - base_energy_bls - base_zpe_bls

                #Different treatment between homolytic scission and normal reaction
                if "hom_sci" in line:
                    barrier = float(ts_energy)
                else:
                    #Save ts energy if there is a ts (eg. not barrierless reaction)
                    ts_energy, ts_zpe = get_energy(jobs, reaction_name, 1, par['high_level'], 
                                               conf=par['conformer_search'])
                    barrier += ts_energy + ts_zpe
                barrier *= constants.AUtoKCAL
                    
                if reactant not in wells:
                    wells.append(reactant)
                    parent[reactant] = reactant
                    do_vdW.append(False)
                
                if len(products) == 1:
                    product = products[0]
                    if product not in wells:
                        if product not in parent:
                            parent[product] = reactant
                        wells.append(product)
                        do_vdW.append(False)
                else:
                    prod_name = '_'.join(sorted(products))
                    if 'vdW' in line:
                        if vdW_well not in wells:
                            wells.append(vdW_well)
                            do_vdW.append(True)
                            parent[vdW_well] = reactant
                            if prod_name not in parent:
                                parent[prod_name] = vdW_well
                    if prod_name not in bimol_products:
                        if prod_name not in parent:
                            parent[prod_name] = reactant
                        bimol_products.append('_'.join(sorted(products)))
                new = 1
                temp = None

                for i, rxn in enumerate(reactions):
                    rxn_prod_name = '_'.join(sorted(rxn[2]))
                    if (reactant == rxn[0] and
                            '_'.join(sorted(products)) == rxn_prod_name):
                        new = 0
                        temp = i
                    if reactant == ''.join(rxn[2]) and ''.join(products) == rxn[0]:
                        new = 0
                        temp = i

                if new :
                    if "vdW" not in line:
                        reactions.append([reactant, reaction_name, products, barrier])
                    else:
                        reactions.append([reactant, reaction_name, products, barrier, vdW_energy, vdW_direction])
                elif not new:
                    if "hom_sci" not in reaction_name:
                        # check if the previous reaction has a lower energy or not
                        if reactions[temp][3] > barrier:
                            reactions.pop(temp)
                            if "vdW" not in line:
                                reactions.append([reactant, reaction_name, products, barrier])
                                if len(reactions[temp]) == 6: #True when reactions[temp] has a vdW well
                                    parent[prod_name] = reactant
                            else:
                                #replace the prod parent by the vdW well
                                parent[prod_name] = vdW_well 
                                reactions.append([reactant, reaction_name, products, barrier, vdW_energy, vdW_direction])
                        elif reactions[temp][3] < barrier and "vdW" in line:
                            wells.pop(-1)
                            do_vdW.pop(-1)
                            parent.pop(vdW_well)
                        elif reactions[temp][3] == barrier:
                            if "vdW" in line:
                                reactions.pop(temp)
                                parent[prod_name] = vdW_well 
                                reactions.append([reactant, reaction_name, products, barrier, vdW_energy, vdW_direction])
                            else:
                                #If reaction exitst with same barrier/products, but with a vdW well, keep it and discard new one.
                                continue
                        elif "hom_sci" in reactions[temp][1]:
                            reactions.pop(temp)


        # copy the xyz files
        copy_from_kinbot(ji, 'xyz')
        # copy the L3 calculations here, whatever was in those directories, inp, out, pbs, etc.
        try:
            copy_from_kinbot(ji, par['single_point_qc'])
        except:
            logger.warning(f'L3 calculations were not copied from {ji}')
    # create a connectivity matrix for all wells and products
    conn, bars = get_connectivity(wells, bimol_products, reactions)
    # create a batch submission for all L3 jobs
    if par['queuing'] == 'pbs':
        cmd = 'qsub'
        ext = 'pbs'
    elif par['queuing'] == 'slurm':
        cmd = 'sbatch'
        ext = 'slurm'
    elif par['queuing'] == 'local':
        cmd = ''
        ext = ''
        pass
    else:
        raise ValueError(f'Unexpected value for queueing: {par["queuing"]}')

    batch_submit = ''

    well_energies = {}
    well_l3energies = {}
    for index, well in enumerate(wells):
        energy, zpe = get_energy(wells, well, do_vdW[index], par['high_level'], 
                            conf=par['conformer_search'])  # from the db
        well_energies[well] = ((energy + zpe) - (base_energy + base_zpe)) * constants.AUtoKCAL
        status, l3energy = get_l3energy(well, par)
        if not status:
            l3done = 0  # not all L3 calculations are done
            batch_submit += f'{cmd} {well}.{ext}\n'
        elif not par['L3_calc']:
            pass
        else:
            well_l3energies[well] = ((l3energy + zpe) - (base_l3energy + base_zpe)) * constants.AUtoKCAL
    prod_energies = {}
    prod_l3energies = {}
    for prods in bimol_products:
        energy = 0. - (base_energy + base_zpe)
        l3energy = 0. - (base_l3energy + base_zpe)
        for pr in prods.split('_'):
            pr_energy, pr_zpe = get_energy(jobs, pr, 0, par['high_level'], 
                                           conf=par['conformer_search'])
            energy += pr_energy + pr_zpe
            status, l3e = get_l3energy(pr, par)
            if not status:
                l3done = 0  # not all L3 calculations are done
                batch_submit += f'{cmd} {pr}.{ext}\n'
            elif not par['L3_calc']:
                pass
            else:
                l3energy += l3e + pr_zpe
        prod_energies[prods] = energy * constants.AUtoKCAL
        prod_l3energies[prods] = l3energy * constants.AUtoKCAL

    ts_l3energies = {}
    for reac in reactions:
        if "hom_sci" not in reac[1]:  # meaning homolytic scission
            if len(reac) > 4:#vdW 
                vdW_well = f"{reac[1]}{reac[5].split('vdW')[1]}"
                reac[4] = well_energies[vdW_well]
            if 'barrierless_saddle' in reac[1]:
                zpe = get_zpe(reac[0], reac[1], 1, par['high_level'])
                status, l3energy = get_l3energy(reac[1], par, bls=1)

                status_prod1, l3energy_prod1 = get_l3energy(reac[2][0], par)
                status_prod2, l3energy_prod2 = get_l3energy(reac[2][1], par)

                status_prod, l3energy_prod = get_l3energy(reac[1] + '_prod', par, bls=1)

                if not status * status_prod:
                    l3done = 0
                    batch_submit += f'{cmd} {reac[1]}.{ext}\n'
                elif not par['L3_calc']:
                    pass
                else:
                    delta1 = l3energy_prod - (l3energy + zpe)  # ZPEs cancel out for fragments
                    delta2 = l3energy_prod1 + l3energy_prod2 - (base_l3energy + base_zpe) 
                    ts_l3energies[reac[1]] = (delta2 - delta1) * constants.AUtoKCAL
            else:
                zpe = get_zpe(reac[0], reac[1], 1, par['high_level'])
                status, l3energy = get_l3energy(reac[1], par)
                if not status:
                    l3done = 0
                    batch_submit += f'{cmd} {reac[1]}.{ext}\n'
                elif not par['L3_calc']:
                    pass
                else:
                    ts_l3energies[reac[1]] = ((l3energy + zpe) - (base_l3energy + base_zpe)) * constants.AUtoKCAL

    logger.info('l3done status {}'.format(l3done))
    # clean duplicates
    batch_submit = list(set(batch_submit.split('\n')))
    batch_submit.reverse()
    batch_submit = '\n'.join(batch_submit)
    batch = f'{par["single_point_qc"]}/batch_L3_{par["queuing"]}.sub'
    if par['queuing'] != 'local':
        with open(batch, 'w') as f:
            f.write(batch_submit)
        os.chmod(batch, stat.S_IRWXU)  # read, write, execute by owner

    if l3done == 1 and par['L3_calc']:
        logger.info('Energies are updated to L3 in ME and PESViewer.')
        well_energies = well_l3energies
        prod_energies = prod_l3energies
        for reac in reactions:  # swap out the barrier
            if 'hom_sci' not in reac[1]:
                reac[3] = ts_l3energies[reac[1]]
            if 'barrierless_saddle' in reac[1]:
                reac[3] = ts_l3energies[reac[1]]
            if len(reac) > 4: #vdW 
                vdW_well = f"{reac[1]}{reac[5].split('vdW')[1]}"
                reac[4] = well_l3energies[vdW_well]

        logger.info('L3 energies in kcal/mol, incl. ZPE')
        for well in wells:
            logger.info('{}   {:.2f}'.format(well, well_l3energies[well]))
        for prod in bimol_products:
            logger.info('{}   {:.2f}'.format(prod, prod_l3energies[prod]))
        for ts in ts_l3energies:
            logger.info('{}   {:.2f}'.format(reaction_name, ts_l3energies[ts]))
    else:
        logger.info(f'Energies used are at the L2 ({par["high_level_method"]}/'
                     f'{par["high_level_basis"]}) level of theory.')

        
    # if L3 was done and requested, everything below is done with that
    # filter according to tasks
    filtered_stpts = filter_stat_points(par, wells, bimol_products, reactions, conn,
                                        bars, well_energies, task, names)
    wells, products, reactions, highlight = filtered_stpts

    create_interactive_graph(wells,
                             bimol_products,
                             reactions,
                             par['title'],
                             well_energies,
                             prod_energies,
                             )

    barrierless = []
    vdW = []
    rxns = []
    for rxn in reactions:
        if 'hom_sci' in rxn[1]:
            #in rxn: [reactant, reaction_name, products, barrier]
            barrierless.append([rxn[0], rxn[1], rxn[2], rxn[3]])
        elif len(rxn) > 4: #vdW for now, find better condition latter?
            #in rxn: [reactant, reaction_name, products, barrier, vdW_energy, vdW_direction]
            vdW.append([rxn[0], rxn[1], rxn[2], rxn[3], rxn[4], rxn[5]])
        else:
            rxns.append([rxn[0], rxn[1], rxn[2], rxn[3]])

    create_rotdPy_inputs(par,
                         barrierless,
                         vdW)


    # write full pesviewer input
    create_pesviewer_input(par,
                           wells,
                           bimol_products,
                           rxns,
                           barrierless,
                           vdW,
                           well_energies,
                           prod_energies,
                           highlight)
    
    if par['me']:
        create_mess_input(par,
                          deepcopy(wells),
                          bimol_products,
                          deepcopy(rxns),
                          deepcopy(barrierless),
                          deepcopy(vdW),
                          well_energies,
                          prod_energies,
                          parent,
                          mass,
                          l3done)

    #if par['single_point_qc'].lower() == 'molpro':
    #    if l3done:
    #        check_l3_l2(par['single_point_key'], parent, reactions)
    #    t1_analysis(par['single_point_key'])


def filter_stat_points(par, wells, products, reactions, conn, bars, well_energies, task,
           names):
    """Filter the wells, products and reactions according to their task and name."""
    # list of reactions to highlight
    highlight = []
    # 1. all: This is the default showing all pathways
    # 2. lowestpath: show the lowest path between the species
    # corresponding to the names
    # 3. allpaths: show all paths between the species
    # corresponding to the names
    # 4. wells: show all reactions of one wells
    # corresponding to the names
    # 5. temperature
    # 6. threshold_reapply: apply the barrier threshold
    # cutoff at the highest level that was done

    # filter the reactions according to the task
    if task == 'all':
        filtered_reactions = reactions
    elif task == 'lowestpath':
        all_rxns = get_all_pathways(wells, products, reactions, names, conn)
        # this is the maximum energy along the minimun energy pathway
        min_energy = None
        min_rxn = None
        for rxn_list in all_rxns:
            barriers = [ri[3] for ri in rxn_list]
            if min_energy is None:
                min_energy = max(barriers)
                min_rxn = rxn_list
            else:
                if max(barriers) < min_energy:
                    min_energy = max(barriers)
                    min_rxn = rxn_list
        filtered_reactions = min_rxn
    elif task == 'allpaths':
        all_rxns = get_all_pathways(wells, products, reactions, names, conn)
        filtered_reactions = []
        for rxnlist in all_rxns:
            for rxn in rxnlist:
                new = 1
                for r in filtered_reactions:
                    if r[1] == rxn[1]:
                        new = 0
                if new:
                    filtered_reactions.append(rxn)
        # this is the maximum energy along the minimun energy pathway
        min_energy = None
        min_rxn = None
        for rxn_list in all_rxns:
            barriers = [ri[3] for ri in rxn_list]
            if min_energy is None:
                min_energy = max(barriers)
                min_rxn = rxn_list
            else:
                if max(barriers) < min_energy:
                    min_energy = max(barriers)
                    min_rxn = rxn_list
        for rxn in min_rxn:
            highlight.append(rxn[1])
    elif task == 'well':
        if len(names) == 1:
            filtered_reactions = []
            for rxn in reactions:
                prod_name = '_'.join(sorted(rxn[2]))
                if names[0] == rxn[0] or names[0] == prod_name:
                    filtered_reactions.append(rxn)
        else:
            logger.error('Only one name should be given for a well filter')
            logger.error('Received: ' + ' '.join(names))
            sys.exit(-1)
    elif task == 'temperature':
        if len(names) == 1:
            try:
                # read the temperature
                temperature = float(names[0])
            except ValueError:
                logger.error('A float is needed for a temperature filter')
                logger.error('Received: ' + ' '.join(names))
                sys.exit(-1)
            filtered_reactions = []
            # iterate the wells
            wells, filtered_reactions = filter_boltzmann(wells[0],
                                                         [wells[0]],
                                                         reactions,
                                                         filtered_reactions,
                                                         well_energies,
                                                         temperature)
        else:
            logger.error('Only one argument should be given for a temperature filter')
            logger.error('Received: ' + ' '.join(names))
            sys.exit(-1)
    elif task == 'l2threshold':
        filtered_reactions = []
        for rxn in reactions:
            if rxn[3] < par['barrier_threshold']:
                filtered_reactions.append(rxn)
    else:
        logger.error('Could not recognize task ' + task)
        sys.exit(-1)

    # filter the wells
    filtered_wells = []
    for well in wells:
        for rxn in filtered_reactions:
            prod_name = '_'.join(sorted(rxn[2]))
            if well == rxn[0] or well == prod_name:
                if well not in filtered_wells:
                    filtered_wells.append(well)
            if "IRC" in well and (len(rxn) == 6):
                vdW_name = f"{rxn[1]}{rxn[5].split('vdW')[1]}"
                if well == vdW_name:
                    if well not in filtered_wells:
                        filtered_wells.append(well)

    # filter the products
    filtered_products = []
    for prod in products:
        for rxn in filtered_reactions:
            prod_name = '_'.join(sorted(rxn[2]))
            if prod == prod_name:
                if prod not in filtered_products:
                    filtered_products.append(prod)

    return filtered_wells, filtered_products, filtered_reactions, highlight


def filter_boltzmann(well, wells, reactions, filtered_reactions,
                     well_energies, temperature):
    """
    Filter the reactions based on branching fractions at a given temperature
    """
    well_energy = well_energies[well]
    # all the reactions that include this well
    one_well_rxns = []
    # iterate all reactions
    for rxn in reactions:
        # check if this reactions belongs to the current well
        if rxn[0] == well or rxn[2][0] == well:
            one_well_rxns.append(rxn)
    # list containing the branching fractions for this well
    branching = []
    for rxn in one_well_rxns:
        # calculate the boltzmann factor
        value = np.exp(-(rxn[3] - well_energy) * 1000 / 1.9872036 / temperature)
        branching.append(value)
    # calculate the actual branching fractions
    br_sum = sum(branching)
    branching = np.array(branching) / br_sum
    for i, rxn in enumerate(one_well_rxns):
        # only add a reaction if the branching fraction is below 1%
        if branching[i] > 0.01:
            new = 1
            for r in filtered_reactions:
                if r[1] == rxn[1]:
                    new = 0
            if new:
                filtered_reactions.append(rxn)
                if not rxn[0] in wells:
                    wells.append(rxn[0])
                    wells, filtered_reactions = filter_boltzmann(rxn[0],
                                                                 wells,
                                                                 reactions,
                                                                 filtered_reactions,
                                                                 well_energies,
                                                                 temperature)
                if len(rxn[2]) == 1:
                    if not rxn[2][0] in wells:
                        wells.append(rxn[2][0])
                        wells, filtered_reactions = filter_boltzmann(rxn[2][0],
                                                                     wells,
                                                                     reactions,
                                                                     filtered_reactions,
                                                                     well_energies,
                                                                     temperature)
    return wells, filtered_reactions


def get_connectivity(wells, products, reactions):
    """
    Create two matrices:
    conn: connectivity (1 or 0) between each pair of stationary points
    bars: barrier height between each pair of stationary points
          0 if not connected
    """
    conn = np.zeros((len(wells) + len(products), len(wells) + len(products)), dtype=int)
    bars = np.zeros((len(wells) + len(products), len(wells) + len(products)))
    for rxn in reactions:
        reac_name = rxn[0]
        prod_name = '_'.join(sorted(rxn[2]))
        i = get_index(wells, products, reac_name)
        j = get_index(wells, products, prod_name)
        conn[i][j] = 1
        conn[j][i] = 1
        barrier = rxn[3]
        bars[i][j] = barrier
        bars[j][i] = barrier
    return conn, bars


def get_all_pathways(wells, products, reactions, names, conn):
    """
    Get all the pathways in which all intermediate species
    are wells and not bimolecular products
    """
    if len(names) == 2:
        # the maximum length between two stationary points
        # is the number of wells+2
        max_length = 5
        n_mol = len(wells) + len(products)
        start = get_index(wells, products, names[0])
        end = get_index(wells, products, names[1])
        # make a graph out of the connectivity
        # nodes of the graph
        nodes = [i for i in range(n_mol)]
        G = nx.Graph()
        G.add_nodes_from(nodes)
        # add the edges of the graph
        for i, ci in enumerate(conn):
            for j, cij in enumerate(ci):
                if cij > 0:
                    G.add_edge(i, j)
        # list of reaction lists for each pathway
        paths = nx.all_simple_paths(G, start, end, cutoff=max_length)
        rxns = []
        for path in paths:
            if is_pathway(wells, products, path, names):
                rxns.append(get_pathway(wells, products, reactions, path, names))
        return rxns
    else:
        logger.error('Cannot find a lowest path if the number of species is not 2')
        logger.error('Found species: ' + ' '.join(names))


def get_index(wells, products, name):
    try:
        i = wells.index(name)
    except ValueError:
        try:
            i = products.index(name) + len(wells)
        except ValueError:
            logger.error('Could not find reactant ' + name)
            sys.exit(-1)
    return i


def get_name(wells, products, i):
    if i < len(wells):
        name = wells[i]
    else:
        name = products[i - len(wells)]
    return name


def get_pathway(wells, products, reactions, ins, names):
    """
    Return the list of reactions between the species in
    the names, according to the instance ins
    """
    # list of reactions
    rxns = []
    for index, i in enumerate(ins[:-1]):
        j = ins[index + 1]
        rxns.append(get_reaction(wells, products, reactions, i, j))
    return rxns


def get_reaction(wells, products, reactions, i, j):
    """
    method to get a reaction in the reactions list
    according to the indices i and j which correspond
    to the index in wells or products
    """
    name_1 = get_name(wells, products, i)
    name_2 = get_name(wells, products, j)
    for rxn in reactions:
        reac_name = rxn[0]
        prod_name = '_'.join(sorted(rxn[2]))
        if ((name_1 == reac_name and name_2 == prod_name) or
                (name_2 == reac_name and name_1 == prod_name)):
            return rxn
    return None


def is_pathway(wells, products, ins, names):
    """
    Method to check if the instance ins
    corresponds to a pathway between the species
    in the names list
    """
    # check of all intermediate species are wells
    if all([insi < len(wells) for insi in ins[1:-1]]):
        name_1 = get_name(wells, products, ins[0])
        name_2 = get_name(wells, products, ins[-1])
        # check if the names correspond
        if ((name_1 == names[0] and name_2 == names[1]) or
                (name_2 == names[0] and name_1 == names[1])):
            return 1
    return 0


def copy_from_kinbot(well, dirname):
    files = os.listdir(f'{well}/{dirname}')
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    for f in files:
        if f.endswith('.out'):
            if not os.path.exists(f'{dirname}/{f}'):
                shutil.copy(f'{well}/{dirname}/{f}', f'{dirname}/{f}')
        else:
            shutil.copy(f'{well}/{dirname}/{f}', f'{dirname}/{f}')


def get_rxn(prods, rxns):
    for rxn in rxns:
        if prods == '_'.join(sorted(rxn[2])):
            return rxn


def create_short_names(wells, products, reactions, barrierless, vdW):
    """
    Create a short name for all the wells, all the
    bimolecular products and all the transition states
    """
    # short names of the wells
    # key: chemid
    # value: short name
    well_short = {}
    # short names of the products
    # key: sorted chemids with underscore
    # value: short name
    pr_short = {}
    # short names of the fragments
    # key: chemid
    # value: short name
    fr_short = {}
    # short name of the ts's
    # key: reaction name (chemid, type and instance)
    # value: short name
    ts_short = {}
    # short name of the barrierless reactions
    # key: reaction number
    # value: reaction number
    nobar_short = {}
    for well in wells:
        if well not in well_short:
            short_name = 'w_' + str(len(well_short) + 1)
            well_short[well] = short_name

    for prod in products:
        if prod not in pr_short:
            short_name = 'pr_' + str(len(pr_short) + 1)
            pr_short[prod] = short_name
        for frag in prod.split('_'):
            if frag not in fr_short:
                short_name = 'fr_' + str(len(fr_short) + 1)
                fr_short[frag] = short_name

    for rxn in reactions:
        if rxn[1] not in ts_short:
            short_name = 'rxn_' + str(len(ts_short) + 1)
            ts_short[rxn[1]] = short_name

    n = 1
    for rxn in barrierless:
        nobar_name = 'nobar_' + str(n)
        key = f'{rxn[0]}_{rxn[2][0]}_{rxn[2][1]}'
        nobar_short[key] = nobar_name
        n = n + 1
    
    for vdw in vdW:
        if f"{vdw[1]}{vdw[5].split('vdW')[1]}" not in well_short:
            short_name = 'w_' + str(len(well_short) + 1)
            well_short[f"{vdw[1]}{vdw[5].split('vdW')[1]}"] = short_name
        nobar_name = 'nobar_' + str(n)
        key = f"{vdw[1]}{vdw[5].split('vdW')[1]}_{vdw[2][0]}_{vdw[2][1]}"
        nobar_short[key] = nobar_name
        n = n + 1
        

    return well_short, pr_short, fr_short, ts_short, nobar_short


def create_mess_input(par, wells, products, reactions, barrierless, vdW,
                      well_energies, prod_energies, parent, mass, l3done):
    """When calculating a full pes, the files from the separate wells
    are read and concatenated into one file
    Two things per file need to be updated
    1. the names of all the wells, bimolecular products and ts's
    2. all the zpe corrected energies
    If it is a multi_conf_tst calculation, then the ZeroEnergies and
    ImaginaryFrequency values need to be also changed for each species in the Union.
    The lowest energy species has L3 energies (if calculated), while the 
    others are corrected with their DeltaL2 energies.
    """

    min_vdW = find_min_vdW(vdW, well_energies)

    for well in wells:
        if well in min_vdW and\
        min_vdW[well] != well:
            wells.pop(wells.index(well))

    logger.info(f"uq value: {par['uq']}")
    for vdw in vdW:
        vdw_name = f"{vdw[1]}{vdw[5].split('vdW')[1]}"
        vdW_well_reac = [vdw[0], vdw[1], [min_vdW[vdw_name]], vdw[3]]
        if vdW_well_reac not in reactions:
            reactions.append(vdW_well_reac)
    short_names = create_short_names(wells, products, reactions, barrierless, vdW)
    well_short, pr_short, fr_short, ts_short, nobar_short = short_names

    # list of the strings to write to mess input file
    s = []

    # Create the header block for MESS
    frame = '######################\n' 
    divider = '! ****************************************\n'
    
    dummy = StationaryPoint('dummy',
                            par['charge'],
                            par['mult'],
                            smiles=par['smiles'],
                            structure=par['structure'])

    mess = MESS(par, dummy)
    uq = UQ(par)

    if l3done:
        lot = 'L3'
    else:
        lot = f'{par["high_level_method"]}/{par["high_level_basis"]}'

    header_file = f'{kb_path}/tpl/mess_header.tpl'
    with open(header_file) as f:
        tpl = f.read()

    for uq_iter in range(par['uq_n']):
        mess_iter = "{0:04d}".format(uq_iter)

        e_well = par['epsilon'] * uq.calc_factor('epsilon', uq_iter)
        s_well = par['sigma'] * uq.calc_factor('sigma', uq_iter)
        enrelfact = par['EnergyRelaxationFactor'] * uq.calc_factor('enrelfact', uq_iter)
        enrelpow = par['EnergyRelaxationPower'] * uq.calc_factor('enrelpow', uq_iter)

        header = tpl.format(TemperatureList=' '.join([str(ti) for ti in par['TemperatureList']]),
                            PressureList=' '.join([str(pi) for pi in par['PressureList']]),
                            EnergyStepOverTemperature=par['EnergyStepOverTemperature'],
                            ExcessEnergyOverTemperature=par['ExcessEnergyOverTemperature'],
                            ModelEnergyLimit=par['ModelEnergyLimit'],
                            CalculationMethod=par['CalculationMethod'],
                            ChemicalEigenvalueMax=par['ChemicalEigenvalueMax'],
                            Reactant=well_short[wells[0]],
                            EnergyRelaxationFactor=round(enrelfact, 2),
                            EnergyRelaxationPower=round(enrelpow, 2),
                            EnergyRelaxationExponentCutoff=par['EnergyRelaxationExponentCutoff'],
                            e_coll=round(constants.epsilon[par['collider']], 2),
                            s_coll=constants.sigma[par['collider']],
                            m_coll=constants.mass[par['collider']],
                            e_well=round(e_well, 2),
                            s_well=round(s_well, 2),
                            m_well=mass,
                            LevelOfTheory=lot
                            )

        s = []  # mess string after header

        well_energies_current = deepcopy(well_energies)
        prod_energies_current = deepcopy(prod_energies)
            
        # write the wells
        s.append(frame + '# WELLS\n' + frame)
        for well in wells:
            name = well_short[well] + ' ! ' + well
            energy = well_energies[well] + uq.calc_factor('energy', uq_iter)
            well_energies_current[well] = energy
            # if 'IRC' not in well:
            with open(parent[well] + '/' + well + '_' + mess_iter + '.mess', 'r') as f:
                s.append(f.read().format(name=name, zeroenergy=round(energy, 2)))
            s.append(divider)
            # else:
                # with open(parent[well] + '/' + parent[well] + '_' + mess_iter + '.mess', 'r') as f:
                    # s.append(f.read().format(name=name, zeroenergy=round(energy, 2)))
                # s.append(divider)

        # write the products and barrierless reactions
        s.append(frame + '# BIMOLECULAR PRODUCTS\n' + frame)
        for prod in products:
            linked_bless_wells = []
            name = pr_short[prod] + ' ! ' + prod
            energy = prod_energies[prod] + uq.calc_factor('energy', uq_iter)
            prod_energies_current[prod] = energy
            fr_names = {}
            for fr in prod.split('_'):
                key = f'fr_name_{fr}'
                value = fr_short[fr] + ' ! ' + fr
                fr_names[key] = value
            # check if is obtained by barrieless reaction
            bless = 0
            for bl in barrierless:
                bl_prod = f'{bl[2][0]}_{bl[2][1]}'
                if prod == bl_prod:
                    if bl[0] in linked_bless_wells:
                        continue
                    linked_bless_wells.append(bl[0])
                    with open(bl[0] + '/' + prod + '_' + mess_iter + '.mess') as f:
                        if bless == 0:
                            s.append(f.read().format(name=name,
                                                     blessname=nobar_short[f'{bl[0]}_{bl_prod}'],
                                                     wellname=well_short[bl[0]],
                                                     prodname=pr_short[prod],
                                                     ground_energy=round(energy, 2),
                                                     **fr_names))
                            bless = 1
                        else:
                            stemp = (f.read().format(name=name,
                                                     blessname=nobar_short[f'{bl[0]}_{bl_prod}'],
                                                     wellname=well_short[bl[0]],
                                                     prodname=pr_short[prod],
                                                     ground_energy=round(energy, 2),
                                                     **fr_names))
                            s.append(stemp[stemp.find('Barrier '):])
            # check if is obtained by vdW dissociation
            for vdw in vdW:
                bl_prod = f'{vdw[2][0]}_{vdw[2][1]}'
                if prod == bl_prod:
                    #Skip this vdW well is bless already added with similar well
                    if min_vdW[f"{vdw[1]}{vdw[5].split('vdW')[1]}"] in linked_bless_wells:
                        continue
                    linked_bless_wells.append(min_vdW[f"{vdw[1]}{vdw[5].split('vdW')[1]}"])
                    with open(vdw[0] + '/' + prod + '_' + mess_iter + '.mess') as f:
                        if bless == 0:
                            s.append(f.read().format(name=name,
                                                     blessname=nobar_short[f"{vdw[1]}{vdw[5].split('vdW')[1]}_{bl_prod}"],
                                                     wellname=well_short[min_vdW[f"{vdw[1]}{vdw[5].split('vdW')[1]}"]],
                                                     prodname=pr_short[prod],
                                                     ground_energy=round(energy, 2),
                                                     **fr_names))
                            bless = 1
                        else:
                            stemp = (f.read().format(name=name,
                                                     blessname=nobar_short[f"{vdw[1]}{vdw[5].split('vdW')[1]}_{bl_prod}"],
                                                     wellname=well_short[min_vdW[f"{vdw[1]}{vdw[5].split('vdW')[1]}"]],
                                                     prodname=pr_short[prod],
                                                     ground_energy=round(energy, 2),
                                                     **fr_names))
                            s.append(stemp[stemp.find('Barrier '):])
            #The product is directly obtained after a reaction with a barrier
            if not bless:
                if 'IRC' not in parent[prod]:
                    try:
                        with open(parent[prod] + '/' + prod + '_' + mess_iter + '.mess') as f:
                            s.append(f.read().format(name=name,
                                                    ground_energy=round(energy, 2),
                                                    **fr_names))
                    except:#When bimolecular template is used both for barrierless and with barrier
                        with open(parent[prod] + '/' + prod + '_' + mess_iter + '.mess') as f:
                            file = f.readlines()
                        f = ''
                        for line in file:
                            if "Barrier" in line:
                                break
                            else:
                                f += line
                        s.append(f.format(name=name,
                                          ground_energy=round(energy, 2),
                                          **fr_names))

                else:
                    with open(parent[parent[prod]] + '/' + prod + '_' + mess_iter + '.mess') as f:
                        s.append(f.read().format(name=name,
                                                ground_energy=round(energy, 2),
                                                **fr_names))
                
            s.append(divider)

        # write the barrier
        s.append(frame + '# BARRIERS\n' + frame)
        for rxn in reactions:
            if rxn[0] == rxn[2][0]:  # Avoid writing identity reactions.
                continue
            with open("reactionList.log", 'a') as f:
                f.write(f'rxn\n')
            name = [ts_short[rxn[1]]]
            name.append(well_short[rxn[0]])
            if len(rxn[2]) == 1:
                name.append(well_short[rxn[2][0]])
            else:
                name.append(pr_short['_'.join(sorted(rxn[2]))])
            name.append('!')
            name.append(rxn[1])
            energy = rxn[3] + uq.calc_factor('barrier', uq_iter)
            welldepth1 = energy - well_energies_current[rxn[0]] 
            if welldepth1 < 0 and par['correct_submerged']== 1:  # submerged, not allowed in MESS
                # tunneling block for submerged is cleaned later
                energy = well_energies_current[rxn[0]]
                logger.warning(f'Submerged barrier corrected for {name}')
            if len(rxn[2]) == 1:
                welldepth2 = energy - well_energies_current[rxn[2][0]] 
                if welldepth2 < 0 and par['correct_submerged'] == 1:  # submerged, not allowed in MESS
                    energy = well_energies_current[rxn[2][0]]
                    logger.warning(f'Submerged barrier corrected for {name}')
            else:
                prodname = '_'.join(sorted(rxn[2]))
                welldepth2 = energy - prod_energies_current[prodname] 
                if welldepth2 < 0 and par['correct_submerged'] == 1:  # submerged, not allowed in MESS
                    energy = prod_energies_current[prodname]
                    logger.warning(f'Submerged barrier corrected for {name}')
            cutoff = min(welldepth1, welldepth2)
            with open(rxn[0] + '/' + rxn[1] + '_' + mess_iter + '.mess') as f:
                s.append(f.read().format(name=' '.join(name), 
                         zeroenergy=round(energy, 2),
                         cutoff=round(cutoff, 2),
                         welldepth1=round(welldepth1, 2),
                         welldepth2=round(welldepth2, 2),
                         ))
            s.append('!****************************************')
        s.append(divider)
        s.append('End ! end kinetics\n')

        if not os.path.exists('me'):
            os.mkdir('me')

        with open(f'me/mess_{mess_iter}.inp', 'w') as f:
            f.write(header)
            f.write('\n'.join(s))

        if par['multi_conf_tst']:
            logger.debug('\tUpdating ZPE and tunneling parameters for multi_conf_tst...')
            with open(f'me/mess_{mess_iter}_corr.inp', 'w') as fcorr:
                with open(f'me/mess_{mess_iter}.inp', 'r') as f:
                    lines = f.read().split('\n')
                    for line in lines:
                        words = line.split()
                        if 'ZeroEnergy' in line:
                            if len(words) == 4:
                                space = ' ' * (len(line) - len(line.lstrip(' ')))
                                ecorr = float(words[3])
                                words[1] = str(round(float(words[1]) + ecorr, 2))
                                newline = ' '.join(words)
                                newline = space + newline
                            else:
                                newline = line
                            fcorr.write(newline)
                        elif 'End ! RRHO' in line:
                            ecorr = None  # reset correction at the end of block
                            fcorr.write(line)
                        elif len(words) == 0:
                            continue
                        elif ('CutoffEnergy' in line or 'WellDepth' in line) and ecorr is not None:
                            space = ' ' * (len(line) - len(line.lstrip(' ')))
                            words[1] = str(round(float(words[1]) + ecorr, 2))
                            newline = ' '.join(words)
                            newline = space + newline
                            fcorr.write(newline)
                        elif 'ImaginaryFrequency' in line:
                            if len(words) == 4:
                                space = ' ' * (len(line) - len(line.lstrip(' ')))
                                words[1] = words[3][1:]  # cutting off - sign
                                newline = ' '.join(words)
                                newline = space + newline
                            else:
                                newline = line
                            fcorr.write(newline)
                        else:
                            fcorr.write(line)
                        fcorr.write('\n')

            shutil.copyfile(f'me/mess_{mess_iter}_corr.inp', f'me/mess_{mess_iter}.inp')
            os.remove(f'me/mess_{mess_iter}_corr.inp')

        # cleaning file from submerged barrier tunneling
        #Tunneling   Eckart
        #  ImaginaryFrequency[1/cm]  2277.51
        #  CutoffEnergy[kcal/mol]    20.56
        #  WellDepth[kcal/mol]       34.0
        #  WellDepth[kcal/mol]       20.56
        #End
        with open(f'me/mess_{mess_iter}_temp.inp', 'w') as ftemp:
            with open(f'me/mess_{mess_iter}.inp', 'r') as f:
                lines = f.read().split('\n')
                submerged = -1
                for ll, line in enumerate(lines):
                    words = line.split()
                    if 'Tunneling' in line:
                        submerged = -1
                        words = lines[ll + 2].split()
                        if float(words[1]) <= 0:
                            submerged = 0
                            ftemp.write('! submerged barrier\n')
                        else:
                            ftemp.write(line)
                            ftemp.write('\n')
                    elif submerged <= 5 and submerged >= 0:
                        submerged += 1
                    else: 
                        ftemp.write(line)
                        ftemp.write('\n')

        shutil.copyfile(f'me/mess_{mess_iter}_temp.inp', f'me/mess_{mess_iter}.inp')
        os.remove(f'me/mess_{mess_iter}_temp.inp')

        if par['me']:
            mess.run()

        #uq.format_uqtk_data() 
    return

def create_rotdPy_inputs(par, bless, vdW):
    """
    Function that create an input file for rotdPy.
    barrierless and vdW are lists of reactions.
    barrierless reactions:
    [reactant, reaction_name, products, barrier]
    vdW reactions:
    [reactant, reaction_name, products, barrier, vdW_energy, vdW_direction]
    """

    barrierless = list(bless) #Avoids modifying barrierless outside of the function
    number_of_barrierless = len(barrierless)

    #format the vdW reactions to be added to the barrierless list.
    for reac in vdW:
        reactant, reaction_name, products, barrier, vdW_energy, vdW_direction = reac
        reaction_name = f"{reaction_name}{vdW_direction.split('vdW')[1]}"
        barrier = vdW_energy
        barrierless.append([reactant, reaction_name, products, barrier])

    for index, reac in enumerate(barrierless):
        reactant, reaction_name, products, barrier = reac
        job_name = f"{reaction_name}_{'_'.join(sorted(products))}"
        rotdpy_inputs = []
        if len(par["rotdPy_inputs"]) != 0:
            for this_input in par["rotdPy_inputs"]:
                inp_products = this_input.split("_")[-2:]
                rotdpy_inputs.append(f"{'_'.join(this_input.split('_')[:-2])}_{'_'.join(sorted(inp_products))}")
        if job_name not in rotdpy_inputs:
            continue
        logger.info(f"Creating rotdPy input for reaction {job_name}")
        parent_chemid = reactant
        if len(products) != 2: #Check if the barrierless reaction has 2 fragments
            logger.warning("The creation of rotdPy inputs currently only work for bimolecular products.")
            logger.warning(f"Skiping rotdPy input creation for reac {job_name}.")
            continue
        tot_frag = len(products)

        vrc_tst_start = 0
        do_correction = True
        corrections = {"1d":{"type": "1d"}}
        for scan_type in ["","_frozen"]:
            if "IRC" in job_name:
                plt_file = f"{job_name.split('IRC')[0]}vrc_tst_scan{scan_type}{job_name.split('prod')[1]}_plt.py"
                if not os.path.isfile(f"{parent_chemid}/{plt_file}"):
                    logger.warning(f"Skipping correction for rotdPy job {job_name}: vrc_tst_scan{scan_type} results not found.")
                    do_correction = False
            elif "hom_sci" in job_name: 
                plt_file = f"{'_'.join(job_name.split('_')[:-2])}_vrc_tst_scan{scan_type}_{inp_products[0]}_{inp_products[1]}_plt.py"
                if not os.path.isfile(f"{parent_chemid}/{plt_file}"):
                    plt_file = f"{'_'.join(job_name.split('_')[:-2])}_vrc_tst_scan{scan_type}_{inp_products[1]}_{inp_products[0]}_plt.py"
                    if not os.path.isfile(f"{parent_chemid}/{plt_file}"):
                        logger.warning(f"Skipping correction for rotdPy job {job_name}: vrc_tst_scan{scan_type} results not found.")
                        do_correction = False
            
            pivot_points_info = []
            
            if do_correction:
                with open(f"{parent_chemid}/{plt_file}") as f:
                    lines = f.readlines()

                for lidx, line in enumerate(lines[4:]):
                    if re.search("^y[0-2] = \[", line):
                        y_data = line
                        if y_data.startswith("y0"):
                            level = "L1"
                        elif y_data.startswith("y1"):
                            level = "L2"
                        elif y_data.startswith("y2"):
                            level = "L3"
                    if re.search("^x =", line):
                        x_data = line
                    if "inf_energy" in line:
                        inf_energy = float(line.split("energy: ")[1])
                    if "scan_ref" in line:
                        corrections["1d"]["scan_ref"] = [[int(line.split('ref = ')[1].split(',')[0][1:]), int(line.split('ref = ')[1].split(',')[1][:-2])]]
                    if "VRC TST Sampling recommended start: " in line :
                        vrc_tst_start = float(line.split("start: ")[1])
                    if 'Number of unique sets of pivot points: ' in line:
                        for pp_number in range(int(line.split()[-1])):
                            indexes0 = [int(nmbr) for nmbr in lines[4:][lidx+4+6*pp_number].split('[[')[1].split('],')[0].split(',')]
                            indexes1 = [int(nmbr) for nmbr in lines[4:][lidx+4+6*pp_number].split('], [')[1].split(']]')[0].split(',')]
                            pivot_points_info.append({
                                '0': int(lines[4:][lidx+1+6*pp_number].split()[-1]),
                                '1': int(lines[4:][lidx+2+6*pp_number].split()[-1]),
                                'type': lines[4:][lidx+3+6*pp_number].split()[-1],
                                'indexes': [indexes0, indexes1],
                                'min': float(lines[4:][lidx+5+6*pp_number].split()[-1]),
                                'max': float(lines[4:][lidx+6+6*pp_number].split()[-1])
                            })

                match scan_type:
                    case "":
                        corrections["1d"]["e_trust"] = [float(value) for value in y_data.split("=")[1][2:-3].split(",")]
                        corrections["1d"]["r_trust"] = [float(value) for value in x_data.split("=")[1][2:-3].split(",")]
                    case "_frozen":
                        corrections["1d"]["e_sample"] = [float(value) for value in y_data.split("=")[1][2:-3].split(",")]
                        corrections["1d"]["r_sample"] = [float(value) for value in x_data.split("=")[1][2:-3].split(",")]
            else:
                corrections = {}
                    
        fragments = []
        for frag_number in range(tot_frag):
            chemid = products[frag_number]

            if par['high_level']:
                #Will read info from L2 structure
                basename = f"{chemid}_well_VTS_high"
            else:
                #Will read info from L1 structure
                basename = f"{chemid}_well_VTS"

            #Create ase.atoms objects for each fragments
            db = connect(f"{parent_chemid}/kinbot.db")
            *_, last_row = db.select(name=f"{basename}")
            atoms = last_row.toatoms() #This is an ase.atoms object
            fragments.append(Fragment.from_ase_atoms(atoms=atoms,
                                                    frag_number=frag_number,
                                                    max_frag=tot_frag,
                                                    chemid=chemid,
                                                    parent_chemid=parent_chemid,
                                                    par=par))
        if par['high_level']:
            #Will read info from L2 structure
            if '_IRC_' in reaction_name:
                basename = f"{reaction_name}_high"
            else:
                basename = f"{parent_chemid}_well_high"
        else:
            #Will read info from L1 structure
            if '_IRC_' in reaction_name:
                basename = f"{reaction_name}"
            else:
                basename = f"{parent_chemid}_well"

        db = connect(f"{parent_chemid}/kinbot.db")
        *_, last_row = db.select(name=f"{basename}", sort="-1")
        parent = StationaryPoint.from_ase_atoms(last_row.toatoms())#This is an ase.atoms object
        parent.characterize()

        fragnames = Fragment.get_fragnames()
        reactive_atoms = []
        reactive_atoms = pp_settings.get_ra(reac, parent_chemid)
        #if index < number_of_barrierless:
            #Map the long range fragments with the index of the parent to find 
        pp_settings.set_parent_map(parent, reactive_atoms, fragments)
            #Make sure the reactive atoms are in the same order as the fragments
            #Shouldn't be needed anymore
            #reactive_atoms = pp_settings.reset_reactive_atoms(reactive_atoms,[frag.map for frag in fragments])
        #endif

        #Set the pivot points on each fragments and create the surfaces
        surfaces = []
        comments = []
        if len(reactive_atoms) == 0:
            logger.warning("No reactive atom detected for this reaction. Pivot points on COMs.")

        for dist in par['vrc_tst_dist_list']:
            if dist < vrc_tst_start:
                logger.info(f"Removing sampling surface {dist} for reaction {reaction_name}")
                continue
            elif dist >= vrc_tst_start:
                surfaces.extend(pp_settings.create_all_surf_for_dist(dist, fragments, par, reactive_atoms))
                
        #Creating the strings to print input file
        #Fragments block:
        Fragments_block = "#Fragments geometries are in Angstroms\n"
        for frag in fragments:
            Fragments_block = Fragments_block + (repr(frag)) + "\n"

        #Surfaces block:
        Surfaces_block = "#Pivot_points and distances are in Bohr\ndivid_surf = [\n"
        for index, surf in enumerate(surfaces):
            Surfaces_block = Surfaces_block + (repr(surf))
            if surf == surfaces[-1]:
                Surfaces_block = Surfaces_block[:-3]
        Surfaces_block += "]\n"

        #Corrections block
        corrections_block = "#Scan energies are in Kcal/mol. The 0 is the assymptotic energy.\n#Scan distances are in Angstrom\ncorrections = {\n"
        for correction in corrections:
            corrections_block += f"'{correction}' : " + "{\n"
            for ckey, cvalue in corrections[correction].items():
                if isinstance(cvalue, str):
                    corrections_block += f"'{ckey}' : '{cvalue}',\n"
                else:
                    corrections_block += f"'{ckey}' : {cvalue},\n"
            corrections_block = corrections_block[:-2]
            corrections_block += "\n}\n"
        corrections_block += "}\n"

        #Calc_block:
        whoami = getpass.getuser()
        Calc_block = "calc = {\n" +\
                    "'code': 'molpro',\n" +\
                    f"'scratch': '/scratch/{whoami}',\n" +\
                    f"'processors': 1,\n" +\
                    f"'queue': '{par['queuing']}',\n" +\
                    "'max_jobs': 2000}"

        #Flux block:
        Flux_block = "flux_parameter = {'pot_smp_max': 6000, 'pot_smp_min': 500, #per facet\n" +\
                        "                  'tot_smp_max': 15000, 'tot_smp_min': 500, \n" +\
                        "                  'flux_rel_err': 5, 'smp_len': 1}\n"


        fname = f"{job_name}.py"
        folder = "rotdPy"
        template_file_path = f'{kb_path}/tpl/rotdPy.tpl'
        with open(template_file_path) as template_file:
            template = template_file.read()
            new_input = template.format(job_name = job_name,
                                   Fragments_block = Fragments_block,
                                   Surfaces_block = Surfaces_block,
                                   frag_names = '[' + ', '.join(fragnames) + ']',
                                   calc_block = Calc_block,
                                   flux_block = Flux_block,
                                   min_dist = par['vrc_tst_dist_list'][0],
                                   corrections_block=corrections_block,
                                   inf_energy=inf_energy)
        if not os.path.exists(folder):
            # Create a new directory because it does not exist
            os.makedirs(folder)

        with open(f"{folder}/{fname}", 'w') as f:
            f.write(new_input)

        #Erase the fragments for this reaction
        Fragment._instances = []

def is_unique_vdW(well, vdW):
    #Return boolean
    if 'prod' not in well:
        return True
    
    #Find well's reaction
    for idx, vdw in enumerate(vdW):
        vdw_name = vdw[1] + vdw[-1].split('vdW')[1]
        if vdw_name == well:
            products = '_'.join(sorted(vdw[2]))
            break

    #Compare products with other reactions
    other_vdW = vdW[:idx]
    if idx+1 < len(vdW):
        other_vdW.extend(vdW[idx+1:])
    for vdw in other_vdW[:idx]:
        other_prod = '_'.join(sorted(vdw[2]))
        if other_prod == products:
            return False
    #No other reaction with a vdW well lead to the same prod.
    return True
    
def find_min_vdW(vdW: list, well_energies: dict) -> dict:
    #Dict linking each vdW well to the lowest equivalent
    min_vdW = {}
    for idx, vdw in enumerate(vdW):
        vdw_name = vdw[1] + vdw[-1].split('vdW')[1]
        products = '_'.join(sorted(vdw[2]))
        vdw_energy = well_energies[vdw_name]

        #Minimum set to itself
        if vdw_name not in min_vdW:
            min_vdW[vdw_name] = vdw_name

        other_vdW = vdW[:idx]
        if idx+1 < len(vdW):
            other_vdW.extend(vdW[idx+1:])
        
        for other_vdw in other_vdW[:idx]:
            other_name = other_vdw[1] + other_vdw[-1].split('vdW')[1]
            other_prod = '_'.join(sorted(other_vdw[2]))
            other_energy = well_energies[other_name]

            #Minimum set to a different well if found
            if other_prod == products:
                if other_energy < well_energies[min_vdW[vdw_name]]:
                    min_vdW[vdw_name] = other_name
    return min_vdW
        
        

def create_pesviewer_input(par, wells, products, reactions, barrierless, vdW,
                           well_energies, prod_energies, highlight):
    """
    highlight: list of reaction names that need a red highlight
    """
    # delete the im_extent and xval files
    try:
        os.remove('{}_xval.txt'.format(par['title']))
    except OSError:
        pass

    if highlight is None:
        highlight = []

    #Dictionary
    #Key: 'vdW_name'
    #Value: 'min_vdW_name' 
    min_vdW = find_min_vdW(vdW, well_energies)

    well_lines = []
    for well in wells:
        if is_unique_vdW(well, vdW):
            energy = well_energies[well]
            well_lines.append('{} {:.2f}'.format(well, energy))
        else:
            energy = well_energies[min_vdW[well]]
            line = '{} {:.2f}'.format(well, energy)
            if line not in well_lines:
                well_lines.append(line)

    bimol_lines = []
    for prods in products:
        energy = prod_energies[prods]
        bimol_lines.append('{} {:.2f}'.format(prods, energy))

    ts_lines = []
    for rxn in reactions:
        high = ''
        if rxn[1] in highlight:
            high = 'red'
        prod_name = '_'.join(sorted(rxn[2]))
        ts_lines.append('{} {:.2f} {} {} {}'.format(rxn[1],
                                                    rxn[3],
                                                    rxn[0],
                                                    prod_name,
                                                    high))
    barrierless_lines = []
    prev_prod = []
    for index, rxn in enumerate(barrierless):
        prod_name = '_'.join(sorted(rxn[2]))
        barrierless_lines.append('{name} {react} {prod}'.format(name='nobar_' + str(index),
                                                                    react=rxn[0],
                                                                    prod=prod_name))
    
    for index, rxn in enumerate(vdW):
        high = ''
        if rxn[1] in highlight:
            high = 'red'
        vdW_name = f"{rxn[1]}{rxn[5].split('vdW')[1]}"
        prod_name = '_'.join(sorted(rxn[2]))
        
        ts_lines.append('{} {:.2f} {} {} {}'.format(rxn[1],
                                                    rxn[3],
                                                    rxn[0],
                                                    min_vdW[vdW_name],
                                                    high))
        line = '{name} {react} {prod}'.format(name='nobar_' + str(index + len(barrierless)),
                                                                react=min_vdW[vdW_name],
                                                                prod=prod_name)
        if line not in barrierless_lines:
            barrierless_lines.append(line)
            

    well_lines = '\n'.join(well_lines)
    bimol_lines = '\n'.join(bimol_lines)
    ts_lines = '\n'.join(ts_lines)
    barrierless_lines = '\n'.join(barrierless_lines)

    # write everything to a file
    fname = 'pesviewer.inp'
    template_file_path = f'{kb_path}/tpl/{fname}.tpl'
    with open(template_file_path) as template_file:
        template = template_file.read()
    template = template.format(id=par['title'],
                               wells=well_lines,
                               bimolecs=bimol_lines,
                               ts=ts_lines,
                               barrierless=barrierless_lines)
    with open(fname, 'w') as f:
        f.write(template)


def create_interactive_graph(wells, products, reactions, title, well_energies, prod_energies):
    """
    Create an interactive plot with pyvis
    """
    if len(wells) < 2:
        return -2
    try:
        from pyvis import network as net
    except ImportError:
        logger.warning('pyvis cannot be imported, no interactive plot is made.')
        return -1

    # For now we are assuming the all of the 2D depictions
    # are in place, which were created with PESViewer
    # Later we can add those in independently, but
    # this is not needed, just requires a quick run of
    # PESViewer

    conn, bars = get_connectivity(wells, products, reactions)

    g = net.Network(height='800px', width='50%',heading='')
    for i, well in enumerate(wells):
        g.add_node(i, label='', borderWidth=3, title=f'{well}: {round(well_energies[well], 1)} kcal/mol', 
                   shape='image', image=f'{os.getcwd()}/{title}/{well}_2d.png')
    for i, prod in enumerate(products):
        g.add_node(i + len(wells), label='', borderWidth=3, title=f'{prod}: {round(prod_energies[prod],1)} kcal/mol', 
                   shape='image', image=f'{os.getcwd()}/{title}/{prod}_2d.png')

    color_min = bars.min()
    color_max = bars.max()
    color_step = (color_max - color_min) / 256.
    for i, ci in enumerate(conn):
        for j, cij in enumerate(ci):
            if cij > 0:
                red = round((bars[i, j] - color_min) / color_step) 
                green = 0
                blue = round((color_max - bars[i, j]) / color_step)
                g.add_edge(i, j, title=f'{round(bars[i, j], 1)} kcal/mol', width=4, color=f'rgb({red},{green},{blue})')
    g.show_buttons(filter_=['physics'])
    g.save_graph(f'{title}.html')
    #display(HTML('example.html'))
    return 0


def get_energy(wells, job, ts, high_level, mp2=0, bls=0, conf=0):
    
    if ts:
        j = job
        wells = [job.split('_')[0]]
    else:
        j = job + '_well'
    if conf and not high_level and not mp2:
        j = f'conf/{job}_low'
    if mp2:
        j += '_mp2'
    if bls:
        j += '_bls'
    if high_level:
        j += '_high'
    energy = np.inf
    zpe = np.inf
    for well in wells:
        if "IRC" in well:
            well = well.split("_")[0]
        db = connect(well + '/kinbot.db')
        rows = db.select(name=j)
        for row in rows:
            try:
                new_energy = row.data.get('energy') * constants.EVtoHARTREE
                new_zpe = row.data.get('zpe')
            except (UnboundLocalError, TypeError):
                continue
            if new_zpe is None:
                continue
            if hasattr(row, 'data') and new_energy + new_zpe < energy + zpe:
                if not ts:
                    # Avoid getting energies from calculations that converged to another structure
                    atoms = row.toatoms()
                    st_pt = StationaryPoint.from_ase_atoms(atoms)
                    st_pt.characterize()
                    chemid_wo_mult = str(st_pt.chemid)[:-1]  # For charged species
                    if chemid_wo_mult != job[:-1]:
                        continue
                energy = new_energy
                zpe = new_zpe
    if np.isinf(energy) or np.isinf(zpe):
        raise ValueError(f'Unable to find an energy for {j}.')

    return energy, zpe


def get_l3energy(job, par, bls=0):
    """
    Get the L3, single-point energies.
    This is not object oriented.
    """

    if bls:
        key = par['barrierless_saddle_single_point_key']
    else:
        key = par['single_point_key']

    if job == '10000000000000000001':  # proton
        return 1, 0.0
    if par['single_point_qc'] == 'molpro':
        if os.path.exists(f'molpro/{job}.out'):
            with open(f'molpro/{job}.out', 'r') as f:
                lines = f.readlines()
            for line in reversed(lines):
            #for line in lines:
                if ('SETTING ' + key) in line:
                    e = float(line.split()[3])
                    logger.info('L3 electronic energy for {} is {} Hartree.'.format(job, e))
                    return 1, e  # energy was found
    elif par['single_point_qc'] == 'orca':
        if os.path.exists(f'orca/{job}_property.txt'):
            with open(f'orca/{job}_property.txt', 'r') as f:
                lines = f.readlines()
                for line in reversed(lines):
                    if (key) in line:
                        e = float(line.split()[-1])
                        logger.info('L3 electronic energy for {} is {} Hartree.'.format(job, e))
                        return 1, e  # energy was found
    elif par['single_point_qc'] == 'gauss':
        if os.path.exists('gauss/' + job + '.log'):
            gaussname = 'gauss/' + job + '.log'
        elif os.path.exists('gauss/' + job + '_high.log'):
            gaussname = 'gauss/' + job + '_high.log'
        elif os.path.exists('gauss/' + job + '_well_high.log'):
            gaussname = 'gauss/' + job + '_well_high.log'
        else:
            logger.info('L3 for {} is missing.'.format(job))
            return 0, -1  # job not yet started to run

        with open(gaussname) as f:
            lines = f.readlines()
            for line in reversed(lines):
                if (key) in line:
                    words = line.split()
                    wi = words.index(key) + 2
                    e = float(words[wi].replace('D', 'E'))
                    logger.info('L3 electronic energy for {} is {} Hartree.'.format(job, e))
                    return 1, e  # energy was found
 
    # if no file or no energy found
    logger.info('L3 for {} is missing.'.format(job))
    return 0, -1  # job not yet started to run or not finished


def get_zpe(jobdir, job, ts, high_level, mp2=0, bls=0):
    if "IRC" in job:
        jobdir = job.split("_")[0]
    db = connect(jobdir + '/kinbot.db')
    if ts:
        j = job
    else:
        j = job + '_well'
    if mp2:
        j += '_mp2'
    if bls:
        j += '_bls'
    if high_level:
        j += '_high'
    rows = db.select(name=j)
    zpe = None 
    for row in rows:
        if hasattr(row, 'data'):
            zpe = row.data.get('zpe')
    if zpe is None: 
        logger.warning('Could not find zpe for {} in directory {}'.format(job, jobdir))
        zpe = 1.  # a large value
    return zpe


def check_status(job, pid):
    command = ['ps', '-u', 'root', '-N', '-o', 'pid,s,user,%cpu,%mem,etime,args']
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stdin=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    out, err = process.communicate()
    out = out.decode()
    lines = out.split('\n')
    for line in lines:
        if len(line) > 0:
            if str(pid) == line.split()[0]:
                return 1
    return 0


def submit_job(chemid, par):
    """
    Submit a kinbot run using subprocess and return the pid
    """
    command = ["kinbot", chemid + ".json", "&"]
    # purge previous summary and monitor files, so that pes doesn't think
    # everything is done
    # relevant if jobs are killed
    try:
        os.system(f'rm -f {chemid}/summary_*.out')
    except OSError:
        pass
    try:
        os.system(f'rm -f {chemid}/kinbot_monitor.out')
    except OSError:
        pass

    for tmpl in ['queue_template', 'q_temp_am1', 'q_temp_mp2', 'q_temp_hi', 
                 'q_temp_l3']:
        if par[tmpl] != '':
            shutil.copyfile('{}'.format(par[tmpl]), 
                            '{}/{}'.format(chemid, par[tmpl]))
    if par['single_point_template'] != '':
        shutil.copyfile('{}'.format(par['single_point_template']), 
                        '{}/{}'.format(chemid, par['single_point_template']))
    if par['barrierless_saddle_single_point_template'] != '':
        shutil.copyfile('{}'.format(par['barrierless_saddle_single_point_template']), '{}/{}'
                        .format(chemid, par['barrierless_saddle_single_point_template']))
        shutil.copyfile('{}'.format(par['barrierless_saddle_prod_single_point_template']), '{}/{}'
                        .format(chemid, par['barrierless_saddle_prod_single_point_template']))
    if par['vrc_tst_scan_parameters']['molpro_tpl'] != '':
        shutil.copyfile('{}'.format(par['vrc_tst_scan_parameters']['molpro_tpl']), 
                        '{}/{}'.format(chemid, par['vrc_tst_scan_parameters']['molpro_tpl']))
    outfile = open(f'{chemid}/kinbot.out', 'w')
    errfile = open(f'{chemid}/kinbot.err', 'w')
    process = subprocess.Popen(command,
                               cwd=chemid,
                               stdout=outfile,
                               stdin=subprocess.PIPE,
                               stderr=errfile)
    time.sleep(1)
    pid = process.pid
    return pid


def write_input(input_file, species, threshold, threshold_L2, root, me):
    # directory for this particular species
    directory = root + '/' + str(species.chemid) + '/'
    if not os.path.exists(directory):
        os.makedirs(directory)

    # make a new parameters instance and overwrite some keys
    input_file = '{}'.format(input_file)
    par2 = Parameters(input_file).par
    # overwrite the title
    par2['title'] = str(species.chemid)
    # make a structure vector and overwrite the par structure
    structure = []
    for at in range(species.natom):
        pos = species.geom[at]
        sym = species.atom[at]
        structure += [sym, pos[0], pos[1], pos[2]]
    par2['structure'] = structure
    # delete the par smiles
    par2['smiles'] = ''
    # overwrite the barrier threshold
    par2['barrier_threshold'] = threshold
    # overwrite the barrier threshold for L2
    par2['barrier_threshold_L2'] = threshold_L2
    # set the pes option to 1
    par2['pes'] = 1
    # don't do ME for these kinbots but write the files
    if me:
        par2['me'] = 2

    file_name = directory + str(species.chemid) + '.json'
    with open(file_name, 'w') as outfile:
        json.dump(par2, outfile, indent=4, sort_keys=True)
    return


def write_input_keep(input_file, keepchemid, root):
    directory = root + '/' + str(keepchemid) + '/'
    if not os.path.exists(directory):
        print(f'Cannot keep a well that is not already explored. Please correct the keep_chemids parameter. Bye!')
        sys.exit(-1)

    par_keep = Parameters(f'{directory}{keepchemid}.json').par
    # make a new parameters instance and overwrite some keys
    input_file = '{}'.format(input_file)
    par_new = Parameters(input_file).par
    # overwrite the title
    par_new['title'] = par_keep['title']
    par_new['structure'] = par_keep['structure']
    par_new['barrier_threshold'] = par_keep['barrier_threshold']
    par_new['barrier_threshold_L2'] = par_keep['barrier_threshold_L2']
    par_new['pes'] = 1
    par_new['me'] = 2

    file_name = directory + str(keepchemid) + '.json'
    with open(file_name, 'w') as outfile:
        json.dump(par_new, outfile, indent=4, sort_keys=True)
    return


def check_l3_l2(l3_key: str, parent_specs: dict, reactions: list) -> None:
    """Perform a check on the difference between L3 and L2 energy differences.

    @param l3_key: Pattern to search for in L3 calculations.
    @param parent_specs: Dictionary of species' chemids with its parent species.
    @param reactions: List of all reaction names in the PES.
    """
    # Get L3 energies
    l3_energies = {}
    logger.info(f'L3-L2 Energy difference Analysis. Energy units: kcal/mol.')
    if not os.path.isdir('molpro'):
        logger.warning("Unable to perform L3-L2 check. The molpro directory "
                        "is missing.")
        return
    if len([f for f in os.listdir('molpro/') if f.endswith('.out')]) == 0:
        logger.warning("Unable to perform L3-L2 check. The molpro directory "
                        "is empty.")
        return

    # Get L3 energies
    for st_pt in list(parent_specs.keys()) + [r[1] for r in reactions]:
        if "_" in st_pt and all([frag.isdigit() for frag in st_pt.split('_')]):
            # Bimolecular species
            l3_energies[st_pt] = 0
            for frag in st_pt.split('_'):
                if not os.path.isfile(f'molpro/{frag}.out'):
                    if st_pt in l3_energies:
                        del l3_energies[st_pt]
                    break
                with open(f'molpro/{frag}.out') as out_fh:
                    for line in out_fh:
                        if l3_key not in line:
                            continue
                        l3_energies[st_pt] += float(line.split()[3])
                        break
        elif os.path.isfile(f'molpro/{st_pt}.out'):
            # Wells and TSs
            with open(f'molpro/{st_pt}.out') as out_fh:
                for line in out_fh:
                    if l3_key not in line:
                        continue
                    l3_energies[st_pt] = float(line.split()[3])
                    break

    # Get L2 Energies and its difference respect L3.
    e_diffs = {}
    for st_pt in l3_energies:
        if any([c.isalpha() for c in st_pt]):  # TSs (have letters in the name)
            db_path = f'{st_pt.split("_")[0]}/kinbot.db'
        else:  # Wells and Bimolecular products
            db_path = f'{parent_specs[st_pt]}/kinbot.db'
        if not os.path.isfile(db_path):
            logger.warning(f"Unable to find L2 energy for {st_pt}.")
            continue
        db = connect(db_path)
        if st_pt.isdigit():
            # Wells
            rows = db.select(name=f'{st_pt}_well_high')
        elif all([fr.isdigit() for fr in st_pt.split('_')]):
            # Bimolecular species.
            l2_energy = 0
            for frag in st_pt.split('_'):
                rows = db.select(name=f'{frag}_well_high')
                try:
                    final_row = next(rows)
                except StopIteration:
                    logger.warning(f"Unable to find L2 energy for {frag}.")
                    l2_energy = np.nan
                    break
                for row in rows:
                    final_row = row
                l2_energy += final_row.data["energy"] * constants.EVtoHARTREE
        else:
            rows = db.select(name=f'{st_pt}_high')

        if "_" not in st_pt or any([c.isalpha() for c in st_pt]):
            # Wells and TSs
            try:
                final_row = next(rows)
            except StopIteration:
                logger.warning(f"Unable to find L2 energy for {st_pt}.")
                continue
            for row in rows:
                final_row = row
            l2_energy = final_row.data["energy"] * constants.EVtoHARTREE
        e_diff = l3_energies[st_pt] - l2_energy
        e_diffs[st_pt] = e_diff / constants.KCALtoHARTREE

    e_diff_avg = np.round(np.average(list(e_diffs.values())), 1)
    e_diff_std = np.round(np.std(list(e_diffs.values())), 1)
    logger.info(f'Avg difference: {e_diff_avg} kcal/mol, '
                 f'Max: {np.round(max(e_diffs.values()), 1)} kcal/mol, '
                 f'Min: {np.round(min(e_diffs.values()), 1)} kcal/mol, '
                 f'STDEV: {e_diff_std} kcal/mol.')
    for st_pt, e_diff in e_diffs.items():
        if not e_diff_avg * 0.9 < e_diff < e_diff_avg * 1.1:
            logger.info(f"Outlying L2-L3 difference found for {st_pt}. "
                         f"Energy difference: {np.round(e_diff, 1)} kcal/mol.")


def t1_analysis(lot='TZ'):
    import os
    try:
        import matplotlib.pyplot as plt
        do_plot = True
    except ModuleNotFoundError:
        do_plot = False

    if 'TZ' in lot.upper():
        lot = 'TZ'
    elif 'DZ' in lot.upper():
        lot = 'DZ'
    else:
        logger.warning('Unable to perform a summary of T1 diagnostics: '
                        'Unrecognized single_point_key.')
        return
    T1s = []
    for f in os.listdir('molpro/'):
        if not f.endswith('.out'):
            continue
        with open(f'molpro/{f}') as fh:
            do_read1 = False
            do_read2 = False
            for line in fh:
                if lot in line:
                    do_read1 = True
                    do_read2 = False
                if "Starting UCCSD calculation" in line or "Starting RCCSD calculation" in line:
                    do_read2 = True
                elif do_read1 and do_read2 and 'T1 diagnostic:' in line:
                    T1s.append(float(line.split()[9]))
                    break
    if T1s:
        counts, bins = np.histogram(T1s)
    else:
        logger.warning('Unable to perform a summary of T1 diagnostics: '
                        'No T1 Diagnostics results found.')
        return

    if do_plot:
        fig1, ax1 = plt.subplots()  # Histogram
        ax1.bar(bins[1:], counts, width=bins[1]*0.7)
        ax1.set_xlabel('T1 Diagnostic')
        ax1.set_ylabel('Counts')
        fig1.savefig('T1_histo.png')

        fig2, ax2 = plt.subplots()  # Cumulative counts.
        ax2.plot(sorted(T1s), range(len(T1s)))
        ax2.set_xlabel('T1 Diagnostic')
        ax2.set_ylabel('Cumulative counts')
        fig2.savefig('T1_cumul_count.png')
    else:
        logger.warning('Matplotlib not found. Unable to plot T1 diagnostics '
                        'summary.')

    logger.info(f"T1 histogram:\nBins: {bins}\nCounts: {counts}.")


if __name__ == "__main__":
    main()
