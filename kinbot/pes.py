"""
This is the main class to run KinBot to explore
a full PES instead of only the reactions of one well
"""
import sys
import os
import stat
import shutil
import logging
import datetime
import time
import subprocess
import json
from distutils.dir_util import copy_tree
import pkg_resources
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from copy import deepcopy

from ase.db import connect

from kinbot import constants
from kinbot import license_message
from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint
from kinbot.mess import MESS
from kinbot.uncertaintyAnalysis import UQ


def main():
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

    # initialize the parameters
    par = Parameters(input_file).par

    # set up the logging environment
    logging.basicConfig(filename='pes.log', level=logging.INFO)

    logging.info(license_message.message)
    msg = 'Starting the PES search at {}'.format(datetime.datetime.now())
    logging.info(msg)

    well0 = StationaryPoint('well0',
                            par['charge'],
                            par['mult'],
                            smiles=par['smiles'],
                            structure=par['structure'])
    well0.characterize()
    write_input(input_file, well0, par['barrier_threshold'], os.getcwd())

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
    # list of all jobs
    jobs = []
    # dict of the pid's for all jobs
    pids = {}
    a = 0
    b = 0
    c = 0
    while 1:
        j = len(jobs)
        if j != a:
            logging.info('{0} {1} {2}'.format("len(jobs): ", j, "\n"))
        if j != a:
            logging.info('{0} {1} {2}'.format("len(jobs): ", j, "\n"))
        a = j
        with open('chemids', 'r') as f:
            jobs = f.read().split('\n')
            jobs = [ji for ji in jobs if ji != '']

        if len(jobs) > j:
            logging.info('\tPicked up new jobs: ' + ' '.join(jobs[j:]))

        k = len(running)
        l = len(finished)
        if b != k:
            logging.info('{0} {1} {2}'.format("len(running): ", len(running), "\n"))
        b = k
        if c != l:
            logging.info('{0} {1} {2}'.format("len(finished): ", len(finished), "\n"))
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
            logging.info('Job: {}'.format(job))
            if 'none' in par['skip_chemids']:
                logging.info('No KinBot runs to be skipped')
            else:
                if job in par['skip_chemids']:
                    kb = 0
            logging.info('kb: {}'.format(kb))
            if kb == 1:
                pid = 0
                if not no_kinbot:
                    pid = submit_job(job, par)  # kinbot is submitted here
                else:
                    get_wells(job)
                pids[job] = pid
                t = datetime.datetime.now()
                logging.info('\tStarted job {} at {}'.format(job, t))
                running.append(job)
            elif kb == 0:
                logging.info('Skipping Kinbot for {}'.format(job))
                finished.append(job)
            else:
                logging.info('kb value not 0 or 1')

        # check if a thread is done
        for job in running:
            if not check_status(job, pids[job]):
                t = datetime.datetime.now()
                logging.info('\tFinished job {} at {}'.format(job, t))
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
            with open('pes_summary.txt', 'w') as f:
                f.write('\n'.join(summary_lines))
            time.sleep(1)

    # delete skipped jobs from the jobs before sending to postprocess
    for skip in par['skip_chemids']:
        try:
            jobs.pop(jobs.index(skip))
        except ValueError:
            pass

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
    logging.info('PES search done!')
    print('PES search done!')


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
            pieces = line.split()
            if line.startswith('SUCCESS'):
                prod = pieces[3:]
            elif line.startswith('HOMOLYTIC_SCISSION'):
                prod = pieces[2:]
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
    jobs: all of the jobs that were run
    temp: this is a temporary output file writing
    """

    l3done = 1  # flag for L3 calculations to be complete

    # base of the energy is the first well, these are L2 energies
    base_energy = get_energy(jobs[0], jobs[0], 0, par['high_level'])
    # L3 energies
    status, base_l3energy = get_l3energy(jobs[0], par)
    if not status:
        l3done = 0
    # L2 ZPE
    base_zpe = get_zpe(jobs[0], jobs[0], 0, par['high_level'])
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
    products = []

    # read all the jobs
    for ji in jobs:
        try:
            summary = open(ji + '/summary_' + ji + '.out', 'r').readlines()
        except:
            failedwells.append(ji)
            continue
        # read the summary file
        for line in summary:
            if line.startswith('SUCCESS'):
                pieces = line.split()
                ts = pieces[2]  # this is the long specific name of the reaction
                if ts in par['skip_reactions']:
                    continue
                reactant = ji
                prod = pieces[3:]  # this is the chemid of the product
                # calculate the barrier based on the new energy base
                barrier = 0. - base_energy - base_zpe

                # overwrite energies with mp2 energy if needed
                mp2_list = ['R_Addition_MultipleBond', 'reac_birad_recombination_R', 
                        'reac_r12_cycloaddition', 'reac_r14_birad_scission']
                if (any([mm in ts for mm in mp2_list]) and not par['high_level']):
                    base_energy_mp2 = get_energy(jobs[0], jobs[0], 0, 
                                                 par['high_level'], mp2=1)
                    base_zpe_mp2 = get_zpe(jobs[0], jobs[0], 0, 
                                           par['high_level'], mp2=1)
                    barrier = 0. - base_energy_mp2 - base_zpe_mp2

                # overwrite energies with bls energy if needed
                if 'barrierless_saddle' in ts and not par[ 'high_level']:
                    base_energy_bls = get_energy(jobs[0], jobs[0], 0,
                                                 par['high_level'], bls=1)
                    base_zpe_bls = get_zpe(jobs[0], jobs[0], 0,
                                           par['high_level'], bls=1)
                    barrier = 0. - base_energy_bls - base_zpe_bls

                ts_energy = get_energy(reactant, ts, 1, par['high_level'])
                ts_zpe = get_zpe(reactant, ts, 1, par['high_level'])
                barrier += ts_energy + ts_zpe
                barrier *= constants.AUtoKCAL
                if reactant not in wells:
                    wells.append(reactant)
                    parent[reactant] = reactant
                if len(prod) == 1:
                    if prod[0] not in wells:
                        if prod[0] not in parent:
                            parent[prod[0]] = reactant
                        wells.append(prod[0])
                else:
                    prod_name = '_'.join(sorted(prod))
                    if prod_name not in products:
                        if prod_name not in parent:
                            parent[prod_name] = reactant
                        products.append('_'.join(sorted(prod)))
                new = 1
                temp = None

                for i, rxn in enumerate(reactions):
                    rxn_prod_name = '_'.join(sorted(rxn[2]))
                    if (reactant == rxn[0] and
                            '_'.join(sorted(prod)) == rxn_prod_name):
                        new = 0
                        temp = i
                    if reactant == ''.join(rxn[2]) and ''.join(prod) == rxn[0]:
                        new = 0
                        temp = i
                if new:
                    reactions.append([reactant, ts, prod, barrier])
                else:
                    # check if the previous reaction has a lower energy or not
                    if reactions[temp][3] > barrier:
                        reactions.pop(temp)
                        reactions.append([reactant, ts, prod, barrier])

            elif line.startswith('HOMOLYTIC_SCISSION'):
                pieces = line.split()
                reactant = ji
                energy = pieces[1]  # energy from summary
                prod = pieces[2:]   # this is the chemid of the product

                if reactant not in wells:
                    wells.append(reactant)
                    parent[reactant] = reactant

                if len(prod) == 1:
                    if prod[0] not in wells:
                        if prod[0] not in parent:
                            parent[prod[0]] = reactant
                        wells.append(prod[0])
                else:
                    prod_name = '_'.join(sorted(prod))
                    if prod_name not in products:
                        if prod_name not in parent:
                            parent[prod_name] = reactant
                        products.append('_'.join(sorted(prod)))
                new = 1
                temp = None
                for i, rxn in enumerate(reactions):
                    rxn_prod_name = '_'.join(sorted(rxn[2]))
                    if (reactant == rxn[0] and
                            '_'.join(sorted(prod)) == rxn_prod_name):
                        new = 0
                        temp = i
                    if reactant == ''.join(rxn[2]) and ''.join(prod) == rxn[0]:
                        new = 0
                        temp = i
                if new:
                    ts = 'barrierless'
                    barrier = energy
                    reactions.append([reactant, ts, prod, barrier])
        # copy the xyz files
        copy_from_kinbot(ji, 'xyz')
        # copy the L3 calculations here, whatever was in those directories, inp, out, pbs, etc.
        copy_from_kinbot(ji, par['single_point_qc'])
    # create a connectivity matrix for all wells and products
    conn, bars = get_connectivity(wells, products, reactions)
    # create a batch submission for all L3 jobs
    # TODO slurm
    if par['queuing'] == 'pbs':
        batch = 'batch_L3_pbs.sub'
        with open(batch, 'w') as f:
            for well in wells:
                f.write('qsub molpro/' + well + '.pbs' + '\n')
            for prod in products:
                for frag in prod.split('_'):
                    f.write('qsub molpro/' + frag + '.pbs' + '\n')
            for reac in reactions:
                f.write('qsub molpro/' + reac[1] + '.pbs' + '\n')
        os.chmod(batch, stat.S_IRWXU)  # read, write, execute by owner

    well_energies = {}
    well_l3energies = {}
    for well in wells:
        energy = get_energy(parent[well], well, 0, par['high_level'])  # from the db
        zpe = get_zpe(parent[well], well, 0, par['high_level'])
        well_energies[well] = ((energy + zpe) - (base_energy + base_zpe)) * constants.AUtoKCAL
        status, l3energy = get_l3energy(well, par)
        if not status:
            l3done = 0  # not all L3 calculations are done
        else:
            well_l3energies[well] = ((l3energy + zpe) - (base_l3energy + base_zpe)) * constants.AUtoKCAL
    prod_energies = {}
    prod_l3energies = {}
    for prods in products:
        energy = 0. - (base_energy + base_zpe)
        l3energy = 0. - (base_l3energy + base_zpe)
        for pr in prods.split('_'):
            energy += get_energy(parent[prods], pr, 0, par['high_level'])
            zpe = get_zpe(parent[prods], pr, 0, par['high_level'])
            energy += zpe
            status, l3e = get_l3energy(pr, par)
            if not status:
                l3done = 0  # not all L3 calculations are done
            else:
                l3energy += l3e + zpe
        prod_energies[prods] = energy * constants.AUtoKCAL
        prod_l3energies[prods] = l3energy * constants.AUtoKCAL

    ts_l3energies = {}
    for reac in reactions:
        if reac[1] != 'barrierless':  # meaning homolytic scission
            if 'barrierless_saddle' in reac[1]:
                zpe = get_zpe(reac[0], reac[1], 1, par['high_level'])
                status, l3energy = get_l3energy(reac[1], par, bls=1)

                status_prod1, l3energy_prod1 = get_l3energy(reac[2][0], par)
                status_prod2, l3energy_prod2 = get_l3energy(reac[2][1], par)

                status_prod, l3energy_prod = get_l3energy(reac[1] + '_prod', par, bls=1)

                if not status * status_prod:
                    l3done = 0
                else:
                    delta1 = l3energy_prod - (l3energy + zpe)  # ZPEs cancel out for fragments
                    delta2 = l3energy_prod1 + l3energy_prod2 - (base_l3energy + base_zpe) 
                    ts_l3energies[reac[1]] = (delta2 - delta1) * constants.AUtoKCAL
            else:
                zpe = get_zpe(reac[0], reac[1], 1, par['high_level'])
                status, l3energy = get_l3energy(reac[1], par)
                if not status:
                    l3done = 0
                else:
                    ts_l3energies[reac[1]] = ((l3energy + zpe) - (base_l3energy + base_zpe)) * constants.AUtoKCAL

    logging.info('l3done status {}'.format(l3done))

    if l3done == 1:
        logging.info('Energies are updated to L3 in ME and PESViewer.')
        well_energies = well_l3energies
        prod_energies = prod_l3energies
        for reac in reactions:  # swap out the barrier
            if 'barrierless' not in reac[1]:
                reac[3] = ts_l3energies[reac[1]]
            if 'barrierless_saddle' in reac[1]:
                reac[3] = ts_l3energies[reac[1]]

        logging.info('L3 energies in kcal/mol, incl. ZPE')
        for well in wells:
            logging.info('{}   {:.2f}'.format(well, well_l3energies[well]))
        for prod in products:
            logging.info('{}   {:.2f}'.format(prod, prod_l3energies[prod]))
        for ts in ts_l3energies:
            logging.info('{}   {:.2f}'.format(ts, ts_l3energies[ts]))


    # if L3 was done, everything below is done with that
    # filter according to tasks
    wells, products, reactions, highlight = filter(par,
                                                   wells,
                                                   products,
                                                   reactions,
                                                   conn,
                                                   bars,
                                                   well_energies,
                                                   task,
                                                   names)

    # draw a graph of the network
    create_graph(wells,
                 products,
                 reactions,
                 well_energies,
                 prod_energies,
                 highlight)

    create_interactive_graph(wells,
                             products,
                             reactions,
                             par['title'],
                             well_energies,
                             prod_energies,
                             )

    barrierless = []
    rxns = []
    for rxn in reactions:
        if rxn[1] == 'barrierless':
            barrierless.append([rxn[0], rxn[1], rxn[2], rxn[3]])
        else:
            rxns.append([rxn[0], rxn[1], rxn[2], rxn[3]])

    # write full pesviewer input
    create_pesviewer_input(par,
                           wells,
                           products,
                           rxns,
                           barrierless,
                           well_energies,
                           prod_energies,
                           highlight)

    create_mess_input(par,
                      wells,
                      products,
                      rxns,
                      barrierless,
                      well_energies,
                      prod_energies,
                      parent,
                      mass)


def filter(par, wells, products, reactions, conn, bars, well_energies, task, names):
    """
    Filter the wells, products and reactions according to the task
    and the names
    """
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
        pass
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
        for list in all_rxns:
            for rxn in list:
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
            logging.error('Only one name should be given for a well filter')
            logging.error('Received: ' + ' '.join(names))
            sys.exit(-1)
    elif task == 'temperature':
        if len(names) == 1:
            try:
                # read the temperature
                temperature = float(names[0])
            except ValueError:
                logging.error('A float is needed for a temperature filter')
                logging.error('Received: ' + ' '.join(names))
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
            logging.error('Only one argument should be given for a temperature filter')
            logging.error('Received: ' + ' '.join(names))
            sys.exit(-1)
    elif task == 'l2threshold':
        filtered_reactions = []
        for rxn in reactions:
            if rxn[3] < par['barrier_threshold']:
                filtered_reactions.append(rxn)
    else:
        logging.error('Could not recognize task ' + task)
        sys.exit(-1)

    # filter the wells
    filtered_wells = []
    for well in wells:
        for rxn in filtered_reactions:
            prod_name = '_'.join(sorted(rxn[2]))
            if well == rxn[0] or well == prod_name:
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
        logging.error('Cannot find a lowest path if the number of species is not 2')
        logging.error('Found species: ' + ' '.join(names))


def get_index(wells, products, name):
    try:
        i = wells.index(name)
    except ValueError:
        try:
            i = products.index(name) + len(wells)
        except ValueError:
            logging.error('Could not find reactant ' + name)
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
    dirname = dirname + '/'
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    copy_tree(well + '/' + dirname, dirname)


def get_rxn(prods, rxns):
    for rxn in rxns:
        if prods == '_'.join(sorted(rxn[2])):
            return rxn


def create_short_names(wells, products, reactions, barrierless):
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
        nobar_short[nobar_name] = nobar_name
        n = n + 1

    return well_short, pr_short, fr_short, ts_short, nobar_short


def create_mess_input(par, wells, products, reactions, barrierless,
                      well_energies, prod_energies, parent, mass):
    """
    When calculating a full pes, the files from the separate wells
    are read and concatenated into one file
    Two things per file need to be updated
    1. the names of all the wells, bimolecular products and ts's
    2. all the zpe corrected energies
    """

    logging.info(f"uq value: {par['uq']}")
    well_short, pr_short, fr_short, ts_short, nobar_short = create_short_names(wells,
                                                                               products,
                                                                               reactions,
                                                                               barrierless)

    # list of the strings to write to mess input file
    s = []

    """
    Create the header block for MESS
    """
    frame = '######################\n' 
    divider = '!****************************************\n'
    
    dummy = StationaryPoint('dummy',
                            par['charge'],
                            par['mult'],
                            smiles=par['smiles'],
                            structure=par['structure'])

    mess = MESS(par, dummy)
    uq = UQ(par)

    header_file = pkg_resources.resource_filename('tpl', 'mess_header.tpl')
    with open(header_file) as f:
        tpl = f.read()

    for uq_iter in range(par['uq_n']):
        mess_iter = "{0:04d}".format(uq_iter)

        e_well = par['epsilon'] * uq.calc_factor('epsilon', '', uq_iter)
        s_well = par['sigma'] * uq.calc_factor('sigma', '', uq_iter)
        enrelfact = par['EnergyRelaxationFactor'] * uq.calc_factor('enrelfact', '', uq_iter)
        enrelpow = par['EnergyRelaxationPower'] * uq.calc_factor('enrelpow', '', uq_iter)

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
                            )

        s = []  # mess string after header

        well_energies_current = deepcopy(well_energies)
        prod_energies_current = deepcopy(prod_energies)
            
        # write the wells
        s.append(frame + '# WELLS\n' + frame)
        for well in wells:
            name = well_short[well] + ' ! ' + well
            energy = well_energies[well] + uq.calc_factor('energy', well_short[well], uq_iter)
            well_energies_current[well] = energy
            with open(parent[well] + '/' + well + '_' + mess_iter + '.mess', 'r') as f:
                s.append(f.read().format(name=name, zeroenergy=round(energy, 2)))
            s.append(divider)

        # write the products
        s.append(frame + '# BIMOLECULAR PRODUCTS\n' + frame)
        for prod in products:
            name = pr_short[prod] + ' ! ' + prod
            energy = prod_energies[prod] + uq.calc_factor('energy', pr_short[prod], uq_iter)
            prod_energies_current[prod] = energy
            fr_names = {}
            for fr in prod.split('_'):
                key = 'fr_name_{}'.format(fr)
                value = fr_short[fr] + ' ! ' + fr
                fr_names[key] = value
            with open(parent[prod] + '/' + prod + '_' + mess_iter + '.mess') as f:
                s.append(f.read().format(name=name,
                                         ground_energy=round(energy, 2),
                                         **fr_names))
            s.append(divider)

        # write the barrier
        s.append(frame + '# BARRIERS\n' + frame)
        for rxn in reactions:
            with open("reactionList.log", 'a') as f:
                f.write('{0} {1}'.format(rxn, "\n"))
            name = [ts_short[rxn[1]]]
            name.append(well_short[rxn[0]])
            if len(rxn[2]) == 1:
                name.append(well_short[rxn[2][0]])
            else:
                name.append(pr_short['_'.join(sorted(rxn[2]))])
            name.append('!')
            name.append(rxn[1])
            energy = rxn[3] + uq.calc_factor('barrier', ts_short[rxn[1]], uq_iter)
            welldepth1 = energy - well_energies_current[rxn[0]] 
            if len(rxn[2]) == 1:
                welldepth2 = energy - well_energies_current[rxn[2][0]] 
            else:
                prodname = '_'.join(sorted(rxn[2]))
                welldepth2 = energy - prod_energies_current[prodname] 
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

        with open('me/' + 'mess_' + mess_iter + '.inp', 'w') as f:
            f.write(header)
            f.write('\n'.join(s))

        if par['me']:
            mess.run()

        #uq.format_uqtk_data() 
    return

def create_pesviewer_input(par, wells, products, reactions, barrierless,
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

    well_lines = []
    for well in wells:
        energy = well_energies[well]
        well_lines.append('{} {:.2f}'.format(well, energy))

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
    index = 0
    new = 1
    prev_prod = []
    for rxn in barrierless:
        prod_name = '_'.join(sorted(rxn[2]))
        for item in prev_prod:
            if prod_name == item:
                new = 0
                break
        if new:
            barrierless_lines.append('{name} {react} {prod}'.format(name='nobar_' + str(index),
                                                                    react=rxn[0],
                                                                    prod=prod_name))
            prev_prod.append(prod_name)
            index = index + 1

    well_lines = '\n'.join(well_lines)
    bimol_lines = '\n'.join(bimol_lines)
    ts_lines = '\n'.join(ts_lines)
    barrierless_lines = '\n'.join(barrierless_lines)

    # write everything to a file
    fname = 'pesviewer.inp'
    template_file_path = pkg_resources.resource_filename('tpl', fname + '.tpl')
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
    try:
        from pyvis import network as net
    except ImportError:
        logging.warning('pyvis cannot be imported, no interactive plot is made.')
        return -1
    try:
        from IPython.core.display import display, HTML
    except ImportError:
        logging.warning('IPython cannot be imported, no interactive plot is made.')
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


def create_graph(wells, products, reactions,
                 well_energies, prod_energies, highlight):
    """
    highlight: list of reaction names that need a red highlight
    """
    if highlight is None:
        highlight = []
    # update the connectivity with the filtered wells, products and reactions
    conn, bars = get_connectivity(wells, products, reactions)

    # get the minimum and maximum well and product energy
    try:
        minimum = min(min(well_energies.values()),
                      min(prod_energies.values()))
        maximum = max(max(well_energies.values()),
                      max(prod_energies.values()))
    except ValueError:
        # list of products can be empty, but list of wells not
        minimum = min(well_energies.values())
        maximum = max(well_energies.values())
    # define the inveresly proportial weights function
    max_size = 400
    min_size = 100
    slope = (min_size - max_size) / (maximum - minimum)
    offset = max_size - minimum * slope
    # define the graph nodes
    nodes = [i for i, wi in enumerate(wells)]
    nodes += [len(wells) + i for i, pi in enumerate(products)]
    # size of the nodes from the weights
    node_size = [slope * well_energies[wi] + offset for wi in wells]
    node_size += [slope * prod_energies[pi] + offset for pi in products]
    # color nodes and wells differently
    node_color = ['lightskyblue' for wi in wells]
    node_color += ['lightcoral' for pi in products]
    # labels of the wells and products
    labels = {}
    name_dict = {}
    for i, wi in enumerate(wells):
        labels[i] = 'w{}'.format(i + 1)
        name_dict[labels[i]] = wi
    for i, pi in enumerate(products):
        labels[i + len(wells)] = 'b{}'.format(i + 1)
        name_dict[labels[i + len(wells)]] = pi
    # write the labels to a file
    with open('species_dict.txt', 'w') as f:
        lines = []
        for name in sorted(name_dict.keys()):
            lines.append('{}  {}'.format(name, name_dict[name]))
        f.write('\n'.join(lines))
    # make a graph object
    G = nx.Graph()
    # add the nodes
    for i, node in enumerate(nodes):
        G.add_node(node, weight=node_size[i])

    # define the inversely proportional weights for the lines
    minimum = min(rxn[3] for rxn in reactions)
    maximum = max(rxn[3] for rxn in reactions)
    max_size = 5
    min_size = 0.5
    try:
        slope = (min_size - max_size) / (maximum - minimum)
    except:
        slope = 1.
    offset = max_size - minimum * slope

    # add the edges
    for i, ci in enumerate(conn):
        for j, cij in enumerate(ci):
            if cij > 0:
                weight = slope * bars[i][j] + offset
                G.add_edge(i, j, weight=weight)
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]

    # position the nodes
    pos = nx.spring_layout(G, scale=1)
    #pos = nx.circular_layout(G)

    
    # make the matplotlib figure
    plt.figure(figsize=(8, 8))
    nx.draw_networkx_edges(G, pos, edgelist=edges, width=weights)
    nx.draw_networkx_nodes(G,
                           pos,
                           nodelist=G.nodes(),
                           node_size=node_size,
                           node_color=node_color)
    nx.draw_networkx_labels(G, pos, labels, font_size=8)
    plt.axis('off')
    plt.savefig('graph.png')
    

def get_energy(directory, job, ts, high_level, mp2=0, bls=0):
    db = connect(directory + '/kinbot.db')
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
    for row in rows:
        if hasattr(row, 'data'):
            energy = row.data.get('energy')
    try:
        # ase energies are always in ev, convert to hartree
        energy *= constants.EVtoHARTREE
    except UnboundLocalError or TypeError:
        # this happens when the job is not found in the database
        logging.error('Could not find {} in directory {} database.'.format(job, directory))
        logging.error('Exiting...')
        sys.exit(-1)
    except TypeError:
        logging.warning('Could not find {} in directory {}'.format(job, directory))
        energy = 0.
    return energy


def get_l3energy(job, par, bls=0):
    """
    Get the L3, single-point energies.
    This is not object oriented.
    """

    if bls:
        key = par['barrierless_saddle_single_point_key']
    else:
        key = par['single_point_key']

    if par['single_point_qc'] == 'molpro':
        if os.path.exists('molpro/' + job + '.out'):
            with open('molpro/' + job + '.out', 'r') as f:
                lines = f.readlines()
                for line in reversed(lines):
                    if ('SETTING ' + key) in line:
                        e = float(line.split()[3])
                        logging.info('L3 electronic energy for {} is {} Hartree.'.format(job, e))
                        return 1, e  # energy was found
    if par['single_point_qc'] == 'gaussian':
        if os.path.exists('gaussian/' + job + '.log'):
            gaussname = 'gaussian/' + job + '.log'
        elif os.path.exists('gaussian/' + job + '_high.log'):
            gaussname = 'gaussian/' + job + '_high.log'
        elif os.path.exists('gaussian/' + job + '_well_high.log'):
            gaussname = 'gaussian/' + job + '_well_high.log'
        else:
            logging.info('L3 for {} is missing.'.format(job))
            return 0, -1  # job not yet started to run

        with open(gaussname) as f:
            lines = f.readlines()
            for line in reversed(lines):
                if (key) in line:
                    words = line.split()
                    wi = words.index(key) + 2
                    e = float(words[wi].replace('D', 'E'))
                    logging.info('L3 electronic energy for {} is {} Hartree.'.format(job, e))
                    return 1, e  # energy was found
 
    # if no file or no energy found
    logging.info('L3 for {} is missing.'.format(job))
    return 0, -1  # job not yet started to run or not finished


def get_zpe(jobdir, job, ts, high_level, mp2=0, bls=0):
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
    if zpe == None: 
        logging.warning('Could not find zpe for {} in directory {}'.format(job, jobdir))
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
        os.system('rm -f {dir}/summary_*.out'.format(dir=chemid))
    except OSError:
        pass
    try:
        os.system('rm -f {dir}/kinbot_monitor.out'.format(dir=chemid))
    except OSError:
        pass

    if par['queue_template'] != '':
        shutil.copyfile('{}'.format(par['queue_template']), '{}/{}'.format(chemid, par['queue_template']))
    if par['single_point_template'] != '':
        shutil.copyfile('{}'.format(par['single_point_template']), '{}/{}'.format(chemid, par['single_point_template']))
    if par['barrierless_saddle_single_point_template'] != '':
        shutil.copyfile('{}'.format(par['barrierless_saddle_single_point_template']), '{}/{}'
                        .format(chemid, par['barrierless_saddle_single_point_template']))
        shutil.copyfile('{}'.format(par['barrierless_saddle_prod_single_point_template']), '{}/{}'
                        .format(chemid, par['barrierless_saddle_prod_single_point_template']))
    outfile = open('{dir}/kinbot.out'.format(dir=chemid), 'w')
    errfile = open('{dir}/kinbot.err'.format(dir=chemid), 'w')
    process = subprocess.Popen(command,
                               cwd=chemid,
                               stdout=outfile,
                               stdin=subprocess.PIPE,
                               stderr=errfile)
    time.sleep(1)
    pid = process.pid
    return pid


def write_input(input_file, species, threshold, root):
    # directory for this particular species
    dir = root + '/' + str(species.chemid) + '/'
    if not os.path.exists(dir):
        os.makedirs(dir)

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
    # set the pes option to 1
    par2['pes'] = 1
    # don't do ME for these kinbots but write the files
    par2['me'] = 2

    file_name = dir + str(species.chemid) + '.json'
    with open(file_name, 'w') as outfile:
        json.dump(par2, outfile, indent=4, sort_keys=True)


if __name__ == "__main__":
    main()
