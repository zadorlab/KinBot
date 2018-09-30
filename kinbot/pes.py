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
import sys
import os
import logging
import datetime
import time
import subprocess

from distutils.dir_util import copy_tree

from ase.db import connect

sys.dont_write_bytecode = True

import kinbot
import license_message
from stationary_pt import *
from constants import *

def main(workdir):
    print license_message.message

    sys.path.append(os.path.expanduser('~/KinBot/code'))
    
    
    #os.chdir(os.path.expanduser(workdir))
    
    # set up the logging environment 
    logging.basicConfig(filename=os.path.expanduser(workdir)+'/pes.log',level=logging.INFO)
    
    logging.info(license_message.message)
    logging.info('Starting KinBot in %s at %s'%(workdir,datetime.datetime.now()))
 
    well0 = stationary_pt('well0')
    well0.geom = par.read_input(os.path.expanduser(workdir)+'/input.inp')
    well0.characterize(par.natom, par.atom, par.mult, par.charge)
    write_input(par,well0,workdir)
    
    f = open(os.path.expanduser(workdir)+'/chemids','w')
    f.write(str(well0.chemid) + '\n')
    f.close()
   
    max_running = 10 
    running = []
    finished = []
    jobs = []
    pids = {}
    while 1:
        j = len(jobs)
        f = open(os.path.expanduser(workdir)+'/chemids','r')
        jobs = f.read().split('\n')
        jobs = [ji for ji in jobs if ji != '']
        f.close()
        
        if len(jobs) > j:
            logging.info('Picked up new jobs: ' + ' '.join([ji for ji in jobs[j:]]))

        if len(finished) == len(jobs):
            break
        
        while len(running) < max_running and len(running) + len(finished) < len(jobs):
            #start a new job
            job = jobs[len(running) + len(finished)]
            pid = submit_job(job, workdir)
            pids[job] = pid
            running.append(job)
        #check if a thread is done
        for job in running:
            if not check_status(job,pids[job]):
                finished.append(job)
        #remove the finished threads
        for job in finished: 
            if job in running:
                running.remove(job)
        f = open(os.path.expanduser(workdir)+'/pes_summary.txt','w+')
        f.write('Total\t\t%i\n'%len(jobs))
        f.write('Running\t\t%i\n'%len(running))
        f.write('Finished\t%i\n\n'%len(finished))
        for job in finished: 
            f.write('\t%s\n'%job)
        
        f.close()
        time.sleep(1)
    postprocess(workdir, jobs)

def postprocess(workdir,jobs):
    #list of lists with four elements
    # reactant chemid
    # reaction name
    # products chemid list
    # reaction barrier height
    reactions = [] 
    wells = []
    products = []
    #read all the jobs
    for ji in jobs:
        summary = open(os.path.expanduser(workdir) + '/' + ji + '/summary_' + ji + '.out','r').readlines()
        for line in summary:
            if line.startswith('SUCCESS'):
                pieces = line.split()
                reactant = ji
                ts = pieces[2]
                prod = pieces[3:]
                barrier = float(pieces[1])
                
                if not reactant in wells:
                    wells.append(reactant)
                if len(prod) == 1:
                    if not prod[0] in wells:
                        wells.append(prod[0])
                else:
                    if not '_'.join(sorted(prod)) in products:
                        products.append('_'.join(sorted(prod)))
                new = 1
                temp = None
                for i,rxn in enumerate(reactions):
                    if reactant == rxn[0] and '_'.join(sorted(prod)) == ' '.join(sorted(rxn[2])):
                        new = 0
                        temp = i
                    if reactant == ''.join(rxn[2]) and ''.join(prod) == rxn[0]:
                        new = 0
                        temp = i
                if new:
                    reactions.append([reactant,ts,prod,barrier])
                else:
                    #check if the previous reaction has a lower energy or not
                    if reactions[i][3] > barrier:
                        reactions.pop(temp)
                        reactions.append([reactant,ts,prod,barrier])
    
    zero_energy = get_energy(workdir,jobs[0], jobs[0],0,par.high_level)
    zero_zpe = get_zpe(workdir,jobs[0], jobs[0],0,par.high_level)
    #copy xyz files
    copy_xyz(workdir,wells)
    
    #write pes input
    create_pesviewer_input(workdir, jobs[0], wells, products, reactions, zero_energy,par.high_level)
    
    #write_mess
    create_mess_input(workdir, jobs[0], wells, products, reactions, zero_energy, zero_zpe, par.high_level)

def copy_xyz(workdir,wells):
    dir_xyz = os.path.expanduser(workdir) + '/xyz/'
    if not os.path.exists(dir_xyz):
        os.mkdir(dir_xyz)
    
    for well in wells:
        copy_tree(  os.path.expanduser(workdir) + '/' + well + '/xyz/', dir_xyz)

    
def get_rxn(prods,rxns):
    for rxn in rxns:
        if prods == '_'.join(sorted(rxn[2])):
            return rxn

def create_mess_input(workdir, well0, wells, products, reactions, zero_energy, zero_zpe, high_level):
    fname = os.path.expanduser(workdir) + '/input.mess'
    f = open(fname, 'w+')
    #todo: add header
    
    s = '##############\n'
    s += '# WELLS \n'
    s += '##############\n'
    
    for well in wells:
        energy = get_energy(workdir,well,well,0,par.high_level)
        zpe = get_zpe(workdir,well,well,0,par.high_level)
        zeroenergy = (  ( energy + zpe )- ( zero_energy + zero_zpe) ) * AUtoKCAL
        s += open(os.path.expanduser(workdir) + '/' + well + '/' + well + '.mess').read().format(zeroenergy = zeroenergy) 
        
    for prods in products:
        energy = 0.
        zpe = 0.
        rxn = get_rxn(prods,reactions)
        for pr in prods.split('_'):
            energy += get_energy(workdir,rxn[0],pr,0,par.high_level)
            zpe += get_zpe(workdir,rxn[0],pr,0,par.high_level)
        zeroenergy = (  ( energy + zpe )- ( zero_energy + zero_zpe) ) * AUtoKCAL
        s += open(os.path.expanduser(workdir) + '/' + rxn[0] + '/' + prods + '.mess').read().format(ground_energy = zeroenergy) 
    f.write('\n')
    
    f.write(s)
    f.close()
    

def create_pesviewer_input(workdir, well0, wells, products, reactions, zero_energy, high_level):
    fname = os.path.expanduser(workdir) + '/pesviewer.inp'
    
    f = open(fname,'w+')
    f.write("> <comments>")
    f.write(license_message.message)
    
    f.write("""This comment is not interpreted, so store any extra info here.
Keywords are case insensitive. Look at the help below.
IMPORTANT: avoid the use of '2d' and '3d' in the names of species, transition states and reactions
(these strings are employed when generating the 2d and 3d files of the molecules)
If you want to use 3D coordinates, store them in a xyz/ directory in the same directory as the python script""")
    f.write('\n\n')
    
    f.write('> <id> Potential_energy_surface\n\n')
    f.write("""> <options> 
units              kcal/mol  #energy units
use_xyz            1         # use xyz, put 0  to switch off
rescale            0         # no rescale , put the well or bimolecular name here to rescale to that value
fh                 9.        # figure height
fw                 18.       # figure width
margin             0.2       # margin fraction on the x and y axis
dpi                120       # dpi of the molecule figures
save               0         # does the plot need to be saved (1) or displayed (0)
write_ts_values    1         # booleans tell if the ts energy values should be written
write_well_values  1         # booleans tell if the well and bimolecular energy values should be written
bimol_color        red       # color of the energy values for the bimolecular products
well_color         blue      # color of the energy values of the wells
ts_color           green     # color or the energy values of the ts, put to 'none' to use same color as line
show_images        1         # boolean tells whether the molecule images should be shown on the graph
rdkit4depict       1         # boolean that specifies which code was used for the 2D depiction""")
    f.write('\n\n')
    
    f.write('> <wells> \n')
    for well in wells:
        energy = (get_energy(workdir,well,well,0,par.high_level) - zero_energy) * AUtoKCAL
        f.write('%s %.2f\n'%(well,energy))
    f.write('\n')
    
    f.write('> <bimolec> \n')
    for prods in products:
        energy = 0. - zero_energy
        rxn = get_rxn(prods,reactions)
        for pr in prods.split('_'):
            energy += get_energy(workdir,rxn[0],pr,0,par.high_level)
        energy = energy * AUtoKCAL
        f.write('%s %.2f\n'%(prods,energy))
    f.write('\n')
    
    f.write('> <ts> \n')
    for rxn in reactions:
        energy = (get_energy(workdir,rxn[0],rxn[1],1,par.high_level) - zero_energy) * AUtoKCAL
        prod_name = '_'.join(sorted(rxn[2]))
        f.write('%s %.2f %s %s\n'%(rxn[1],energy,rxn[0],prod_name))
    f.write('\n')
    
    f.write('> <barrierless> \n\n')
    
    f.write("""> <help>
File follows the rules of SD file format for keywords. Keywords are case
insensitive when parsed.
Keywords:
units: units of the energies supplied above

usexyz: use the xyz coordinates of all the species and render a 2D/3D depiction

rescale: energies are rescaled relative to the energy of the species given here 

wells: all the wells of the PES, separated by lines
each line contains the name, the energy, and optionally the smiles

bimolec: all the bimolecular products of the PES, separated by lines
each line contains the name, the energy, and optionally the smiles of both bimolecular products

ts: all the transition states of the PES, separated by lines
each line contains the name, the energy, and the names of the reactant and product

barrierless: all the barrierless reactions of the PES, separated by lines
each line contains the name and the names of the reactant and product""")
    f.close()

def get_energy(workdir,dir,job,ts,high_level):
    db = connect(os.path.expanduser(workdir) + '/' + dir + '/kinbot.db')
    if ts:
        j = job
    else:
        j = job + '_well'
    if high_level:
        j += '_high'
    
    rows = db.select(name = j)
    for row in rows:
        if hasattr(row, 'data'):
            energy = row.data.get('energy')
    #ase energies are always in ev, convert to hartree
    energy *= EVtoHARTREE
    return energy

def get_zpe(workdir,dir,job,ts,high_level):
    db = connect(os.path.expanduser(workdir) + '/' + dir + '/kinbot.db')
    if ts:
        j = job
    else:
        j = job + '_fr'
    if high_level:
        j += '_high'
    
    rows = db.select(name = j)
    for row in rows:
        if hasattr(row, 'data'):
            zpe = row.data.get('zpe')

    return zpe

def check_status(job,pid):
    command = ['ps','-u','root','-N','-o','pid,s,user,%cpu,%mem,etime,args']
    process = subprocess.Popen(command,shell=False,stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    out,err = process.communicate()
    lines = out.split('\n')
    for line in lines:
        if len(line)> 0:
            if '%i'%pid == line.split()[0]:
                return 1
    return 0


def submit_job(chemid, workdir):
    """
    Submit a kinbot run usung subprocess and return the pid
    """
    if '~' in workdir:
        workdir = workdir.replace('~/','')
    job = workdir + '/' + chemid
    command = ["python","kinbot.py",job,"&"]
    process = subprocess.Popen(command,stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    time.sleep(1)
    pid = process.pid
    return pid 

def write_input(par,species,workdir):
    if '~' in workdir:
        dir = os.path.expanduser(workdir) + '/' + str(species.chemid) + '/'
    else:
        dir = workdir + '/' + str(species.chemid) + '/'
    if not os.path.exists(dir):
        os.mkdir(dir)
        
    file_name = dir + 'input.inp'
    f = open(file_name,'w')
    f.write('title = \'%s\'\n\n'%str(species.chemid))
    f.write('method = \'%s\'\n'%par.method)
    f.write('basis = \'%s\'\n'%par.basis)
    f.write('high_level = %i\n'%par.high_level)
    f.write('high_level_method = \'%s\'\n'%par.high_level_method)
    f.write('high_level_basis = \'%s\'\n'%par.high_level_basis)
    
    f.write('qc = \'%s\'\n'%par.qc)
    f.write('queuing = \'%s\'\n'%par.queuing)
    f.write('username = \'%s\'\n'%par.username)
    f.write('queue_name = \'%s\'\n'%par.queue_name)
    f.write('conformer_search = %i\n'%par.conformer_search)
    f.write('reaction_search = %i\n'%par.reaction_search)
    f.write('rotor_scan = %i\n'%par.rotor_scan)
    f.write('barrier_threshold = %.1f\n'%par.barrier_threshold)
    f.write('pes = %i\n'%par.pes)
    f.write('ga = %i\n'%par.ga)
    f.write('ngen = %i\n'%par.ngen)
    f.write('ppn = %i\n\n'%par.ppn)
    
    f.write('charge = %i\n'%par.charge)
    f.write('mult = %i\n'%par.mult)
    f.write('natom = %i\n'%par.natom)
    
    structure = []
    for at in range(par.natom):
        pos = species.geom[at]
        sym = par.atom[at]
        structure += [sym,pos[0],pos[1],pos[2]]

    f.write('structure = %s\n\n'%structure)
    f.close()
    
if __name__ == "__main__":
    workdir = '~/' + sys.argv[1]
    main(workdir)
