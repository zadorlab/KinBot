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
import os, sys
import subprocess
import logging
import numpy as np
import re
import time
import time
import copy

from ase.db import connect

from constants import *
from geom import *
from jobs import *
import par

def get_qc_arguments(job,mult,ts = 0, step = 0, max_step = 0, irc = None,scan = 0,high_level=0, hir = 0):
    """
    Method to get the argument to pass to ase, which are then passed to the qc codes.
    
    Job: name of the job
    
    mult: multiplicity of the job
    
    ts: 1 for transition state searches, 0 for wells and products
    
    step: Step number in the transition state search
    
    max_step: total number of necessary steps to get to the final transition state structure
        This is different for every reaction family
    
    irc: direction of the irc, None if this is not an irc job
    
    scan: is this calculation part of a scan of a bond length to find a maximum energy
    """
    if par.qc == 'gauss':
        # arguments for Gaussian
        kwargs = {
        'method': par.method, 
        'basis': par.basis, 
        'nprocshared' : par.ppn,
        'mem' : '1000MW',
        'chk' : job,
        'label': job, 
        'NoSymm' : 'NoSymm',
        'multiplicity': mult,
        'scf' : 'xqc',
        }
        if ts:
            # arguments for transition state searches
            kwargs['method'] = 'am1'
            kwargs['basis'] = ''
            
            if step == 0:
                kwargs['opt'] = 'ModRedun,Tight,CalcFC,MaxCycle=999'
            elif step < max_step: 
                kwargs['opt'] = 'ModRedun,Tight,CalcFC,MaxCycle=999'
                #kwargs['geom'] = 'AllCheck'
                kwargs['guess'] = 'Read'
            else:
                kwargs['method'] = par.method
                kwargs['basis'] = par.basis
                kwargs['opt'] = 'NoFreeze,TS,CalcFC,NoEigentest,MaxCycle=999'
                kwargs['freq'] = 'freq'
                #kwargs['geom'] = 'AllCheck,NoKeepConstants'
                #kwargs['guess'] = 'Read'
            if scan:
                kwargs['method'] = 'mp2'
                kwargs['basis'] = par.basis
        if irc is not None: 
            #arguments for the irc calculations
            kwargs['geom'] = 'AllCheck,NoKeepConstants'
            kwargs['guess'] = 'Read'
            kwargs['irc'] = 'RCFC,{},MaxPoints=30,StepSize=20'.format(irc)
        if 'R_Addition_MultipleBond' in job:
            #todo: make this generally applicable
            kwargs['method'] = 'mp2'
        if high_level:
            kwargs['method'] = par.high_level_method
            kwargs['basis'] = par.high_level_basis
        if hir:
            kwargs['opt'] = 'ModRedun,CalcFC'
            if ts:
                del kwargs['freq']
                kwargs['opt'] = 'ModRedun,CalcFC,TS,NoEigentest,MaxCycle=999'
        return kwargs
        
    if par.qc == 'nwchem':
        # arguments for NWChem
        odft = mult > 1
        kwargs = {
        'xc':par.method, 
        'basis':par.basis, 
        'scratch_dir' : 'scratch',
        'permanent_dir' : './perm', 
        'label' : job,
        'mult' : mult,
        'odft' : odft,
        'task' : 'optimize',
        'driver' : 'tight',
        'maxiter' : 100,
        'clear' : 1,
        'inhess' : 1
        }
        if ts:
            # arguments for transition state searches
            if step == max_step:
                kwargs['task'] = 'saddle'
            else:
                kwargs['driver'] = 'default'
        if irc is not None: 
            # arguments for the irc calculations
            irc_kwargs = {
            'task' : 'mepgs',
            'mepgs' : 1,
            'maxmep' : 30,
            'maxiter' : 20,
            'inhess' : 2,
            'xyz' : 1,
            'evib' : 0.0005,
            'stride' : 0.1,
            'direction' : irc
            }
            kwargs.update(irc_kwargs)
        return kwargs


def qc_singlept(species, geom, natom, atom, mult, charge, index=-1):
    """ 
    Creates a geometry optimization input and runs it. 
    qc: 'gauss' or 'nwchem'
    wellorts: 0 for wells and 1 for saddle points
    index: >=0 for sampling, each job will get numbered with index
    """
    
    
    if par.qc == 'gauss':
        with open(par.tpldir + par.qc + '_singlept.tpl') as f:
            lines = f.readlines()
    elif par.qc == 'nwchem':
        with open(par.tpldir + par.qc + '.tpl') as f:
            lines = f.readlines()
    
    if index == -1:
        job = str(species.chemid) + '_well'
    elif par.ga == 1:
        job = str(species.chemid) + '_' + str(par.ngen) + '_' + str(index).zfill(par.zf)
    else:
        job = str(species.chemid) + '_' + str(index).zfill(par.zf)

    tpl_replace(job, lines, species, charge, mult, atom, natom, geom, 3)
        
    submit_qc(job)

    return 0

def qc_hir(species, geom, wellorts, natom, atom, mult, charge, rot_index, ang_index,fix):
    """ 
    Creates a constrained geometry optimization input and runs it. 
    wellorts: 0 for wells and 1 for saddle points
    rot_index: index of the rotor in the molecule
    ang_index: index for the current size of the angle
    fix: four atoms of the dihedral that is currently fixed
    """
    if wellorts:
        job = 'hir/' + species.name + '_hir_' + str(rot_index) + '_' + str(ang_index).zfill(2)
    else:
        job = 'hir/' + str(species.chemid) + '_hir_' + str(rot_index) + '_' + str(ang_index).zfill(2)

    kwargs = get_qc_arguments(job,mult, ts = wellorts, step = 1, max_step = 1, high_level = 1, hir = 1)
    kwargs['fix'] = fix
    del kwargs['chk']
    
    dummy = is_linear(geom,species.bond)
    if len(dummy) > 0: # add a dummy atom for each close to linear angle
        for d in dummy:
            atom = np.append(atom,['X'])
            geom = np.concatenate((geom, [d]), axis=0)
            natom += 1
    dummy = [d.tolist() for d in dummy]
    
    template = open(par.tpldir + 'ase_{qc}_hir.py.tpl'.format(qc = par.qc),'r').read()
    template = template.format(label = job, kwargs = kwargs, atom = list(atom), geom = list([list(gi) for gi in geom]), ppn = par.ppn, dummy = dummy)

    f_out = open('{}.py'.format(job),'w')
    f_out.write(template)
    f_out.close()
    
    submit_qc(job)

    return 0


def qc_ring_conf(species, geom, wellorts, natom, atom, mult, charge, relaxed_scan, fix, index=-1):
    """ 
    Creates a constrained geometry optimization input for the conformational search of cyclic structures and runs it.
    qc: 'gauss' or 'nwchem'
    scan: list of dihedrals to be scanned and their values
    wellorts: 0 for wells and 1 for saddle points
    index: >=0 for sampling, each job will get numbered with index
    """
    if index == -1:
        job = 'conf/' + str(species.chemid) + '_well'
    else:
        if wellorts:
            job = 'conf/' + species.name + '_r' + str(index).zfill(par.zf)
        else:
            job = 'conf/' + str(species.chemid) + '_r' + str(index).zfill(par.zf)
    
    kwargs = get_qc_arguments(job,mult, ts = wellorts, step = 1, max_step = 1, hir = 1)
    
    kwargs['relaxed_scan'] = relaxed_scan
    kwargs['fix'] = fix
    kwargs['method'] = 'am1'
    kwargs['basis'] = ''
    del kwargs['chk']
    
    dummy = is_linear(geom,species.bond)
    if len(dummy) > 0: # add a dummy atom for each close to linear angle
        for d in dummy:
            atom = np.append(atom,['X'])
            geom = np.concatenate((geom, [d]), axis=0)
            natom += 1
    dummy = [d.tolist() for d in dummy]
    
    template = open(par.tpldir + 'ase_{qc}_opt_well.py.tpl'.format(qc = par.qc),'r').read()
    template = template.format(label = job, kwargs = kwargs, atom = list(atom), geom = list([list(gi) for gi in geom]), ppn = par.ppn, dummy = dummy)

    f_out = open('{}.py'.format(job),'w')
    f_out.write(template)
    f_out.close()
    
    submit_qc(job)

    return 0

def qc_conf(species, geom, wellorts, natom, atom, mult, charge, index=-1, ring = 0):
    """ 
    Creates a geometry optimization input for the conformational search and runs it.
    qc: 'gauss' or 'nwchem'
    wellorts: 0 for wells and 1 for saddle points
    index: >=0 for sampling, each job will get numbered with index
    """
    if index == -1:
        job = 'conf/' + str(species.chemid) + '_well'
    else:
        r = ''
        if ring: r = 'r'
        if wellorts:
            job = 'conf/' + species.name + '_' + r + str(index).zfill(par.zf)
        else:
            job = 'conf/' + str(species.chemid) + '_' + r + str(index).zfill(par.zf)
    
    if wellorts:
        kwargs = get_qc_arguments(job,mult, ts = 1, step = 1, max_step = 1)
    else:
        kwargs = get_qc_arguments(job,mult)
        if par.qc == 'gauss':
            kwargs['opt'] = 'CalcFC, Tight'
    
    del kwargs['chk']
    
    dummy = is_linear(geom,species.bond)
    if len(dummy) > 0: # add a dummy atom for each close to linear angle
        for d in dummy:
            atom = np.append(atom,['X'])
            geom = np.concatenate((geom, [d]), axis=0)
            natom += 1
    dummy = [d.tolist() for d in dummy]
    
    template = open(par.tpldir + 'ase_{qc}_opt_well.py.tpl'.format(qc = par.qc),'r').read()
    template = template.format(label = job, kwargs = kwargs, atom = list(atom), geom = list([list(gi) for gi in geom]), ppn = par.ppn, dummy = dummy)

    f_out = open('{}.py'.format(job),'w')
    f_out.write(template)
    f_out.close()
    
    submit_qc(job)

    return 0

def qc_opt(species, geom, wellorts, natom, atom, mult, charge, index=-1, high_level = 0, mp2 = 0):
    """ 
    Creates a geometry optimization input and runs it. 
    qc: 'gauss' or 'nwchem'
    wellorts: 0 for wells and 1 for saddle points
    index: >=0 for sampling, each job will get numbered with index
    """
    
    if index == -1:
        job = str(species.chemid) + '_well'
    else:
        job = str(species.chemid) + '_' + str(index).zfill(par.zf)
    if high_level:
        job = str(species.chemid) + '_well_high'
    if mp2:
        job = str(species.chemid) + '_well_mp2'
    
    if wellorts:
        job = species.name + '_' + str(index).zfill(par.zf)
        kwargs = get_qc_arguments(job,mult, ts = 1, step = 1, max_step = 1, high_level = high_level)
    else:
        kwargs = get_qc_arguments(job,mult, high_level = high_level)
        if par.qc == 'gauss':
            kwargs['opt'] = 'CalcFC, Tight'
    
    if mp2:
        kwargs['method'] = 'mp2'
        
    dummy = is_linear(geom,species.bond)
    if len(dummy) > 0: # add a dummy atom for each close to linear angle
        for d in dummy:
            atom = np.append(atom,['X'])
            geom = np.concatenate((geom, [d]), axis=0)
            natom += 1
    dummy = [d.tolist() for d in dummy]
    
    template = open(par.tpldir + 'ase_{qc}_opt_well.py.tpl'.format(qc = par.qc),'r').read()
    template = template.format(label = job, kwargs = kwargs, atom = list(atom), geom = list([list(gi) for gi in geom]), ppn = par.ppn, dummy = dummy)

    f_out = open('{}.py'.format(job),'w')
    f_out.write(template)
    f_out.close()
    
    submit_qc(job)

    return 0
        


def qc_freq(species, geom, natom, atom, mult, charge, index=-1, high_level = 0):
    """ 
    Creates a frequency input and runs it. 
    index: >=0 for sampling, each job will get numbered with index
    """

    if index == -1:           
        job = str(species.chemid) + '_fr'
    else:
        job = str(species.chemid) + '_fr_' + str(index).zfill(par.zf)
    if high_level:
        job = str(species.chemid) + '_fr_high'

    kwargs = get_qc_arguments(job,mult, high_level = high_level)
    if par.qc == 'gauss':
        kwargs['freq'] = 'freq'
        kwargs['ioplist'] = ['7/33=1']
    elif par.qc == 'nwchem':
        kwargs['task'] = 'frequencies'
    
    dummy = is_linear(geom,species.bond)
    if len(dummy) > 0: # add a dummy atom for each close to linear angle
        for d in dummy:
            atom = np.append(atom,['X'])
            geom = np.concatenate((geom, [d]), axis=0)
            natom += 1
        #switch on the symmetry of gaussian
        if 'NoSymm' in kwargs:
            del kwargs['NoSymm']
    dummy = [d.tolist() for d in dummy]
    
    template = open(par.tpldir + 'ase_{qc}_freq_well.py.tpl'.format(qc = par.qc),'r').read()
    template = template.format(label = job, kwargs = kwargs, atom = list(atom), geom = list([list(gi) for gi in geom]), 
                               ppn = par.ppn, dummy = dummy)

    f_out = open('{}.py'.format(job),'w')
    f_out.write(template)
    f_out.close()
    
    submit_qc(job)
    
    return 0



def submit_qc(job, singlejob=1):
    """
    Submit a job to the queue, unless the job:
        * has finished with Normal termination
        * has finished with Error termination
        * is currently running
    However, if the optional parameter singlejob is set to zero, then 
    the job is run only if it has finished earlier with normal termination.
    This is for continuations, when the continuing jobs overwrite each other.
    """

    if singlejob == 1:
        if check_qc(job) != 0: return 0
    else:
        if check_qc(job) == 'running': return 0
        #if check_qc(job) == 'error': return 0


    """#check for existing log or out file and rename those
    ext = ['.out','.log','.com','.nw']
    for e in ext:
        fname = '{}{}'.format(job,e)
        if os.path.exists(fname):
            it = 0
            new_name = '{}_{}{}'.format(job,it,e)
            while os.path.exists(new_name):
                it += 1
                new_name = '{}_{}{}'.format(job,it,e)
            os.rename(fname,new_name)"""
    
    if par.queuing == 'pbs':
        pbs_file = '{}.pbs'.format(job)
        python_file = '{}.py'.format(job)
        python_template = open(par.tpldir + 'pbs_python.tpl','r').read()
        python_template = python_template.format(   name = job, ppn = par.ppn, queue_name = par.queue_name, 
                                                    dir = 'perm', python_file = python_file, arguments = '' )
        f_out_pbs = open(pbs_file,'w')
        f_out_pbs.write(python_template)
        f_out_pbs.close()

        command = ['qsub',job + '.pbs']
        process = subprocess.Popen(command,shell=False,stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err = process.communicate()
        pid = out.split('\n')[0].split('.')[0]
        job_ids[job] = pid
    elif par.queuing == 'slurm':
        slurm_file = '{}.sbatch'.format(job)
        python_file = '{}.py'.format(job)
        python_template = open(par.tpldir + 'slurm_python.tpl','r').read()
        python_template = python_template.format(   name = job, ppn = par.ppn, queue_name = par.queue_name, dir = 'perm', 
                                                    slurm_feature = slurm_feature, python_file = python_file, arguments = '' )
        f_out_slurm = open(slurm_file,'w')
        f_out_slurm.write(python_template)
        f_out_slurm.close()

        command = ['sbatch',job + '.sbatch']
        process = subprocess.Popen(command,shell=False,stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err = process.communicate()
        pid = out.split('\n')[0].split()[-1]
        job_ids[job] = pid
    else:
        logging.error('KinBot does not recognize queuing system {}.'.format(par.queuing))
        logging.error('Exiting')
        sys.exit()
    
    return 1 # important to keep it 1, this is the natural counter of jobs submitted


def get_qc_geom(job, natom, wait=0, allow_error = 0):
    """
    Get the geometry from the ase database file.
    Returns it, with the following conditions about the status of the job.
    If wait = 0, return an (1, empty array) when the job is still running (instead of the geometry).
    If wait = 1, wait for the job to finish.
    If wait = 2, return the last geometry while the job is still running.
        This option is to monitor the formation of bimolecular products.
    If allow_error = 0, do not read the final geometry if the status is not "normal"
    if allow_error = 1, read the geometry even though there is an error in the output file
        This option is to read the final IRC geometry when it did not converge
    """ 
    
   
    geom = np.zeros((natom,3))    
    
    if check_qc(job) == 'error' and not allow_error: return -1, geom
    status = 0
    while 1:
        if check_qc(job) == 'running': 
            if wait == 1:
                time.sleep(1)
            elif wait == 2:
                status = 2
                break
            else: 
                return 1, geom 
        else:
            break
    if check_qc(job) != 'normal':
        if not allow_error:
            if wait != 2: return -1, geom

    #open the database
    rows = par.db.select(name = job)
    
    found_entry = 0
    #take the last entry
    for row in rows:
        mol = row.toatoms()
        geom = mol.positions
        found_entry = 1
    
    if found_entry:
        return status, geom
    else:
        return -1, geom

def get_second_to_last_geom(job, natom, wait=0, allow_error = 0):
    """
    Get the geometry from the ase database file.
    Returns it, with the following conditions about the status of the job.
    If wait = 0, return an (1, empty array) when the job is still running (instead of the geometry).
    If wait = 1, wait for the job to finish.
    If wait = 2, return the last geometry while the job is still running.
        This option is to monitor the formation of bimolecular products.
    If allow_error = 0, do not read the final geometry if the status is not "normal"
    if allow_error = 1, read the geometry even though there is an error in the output file
        This option is to read the final IRC geometry when it did not converge
    """ 
    
   
    geom = np.zeros((natom,3))    
    
    if check_qc(job) == 'error' and not allow_error: return -1, geom
    status = 0
    while 1:
        if check_qc(job) == 'running': 
            if wait == 1:
                time.sleep(1)
            elif wait == 2:
                status = 2
                break
            else: 
                return 1, geom 
        else:
            break
    if check_qc(job) != 'normal':
        if not allow_error:
            if wait != 2: return -1, geom

    #open the database
    rows = par.db.select(name = job)
    
    found_entry = 0
    geoms = []
    #take the last entry
    for row in rows:
        mol = row.toatoms()
        geoms.append(mol.positions)
        found_entry = 1
    
    if found_entry:
        return status, geoms[-2]
    else:
        return -1, geom

def get_qc_freq(job, natom, wait=0, allow_error = 0):
    """
    Get the frequencies from the ase database file
    If wait is set to 1, it will wait for the job to finish.
    """ 
    
    if check_qc(job) == 'error': return -1, [0]
    while 1:
        if check_qc(job) == 'running':
            if wait == 1:
                time.sleep(1)
            else:
                return 1, []
        else:
            break
    
    if check_qc(job) != 'normal': return -1, [0]

    freq = []
    
    #open the database
    rows = par.db.select(name = job)
    
    #take the last entry
    for row in rows:
        if hasattr(row, 'data'):
            if not row.data.get('frequencies') is None:
                freq = list(row.data.get('frequencies'))
    
    if len(freq) == 0 and natom > 1:
        return -1,freq

    return 0, freq

def get_qc_energy(job, wait=0):
    """
    Read the last energy from a job. 
    For Gaussian currently works for DFT and HF only.
    For NWChem it works for optimization jobs, using the @ handle.
    If wait is set to 1, it will wait for the job to finish, otherwise
    just reads the last energy in the file.
    Returns the error code and the energy
    Error codes:
    -1: error
     0: normal and done
     1: running
    """
    


    if check_qc(job) == 'error': return -1, 0.
    while 1:
        if check_qc(job) == 'running': 
            if wait == 1:
                time.sleep(1)
            else:
                return 1, 0.
        else:
            break

    #open the database
    rows = par.db.select(name = job)
    
    #take the last entry
    for row in rows:
        if hasattr(row, 'data'):
            if not row.data.get('energy') is None:
                energy = row.data.get('energy')
    
    #ase energies are always in ev, convert to hartree
    energy *= EVtoHARTREE
    
    return 0, energy



def get_qc_zpe(job, wait=1):
    """
    Read the zero point energy. 
    If wait is set to 1 (default), it will wait for the job to finish.
    """
    if check_qc(job) == 'error': return -1, 0.
    while 1:
        if check_qc(job) == 'running': 
            if wait == 1:
                time.sleep(1)
            else:
                return 0, 0.
        else:
            break

    #open the database
    rows = par.db.select(name = job)

    #take the last entry
    for row in rows:
        if hasattr(row, 'data'):
            zpe = row.data.get('zpe')
    
    return 0, zpe

def read_qc_hess(job, natom):
    """
    Read the hessian of a gaussian chk file
    """
    if check_qc(job) != 'normal': 
        return []
    
    #initialize hessian matrix
    hess = np.zeros((3*natom,3*natom))
    
    if par.qc == 'gauss':
        
        fchk = str(job) + '.fchk'
        #if not os.path.exists(fchk):
        #create the fchk file using formchk
        os.system('formchk ' + job + '.chk > /dev/null')
        
        with open(fchk) as f:
            lines = f.read().split('\n')
        
        nvals = 3*natom * (3 * natom +1) / 2

        for index, line in enumerate(reversed(lines)):
            if re.search('Cartesian Force Constants', line) != None:
                hess_flat = []
                n = 0
                while len(hess_flat) < nvals:
                    hess_flat.extend([float(val) for val in lines[-index + n].split()])
                    n += 1
                n = 0
                for i in range(3*natom):
                    for j in range(i+1):
                        hess[i][j] = hess_flat[n]
                        hess[j][i] = hess_flat[n]
                        n += 1
                break
    return hess

    
def read_qc_geom(job, natom, wait=0, allow_error = 0):
    """
    Reads the last geometry from an output file.
    Returns it, with the following conditions about the status of the job.
    If wait = 0, return an (1, empty array) when the job is still running (instead of the geometry).
    If wait = 1, wait for the job to finish.
    If wait = 2, return the last geometry while the job is still running.
        This option is to monitor the formation of bimolecular products.
    If allow_error = 0, do not read the final geometry if the status is not "normal"
    if allow_error = 1, read the geometry even though there is an error in the output file
        This option is to read the final IRC geometry when it did not converge
    """ 
    
    
    geom = np.zeros((natom,3))    
    
    if check_qc(job) == 'error' and not allow_error: return -1, geom
    status = 0
    while 1:
        if check_qc(job) == 'running': 
            if wait == 1:
                time.sleep(1)
            elif wait == 2:
                status = 2
                break
            else: 
                return 1, geom 
        else:
            break

    if check_qc(job) != 'normal':
        if not allow_error:
            if wait != 2: return -1, geom

    if par.qc == 'gauss': 
        outfile = str(job) + '.log'
    elif par.qc == 'nwchem':
        outfile = str(job) + '.out'
    
    if not os.path.isfile(outfile) and wait == 2:
        return 1, geom
    
    with open(outfile) as f:
        lines = f.readlines()

    if par.qc == 'gauss':        
        for index, line in enumerate(reversed(lines)):
            if re.search('Input orientation:', line) != None:
                for n in range(natom):
                    geom[n][0:3] = np.array(lines[-index+4+n].split()[3:6]).astype(float)
                break
    elif par.qc == 'nwchem':
        for index, line in enumerate(reversed(lines)):
            if re.search('Output coordinates in angstroms', line) != None:
                for n in range(natom):
                    geom[n][0:3] = np.array(lines[-index+3+n].split()[3:6]).astype(float)
                break        
                            
    return status, geom
 
        

def read_qc_grad(job, natom, wait=1):
    """
    Reads the last geometry from an output file.
    Returns it, with the following conditions about the status of the job.
    If wait = 0, return an empty array when the job is still running (instead of the gradients).
    If wait = 1, wait for the job to finish.
    If wait = 2, return the last gradient while the job is still running.
    """ 
    

    grad = np.zeros((natom,3))    
    
    if check_qc(job) == 'error': return -1, geom
    status = 0
    while 1:
        if check_qc(job) == 'running': 
            if wait == 1:
                time.sleep(1)
            elif wait == 2:
                status = 2
                break
            else: 
                status = 1
                return status, geom 
        else:
            break

    if check_qc(job) != 'normal' and wait != 2: return -1, geom

    
    if par.qc == 'gauss': 
        outfile = str(job) + '.log'
    elif par.qc == 'nwchem':
        outfile = str(job) + '.out'
        
    with open(outfile) as f:
        lines = f.readlines()

    if par.qc == 'gauss':  # FIXME, the Gaussian version is not yet done      
        for index, line in enumerate(reversed(lines)):
            if re.search('Input orientation:', line) != None:
                for n in range(natom):
                    geom[n][0:3] = np.array(lines[-index+4+n].split()[3:6]).astype(float)
                break
    elif par.qc == 'nwchem':
        for index, line in enumerate(reversed(lines)):
            if re.search('GRADIENTS', line) != None:
                for n in range(natom):
                    grad[n][0:3] = np.array(lines[-index+3+n].split()[5:8]).astype(float)
                break        
    return status, grad


def read_qc_freq(job, natom, wait=0):
    """
    Reads the frequencies.
    If wait is set to 1, it will wait for the job to finish.
    """ 
    
    if check_qc(job) == 'error': return -1, [0]
    while 1:
        if check_qc(job) == 'running':
            if wait == 1:
                time.sleep(1)
            else:
                return 1, []
        else:
            break

    if check_qc(job) != 'normal': return -1, [0]

    if natom == 1:
        return 0, np.array([])
    elif natom == 2:
        freq = np.array([0.])
    else:
        freq = np.zeros(3 * natom - 6) 

    freq_all = np.zeros(3 * natom)    
            
    if par.qc == 'gauss': 
        outfile = str(job) + '.log'
    elif par.qc == 'nwchem':
        outfile = str(job) + '.out'
          
    with open(outfile) as f:
        lines = f.readlines()
    
    nfreq = 0

    if par.qc == 'gauss':
        for line in lines:
            if re.search('Frequencies', line) != None:
                if natom == 2:
                    freq[0] = np.array(line.split()[2]).astype(float)
                    break
                else:
                    freq[nfreq:nfreq+3] = np.array(line.split()[2:5]).astype(float)
                    nfreq = nfreq + 3
                    
    elif par.qc == 'nwchem': 
        for line in lines:
            if re.search('P.Frequency', line) != None:
                nval = len(line.split()) - 1
                freq_all[nfreq:nval+nfreq] = np.array(line.split()[1:nval+1]).astype(float)
                nfreq += nval
        i = 0
        for fr in range(nfreq):
            if abs(freq_all[fr]) > 5:
                freq[i] = freq_all[fr]
                i += 1       
                         
    return 0, freq



def read_qc_imode(job, natom):
    """
    Read the normal mode for the imaginary frequency.
    """ 
    
    
    imode = np.zeros((natom, 3))

    if check_qc(job) == 'error': return -1, imode
    if check_qc(job) == 'running': return 1, imode
    if natom == 1:
        return 0, imode
         
    if par.qc == 'gauss': 
        outfile = str(job) + '.log'
    elif par.qc == 'nwchem':
        outfile = str(job) + '.out'
           
    with open(outfile) as f:
        lines = f.readlines()
    
    if par.qc == 'gauss':
        for index, line in enumerate(lines):
            if re.search('Frequencies', line) != None:
                for atom in range(natom):
                    imode[atom][0:3] = np.array(lines[atom+index+5].split()[2:5]).astype(float)
                break
    elif par.qc == 'nwchem': 
        for index, line in enumerate(lines):
            if re.search('P.Frequency', line) != None:
                for atom in range(natom): 
                    imode[atom][0] = np.array(lines[3*atom+index+2].split()[1]).astype(float)
                    imode[atom][1] = np.array(lines[3*atom+index+3].split()[1]).astype(float)
                    imode[atom][2] = np.array(lines[3*atom+index+4].split()[1]).astype(float)
                break                       
                                                                     
    return 0, imode



def read_qc_nmode(job, natom):
    """
    Read all normal modes.
    """ 
    
    
    if natom == 1:
        return 0, np.array([])
    elif natom == 2:
        nm = 1
    else:
        nm = 3 * natom - 6
        
    nmode = np.zeros((nm, natom, 3)) 
    nmode_all = np.zeros((3 * natom, natom, 3)) 
    freq_all = np.zeros(3 * natom) 
        
    if check_qc(job) == 'error': return -1, nmode
    if check_qc(job) == 'running': return 1, nmode
         
    if par.qc == 'gauss': 
        outfile = str(job) + '.log'
    elif par.qc == 'nwchem':
        outfile = str(job) + '.out'
           
    with open(outfile) as f:
        lines = f.readlines()
    
    if par.qc == 'gauss':
        mode = 0
        for index, line in enumerate(lines):
            if re.search('Frequencies', line) != None and nm > 1:
                for atom in range(natom):
                    nmode[mode,atom,0:3] = np.array(lines[atom+index+5].split()[2:5]).astype(float)
                    nmode[mode+1,atom,0:3] = np.array(lines[atom+index+5].split()[5:8]).astype(float)                    
                    nmode[mode+2,atom,0:3] = np.array(lines[atom+index+5].split()[8:11]).astype(float)
                mode += 3
            elif re.search('Frequencies', line) != None:
                for atom in range(natom):
                    nmode[0,atom,0:3] = np.array(lines[atom+index+5].split()[2:5]).astype(float)
                                        
    elif par.qc == 'nwchem':
        mode = 0 
        for index, line in enumerate(lines):
            if re.search('P.Frequency', line) != None:
                nval = len(line.split())-1
                freq_all[mode:mode+nval] = np.array(line.split()[1:nval+1]).astype(float)
                for atom in range(natom):
                    for coo in range(3):
                        nmode_all[mode:mode+nval,atom,coo] = np.array(lines[3*atom+index+2+coo].split()[1:2+nval]).astype(float)
                mode += nval

        i = 0
        for m in range(mode):
            if abs(freq_all[m]) > 5:
                nmode[i] = nmode_all[m]
                i += 1       
                
                                                                     
    return 0, nmode



def read_qc_energy(job, wait=0):
    """
    Read the last energy from a job. 
    For Gaussian currently works for DFT and HF only.
    For NWChem it works for optimization jobs, using the @ handle.
    If wait is set to 1, it will wait for the job to finish, otherwise
    just reads the last energy in the file.
    Returns the error code and the energy
    Error codes:
    -1: error
     0: normal and done
     1: running
    """


    if check_qc(job) == 'error': return -1, 0.
    while 1:
        if check_qc(job) == 'running': 
            if wait == 1:
                time.sleep(1)
            else:
                return 1, 0.
        else:
            break

    if par.qc == 'gauss': 
        outfile = str(job) + '.log'
    elif par.qc == 'nwchem':
        outfile = str(job) + '.out'
        
    with open(outfile) as f:
        lines = f.readlines()
    
    if par.qc == 'gauss':
        for line in reversed(lines):
            if re.search('SCF Done', line) != None:
                energy = float(line.split()[4])
                break
    elif par.qc == 'nwchem':
        for line in reversed(lines):  
            if re.search('@', line) != None:
                energy = float(line.split()[2])
                break          
            if re.search('Total', line) != None and re.search('energy = ', line) != None:
                energy = float(line.split()[4])
                break
                       
    return 0, energy



def read_qc_zpe(job, wait=1):
    """
    Read the zero point energy. 
    If wait is set to 1 (default), it will wait for the job to finish.
    """


    if check_qc(job) == 'error': return -1, 0.
    while 1:
        if check_qc(job) == 'running': 
            if wait == 1:
                time.sleep(1)
            else:
                return 0, 0.
        else:
            break

    if par.qc == 'gauss': 
        outfile = str(job) + '.log'
    elif par.qc == 'nwchem':
        outfile = str(job) + '.out'
        
    with open(outfile) as f:
        lines = f.readlines()
    
    if par.qc == 'gauss':
        for line in reversed(lines):
            if re.search('Zero-point correction=', line) != None:
                zpe = float(line.split()[2])
                break
    elif par.qc == 'nwchem':
        for line in reversed(lines):  
            if re.search('Zero-Point correction to Energy', line) != None:
                zpe = float(line.split()[8])
                break          
                      
    return 0, zpe
    
def is_in_database(job):
    """
    Checks if the current job is in the database:
    """
    #open the database
    rows = par.db.select(name = job)
    
    mol = None
    
    #take the last entry
    for row in rows:
        mol = row.toatoms()
    
    if mol is None:
        return 0
    
    return 1

def check_qc(job):
    """
    Checks the status of a Gaussian or NWChem job.
    """
    if par.qc == 'gauss':
        log_file = job + '.log'
    elif par.qc == 'nwchem':
        log_file = job + '.out'
    log_file_exists = os.path.exists(log_file)
    
    devnull = open(os.devnull, 'w')
    if par.queuing == 'pbs':
        command = 'qstat -f | grep ' + '"Job Id: ' + job_ids.get(job,'-1') + '"' + ' > /dev/null'
    elif par.queuing == 'slurm':
        command = 'scontrol show job ' + job_ids.get(job,'-1') + ' | grep "JobId=' + job_ids.get(job,'-1') + '"' + ' > /dev/null'
    else:
        logging.error('KinBot does not recognize queuing system {}.'.format(par.queuing))
        logging.error('Exiting')
        sys.exit()
    if int(subprocess.call(command, shell = True, stdout=devnull, stderr=devnull)) == 0: 
        return 'running' 
    elif is_in_database(job) and log_file_exists: #by deleting a log file, you allow restarting a job
        #open the database
        rows = par.db.select(name = job)
        data = None
        #take the last entry
        for row in rows:
            if hasattr(row,'data'):
                data = row.data
        if data is None:
            return 0
        else:
            return data['status']
    else: 
        return 0
        


def kill_job(job):
    """ 
    Deletes job from the queue.
    First checks if the job is running.
    """
    

    job = str(job)
    if check_qc(job) != 'running':
        return -1
             
    command = 'qstat -f > running.jobs'
    os.system(command)
    
    with open('running.jobs') as f:
        lines = f.readlines()
        
    jobid = 0
    for line in lines:
        if re.search('Job Id:', line) != None:
            jobid = line.split()[2].split(".")[0]
        if re.search(job, line) != None:
            break
            
    if jobid != 0:
        command = 'qdel ' + jobid
        os.system(command)
    
    return 0



def tpl_replace(fname, lines, species, charge, mult, atom, natom, geom, task):
    """ 
    Replace strings in template files.
    task is for NWChem files
        0: optimize
        1: saddle
        2: freqency (and Hessian)
        3: gradient
    """ 
    
    
    if par.qc == 'gauss':
        f = open(fname + '.com', 'w')    
    elif par.qc == 'nwchem':
        f = open(fname + '.inp', 'w')
        
    for line in lines:
        line = line.replace('_chk_', str(species.chemid))
        line = line.replace('_charge_', str(charge))
        if mult > 1:
            line = line.replace('_odft_', 'odft')
        else:
            line = line.replace('_odft_', '')  

        line = line.replace('_scfmultiplicity_', scf_mult(mult))
        line = line.replace('_multiplicity_', str(mult))
        line = line.replace('_method_', str(par.method))
        line = line.replace('_basis_', str(par.basis))
        line = line.replace('_methodclass_', str(par.methodclass))
        line = line.replace('_name_', str(fname))
        line = line.replace('_ppn_', str(par.ppn))
        line = line.replace('_scratch_', str(par.scratch))
        line = line.replace('_level_', str(par.level))
        line = line.replace('_convergence_', str(par.convergence))
        line = line.replace('_maxiter_', str(par.maxiter))
        if task == 0:
            line = line.replace('_task_', 'optimize')
        elif task == 1:
            line = line.replace('_task_', 'saddle')
        elif task == 2:
            line = line.replace('_task_', 'frequencies')
        elif task == 3:
            line = line.replace('_task_', 'gradient')
        if line.split():
            if line.split()[0] == '_geom_':
                for atomi in range(natom):
                    f.write(atom[atomi] + '  ')
                    #[f.write(str(geom[atomi][i]) + '  ') for i in range(3)]
                    [f.write("{g:12.8f}".format(g=geom[atomi][i]) ) for i in range(3)]
                    f.write('\n')
            else:
                f.write(line)
        else:
            f.write(line)

    f.close()



def scf_mult(mult):
    """
    Creates string for NWChem multiplicity in the SCF part of the deck.
    """
    
    
    if mult == 1:
        return 'singlet'
    elif mult == 2:
        return 'doublet'
    elif mult == 3:
        return 'triplet'



def main():
    """
    This module contains reading and submission of Gaussian and NWChem jobs. 
    """



if __name__ == "__main__":
    main()
