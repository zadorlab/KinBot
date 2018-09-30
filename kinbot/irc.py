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
import numpy as np
import os
import copy
import time
import logging
from shutil import copyfile

from vector import *
from qc import *
from constants import *
from kinbot import *
from stationary_pt import *
from geom import *
from qc import *
import par

def irc2stationary_pt(species,index):
    """
    Read the irc files
    There are three possible scenarios: 
    1. One of the ircs leads the initial well and the other to another well or bimolecular product
    2. Neither of the ircs lead to the inital well, transition state structure is not the one kinbot was looking for
    3. Both the ircs lead to the initial well, KinBot found either an identical reaction or the ts is not correct
    """
    instance_name = species.reac_name[index]
    
    directions = ['Forward','Reverse']
    
    ini_well_hits = 0
    prod_hit = -1
    st_pts = [-1,-1]
    for i,direction in enumerate(directions):
        irc_name = '%s_IRC_%s_prod'%(instance_name,direction[0])
        
        err, geom = get_qc_geom(irc_name, par.natom, allow_error = 1)
        if err == -1:
            return 0
        if problem_in_geom(geom):
            #this happens seldomly that all the atoms are very close to one another (problem in Gaussian)
            logging.info('\tProblem with product geometry for %s'%species.reac_name[index])
            return 0

        temp = stationary_pt(irc_name)
        temp.geom = geom

        #temp.characterize(par.natom,par.atom,par.mult,par.charge)
        temp.calc_chemid(par.natom,par.atom,par.mult)
        
        st_pts[i] = temp
        if temp.chemid == species.chemid:
            ini_well_hits += 1
        else:
            prod_hit = i

    if ini_well_hits == 0:
        logging.info('\tNeither IRC leads to the well for %s'%species.reac_name[index])
        return 0
    elif ini_well_hits == 2:
        logging.info('\tBoth IRCs lead to the well, identical reaction found: %s'%species.reac_name[index])
        return 0
    else:
        # ircs OK: well and product found
        logging.info('\tIRCs succesful for %s'%species.reac_name[index])
        return st_pts[prod_hit]

def problem_in_geom(geom):
    #check if interatomic distances are closer than 0.3 Angstrom
    for i in range(len(geom)):
        for j in range(i+1, len(geom)):
            dist = np.linalg.norm(geom[i] - geom[j])
            if dist < 0.3:
                return 1

    return 0

def check_irc(species,index):
    instance_name = species.reac_name[index]
    
    directions = ['Forward','Reverse']
    
    status = [-1,-1]
    
    for i,direction in enumerate(directions):
        irc_name = '%s_IRC_%s'%(instance_name,direction[0])
        status[i] = check_qc(irc_name)
        
    return status

def do_irc_calculations(species, index):
    """
    Carry out the IRC calculation.
    """
    
    instance_name = species.reac_name[index]
    err, geom = get_qc_geom(instance_name, par.natom)

    directions = ['Forward','Reverse']
    
    for i,direction in enumerate(directions):
        irc_name = '%s_IRC_%s'%(instance_name,direction[0])
        
        if par.qc == 'gauss':
            #copy the chk file
            if os.path.exists(instance_name + '.chk'):
                copyfile(instance_name + '.chk',irc_name + '.chk')

        if par.qc == 'nwchem' and direction == 'Reverse':
            direction = 'Backward'
        
        odft = par.mult > 1
        kwargs = get_qc_arguments(irc_name, par.mult, irc = direction.lower())
        prod_kwargs = get_qc_arguments(irc_name + '_prod', par.mult)
        if par.qc == 'gauss':
            prod_kwargs['opt'] = 'CalcFC, Tight, MaxCycle=10'
        
        template = open('{dir}/ase_{qc}_irc.py.tpl'.format(dir = par.tpldir,qc = par.qc),'r').read()
        template = template.format(label = irc_name, kwargs = kwargs, prod_kwargs = prod_kwargs, atom = list(par.atom), 
                               geom = list([list(gi) for gi in geom]), ppn = par.ppn)
        
        f_out = open('{}.py'.format(irc_name),'w')
        f_out.write(template)
        f_out.close()

        submit_qc(irc_name, 0)

    return 0
