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
Generic methods for the reaction families
"""

import numpy as np
import copy
import time

from vector import *
from qc import *
from constants import *
from kinbot import *
from stationary_pt import *
from geom import *
from qc import *
from modify_geom import *
import par

def initialize_reaction(species, instance, step, instance_name,max_step,skip = 0,scan = 0):
    """
    Verify what has been done and what needs to be done
    
    skip: boolean which tells to skip the first 12 steps in case of an instance shorter than 4
    
    scan: boolean which tells if this is part of an energy scan along a bond length coordinate
    """
    kwargs = get_qc_arguments(instance_name,par.mult,ts = 1,step = step,max_step=max_step,scan = scan)
    
    if step == 0:
        if is_in_database(instance_name):
            if check_qc(instance_name) == 'normal': 
                err, freq = get_qc_freq(instance_name, par.natom)
                if err == 0 and len(freq) > 0.:
                    err, geom = get_qc_geom(instance_name, par.natom)
                    step = max_step + 1
                    return step, geom, kwargs
        if skip and len(instance) < 4: 
            step = 12
        geom = species.geom
    else:
        err, geom = get_qc_geom(instance_name, par.natom, allow_error = 1)
        
    
    return step, geom, kwargs

def carry_out_reaction(species, instance, step, instance_name, max_step, geom, kwargs, fix, change, release):
    if step > max_step:
        return step
    
    
    #apply the geometry changes here and fix the coordinates that changed
    change_starting_zero = []
    for c in change:
        c_new = [ci - 1 for ci in c[:-1]]
        c_new.append(c[-1])
        change_starting_zero.append(c_new)
    if len(change_starting_zero) >0 :
        geom = modify_coordinates(species,instance_name,geom,change_starting_zero,species.bond,par.natom,par.atom)
        for c in change:
            fix.append(c[:-1])
        change = []
    
    
    kwargs['fix'] = fix
    kwargs['change'] = change
    kwargs['release'] = release

    if step < max_step:
        template = open(par.tpldir + 'ase_{qc}_ts_search.py.tpl'.format(qc = par.qc),'r').read()
    else:
        template = open(par.tpldir + 'ase_{qc}_ts_end.py.tpl'.format(qc = par.qc),'r').read()
    
    template = template.format(label = instance_name, kwargs = kwargs, atom = list(par.atom), 
                               geom = list([list(gi) for gi in geom]), ppn = par.ppn)
    
    f_out = open('{}.py'.format(instance_name),'w')
    f_out.write(template)
    f_out.close()
    
    step += submit_qc(instance_name, 0)
    
    return step
    
