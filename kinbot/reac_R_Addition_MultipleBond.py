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
import copy
import time

from vector import *
from qc import *
from constants import *
from kinbot import *
from stationary_pt import *
from geom import *
from qc import *
from reac_family import *
import par



def do_R_Addition_MultipleBond(species, instance, step, instance_name):
    """
    Carry out the reaction in reverse, which is a beta scission.
    The number of steps of 0.1 Angstrom is controlled by a parameter.
    """
    
    if step > 0:
        if check_qc(instance_name) != 'normal' and check_qc(instance_name) != 'error': return step
    
    #maximum number of steps for this reaction family
    max_step = 12

    # have a look at what has been done and what needs to be done
    step,geom, kwargs = initialize_reaction(species, instance, step, instance_name, max_step,scan = 1)

    #the the constraints for this step
    step, fix, change, release = get_constraints_R_Addition_MultipleBond(step, species, instance,geom,max_step)
    
    #carry out the reaction and return the new step
    return carry_out_reaction(species, instance, step, instance_name, max_step, geom, kwargs, fix, change, release)

def get_constraints_R_Addition_MultipleBond(step,species,instance,geom,max_step):
    """
    There are three types of constraints:
    1. fix a coordinate to the current value
    2. change a coordinate and fix is to the new value
    3. release a coordinate (only for gaussian)
    """
    fix = []
    change = []
    release = []
    #fix all the bond lengths
    #~ if step < max_step:
        #~ for i in range(par.natom - 1):
            #~ for j in range(i+1, par.natom):
                #~ if species.bond[i][j] > 0:
                    #~ fix.append([i+1,j+1])
    if step == 0:
        # verify if the radical atom has more than one neighbor and 
        # change the dihedral to 90 degrees in that case
        neigh = [i for i,ni in enumerate(species.bond[instance[0]]) if ni > 0]
        if len(neigh) > 1:
            for ni in neigh:
                if not ni in instance:
                    change.append([ni + 1, instance[0] + 1, instance[1] + 1, instance[2] + 1, 90.])
                    break
            
    if step < 12:
        final_dist = 1.42
        if par.atom[instance[0]] == 'C' and par.atom[instance[1]] == 'C' and par.atom[instance[2]] == 'C':
            final_dist = 2.20
        if par.atom[instance[0]] == 'C' and par.atom[instance[1]] == 'C' and par.atom[instance[2]] == 'H':
            final_dist = 1.79
        if par.atom[instance[0]] == 'C' and par.atom[instance[1]] == 'C' and par.atom[instance[2]] == 'O':
            final_dist = 2.04
        if par.atom[instance[0]] == 'O' and par.atom[instance[1]] == 'C' and par.atom[instance[2]] == 'C':
            final_dist = 2.12
        if par.atom[instance[0]] == 'O' and par.atom[instance[1]] == 'C' and par.atom[instance[2]] == 'H':
            final_dist = 1.84
        if par.atom[instance[0]] == 'O' and par.atom[instance[1]] == 'C' and par.atom[instance[2]] == 'O':
            final_dist = 2.04 #TODO: verify if this value is OK
        if par.atom[instance[0]] == 'C' and par.atom[instance[1]] == 'O' and par.atom[instance[2]] == 'C':
            final_dist = 2.04 #TODO: verify if this value is OK
        if par.atom[instance[0]] == 'C' and par.atom[instance[1]] == 'O' and par.atom[instance[2]] == 'H':
            final_dist = 1.42
        if par.atom[instance[0]] == 'C' and par.atom[instance[1]] == 'O' and par.atom[instance[2]] == 'O':
            final_dist = 2.04 #TODO: verify if this value is OK
        if par.atom[instance[0]] == 'O' and par.atom[instance[1]] == 'O' and par.atom[instance[2]] == 'C':
            final_dist = 2.04 #TODO: verify if this value is OK
        
        val = new_bond_length(species,instance[1],instance[2],step,12,final_dist,geom)
        constraint = [instance[1] + 1,instance[2] + 1,val]
        change.append(constraint)

    
    #remove the bonds from the fix if they are in another constaint
    for c in change:
        if len(c) == 3:
            index = -1
            for i,fi in enumerate(fix):
                if len(fi) == 2:
                    if sorted(fi) == sorted(c[:2]):
                        index = i
            if index > -1:
                del fix[index]
    
    return step, fix, change, release


