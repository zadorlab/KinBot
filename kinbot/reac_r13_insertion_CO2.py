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



def do_r13_insertion_CO2(species, instance, step, instance_name):
    """
    Carry out the reaction.
    """
    
    if step > 0:
        if check_qc(instance_name) != 'normal' and check_qc(instance_name) != 'error': return step
    
    #maximum number of steps for this reaction family
    max_step = 22

    # have a look at what has been done and what needs to be done
    step,geom, kwargs = initialize_reaction(species, instance, step, instance_name, max_step)

    #the the constraints for this step
    step, fix, change, release = get_constraints_r13_insertion_CO2(step, species, instance,geom)
    
    #carry out the reaction and return the new step
    return carry_out_reaction(species, instance, step, instance_name, max_step, geom, kwargs, fix, change, release)

def get_constraints_r13_insertion_CO2(step,species,instance,geom):
    """
    There are three types of constraints:
    1. fix a coordinate to the current value
    2. change a coordinate and fix is to the new value
    3. release a coordinate (only for gaussian)
    """
    fix = []
    change = []
    release = []
    #if step < 12:
    #    #fix all the bond lengths
    #    for i in range(par.natom - 1):
    #        for j in range(i+1, par.natom):
    #            if species.bond[i][j] > 0:
    #                fix.append([i+1,j+1])
    if step < 12:
        new_dihs = new_ring_dihedrals(species, instance, step, 12)
        for dih in range(len(instance)-3):
            constraint = []
            for i in range(4):
                constraint.append(instance[dih+i] + 1)
            constraint.append(new_dihs[dih])
            change.append(constraint)
    elif step < 22:
        for dih in range(len(instance)-3):  
            constraint = []
            for i in range(4):
                constraint.append(instance[dih+i] + 1)
            release.append(constraint)
            
        fval = [2.0,1.3,2.0,2.0]
        
        if par.atom[instance[-1]] == 'H':
            fval[2] = 1.35
            fval[3] = 1.35
        
        val = new_bond_length(species,instance[0],instance[1],step - 11,10,fval[0],geom)
        constraint = [instance[0] + 1,instance[1] + 1,val]
        change.append(constraint)
        
        val = new_bond_length(species,instance[1],instance[2],step - 11,10,fval[1],geom)
        constraint = [instance[1] + 1,instance[2] + 1,val]
        change.append(constraint)
        
        val = new_bond_length(species,instance[2],instance[3],step - 11,10,fval[2],geom)
        constraint = [instance[2] + 1,instance[3] + 1,val]
        change.append(constraint)
        
        val = new_bond_length(species,instance[3],instance[0],step - 11,10,fval[3],geom)
        constraint = [instance[3] + 1,instance[0] + 1,val]
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

