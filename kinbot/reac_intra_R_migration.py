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


def do_intra_R_migration(species, instance, step, instance_name):
    """
    Carry out the reaction.
    """
    
    if step > 0:
        if check_qc(instance_name) != 'normal' and check_qc(instance_name) != 'error': return step
    
    #maximum number of steps for this reaction family
    max_step = 14

    # have a look at what has been done and what needs to be done
    skip = 1
    step,geom, kwargs = initialize_reaction(species, instance, step, instance_name, max_step, skip)

    #the the constraints for this step
    step, fix, change, release = get_constraints_intra_R_migration(step, species, instance,geom)
    
    #carry out the reaction and return the new step
    return carry_out_reaction(species, instance, step, instance_name, max_step, geom, kwargs, fix, change, release)

def get_constraints_intra_R_migration(step,species,instance,geom):
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
    elif step == 12:
        if len(instance) > 3:
            for dih in range(len(instance)-3):
                f = []
                for i in range(4):
                    f.append(instance[dih+i] + 1)
                fix.append(f)
            for angle in range(len(instance)-2):
                constraint = []
                for i in range(3):
                    constraint.append(instance[angle+i] + 1)
                constraint.append(180. * (len(instance)-2.) / len(instance))
                change.append(constraint)
        else:
            #val = 60.0
            #constraint = [instance[-2] + 1,instance[-1] + 1,instance[0] + 1,val]
            #change.append(constraint)
            
            fval = 1.35
            if par.atom[instance[0]] == 'O': fval = 1.2
            constraint = [instance[0] + 1,instance[-1] + 1,fval]
            change.append(constraint)
            
            fval = 1.35
            if par.atom[instance[-2]] == 'O': fval = 1.2
            constraint = [instance[-2] + 1,instance[-1] + 1,fval]
            change.append(constraint)
            step = 13
    elif step == 13:
        for angle in range(len(instance)-2):
            constraint = []
            for i in range(3):
                constraint.append(instance[angle+i] + 1)
            release.append(constraint)
        for dih in range(len(instance)-3):  
            constraint = []
            for i in range(4):
                constraint.append(instance[dih+i] + 1)
            release.append(constraint)
            
        fval = st_bond[''.join(sorted(par.atom[instance[0]]+par.atom[instance[-1]]))]*1.0
        constraint = [instance[0] + 1,instance[-1] + 1,fval]
        change.append(constraint)
        
        #if len(instance) > 3 and np.sum(species.rad) == 0:
        #    fval = 1.35
        #    constraint = [instance[1] + 1,instance[-1] + 1,fval]
        #    change.append(constraint)
        
        #fval = 1.0
        #constraint = [instance[1] + 1,instance[-1] + 1,fval]
        #change.append(constraint)
        
        fval = st_bond[''.join(sorted(par.atom[instance[-2]]+par.atom[instance[-1]]))]*1.0
        constraint = [instance[-2] + 1,instance[-1] + 1,fval]
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
