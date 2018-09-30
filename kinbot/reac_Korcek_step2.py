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
import numpy as np

from vector import *
from qc import *
from constants import *
from kinbot import *
from stationary_pt import *
from geom import *
from qc import *
from reac_family import *
import par



def do_Korcek_step2(species, instance, step, instance_name):
    """
    Carry out the reaction.
    """
    
    if step > 0:
        if check_qc(instance_name) != 'normal' and check_qc(instance_name) != 'error': return step

    #maximum number of steps for this reaction family
    max_step = 12

    # have a look at what has been done and what needs to be done
    step,geom, kwargs = initialize_reaction(species, instance, step, instance_name, max_step)

    #the the constraints for this step
    step, fix, change, release = get_constraints_Korcek_step2(step, species, instance,geom)
    
    #carry out the reaction and return the new step
    return carry_out_reaction(species, instance, step, instance_name, max_step, geom, kwargs, fix, change, release)

def get_constraints_Korcek_step2(step,species,instance,geom):
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
        
        if step < 10:
            dih = calc_dihedral(geom[instance[-4]], geom[instance[-3]], geom[instance[-2]], geom[instance[-1]])[0]
            if np.abs(dih) < 160:
                #move the dihedral to 160 degrees in 10 steps
                frac = 1. / (10 - step + 0.)
                new_dih = dih + frac * (160. - dih)
                constraint = [instance[-4] + 1,instance[-3] + 1,instance[-2] + 1,instance[-1] + 1,new_dih]
                change.append(constraint)
        
        fval = [2.0,2.0,1.8,1.8]
        if par.atom[instance[-1]] == 'H':
            fval[2] = 1.35
            fval[3] = 1.35
        
        val = new_bond_length(species,instance[0],instance[1],step+1,12,fval[0],geom)
        constraint = [instance[0] + 1,instance[1] + 1,val]
        change.append(constraint)
        
        val = new_bond_length(species,instance[2],instance[3],step+1,12,fval[1],geom)
        constraint = [instance[2] + 1,instance[3] + 1,val]
        change.append(constraint)
        
        if species.bond[instance[-1]][instance[-2]] == 1:
            val = new_bond_length(species,instance[-1],instance[-2],step+1,12,fval[2],geom)
            constraint = [instance[-1] + 1,instance[-2] + 1,val]
            change.append(constraint)
        #else do not change this bond length, the bond needs to stay and just change in order
        
        val = new_bond_length(species,instance[-1],instance[3],step+1,12,fval[3],geom)
        constraint = [instance[-1] + 1,instance[3] + 1,val] #todo: larger rings, this only work for 5 membered rings
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
