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
from motif import *
from geom import *
from qc import *
from cheminfo import *
import par


def do_combinatorial(species, instance, step, instance_name):
    """
    Carry out the reaction.
    first 12 steps: curl up to form first bond
    steps 12 to 23: curl up to form second bond
    steps 24 to 35: curl up to form third bond
    step 36: relax the dihedrals, fix forming bonds
    steps 37 to 45: change the bond lengths
    step 46: saddle point optimization
    step 47: Done.
    
    The instance looks like: [[pr_bonds],[re_bonds]] with a mol_bonds [bond],[],[] and bond being two indices
    """
    
    if step > 0:
        if check_qc(instance_name) != 'normal' and check_qc(instance_name) != 'error': return step
    
    #maximum number of steps for this reaction family
    max_step = 46

    # have a look at what has been done and what needs to be done
    step,geom, kwargs = initialize_reaction(species, instance, step, instance_name, max_step, skip)

    #the the constraints for this step
    step, fix, change, release = get_combinatorial(step, species, instance,geom)
    
    #carry out the reaction and return the new step
    return carry_out_reaction(species, instance, step, instance_name, max_step, geom, kwargs, fix, change, release)
    

def get_combinatorial(step,species,instance,geom):
    """
    There are three types of constraints:
    1. fix a coordinate to the current value
    2. change a coordinate and fix is to the new value
    3. release a coordinate (only for gaussian)
    """
    fix = []
    change = []
    release = []

    if step < 12:
        bool,inst = find_inst(species, instance[0][0])
        
        if bool and len(inst) > 3:
            dihedral_diff = init_ring_dihedral(species, inst)
            for dih in range(len(inst)-3):
                for i in range(4):
                    f.write(str(inst[dih+i] + 1) + ' ') 
                f.write('+= ' + str(dihedral_diff[dih] / 12.) + ' F\n')
        else:
            step = 12
    if step > 11 and step < 24:
        additional_bond_fix = [[instance[0][0][0] + 1, instance[0][0][1] + 1]]
        err, geom = read_qc_geom(instance_name, par.natom)

        if all([all([geom[i][j] == 0 for j in range(3)]) for i in range(par.natom)]):
            geom = species.geom
        
        #relax the previous dihedral fixing
        bool_old,inst_old = find_inst(species, instance[0][0])
        if bool_old and len(inst_old) > 3:
            for dih in range(len(inst_old)-3):
                for i in range(4):
                    f.write(str(inst_old[dih+i] + 1) + ' ') 
                f.write(' A\n')
        
        #generate the new dihedral fixing
        bool,inst = find_inst(species, instance[0][1])
        if bool and len(inst) > 3:
            dihedral_diff = init_ring_dihedral(species, inst, geom)
            for dih in range(len(inst)-3):
                for i in range(4):
                    f.write(str(inst[dih+i] + 1) + ' ') 
                f.write('+= ' + str(dihedral_diff[dih] / (23. - step + 1.) ) + ' F\n')
            for b in additional_bond_fix:
                f.write(str(b[0]) + ' ' + str(b[1]) + ' F\n')
        else:
            step = 24
    if step > 23 and step < 36:
        additional_bond_fix = [[instance[0][0][0] + 1, instance[0][0][1] + 1]]
        additional_bond_fix.append([instance[0][1][0] + 1, instance[0][1][1] + 1])
        err, geom = read_qc_geom(instance_name, par.natom)
        
        if all([all([geom[i][j] == 0 for j in range(3)]) for i in range(par.natom)]):
            geom = species.geom
            
        #relax the previous dihedral fixing
        bool_old,inst_old = find_inst(species, instance[0][1])
        if bool_old and len(inst_old) > 3:
            for dih in range(len(inst_old)-3):
                for i in range(4):
                    f.write(str(inst_old[dih+i] + 1) + ' ') 
                f.write(' A\n')

        #generate the new dihedral fixing
        bool,inst = find_inst(species, instance[0][2])
        if bool and len(inst) > 3:
            dihedral_diff = init_ring_dihedral(species, inst, geom)
            for dih in range(len(inst)-3):
                for i in range(4):
                    f.write(str(inst[dih+i] + 1) + ' ') 
                f.write('+= ' + str(dihedral_diff[dih] / (35. - step + 1.) ) + ' F\n')
            for b in additional_bond_fix:
                f.write(str(b[0]) + ' ' + str(b[1]) + ' F\n')
        else:
            step = 36
    if step == 36: #relax the dihedrals and keep de bond lengths fixed
        #relax the previous dihedral fixing
        bool_old,inst_old = find_inst(species, instance[0][2])
        if bool_old and len(inst_old) > 3:
            for dih in range(len(inst_old)-3):
                for i in range(4):
                    f.write(str(inst_old[dih+i] + 1) + ' ') 
                f.write(' A\n')
        
        if len(instance[0][0]) == 2:
            additional_bond_fix = [[instance[0][0][0] + 1, instance[0][0][1] + 1]]
        if len(instance[0][1]) == 2:
            additional_bond_fix.append([instance[0][1][0] + 1, instance[0][1][1] + 1])
        if len(instance[0][2]) == 2:
            additional_bond_fix.append([instance[0][2][0] + 1, instance[0][2][1] + 1])
        
        for b in additional_bond_fix:
            if not b[0] == b[1]:
                if species.bond[b[0]-1][b[1]-1] == 0:
                    f.write(str(b[0]) + ' ' + str(b[1]) + ' F\n')
    if step > 36 and step < 46: # iterativelly change the bond lengths to their desired value
        err, geom = read_qc_geom(instance_name, par.natom)
        if all([all([geom[i][j] == 0 for j in range(3)]) for i in range(par.natom)]):
            geom = species.geom
        
        if step == 37:
            f_geom = open(instance_name + '_reac_geom.xyz','w')
            f_geom.write('%i\n\n'%par.natom)
            for i,at in enumerate(par.atom):
                f_geom.write('%s\t%.8f\t%.8f\t%.8f\n'%(at,geom[i][0],geom[i][1],geom[i][2]))
            f_geom.write('\n\n')
            f_geom.close()

            
        final_dist = []
        for form_bond in instance[0]:
            if len(form_bond) == 2 and form_bond[0] != form_bond[1]:
                if species.bond[form_bond[0],form_bond[1]] == 0: #else do not change that bond
                    fdist = st_bond[''.join(sorted(par.atom[form_bond[0]]+par.atom[form_bond[1]]))] #final distance
                    cdist = np.linalg.norm(geom[form_bond[0]] - geom[form_bond[1]]) #current distance
                    dist_diff = (fdist - cdist )/(45 - step + 1) #distance increment
                    f.write(str(form_bond[0] + 1) + ' ' + str(form_bond[1] + 1) + ' += ' + str(dist_diff) + ' F\n')
        for break_bond in instance[1]:
            if len(break_bond) == 2 and break_bond[0] != break_bond[1]:
                if species.bond[break_bond[0],break_bond[1]] == 1: #else do not change that bond
                    #check if bond should be breaking
                    br = 1
                    for fb in instance[0]:
                        if sorted(fb) == sorted(break_bond):
                            br = 0
                    if br:
                        fdist = st_bond[''.join(sorted(par.atom[break_bond[0]]+par.atom[break_bond[1]]))] #final distance
                        cdist = np.linalg.norm(geom[break_bond[0]] - geom[break_bond[1]]) #current distance
                        dist_diff = (fdist - cdist)/(45 - step + 1) #distance increment
                        f.write(str(break_bond[0] + 1) + ' ' + str(break_bond[1] + 1) + ' += ' + str(dist_diff) + ' F\n')
        
    f.write('\n\n')                                    
    f.close()
   
    step += submit_qc(instance_name, 0)

    return step

def find_inst(species, bond):
    """
    Find the shortest motif from the first atom to the last atom of the forming bond
    """
    
    if len(bond) == 2:
        if bond[0] == bond[1]:
            #forming a lone electron pair
            return 0,0
        else: 
            size = 2
            while 1: 
                motif = ['X' for i in range(size)]
                instances = start_motif(motif, par.natom, species.bond, par.atom, -1, [[i] for i in range(par.natom)])
                for ins in instances: 
                    if bond[0] == ins[0] and bond[1] == ins[-1]:
                        return 1, ins
                    if bond[1] == ins[0] and bond[0] == ins[-1]:
                        return 1, ins
                size += 1
                if size > par.natom:
                    logging.error('Could not find path between atoms %i and %i'%(bond[0],bond[1]))
                    return 0,0
    else:
        # forming a radical
        return 0, 0
    