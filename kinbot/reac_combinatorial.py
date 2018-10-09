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

import constants
import reac_family
import geometry

class  Combinatorial:
    """
    Carry out the reaction.
    first 12 steps: curl up to form first bond
    steps 12 to 23: curl up to form second bond
    steps 24 to 35: curl up to form third bond
    step 36: relax the dihedrals, fix forming bonds
    steps 37 to 45: change the bond lengths
    step 46: saddle point optimization
    step 47: Done.
    
    The self.instance looks like: [[pr_bonds],[re_bonds]] with a mol_bonds [bond],[],[] and bond being two indices
    """

    def __init__(self,species,qc,par,instance,instance_name):
        #st_pt of the reactant
        self.species = species
        #st_pt of the ts
        self.ts = None
        #st_pt of the product(s)
        self.products = []
        #bond matrix of the products
        self.product_bonds = [] 
        
        #optimization objects
        self.ts_opt = None
        self.prod_opt = []
        
        self.qc = qc
        self.par = par
        
        #indices of the reactive atoms
        self.instance = instance
        #name of the reaction
        self.instance_name = instance_name
        
        #maximum number of steps for this reaction family
        self.max_step = 46
        #do a scan?
        self.scan = 0
        #skip the first 12 steps in case the instance has a length of 3?
        self.skip = 0

    def get_constraints(self,step, geom):
        """
        There are three types of constraints:
        1. fix a coordinate to the current value
        2. change a coordinate and fix is to the new value
        3. release a coordinate (only for gaussian)
        """
        fix = []
        change = []
        release = []
        if step < self.max_step:
            #fix all the bond lengths
            for i in range(self.species.natom - 1):
                for j in range(i+1, self.species.natom):
                    if self.species.bond[i][j] > 0:
                        fix.append([i+1,j+1])
        else:
            step = 12
        if step > 11 and step < 24:
            additional_bond_fix = [[self.instance[0][0][0] + 1, self.instance[0][0][1] + 1]]
            err, geom = read_qc_geom(self.instance_name, self.species.natom)

            if all([all([geom[i][j] == 0 for j in range(3)]) for i in range(self.species.natom)]):
                geom = self.species.geom
            
            #relax the previous dihedral fixing
            bool_old,inst_old = find_inst(self.species, self.instance[0][0])
            if bool_old and len(inst_old) > 3:
                for dih in range(len(inst_old)-3):
                    for i in range(4):
                        f.write(str(inst_old[dih+i] + 1) + ' ') 
                    f.write(' A\n')
            
            #generate the new dihedral fixing
            bool,inst = find_inst(self.species, self.instance[0][1])
            if bool and len(inst) > 3:
                dihedral_diff = init_ring_dihedral(self.species, inst, geom)
                for dih in range(len(inst)-3):
                    for i in range(4):
                        f.write(str(inst[dih+i] + 1) + ' ') 
                    f.write('+= ' + str(dihedral_diff[dih] / (23. - step + 1.) ) + ' F\n')
                for b in additional_bond_fix:
                    f.write(str(b[0]) + ' ' + str(b[1]) + ' F\n')
            else:
                step = 24
        if step > 23 and step < 36:
            additional_bond_fix = [[self.instance[0][0][0] + 1, self.instance[0][0][1] + 1]]
            additional_bond_fix.append([self.instance[0][1][0] + 1, self.instance[0][1][1] + 1])
            err, geom = read_qc_geom(self.instance_name, self.species.natom)
            
            if all([all([geom[i][j] == 0 for j in range(3)]) for i in range(self.species.natom)]):
                geom = self.species.geom
                
            #relax the previous dihedral fixing
            bool_old,inst_old = find_inst(self.species, self.instance[0][1])
            if bool_old and len(inst_old) > 3:
                for dih in range(len(inst_old)-3):
                    for i in range(4):
                        f.write(str(inst_old[dih+i] + 1) + ' ') 
                    f.write(' A\n')

            #generate the new dihedral fixing
            bool,inst = find_inst(self.species, self.instance[0][2])
            if bool and len(inst) > 3:
                dihedral_diff = init_ring_dihedral(self.species, inst, geom)
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
            bool_old,inst_old = find_inst(self.species, self.instance[0][2])
            if bool_old and len(inst_old) > 3:
                for dih in range(len(inst_old)-3):
                    for i in range(4):
                        f.write(str(inst_old[dih+i] + 1) + ' ') 
                    f.write(' A\n')
            
            if len(self.instance[0][0]) == 2:
                additional_bond_fix = [[self.instance[0][0][0] + 1, self.instance[0][0][1] + 1]]
            if len(self.instance[0][1]) == 2:
                additional_bond_fix.append([self.instance[0][1][0] + 1, self.instance[0][1][1] + 1])
            if len(self.instance[0][2]) == 2:
                additional_bond_fix.append([self.instance[0][2][0] + 1, self.instance[0][2][1] + 1])
            
            for b in additional_bond_fix:
                if not b[0] == b[1]:
                    if self.species.bond[b[0]-1][b[1]-1] == 0:
                        f.write(str(b[0]) + ' ' + str(b[1]) + ' F\n')
        if step > 36 and step < 46: # iterativelly change the bond lengths to their desired value
            err, geom = read_qc_geom(self.instance_name, self.species.natom)
            if all([all([geom[i][j] == 0 for j in range(3)]) for i in range(self.species.natom)]):
                geom = self.species.geom
            
            if step == 37:
                f_geom = open(self.instance_name + '_reac_geom.xyz','w')
                f_geom.write('%i\n\n'%self.species.natom)
                for i,at in enumerate(self.species.atom):
                    f_geom.write('%s\t%.8f\t%.8f\t%.8f\n'%(at,geom[i][0],geom[i][1],geom[i][2]))
                f_geom.write('\n\n')
                f_geom.close()

                
            final_dist = []
            for form_bond in self.instance[0]:
                if len(form_bond) == 2 and form_bond[0] != form_bond[1]:
                    if self.species.bond[form_bond[0],form_bond[1]] == 0: #else do not change that bond
                        fdist = constants.st_bond[''.join(sorted(self.species.atom[form_bond[0]]+self.species.atom[form_bond[1]]))] #final distance
                        cdist = np.linalg.norm(geom[form_bond[0]] - geom[form_bond[1]]) #current distance
                        dist_diff = (fdist - cdist )/(45 - step + 1) #distance increment
                        f.write(str(form_bond[0] + 1) + ' ' + str(form_bond[1] + 1) + ' += ' + str(dist_diff) + ' F\n')
            for break_bond in self.instance[1]:
                if len(break_bond) == 2 and break_bond[0] != break_bond[1]:
                    if self.species.bond[break_bond[0],break_bond[1]] == 1: #else do not change that bond
                        #check if bond should be breaking
                        br = 1
                        for fb in self.instance[0]:
                            if sorted(fb) == sorted(break_bond):
                                br = 0
                        if br:
                            fdist = constants.st_bond[''.join(sorted(self.species.atom[break_bond[0]]+self.species.atom[break_bond[1]]))] #final distance
                            cdist = np.linalg.norm(geom[break_bond[0]] - geom[break_bond[1]]) #current distance
                            dist_diff = (fdist - cdist)/(45 - step + 1) #distance increment
                            f.write(str(break_bond[0] + 1) + ' ' + str(break_bond[1] + 1) + ' += ' + str(dist_diff) + ' F\n')
            
        f.write('\n\n')                                    
        f.close()
       
        step += submit_qc(self.instance_name, 0)

        return step

    def find_inst(self, bond):
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
                    self.instances = start_motif(motif, self.species.natom, self.species.bond, self.species.atom, -1, [[i] for i in range(self.species.natom)])
                    for ins in self.instances: 
                        if bond[0] == ins[0] and bond[1] == ins[-1]:
                            return 1, ins
                        if bond[1] == ins[0] and bond[0] == ins[-1]:
                            return 1, ins
                    size += 1
                    if size > self.species.natom:
                        logging.error('Could not find path between atoms %i and %i'%(bond[0],bond[1]))
                        return 0,0
        else:
            # forming a radical
            return 0, 0
        