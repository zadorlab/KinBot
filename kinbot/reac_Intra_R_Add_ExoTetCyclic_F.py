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

class IntraRAddExoTetCyclicF:
    
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
        self.max_step = 22
        #do a scan?
        self.scan = 0
        #skip the first 12 steps in case the instance has a length of 3?
        self.skip = 1

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
        if step < 12:
            new_dihs = geometry.new_ring_dihedrals(self.species, self.instance, step, 12,geom)
            for dih in range(len(self.instance)-4): # do not include last atom
                constraint = []
                for i in range(4):
                    constraint.append(self.instance[dih+i] + 1)
                constraint.append(new_dihs[dih])
                change.append(constraint)
            ldih = [] # constraint for the last dihedral, which needs to be 180 degrees    
            for i in range(4):
                ldih.append(self.instance[len(self.instance)-4+i] + 1)
            dih = geometry.calc_dihedral(geom[ldih[0] - 1], geom[ldih[1] - 1], geom[ldih[2] - 1], geom[ldih[3] - 1])[0]
            frac = 1./(12. - step)
            if dih < 0:
                new_dih = dih - frac * (180. + dih) 
                ldih.append(new_dih)
            else:
                new_dih = dih + frac * (180. - dih) 
                ldih.append(new_dih)
            change.append(ldih)
        elif step < 22:
            for dih in range(len(self.instance)-3):  
                constraint = []
                for i in range(4):
                    constraint.append(self.instance[dih+i] + 1)
                release.append(constraint)

            fdist1 = constants.st_bond[''.join(sorted(self.species.atom[self.instance[0]]+self.species.atom[self.instance[-2]]))]*1.0
            if ''.join(sorted(self.species.atom[self.instance[0]]+self.species.atom[self.instance[-2]])) == 'CO':
                fdist1 = 1.68
            ndist1 = geometry.new_bond_length(self.species,self.instance[0],self.instance[-2],step - 11 ,10,fdist1,geom)
            constraint = [self.instance[0] + 1,self.instance[-2] + 1,ndist1]
            change.append(constraint)
            
            fdist2 = constants.st_bond[''.join(sorted(self.species.atom[self.instance[-1]]+self.species.atom[self.instance[-2]]))]*1.0
            if ''.join(sorted(self.species.atom[self.instance[-1]]+self.species.atom[self.instance[-2]])) == 'CO':
                fdist2 = 1.68
            ndist2 = geometry.new_bond_length(self.species,self.instance[-1],self.instance[-2],step - 11 ,10,fdist2,geom)
            constraint = [self.instance[-1] + 1,self.instance[-2] + 1,ndist2]
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


