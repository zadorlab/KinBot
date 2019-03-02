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

from kinbot import reac_family
from kinbot import geometry

class KorcekStep2:
    
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
        self.max_step = 12
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
        if step < 12:
            
            if step < 10:
                dih = geometry.calc_dihedral(geom[self.instance[-4]], geom[self.instance[-3]], geom[self.instance[-2]], geom[self.instance[-1]])[0]
                if np.abs(dih) < 160:
                    #move the dihedral to 160 degrees in 10 steps
                    frac = 1. / (10 - step + 0.)
                    new_dih = dih + frac * (160. - dih)
                    constraint = [self.instance[-4] + 1,self.instance[-3] + 1,self.instance[-2] + 1,self.instance[-1] + 1,new_dih]
                    change.append(constraint)
            
            fval = [2.0,2.0,1.8,1.8]
            if self.species.atom[self.instance[-1]] == 'H':
                fval[2] = 1.35
                fval[3] = 1.35
            
            val = geometry.new_bond_length(self.species,self.instance[0],self.instance[1],step+1,12,fval[0],geom)
            constraint = [self.instance[0] + 1,self.instance[1] + 1,val]
            change.append(constraint)
            
            val = geometry.new_bond_length(self.species,self.instance[2],self.instance[3],step+1,12,fval[1],geom)
            constraint = [self.instance[2] + 1,self.instance[3] + 1,val]
            change.append(constraint)
            
            if self.species.bond[self.instance[-1]][self.instance[-2]] == 1:
                val = geometry.new_bond_length(self.species,self.instance[-1],self.instance[-2],step+1,12,fval[2],geom)
                constraint = [self.instance[-1] + 1,self.instance[-2] + 1,val]
                change.append(constraint)
            #else do not change this bond length, the bond needs to stay and just change in order
            
            val = geometry.new_bond_length(self.species,self.instance[-1],self.instance[3],step+1,12,fval[3],geom)
            constraint = [self.instance[-1] + 1,self.instance[3] + 1,val] #todo: larger rings, this only work for 5 membered rings
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
