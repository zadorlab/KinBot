import numpy as np
import copy
import time
import sys

from kinbot import reac_family
from kinbot import geometry

class IntraRHAddEndoR:
    
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
            final_dist = [2.0, 1.3, 1.35, 1.3]
            for i in range(len(self.instance)):
                if i == len(self.instance) - 1:
                    j = 0
                else:
                    j = i+1
                val = geometry.new_bond_length(self.species,self.instance[i],self.instance[j],step,12,final_dist[i],geom)
                constraint = [self.instance[i] + 1,self.instance[j] + 1,val]
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
