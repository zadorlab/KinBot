import numpy as np
import copy
import time

from kinbot import cheminfo
from kinbot import constants
from kinbot import reac_family
from kinbot import geometry

class Combinatorial:
    """
    Carry out the reaction.
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
        self.reac = instance[0]
        self.prod = instance[1]
        self.ts_bond = instance[2]
        # the position of the ts can be either
        # 0, 1 or 2, corresponding to early, mid
        # and late transition states
        self.position = instance[3]
        self.fvals = []
        self.get_final_vals()
        
        #name of the reaction
        self.instance_name = instance_name
        
        #maximum number of steps for this reaction family
        self.max_step = 20
        #do a scan?
        self.scan = 0
        #skip the first 12 steps in case the instance has a length of 3?
        self.skip = 0
        
        self.get_expected_products()

    def get_final_vals(self):
        """
        Method to get the final values of the initial ts geometry
        1. All bonds that need to be broken
        2. All bonds that need to be formed
        """
        fdists = {}
        fdists['CC'] = [1.60, 1.80, 2.16]
        fdists['CO'] = [1.50, 1.96, 2.62]
        fdists['CH'] = [1.20, 1.31, 1.50]
        fdists['OO'] = [1.69, 1.78, 1.92]
        fdists['HO'] = [0.99, 1.14, 1.38]
        fdists['HH'] = [0.75, 0.90, 1.10]
        fdists['CS'] = [2.18, 2.18, 2.18]
        fdists['OS'] = [2.18, 2.18, 2.18]
        fdists['HS'] = [1.60, 1.60, 1.60]
        fdists['SS'] = [2.48, 2.48, 2.48]

        if self.prod[0]:
            for pi in self.prod:
                i = pi[0]
                j = pi[1]
                syms = ''.join(sorted(self.species.atom[i]+self.species.atom[j]))
                if self.species.bond[i][j] == 0:
                    fdist = constants.st_bond[syms]
                    if syms in fdists:
                        fdist = list(reversed(fdists[syms]))[self.position]
                    self.fvals.append([i, j, fdist])

        if self.reac[0]:
            for ri in self.reac:
                i = ri[0]
                j = ri[1]
                syms = ''.join(sorted(self.species.atom[i]+self.species.atom[j]))
                if self.species.bond[i][j] == 1:
                    fdist = constants.st_bond[syms]
                    if syms in fdists:
                        fdist = fdists[syms][self.position]
                    self.fvals.append([i, j, fdist])

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
        #~ if step < self.max_step - 1:
            #~ # fix all the bond lengths
            #~ for i in range(self.species.natom - 1):
                #~ for j in range(i+1, self.species.natom):
                    #~ if self.species.bond[i][j] > 0:
                        #~ fix.append([i+1,j+1])
        if step < self.max_step:
            for fval in self.fvals:
                ndist = geometry.new_bond_length(self.species,
                                                 fval[0],
                                                 fval[1],
                                                 step,
                                                 self.max_step-1,
                                                 fval[2],
                                                 geom)
                constraint = [fval[0] + 1, fval[1] + 1, ndist]
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

    def get_expected_products(self):
        """
        Make smiles and inchis of the expected products
        """
        self.product_bond = np.zeros((self.species.natom, self.species.natom))
        for i in range(self.species.natom - 1):
            for j in range(i + 1, self.species.natom):
                if [i, j] in self.reac or [j, i] in self.reac:
                    self.product_bond[i][j] = self.species.bond[i][j] - 1
                elif [i, j] in self.prod or [j, i] in self.prod:
                    self.product_bond[i][j] = self.species.bond[i][j] + 1
                else:
                    self.product_bond[i][j] = self.species.bond[i][j]
                self.product_bond[j][i] = self.product_bond[i][j]
        rdmol, smi = cheminfo.create_rdkit_mol(self.product_bond, self.species.atom)
        self.prod_smi = smi.split('.')
        self.prod_inchi = []
        for smi in self.prod_smi:
            self.prod_inchi.append(cheminfo.create_inchi_from_smi(smi))

    def get_final_inchis(self):
        inchis = []
        for opt in self.prod_opt:
            species = opt.species
            inchi = cheminfo.create_inchi('', '', 'xyz/{}.xyz'.format(species.chemid))
            inchis.append(inchi)
        return inchis
