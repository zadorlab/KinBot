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
from __future__ import print_function, division
import sys
import os
import numpy as np
import re
import logging
import copy
import math

import cheminfo
import constants
import find_motif
import geometry

class StationaryPoint:
    """
    This object contains the properties of wells.
    """

    def __init__(self, name, charge, mult, smiles = '', structure = None, natom = 0, atom = None,geom = None, wellorts = 0):
        self.name = name
        self.mult = mult
        self.charge = charge
        self.short_name = '' #name for the MESS calculations needs to be shorter
        if geom is None:
            self.geom = []
        else:
            self.geom = geom
        if structure is None:
            self.structure = []
        else:
            self.structure = structure
        self.smiles = smiles
        self.obmol = None
        self.inchi = ''
        self.natom = natom
        if atom is None:
            self.atom = []
        else:
            self.atom = atom
        self.wellorts = wellorts

        self.energy = 0.
        self.zpe = 0.
        self.elec = 1.
        self.freq = [] # frequencies calculated by the qc program
 
        self.symm = 1.
        self.rot = []
        self.rads = [] # unique list of radical centers in case of resonance
        self.bonds = [] # unique list of bond matrices in case of resonance

        self.reac_type = []
        self.reac_inst = [] # holds the key atoms
        self.reac_obj = [] #instances of the reac_name objects
        self.reac_name = [] # holds the file name
        self.reac_step = []
        self.reac_ts_done = []
        self.reac_ts_geom = []
        self.reac_ts_freq = []
        self.reac_scan_energy = []
        
        #Instance of HomolyticScissions class
        self.homolytic_scissions = None
        
        #Instance of HIR class
        self.hir = None
        
        #Instance of the Conformers class
        self.confs = None
        
        #symmetry numbers
        self.sigma_ext = -1 #extermal symmetry number
        self.sigma_int = [] #internal summetry number around each atom
        self.nopt = -1 #number of optical isomers
        
        # frequencies calculated by kinbot
        self.kinbot_freqs = []
        self.reduced_freqs = []
        
        if len(self.geom) == 0:
            self.get_geom()
        if self.natom == 0:
            self.natom = len(atom)

    def get_geom(self):
        """
        Method reads the structure and converts it to a geometry
        or converts the smiles to a geometry using open babel
        """
        if len(self.structure) == 0:
            self.obmol,self.structure,self.bond = cheminfo.generate_3d_structure(self.smiles)
        self.natom = len(self.structure) // 4
        self.structure = np.reshape(self.structure, (self.natom,4))
        self.atom = self.structure[:,0]
        self.geom = self.structure[:,1:4].astype(float)

    def characterize(self):
        """
        With one call undertake a typical set of structural characterizations.
        """
        self.bond_mx()
        self.find_conf_dihedral()
        self.find_atom_eqv()

    def distance_mx(self):
        """ 
        Create a distance matrix 
        """
        self.dist = np.zeros((self.natom, self.natom))    
        for i in range(self.natom):
            for j in range(self.natom):
                self.dist[i][j] = np.linalg.norm(self.geom[i] - self.geom[j])
        return 0 

        
    def bond_mx(self):
        """ 
        Create bond matrix 
        """
        if not hasattr(self,'bond'): 
            self.distance_mx()
            self.bond = np.zeros((self.natom, self.natom), dtype=int)
        
            for i in range(self.natom):
                for j in range(self.natom):
                    if i == j: continue
                    if self.dist[i][j] < constants.st_bond[''.join(sorted(self.atom[i]+self.atom[j]))]:
                        self.bond[i][j] = 1
        else:
            for i in range(self.natom):
                for j in range(self.natom):
                    if i == j: continue
                    if self.bond[i][j] > 0:
                        self.bond[i][j] = 1
        
        max_bond = [constants.st_bond[self.atom[i]] for i in range(self.natom)]

        n_bond = np.sum(self.bond, axis=0)

        save_n_bond = n_bond
        self.rad = max_bond - n_bond

        # create all the permutations of the heavy atoms
        rad_atoms = [i for i in range(self.natom) if self.rad[i] > 0]
        all_permutations = False 
        if all_permutations: #use all the permutations (slow for more than 6 atoms in conjugated system)
            perms = list(itertools.permutations(rad_atoms))
        else: # use the same atom ordering but a different starting atoms and searching directions
            new_algo = False
            if new_algo:
                perms = []
                for i in range(1000):
                    perm = np.ndarray.tolist(np.random.permutation(rad_atoms))
                    perms.append(perm)
            else:
                perms = []
                if len(rad_atoms) > 0:
                    for index in range(len(rad_atoms)):
                        list1 = np.ndarray.tolist(np.roll(rad_atoms, index)) #forward search
                        list2 = np.ndarray.tolist(np.roll(list1[::-1], 1)) #reverse search
                        perms.append(list1)
                        perms.append(list2)
                else:
                    perms.append([0])
        
        # create lists to save the results of all the searches
        perm_bond = []
        perm_rad = []
        
        for index,perm in enumerate(perms): # iterate the permutations
            # copy the objects of the molecule into temporary objects for this search
            perm_bond.append(copy.deepcopy(self.bond))
            perm_rad.append(np.copy(self.rad))
            perm_n_bond= np.copy(n_bond)
            perm_save_n_bond = np.copy(save_n_bond)
            while np.sum(perm_rad[index]) > 1:
                for ind1 in range(len(perm)):
                    i = perm[ind1]
                    if self.atom[i] == 'S':
                        if perm_rad[index][i]%2 == 0:
                            perm_rad[index][i] = 0
                    for ind2 in range(ind1, len(perm)):   
                        j = perm[ind2]
                        if perm_rad[index][i] > 0 and perm_rad[index][j] > 0 and perm_bond[index][i][j] > 0:
                            incr = 1
                            if perm_rad[index][i] == 2 and perm_rad[index][j] == 2:
                                incr = 2
                            perm_bond[index][i][j] += incr
                            perm_bond[index][j][i] += incr
                            perm_n_bond[i] += incr
                            perm_n_bond[j] += incr
                            perm_rad[index][i] -= incr
                            perm_rad[index][j] -= incr

                if perm_n_bond.all == perm_save_n_bond.all: 
                    # bond orders do not change anymore
                    # check for sulfur atoms, if rad == 2 or 4, bring it back to zero
                    for i,at in enumerate(self.atom):
                        if at == 'S':
                            if perm_rad[index][i] > 1:
                                if perm_rad[index][i]%2 == 0:
                                    perm_rad[index][i] = 0
                                else:
                                    perm_rad[index][i] = 1
                    break 
                else: 
                    perm_save_n_bond = perm_n_bond
        
        #get the permutation that leads to the lowest number of radicals
        if len(perm_rad) > 0:
            tot_rad_sum = [np.sum(x) for x in perm_rad] # total number of radicals in the molecule
            value,idx = min((val,i) for (i,val) in enumerate(tot_rad_sum)) 
            
            #take a "random" bond matrix, corresponding to the lowest number of radical centers
            #as the standard bond matrix for this stationary point
            self.bond = perm_bond[idx] 
            self.rad = perm_rad[idx]
        
        # collect all the resonance isomers 
        for i,perm_b in enumerate(perm_bond):
            #only consider the resonance structure with the minimum number of radical centers
            if np.sum(perm_rad[i]) == value: 
                #check the uniqueness of the rad vector
                is_unique = 1
                for r in self.rads:
                    if all([perm_rad[i][j] == r[j] for j in range(self.natom)]):
                        is_unique = 0
                if is_unique:
                    self.rads.append(perm_rad[i])
                    self.bonds.append(perm_bond[i])
                else:
                    #check the uniqueness of the bond matrix
                    is_unique = 1
                    for b in self.bonds:
                        if all([all([b[j][k] == perm_b[j][k] for k in range(self.natom)]) for j in range(self.natom)]):
                            is_unique = 0
                    if is_unique:
                        self.rads.append(perm_rad[i])
                        self.bonds.append(perm_bond[i])
        
        return 0


    def calc_multiplicity(self,atomlist):
        """ 
        Calculate the multiplicity based on atom types.
        Returns the lowest multiplicity possible, i.e., singlet or dublet,
        Gaussian style.
        """

        if all([element == 'O' for element in atomlist]):
            return 3 # O and O2 are triplet
        if len(atomlist) == 1 and atomlist[0] == 'C':
            return 3 # C atom is triplet?

        mult = 0
        for element in atomlist:
            mult += constants.st_bond[element]
            
        return 1 + mult % 2
        

    def start_multi_molecular(self):
        """
        Iterative method to find all the separate products from a bond matrix
        """
        self.bond_mx()
        
        bond = copy.deepcopy(self.bond)
        
        max_step = 1000
        status = [0 for i in range(self.natom)] # 1: part of a molecule, 0: not part of a molecule
        atoms = [i for i in range(self.natom)]
        mols = [] # list of lists with atom indices of the separate molecules
        atomlist = np.asarray(self.atom)
        
        while 1:
            if any([status[i] == 0 for i in range(len(status))]):
                # reduce the bond matrix to the atoms that have a 0 as status
                bondi = [[bond[i][j] for j in range(len(status)) if status[j] == 0] for i in range(len(status)) if status[i] == 0]
                at = [atoms[i] for i in range(len(status)) if status[i] == 0]
                natomi = len(at)
                fragi = [0 for i in range(self.natom)]
                if natomi == 1:
                    #this is a molecule containing only one atom
                    fragi[at[0]] = 1
                    atomi = [at[0]]
                    bool = 0
                else:
                    natomi = len(at)
                    bool, sta = self.extract_next_mol(natomi,bondi)
                    atomi = [at[i] for i in range(natomi) if sta[i] == 1]
                    for i in range(len(sta)):
                        if sta[i] == 1:
                            status[at[i]] = 1
                            fragi[at[i]] = 1
                    
                if not bool and len(mols) == 0:
                    #the bond matrix corresponds to one molecule only
                    self.calc_chemid()
                    self.name = str(self.chemid)
                    mols.append(self)
                    break
                geomi = np.asarray(self.geom)[np.where(np.asarray(fragi) == 1)]
                natomi = np.sum(fragi)
                atomi = atomlist[np.where(np.asarray(fragi) == 1)]
                multi = self.calc_multiplicity(atomi)
                chargei = self.charge # todo
                moli = StationaryPoint('prod_%i'%(len(mols)+1),chargei,multi,atom=atomi,natom=natomi,geom=geomi)
                moli.characterize()
                moli.calc_chemid()
                moli.name = str(moli.chemid)

                mols.append(moli)
 
                if bool:
                    continue 
                else:
                    #reached the end, return the molecules
                    break
        return mols

   

    def extract_next_mol(self, natom, bond):
        """
        Test if a structure is bimolecular or one complex.
        Strategy: start walking from an arbitrary atom, and if not all are visited when all
        bonds available are walked, it is two or more separate fragments.
        It is a similar algorithm to which finds cycles.
        """
        
        max_step = 1000
        status = [0 for i in range(natom)] # 1: parent, -1: child, 0: not yet checked
        edge = [0 for i in range(natom*(natom-1)//2)] # the ordered list of edges in the graph, 0-1, 0-2, ..., 1-2, 1-3, ...
        chain = [] # steps made during the search

        n_connect = np.count_nonzero(bond) // 2 # number of atom connections

        visited_all = 0 # 1 if all are visited

        k = 0 # number of steps we made in the tree
        i = 0 # starting atom, arbitrary 
        chain = chain + [i]
        
        status[i] = -1
        
        while 1:
            if visited_all == n_connect or k == max_step: break
            j = 0
            while j < natom:
                if visited_all == n_connect or k == max_step: break
                # position of edge i, j in the ordered list
                if i < j: n = (natom - 1) * i - i * (i + 1) // 2 + j - 1 
                if i > j: n = (natom - 1) * j - j * (j + 1) // 2 + i - 1;   
                if i == j: 
                    j += 1
                else:
                    # they are connected and we did not go along that edge before 
                    if bond[i][j] > 0 and edge[n] == 0:
                        visited_all = visited_all + 1 
                        edge[n] = 1; 
                        status[i] = 1 # parent, the first atom, even if leaf, will become a parent.
                        status[j] = -1 # child
                        k = k + 1 
                        chain = chain + [j] 
                        i = j # now we'll start from j
                        j = 0 #restart the search for neighbors from 0
                    else:
                        j += 1

            # retract
            chain = chain[:-1]
            if len(chain) == 0:
                break
            i = chain[-1] # step back on the chain
            k = k + 1
        
        if any([status[i] == 0 for i in range(len(status))]):
            return 1, [abs(si) for si in status]
        else:
            return 0, [abs(si) for si in status]

    def find_cycle(self):
        """
        Find all the cycles in a molecule, if any
        This is done by searching from motifs ['X','X', ..., 'X']
        with length 3 to natom, and the cycles are defined by 
        the motif instances of which the first and last atom are bonded
        
        The search is halted before reaching natoms if a certain morif length 
        does not give any hit
        
        
        TODO: leave all the leaves of the graph out for the search, i.e.
        the atoms that only have neighbor, as they never participate in a cycle
        
        The cycles are keps in the cycle_chain list, which is a list of lists
        Its lists contain the atom indices participating in each cycle.
        
        In the case of fused cycles, keep all the possible cycles (e.g. two fused
        rings lead to three cycles, and they are all defined in the cycle_chain
        """
        
        self.cycle_chain = [] #list of the cycles
        self.cycle = [0 for i in range(self.natom)] # 0 if atom is not in cycle, 1 otherwise
        
        for cycle_size in range(3,self.natom+1):
            motif = ['X' for i in range(cycle_size)]
            instances = find_motif.start_motif(motif, self.natom, self.bond, self.atom, -1, [[k] for k in range(self.natom)])
            if len(instances) == 0:
                break
            for ins in instances:
                if self.bond[ins[0]][ins[-1]]:
                    #cycle found, check if it is new
                    new = 1
                    for cyc in self.cycle_chain:
                        if sorted(cyc) == sorted(ins):
                            new = 0
                    if new:
                        self.cycle_chain.append(ins)
                        for at in ins:
                            self.cycle[at] = 1
        return 0

    def calc_chemid(self):
        """ 
        The total id for a species.
        It is the sum of the atomids, plus a number for the multiplicity, Gaussian style.
        """
        self.chemid = int(0)
        #self.atomid = np.zeros(natom, dtype=long)
        self.atomid = [int(0) for i in range(self.natom)]
                      
        for i in range(self.natom):
            self.start_id(i) 
        
        for i in range(self.natom):
            self.chemid += self.atomid[i]
        self.chemid *= 10
        self.chemid += self.mult

        return 0
        
                                                                                                
                                                                                                                        
    def start_id(self, i):
        """ 
        Initialize recursive loop for id.
        i is the index for the atom to start at.
        """
        
        
        visit = [0 for k in range(self.natom)]
        depth = 0
        atomid = int(0)
        
        self.atomid[i], visit = self.calc_atomid(visit, depth, i, atomid)
        #a, visit = self.calc_atomid(visit, depth, i, atomid, natom, atom)

        #self.atomid[i] = a

        
        return 0
                                
                                
                                
    def calc_atomid(self, visit, depth, i, atomid):
        """ Caclulate chemical ID for a given atom. """        
        
        if not hasattr(self,'bond'): 
            #recalculate the bond matrix only if it is not there yet
            self.bond_mx()
            
        maxdepth = 7
        digit = 3
        if depth == maxdepth: return atomid, visit

        atomid += constants.mass[self.atom[i]] * int(math.pow(10, digit * (maxdepth - 1 -depth)))
        
        visit[i] = 1
        
        for j in range(self.natom):
            if self.bond[i][j] > 0:
                if visit[j] == 0:
                    depth += 1
                    atomid, visit = self.calc_atomid(visit, depth, j, atomid)
                    depth -= 1
                    visit[j] = 0

        return atomid, visit
    
    def find_dihedral(self): 
        """ 
        Identify unique rotatable bonds in the structure 
        No rotation around ring bonds and double and triple bonds.
        """
        
        if not hasattr(self,'chemid'):
            #only calculate this if it has not been calculated before
            self.calc_chemid()
        if not hasattr(self,'cycle_chain'):
            self.find_cycle()
        if len(self.bonds) == 0:
            self.bonds = [self.bond]
        
        self.dihed = []
        hit = 0

        if self.natom < 4: return 0

        # a-b-c-d, rotation around b-c
        for b in range(self.natom):
            if hit == 1: hit = 0
            for c in range(b, self.natom):
                if hit == 1: hit = 0
                if all([bi[b][c]==1 for bi in self.bonds]) and self.cycle[b] * self.cycle[c] == 0:
                    for a in range(self.natom):
                        if hit == 1: break
                        if self.bond[a][b] == 1 and a != c:
                            for d in range(self.natom):
                                if hit == 1: break
                                if self.bond[c][d] == 1 and d != b:
                                    dihedral_angle, warning = geometry.calc_dihedral(self.geom[a], self.geom[b], self.geom[c], self.geom[d])
                                    if warning == 0:
                                        self.dihed.append([a, b, c, d])
                                        hit = 1 

        return 0

    def find_conf_dihedral(self):
        """
        Just keep those rotatable bonds that are needed for conformer search.
        This way we exclude things like methyl groups, t-butyl groups, etc.            
        The result is stored in self.conf_dihed.
        """
        
        
        self.find_dihedral()
        self.conf_dihed = []
        dihed_sideb = []
        dihed_sidec = []
        
        for rotbond in range(len(self.dihed)):
            start = 0
            for i in range(self.natom):
                if np.count_nonzero(self.bond[self.dihed[rotbond][1]]) == 2:
                    dihed_sideb.append(self.dihed[rotbond][:])
                    break
                if i != self.dihed[rotbond][2] and self.bond[self.dihed[rotbond][1]][i] > 0:
                    if start == 0: 
                        base = self.atomid[i]
                        start = 1
                    elif self.atomid[i] != base: 
                        dihed_sideb.append(self.dihed[rotbond][:])
                        break
                
                        
        for rotbond in range(len(self.dihed)):
            start = 0
            for i in range(self.natom):
                if np.count_nonzero(self.bond[self.dihed[rotbond][2]]) == 2:
                    dihed_sidec.append(self.dihed[rotbond][:])
                    break
                if i != self.dihed[rotbond][1] and self.bond[self.dihed[rotbond][2]][i] > 0:
                    if start == 0: 
                        base = self.atomid[i]
                        start = 1
                    elif self.atomid[i] != base: 
                        dihed_sidec.append(self.dihed[rotbond][:])
                        break
        
        for b in range(len(dihed_sideb)):
            for c in range(len(dihed_sidec)):
                if dihed_sideb[b] == dihed_sidec[c]: self.conf_dihed.append(dihed_sideb[b][:])
                
        return 0
                
                
                
    def find_atom_eqv(self):
        """
        Determines which atoms are equivalent.
        """

        #this list contains a list of each set of equivalent atoms
        self.atom_eqv = []
        for atomi in range(self.natom):
            new_list = 1
            for list in self.atom_eqv:
                for atomj in list:
                    if self.atomid[atomi] == self.atomid[atomj] and not atomi in list:
                        if self.rigid_along_path(atomi,atomj):
                            break
                        else:
                            list.append(atomi)
                            new_list = 0
                            break
            if new_list:
                self.atom_eqv.append([atomi])
        
        return 0


    def rigid_along_path(self,atomi, atomj):
        """
        Method finds the shortest path between two atoms and checks if any atom along that
        pathway is rigid. An atom is rigid if it is in a cycle or is doubly bonded to another atom
        which has more than one neighbor. 
        """
        
        if self.bond[atomi][atomj] > 0:
            if self.bond[atomi][atomj] > 1: #atoms are doubly bonded
                return 1
            elif self.cycle[atomi] == 1: #atoms are in a cycle
                return 1
            else:
                return 0
        
        for chain_length in range(3,self.natom):
            motif = ['X' for i in range(chain_length)]
            instances = find_motif.start_motif(motif, self.natom, self.bond, self.atom, -1, [[k] for k in range(self.natom)])
            if len(instances) == 0:
                break
            for ins in instances:
                if (ins[0] == atomi and ins[-1] == atomj) or (ins[0] == atomj and ins[-1] == atomi):
                    for at in ins[1:-1]:
                        if self.cycle[at] == 1:
                            return 1
                        elif 2 in self.bond[at]:
                            double_neigh = [i for i, x in enumerate(self.bond[at]) if x == 2]
                            for neigh in double_neigh:
                                if sum(self.bond[neigh]) > 2: # atom has at least on other neighbor
                                    return 1
                    return 0
        return 0

def main():
    """
    This is the main object of the code, the stationary point
    """

if __name__ == "__main__":
    main()
