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
import sys
import os
import numpy as np
import re
import subprocess
import logging
import copy
import math

from qc import *
from vector import *
from motif import *
from reaction import *
from write_mess import *
from constants import *
from par import *
from bond_combinations import *


class stationary_pt(object):
    """
    This object contains the properties of wells.
    """

    def __init__(self, name):
        self.name = name
        self.short_name = '' #name for the MESS calculations needs to be shorter
        self.geom = np.empty((0,0))
        self.energy = 0.
        self.zpe = 0.
        self.symm = 1.
        self.elec = 1.
        self.brot = np.zeros((3))
        self.freq = [] # frequencies calculated by the qc program
        self.rot = []
        self.rads = [] # unique list of radical centers in case of resonance
        self.bonds = [] # unique list of bond matrices in case of resonance
        self.reac_type = []
        self.reac_inst = [] # holds the key atoms
        self.reac_name = [] # holds the file name
        self.reac_step = []
        self.reac_ts_done = []
        self.reac_ts_geom = []
        self.reac_ts_freq = []
        self.reac_scan_energy = []
        self.natom = -1
        self.mult = -1
        self.charge = 0
        self.atom = []
        # list -1 (not finished), 0 (successful) or 1 (failed) for each hir calculation
        # each element corresponds to one point along the scan
        self.hir_status = []
        self.hir_energies = []
        self.hir_fourier = []
        self.hir_geoms = []
        
        # -1 (not finished), 0 (successful) or 1 (failed) for each cyclic conformer
        self.cyc_conf_status = []
        # -1 (not finished), 0 (successful) or 1 (failed) for each open chain conformer
        self.conf_status = []
        
        #symmetry numbers
        self.sigma_ext = -1 #extermal symmetry number
        self.sigma_int = [] #internal summetry number around each atom
        self.nopt = -1 #number of optical isomers
        
        # frequencies calculated by kinbot
        self.kinbot_freqs = []
        self.reduced_freqs = []



    def characterize(self, natom, atom, mult, charge):
        """
        With one call undertake a typical set of structural characterizations.
        """


        self.find_conf_dihedral(natom, atom, mult)
        self.find_atom_eqv(natom, atom)


        
    def distance_mx(self, natom):
        """ Create a distance matrix """


        self.dist = np.zeros((natom, natom))    
        for i in range(natom):
            for j in range(natom):
                self.dist[i][j] = np.linalg.norm(self.geom[i] - self.geom[j])
  
        return 0 
        


    def bond_mx(self, natom, atom):
        """ Create bond matrix """


        self.distance_mx(natom)
        self.bond = np.zeros((natom, natom), dtype=int)
    
        for i in range(natom):
            for j in range(natom):
                if i == j: continue
                if self.dist[i][j] < st_bond[''.join(sorted(atom[i]+atom[j]))]:
                    self.bond[i][j] = 1
        
        max_bond = [st_bond[atom[i]] for i in range(natom)]

        n_bond = np.sum(self.bond, axis=0)

        save_n_bond = n_bond
        self.rad = max_bond - n_bond

        # create all the permutations of the heavy atoms
        rad_atoms = [i for i in range(natom) if self.rad[i] > 0]
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
                    if atom[i] == 'S':
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
                    for i,at in enumerate(atom):
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
                    if all([perm_rad[i][j] == r[j] for j in range(natom)]):
                        is_unique = 0
                if is_unique:
                    self.rads.append(perm_rad[i])
                    self.bonds.append(perm_bond[i])
                else:
                    #check the uniqueness of the bond matrix
                    is_unique = 1
                    for b in self.bonds:
                        if all([all([b[j][k] == perm_b[j][k] for k in range(natom)]) for j in range(natom)]):
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
            mult += st_bond[element]
            
        return 1 + mult % 2
        

    def start_multi_molecular(self, natom, atom):
        """
        Iterative method to find all the separate products from a bond matrix
        """
        self.bond_mx(natom, atom)
        
        bond = copy.deepcopy(self.bond)
        
        max_step = 1000
        status = [0 for i in range(natom)] # 1: part of a molecule, 0: not part of a molecule
        atoms = [i for i in range(natom)]
        mols = [] # list of lists with atom indices of the separate molecules
        atomlist = np.asarray(atom)
        
        while 1:
            if any([status[i] == 0 for i in range(len(status))]):
                # reduce the bond matrix to the atoms that have a 0 as status
                bondi = [[bond[i][j] for j in range(len(status)) if status[j] == 0] for i in range(len(status)) if status[i] == 0]
                at = [atoms[i] for i in range(len(status)) if status[i] == 0]
                natomi = len(at)
                fragi = [0 for i in range(natom)]
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
                    self.mult = par.mult
                    self.charge = par.charge
                    self.natom = par.natom
                    self.atom = par.atom
                    mols.append(self)
                    break
                geomi = np.asarray(self.geom)[np.where(np.asarray(fragi) == 1)]
                natomi = np.sum(fragi)
                atomi = atomlist[np.where(np.asarray(fragi) == 1)]
                multi = self.calc_multiplicity(atomi)
                chargei = par.charge # todo
                
                moli = stationary_pt('prod_%i'%(len(mols)+1))
                moli.geom = copy.copy(geomi)
                moli.characterize(natomi,atomi,multi,chargei)
                
                moli.mult = multi
                moli.charge = chargei
                moli.natom = natomi
                moli.atom = atomi
                
                moli.calc_chemid(natomi,atomi,multi)
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
        edge = [0 for i in range(natom*(natom-1)/2)] # the ordered list of edges in the graph, 0-1, 0-2, ..., 1-2, 1-3, ...
        chain = [] # steps made during the search

        n_connect = np.count_nonzero(bond) / 2 # number of atom connections

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
                if i < j: n = (natom - 1) * i - i * (i + 1) / 2 + j - 1 
                if i > j: n = (natom - 1) * j - j * (j + 1) / 2 + i - 1;   
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
            return 1, map(abs, status)
        else:
            return 0, map(abs, status)

    def find_cycle(self, natom, atom):
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
        self.cycle = [0 for i in range(natom)] # 0 if atom is not in cycle, 1 otherwise
        
        for cycle_size in range(3,natom):
            motif = ['X' for i in range(cycle_size)]
            instances = start_motif(motif, natom, self.bond, atom, -1, [[k] for k in range(natom)])
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
    
    def find_cycle_old(self, natom, atom):
        """ 
        Find cycles in the structure, if any.
        It the sysmplest version it just deals with one cycle.
        """         
        
        
        self.bond_mx(natom, atom)
        max_step = 10000
        status = [0 for i in range(natom)] # 1: parent, -1: child, 0: not yet checked
        edge = [0 for i in range(natom*(natom-1)/2)] # the ordered list of edges in the graph, 0-1, 0-2, ..., 1-2, 1-3, ...
        chain = [-999 for i in range(max_step)] # steps made during the search
        cyc = [-999 for i in range(max_step)] # untrimmed cycle variable        
        self.cycle = [0 for i in range(natom)] # final cycle variable        

        n_connect = np.count_nonzero(self.bond) / 2 # number of atom connections

        visited_all = 0 # 1 if all are visited
        n_cycle = 0 # number of cycles

        k = 0 # number of steps we made in the tree
        i = 0 # starting atom
        m = 0
        l = 0 
        chain[k] = i  
        notfound = 0 

        while 1:
            if visited_all == n_connect or k == max_step: break
            j = 0
            while j < natom:
                if visited_all == n_connect or k == max_step: break
                # position of edge i, j (i > j) in the ordered list
                if i < j: n = (natom - 1) * i - i * (i + 1) / 2 + j - 1 
                if i > j: n = (natom - 1) * j - j * (j + 1) / 2 + i - 1;   
                if i == j: 
                    j += 1
                else:
                    # they are connected and we did not go along that edge before 
                    if self.bond[i][j] > 0 and edge[n] == 0:
                        notfound = 0; 
                        visited_all = visited_all + 1 
                        edge[n] = 1; 
                        if status[j] == 1: # cycle found!
                            n_cycle = n_cycle + 1 
                            # need to enumerate the atoms in the cycle
                            m = k # number of steps so far - 1
                            if n_cycle == 1: cyc[0] = j #closing point in the cycle
                            if n_cycle == 1: l = 1 
                            while chain[m] != cyc[0]: # we are not back to the starting point
                                cyc[l] = chain[m]
                                l = l + 1
                                m = m - 1
                        status[i] = 1 # parent, the first atom, even if leaf, will become a parent.
                        status[j] = -1 # child
                        k = k + 1 
                        chain[k] = j 
                        i = j # now we'll start from j
                        j = 0
                    else:
                        j += 1

            # retract
            notfound = notfound + 1 
            i = chain[k - (2*notfound - 1)] # step back on the chain
            k = k + 1
            chain[k] = i
 
        # clean up cyc variable from side branches and put it in cycle variable
        for i in range(l+1): # this was changed to i<=l from i<l, now all atoms in cycle are enumerated
            for j in range(i+1, l):
                if cyc[i] == cyc[j]:
                    for m in range(i+1, j+1):
                        cyc[m] = -999; 
 
        self.cycle_chain = []
        for i in range(max_step): 
            if cyc[i] != -999: 
                self.cycle[cyc[i]] = 1 # the cyc[i]th element will be 1, and not the ith
                self.cycle_chain.append(cyc[i]) # atoms in cycle in order
        self.cycle_size = np.sum(self.cycle)
        
        return 0 


        
    def calc_chemid(self, natom, atom, mult):
        """ 
        The total id for a species.
        It is the sum of the atomids, plus a number for the multiplicity, Gaussian style.
        """
        
        
        self.chemid = long(0)
        #self.atomid = np.zeros(natom, dtype=long)
        self.atomid = [long(0) for i in range(natom)]
                      
        for i in range(natom):
            self.start_id(i, natom, atom) 
        
        for i in range(natom):
            self.chemid += self.atomid[i]
        self.chemid *= 10
        self.chemid += mult

        return 0
        
                                                                                                
                                                                                                                        
    def start_id(self, i, natom, atom):
        """ 
        Initialize recursive loop for id.
        i is the index for the atom to start at.
        """
        
        
        visit = [0 for k in range(natom)]
        depth = 0
        atomid = long(0)
        
        self.atomid[i], visit = self.calc_atomid(visit, depth, i, atomid, natom, atom)
        #a, visit = self.calc_atomid(visit, depth, i, atomid, natom, atom)

        #self.atomid[i] = a

        
        return 0
                                
                                
                                
    def calc_atomid(self, visit, depth, i, atomid, natom, atom):
        """ Caclulate chemical ID for a given atom. """        
        
        if not hasattr(self,'bond'): 
            #recalculate the bond matrix only if it is not there yet
            self.bond_mx(natom, atom)
            
        maxdepth = 7
        digit = 3
        if depth == maxdepth: return atomid, visit

        atomid += mass[atom[i]] * long(math.pow(10, digit * (maxdepth - 1 -depth)))
        
        visit[i] = 1
        
        for j in range(natom):
            if self.bond[i][j] > 0:
                if visit[j] == 0:
                    depth += 1
                    atomid, visit = self.calc_atomid(visit, depth, j, atomid, natom, atom)
                    depth -= 1
                    visit[j] = 0

        return atomid, visit
        
    
    
    def find_dihedral(self, natom, atom, mult): 
        """ 
        Identify unique rotatable bonds in the structure 
        No rotation around ring bonds and double and triple bonds.
        """
        
        if not hasattr(self,'chemid'):
            #only calculate this if it has not been calculated before
            self.calc_chemid(natom, atom, mult)
        if not hasattr(self,'cycle_chain'):
            self.find_cycle(natom, atom)
        
        self.dihed = []
        hit = 0

        if natom < 4: return 0

        # a-b-c-d, rotation around b-c
        for b in range(natom):
            if hit == 1: hit = 0
            for c in range(b, natom):
                if hit == 1: hit = 0
                if self.bond[b][c] == 1 and self.cycle[b] * self.cycle[c] == 0:
                    for a in range(natom):
                        if hit == 1: break
                        if self.bond[a][b] == 1 and a != c:
                            for d in range(natom):
                                if hit == 1: break
                                if self.bond[c][d] == 1 and d != b:
                                    dihedral_angle, warning = dihedral(self.geom[a], self.geom[b], self.geom[c], self.geom[d])
                                    if warning == 0:
                                        self.dihed.append([a, b, c, d])
                                        hit = 1 
                                    
        return 0
        
    
    
    def find_conf_dihedral(self, natom, atom, mult):
        """
        Just keep those rotatable bonds that are needed for conformer search.
        This way we exclude things like methyl groups, t-butyl groups, etc.            
        The result is stored in self.conf_dihed.
        """
        
        
        self.find_dihedral(natom, atom, mult)
        self.conf_dihed = []
        dihed_sideb = []
        dihed_sidec = []
        
        for rotbond in range(len(self.dihed)):
            start = 0
            for i in range(natom):
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
            for i in range(natom):
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
                
                
                
    def find_atom_eqv(self, natom, atom):
        """
        Determines which atoms are equivalent.
        """

        #this list contains a list of each set of equivalent atoms
        self.atom_eqv = []
        for atomi in range(natom):
            new_list = 1
            for list in self.atom_eqv:
                for atomj in list:
                    if self.atomid[atomi] == self.atomid[atomj] and not atomi in list:
                        if self.rigid_along_path(atomi,atomj,natom, atom):
                            break
                        else:
                            list.append(atomi)
                            new_list = 0
                            break
            if new_list:
                self.atom_eqv.append([atomi])
        
        return 0


    def rigid_along_path(self,atomi, atomj, natom, atom):
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
        
        for chain_length in range(3,natom):
            motif = ['X' for i in range(chain_length)]
            instances = start_motif(motif, natom, self.bond, atom, -1, [[k] for k in range(natom)])
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

    def find_reactions(self, natom, atom):
        """
        List all reaction types available, and find the key atoms for them 
        for the current structure.
        """
        
        #keys: names of the families
        #values: list of instances
        #this dict is used to keep track of the unique reactions found,
        #and to verify whether a new reaction is indeed unique 
        self.reactions = {}
        
       
        for i, bond in enumerate(self.bonds):
            rad = self.rads[i]
            
            
            if 'intra_H_migration' in par.families or 'all' in par.families:
                self.search_intra_H_migration(natom,atom,bond,rad)
                
            if 'intra_R_migration' in par.families or 'all' in par.families:
                self.search_intra_R_migration(natom,atom,bond,rad)
                
            if 'intra_OH_migration' in par.families or 'all' in par.families:
                self.search_intra_OH_migration(natom,atom,bond,rad)
                
            if 'cpd_H_migration' in par.families or 'all' in par.families:
                self.search_cpd_H_migration(natom,atom,bond,rad)
                
            if 'Intra_RH_Add_Endocyclic_F' in par.families or 'all' in par.families:
                self.search_Intra_RH_Add_Endocyclic_F(natom,atom,bond,rad)
                
            if 'Intra_RH_Add_Endocyclic_R' in par.families or 'all' in par.families:
                self.search_Intra_RH_Add_Endocyclic_R(natom,atom,bond,rad)
                
            if 'Cyclic_Ether_Formation' in par.families or 'all' in par.families:
                self.search_Cyclic_Ether_Formation(natom,atom,bond,rad)
                
            if 'Intra_RH_Add_Exocyclic_F' in par.families or 'all' in par.families:
                self.search_Intra_RH_Add_Exocyclic_F(natom,atom,bond,rad)
                
            if 'Intra_RH_Add_Exocyclic_R' in par.families or 'all' in par.families:
                self.search_Intra_RH_Add_Exocyclic_R(natom,atom,bond,rad)
                
            if 'Retro_Ene' in par.families or 'all' in par.families:
                self.search_Retro_Ene(natom,atom,bond,rad)
                
            if 'Intra_R_Add_Endocyclic_F' in par.families or 'all' in par.families:
                self.search_Intra_R_Add_Endocyclic_F(natom,atom,bond,rad)
                
            if 'Intra_R_Add_ExoTetCyclic_F' in par.families or 'all' in par.families:
                self.search_Intra_R_Add_ExoTetCyclic_F(natom,atom,bond,rad)
                
            if 'Intra_R_Add_Exocyclic_F' in par.families or 'all' in par.families:
                self.search_Intra_R_Add_Exocyclic_F(natom,atom,bond,rad)
                
            if 'Korcek_step2' in par.families or 'all' in par.families:
                self.search_Korcek_step2(natom,atom,bond,rad)
                
            if 'r22_cycloaddition' in par.families or 'all' in par.families:
                self.search_r22_cycloaddition(natom,atom,bond,rad)
                
            if 'r12_cycloaddition' in par.families or 'all' in par.families:
                self.search_r12_cycloaddition(natom,atom,bond,rad)
                
            if 'r12_insertion_R' in par.families or 'all' in par.families:
                self.search_r12_insertion_R(natom,atom,bond,rad)
                
            if 'r13_insertion_CO2' in par.families or 'all' in par.families:
                self.search_r13_insertion_CO2(natom,atom,bond,rad)
                
            if 'r13_insertion_ROR' in par.families or 'all' in par.families:
                self.search_r13_insertion_ROR(natom,atom,bond,rad)
                
            if 'Diels_alder_addition' in par.families or 'all' in par.families:
                self.search_Diels_alder_addition(natom,atom,bond,rad)
                
            if 'Intra_Diels_alder_R' in par.families or 'all' in par.families:
                self.search_Intra_Diels_alder_R(natom,atom,bond,rad)
                
            if 'ketoenol' in par.families or 'all' in par.families:
                self.search_ketoenol(natom,atom,bond,rad)
                
            if 'HO2_Elimination_from_PeroxyRadical' in par.families or 'all' in par.families:
                self.search_HO2_Elimination_from_PeroxyRadical(natom,atom,bond,rad)
                
            if 'R_Addition_COm3_R' in par.families or 'all' in par.families:
                self.search_R_Addition_COm3_R(natom,atom,bond,rad)
                
            if 'R_Addition_MultipleBond' in par.families or 'all' in par.families:
                self.search_R_Addition_MultipleBond(natom,atom,bond,rad)
                
            if '12_shift_S_F' in par.families or 'all' in par.families:
                self.search_12_shift_S_F(natom,atom,bond,rad)
                
            if '12_shift_S_R' in par.families or 'all' in par.families:
                self.search_12_shift_S_R(natom,atom,bond,rad)
                
            if 'R_Addition_CSm_R' in par.families or 'all' in par.families:
                self.search_R_Addition_CSm_R(natom,atom,bond,rad)
                
            if 'r13_insertion_RSR' in par.families or 'all' in par.families:
                self.search_r13_insertion_RSR(natom,atom,bond,rad)
            
            
            #if 'birad_recombination_F' in par.families or 'all' in par.families:
            #    self.search_birad_recombination_F(natom,atom,bond,rad)
            #if 'birad_recombination_R' in par.families or 'all' in par.families:
            #    self.search_birad_recombination_R(natom,atom,bond,rad)
            #if 'Intra_disproportionation_F' in par.families or 'all' in par.families:
            #    self.search_Intra_disproportionation_F(natom,atom,bond,rad)
            #if 'Intra_disproportionation_R' in par.families or 'all' in par.families:
            #    self.search_Intra_disproportionation_R(natom,atom,bond,rad)
            #if 'r14_birad_scission' in par.families or 'all' in par.families:
            #    self.search_r14_birad_scission(natom,atom,bond,rad)
            #if 'r14_cyclic_birad_scission_R' in par.families or 'all' in par.families:
            #    self.search_r14_cyclic_birad_scission_R(natom,atom,bond,rad)

        
        for name in self.reactions:
            self.reaction_matrix(self.reactions[name], name) 
        
        #verify if every name is unique
        for index in range(len(self.reac_name)-1):
            if self.reac_name[index] in self.reac_name[index+1:]:
                logging.error('Found reaction name "%s" more than once'%self.reac_name[index])
                logging.error('Exiting')
                sys.exit()

        
        #write the reactions that were found to the log
        logging.info('\tFound the following reactions:')
        for rxn in self.reac_name:
            logging.info('\t\t%s'%rxn)

        #sys.exit()
        
        return 0  
            
        
    
    def search_combinatorial(self, natom, atom, bond, rad):
        """ 
        This is a method to create all possible combinations of maximum 3 bond breakings 
        and maximum 3 bond formations.
        
        TODO: allow bond breaking without the atoms forming new bond (only for radicals)

        """
        
        name = 'combinatorial'
        
        if not name in self.reactions:
            self.reactions[name] = []

        instances = generate_product_bond_matrices(self,atom,natom,bond,rad)
        for inst in instances:
            self.reactions[name].append(inst)
        return 0


    def search_intra_H_migration(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        H-R~~~~~~~R* <==> R*~~~~~~~R-H

        Find all unique cases for ring sizes between 3 and 9. Works in both directions.
        """
        
        name = 'intra_H_migration'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        if np.sum(rad) == 0: 
            #find H-migrations over double bonds
            
            if 0:
                #special case of H migrations similar to keto-enol
                motif = ['X' for i in range(4)]
                motif[-1] = 'H'
                instances = start_motif(motif, natom, bond, atom, -1, [[i] for i in range(natom)])
           
                for instance in instances:
                    if bond[instance[0]][instance[1]] > 1:
                        rxns += [instance] 
            else:
                for ringsize in range(3, 9):
                    motif = ['X' for i in range(ringsize)]
                    motif[-1] = 'H'
                    instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)
               
                    for instance in instances:
                        if any([bi > 1 for bi in bond[instance[0]]]):
                            rxns += [instance] 
            
        else:
            instances = []
            for ringsize in range(3, 9):
                motif = ['X' for i in range(ringsize)]
                motif[-1] = 'H'
                instances += start_motif(motif, natom, bond, atom, np.nonzero(rad)[0][0], self.atom_eqv)
            for instance in instances: 
                rxns.append(instance)
        
        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0
        
    
    def search_intra_R_migration(self, natom, atom, bond, rad):
        """ 
        This is an class that covers several RMG classes.
        
        R cannot be an H, this is already taken care of in the intra_H_migration
        
        TODO: merge this with intra H migration families?
        """
        
        if np.sum(rad) != 1: return 

        name = 'intra_R_migration'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        instances = []
        for ringsize in range(3, 9):
            motif = ['X' for i in range(ringsize)]
            instances += start_motif(motif, natom, bond, atom, np.nonzero(rad)[0][0], self.atom_eqv)

        for instance in instances: 
            if not atom[instance[-1]] == 'H':
                rxns.append(instance)
        
        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0

    def search_cpd_H_migration(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        H-C1-C=C-C=C-1 <==> C1=C-C=C-C(-H)-1

        """
        
        if not any([len(ci) == 5 for ci in self.cycle_chain]) : return
        
        name = 'cpd_H_migration'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        bondsum = 0
        
        for cycle in self.cycle_chain:
            if len(cycle) == 5:
                for index, atomi in enumerate(cycle):
                    if index < 4:
                        atomj = cycle[index + 1]
                    else:
                        atomj = cycle[0]
                    if index == 0:
                        atomk = cycle[-1]
                    else:
                        atomk = cycle[index - 1]
                    bondsum += bond[atomi][atomj]
                    if bond[atomi][atomj] == 1 and bond[atomi][atomk] == 1:
                        start = atomi
                        startindex = index
                if bondsum != 7: return # exactly two double bonds
                ring_forw = np.ndarray.tolist(np.roll(cycle, 5 - startindex))
                ring_rev = ring_forw[::-1] # look at the ring in the reverse direction for an H-shift to the other side
                ring_rev = np.ndarray.tolist(np.roll(ring_rev, 1))
                rings = [ring_forw,ring_rev]
                
                Hatomi = -1
                for atomi in range(natom):
                    if atom[atomi] == 'H':
                        if bond[atomi][start] == 1:
                            Hatomi = atomi
                if Hatomi > -1:
                    for ring in rings:
                        instance = ring[:]
                        instance.append(Hatomi)
                        rxns += [instance]

        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)

        return 0
        
    

    def search_intra_OH_migration(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        R*~~~~~~~O-OH <==> HOR~~~~~~~O*

        Find all unique cases for ring sizes between 3 and 9. The H atom is not counted in the cycle size but has to be there.
        """
        
        
        if np.sum(rad) == 0: return
        
        name = 'intra_OH_migration'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(3, 9):
            # forward direction
            motif = ['X' for i in range(ringsize+1)]
            motif[-1] = 'H'
            motif[-2] = 'O'
            motif[-3] = 'O'
            instances = start_motif(motif, natom, bond, atom, np.nonzero(rad)[0][0], self.atom_eqv)
            # reverse direction
            motif = ['X' for i in range(ringsize+1)]
            motif[-1] = 'H'
            motif[-2] = 'O'
            motif[0] = 'O'
            instances += start_motif(motif, natom, bond, atom, np.nonzero(rad)[0][0], self.atom_eqv)
            for ins in instances:
                rxns.append(ins)

        for case in range(len(rxns)):
            rxns[case] = rxns[case][:-1] #cut off H
            
        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)

        return 0



    def search_Intra_RH_Add_Endocyclic_F(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

                                   H
                                   | 
        H-R~~~~~~~R=R ==> R~~~~~~~R-R
                          |         |
                           ---------

        Find all unique cases for ring sizes between 3 and 9. This is for the forward direction.
        """
        
        if np.sum(rad) != 0: return
        if len(self.cycle_chain) > 0: return
        
        name = 'Intra_RH_Add_Endocyclic_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in range(5, 9):
            motif = ['X' for i in range(ringsize + 1)]
            motif[-1] = 'H'
            instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)

            bondpattern = ['X' for i in range(ringsize)]
            bondpattern[0] = 2
            for instance in instances:
                if bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance] 
            

        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-2] == instance[-2] and len(inst) == len(instance):
                    new = 0
            if new:
                self.reactions[name].append(inst)
                
        return 0
        


    def search_Intra_RH_Add_Endocyclic_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

                                   H
                                   | 
        H-R~~~~~~~R=R ==> R~~~~~~~R-R
                          |         |
                           ---------

        Find all unique cases for ring sizes between 3 and 9. This is for the reverse direction.
        """
        
        if len(self.cycle_chain) == 0: return
        if np.sum(rad) != 0: return
        
        name = 'Intra_RH_Add_Endocyclic_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ci in self.cycle_chain:
            motif = ['X' for i in range(len(ci) + 1)]
            motif[-1] = 'H'
            
            instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)
            
            # check if there is a bond between the first and second to last atom
            for instance in instances:
                if bond[instance[0]][instance[-2]] > 0:
                    rxns += [instance[-4:]]
        
        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        self.reactions[name]

        return 0



    def search_Cyclic_Ether_Formation(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        R*~~~~~~~O-OR ==> R~~~~~~~O + OR
                          |_______|

        Find all unique cases for ring sizes between 3 and 9. The OR groups are not counted in the cycle size but have to be there.
        Only the forward direction is included.
        """
        
        
        if np.sum(rad) == 0: return
        
        name = 'Cyclic_Ether_Formation'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(4, 10):
            motif = ['X' for i in range(ringsize)]
            motif[-2] = 'O'
            motif[-3] = 'O'
            motif[0] = 'C'
            rxns += start_motif(motif, natom, bond, atom, np.nonzero(rad)[0][0], self.atom_eqv)

        for instance in range(len(rxns)):
            rxns[instance] = rxns[instance][:-2] #cut off OR
            
        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_Intra_R_Add_Endocyclic_F(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.
        """
        
        if np.sum(rad) == 0: return

        name = 'Intra_R_Add_Endocyclic_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in range(3, 9):
            motif = ['X' for i in range(ringsize)]

            instances = start_motif(motif, natom, bond, atom, np.nonzero(rad)[0][0], self.atom_eqv)
            bondpattern = ['X' for i in range(ringsize-1)]
            bondpattern[-1] = 2
            for instance in instances:
                if bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance]
            
            bondpattern[-1] = 3
            for instance in instances:
                if bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance]



        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
                
        return 0


    def search_Intra_R_Add_ExoTetCyclic_F(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.
        """

        if np.sum(rad) == 0: return

        name = 'Intra_R_Add_ExoTetCyclic_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(3, 9):
            motif = ['X' for i in range(ringsize + 1)]
            rxns += start_motif(motif, natom, bond, atom, np.nonzero(rad)[0][0], self.atom_eqv)

        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst) 

        return 0


    def search_Intra_R_Add_Exocyclic_F(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.
        """
        
        if np.sum(rad) == 0: return

        name = 'Intra_R_Add_Exocyclic_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(3, 9):
            motif = ['X' for i in range(ringsize + 1)]

            instances = start_motif(motif, natom, bond, atom, np.nonzero(rad)[0][0], self.atom_eqv)
            bondpattern = ['X' for i in range(ringsize)]
            bondpattern[-1] = 2
            for instance in instances:
                if bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance]
                    
            bondpattern[-1] = 3
            for instance in instances:
                if bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance]


        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst) 

        return 0



    def search_Intra_RH_Add_Exocyclic_F(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        The general scheme is:

        H-R~~~~~~R=R ==> R~~~~~~R-R-H
                         |      |
                          ------

        The special case of this reaction is Korcel_step1:
              
                R        R   OH 
        R      /          \ /
         \ / \C=O          C
          |       ==>     / \
          O   H          |   O
           \ /          / \ /
            O          R   O

        Implemented as:

                                   --O--O--
                                  |        |
        O=C~~~~~~~~C-O-O-H ==> HO-C~~~~~~~~C
          |                       |
          R                       R

        Find all unique cases for final ring sizes between 3 and 9. The carbonyl dangling R and the
        tail H are included, but are not counted as the ring size, but these two atoms are kept
        because they are needed in the geometry manipulation step.
        """
        
        if len(self.cycle_chain) > 0: return
        
        name = 'Intra_RH_Add_Exocyclic_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(3, 9):
            motif = ['X' for i in range(ringsize+2)]
            motif[-1] = 'H'

            instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)
            bondpattern = ['X' for i in range(ringsize+1)]
            bondpattern[0] = 2
            for instance in instances:
                if bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance]

        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst) 
        
        return 0


    def search_Intra_RH_Add_Exocyclic_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

                                   H
                                   | 
        H-R~~~~~~~R=R ==> R~~~~~~~R-R
                          |         |
                           ---------

        Find all unique cases for ring sizes between 3 and 9. This is for the reverse direction.
        """
        
        name = 'Intra_RH_Add_Exocyclic_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer


        if len(self.cycle_chain) == 0: return
        if np.sum(rad) != 0: return

        for ci in self.cycle_chain:
            motif = ['X' for i in range(len(ci) + 2)]
            motif[-1] = 'H'
            instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)

            # check if there is a bond between the first and second to last atom
            for instance in instances:
                if bond[instance[0]][instance[-3]] > 0:
                    rxns += [instance[-4:]]


        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)

        return 0


    def search_Retro_Ene(self, natom, atom, bond, rad):
        """ 
        This is not an RMG class.
        """
        
        
        if np.sum(rad) != 0: return
        
        name = 'Retro_Ene'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X' for i in range(6)]
        motif[-1] = 'H'
        instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)

        bondpattern = ['X' for i in range(5)]
        bondpattern[0] = 2
        for instance in instances:
            if bondfilter(instance, bond, bondpattern) == 0:
                rxns += [instance] 
        
        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0

    def search_Korcek_step2_judit(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.
              
        Implemented for 4, 5, and 6-membered rings (the R groups are not necessarily identical):

        4-membered ring:

            --O--O--
           |        |
        HO-C--------C-R  ==> RCOOH + R2CO
           |        |
           R        R

        5-membered ring:

            --O--O--
           |        |
        HO-C---C----C-R  ==> RCOOH + R3CC(R)O
           |  / \   |
           R R   R  R

            --O--O--
           |        |
        HO-C---C----C-R  ==> RCOOH + R3CC(R)O
           |  / \   |
           R R   R  R

        6-membered ring:

            ----O---O----
           |             |
        HO-C---C- ---C---C-R  ==> RCOOH + C2R2 + R2CO 
           |  / \   / \  |
           R R   R R   R R

        FIXME: need to generalize to larger structures
        Only the forward direction is included.

        """
        
        
        name = 'Korcek_step2'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(6, 8):
            motif = ['C' for i in range(ringsize)]
            motif[0] = 'O'
            motif[1] = 'O'
            motif[-1] = 'H'
            motif[-2] = 'O'

            for atomi in range(natom):
                if atom[atomi] == 'O':
                    for atomj in range(natom):
                        if atom[atomj] == 'O':
                            if bond[atomi][atomj] == 1:
                                korcek_chain =  start_motif(motif, natom, bond, atom, atomi, self.atom_eqv)
                                for case in range(len(korcek_chain)):
                                    if bond[korcek_chain[case][0]][korcek_chain[case][-3]] == 1:
                                        for ringbond in range(len(korcek_chain[0]) - 2 - 3): # FIXME, assuming just one Korcek hit
                                            marked_chain = korcek_chain[case] + [ringbond + 2] # mark the first atom of the second bond fission
                                            rxns += [marked_chain]

        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst) 
        
        return 0


    def search_Korcek_step2(self, natom, atom, bond, rad):
        """ 
        Generalized Korcek step 
        
        The 4 membered ring equals a 2,2 cycloaddition and is not considered here (no H shift involved)
        
        The 5 membered ring proceeds through a 6 membered transition state (including a 1,2 H migration):

            --O--O--
           |        |
        HO-C---C----C-R  ==> RCOOH + R3CC(R)O
           |  / \   |
           R R   R  R

        6-membered ring: TODO

        Only the forward direction is included.

        """
        
        
        name = 'Korcek_step2'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(5, 6):
            motif = ['X' for i in range(ringsize + 1)]
            #motif[-1] = 'H'
            korcek_chain =  start_motif(motif, natom, bond, atom, -1, self.atom_eqv)
            for ins in korcek_chain:
                if bond[ins[0]][ins[-2]] == 1:
                    rxns += [ins]

        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0

    def search_r22_cycloaddition(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        R      R         R---R
        ||  +  ||  <==   |   |
        R      R         R---R

        N.B.: only the reverse direction is available. Also, the 3 related RMG classes are treated as one.

        """
        
        
        name = 'r22_cycloaddition'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        if not any([len(ci) == 4 for ci in self.cycle_chain]): return 
        
        for ci in self.cycle_chain:
            if len(ci) == 4:
                # there are two ways to slice a 4-mem ring
                ring1 = ci
                ring2 = np.ndarray.tolist(np.roll(ring1, 1))

                # FIXME only works for 1 cycle
                rxns += [ring1]
                rxns += [ring2]

        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0



    def search_r12_cycloaddition(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

                       R--R
        R=R + R: <==   \  /
                        R 

        N.B.: only the reverse direction is available. 

        """
        
        
        name = 'r12_cycloaddition'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        if not any([len(ci) == 3 for ci in self.cycle_chain]): return 
        
        for ci in self.cycle_chain:
            if len(ci) == 3:
                # there are three ways to slice a 3-mem ring
                ring1 = self.cycle_chain
                ring2 = np.ndarray.tolist(np.roll(ring1, 1))
                ring3 = np.ndarray.tolist(np.roll(ring1, 2))

                # FIXME only works for 1 cycle
                rxns += ring1
                rxns += ring2
                rxns += ring3 

        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0



    def search_r12_insertion_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.
        """
        
        #if np.sum(rad) != 0: return
        
        name = 'r12_insertion_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        motif = ['X','X','X']
        instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)
        
        for instance in instances:
            #if all([atom[atomi] != 'H' for atomi in instance]):
            rxns += [instance]

        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1] and inst[2] == instance[2]:
                    new = 0
            if new:
                self.reactions[name].append(inst) 
        
        return 0


    def search_r13_insertion_CO2(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.
        """
        
        #if np.sum(rad) != 0: return
        
        name = 'r13_insertion_CO2'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        motif = ['X','C','O','X']
        instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)
        for instance in instances:
            for atomi in range(natom):
                if not atomi in instance:
                    if atom[atomi] == 'O':
                        if bond[atomi][instance[1]] == 2:
                            rxns += [instance]
                    


        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_r13_insertion_ROR(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.
        """
        
        #if np.sum(rad) != 0: return
        name = 'r13_insertion_ROR'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X','X','X','O']
        rxns = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)
        

        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst) 
        
        return 0


    def search_Diels_alder_addition(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

          R                  R
        //                 /   \
        R       R         R     R
        |  +    ||  <==   ||    |
        R       R         R     R
         \\                \   /
           R                 R

        N.B.: only the reverse direction is available. 

        """
        
        
        name = 'Diels_alder_addition'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        if not any([len(ci) == 6 for ci in self.cycle_chain]): return 

        for ci in self.cycle_chain:
            if len(ci) == 6:
                bondsum = 0
                for index, atomi in enumerate(ci):
                    if index < 5:
                        atomj = ci[index + 1]
                    else:
                        atomj = ci[0]
                    bondsum += bond[atomi][atomj]
                    if bond[atomi][atomj] == 2:
                        start = atomi
                        startindex = index
                if bondsum != 7: return # exactly one double bond
                ring = np.ndarray.tolist(np.roll(ci, 6 - startindex))

                rxns += [ring] # FIXME only works for 1 cycle

        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0



    def search_Intra_Diels_alder_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.
        """
        
        name = 'Intra_Diels_alder_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(3, 9):
            motif = ['X' for i in range(ringsize + 4)]
            instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)

            bondpattern = ['X' for i in range(ringsize + 3)]
            bondpattern[0] = 2
            bondpattern[2] = 2
            bondpattern[-1] = 2

            for instance in instances:
                if bondfilter(instance, bond, bondpattern) == 0:
                    #inst = instance[:4] + instance[-2:]
                    rxns += [instance] 
     

        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
                
        return 0



    def search_ketoenol(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        R=R-O-R <==> R-R-R=O
        """
        
        name = 'ketoenol'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        # enol to keto
        motif = ['C', 'C', 'O', 'X']
        instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)

        # keto to enol
        motif = ['O', 'C', 'C', 'X']
        instances += start_motif(motif, natom, bond, atom, -1, self.atom_eqv)
        bondpattern = [2, 'X', 'X', 'X']
        for instance in instances:
            if bondfilter(instance, bond, bondpattern) == 0:
                rxns += [instance]
            
        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if (inst[0] == instance[0] and inst[1] == instance[1]
                    and inst[2] == instance[2] and inst[3] == instance[3]):
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0
 


    def search_HO2_Elimination_from_PeroxyRadical(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        H-R-R-O-O* ==> R=R + HO2

        N.B.: only the forward direction is available.
        """
        
        
        if np.sum(rad) == 0: return
        
        name = 'HO2_Elimination_from_PeroxyRadical'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['H', 'X', 'X', 'O', 'O']
        rxns += start_motif(motif, natom, bond, atom, -1, self.atom_eqv)
            
        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0
        

    def search_R_Addition_COm3_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        C#O + R* <== R-C*=O

        N.B.: only the reverse direction is available. 
        """
        
        
        if np.sum(rad) == 0: return
        
        name = 'R_Addition_COm3_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X', 'C', 'O']

        instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)

        for instance in instances:
            bondpattern = [1, 2]
            if bondfilter(instance, bond, bondpattern) == 0:
                if rad[instance[1]] == 1:
                    rxns += [instance]
        
        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1] and inst[2] == instance[2]:
                    new = 0
            if new:
                self.reactions[name].append(inst)

        return 0


        
    def search_R_Addition_MultipleBond(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        R=R + R* <== R*-R-R

        N.B.: only the reverse direction is available. 
        """
        
        
        if np.sum(rad) == 0: return
        
        name = 'R_Addition_MultipleBond'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X', 'X', 'X']
        rxns += start_motif(motif, natom, bond, atom, np.nonzero(rad)[0][0], self.atom_eqv)
        
        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1] and inst[2] == instance[2]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_12_shift_S_F(self,natom, atom, bond, rad):
        """
        This is an RMG class.
        """

        if np.sum(rad) != 1: return
        
        name = '12_shift_S_F'
        
        if not name in self.reactions:
            self.reactions[name] = []
        
        motif = ['X','S','X']
        rxns = start_motif(motif, natom, bond, atom, np.nonzero(rad)[0][0], self.atom_eqv)

        #filter for identical reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1] and inst[2] == instance[2]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        return 0


    def search_12_shift_S_R(self,natom, atom, bond, rad):
        """
        This is an RMG class.
        """
        
        if np.sum(rad) != 1: return
        
        name = '12_shift_S_R'
        
        if not name in self.reactions:
            self.reactions[name] = []
        
        motif = ['S','X','X']
        rxns = start_motif(motif, natom, bond, atom, np.nonzero(rad)[0][0], self.atom_eqv)
        
        #filter for identical reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1] and inst[2] == instance[2]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        return 0


    def search_r13_insertion_RSR(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.
        """
        
        #if np.sum(rad) != 0: return
        name = 'r13_insertion_RSR'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X','X','X','S']
        rxns = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)
        

        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst) 
        
        return 0


    def search_R_Addition_CSm_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        C#S + R* <== R-C*=S

        N.B.: only the reverse direction is available. 
        """
        
        
        if np.sum(rad) == 0: return
        
        name = 'R_Addition_CSm_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X', 'C', 'S']

        instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)

        for instance in instances:
            bondpattern = [1, 2]
            if bondfilter(instance, bond, bondpattern) == 0:
                if rad[instance[1]] == 1:
                    rxns += [instance]
        
        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1] and inst[2] == instance[2]:
                    new = 0
            if new:
                self.reactions[name].append(inst)

        return 0


    def search_r14_birad_scission(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.
        """

        if np.sum(rad) != 2: return
        
        
        name = 'r14_birad_scission'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        
        motif = ['X','X','X','X']
        instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)
        for instance in instances: 
            if rad[instance[0]] == 1 and rad[instance[-1]] == 1:
                rxns += [instance]

        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[1] == instance[1] and inst[2] == instance[2]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_r14_cyclic_birad_scission_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.
        """

        if np.sum(rad) != 0: return
        
        name = 'r14_cyclic_birad_scission_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in range(5, 9):
            motif = ['X' for i in range(ringsize)]
            instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)
            
            bondpattern = ['X' for i in range(ringsize - 1)]
            bondpattern[0] = 2
            bondpattern[-1] = 2
            
            for instance in instances:
                if bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance] 
                    


        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_birad_recombination_F(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.
        """

        if np.sum(rad) != 2: return
        
        name = 'birad_recombination_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in range(3, 9):
            motif = ['X' for i in range(ringsize)]
            instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)
           
            for instance in instances: 
                if rad[instance[0]] == 1 and rad[instance[-1]] == 1:
                    rxns += [instance]


        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_birad_recombination_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.
        """

        if np.sum(rad) != 0: return
        if len(self.cycle_chain) == 0: return
        
        name = 'birad_recombination_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X','X']
        instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)
        for instance in instances: 
            if instance[0] in self.cycle and instance[1] in self.cycle :
                rxns += [instance]


        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_Intra_disproportionation_F(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.
        """

        if np.sum(rad) != 2: return
        
        name = 'Intra_disproportionation_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in range(5, 9):
            motif = ['X' for i in range(ringsize)]
            motif[-1] = 'H'
            instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)
           
            for instance in instances: 
                if rad[instance[0]] == 1 and rad[instance[-3]] == 1:
                    rxns += [instance]


        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_Intra_disproportionation_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.
        """

        if np.sum(rad) != 0: return
        
        name = 'Intra_disproportionation_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(5, 9):
            motif = ['X' for i in range(ringsize)]
            motif[-1] = 'H'
            
            instances = start_motif(motif, natom, bond, atom, -1, self.atom_eqv)
            
            bondpattern = ['X' for i in range(ringsize - 1)]
            bondpattern[0] = 2
            
            for instance in instances:
                if bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance] 


        #filter for the same reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def reaction_matrix(self, reac_list, reac_id):
        """ 
        Create arrays to store all reactions for species.
        input: 
        reac_list: atom motifs from individual searches
        reac_id: reaction name (e.g., HO2_Elimination_from_PeroxyRadical) from individual searc functions
        Every reaction type just makes the below arrays longer, generated as reactions are found.

        generated:
        reac_type: reaction class identifier
        reac_inst: reaction instance defined by the important atoms
        reac_step: the step at which the search is at
        reac_scan_energy: for each reaction the energy as a function of steps, only used for scanning type searches, e.g. R_Addition_MultipleBond
        rec_ts_done: the last calculations is submitted in the sequence
        reac_ts_geom: the geometry of the TS
        reac_ts_freq: the freqencies of the TS
        reac_name: the base name of the file to run - created for each reaction later
        """
        
        
        self.reac_type += [reac_id for i in range(len(reac_list))]
        self.reac_inst += reac_list
        self.reac_step += [0 for i in range(len(reac_list))]
        self.reac_scan_energy += [[] for i in range(len(reac_list))]
        self.reac_ts_done += [0 for i in range(len(reac_list))] 
        self.reac_ts_geom += [0 for i in range(len(reac_list))]
        self.reac_ts_freq += [0 for i in range(len(reac_list))]
        
        if reac_id == 'intra_H_migration':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'intra_R_migration':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'intra_OH_migration':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'cpd_H_migration':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) + '_' + str(reac_list[i][-2] + 1) for i in range(len(reac_list))]
        elif reac_id == 'Intra_RH_Add_Endocyclic_F':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(len(reac_list[i])) + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-2] + 1) for i in range(len(reac_list))]
        elif reac_id == 'Intra_RH_Add_Endocyclic_R':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'Cyclic_Ether_Formation':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'Intra_RH_Add_Exocyclic_F':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'Intra_RH_Add_Exocyclic_R':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'Retro_Ene':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'Intra_R_Add_Endocyclic_F':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'Intra_R_Add_ExoTetCyclic_F':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-2] + 1) + '_' + str(reac_list[i][-1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'Intra_R_Add_Exocyclic_F':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-2] + 1) for i in range(len(reac_list))]
        elif reac_id == 'Korcek_step2':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'r22_cycloaddition':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'r12_cycloaddition':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'r12_insertion_R':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1) for i in range(len(reac_list))]
        elif reac_id == 'r13_insertion_CO2':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'r13_insertion_ROR':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1) + '_' + str(reac_list[i][3] + 1) for i in range(len(reac_list))]
        elif reac_id == 'r14_birad_scission':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1) for i in range(len(reac_list))]
        elif reac_id == 'r14_cyclic_birad_scission_R':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'birad_recombination_F':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'birad_recombination_R':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'Intra_disproportionation_F':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'Intra_disproportionation_R':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'Diels_alder_addition':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'Intra_Diels_alder_R':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'ketoenol':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1)  + '_' + str(reac_list[i][2] + 1)  + '_' + str(reac_list[i][3] + 1) for i in range(len(reac_list))]
        elif reac_id == 'HO2_Elimination_from_PeroxyRadical':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) for i in range(len(reac_list))]
        elif reac_id == 'R_Addition_COm3_R':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1) for i in range(len(reac_list))]
        elif reac_id == 'R_Addition_MultipleBond':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1) for i in range(len(reac_list))]
        elif reac_id == '12_shift_S_F':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1) for i in range(len(reac_list))]
        elif reac_id == '12_shift_S_R':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1) for i in range(len(reac_list))]
        elif reac_id == 'R_Addition_CSm_R':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1) for i in range(len(reac_list))]
        elif reac_id == 'r13_insertion_RSR':
            self.reac_name += [str(self.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1) + '_' + str(reac_list[i][3] + 1) for i in range(len(reac_list))]
        else:
            self.reac_name += [0 for i in range(len(reac_list))] 
        
        return 0


def main():
    """
    This is the main object of the code, the stationary point
    """



if __name__ == "__main__":
    main()
