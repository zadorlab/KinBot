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
import os,sys
import logging
import numpy as np

from geom import *
from stationary_pt import *
import par


"""

PROBLEMS with external symmetry; 

** Consecutive double bonds
** C-C#C-C
** R1-C(-R1) = C(-R2)-R2

special cases:
** ethane: external vs internal
** ethylene


"""

def calculate_symmetry(species, natom, atom):
    """
    Calculate the symmetry numbers (external and internal) and 
    the number of optical isomers of a molecule based on some
    heuristics for atom-centered symmetry, bond-centered 
    symmetry and ring-centered symmetry. 
    
    TODO: 
    * Symmetry along consecutive double bonds is not well perceived
    
    """
    sigma_ext = 1
    nopt = 1
    sigma_int = [[1 for i in range(natom)] for i in range(natom)]
    
    #get the number of optical isomers
    for at in range(natom):
        nei = get_neighbors(species,at)
        if len(nei) == 4:
            nei_atomid = sorted([species.atomid[ni] for ni in nei])
            if len(nei_atomid) == len(set(nei_atomid)):
                nopt = 2

    lin = start_linear(species,natom)
    
    #get all atom-centered symmetries
    for at in range(natom):
        if species.cycle[at] == 0:
            nei = get_neighbors(species,at)
            if len(nei) == 1: #no symmetry contributions
                continue
            elif len(nei) == 2:
                if species.atomid[nei[0]] == species.atomid[nei[1]]:
                    linear = 0
                    for li in lin:
                        if at in li[1:-1]:
                            linear = 1
                    if not linear:
                        sigma_ext *= 2
            elif len(nei) == 3:
                if (species.atomid[nei[0]] == species.atomid[nei[1]] and 
                species.atomid[nei[1]] == species.atomid[nei[2]]):
                    sigma_ext *= 6
                """elif species.atomid[nei[0]] == species.atomid[nei[1]]:
                    if species.bond[at][nei[2]] == 1:
                        sigma_int[at] *= 2
                elif species.atomid[nei[1]] == species.atomid[nei[2]]:
                    if species.bond[at][nei[0]] == 1:
                        sigma_int[at] *= 2
                elif species.atomid[nei[0]] == species.atomid[nei[2]]:
                    if species.bond[at][nei[1]] == 1:
                        sigma_int[at] *= 2"""
            elif len(nei) == 4:
                nei_atomid = sorted([species.atomid[ni] for ni in nei])
                if all([ati == nei_atomid[0] for ati in nei_atomid]):
                    sigma_ext *= 12
                elif any([nei_atomid.count(ati) == 3 for ati in nei_atomid]):
                    #sigma_int[at] *= 3
                    continue
                elif all([nei_atomid.count(ati) == 2 for ati in nei_atomid]):
                    sigma_ext *= 2
                elif any([nei_atomid.count(ati) == 2 for ati in nei_atomid]):
                    #not symmetric nor optically active
                    continue
                else:
                    continue
    
    #internal rotational symmetry
    #both atoms can be in a cycle, but not in the same cycle
    for li in lin:
        i = li[0]
        cycle = []
        for cyc in species.cycle_chain:
            if i in cyc:
                cycle.extend(cyc)
        j = li[-1]
        if not j in cycle:
            if any([species.bond[li[ai]][li[ai+1]] == 1 for ai in range(len(li)-1)]):
                nei1 = get_neighbors(species,i)
                nei1 = [species.atomid[ni] for ni in nei1 if ni != li[1]]
                nei2 = get_neighbors(species,j)
                nei2 = [species.atomid[ni] for ni in nei2 if ni != li[-2]]
                if len(nei1) > 0 and len(nei2) > 0:
                    s1 = 1
                    if all([ati == nei1[0] for ati in nei1]):
                        s1 = len(nei1)
                    s2 = 1
                    if all([ati == nei2[0] for ati in nei2]):
                        s2 = len(nei2)
                    sigma_int[i][j] = lcm(s1,s2)

    #get all bond-centered symmetries
    for li in lin:
        i = li[0]
        j = li[-1]
        if species.atomid[i] == species.atomid[j] and species.cycle[i] == 0:
            sigma_ext *= 2

        if species.cycle[i] == 0 and species.cycle[j] == 0:
            nei1 = get_neighbors(species,i)
            nei1 = [species.atomid[ni] for ni in nei1 if ni != li[1]]
            nei2 = get_neighbors(species,j)
            nei2 = [species.atomid[ni] for ni in nei2 if ni != li[-2]]
            
            if len(nei1) == 0 and len(nei2) > 1 and all([ati == nei2[0] for ati in nei2]):
                #if all the neighbors of j are the same, it should have been taken into account earlier
                nei = get_neighbors(species,j)
                nei = [species.atomid[ni] for ni in nei]
                if not all([ati == nei[0] for ati in nei]):
                    sigma_ext *= len(nei2)
            elif len(nei2) == 0 and len(nei1) > 1 and all([ati == nei1[0] for ati in nei1]):
                #if all the neighbors of i are the same, it should have been taken into account earlier
                nei = get_neighbors(species,i)
                nei = [species.atomid[ni] for ni in nei]
                if not all([ati == nei[0] for ati in nei]):
                    sigma_ext *= len(nei1)
            elif all([ati == nei1[0] for ati in nei1]) and all([ati == nei2[0] for ati in nei2]):
                if len(nei1) == len(nei2) and len(nei1) > 0:
                    sigma_ext *= len(nei1)

    #get all ring-centered symmetries
    for cyc in species.cycle_chain:
        cyc_atomid = [species.atomid[ci] for ci in cyc]
        symm = 0
        for i in range(len(cyc)):
            new_order = np.roll(np.array(cyc_atomid),-i)
            if all([cyc_atomid[at] == new_order[at] for at in range(len(cyc))]):
                symm += 1
            new_order_reversed = np.roll(np.array(cyc_atomid[::-1]),-i)
            if all([cyc_atomid[at] == new_order_reversed[at] for at in range(len(cyc))]):
                symm += 1
        #additional patch: if an atom has two identical neighbors 
        #along the ring, but two distinct neighbors outside the ring
        #the symmetry number needs to be divided by 2
        if symm > 1:
            divide = 1
            for at in cyc:
                nei = get_neighbors(species,at)
                cyc_nei = [ni for ni in nei if ni in cyc]
                other_nei = [ni for ni in nei if ni not in cyc]
                if len(other_nei) > 1:
                    if species.atomid[cyc_nei[0]] == species.atomid[cyc_nei[1]]:
                        if not all([species.atomid[other_nei[0]] == species.atomid[oi] for oi in other_nei]):
                            divide = 2
            symm /= divide
        if symm > 0:
            sigma_ext *= symm

    species.sigma_ext = sigma_ext
    species.sigma_int = sigma_int
    species.nopt = nopt

def start_linear(species,natom):
    """
    Get all the 'neighbors' of i which are connected to i via
    one or a set of linear bonds. 
    
    This includes 
    * all direct bonds if they are not part of a linear system
    * all sets of consecutive double bonds
    * all sets of alternating triple and single bonds
    
    Only look at indices larger than i, as the other bonds have been considered with lower 
    indices
    """
    lin = []
    
    for i in range(natom - 1):
        for j in range(i+1,natom):
            if species.bond[i][j] > 0:
                if len(get_neighbors(species,j)) == 2:
                    k = [ni for ni in get_neighbors(species,j) if ni != i][0]
                    if calc_angle(species.geom[i],species.geom[j],species.geom[k]) > np.pi * 175. / 180.:
                        new_lin = get_linear(species,[i,j,k],natom)
                        new = 1
                        for li in lin:
                            if all([ni in li for ni in new_lin]):
                                new = 0
                        if new:
                            lin.append(new_lin)
    
    for i in range(natom - 1):
        for j in range(i+1,natom):
            if species.bond[i][j] > 0:
                new = 1
                for li in lin:
                    if i in li and j in li:
                        new = 0
                if new:
                    lin.append([i,j])
    return lin

def get_linear(species,visited,natom):
    for j in range(natom):
        if j not in visited:
            if species.bond[visited[-1]][j] > 0:
                if calc_angle(species.geom[visited[-2]],species.geom[visited[-1]],species.geom[j]) > np.pi * 175. / 180.:
                    visited.append(j)
                    return get_linear(species,visited,natom)
                else:
                    return visited
    return visited
    
def lcm(x, y):
   if x > y:
       z = x
   else:
       z = y

   while(True):
       if((z % x == 0) and (z % y == 0)):
           lcm = z
           break
       z += 1

   return lcm

def get_neighbors(species,at):
    """
    Get the neighbors of atom at in species
    """
    return [i for i,b in enumerate(species.bond[at]) if b > 0]
    
def main():
    """
    Calculate the total rotational symmetry number 
    Calculate the number of optical isomers
    """

if __name__ == "__main__":
    main()

