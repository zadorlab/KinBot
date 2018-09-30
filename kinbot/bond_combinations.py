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
import os, sys, copy
import itertools

sys.dont_write_bytecode = True

from stationary_pt import *
from cheminfo import *

smi_list = []

def equivalent_bond(bonds,b1,well):
    """
    Method checks if the current bond is equivalent to one of the bonds in the bonds list
    """
    for b2 in bonds:
        # check if the first atoms of both bonds are in the same equivalency list
        # and if the second atoms of both bonds are in the same equivalency list
        b1_e1 = -1 # first bond first atom equivalence index
        b1_e2 = -1 # first bond second atom equivalence index
        b2_e1 = -1 # second bond first atom equivalence index
        b2_e2 = -1 # second bond second atom equivalence index
        for i,eq in enumerate(well.atom_eqv):
            if b1[0] in eq:
                b1_e1 = i
            if b1[1] in eq:
                b1_e2 = i
            if b2[0] in eq:
                b2_e1 = i
            if b2[1] in eq:
                b2_e2 = i
        if b1_e1 == b2_e1 and b1_e2 == b2_e2:
            return 1
    return 0

def generate_product_bond_matrices(mol,atom,natom,bond,rad):
    """
    This method does the following:
    1. Generate all possible combinations of three bonds in the molecule
    2. For each combination, create all the possible permutation of the 6 atoms involved
    3. Filter the combinations that lead to identical atom rearrangements
    4. Create (and return) a list of bond matrices of the product and the transition state
    """

    # assume up to three bonds can break
    nbonds = 3

    bonds = []
    for i in range(len(atom)-1):
        for j in range(i+1,len(atom)):
            if bond[i,j] > 0:
                if not equivalent_bond(bonds,[i,j],mol):
                    bonds.append([i,j])

    #contains all possibilities of combining three different bonds (six atoms)
    #from which the permutations are made to define reactions
    reactive_atoms = [] 
    for i in range(len(bonds)-2):
        for j in range(i+1,len(bonds)-1):
            for k in range(j+1,len(bonds)):
                # an atom can be in maximum two of the three bonds 
                max = 0
                comb = [bonds[i],bonds[j],bonds[k]]
                for at in range(len(atom)):
                    value = sum([(at in c) for c in comb])
                    if value > max:
                        max = value
                react = bonds[i] + bonds[j] + bonds[k]
                if max < 3 :
                    reactive_atoms.append(react)


    #contains product bonds
    all_prods = []
    #contains react bonds
    all_reacts = []

    for comb in reactive_atoms:
        # instead of directly adding the reactions to the all_prods and all_reacts list, 
        # first add them to this rxns list. Reactions from one comb are compared to one 
        # another to assure uniqueness. Reactions from different comb lists are always different 
        # form one another. 
        # this significantly speeds up the generation of the reactions 
        rxns = []
        #create all the permutation for the current set of three bonds
        perms = list(itertools.permutations(comb,6)) 
        
        r = [] # put the reactant in a set of three bonds
        for i in range(nbonds):
            b = sorted([comb[2*i], comb[2*i+1]])
            r.append(b)
        
        for perm in perms: #iterate the permutations
            p = [] # create a list of three bonds from the permutation
            for i in range(nbonds):
                b = sorted([perm[2*i], perm[2*i+1]])
                p.append(b)
            new = 1 # is this reaction new?
            for re in rxns: # check if this list of bonds has been created before
                if sorted(re) == sorted(p):
                    new = 0
            if sorted(r) == sorted(p): #check if this is the identity permutation
                new = 0
            if new: # add the reaction to the reactions of the current bond set if it is new
                rxns.append(p)
        for rxn in rxns: # add the reaction to the list of all reactions
            all_prods.append(rxn)
            
            all_reacts.append(r)
    
    
    reactions = []

    # generate the ts and product bond matrices
    for i,prod in enumerate(all_prods):
        reac = all_reacts[i]
        reactions.append([prod,reac])
        
        #initialize the bond matrices to the one of the reactant
        prod_bond = copy.deepcopy(bond) 
        ts_bond = [[float(bond[i][j]) for j in range(len(atom))] for i in range(len(atom))]
        
        for i in range(len(atom)-1):
            for j in range(i+1,len(atom)):
                if [i,j] in reac:
                    if not [i,j] in prod:
                        # a bond in the reactant and not in the product is a bond that is 
                        # either broken or decreases in order
                        prod_bond[i][j] -= 1
                        ts_bond[i][j] = float(ts_bond[i][j]) - .5
                elif [i,j] in prod:
                    # a bond in the product and not in the reactant is a bond that is 
                    # either formed or increases in order
                    prod_bond[i][j] += 1
                    ts_bond[i][j] = float(ts_bond[i][j]) + .5
        #reactions.append([prod_bond,ts_bond])
        rdmol, smi,struc = create_rdkit_mol(prod_bond,atom)
        smi_list.append(smi)
    
    return reactions

    
    
    

def main():
    mol = stationary_pt('well0')
    smi = 'O=CCCOO'
    mult = 1
    charge = 0

    obmol, structure = generate_3d_structure(smi)
    natom = len(obmol.atoms)
    structure = np.reshape(structure, (natom,4))

    atom = structure[:,0]
    
    mol.geom = structure[:,1:4].astype(float)
    mol.characterize(natom,atom,mult,charge)
    
    #make and rdkit molecule from the mol.bond matrix
    rdmol, smi, struc = create_rdkit_mol(mol.bond,atom)
    smi_list.append(smi)
    
    generate_product_bond_matrices(mol,atom,natom,mol.bond,[])

if __name__ == "__main__":
    main()



    
