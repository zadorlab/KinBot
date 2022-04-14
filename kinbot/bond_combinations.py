import itertools
import numpy as np

from kinbot import find_motif
from kinbot.stationary_pt import StationaryPoint


def equivalent_bond(b1, b2, species):
    """
    Method checks if the current bond is equivalent to
    one of the bonds in the bonds list.
    Bonds are always defined by the smallest atom index
    first and the largest second.
    """
    # check if the first atoms of both bonds are in the
    # same equivalency list and if the second atoms of
    # both bonds are in the same equivalency list
    b1_e1 = -1  # first bond first atom equivalence index
    b1_e2 = -1  # first bond second atom equivalence index
    b2_e1 = -1  # second bond first atom equivalence index
    b2_e2 = -1  # second bond second atom equivalence index
    for i, eq in enumerate(species.atom_eqv):
        if b1[0] in eq:
            b1_e1 = i
        if b1[1] in eq:
            b1_e2 = i
        if b2[0] in eq:
            b2_e1 = i
        if b2[1] in eq:
            b2_e2 = i
    return b1_e1 == b2_e1 and b1_e2 == b2_e2


def generate_all_product_bond_matrices(mol, par):
    """
    Generate all product bond matrices with the maximum number of bonds
    being 2 or 3
    """
    reactions = []
    # generate the reactions for closed shell species
    # and for the non-radical part of open shell species
    if par['comb_molec']:
        nbonds_list = range(par['min_bond_break'], par['max_bond_break'] + 1)
        for bond in mol.bonds:
            for nbonds in nbonds_list:
                rxns = generate_product_bond_matrices(mol, bond, nbonds, par, rad=-1)
                reactions.extend(rxns)
    # generate the reactions in which radicals participate
    nbonds_list = range(par['min_bond_break'] - 1, par['max_bond_break'])
    for i, bond in enumerate(mol.bonds):
        rads = np.nonzero(mol.rads[i])[0]
        for nbonds in nbonds_list:
            # reactions involving a radical atom
            # a bond is formed with that atom, no bond is broken
            # and another atom only has bond breaking, no forming
            if par['comb_rad']:
                for rad in rads:
                    rxns = generate_product_bond_matrices(mol, bond,
                                                          nbonds, par, rad=rad)
                    reactions.extend(rxns)
            # reactions involving a pi electron leading to a new lone pair
            if par['comb_pi']:
                for i in range(mol.natom):
                    if any(bij == 2 for bij in mol.bond[i]):
                        rxns = generate_product_bond_matrices(mol, bond,
                                                              nbonds, par, rad=i)
                        reactions.extend(rxns)
            # TODO: reactions with lone electron pairs
            if par['comb_lone']:
                for i, ai in enumerate(mol.atom):
                    if ai == 'O':
                        rxns = generate_product_bond_matrices(mol, bond,
                                                              nbonds, par, rad=i)
                        reactions.extend(rxns)
    return reactions


def get_product_bonds(bonds, par, rad=-1):
    """
    This method creates a list of new atom pairs
    which are all different than the atom pairs in
    the bonds provided as argument.
    """
    # list of all the atoms in the bonds list
    atoms = []
    for bi in bonds:
        atoms.extend(bi)
    if rad > -1:
        for bond in bonds:
            if rad in bond:
                return []
        atoms.append(rad)
    # generate all the possible (new) atom pairs
    pairs = []
    for i, at1 in enumerate(atoms[:-1]):
        for at2 in atoms[i+1:]:
            if not sorted([at1, at2]) in bonds:
                pairs.append([at1, at2])
    # generate all the lists of pairs which are
    # of the same length as the bonds list
    all_prods = itertools.combinations(pairs, len(bonds))
    # final list of all products
    prods = []
    for prod in all_prods:
        prod_atoms = []
        for pi in prod:
            prod_atoms.extend(pi)
        if rad == -1:
            # to verify if a list of product bonds is meaningful
            # the list of product atoms needs to be identical to
            # the list of reactant atoms
            if sorted(atoms) == sorted(prod_atoms):
                prods.append(sorted(prod))
                # also make a list in which one of the bonds is not formed,
                # this is needed to break the valence of atoms and form
                # for example zwitterionic species.
                if par['break_valence']:
                    for i in range(len(prod)):
                        prods.append(prod[:i]+prod[i+1:])
        else:
            # to verify if a list of product bonds is meaningful
            # the list of product atoms should be the list of
            # reactant atoms minus one atoms, which cannot be the
            # radical atom
            count = 0  # number of atoms for which there is a discrepancy
            atom = []  # list of atoms for which there is a discrepancy
            for at in atoms:
                if atoms.count(at) != prod_atoms.count(at):
                    count += np.abs(atoms.count(at) - prod_atoms.count(at))
                    atom.append(at)
            if count == 1:
                if atom[0] != rad:
                    prods.append(sorted(prod))
    return prods


def is_identical(prod1, prod2, mol):
    """
    Check if two products are unique by verifying if the
    shorted chain between the new bonds is equivalent
    Idential chains are chains of the same length of which
    all elements are pairwise equivalent.
    """
    # Check if all bonds of one product are equivalent to a bond
    # of the other product and put this in a dictionary
    eq = {}
    for i, pi in enumerate(prod1):
        if i not in eq.values():
            for j, pj in enumerate(prod2):
                if j not in eq:
                    if equivalent_bond(pi, pj, mol):
                        eq[j] = i
    if len(eq) == len(prod1):
        # this means all bonds have one equivalent bond in the
        # other bonds list

        # now, create chains between the two members of each new bond
        # and check if the full chain is of identical length
        for i in eq:
            b1 = prod1[eq[i]]
            b2 = prod2[i]
            chain1 = get_chain(b1[0], b1[1], mol)
            chain2 = get_chain(b2[0], b2[1], mol)
            if len(chain1) != len(chain2):
                return 0
        return 1
    else:
        return 0


def get_chain(a1, a2, mol):
    """
    Get the shortest chain between two atoms
    """
    for i in range(1, mol.natom):
        motif = ['X' for j in range(i)]
        instances = find_motif.start_motif(motif,
                                           mol.natom,
                                           mol.bond,
                                           mol.atom,
                                           a1,
                                           [[k] for k in range(mol.natom)])
        for ins in instances:
            if ins[-1] == a2:
                return ins
    return []


def generate_ts(reac, prod, bond):
    """
    Method to generate the bond matrix of the transition state
    This is needed to know whether to break (form) a bond instead
    of decreasing (increasing) its bond order
    """
    ts_bond = [[float(bij) for bij in bi] for bi in bond]
    if reac[0]:
        for ri in reac:
            i = ri[0]
            j = ri[1]
            ts_bond[i][j] -= 0.5
    if prod[0]:
        for pi in prod:
            i = pi[0]
            j = pi[1]
            ts_bond[i][j] += 0.5
    return ts_bond


def generate_product_bond_matrices(mol, bond, nbonds, par, rad=-1):
    """
    This method does the following:
    1. Generate all possible combinations of three bonds in the molecule
    2. For each combination, create all the possible
       permutation of the 6 atoms involved
    3. Filter the combinations that lead to identical atom rearrangements
    4. Create (and return) a list of bond matrices of the product
       and the transition state
    rad: index of the radical site to consider, -1 if this search is applied
    to closed shell (or the non-radical part of a) molecule.
    """
    # list of all reactions to explore
    reactions = []

    bonds = []
    for i in range(len(mol.atom)-1):
        for j in range(i+1, len(mol.atom)):
            if bond[i, j] > 0:
                bonds.append([i, j])

    reactive_atoms = itertools.combinations(bonds, nbonds)

    # contains product bonds
    all_prods = []
    # contains react bonds
    all_reacts = []

    for comb in reactive_atoms:
        # instead of directly adding the reactions to the all_prods
        # and all_reacts list, first add them to this rxns list.
        # Reactions from one comb are compared to one another to
        # assure uniqueness. Reactions from different comb lists
        # are always different form one another.
        # This significantly speeds up the generation of the reactions
        rxns = []

        # look for all possibilies in which all of the nbonds break
        prod_bonds = get_product_bonds(comb, par, rad=rad)

        for prod in prod_bonds:
            # verify if this reaction is unique
            new = 1
            for rxn in reactions:
                if (is_identical(rxn[1], prod, mol) and
                        is_identical(rxn[0], comb, mol)):
                    new = 0
            if new:
                ts = generate_ts(comb, prod, bond)
                # add the reaction three times, with an early, a mid
                # and a late ts
                # reactions.append([comb, prod, ts, 0])
                reactions.append([comb, prod, ts, 1])
                # reactions.append([comb, prod, ts, 2])
    return reactions


def main():
    smi = '[CH2]CC'
    mult = 1
    charge = 0
    mol = StationaryPoint('well0', charge, mult, smiles=smi)
    mol.characterize()
    reactions = generate_all_product_bond_matrices(mol)


if __name__ == "__main__":
    main()
