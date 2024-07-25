import os
import logging
import numpy as np
import copy

from kinbot import geometry

logger = logging.getLogger('KinBot')


def make_zmat_from_cart(species, rotor, cart, mode):
    """
    Rearrange geometry defined in Cartesian into a Z-matrix,
    with references suitable for a 1-D hindered rotor scan.
    If mode = 0: all rotatable bonds
    If mode = 1: only those bonds, which generate conformers
    If mode = 2: suply your one rotor in rotor as a list of atom indices
    """
    natom = species.natom
    atom = species.atom
    if mode == 0:
        a = species.dihed[rotor][0]
        b = species.dihed[rotor][1]
        c = species.dihed[rotor][2]
        d = species.dihed[rotor][3]
    elif mode == 1:
        a = species.conf_dihed[rotor][0]
        b = species.conf_dihed[rotor][1]
        c = species.conf_dihed[rotor][2]
        d = species.conf_dihed[rotor][3]
    elif mode == 2:
        a = rotor[0]
        b = rotor[1]
        c = rotor[2]
        d = rotor[3]

    groupA = np.zeros(natom, dtype=int)
    groupB = np.zeros(natom, dtype=int)
    groupC = np.zeros(natom, dtype=int)
    groupD = np.zeros(natom, dtype=int)

    groupA[a] = 1
    groupB[b] = 1
    groupC[c] = 1
    groupD[d] = 1

    # get all immediate neighbors of a A, B, C, and D
    for i in range(natom):
        if species.bond[a][i] > 0 and i != b:
            groupA[i] = 1
        elif species.bond[b][i] > 0 and i != a and i != c:
            groupB[i] = 1
        elif species.bond[c][i] > 0 and i != b and i != d:
            groupC[i] = 1
        elif species.bond[d][i] > 0 and i != c:
            groupD[i] = 1

    found = 1
    while found > 0:
        found = 0
        for i in range(natom):
            if groupA[i] == 1:
                for j in range(natom):
                    if species.bond[i][j] > 0 and j != b and groupA[j] != 1:
                        groupA[j] = 1
                        found = 1
            elif groupB[i] == 1:
                for j in range(natom):
                    if species.bond[i][j] > 0 and j != a and j != c and groupB[j] != 1:
                        groupB[j] = 1
                        found = 1
            elif groupC[i] == 1:
                for j in range(natom):
                    if species.bond[i][j] > 0 and j != b and j != d and groupC[j] != 1:
                        groupC[j] = 1
                        found = 1
            elif groupD[i] == 1:
                for j in range(natom):
                    if species.bond[i][j] > 0 and j != c and groupD[j] != 1:
                        groupD[j] = 1
                        found = 1

    zmat_atom = ['X' for i in range(natom)]
    zmat_ref = np.zeros((natom, 3), dtype=int) - 1
    zmat = np.zeros((natom, 3)) - 1

    # FIXME need to take care about TS structures maybe
    zmatorder = [-1 for i in range(natom)]

    zmat_atom[0] = atom[a]
    zmatorder[0] = a

    zmat_atom[1] = atom[b]
    zmatorder[1] = b
    zmat_ref[1][0] = 1
    zmat[1][0] = np.linalg.norm(cart[b] - cart[a])

    zmat_atom[2] = atom[c]
    zmatorder[2] = c
    zmat_ref[2][0] = 2
    zmat[2][0] = np.linalg.norm(cart[c] - cart[b])
    zmat_ref[2][1] = 1
    zmat[2][1] = np.degrees(geometry.calc_angle(cart[c], cart[b], cart[a]))

    zmat_atom[3] = atom[d]
    zmatorder[3] = d
    zmat_ref[3][0] = 3
    zmat[3][0] = np.linalg.norm(cart[d] - cart[c])
    zmat_ref[3][1] = 2
    zmat[3][1] = np.degrees(geometry.calc_angle(cart[d], cart[c], cart[b]))
    zmat_ref[3][2] = 1
    zmat[3][2], collin = geometry.calc_dihedral(cart[d], cart[c], cart[b], cart[a])

    j = 4
    for i in range(natom):
        if i == a or i == b or i == c or i == d:
            continue
        zmat_atom[j] = atom[i]
        zmatorder[j] = i
        if groupA[i] == 1:
            zmat_ref[j][0] = 1
            zmat[j][0] = np.linalg.norm(cart[i] - cart[a])
            zmat_ref[j][1] = 2
            zmat[j][1] = np.degrees(geometry.calc_angle(cart[i], cart[a], cart[b]))
            zmat_ref[j][2] = 3
            zmat[j][2], collin = geometry.calc_dihedral(cart[i], cart[a], cart[b], cart[c])
        elif groupB[i] == 1:
            zmat_ref[j][0] = 2
            zmat[j][0] = np.linalg.norm(cart[i] - cart[b])
            zmat_ref[j][1] = 3
            zmat[j][1] = np.degrees(geometry.calc_angle(cart[i], cart[b], cart[c]))
            zmat_ref[j][2] = 4
            zmat[j][2], collin = geometry.calc_dihedral(cart[i], cart[b], cart[c], cart[d])
        elif groupC[i] == 1:
            zmat_ref[j][0] = 3
            zmat[j][0] = np.linalg.norm(cart[i] - cart[c])
            zmat_ref[j][1] = 2
            zmat[j][1] = np.degrees(geometry.calc_angle(cart[i], cart[c], cart[b]))
            zmat_ref[j][2] = 1
            zmat[j][2], collin = geometry.calc_dihedral(cart[i], cart[c], cart[b], cart[a])
        elif groupD[i] == 1:
            zmat_ref[j][0] = 4
            zmat[j][0] = np.linalg.norm(cart[i] - cart[d])
            zmat_ref[j][1] = 3
            zmat[j][1] = np.degrees(geometry.calc_angle(cart[i], cart[d], cart[c]))
            zmat_ref[j][2] = 2
            zmat[j][2], collin = geometry.calc_dihedral(cart[i], cart[d], cart[c], cart[b])
        j += 1

    return zmat_atom, zmat_ref, zmat, zmatorder

def make_zmat_from_cart_all_dihedrals(bond, cycle, dihed, conf_dihed, natom, atom, cart, mode):
    """
    Rearrange geometry defined in Cartesian into a Z-matrix,
    with references suitable for a 1-D hindered rotor scan.
    Include all rotors which are:
    If mode = 0: all rotatable bonds
    If mode = 1: only those bonds, which generate conformers,
    i.e. no methyl bonds, t-butyl bonds, etc

    bond: bond matrix of the species
    cycle: atom list of species with 0 if atom is not in cycle and 1 otherwise
    dihed: total list of dihedrals
    conf_dihed: list of dihedrals without the symmetrical ones
    natom: number of atoms
    atom: symbols of all the atoms
    cart: cartesian coordinates of the species
    """
    rotors = []  # rotors to consider

    if mode == 0:
        for rotor in dihed:
            rotors.append(rotor[:])
    if mode == 1:
        for rotor in conf_dihed:
            rotors.append(rotor[:])

    if len(rotors) > 0:
        # order in which the atoms will be added:
        # 1. First rotor
        # 2. Smallest path between first and second rotor (if any)
        # 3. Second rotor
        # 4. Smallest path between second and third rotor (if any)
        # ...
        # n. Last rotor
        # n+1. All atoms that are not part of a rotor,
        #      starting by the neighbors of all rotor atoms

        # order the rotors as such that the path between subsequent rotors does not
        # cross another rotor
        rotors, connecting_list, connected_rotor, path_length = order_rotors(rotors, bond, natom, atom)
        # Add first rotor:
        zmat_atom = ['X' for i in range(natom)]
        zmat_ref = np.zeros((natom, 3), dtype=int) - 1
        zmat = np.zeros((natom, 3)) - 1
        zmatorder = [-1 for i in range(natom)]

        zmat_atom[0] = atom[rotors[0][0]]
        zmatorder[0] = rotors[0][0]

        zmat_atom[1] = atom[rotors[0][1]]
        zmatorder[1] = rotors[0][1]
        zmat_ref[1][0] = 1
        zmat[1][0] = np.linalg.norm(cart[rotors[0][1]] - cart[rotors[0][0]])

        zmat_atom[2] = atom[rotors[0][2]]
        zmatorder[2] = rotors[0][2]
        zmat_ref[2][0] = 2
        zmat[2][0] = np.linalg.norm(cart[rotors[0][2]] - cart[rotors[0][1]])
        zmat_ref[2][1] = 1
        zmat[2][1] = np.degrees(geometry.calc_angle(cart[rotors[0][2]], cart[rotors[0][1]], cart[rotors[0][0]]))

        zmat_atom[3] = atom[rotors[0][3]]
        zmatorder[3] = rotors[0][3]
        zmat_ref[3][0] = 3
        zmat[3][0] = np.linalg.norm(cart[rotors[0][3]] - cart[rotors[0][2]])
        zmat_ref[3][1] = 2
        zmat[3][1] = np.degrees(geometry.calc_angle(cart[rotors[0][3]], cart[rotors[0][2]], cart[rotors[0][1]]))
        zmat_ref[3][2] = 1
        zmat[3][2], collin = geometry.calc_dihedral(cart[rotors[0][3]], cart[rotors[0][2]], cart[rotors[0][1]], cart[rotors[0][0]])

        # Add subsequent rotors:
        rot_index = 1
        j = 4

        while rot_index < len(rotors):
            rotor = rotors[rot_index]
            chain_length = path_length[rot_index - 1]
            chain = connecting_list[rot_index - 1]
            conn_rotor = connected_rotor[rot_index - 1]
            to_add = []
            added = []
            for at in rotor:
                if at not in zmatorder:
                    to_add.append(at)
                else:
                    added.append(at)

            if len(to_add) == 0:
                logger.error("error, all atoms of a rotor have been added, without adding this dihedral explicitly")
            elif len(to_add) == 1:
                at = to_add[0]
                if rotor.index(at) == 1 or rotor.index(at) == 2:
                    logger.error("error, all atoms except a middle atom have been added, strange...")
                elif rotor.index(at) == 0:
                    add(j, [at, rotor[1], rotor[2], rotor[3]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    j += 1
                elif rotor.index(at) == 3:
                    add(j, [at, rotor[2], rotor[1], rotor[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    j += 1
            elif len(to_add) == 2:
                at1 = to_add[0]
                at2 = to_add[1]
                if sorted([rotor.index(at1), rotor.index(at2)]) == [0, 1]:
                    neighbor_list = get_neighbors(rotor[2], bond, natom, rotors, zmatorder)
                    # add rotor[1]
                    add(j, [rotor[1], rotor[2], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    j += 1
                    # add rotor[0]
                    add(j, [rotor[0], rotor[1], rotor[2], rotor[3]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    j += 1
                elif sorted([rotor.index(at1), rotor.index(at2)]) == [2, 3]:
                    neighbor_list = get_neighbors(rotor[1], bond, natom, rotors, zmatorder)
                    # add rotor[2]
                    add(j, [rotor[2], rotor[1], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    j += 1
                    # add rotor[3]
                    add(j, [rotor[3], rotor[2], rotor[1], rotor[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    j += 1
                elif sorted([rotor.index(at1), rotor.index(at2)]) == [0, 3]:
                    r1, pos1 = middle_atom(rotors, rotor[1])
                    r2, pos2 = middle_atom(rotors, rotor[2])
                    if r1 and pos1[0] in zmatorder and pos1[1] in zmatorder:
                        # add rotor[0]
                        add(j, [rotor[0], rotor[1], pos1[0], pos1[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        # add rotor[3]
                        add(j, [rotor[3], rotor[2], rotor[1], rotor[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                    elif r2 and pos2[0] in zmatorder and pos2[1] in zmatorder:
                        # add rotor[3]
                        add(j, [rotor[3], rotor[2], pos2[0], pos2[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        # add rotor[0]
                        add(j, [rotor[0], rotor[1], rotor[2], rotor[3]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                    else:
                        if cycle[rotor[2]] == 1:
                            neighbor_list = get_neighbors(rotor[2], bond, natom, rotors, zmatorder)
                            # add rotor[3]
                            add(j, [rotor[3], rotor[2], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            # add rotor[0]
                            add(j, [rotor[0], rotor[1], rotor[2], rotor[3]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                        else:
                            neighbor_list = get_neighbors(rotor[1], bond, natom, rotors, zmatorder)
                            # add rotor[0]
                            add(j, [rotor[0], rotor[1], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            # add rotor[3]
                            add(j, [rotor[3], rotor[2], rotor[1], rotor[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                elif sorted([rotor.index(at1), rotor.index(at2)]) == [0, 2]:
                    neighbor_list = get_neighbors(rotor[1], bond, natom, rotors, zmatorder)
                    # add rotor[2]
                    add(j, [rotor[2], rotor[1], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    # add rotor[0]
                    add(j, [rotor[0], rotor[1], rotor[2], rotor[3]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    j += 1
                elif sorted([rotor.index(at1), rotor.index(at2)]) == [1, 3]:
                    neighbor_list = get_neighbors(rotor[0], bond, natom, rotors, zmatorder)
                    # add rotor[1]
                    add(j, [rotor[1], rotor[0], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    # add rotor[3]
                    add(j, [rotor[3], rotor[2], rotor[1], rotor[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    j += 1
                else:
                    logger.error("error, two atoms that need to be added are not the outer atoms of a rotor, strange...")
            elif len(to_add) == 3:
                if added[0] == rotor[0]:
                    neighbor_list = get_neighbors(rotor[0], bond, natom, rotors, zmatorder)
                    add(j, [rotor[1], rotor[0], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    j += 1
                    add(j, [rotor[2], rotor[1], rotor[0], neighbor_list[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    j += 1
                    add(j, [rotor[3], rotor[2], rotor[1], rotor[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    j += 1
                elif added[0] == rotor[1]:
                    neighbor_list = get_neighbors(rotor[1], bond, natom, rotors, zmatorder)
                    if cycle[rotor[3]] == 1 or cycle[rotor[2]] == 1:
                        add(j, [rotor[2], rotor[1], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[3], rotor[2], rotor[1], neighbor_list[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[0], rotor[1], rotor[2], rotor[3]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                    else:
                        add(j, [rotor[0], rotor[1], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[2], rotor[1], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[3], rotor[2], rotor[1], rotor[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                elif added[0] == rotor[2]:
                    neighbor_list = get_neighbors(rotor[2], bond, natom, rotors, zmatorder)
                    if cycle[rotor[3]] == 1 or cycle[rotor[2]] == 1:
                        add(j, [rotor[3], rotor[2], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[1], rotor[2], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[0], rotor[1], rotor[2], rotor[3]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                    else:
                        add(j, [rotor[1], rotor[2], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[0], rotor[1], rotor[2], neighbor_list[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[3], rotor[2], rotor[1], rotor[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                elif added[0] == rotor[3]:
                    neighbor_list = get_neighbors(rotor[3], bond, natom, rotors, zmatorder)
                    add(j, [rotor[2], rotor[3], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    j += 1
                    add(j, [rotor[1], rotor[2], rotor[3], neighbor_list[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    j += 1
                    add(j, [rotor[0], rotor[1], rotor[2], rotor[3]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    j += 1
                else:
                    logger.error('error')
            else:
                if len(chain) == 2:  # this is the minimum chain length
                    neighbor_list = get_neighbors(chain[1], bond, natom, rotors, zmatorder)
                    if chain[0] == rotor[0]:
                        add(j, [rotor[0], chain[1], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[1], rotor[0], chain[1], neighbor_list[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[2], rotor[1], rotor[0], chain[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[3], rotor[2], rotor[1], rotor[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                    elif chain[0] == rotor[1]:
                        if cycle[rotor[3]] == 1 or cycle[rotor[2]] == 1:
                            add(j, [rotor[1], chain[1], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[2], rotor[1], chain[1], neighbor_list[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[3], rotor[2], rotor[1], chain[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[0], rotor[1], rotor[2], rotor[3]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                        else:
                            add(j, [rotor[1], chain[1], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[0], rotor[1], chain[1], neighbor_list[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[2], rotor[1], chain[1], neighbor_list[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[3], rotor[2], rotor[1], rotor[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                    elif chain[0] == rotor[2]:
                        if cycle[rotor[3]] == 1 or cycle[rotor[2]] == 1:
                            add(j, [rotor[2], chain[1], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[1], rotor[2], chain[1], neighbor_list[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[3], rotor[2], chain[1], neighbor_list[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[0], rotor[1], rotor[2], rotor[3]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                        else:
                            add(j, [rotor[2], chain[1], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[1], rotor[2], chain[1], neighbor_list[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[0], rotor[1], rotor[2], chain[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[3], rotor[2], rotor[1], rotor[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                    elif chain[0] == rotor[3]:
                        add(j, [rotor[3], chain[1], neighbor_list[0], neighbor_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[2], rotor[3], chain[1], neighbor_list[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[1], rotor[2], rotor[3], chain[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[0], rotor[1], rotor[2], rotor[3]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                elif len(chain) > 2:
                    neighbor_list = get_neighbors(chain[-1], bond, natom, rotors, zmatorder)
                    # add the chain of atoms between this rotor and the closest rotor that is already in the zmatrix
                    for i in range(1, len(chain) - 1):  # do not consider two outer atoms, which are part of a rotor
                        ref_list = []
                        if i == 1:
                            ref_list = [chain[-1], neighbor_list[0], neighbor_list[1]]
                        elif i == 2:
                            ref_list = [chain[-2], chain[-1], neighbor_list[0]]
                        else:
                            ref_list = [chain[-i], chain[-(i - 1)], chain[-(i - 2)]]
                        add(j, [chain[-(i + 1)], ref_list[0], ref_list[1], ref_list[2]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                    # add the rotor itself:
                    ref_list = []

                    if len(chain) == 3:  # only one atom connects the two rotors:
                        ref_list = [chain[1], chain[2], neighbor_list[0]]
                    else:
                        ref_list = [chain[1], chain[2], chain[3]]

                    if chain[0] == rotor[0]:
                        add(j, [rotor[0], ref_list[0], ref_list[1], ref_list[2]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[1], rotor[0], ref_list[0], ref_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[2], rotor[1], rotor[0], ref_list[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[3], rotor[2], rotor[1], rotor[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                    elif chain[0] == rotor[1]:
                        if cycle[rotor[3]] == 1 or cycle[rotor[2]] == 1:
                            add(j, [rotor[1], ref_list[0], ref_list[1], ref_list[2]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[2], rotor[1], ref_list[0], ref_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[3], rotor[2], rotor[1], ref_list[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[0], rotor[1], rotor[2], rotor[3]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                        else:
                            add(j, [rotor[1], ref_list[0], ref_list[1], ref_list[2]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[0], rotor[1], ref_list[0], ref_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[2], rotor[1], ref_list[0], ref_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[3], rotor[2], rotor[1], rotor[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                    elif chain[0] == rotor[2]:
                        if cycle[rotor[3]] == 1 or cycle[rotor[2]] == 1:
                            add(j, [rotor[2], ref_list[0], ref_list[1], ref_list[2]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[1], rotor[2], ref_list[0], ref_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[3], rotor[2], ref_list[0], ref_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[0], rotor[1], rotor[2], rotor[3]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                        else:
                            add(j, [rotor[2], ref_list[0], ref_list[1], ref_list[2]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[1], rotor[2], ref_list[0], ref_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[0], rotor[1], rotor[2], ref_list[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                            add(j, [rotor[3], rotor[2], rotor[1], rotor[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                            j += 1
                    elif chain[0] == rotor[3]:
                        add(j, [rotor[3], ref_list[0], ref_list[1], ref_list[2]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[2], rotor[3], ref_list[0], ref_list[1]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[1], rotor[2], rotor[3], ref_list[0]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                        add(j, [rotor[0], rotor[1], rotor[2], rotor[3]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                        j += 1
                else:
                    logger.error('error, chain length too small')
            rot_index += 1
    else:
        # no rotors here, add the first four atoms of the molecule
        zmat_atom = ['X' for i in range(natom)]
        zmat_ref = np.zeros((natom, 3), dtype=int) - 1
        zmat = np.zeros((natom, 3)) - 1
        zmatorder = [-1 for i in range(natom)]

        if len(atom) == 1:
            zmat_atom[0] = atom[0]
            zmatorder[0] = 0
        elif len(atom) == 2:
            zmat_atom[0] = atom[0]
            zmatorder[0] = 0

            zmat_atom[1] = atom[1]
            zmatorder[1] = 1
            zmat_ref[1][0] = 1
            zmat[1][0] = np.linalg.norm(cart[1] - cart[0])
        else:
            # get three atoms that are bonded to one another
            motif = ['X', 'X', 'X']
            ins = start_motif(motif, natom, bond, atom, -1, [[k] for k in range(natom)])[0]

            zmat_atom[0] = atom[ins[0]]
            zmatorder[0] = ins[0]

            zmat_atom[1] = atom[ins[1]]
            zmatorder[1] = ins[1]
            zmat_ref[1][0] = 1
            zmat[1][0] = np.linalg.norm(cart[ins[1]] - cart[ins[0]])

            zmat_atom[2] = atom[ins[2]]
            zmatorder[2] = ins[2]
            zmat_ref[2][0] = 2
            zmat[2][0] = np.linalg.norm(cart[ins[2]] - cart[ins[1]])
            zmat_ref[2][1] = 1
            zmat[2][1] = np.degrees(geometry.calc_angle(cart[ins[2]], cart[ins[1]], cart[ins[0]]))

            j = 3

    # add all the remaining atoms
    while -1 in zmatorder:
        for i in range(natom):
            if i not in zmatorder:
                neighbor_list = get_three_neighbors(i, bond, natom, rotors, zmatorder)
                if len(neighbor_list) > 2:
                    add(j, [i, neighbor_list[0], neighbor_list[1], neighbor_list[2]], zmat, zmat_atom, zmatorder, zmat_ref, atom, cart)
                    j += 1

    return zmat_atom, zmat_ref, zmat, zmatorder


def get_three_neighbors(at, bond, natom, rotors, zmatorder):
    """
    Generate a list of three neighbors for atom at
    (only works if at is not in the zmat yet!!
    """
    # get three atoms that are in zmat_order
    list = []
    for k in range(natom):
        if k in zmatorder:
            if bond[at][k] > 0:
                rotor, pos = middle_atom(rotors, k)
                if rotor:
                    list.append(k)
                    list.append(pos[0])
                    list.append(pos[1])
                    return list
    for k in range(natom):
        if k in zmatorder:
            if bond[at][k] > 0:
                for l in range(natom):
                    if l in zmatorder:
                        if bond[k][l] > 0:
                            if l != at:
                                # check for neighbors of l
                                for m in range(natom):
                                    if m in zmatorder > -1:
                                        if bond[l][m] > 0:
                                            if m != k and m != at:
                                                list.append(k)
                                                list.append(l)
                                                list.append(m)
                                                return list
                                # check for another neighbor of k
                                for m in range(natom):
                                    if m in zmatorder:
                                        if bond[k][m] > 0:
                                            if m != l and m != at:
                                                list.append(k)
                                                list.append(l)
                                                list.append(m)
                                                return list
    return list


def get_neighbors(at, bond, natom, rotors, zmatorder):
    """
    Get two neighbors or at
    """
    # get three atoms that are in zmat_order
    list = []
    rotor, pos = middle_atom(rotors, at)
    if rotor:
        if pos[0] in zmatorder and pos[1] in zmatorder:
            list.append(pos[0])
            list.append(pos[1])
            return list
    for k in range(natom):
        if k in zmatorder:
            if bond[at][k] > 0:
                rotor, pos = middle_atom(rotors, k)
                if rotor:
                    if pos[0] in zmatorder:
                        if at != pos[0]:
                            list.append(k)
                            list.append(pos[0])
                            return list
                    if pos[1] in zmatorder:
                        if at != pos[1]:
                            list.append(k)
                            list.append(pos[1])
                            return list
    for k in range(natom):
        if k in zmatorder:
            if bond[at][k] > 0:
                for l in range(natom):
                    if l in zmatorder:
                        if bond[k][l] > 0:
                            if l != at:
                                list.append(k)
                                list.append(l)
                                return list
    return list


def middle_atom(rotors, k):
    """
    Verify if the kth atom is in a rotor and return two appropriate zmat neighbors
    return its neighbors
    """

    for rot in rotors:
        if k == rot[1]:
            return 1, [rot[0], rot[2]]
        if k == rot[2]:
            return 1, [rot[1], rot[3]]

    for rot in rotors:
        if k == rot[0]:
            return 1, [rot[1], rot[2]]
        if k == rot[3]:
            return 1, [rot[2], rot[1]]

    return 0, []


def add(j, list, zmat, zmat_atom, zmatorder, zmat_ref, atom, cart):
    """
    add an atom on the jth position of the zmat according to the four atoms in the list
    """
    zmat_atom[j] = atom[list[0]]
    zmatorder[j] = list[0]
    zmat_ref[j][0] = zmatorder.index(list[1]) + 1
    zmat[j][0] = np.linalg.norm(cart[list[0]] - cart[list[1]])
    zmat_ref[j][1] = zmatorder.index(list[2]) + 1
    zmat[j][1] = np.degrees(geometry.calc_angle(cart[list[0]], cart[list[1]], cart[list[2]]))
    zmat_ref[j][2] = zmatorder.index(list[3]) + 1
    zmat[j][2], collin = geometry.calc_dihedral(cart[list[0]], cart[list[1]], cart[list[2]], cart[list[3]])


def order_rotors(rotors, bond, natom, atom):
    """
    Order rotors such that the shortest path between two subsequent rotors
    does not involve atoms of other rotors
    path_length: negative in case of overlap, 0 in case of directly bonded
    and positive in case of a pathway between both subsequent rotors
    connecting list: empty list in case of overlap or if directly bonded,
    chain atoms between to rotors in all other cases.
    """

    ordered = [rotors[0]]
    path_length = []
    connecting_list = []
    connected_rotor = []

    while len(ordered) < len(rotors):
        min = 999
        list = []
        rotor = []
        for r1 in rotors:
            new = 1
            for r2 in ordered:
                if r1 == r2:
                    new = 0
            if new:
                for r2 in ordered:
                    m, l = get_minimum_pathway(bond, r1, r2, natom, atom)
                    if m < min:
                        rotor = r2
                        next = r1
                        min = m
                        list = l
        connected_rotor.append(rotor)
        path_length.append(m)
        ordered.append(next)
        connecting_list.append(list)

    return ordered, connecting_list, connected_rotor, path_length


def get_minimum_pathway(bond, r1, r2, natom, atom):
    """
    This method gets the minimum pathway between r1 and r2
    if r1 and r2 overlaps, it returns a negative number with the number of overlapping atoms
    if the two rotors are bonded with no overlap, it returns the value 0
    if the rotors are not bonded, it returns the length of the minimum pathway between both,
    and the pathway atoms
    """

    overlap = []
    for a1 in r1:
        for a2 in r2:
            if a1 == a2:
                overlap.append(a1)
    if len(overlap) > 0:
        return -len(overlap), overlap

    # the rotors do not overlap, check for direct bonds
    for a1 in r1:
        for a2 in r2:
            if bond[a1][a2] > 0:
                return 0, [a1, a2]

    # no bonds, check for the shortest pathway
    for i in range(3, natom):
        motif = ['X' for j in range(i)]
        instances = start_motif(motif, natom, bond, atom, -1, [[k] for k in range(natom)])
        for ins in instances:
            if ins[0] in r1 and ins[-1] in r2:
                return i - 2, ins
            if ins[0] in r2 and ins[-1] in r1:
                ins.reverse()
                return i, ins


def estoktp_zmat(species, fname):
    """
    Write a zmat to a file as such that all the rotors are independently defined
    This can be used by EStoKTp to get to a rate coefficient
    Also, visualizations of the rotors are created to verify their correctness
    """

    natom = species.natom
    atom = species.atom

    if not os.path.exists('rotors'):
        os.mkdir('rotors')
    zmat_atom, zmat_ref, zmat, zmatorder = make_zmat_from_cart_all_dihedrals(species.bond, species.cycle, species.dihed, species.conf_dihed, natom, atom, species.geom, 0)
    with open('rotors/' + fname + '.zmat', 'w') as zfile:
        write_zmat_molden(zfile, zmat_atom, zmat_ref, zmat)
        zfile.write('// All dihedrals\n')
        for dih in species.dihed:
            zfile.write('// {}\n'.format(' '.join([str(zmatorder.index(di) + 1) for di in dih])))
        zfile.write('\n// Dihedrals for conformational scan\n')
        for dih in species.conf_dihed:
            zfile.write('// {}\n'.format(' '.join([str(zmatorder.index(di) + 1) for di in dih])))
        zfile.write('\n\n')

    for i, rotor in enumerate(species.dihed):
        # write a geometry file
        # get the appropriate dihedral position
        for j in range(3, natom):
            indices = [zmatorder[j], zmatorder[zmat_ref[j][0] - 1], zmatorder[zmat_ref[j][1] - 1], zmatorder[zmat_ref[j][2] - 1]]
            if indices == rotor or indices == rotor[::-1]:
                break

        fname = 'rotors/' + fname + '_' + str(i) + '.xyz'
        with open(fname, 'w') as f:
            # do 36 iterations of 10 degrees
            for k in range(36):
                zmat_copy = copy.deepcopy(zmat)
                zmat_copy[j][2] += float(k * 10)
                cart = make_cart_from_zmat(zmat_copy, zmat_atom, zmat_ref, natom, atom, zmatorder)
                cart = translate_and_rotate(cart, rotor[1], rotor[2])

                f.write(write_cart(cart, atom))

    return 0


def write_zmat(zfile, zmat_atom, zmat_ref, zmat):
    """
    Write a Z-matrix in Gaussian format into zfile.
    File is open already and is writeable.
    """
    zfile.write(str(zmat_atom[0]) + '\n')
    zfile.write(str(zmat_atom[1]) + '  ' + str(zmat_ref[1][0]) + '  ' + str(zmat[1][0]) + '\n')
    zfile.write(str(zmat_atom[2]) + '  ' + str(zmat_ref[2][0]) + '  ' + str(zmat[2][0]) + '  ' + str(zmat_ref[2][1]) + '  ' + str(zmat[2][1]) + '\n')
    for i in range(3, len(zmat_atom)):
        zfile.write(str(zmat_atom[i]) + '  ' + str(zmat_ref[i][0]) + '  ' + str(zmat[i][0]) + '  ' + str(zmat_ref[i][1]) + '  ' + str(zmat[i][1]) + '  ' + str(zmat_ref[i][2]) + '  ' + str(zmat[i][2]) + '\n')

    return 0


def write_zmat_molden(zfile, zmat_atom, zmat_ref, zmat):
    """
    Write a Z-matrix in Molden format into zfile.
    File is open already and is writeable.
    """

    zfile.write('zmat angstroms\n')
    zfile.write(str(zmat_atom[0]) + '\n')

    if len(zmat_atom) > 1:
        zfile.write(str(zmat_atom[1]) + '  ' + str(zmat_ref[1][0]) + '  dist1\n')

    if len(zmat_atom) > 2:
        zfile.write(str(zmat_atom[2]) + '  ' + str(zmat_ref[2][0]) + '  dist2  ' + str(zmat_ref[2][1]) + '  ang2\n')

    for i in range(3, len(zmat_atom)):
        zfile.write(str(zmat_atom[i]) + '  ' + str(zmat_ref[i][0]) + '  dist' + str(i) + '  ' + str(zmat_ref[i][1]) + '  ang' + str(i) + '  ' + str(zmat_ref[i][2]) + '  dih' + str(i) + '\n')

    zfile.write('variables\n')
    if len(zmat_atom) > 1:
        zfile.write('dist1     ' + str(zmat[1][0]) + '\n')
    if len(zmat_atom) > 2:
        zfile.write('dist2     ' + str(zmat[2][0]) + '\n')
        zfile.write('ang2      ' + str(zmat[2][1]) + '\n')
    for i in range(3, len(zmat_atom)):
        zfile.write('dist' + str(i) + '     ' + str(zmat[i][0]) + '\n')
        zfile.write('ang' + str(i) + '      ' + str(zmat[i][1]) + '\n')
        zfile.write('dih' + str(i) + '      ' + str(zmat[i][2]) + '\n')

    zfile.write('constants\n')
    zfile.write('end\n\n')

    return 0


def write_cart(geom, atom):
    s = '%i\n' % len(geom)
    s += 'Energy = 1.0\n'

    for index in range(len(geom)):
        s += '%s %.6f %.6f %.6f\n' % (atom[index], geom[index][0], geom[index][1], geom[index][2])
    return s


def make_cart_from_zmat(zmat, zmat_atom, zmat_ref, natom, atom, zmatorder):
    """
    Create Cartesian coordinates from a Z-matrix representation.
    Let's assume that we are dealing with atom D, whose distance r is defined relative to C,
    angle theta relative to B, and dihedral phi relative to A.
    First we place D on the B-C axis at r distance from C.
    Then D is rotated in the ABC plance by theta, and
    then around the B-C axis by phi.
    All this is done by a general rotation matrix formalism.
    The first three atoms are special.
    """

    cart = np.zeros((natom, 3))
    zm = copy.deepcopy(zmat)

    for i in range(natom):
        zm[i][1] = np.radians(zm[i][1])
        zm[i][2] = np.radians(zm[i][2])

    # A
    cart[0] = [0., 0., 0.]
    if natom == 1:
        return 1

    if len(zmat_atom) > 1:
        # B
        cart[1] = [zm[1][0], 0., 0.]
        if natom == 2:
            return 2

        if len(zmat_atom) > 2:
            # C
            cart[2][0] = np.sign(np.pi / 2 - zm[2][1]) * zm[2][0] * np.cos(zm[2][1]) + cart[1][0]
            cart[2][1] = zm[2][0] * np.sin(zm[2][1])
            if natom == 3:
                return 3

            # D
            for i in range(3, len(zmat_atom)):
                c = zmat_ref[i][0] - 1  # distance
                b = zmat_ref[i][1] - 1  # angle
                a = zmat_ref[i][2] - 1  # dihedral

                # D is placed parallel to the B-C axis, relative to A (the origin)
                # B->C vector
                bc = [cart[c][j] - cart[b][j] for j in range(3)]
                bc = bc / np.linalg.norm(bc)
                # y/x = p
                # z/x = q
                # x^2 + y^2 + z^2 = r^2
                # therefore x = sqrt(r ^ 2 / (p ^ 2 + q ^ 2 + 1))
                p = bc[1] / bc[0]
                q = bc[2] / bc[0]
                r = zm[i][0]
                if cart[c][0] > cart[b][0]:
                    x = np.sqrt(r * r / (p * p + q * q + 1))
                else:
                    x = -np.sqrt(r * r / (p * p + q * q + 1))
                cart[i] = [x, x * p, x * q]

                # A->B vector
                ab = [cart[b][j] - cart[a][j] for j in range(3)]

                # |AB x BC|
                n = np.cross(ab, bc) / np.linalg.norm(np.cross(ab, bc))

                # rotation around the normal to ABC by the angle
                th = zm[i][1] - np.pi
                cart[i] = geometry.rotate_atom(cart[i], n, th)

                # rotation around the BC axis by the dihedral angle
                th = zm[i][2] - np.pi
                cart[i] = geometry.rotate_atom(cart[i], bc, th)

                # shift the vector to C
                cart[i] = [cart[i][j] + cart[c][j] for j in range(3)]

    cart_reordered = []
    for i in range(natom):
        cart_reordered.append(cart[zmatorder.index(i)])

    return cart_reordered


def main():
    """
    zmat to cart conversion, and vice versa.
    """


if __name__ == "__main__":
    main()
