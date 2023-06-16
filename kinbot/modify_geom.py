"""
This class modifies a given geometry according to a set of coordinates that
need to have a new value assigned.

The optimization is done based on interatomic distances only
The deviations of the the distances are weighted by the inverse of the distance itself
"""
import os
import copy
import logging
import time
import numpy as np

from ase import Atoms
from ase.io import write
from ase.calculators.singlepoint import SinglePointCalculator

from kinbot import bfgs
from kinbot import constants
from kinbot import find_motif
from kinbot import geometry
from kinbot import zmatrix

logger = logging.getLogger('KinBot')


class cost_function():
    def __init__(self, coords):
        self.coords = coords

    def eval(self, x):
        """
        x is a vector of length 3N with N the number of atoms
        containing the cartesian coordinates [x1, y1, z1, x2, ..., xN, yN, zN]
        """
        e = 0
        for coord in self.coords:
            i = coord[0]
            j = coord[1]
            d = coord[2]
            weight = coord[3]
            dc = (x[3 * i] - x[3 * j]) ** 2 + (x[3 * i + 1] - x[3 * j + 1]) ** 2 + (x[3 * i + 2] - x[3 * j + 2]) ** 2
            if len(coord) == 5:
                if dc < d:
                    # add the one-sided potential
                    e += ((dc - d) * weight) ** 2
            else:
                e += ((dc - d) * weight) ** 2
        return e

    def gradient(self, x):
        grad = np.zeros(len(x))
        for coord in self.coords:
            i = coord[0]
            j = coord[1]
            d = coord[2]
            weight = coord[3]
            dc = (x[3 * i] - x[3 * j]) ** 2 + (x[3 * i + 1] - x[3 * j + 1]) ** 2 + (x[3 * i + 2] - x[3 * j + 2]) ** 2
            if len(coord) == 5:
                if dc < d:
                    grad[3 * i] += 2 * ((dc - d) * weight) * 2 * (x[3 * i] - x[3 * j])
                    grad[3 * i + 1] += 2 * ((dc - d) * weight) * 2 * (x[3 * i + 1] - x[3 * j + 1])
                    grad[3 * i + 2] += 2 * ((dc - d) * weight) * 2 * (x[3 * i + 2] - x[3 * j + 2])

                    grad[3 * j] += 2 * ((dc - d) * weight) * 2 * (x[3 * i] - x[3 * j]) * -1
                    grad[3 * j + 1] += 2 * ((dc - d) * weight) * 2 * (x[3 * i + 1] - x[3 * j + 1]) * -1
                    grad[3 * j + 2] += 2 * ((dc - d) * weight) * 2 * (x[3 * i + 2] - x[3 * j + 2]) * -1
            else:
                grad[3 * i] += 2 * ((dc - d) * weight) * 2 * (x[3 * i] - x[3 * j])
                grad[3 * i + 1] += 2 * ((dc - d) * weight) * 2 * (x[3 * i + 1] - x[3 * j + 1])
                grad[3 * i + 2] += 2 * ((dc - d) * weight) * 2 * (x[3 * i + 2] - x[3 * j + 2])

                grad[3 * j] += 2 * ((dc - d) * weight) * 2 * (x[3 * i] - x[3 * j]) * -1
                grad[3 * j + 1] += 2 * ((dc - d) * weight) * 2 * (x[3 * i + 1] - x[3 * j + 1]) * -1
                grad[3 * j + 2] += 2 * ((dc - d) * weight) * 2 * (x[3 * i + 2] - x[3 * j + 2]) * -1

        return grad


def append_geom(natom, step, new_e, atom, x_new, grad, atoms_list, f_out=None):
    if f_out is not None:
        f_out.write('{}\nPoint  {} Energy=  {}\n'.format(natom, step, new_e))
        for at in range(natom):
            f_out.write(atom[at] + ' ')
            [f_out.write(str(np.reshape(x_new, (natom, 3))[at][i]) + '  ') for i in range(3)]
            f_out.write('\n')
        step += 1

    atoms = Atoms(symbols=atom, positions=np.reshape(x_new, (natom, 3)))
    calc = SinglePointCalculator(atoms, energy=new_e, forces=10. * np.reshape(grad, (natom, 3)))
    atoms.set_calculator(calc)
    atoms_list.append(atoms)
    return step


def modify_coordinates(species, name, geom, changes, bond, write_files=0):
    """
    Geom is the geometry (n x 3 matrix with n the number of atoms)
    in cartesian coordinates
    Changes is a list of lists, each list containing the coordinates
    and their new value (atom indices start at 0):
    To change a bond length: [atom1, atom2, bond_length]
    To change a bond angle: [neighbor1, central_atom, neighbor2,
                             angle_in_degrees]
    To change a dihedral angle: [atom1, atom2, atom3, atom4,
                                 dihedarl_angle_in_degrees]

    Bond is the bond matrix of the molecule
    """

    start_time = time.time()
    logger.debug('Starting coordinate modification for {}'.format(name))
    logger.debug('Changes:')
    for c in changes:
        logger.debug('\t{}'.format('\t'.join([str(ci) for ci in c])))

    step = 1
    atoms_list = []

    count = 0
    fname = '{}_{}.xyz'.format(name, count)
    while os.path.exists(fname):
        count += 1
        fname = '{}_{}.xyz'.format(name, count)
    f_out = None
    if write_files:
        f_out = open(fname, 'w')

    new_geom = copy.deepcopy(geom)
    step = append_geom(species.natom, step, 0., species.atom, new_geom, np.zeros((species.natom * 3)), atoms_list, f_out=f_out)

    # change dihedrals, if necessary
    for ci in changes:
        if len(ci) == 5:
            zmat_atom, zmat_ref, zmat, zmatorder = zmatrix.make_zmat_from_cart(species, ci[:-1], new_geom, 2)
            # write_zmat(zmat_atom, zmat_ref, zmat, new_geom, species.atom)
            orig_dih = zmat[3][2]
            new_dih = ci[-1]
            dih_diff = new_dih - orig_dih
            zmat[3][2] += dih_diff
            for i in range(4, species.natom):
                if zmat_ref[i][2] == 4:
                    zmat[i][2] += dih_diff
                if zmat_ref[i][2] == 1:
                    zmat[i][2] += dih_diff
            new_geom = zmatrix.make_cart_from_zmat(zmat, zmat_atom, zmat_ref, species.natom, species.atom, zmatorder)
            # write_zmat(zmat_atom, zmat_ref, zmat, new_geom, species.atom)
            step = append_geom(species.natom, step, 0., species.atom, new_geom, np.zeros((species.natom * 3)), atoms_list, f_out=f_out)
        # change angles, if necessary
        if len(ci) == 4:
            # original angle in radians
            orig_angle = geometry.calc_angle(new_geom[ci[0]], new_geom[ci[1]], new_geom[ci[2]])
            new_angle = np.radians(ci[-1])  # new angle in radians

            v1 = new_geom[ci[0]] - new_geom[ci[1]]
            v2 = new_geom[ci[2]] - new_geom[ci[1]]
            rot_ax = [0., 0., 0.]

            # create a vector perpendicular to v1 and v2
            # verify if points are collinear
            if np.linalg.norm(np.cross(v1, v2)) == 0:
                # rotate around any axis perpendicular to the axis along the three points:
                if v1[0] != 0 or v1[1] != 0:
                    rot_ax = [v1[1], -v1[0], 0.]
                elif v1[0] != 0 or v1[2] != 0:
                    rot_ax = [v1[2], 0., -v1[0]]
                else:
                    rot_ax = [1., 0., 0.]
            else:
                rot_ax = np.cross(v1, v2)

            rot_ax = rot_ax / np.linalg.norm(rot_ax)
            # rotate all the atoms on the side of the last atom
            st, ats, ats2 = divide_atoms(ci[2], ci[1], bond, species.natom, species.atom)
            if not st:
                break
            for atj in ats:
                new_geom[atj] = perform_rotation(new_geom[atj], new_geom[ci[1]], rot_ax, new_angle - orig_angle)
                step = append_geom(species.natom, step, 1., species.atom, new_geom, np.zeros((species.natom * 3)), atoms_list, f_out=f_out)

    coords = get_coords(species, bond, new_geom, changes, 0)
    # optimize the geometry to meet the coords list
    x0 = np.reshape(new_geom, 3 * species.natom)
    cost_fct = cost_function(coords)
    logger.debug('Starting BFGS')
    gs = ''  # initial geomtry string
    for i, at in enumerate(species.atom):
        x, y, z = new_geom[i]
        gs += '{}, {:.8f}, {:.8f}, {:.8f}, \n'.format(at, x, y, z)
    logger.debug("For the following initial geometry:\n" + gs)

    opt = bfgs.BFGS()
    x_opt, x_i, g_i = opt.optimize(cost_fct, x0)

    new_geom = np.reshape(x_opt, (species.natom, 3))
    for i, xi in enumerate(x_i):
        geomi = np.reshape(xi, (species.natom, 3))
        gradi = np.reshape(g_i[i], (species.natom, 3))
        step = append_geom(species.natom, step, 2., species.atom, geomi, gradi, atoms_list, f_out=f_out)

    if write_files:
        write(fname.replace('.xyz', '.traj'), atoms_list)
        f_out.close()

    success = control_changes(species, name, geom, new_geom, changes, bond)

    end_time = time.time()
    logger.debug('Finished coordinate changes after {:.2f} seconds'.format(end_time - start_time))
    return success, new_geom


def control_changes(species, name, geom, new_geom, changes, bond):
    """
    This method controls three things:
    1. the changes are all correct
    2. the bonded atom pairs that are not part of the changes have the same bond length
    3. the non-bonded atom pairs are far enough from each other
    """
    success = 1
    gs = ''  # initial geomtry string
    for i, at in enumerate(species.atom):
        x, y, z = geom[i]
        gs += '{}, {:.8f}, {:.8f}, {:.8f}, \n'.format(at, x, y, z)

    cs = ''  # changes list
    for ci in changes:
        cs += '[{}]\n'.format(', '.join([str(cij) for cij in ci]))

    # check if the changes are good
    for ci in changes:
        if len(ci) == 3:
            bond_length = np.linalg.norm(new_geom[ci[0]] - new_geom[ci[1]])
            if np.abs(bond_length - ci[2]) > 0.05:  # use a 0.05 Angstrom cutoff
                logger.debug("The modified bond length is not correct")
                logger.debug("Expected {}, got {}".format(ci[2], bond_length))
                logger.debug("Species name: " + name)
                logger.debug("For the following initial geometry:\n" + gs)
                logger.debug("And the following change list:\n" + cs)
                success = 0
        if len(ci) == 4:
            angle = geometry.calc_angle(new_geom[ci[0]], new_geom[ci[1]], new_geom[ci[2]])
            change_angle = np.radians(ci[3])
            if np.abs(angle - change_angle) > 0.05:  # use a 0.05 radians cutoff
                logger.debug("The modified angle is not correct")
                logger.debug("Expected {}, got {}".format(change_angle, angle))
                logger.debug("Species name: " + name)
                logger.debug("For the following initial geometry:\n" + gs)
                logger.debug("And the following change list:\n" + cs)
                success = 0

    # check if the bond lengths are good:
    for i in range(species.natom - 1):
        for j in range(i + 1, species.natom):
            if bond[i][j] > 0:
                is_change = 0
                for ci in changes:
                    if i in ci and j in ci:
                        is_change = 1
                if not is_change:
                    new_bond = np.linalg.norm(new_geom[i] - new_geom[j])
                    orig_bond = np.linalg.norm(geom[i] - geom[j])
                    if np.abs(new_bond - orig_bond) > 0.1:  # use a 0.1 Angstrom cutoff
                        logger.debug("The bond length {} {} is not correct after modifications".format(i, j))
                        logger.debug("Expected {}, got {}".format(orig_bond, new_bond))
                        logger.debug("Species name: " + name)
                        logger.debug("For the following initial geometry:\n" + gs)
                        logger.debug("And the following change list:\n" + cs)
                        success = 0

    # check if non-bonded atoms are too close:
    for i in range(species.natom - 1):
        for j in range(i + 1, species.natom):
            if bond[i][j] == 0:
                dist = np.linalg.norm(new_geom[i] - new_geom[j])
                minb = constants.st_bond[''.join(sorted([species.atom[i], species.atom[j]]))]
                if dist < minb:
                    logger.debug("The atoms {} {} are too close after modifications".format(i, j))
                    logger.debug("Species name: " + name)
                    logger.debug("For the following initial geometry:\n" + gs)
                    logger.debug("And the following change list:\n" + cs)
                    success = 0

    return success


def write_zmat(zmat_atom, zmat_ref, zmat, geom, atom):
    i = 0
    zmat_file = 'zmatrix_' + str(i) + '.zmat'
    while os.path.exists(zmat_file):
        i += 1
        zmat_file = 'zmatrix_' + str(i) + '.zmat'

    with open(zmat_file, 'w') as ff:
        zmatrix.write_zmat_molden(ff, zmat_atom, zmat_ref, zmat)

    cart_file = zmat_file.replace('.zmat', '.xyz')
    with open(cart_file, 'w') as ff:
        ff.write(str(len(atom)) + '\n\n')
        for i, at in enumerate(atom):
            x, y, z = geom[i]
            ff.write('{} {:.8f} {:.8f} {:.8f}\n'.format(at, x, y, z))
        ff.write('\n')


def get_coords(species, bond, geom, changes, mode):
    """
    list the (N*(N-1)) / 2 possible bond lengths and the value we optimize to
    """
    natom = species.natom
    atom = species.atom
    coords = []
    for i in range(natom - 1):
        for j in range(i + 1, natom):
            is_change = 0  # the distance from i to j needs to be changed
            change = []
            for ci in changes:
                if [i, j] == sorted([ci[0], ci[-2]]):
                    is_change = 1
                    change = ci
            if is_change:
                if len(change) == 3:
                    coords.append([i, j, change[-1] ** 2, 1.])
                elif len(change) == 4:
                    # calculate the bond length that corresponds to the new angle
                    b1 = np.linalg.norm(geom[i] - geom[change[1]])
                    b2 = np.linalg.norm(geom[j] - geom[change[1]])
                    a = np.radians(change[-1])
                    d = b1 ** 2 + b2 ** 2 - 2 * b1 * b2 * np.cos(a)
                    coords.append([i, j, d, 10])
                elif len(change) == 5:
                    # take the current interatomic distance
                    d = np.linalg.norm(geom[i] - geom[j]) ** 2
                    coords.append([i, j, d, 10])
            else:
                if mode == 0:
                    d = np.linalg.norm(geom[i] - geom[j]) ** 2
                    d_min = (constants.st_bond[''.join(sorted([atom[i], atom[j]]))] * 1.2) ** 2
                    if np.sqrt(d) < 4.:  # use a cutoff of 4 angstroms
                        if bond[i][j] > 0:
                            coords.append([i, j, d, 1. / d])  # use a larger weight for bonds
                        else:
                            # check if i and j have the same neighbor
                            same_neighbor = 0
                            for k in range(natom):
                                if k != i and k != j:
                                    if bond[i][k] > 0 and bond[j][k] > 0:
                                        same_neighbor = 1
                            if same_neighbor:
                                coords.append([i, j, d, .5 / d])  # this is a shallow two-sided potential
                            else:
                                coords.append([i, j, d_min, 1. / d_min, 1])  # this is stronger one-sided potential

    return coords


def divide_atoms(ati, atj, bond, natom, atom):
    """
    This method divides the atoms in a molecule in two sets,
    which are separated by a bond
    In the case of rings, the atoms are equally divided in the two sets,
    which will change the bond length of the bond furthest away from
    the given bond.
    Be careful when using this method for cyclic structures!
    """
    status = 1
    if bond[ati, atj] == 0:
        return 0, [ati], []

    # Get all the atoms on the side of ati
    visited = [ati]
    forbidden = [atj]
    division1 = [ati]

    # check for cycles and cut them in half
    for ring_size in range(3, natom + 1):
        motif = ['X' for at in range(ring_size)]
        inst = find_motif.start_motif(motif, natom, bond, atom, -1, [])
        for ins in inst:
            if bond[ins[0]][ins[-1]] > 0:
                # cycle found
                if ins[0] == ati and ins[-1] == atj:
                    forbidden.append(ins[ring_size // 2])
                if ins[0] == atj and ins[-1] == ati:
                    forbidden.append(ins[- ring_size // 2 - 1])
        if len(inst) == 0:
            break

    get_neighbors(ati, visited, forbidden, division1, bond, natom)
    division2 = [x for x in range(natom) if x not in division1]

    return status, division1, division2


def get_neighbors(ati, visited, forbidden, division1, bond, natom):
    for atj in range(natom):
        if atj not in visited and atj not in forbidden:
            if bond[atj, ati] > 0:
                division1.append(atj)
                visited.append(atj)
                get_neighbors(atj, visited, forbidden, division1, bond, natom)


def perform_rotation(at, center, axis, angle):
    # move center to origin:
    at -= center

    # create rotation matrix:
    a = np.cos(angle / 2)
    b, c, d = -axis * np.sin(angle / 2)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    rot_matrix = ([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                   [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                   [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    # perform the rotation:
    at = np.dot(rot_matrix, at)

    # put the center back to its original coordinates:
    at += center

    return at
