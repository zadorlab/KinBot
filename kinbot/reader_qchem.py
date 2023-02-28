"""
Functions to read QChem output files.
"""
import logging

import numpy as np
from ase import Atoms

logger = logging.getLogger('KinBot')


def read_geom(outfile, mol=Atoms(), irc=False):
    """Read the final geometry from a QChem output file.
    """
    do_read, done = (False, False)
    with open(outfile) as f:
        for line in f:
            if 'I     Atom' in line:
                do_read = True
                geom = np.zeros((len(mol), 3))
            elif irc and 'IRC -- convergence criterion reached.' in line:
                break
            elif not do_read or '---------' in line:
                continue
            elif do_read and line.split()[0].isdecimal():
                geom[int(line.split()[0]) - 1] = np.array([float(c) for c in line.split()[2:5]])
                if int(line.split()[0]) == len(mol):
                    done = True
                    do_read = False
            else:
                raise IndexError('Unrecognised QChem output format')
    if not done:
        # Coordinates not found, assuming it is a SPE calculation.
        raise ValueError(f'Coordinates not found on {outfile}')
    return geom


def read_zpe(outfile):
    """Read the zpe
    """
    from kinbot import constants

    with open(outfile) as f:
        for line in f:
            if 'Zero point vibrational energy:' in line:
                zpe = line.split()[4]
    try:
        return float(zpe) * constants.KCALtoHARTREE
    except ValueError:
        logger.warning(f'Non-numeric ZPE: {zpe}')
        return np.NAN


def read_freq(outfile, atoms):
    """
    Read the frequencies
    """
    freqs = []
    natoms = len(atoms)  # filter out the dummy atoms
    if natoms == 1:
        return []

    with open(outfile) as f:
        for line in f:
            if 'Frequency:' in line:
                freqs.extend([float(fr) for fr in line.split()[1:]])
    return freqs


# def read_convergence(outfile):
#     """
#     Check for the four YES.
#     0: did not converge
#     1: forces and displacements converged
#     2: forces converged
#     """
#
#     with open(outfile) as f:
#         lines = f.readlines()
#
#     for n, line in enumerate(lines):
#         if 'Item               Value     Threshold  Converged?' in line:
#             if 'YES' in lines[n + 1]:
#                 if 'YES' in lines[n + 2]:
#                     if 'YES' in lines[n + 3]:
#                         if 'YES' in lines[n + 4]:
#                             return 1
#                     else:
#                         return 2
#
#     return 0  # will look through the whole file


def constraint(mol, fix, change):  # TODO
    """
    Convert constraints into PCBFGS constraints.
    """

    bonds = []
    angles = []
    dihedrals = []
    for fi in fix:
        if len(fi) == 2:
            # careful: atom indices in the fix lists start at 1
            bondlength = mol.get_distance(fi[0] - 1, fi[1] - 1)
            bonds.append([bondlength, [fi[0] - 1, fi[1] - 1]])
        if len(fi) == 3:
            # careful: atom indices in the fix lists start at 1
            angle = mol.get_angle(fi[0] - 1, fi[1] - 1,
                                  fi[2] - 1) * np.pi / 180.
            angles.append([angle, [fi[0] - 1, fi[1] - 1, fi[2] - 1]])
        if len(fi) == 4:
            # careful: atom indices in the fix lists start at 1
            dihed = mol.get_dihedral(fi[0] - 1, fi[1] - 1, fi[2] - 1,
                                     fi[3] - 1) * np.pi / 180.
            dihedrals.append(
                [dihed, [fi[0] - 1, fi[1] - 1, fi[2] - 1, fi[3] - 1]])
    for ci in change:
        if len(ci) == 3:
            # careful: atom indices in the fix lists start at 1
            bondlength = ci[2]
            bonds.append([bondlength, [ci[0] - 1, ci[1] - 1]])
        if len(ci) == 4:
            # careful: atom indices in the fix lists start at 1
            angle = ci[3] * np.pi / 180.
            angles.append([angle, [ci[0] - 1, ci[1] - 1, ci[2] - 1]])
        if len(ci) == 5:
            # careful: atom indices in the fix lists start at 1
            dihed = ci[4] * np.pi / 180.
            dihedrals.append(
                [dihed, [ci[0] - 1, ci[1] - 1, ci[2] - 1, ci[3] - 1]])

    return bonds, angles, dihedrals


def read_hess(job, natom):
    """
    Read the hessian of a QChem chk file
    """

    # initialize Hessian
    hess = np.zeros((3 * natom, 3 * natom))

    fchk = str(job) + '.fchk'
    chk = str(job) + '.chk'
    if os.path.exists(chk):
        # create the fchk file using formchk
        os.system('formchk ' + job + '.chk > /dev/null')

    with open(fchk) as f:
        lines = f.read().split('\n')

    nvals = 3 * natom * (3 * natom + 1) / 2

    for index, line in enumerate(reversed(lines)):
        if re.search('Cartesian Force Constants', line) != None:
            hess_flat = []
            n = 0
            while len(hess_flat) < nvals:
                hess_flat.extend(
                    [float(val) for val in lines[-index + n].split()])
                n += 1
            n = 0
            for i in range(3 * natom):
                for j in range(i + 1):
                    hess[i][j] = hess_flat[n]
                    hess[j][i] = hess_flat[n]
                    n += 1
            break
    return hess


def read_imag_mode(job, natom):
    """
    Read the imaginary normal mode displacements from a log file.
    Only for saddle points! It will read the firs normal mode
    for a well, but that's not very useful.
    """
    joblog = '{}_freq.out'.format(job)
    do_read = False
    nmode = []
    with open(joblog) as f:
        for line in f:
            if '               X' in line:
                do_read = True
            elif not do_read:
                continue
            elif do_read:
                if 'TransDip' in line:
                    break
                nmode.append([float(val) for val in line.split()[1:4]])

    return np.array(nmode)


# def read_all_irc_geoms(outfile):
#     """
#     Read the IRC geometries from a QChem 16 log file.
#     Used in sampler code.
#     """
#
#     with open(outfile) as f:
#         lines = f.readlines()
#
#     start = True
#     all_geoms = None
#     atom = None
#     for index, line in enumerate(lines):
#         if 'Charge = ' in line:
#             charge = line.split()[2]
#             mult = line.split()[5]
#         if 'CURRENT STRUCTURE' in line:
#             geom = np.array([])
#             atom = np.array([])
#             natom = 0
#             while True:
#                 current_line = lines[index + 6 + natom]
#                 if '-------' in current_line:
#                     geom = np.reshape(geom, (-1, 3))
#                     if start:
#                         all_geoms = np.array([copy.deepcopy(geom)])
#                         start = False
#                     else:
#                         all_geoms = np.vstack((all_geoms, geom[None]))
#                     break
#                 atom = np.append(atom, int(current_line.split()[1]))
#                 g = np.array(current_line.split()[2:5]).astype(float)
#                 geom = np.append(geom, g)
#                 natom += 1
#     return atom, all_geoms, charge, mult


def write_constraints(inp_file):
    pass
