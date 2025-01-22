import os
import re
import numpy as np
import copy
from kinbot import constants

"""
Functions to read Gaussian output files.
"""


def read_energy(outfile, from_line=0):
    """
    Read the last SCF Done line.
    """

    with open(outfile) as f:
        lines = f.readlines()

    if from_line:
        lines = lines[:-from_line]

    energy = np.nan
    for line in reversed(lines):
        if 'SCF Done' in line:
            energy = float(line.split()[4]) / constants.EVtoHARTREE
            break

    return energy

def read_all_energies(outfile, from_line=0):
    """
    Read all SCF Done lines.
    """

    with open(outfile) as f:
        lines = f.readlines()

    if from_line:
        lines = lines[:-from_line]

    energy = []
    for line in reversed(lines): 
        if 'SCF Done' in line:
            energy.append(float(line.split()[4]) / constants.EVtoHARTREE)

    return energy

def read_geom(outfile, mol, max2frag=False, charge=None, mult=None, from_line=0):
    """
    Read the final geometry from a Gaussian file.
    max2frag: if set to true, the last geometry is taken where
              the structure has max 2 fragments
              This is to avoid "exploding" IRC_prod trajectories.
    """

    with open(outfile) as f:
        lines = f.readlines()

    start = from_line  # start reading here
    while 1:
        geom = np.zeros((len(mol), 3))
        if start == 0:
            data = lines
        else:
            data = lines[:-start]
        for index, line in enumerate(reversed(data)):
            if 'Input orientation:' in line:
                for n in range(len(mol)):
                    geom[n][0:3] = np.array(data[-index+4+n].split()[3:6]).astype(float)
                start += index  # mark end for next loop if there is
                break

        if not max2frag:
            break
        else:
            from kinbot.stationary_pt import StationaryPoint
            temp = StationaryPoint('dummy',
                                   charge,
                                   mult,
                                   atom=[i for i in mol.symbols],
                                   geom=geom)
            temp.characterize()
            fragments, maps = temp.start_multi_molecular()
            if len(fragments) <= 2:
                break

    return geom

def read_all_geoms(outfile, mol, charge=None, mult=None, from_line=0):
    """
    Read all geometries from a Gaussian file.
    """

    with open(outfile) as f:
        lines = f.readlines()

    geoms = []
    geom = np.zeros((len(mol), 3))
    if not from_line:
        data = lines
    else:
        data = lines[:-from_line]
    for index, line in enumerate(reversed(data)):
        if 'Input orientation:' in line:
            for n in range(len(mol)):
                geom[n][0:3] = np.array(data[-index+4+n].split()[3:6]).astype(float)
            geoms.append(copy.deepcopy(geom)) 

    return geoms

def read_zpe(outfile):
    """
    Read the zpe
    """

    with open(outfile) as f:
        lines = f.readlines()

    zpe = np.nan
    for line in reversed(lines):
        if 'Zero-point correction=' in line:
            zpe = float(line.split()[2])

    return zpe


def read_freq(outfile, atom):
    """
    Read the frequencies
    """

    with open(outfile) as f:
        lines = f.readlines()

    natom = len(atom)

    freq = []
    if natom == 1:
        return freq
    else:
        for line in lines:
            if 'Frequencies' in line:
                if natom == 2:
                    freq.append(np.array(line.split()[2]).astype(float))
                    break
                else:
                    f = np.array(line.split()[2:5]).astype(float)
                    freq.extend(f)

    return freq


def read_convergence(outfile):
    """
    Check for the four YES.
    0: did not converge
    1: forces and displacements converged
    2: forces converged
    """

    with open(outfile) as f:
        lines = f.readlines()
        
    for n, line in enumerate(lines):
        if 'Item               Value     Threshold  Converged?' in line:
            if 'YES' in lines[n+1]:
                if 'YES' in lines[n+2]:
                    if 'YES' in lines[n+3]:
                        if 'YES' in lines[n+4]:
                            return 1
                    else:
                        return 2

    return 0  # will look through the whole file

def read_converged_geom_energy(outfile, mol):
    """
    Return the lowest energy geometry with converged forces during an optimization, with its energy
    """

    with open(outfile) as f:
        lines = f.readlines()
    
    msg1 = "Lowest energy point so far.  Saving SCF results."
    msg2 = "Optimization completed"
    

    for n, line in enumerate(reversed(lines)):
        if 'Item               Value     Threshold  Converged?' in line:
            if 'YES' in lines[len(lines)-n]:
                if 'YES' in lines[len(lines)-n+1]:
                    forces_converged = True
                else:
                    forces_converged = False
            else:
                forces_converged = False

            if forces_converged:
                if msg1 in lines[len(lines)-n+5] or\
                msg1 in lines[len(lines)-n+6] or\
                msg2 in lines[len(lines)-n+5] or\
                msg2 in lines[len(lines)-n+6] or\
                read_convergence(outfile) == 1:
                    geom = read_geom(outfile, mol, from_line=n)
                    e = read_energy(outfile, from_line=n)

                    return e, geom
            else:
                raise Exception('Forces are not converged for the lowest energy geometry.')
                
def read_lowest_geom_energy(outfile, mol):
    """
    Return the lowest energy geometry during an optimization, with its energy
    """

    with open(outfile) as f:
        lines = f.readlines()
    
    msg = "Lowest energy point so far.  Saving SCF results."

    geom = mol.positions

    for n, line in enumerate(reversed(lines)):
        if 'Item               Value     Threshold  Converged?' in line:
            if msg in lines[len(lines)-n+5] or msg in lines[len(lines)-n+6]:
                geom = read_geom(outfile, mol, from_line=n)
                e = read_energy(outfile, from_line=n)

                return e, geom

def constraint(mol, fix, change):
    """
    Convert constraints into ASE constraints.
    Angles are in degrees according to angles_deg and dihedrals_deg.
    Careful: atom indices in the fix and change lists start at 1
    """

    bonds = []
    angles = []
    dihedrals = []
    for fi in fix:
        if len(fi) == 2:
            bondlength = mol.get_distance(fi[0]-1, fi[1]-1)
            bonds.append([bondlength, [fi[0]-1, fi[1]-1]])
        if len(fi) == 3:
            angle = mol.get_angle(fi[0]-1, fi[1]-1, fi[2]-1)
            angles.append([angle, [fi[0]-1, fi[1]-1, fi[2]-1]])
        if len(fi) == 4:
            dihed = mol.get_dihedral(fi[0]-1, fi[1]-1, fi[2]-1, fi[3]-1)
            dihedrals.append([dihed, [fi[0]-1, fi[1]-1, fi[2]-1, fi[3]-1]])
    for ci in change:
        if len(ci) == 3:
            bondlength = ci[2]
            bonds.append([bondlength, [ci[0]-1, ci[1]-1]])
        if len(ci) == 4:
            angle = ci[3]
            angles.append([angle, [ci[0]-1, ci[1]-1, ci[2]-1]])
        if len(ci) == 5:
            dihed = ci[4]
            dihedrals.append([dihed, [ci[0]-1, ci[1]-1, ci[2]-1, ci[3]-1]])

    return bonds, angles, dihedrals


def read_hess(job, natom):
    """
    Read the hessian of a gaussian chk file
    """

    # initialize Hessian
    hess = np.zeros((3*natom, 3*natom))
    
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
                hess_flat.extend([float(val) for val in lines[-index + n].split()])
                n += 1
            n = 0
            for i in range(3*natom):
                for j in range(i+1):
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

    nmode = np.zeros([natom, 3])
    joblog = '{}.log'.format(job)
    with open(joblog) as f:
        lines = f.read().split('\n')
    
    for l, line in enumerate(lines):
        if line[:10] == '  Atom  AN':
            for n in range(natom):
                mm = lines[l + n + 1].split() 
                nmode[n][0] = float(mm[2])
                nmode[n][1] = float(mm[3])
                nmode[n][2] = float(mm[4])
            break

    return nmode


def read_all_irc_geoms(outfile):
    """
    Read the IRC geometries from a Gaussian 16 log file.
    Used in sampler code.
    """

    with open(outfile) as f:
        lines = f.readlines()

    start = True
    all_geoms = None
    atom = None
    for index, line in enumerate(lines):
        if 'Charge = ' in line:
            charge = line.split()[2]
            mult = line.split()[5]
        if 'CURRENT STRUCTURE' in line:
            geom = np.array([])
            atom = np.array([])
            natom = 0
            while True:
                current_line = lines[index+6+natom]
                if '-------' in current_line:
                    geom = np.reshape(geom, (-1, 3))
                    if start:
                        all_geoms = np.array([copy.deepcopy(geom)])
                        start = False
                    else:
                        all_geoms = np.vstack((all_geoms, geom[None]))
                    break
                atom = np.append(atom, int(current_line.split()[1]))
                g = np.array(current_line.split()[2:5]).astype(float)
                geom = np.append(geom, g)
                natom += 1
    return atom, all_geoms, charge, mult


def correct_kwargs(outfile, kwargs):
    """Function to correct errors by modifying the kwargs of the ase calculator.

    @param outfile: Output file of gaussian (.log file)
    @param kwargs: Original keyword arguments.
    @return kwargs: New keyword arguments.
    """
    if 'opt' not in kwargs:
        return kwargs

    from kinbot.utils import tail
    outf_end = tail(outfile, 10)
    # Use cartesian coordinates when internal ones fail.
    if 'Error in internal coordinate system' in outf_end \
            and 'cartesian' not in kwargs['opt']:
        kwargs['opt'] += ', cartesian'
    elif 'Error termination request processed by link 9999.' in outf_end:
        kwargs['opt'] = kwargs['opt'].replace('CalcFC', 'CalcAll')
        kwargs.pop('freq', None)
        if 'cartesian' not in kwargs['opt']:
            kwargs['opt'] += ',Cartesian'

    return kwargs
