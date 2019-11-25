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

import os
import re
import numpy as np

"""
Functions to read quantum chemistry output files.
"""

def read_geom(outfile, mol, dummy):
    """
    Read the final geometry from a Gaussian file.
    """
    with open(outfile) as f:
        lines = f.readlines()

    geom = np.zeros((len(mol),3))
    for index, line in enumerate(reversed(lines)):
        if 'Input orientation:' in line:
            for n in range(len(mol)):
                geom[n][0:3] = np.array(lines[-index+4+n].split()[3:6]).astype(float)
            break

    for i,d in enumerate(dummy):
        geom[-(i+1)][0:3] = d[0:3]

    return geom


def read_zpe(outfile):
    """
    Read the zpe
    """

    with open(outfile) as f:
        lines = f.readlines()

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

    natom = len([at for at in atom if at !='X']) #filter out the dummy atoms
    if natom == 1:
        freq = []
    else:
        freq = []
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


def constraint(mol, fix, change):
    """
    Convert constraints into PCBFGS constraints.
    """

    bonds = []
    angles = []
    dihedrals = []
    for fi in fix:
        if len(fi) == 2:
            #careful: atom indices in the fix lists start at 1
            bondlength = mol.get_distance(fi[0]-1, fi[1]-1)
            bonds.append([bondlength,[fi[0]-1, fi[1]-1]])
        if len(fi) == 3:
            #careful: atom indices in the fix lists start at 1
            angle = mol.get_angle(fi[0]-1, fi[1]-1, fi[2]-1) * np.pi / 180.
            angles.append([angle,[fi[0]-1, fi[1]-1, fi[2]-1]])
        if len(fi) == 4:
            #careful: atom indices in the fix lists start at 1
            dihed = mol.get_dihedral(fi[0]-1, fi[1]-1, fi[2]-1, fi[3]-1) * np.pi / 180.
            dihedrals.append([dihed,[fi[0]-1, fi[1]-1, fi[2]-1, fi[3]-1]])
    for ci in change:
        if len(ci) == 3:
            #careful: atom indices in the fix lists start at 1
            bondlength = ci[2]
            bonds.append([bondlength,[ci[0]-1, ci[1]-1]])
        if len(ci) == 4:
            #careful: atom indices in the fix lists start at 1
            angle = ci[3] * np.pi / 180.
            angles.append([angle, [ci[0]-1, ci[1]-1, ci[2]-1]])
        if len(ci) == 5:
            #careful: atom indices in the fix lists start at 1
            dihed = ci[4] * np.pi / 180.
            dihedrals.append([dihed,[ci[0]-1, ci[1]-1, ci[2]-1, ci[3]-1]])

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
    #create the fchk file using formchk
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
                nmode[n][0]= float(mm[2])
                nmode[n][1]= float(mm[3])
                nmode[n][2]= float(mm[4])
            break

    return(nmode)


