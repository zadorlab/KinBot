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
import numpy as np

from kinbot import constants
from kinbot import geometry


def get_frequencies(species, hess, geom):
    """"
    Calculates three sets of frequencies:

    1: all the frequencies including translations and rotations
    2: frequencies when translation and external rotations are projected out
        these should be identical to the frequencies supplied by Gaussian
    3: frequencies when internal rotations are also projected out

    The units of the hessian should be: Hartree/Bohr^2
    """
    atom = species.atom
    natom = species.natom

    masses = []
    for at in atom:
        masses += [constants.exact_mass[at]] * 3
    masses = np.array(masses)

    # Translate molecule's center of mass to (0, 0, 0)
    geom = geom - geometry.get_center_of_mass(geom, atom)

    # Mass-weight the hessian
    hess /= np.sqrt(np.outer(masses, masses))

    # STEP 1: calculate the initial frequencies
    all_eigvals, all_eigvecs = np.linalg.eig(hess)

    all_modes = all_eigvecs.T
    all_modes /= np.sqrt(masses[np.newaxis, :])
    for mode in all_modes:
        mode /= np.linalg.norm(mode)

    # STEP 2: project out translation and rotation
    # Build set of translation vectors to project out
    tvecs = np.zeros((3, 3 * natom))
    for i in range(3):
        ar = [np.array([1., 0., 0.] * natom) * np.sqrt(masses)]
        tvecs[i] = np.roll(ar, i)
        tvecs[i] /= np.linalg.norm(tvecs[i])

    # Start to build rotation vectors
    I = geometry.get_moments_of_inertia(geom, atom)[1]
    X = I.T
    P = np.dot(geom, X)
    D = np.zeros((natom, 3, 3))
    for i, Pi in enumerate(P):
        D[i] = np.cross(Pi, I.T) * np.sqrt(constants.exact_mass[atom[i]])
    D4 = D[:, :, 0].ravel()
    D5 = D[:, :, 1].ravel()
    D6 = D[:, :, 2].ravel()

    # Small rotation vector magnitudes mean it's not a real rotation
    # (this can happen for linear molecules), so don't include those
    rvecs = []
    D4_norm = np.linalg.norm(D4)
    if D4_norm > 1e-5:
        rvecs.append(D4 / D4_norm)

    D5_norm = np.linalg.norm(D5)
    if D5_norm > 1e-5:
        rvecs.append(D5 / D5_norm)

    D6_norm = np.linalg.norm(D6)
    if D6_norm > 1e-5:
        rvecs.append(D6 / D6_norm)
    rvecs = np.array(rvecs)

    nvecs = 3 * natom - 3 - len(rvecs)
    vecs = np.zeros((nvecs, 3 * natom))
    n = 0

    # Use Gram-Schmidt orthonormalization to build new projected basis.
    # The details of this basis don't matter, because we're going to
    # convert the eigenvectors of the projected Hessian back to Cartesian
    # coordinates.
    while n < nvecs:
        vec = np.random.random(3 * natom)
        vec /= np.linalg.norm(vec)
        for tvec in tvecs:
            vec -= np.dot(vec, tvec) * tvec
        for rvec in rvecs:
            vec -= np.dot(vec, rvec) * rvec
        for i in range(n):
            vec -= np.dot(vec, vecs[i]) * vecs[i]

        if np.linalg.norm(vec) > 1e-4:
            vecs[n] = vec / np.linalg.norm(vec)
            n += 1

    # Projected Hessian
    Phess = np.dot(np.dot(vecs, hess), vecs.T)
    eigvals, eigvecs = np.linalg.eig(Phess)

    modes = np.dot(eigvecs.T, vecs)
    modes /= np.sqrt(masses[np.newaxis, :])
    for mode in modes:
        mode /= np.linalg.norm(mode)

    freqs = [convert_to_wavenumbers(ei) for ei in sorted(eigvals)]

    # STEP 3: project out internal rotations

    # Build set of internal rotation vectors to project out
    R = []
    for rot in species.dihed:
        # mass weight the cartesian coordinates
        mgeom = np.zeros((natom, 3))
        for i in range(natom):
            mgeom[i][0:3] += geom[i] * np.sqrt(constants.exact_mass[atom[i]])
        # partition the molecule in two parts divided by the rotor bond
        Ri = np.zeros(3*natom)
        l1, l2 = partition(species, rot, natom)

        axis = mgeom[rot[1]] - mgeom[rot[2]]
        axis = axis / np.linalg.norm(axis)
        for at in range(natom):
            vect = mgeom[at] - mgeom[rot[2]]
            proj = np.dot(vect, axis)*axis/np.linalg.norm(axis)**2
            # vector perpendicular to the rotational axis through the atom
            per = vect - proj
            dist = np.linalg.norm(per)
            if dist > 1e-6:
                sign = 1
                if at in l2:
                    sign = -1
                rot_vect = sign * np.cross(per, axis)
                Ri[3*at:3*at+3] = rot_vect

        # project the translational, external rotational and previous
        # internal rotations out of the current vector
        Ri = Ri / np.linalg.norm(Ri)
        for tvec in tvecs:
            Ri -= np.dot(Ri, tvec) * tvec
        for rvec in rvecs:
            Ri -= np.dot(Ri, rvec) * rvec
        for Rj in R:
            Ri -= np.dot(Ri, Rj) * Rj
        Ri = Ri / np.linalg.norm(Ri)

        R.append(Ri / np.linalg.norm(Ri))

    nvecs = 3 * natom - 3 - len(rvecs) - len(R)
    vecs = np.zeros((nvecs, 3 * natom))
    n = 0

    # Use Gram-Schmidt orthonormalization to build new projected basis.
    # The details of this basis don't matter, because we're going to
    # convert the eigenvectors of the projected Hessian back to Cartesian
    # coordinates.
    while n < nvecs:
        vec = np.random.random(3 * natom)
        vec /= np.linalg.norm(vec)
        for tvec in tvecs:
            vec -= np.dot(vec, tvec) * tvec
        for rvec in rvecs:
            vec -= np.dot(vec, rvec) * rvec
        for Rvec in R:
            vec -= np.dot(vec, Rvec) * Rvec
        for i in range(n):
            vec -= np.dot(vec, vecs[i]) * vecs[i]

        if np.linalg.norm(vec) > 1e-4:
            vecs[n] = vec / np.linalg.norm(vec)
            n += 1

    # Projected Hessian
    Phess = np.dot(np.dot(vecs, hess), vecs.T)
    eigvals, eigvecs = np.linalg.eig(Phess)

    modes = np.dot(eigvecs.T, vecs)
    modes /= np.sqrt(masses[np.newaxis, :])
    for mode in modes:
        mode /= np.linalg.norm(mode)

    reduced_freqs = [convert_to_wavenumbers(ei) for ei in sorted(eigvals)]

    return freqs, reduced_freqs


def convert_to_wavenumbers(val):
    """
    Convert an eigenvalue in Hartree * bohr^-2 * amu^-1
    to a frequency in wavenumbers
    """
    denominator = 2. * np.pi * constants.SPEEDofLIGHT * constants.BOHRtoCM
    if val < 0.:
        fr = -np.sqrt(-val / constants.MEtoAMU) / denominator
    else:
        fr = np.sqrt(val / constants.MEtoAMU) / denominator

    return fr


def partition(species, rotor, natom):
    l1 = [rotor[1]]
    forbidden = [rotor[2]]
    visited = [rotor[1], rotor[2]]
    get_neighbors(rotor[1], visited, forbidden, l1, species.bond, natom)
    return l1, [x for x in range(natom) if x not in l1]


def get_neighbors(ati, visited, forbidden, division, bond, natom):
    for atj in range(natom):
        if atj not in visited and atj not in forbidden:
            if bond[atj, ati] > 0:
                division.append(atj)
                visited.append(atj)
                get_neighbors(atj, visited, forbidden, division, bond, natom)
