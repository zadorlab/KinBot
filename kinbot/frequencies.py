import os
import shutil
import logging
import numpy as np
from ase.data import atomic_numbers, covalent_radii
from ase import Atoms
from ase.vibrations import Vibrations

from kinbot import constants
from kinbot import geometry
from kinbot.stationary_pt import StationaryPoint
from kinbot.constants import EVtoHARTREE

logger = logging.getLogger('KinBot')

# Linearity is decided in get_frequencies by three gates (see assess_linearity):
#  1. every bond angle lies within LINEAR_ANGLE_TOLERANCE degrees of 180, which
#     is what an optimiser without symmetry constraints leaves behind;
#  2. the curvature of the Hessian along the rigid rotation about the
#     molecular axis, expressed as a wavenumber, exceeds LINEAR_AXIS_MODE_MIN.
#     For a molecule whose minimum is linear that curvature is the second
#     component of the degenerate bend (hundreds of cm-1); for a molecule whose
#     minimum is genuinely bent it is a free rotation (zero);
#  3. it also exceeds LINEAR_AXIS_MODE_NOISE_FACTOR times the residual
#     curvature of the other five rigid-body motions, i.e. the noise of that
#     particular Hessian.
# A geometry is never modified; a linearised copy is used for the rigid-body
# vectors here and for the MESS geometry block, so MESS's own linearity test
# reaches the same conclusion.
LINEAR_ANGLE_TOLERANCE = 2.0
# MESS decides the rotational dimension itself from the geometry it is given:
# linear when I_min/I_mid < 1e-5 (libmess model.cc, RigidRotor). KinBot's
# decision above is stricter about a bent minimum, so a species can be
# non-linear for KinBot yet linear for MESS; MESS.rotor_core_line handles it.
MESS_LINEAR_MOMENT_RATIO = 1.e-5
# h / (8 pi^2 c) in cm-1 amu Angstrom^2: B = ROTATIONAL_CONSTANT_FACTOR / I
ROTATIONAL_CONSTANT_FACTOR = 16.857629
LINEAR_AXIS_MODE_MIN = 50.
LINEAR_AXIS_MODE_NOISE_FACTOR = 3.


def max_bend_deviation(geom, bond):
    """Largest deviation (degrees) of any bond angle a-b-c from 180."""
    geom = np.asarray(geom, dtype=float)
    natom = len(geom)
    worst = 0.
    for b in range(natom):
        neighbours = [i for i in range(natom) if i != b and bond[b][i] > 0]
        for k, a in enumerate(neighbours):
            for c in neighbours[k + 1:]:
                u = geom[a] - geom[b]
                v = geom[c] - geom[b]
                cosine = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
                worst = max(worst, 180. - np.degrees(np.arccos(np.clip(cosine, -1., 1.))))
    return worst


def moment_ratio(geom, atom):
    """I_min / I_mid of the principal moments of inertia; MESS's linearity measure."""
    moments = geometry.get_moments_of_inertia(np.asarray(geom, dtype=float), atom)[0]
    return moments[0] / moments[1] if moments[1] > 0. else 0.


def rotational_constants(geom, atom):
    """Rotational constants in cm-1 from the principal moments (amu Angstrom^2)."""
    moments = geometry.get_moments_of_inertia(np.asarray(geom, dtype=float), atom)[0]
    return [ROTATIONAL_CONSTANT_FACTOR / moment for moment in moments]


def linearize(geom, atom):
    """Copy of geom with every atom projected onto the principal axis with the
    smallest moment of inertia. The input array is not modified."""
    geom = np.asarray(geom, dtype=float)
    com = geometry.get_center_of_mass(geom, atom)
    axis = geometry.get_moments_of_inertia(geom - com, atom)[1][0]
    axis = axis / np.linalg.norm(axis)
    return com + np.outer(np.dot(geom - com, axis), axis)


def rigid_body_vectors(geom, atom):
    """Mass-weighted translation vectors (normalised) and rotation vectors
    (unnormalised; their norms are the square roots of the principal moments)
    for a geometry centred on its centre of mass."""
    natom = len(atom)
    masses = np.repeat([constants.exact_mass[at] for at in atom], 3)
    tvecs = np.zeros((3, 3 * natom))
    for i in range(3):
        ar = [np.array([1., 0., 0.] * natom) * np.sqrt(masses)]
        tvecs[i] = np.roll(ar, i)
        tvecs[i] /= np.linalg.norm(tvecs[i])
    I = geometry.get_moments_of_inertia(geom, atom)[1]
    P = np.dot(geom, I.T)
    D = np.zeros((natom, 3, 3))
    for i, Pi in enumerate(P):
        D[i] = np.cross(Pi, I.T) * np.sqrt(constants.exact_mass[atom[i]])
    return tvecs, [D[:, :, k].ravel() for k in range(3)]


def assess_linearity(species, hess_mw, geom):
    """Three-gate linearity decision; see the module constants.

    geom must be centred on the centre of mass and hess_mw mass-weighted and
    symmetric, as prepared in get_frequencies. Returns a dict with 'linear'
    and the diagnostics behind the decision.
    """
    natom = len(species.atom)
    result = {'linear': False, 'angle_deviation': 0., 'axis_mode': None,
              'noise_floor': None, 'reason': ''}
    if natom < 3:
        result.update(linear=(natom == 2), reason='diatomic' if natom == 2 else 'atom')
        return result
    deviation = max_bend_deviation(geom, species.bond)
    result['angle_deviation'] = deviation
    if deviation > LINEAR_ANGLE_TOLERANCE:
        result['reason'] = f'bond angles deviate up to {deviation:.1f} degrees from 180'
        return result
    tvecs, rotations = rigid_body_vectors(geom, species.atom)
    norms = [np.linalg.norm(r) for r in rotations]
    k_axis = int(np.argmin(norms))
    if norms[k_axis] == 0.:
        result.update(linear=True, axis_mode=np.inf, noise_floor=0.,
                      reason='exactly linear geometry')
        return result
    others = list(tvecs) + [rotations[k] / norms[k] for k in range(3) if k != k_axis]
    v = rotations[k_axis] / norms[k_axis]
    for b in others:
        v = v - np.dot(v, b) * b
    v /= np.linalg.norm(v)
    axis_mode = convert_to_wavenumbers(v @ hess_mw @ v)
    noise = max(abs(convert_to_wavenumbers(b @ hess_mw @ b)) for b in others)
    result.update(axis_mode=axis_mode, noise_floor=noise)
    if axis_mode > max(LINEAR_AXIS_MODE_MIN, LINEAR_AXIS_MODE_NOISE_FACTOR * noise):
        result.update(linear=True, reason=(
            f'curvature along the axis rotation is a vibration ({axis_mode:.0f} cm-1, '
            f'rigid-body noise {noise:.0f} cm-1)'))
    else:
        result['reason'] = (f'curvature along the axis rotation is {axis_mode:.0f} cm-1 '
                            f'(rigid-body noise {noise:.0f} cm-1): a free rotation, so the '
                            'minimum is bent')
    return result


def thermochemical_frequencies(raw, wellorts=0, imagfreq_threshold=50.):
    """Copy modes and apply KinBot's accepted small-imaginary correction.

    The first saddle mode remains the reaction coordinate. This prepares
    thermochemical data; it does not validate or change the raw calculation.
    """
    result = [float(value) for value in raw]
    for index in range(int(bool(wellorts)), len(result)):
        if -imagfreq_threshold <= result[index] < 0.:
            result[index] *= -1.
    return result


def get_frequencies(species, hess, geom, checkdist=0, massweighted=False):
    """"Calculates three sets of frequencies:

    1: all the frequencies including translations and rotations
    2: frequencies when translation and external rotations are projected out
        these should be identical to the frequencies supplied by Gaussian
    3: frequencies when internal rotations are also projected out

    The units of the hessian should be: Hartree/Bohr^2

    checkdist: if set to 1, then in the partitioning of the rotating fragments
        only strongly bonded atoms are included
    massweighted: whether the hessian is already mass-weighted or not.
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
    if not massweighted:
        # Cannot use /= on immutable arrays read from db. (Sella)
        hess = hess / np.sqrt(np.outer(masses, masses))
    # The Hessian is symmetric up to numerical noise; enforce it so that the
    # symmetric eigensolver below is exact and never returns complex modes.
    hess = 0.5 * (hess + hess.T)

    # STEP 1: calculate the initial frequencies
    all_eigvals, all_eigvecs = np.linalg.eigh(hess)

    all_modes = all_eigvecs.T
    all_modes /= np.sqrt(masses[np.newaxis, :])
    for mode in all_modes:
        mode /= np.linalg.norm(mode)

    # STEP 2: project out translation and rotation. For a linear molecule the
    # rotation about the molecular axis is part of the degenerate bend and must
    # stay in the spectrum; assess_linearity decides, and a linearised copy of
    # the geometry then makes that rotation vanish exactly.
    verdict = assess_linearity(species, hess, geom)
    if natom >= 3:
        name = getattr(species, 'name', '')
        if verdict['linear']:
            logger.info(f'{name}: treated as a linear rotor, 3N-5 frequencies '
                        f'(max bend deviation {verdict["angle_deviation"]:.2f} deg; '
                        f'{verdict["reason"]}).')
        elif verdict['angle_deviation'] <= 10.:
            logger.warning(f'{name}: geometry is within {verdict["angle_deviation"]:.1f} deg '
                           f'of linear but is treated as non-linear ({verdict["reason"]}). '
                           'If this species is linear, tighten the optimisation.')
    if verdict['linear']:
        geom = linearize(geom, atom)
    tvecs, rotations = rigid_body_vectors(geom, atom)
    norms = [np.linalg.norm(r) for r in rotations]
    rvecs = np.array([r / n for r, n in zip(rotations, norms) if n > 1e-8 * max(norms)])

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
    eigvals, eigvecs = np.linalg.eigh(Phess)

    modes = np.dot(eigvecs.T, vecs)
    modes /= np.sqrt(masses[np.newaxis, :])
    for mode in modes:
        mode /= np.linalg.norm(mode)

    freqs = [convert_to_wavenumbers(ei) for ei in sorted(eigvals)]

    # STEP 3: project out internal rotations

    # Build set of internal rotation vectors to project out
    R = []
    for rotor, rot in enumerate(species.dihed):
        hir = getattr(species, 'hir', None)
        if hir is not None and not hir.is_valid_rotor(rotor):
            continue
        if skip_rotor(species.name, rot) == 1:
            continue
            
        # partition the molecule in two parts divided by the rotor bond
        Ri = np.zeros(3 * natom)
        l1, l2 = partition(species, rot, natom, checkdist=checkdist)
        # Construct the physical rotation in Cartesian coordinates. Weight
        # the displacement afterwards: differences of mass-weighted atom
        # positions do not describe rotation about the actual bond axis.
        axis = geom[rot[1]] - geom[rot[2]]
        axis = axis / np.linalg.norm(axis)
        for at in range(natom):
            vect = geom[at] - geom[rot[2]]
            proj = np.dot(vect, axis)*axis/np.linalg.norm(axis)**2
            # vector perpendicular to the rotational axis through the atom
            per = vect - proj
            dist = np.linalg.norm(per)
            if dist > 1e-6:
                if at in l1:
                    sign = 1
                elif at in l2:
                    sign = -1
                else:
                    sign = 0
                rot_vect = (sign * np.cross(per, axis)
                            * np.sqrt(constants.exact_mass[atom[at]]))
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
    eigvals, eigvecs = np.linalg.eigh(Phess)

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


def partition(species, rotor, natom, checkdist=0):
    l1 = [rotor[1]]
    forbidden = [rotor[2]]
    visited = [rotor[1], rotor[2]]
    get_neighbors(rotor[1], visited, forbidden, l1, species, natom, checkdist)
    if checkdist == 0:
        return l1, [x for x in range(natom) if x not in l1]
    else:
        l2 =[rotor[2]]
        forbidden = [rotor[1]]
        visited = [rotor[1], rotor[2]]
        get_neighbors(rotor[2], visited, forbidden, l2, species, natom, checkdist)
        return l1, l2


def get_neighbors(ati, visited, forbidden, division, species, natom, checkdist):
    for atj in range(natom):
        if atj not in visited and atj not in forbidden:
            if species.bond[atj, ati] > 0:
                if checkdist == 0:
                    division.append(atj)
                    visited.append(atj)
                    get_neighbors(atj, visited, forbidden, division, species, natom, checkdist)
                else:
                    try:
                        cutoff = constants.st_bond[''.join(sorted(species.atom[atj] + species.atom[ati]))]
                    except KeyError:
                        cutoff = 1.2 * (covalent_radii[atomic_numbers[species.atom[ati]]] 
                                        + covalent_radii[atomic_numbers[species.atom[atj]]])
                    if species.dist[atj, ati] < cutoff:
                        division.append(atj)
                        visited.append(atj)
                        get_neighbors(atj, visited, forbidden, division, species, natom, checkdist)


def skip_rotor(name, rot):
    if 'barrierless_saddle' in name:
        return 0
        l0 = name.split('_')
        l = [int(l0[3]) - 1, int(l0[4]) - 1] 
        if any(rot[i:i+2] == l for i in range(3)):
            return 1
        if any(rot[i:i+2] == l[::-1] for i in range(3)):
            return 1
        return 0
    elif 'prod' in name:
        return 1


def calc_vibrations(mol, label):
    # this is a frequency calculator when ASE is used
    mol.calc.label = f'{label}_vib'
    if 'chk' in mol.calc.parameters:
        del mol.calc.parameters['chk']
    # Compute frequencies in a separate temporary directory to avoid 
    # conflicts accessing the cache in parallel calculations.
    if not os.path.isdir(f'{label}_vib'):
        os.mkdir(f'{label}_vib')
    init_dir = os.getcwd()
    os.chdir(f'{label}_vib')
    try:
        if os.path.isdir('vib'):
            shutil.rmtree('vib')
        vib = Vibrations(mol)
        vib.run()
        vib.write_jmol()
        # Use kinbot frequencies to avoid mixing low vib frequencies with 
        # the values associated with external rotations.
        _ = vib.get_frequencies()
        zpe = vib.get_zero_point_energy() * EVtoHARTREE
        hessian = vib.H / 97.17370087
        st_pt = StationaryPoint.from_ase_atoms(mol)
        st_pt.characterize()
        freqs, _ = get_frequencies(st_pt, hessian, st_pt.geom)
        return freqs, zpe, hessian
    except Exception:
        return None, None, None
    finally:
        os.chdir(init_dir)
