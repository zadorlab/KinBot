"""Check rotor projection against finite physical rotations, without QC jobs."""

from copy import deepcopy
from types import SimpleNamespace
import unittest

from ase.build import molecule
import numpy as np
from scipy.linalg import null_space
from scipy.spatial.transform import Rotation

from kinbot import constants, frequencies


# ASE reference geometries; bonds and the rotated top are specified explicitly
# so the reference motion does not depend on KinBot's rotor enumeration.
ROTORS = {
    'methanol_co': ('CH3OH', [(0, 1), (0, 2), (0, 4), (0, 5), (1, 3)],
                    [2, 0, 1, 3], [0, 2, 4, 5]),
    'ethanol_cc': ('CH3CH2OH', [(0, 1), (0, 6), (0, 7), (0, 8),
                               (1, 2), (1, 4), (1, 5), (2, 3)],
                   [6, 0, 1, 2], [0, 6, 7, 8]),
    'ethanol_co': ('CH3CH2OH', [(0, 1), (0, 6), (0, 7), (0, 8),
                               (1, 2), (1, 4), (1, 5), (2, 3)],
                   [0, 1, 2, 3], [0, 1, 4, 5, 6, 7, 8]),
    'ethane_cc': ('C2H6', [(0, 1), (0, 2), (0, 3), (0, 4),
                          (1, 5), (1, 6), (1, 7)],
                  [2, 0, 1, 5], [0, 2, 3, 4]),
}


def rotor_fixture(name, saddle=False):
    """Construct a Hessian whose torsional eigenvector is known independently."""
    formula, bonds, rotor, top = ROTORS[name]
    atoms = molecule(formula)
    bond = np.zeros((len(atoms), len(atoms)))
    for i, j in bonds:
        bond[i, j] = bond[j, i] = 1
    species = SimpleNamespace(name=name, atom=atoms.get_chemical_symbols(),
                              natom=len(atoms), geom=atoms.positions.copy(),
                              bond=bond, dihed=[rotor.copy()], wellorts=int(saddle))
    mass = np.array([constants.exact_mass[a] for a in species.atom])
    xyz = species.geom - np.average(species.geom, axis=0, weights=mass)
    roots = np.sqrt(mass[:, None])
    rigid = [np.tile(axis, (species.natom, 1)) * roots for axis in np.eye(3)]
    rigid += [np.cross(xyz, axis) * roots for axis in np.eye(3)]
    vibrational = null_space(np.array(rigid).reshape(6, -1))

    # Differentiate actual rotations of one whole top about the Cartesian bond.
    # This supplies an independent reference, rather than reimplementing the
    # analytical cross product used by get_frequencies.
    axis = xyz[rotor[2]] - xyz[rotor[1]]
    axis /= np.linalg.norm(axis)
    displaced = []
    for angle in (-1.e-5, 1.e-5):
        positions = xyz.copy()
        positions[top] = Rotation.from_rotvec(angle * axis).apply(
            xyz[top] - xyz[rotor[1]]) + xyz[rotor[1]]
        displaced.append(positions)
    derivative = ((displaced[1] - displaced[0]) / 2.e-5 * roots).ravel()
    torsion = vibrational @ (vibrational.T @ derivative)
    torsion /= np.linalg.norm(torsion)
    other = vibrational @ null_space((torsion @ vibrational)[None, :])
    basis = np.column_stack([torsion, other])
    eigenvalues = np.r_[.003, np.arange(1, other.shape[1] + 1) * .02]
    if saddle:
        eigenvalues[1] = -.01  # A reaction coordinate orthogonal to the torsion.
    weighted = (basis * eigenvalues) @ basis.T
    masses = np.repeat(mass, 3)
    hessian = weighted * np.sqrt(np.outer(masses, masses))
    return species, hessian, eigenvalues


class TestPhysicalRotorProjection(unittest.TestCase):
    def test_only_the_known_torsional_mode_is_removed(self):
        for name in ROTORS:
            for saddle in (False, True):
                with self.subTest(rotor=name, saddle=saddle):
                    species, hessian, values = rotor_fixture(name, saddle)
                    raw, reduced = frequencies.get_frequencies(species, hessian, species.geom)
                    expected = sorted(map(frequencies.convert_to_wavenumbers, values))
                    expected_reduced = sorted(map(frequencies.convert_to_wavenumbers, values[1:]))
                    np.testing.assert_allclose(raw, expected, rtol=1.e-9, atol=1.e-6)
                    np.testing.assert_allclose(reduced, expected_reduced, rtol=1.e-9, atol=1.e-6)
                    self.assertEqual(sum(f < 0 for f in reduced), int(saddle))

    def test_rigid_motion_and_axis_reversal_preserve_both_spectra(self):
        rotation = Rotation.from_rotvec([.3, -.8, .4]).as_matrix()
        for name in ROTORS:
            with self.subTest(rotor=name):
                species, hessian, _ = rotor_fixture(name)
                expected = frequencies.get_frequencies(species, hessian, species.geom)
                transform = np.kron(np.eye(species.natom), rotation)
                species.dihed = [species.dihed[0][::-1]]
                actual = frequencies.get_frequencies(
                    species, transform @ hessian @ transform.T,
                    species.geom @ rotation.T + [20., -13., 4.])
                for reference, result in zip(expected, actual):
                    np.testing.assert_allclose(result, reference, rtol=1.e-9, atol=1.e-6)

    def test_atom_permutation_and_preweighted_hessian_preserve_results(self):
        for name in ROTORS:
            with self.subTest(rotor=name):
                species, hessian, _ = rotor_fixture(name)
                expected = frequencies.get_frequencies(species, hessian, species.geom)
                order = np.random.default_rng(7).permutation(species.natom)
                inverse = np.argsort(order)
                permuted = deepcopy(species)
                permuted.atom = np.asarray(species.atom)[order]
                permuted.bond = species.bond[np.ix_(order, order)]
                permuted.dihed = [inverse[species.dihed[0]].tolist()]
                indices = (order[:, None] * 3 + np.arange(3)).ravel()
                masses = np.repeat([constants.exact_mass[a] for a in permuted.atom], 3)
                weighted = hessian[np.ix_(indices, indices)] / np.sqrt(np.outer(masses, masses))
                actual = frequencies.get_frequencies(
                    permuted, weighted, species.geom[order], massweighted=True)
                for reference, result in zip(expected, actual):
                    np.testing.assert_allclose(result, reference, rtol=1.e-9, atol=1.e-6)

    def test_no_rotor_preserves_the_complete_harmonic_spectrum(self):
        species, hessian, values = rotor_fixture('ethanol_co')
        species.dihed = []
        expected = sorted(map(frequencies.convert_to_wavenumbers, values))
        for result in frequencies.get_frequencies(species, hessian, species.geom):
            np.testing.assert_allclose(result, expected, rtol=1.e-9, atol=1.e-6)


if __name__ == '__main__':
    unittest.main()
