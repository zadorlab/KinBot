"""Rotating bonds stay single; torsion reference bonds need not be single."""
from types import SimpleNamespace
import unittest

import numpy as np
from ase import Atoms
from ase.build import molecule
from scipy.linalg import null_space

from kinbot import constants, frequencies, symmetry
from kinbot.hindered_rotors import HIR
from kinbot.stationary_pt import StationaryPoint


def point(atoms, mult=1, order=None):
    if order is not None:
        atoms = atoms[order]
    species = StationaryPoint('rotor_detection', 0, mult,
        atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    species.characterize()
    return species


class TestRotorDetection(unittest.TestCase):
    def test_nitro_nitroso_and_acyl_rotors_in_both_axis_orders(self):
        for name, mult in [('CH3NO2', 1), ('CH3NO', 1), ('CH3CO', 2)]:
            # Detection fixture only: remove one O for the nitroso group.
            # This is not an optimized nitrosomethane geometry.
            atoms = (molecule('CH3NO2')[:-1]
                     if name == 'CH3NO' else molecule(name))
            rng = np.random.default_rng(173)
            orders = [np.arange(len(atoms)), np.arange(len(atoms))[::-1]]
            orders += [rng.permutation(len(atoms)) for _ in range(12)]
            for order in orders:
                with self.subTest(name=name, order=order.tolist()):
                    species = point(atoms, mult, order)
                    self.assertEqual(len(species.dihed), 1)
                    a, b, c, d = species.dihed[0]
                    self.assertEqual(sorted([species.atom[b], species.atom[c]]),
                                     ['C', 'C'] if name == 'CH3CO' else ['C', 'N'])
                    self.assertEqual(sorted([species.bond[a, b], species.bond[c, d]]), [1, 2])
                    self.assertTrue(all(bond[b, c] == 1 for bond in species.bonds))
                    # Finding a HIR rotor must not sample symmetry-equivalent
                    # methyl orientations as extra conformers.
                    self.assertEqual(species.conf_dihed, [])

    def test_redundant_nitro_references_do_not_duplicate_the_axis(self):
        species = point(molecule('CH3NO2'))
        graph, identity = species.bond.copy(), species.chemid
        symmetry.calculate_symmetry(species)
        self.assertEqual(species.sigma_int[0][1], 6)
        species.find_dihedral(findall=1)
        self.assertEqual(len(species.dihed_allrot), 6)  # three H times two O
        self.assertEqual({tuple(r[1:3]) for r in species.dihed_allrot}, {(0, 1)})
        species.find_dihedral()
        self.assertEqual(len(species.dihed), 1)
        np.testing.assert_array_equal(species.bond, graph)
        self.assertEqual(species.chemid, identity)

    def test_existing_rotor_references_remain_unchanged(self):
        expected = {
            'CH3OH': [[2, 0, 1, 3]],
            'CH3CH2OH': [[6, 0, 1, 2], [0, 1, 2, 3]],
            'CH3COCH3': [[3, 1, 2, 4], [2, 1, 3, 5]],
            'CH3CHO': [[2, 1, 3, 4]],
            'C2H6': [[2, 0, 1, 5]],
        }
        for name, rotors in expected.items():
            with self.subTest(name=name):
                self.assertEqual(point(molecule(name)).dihed, rotors)

    def test_new_axis_is_appended_without_renumbering_existing_rotors(self):
        # Fixed nitroethane geometry (RDKit/MMFF, seed 71); no RDKit runtime
        # dependency. Reverse the atoms so the missing C-N axis sorts first.
        atoms = Atoms('CCNOOHHHHH', positions=[
            [1.06018423, -.40985387, .02750706],
            [-.03570052, .62843464, -.08975162],
            [-1.35865392, -.05859812, .04912401],
            [-1.82204219, -.16237203, 1.19162857],
            [-1.88726904, -.47739185, -.98817595],
            [.97245909, -1.17262515, -.75381576],
            [2.04182915, .06402358, -.07034597],
            [1.02482634, -.9196928, .99642097],
            [.02859292, 1.38158493, .70127792],
            [-.02422606, 1.12649068, -1.06386922]])
        species = point(atoms, order=np.arange(len(atoms))[::-1])
        self.assertEqual([r[1:3] for r in species.dihed], [[8, 9], [7, 8]])
        self.assertEqual(species.dihed[0], [0, 8, 9, 2])

    def test_multiple_ring_and_linear_axes_stay_excluded(self):
        for name in ['C2H4', 'C2H2', 'C6H6', 'C3H6_D3h', 'CH3CN']:
            with self.subTest(name=name):
                self.assertEqual(point(molecule(name)).dihed, [])
        species = point(molecule('CH3NO2'))
        other_resonance = species.bond.copy()
        other_resonance[0, 1] = other_resonance[1, 0] = 2
        species.bonds.append(other_resonance)
        species.find_dihedral()
        self.assertEqual(species.dihed, [])

    def test_nitro_hir_geometries_preserve_both_fragments(self):
        species = point(molecule('CH3NO2'))
        self.assertEqual(len(species.dihed), 1)
        captures = []

        def capture(species, xyz, rotor, angle, fix, rigid):
            captures.append(np.asarray(xyz).copy())
            return f'captured_{rotor}_{angle}'

        hir = HIR(species, SimpleNamespace(qc='gauss', qc_hir=capture),
                  {'nrotation': 12, 'plot_hir_profiles': False, 'rotor_0_test': True})
        hir.generate_hir_geoms(species.geom, False)
        self.assertEqual(len(captures), 12)
        groups = frequencies.partition(species, species.dihed[0], species.natom)
        self.assertEqual(sorted(map(len, groups)), [3, 4])
        angles = []
        for xyz in captures:
            self.assertTrue(np.isfinite(xyz).all())
            bonded = np.argwhere(np.triu(species.bond) > 0)
            np.testing.assert_allclose(
                [np.linalg.norm(xyz[a] - xyz[b]) for a, b in bonded],
                [np.linalg.norm(species.geom[a] - species.geom[b]) for a, b in bonded],
                atol=1e-7)
            for group in groups:
                old, new = species.geom[group], xyz[group]
                np.testing.assert_allclose(np.linalg.norm(old[:, None]-old[None, :], axis=-1),
                                           np.linalg.norm(new[:, None]-new[None, :], axis=-1), atol=1e-7)
            angles.append(Atoms(species.atom, positions=xyz).get_dihedral(*species.dihed[0]))
        np.testing.assert_allclose(np.abs(np.diff(np.unwrap(np.radians(angles))))*180/np.pi,
                                   30., atol=1e-7)

    def test_nitro_projection_removes_one_mode(self):
        species = point(molecule('CH3NO2'))
        mass = np.array([constants.exact_mass[a] for a in species.atom])
        xyz = species.geom - np.average(species.geom, axis=0, weights=mass)
        rigid = [np.tile(v, (species.natom, 1))*np.sqrt(mass[:, None]) for v in np.eye(3)]
        rigid += [np.cross(xyz, v)*np.sqrt(mass[:, None]) for v in np.eye(3)]
        basis = null_space(np.array(rigid).reshape(6, -1))
        # Analytic Hessian with known positive eigenvalues in the physical
        # vibrational space; this is not a computed nitromethane spectrum.
        weighted = .01 * (basis @ basis.T)
        hess = weighted*np.sqrt(np.outer(np.repeat(mass, 3), np.repeat(mass, 3)))
        raw, reduced = frequencies.get_frequencies(species, hess, species.geom)
        self.assertEqual((len(raw), len(reduced)), (15, 14))


if __name__ == '__main__':
    unittest.main()
