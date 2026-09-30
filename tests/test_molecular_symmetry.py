"""Legacy rotational rules remain separate from spatial duplicate checks."""
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
import numpy as np
from ase.build import molecule
from kinbot.molecular_symmetry import equivalent_geometry, geometric_mirror_states
from kinbot.symmetry import calculate_symmetry, conformer_symmetry_numbers
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereochemistry import conformer_symmetry


class TestLegacyConformerSymmetry(unittest.TestCase):
    def test_rotational_benchmark_including_pyramidal_correction_for_members(self):
        with Path(__file__).with_name('symmetry_data.json').open() as handle:
            cases = json.load(handle)
        for name, case in cases.items():
            with self.subTest(name=name):
                # Match the original benchmark's SMILES-based construction.
                point = StationaryPoint(name, 0, 1, smiles=name)
                point.characterize()
                numbers = conformer_symmetry_numbers(point)
                self.assertEqual([numbers['sigma_ext'], numbers['nopt']],
                                 case['expected_values'][::2])
                self.assertEqual(point.sigma_ext, -1)

    def test_saved_methoxy_keeps_three_in_member_metadata(self):
        point = methoxy()
        calculate_symmetry(point)
        original = point.geom.copy()
        bond = point.bond.copy()
        point.conformer_index = [4, -999, 8]
        point.conformer_geom = [point.geom, None, point.geom * [-1, 1, 1]]
        for tolerance in (.03, .05, .1):
            records = conformer_symmetry(point, tolerance)
            self.assertEqual([r['sigma_ext'] for r in records], [3, 3])
            self.assertEqual([r['index'] for r in records], [4, 8])
        np.testing.assert_array_equal(point.geom, original)
        np.testing.assert_array_equal(point.bond, bond)
        self.assertEqual(point.sigma_ext, 3)

    def test_member_geometry_reaches_existing_linear_axis_rule(self):
        atoms = molecule('CO2')
        point = StationaryPoint('CO2', 0, 1, atom=atoms.get_chemical_symbols(), geom=atoms.positions)
        point.characterize()
        bent = point.geom.copy()
        bent[2, 0] += .3
        self.assertEqual(conformer_symmetry_numbers(point, bent)['sigma_ext'], 2)
        np.testing.assert_array_equal(point.geom, atoms.positions)

    def test_ch2oh_keeps_existing_carbon_oxygen_rotor_and_symmetry(self):
        with (Path(__file__).parent / 'reference' / 'ch2oh_continuations.json').open() as handle:
            data = json.load(handle)
        geom = data['scans']['fc']['points'][0]['geometry']
        point = StationaryPoint('CH2OH', 0, 2,
            atom=data['provenance']['atoms'], geom=np.array(geom))
        point.characterize()
        calculate_symmetry(point)
        self.assertTrue(any(set(rotor[1:3]) == {0, 1} for rotor in point.dihed))
        self.assertEqual(point.sigma_int[0][1], 2)
        self.assertEqual(conformer_symmetry_numbers(point)['sigma_int'][0][1], 2)


class TestSpatialMappings(unittest.TestCase):
    def test_chiral_geometry_retains_optical_pair_detection(self):
        atoms = molecule('CH4')
        point = SimpleNamespace(atom=['C', 'H', 'F', 'Cl', 'Br'], geom=atoms.positions)
        point.bond = np.zeros((5, 5), dtype=int)
        point.bond[0, 1:] = point.bond[1:, 0] = 1
        self.assertEqual(geometric_mirror_states(point), 2)
        self.assertFalse(equivalent_geometry(point, point.geom, point.geom * [-1, 1, 1]))

    def test_resonance_mappings_and_rigid_transform_duplicates_are_preserved(self):
        for name in ('C6H6', 'NO2'):
            atoms = molecule(name)
            point = StationaryPoint(name, 0, 1 if name == 'C6H6' else 2,
                atom=atoms.get_chemical_symbols(), geom=atoms.positions)
            point.characterize()
            self.assertGreater(len(point.bonds), 1)
            rotated = point.geom @ np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]) + 7
            for bond in point.bonds:
                point.bond = bond
                self.assertTrue(equivalent_geometry(point, point.geom, rotated))
                self.assertEqual(geometric_mirror_states(point), 1)


def methoxy():
    """Saved MeOH_H product geometry; no calculation-directory dependency."""
    with (Path(__file__).parent / 'reference' / 'methoxy_symmetry.json').open() as handle:
        data = json.load(handle)
    point = StationaryPoint('methoxy', 0, 2, atom=data['atom'], geom=np.array(data['geometry']))
    point.characterize()
    point.freq = point.reduced_freqs = data['frequencies']
    point.energy, point.zpe = data['energy'], data['zpe']
    return point


if __name__ == '__main__':
    unittest.main()
