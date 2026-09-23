import json
import os
from pathlib import Path
from dataclasses import replace
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
from ase import Atoms
from ase.build import molecule
from kinbot.conformer_records import inventory, retain
from kinbot.conformer_counting import evaluate_members, writer_members, CountingError
from kinbot.conformers import Conformers
from kinbot.symmetry import conformer_symmetry_numbers
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import optical_scope


def peroxide():
    atoms = Atoms('OOHH', positions=[[0, 0, 0], [1.45, 0, 0],
                                    [-.3, .9, 0], [1.75, 0, .9]])
    point = StationaryPoint('peroxide', 0, 1, atom=atoms.get_chemical_symbols(),
                            geom=atoms.positions)
    point.characterize()
    return point


class TestConformerCounting(unittest.TestCase):
    def records(self, species, geoms):
        return inventory(species, geoms, [-1.] * len(geoms),
                         [[100., 200., 300., 400., 500., 600.]] * len(geoms),
                         [0] * len(geoms))

    def test_explicit_and_representative_mirror_partition_weights_agree(self):
        species = peroxide()
        one, _ = evaluate_members(species, self.records(species, [species.geom]))
        two, groups = evaluate_members(species, self.records(
            species, [species.geom, species.geom * [-1, 1, 1]]))
        self.assertEqual(one[0].mirror_states, 2)
        self.assertEqual(one[0].remaining_optical_weight, 2)
        self.assertEqual(groups, [[0, 1]])
        self.assertEqual(one[0].remaining_optical_weight / one[0].sigma_ext,
                         sum(record.remaining_optical_weight / record.sigma_ext for record in two))

    def test_supported_reference_cannot_license_an_unsupported_stable_member(self):
        species = peroxide()
        optical_scope(species)
        species.conformer_index = [7]
        species.conformer_geom = [species.geom.copy()]
        species.conformer_zeroenergy = [-1.]
        species.conformer_freq = [[100.] * 6]
        with TemporaryDirectory() as directory:
            previous = Path.cwd()
            try:
                os.chdir(directory)
                with patch('kinbot.conformer_counting.canonical_identity', return_value={
                        'status': 'unsupported', 'reason': 'unassigned observed configuration'}):
                    with self.assertLogs('KinBot', level='WARNING'):
                        self.assertEqual(writer_members(species), {})
                self.assertEqual(species.conformer_inventory[0].index, 7)
                self.assertEqual(species.conformer_records, {})
                self.assertEqual(species.conformer_inventory[0].exclusion_reason,
                                 'unassigned configured stereoisomer')
            finally:
                os.chdir(previous)
        # A TS union graph intentionally uses its supported reactant reference.
        species.wellorts = 1
        records = self.records(species, [species.geom])
        result, groups = evaluate_members(species, records)
        self.assertEqual(groups, [[0]])
        self.assertIsNone(result[0].stereo_identity)

    def test_symmetry_and_optical_weights_precede_population_filter(self):
        species = peroxide()
        search = object.__new__(Conformers)
        search.species = species
        search.optical_population = 'specified'
        reflected = species.geom * [-1, 1, 1]
        search.find_unique([species.geom, reflected], [-1., -1.],
                           [[100., 200., 300., 400., 500., 600.]] * 2,
                           [0, 0], temp=300., boltz=.05)
        records = species.conformer_records
        self.assertEqual(len(records), 2)
        self.assertTrue(all(r.sigma_ext == 2 and r.remaining_optical_weight == 1
                            and r.population_ratio == 1 for r in records.values()))

    def test_proper_atom_permutations_remove_only_duplicates(self):
        species = peroxide()
        permuted = species.geom[[1, 0, 3, 2]]
        records, groups = evaluate_members(species, self.records(species, [species.geom, permuted]))
        self.assertEqual(groups, [[0]])
        self.assertEqual(records[1].duplicate_of, 0)

    def test_inconsistent_duplicate_and_mirror_energies_warn_without_mixing_properties(self):
        species = peroxide()
        for second in (species.geom.copy(), species.geom * [-1, 1, 1]):
            records = list(self.records(species, [species.geom, second]))
            records[1] = replace(records[1], zero_energy_hartree=-1.01,
                                 frequencies_cm1=(111.,)*6, source_job='lower_calculation')
            with self.assertLogs('KinBot', level='WARNING'):
                result, groups = evaluate_members(species, records)
            self.assertEqual(result[1].source_job, 'lower_calculation')
            self.assertEqual(result[1].frequencies_cm1, (111.,)*6)
            if np.array_equal(second, species.geom):
                self.assertEqual(groups, [[1]])
                self.assertEqual(result[0].duplicate_of, 1)
                retain(species, result, [1])
                species.conformer_index = [1]
                species.conformer_geom = [result[1].geometry]
                species.conformer_zeroenergy = [result[1].zero_energy_hartree]
                species.conformer_freq = [result[1].frequencies_cm1]
                self.assertTrue(writer_members(species)[0].optical_evidence['warnings'])
            else:
                self.assertEqual(groups, [[1, 0]])
                self.assertEqual([r.remaining_optical_weight for r in result], [1., 1.])
            # Original observations remain intact, not converted into states.
            self.assertEqual(records[1].zero_energy_hartree, -1.01)

    def test_mirror_energy_tolerance_does_not_relax_duplicate_tolerance(self):
        from kinbot.constants import AUtoKCAL
        species = peroxide()
        records = list(self.records(species, [species.geom, species.geom * [-1, 1, 1]]))
        for delta in (.3644, .99, 1.01):
            with self.subTest(delta_kcal_mol=delta):
                records[1] = replace(records[1], zero_energy_hartree=-1. + delta / AUtoKCAL)
                result, groups = evaluate_members(species, records)
                self.assertEqual(groups, [[0, 1]])
                self.assertEqual([r.remaining_optical_weight for r in result], [1., 1.])
                self.assertEqual(result[1].zero_energy_hartree, records[1].zero_energy_hartree)
                self.assertEqual(bool(result[0].optical_evidence.get('warnings')), delta > 1.)
        records[1] = replace(records[1], geometry=records[0].geometry,
                             zero_energy_hartree=-1. + .3644 / AUtoKCAL)
        with self.assertLogs('KinBot', level='WARNING'):
            result, groups = evaluate_members(species, records)
        self.assertEqual(groups, [[0]])
        self.assertEqual(result[1].duplicate_of, 0)

    def test_optical_tolerance_does_not_merge_distinct_proper_geometries(self):
        from kinbot.molecular_symmetry import equivalent_geometry
        species = peroxide()
        second = species.geom.copy()
        second[3, 2] += .2
        self.assertFalse(equivalent_geometry(species, species.geom, second, .05))
        self.assertTrue(equivalent_geometry(species, species.geom, second, .1))
        records, groups = evaluate_members(species, self.records(species, [species.geom, second]))
        self.assertEqual(groups, [[0], [1]])
        self.assertTrue(all(record.duplicate_of is None for record in records))

    def test_duplicate_winner_is_independent_of_observation_order(self):
        species = peroxide()
        records = list(self.records(species, [species.geom, species.geom.copy()]))
        records[1] = replace(records[1], zero_energy_hartree=-1.0001)
        for observations in (records, records[::-1]):
            result, groups = evaluate_members(species, observations)
            self.assertEqual([result[i].index for group in groups for i in group], [1])

    def test_nontransitive_mirror_matches_cannot_make_a_three_state_orbit(self):
        species = peroxide()
        left = species.geom * [-1, 1, 1]
        right = left.copy()
        left[3, 2] += .12
        right[3, 2] -= .12
        with self.assertLogs('KinBot', level='WARNING'):
            records, groups = evaluate_members(species, self.records(species, [species.geom, left, right]))
        self.assertEqual(groups, [[0, 1, 2]])
        self.assertEqual([r.remaining_optical_weight for r in records], [1., 1., 1.])
        self.assertTrue(all(r.optical_evidence['states_covered_by_mc_tst'] is None for r in records))
        self.assertTrue(all(r.optical_evidence['status'] == 'assumed' for r in records))

    def test_distortion_does_not_add_a_symmetry_penalty_to_population_cutoff(self):
        species = peroxide()
        second = species.geom.copy()
        second[3, 2] += .5
        search = object.__new__(Conformers)
        search.species = species
        geoms = [species.geom, second]
        moments = [Atoms(species.atom, positions=geom).get_moments_of_inertia()
                   for geom in geoms]
        # Same electronic, vibrational and translational factors; rotational
        # contributions differ by sqrt(Ia Ib Ic)/sigma. Both omit a mirror.
        expected = np.sqrt(np.prod(moments[0]) / np.prod(moments[1]))
        result = search.find_unique(geoms, [-1., -1.], [[100.] * 6] * 2,
                                    [0, 0], temp=300., boltz=.6)
        self.assertEqual(result[-1], [0, 1])
        self.assertAlmostEqual(species.conformer_inventory[0].population_ratio, expected)

    def test_canonical_evidence_refreshes_final_member_weights(self):
        from kinbot.thermochemistry import thermochemistry_evidence
        species = peroxide()
        species.energy, species.zpe = -1.01, .01
        species.freq = species.reduced_freqs = [100.] * 6
        species.conformer_geom = [species.geom, None]
        species.conformer_index = [7, -999]
        species.conformer_freq = [[100.] * 6, None]
        species.conformer_zeroenergy = [-1., -999.]
        species.conformer_representation = 'MC-RRHO'
        evidence = thermochemistry_evidence(species)
        record = evidence['conformer_inventory'][0]
        self.assertEqual(record['index'], 7)
        self.assertEqual(record['remaining_optical_weight'], 2.)
        self.assertFalse(evidence['rate_model_ready'])
        json.dumps(evidence, allow_nan=False)

    def test_unknown_scope_preserves_observations_with_warned_legacy_weights(self):
        from kinbot.thermochemistry import thermochemistry_evidence
        species = peroxide()
        species.optical_reference = {'status': 'unsupported', 'reason': 'axial configuration'}
        species.energy, species.zpe = -1.01, .01
        species.freq = species.reduced_freqs = [100.] * 6
        species.conformer_index = [7]
        species.conformer_geom = [species.geom]
        species.conformer_zeroenergy = [-1.]
        species.conformer_freq = [[100.] * 6]
        species.conformer_representation = 'MC-RRHO'
        with self.assertLogs('KinBot', level='WARNING'):
            members = writer_members(species)
        self.assertEqual(members[0].optical_evidence['status'], 'legacy_unverified')
        self.assertEqual(members[0].remaining_optical_weight, 1.)
        self.assertEqual(species.conformer_inventory[0].index, 7)
        self.assertIsNotNone(species.conformer_inventory[0].geometry)
        result = thermochemistry_evidence(species)
        self.assertFalse(result['rate_model_ready'])
        json.dumps(result, allow_nan=False)
        records, groups = evaluate_members(species, self.records(species, [species.geom]), strict=False)
        self.assertEqual(groups, [[0]])
        self.assertEqual(records[0].remaining_optical_weight, 1.)

    def test_methane_symmetry_ratio_already_contains_four_equivalent_sites(self):
        methane = molecule('CH4')
        reactant = StationaryPoint('methane', 0, 1,
            atom=methane.get_chemical_symbols(), geom=methane.positions)
        reactant.characterize()
        # Collinear abstraction axis: C--H--H, with three spectator H atoms.
        xyz = [[0, 0, 0], [0, 0, 1.3], [0, 0, 2.3]] + [
            [np.cos(a), np.sin(a), -.3] for a in np.arange(3) * 2*np.pi/3]
        ts = StationaryPoint('ts', 0, 2, atom=['C'] + ['H'] * 5,
                             geom=np.array(xyz), wellorts=1)
        bonds = np.zeros((6, 6), int)
        for i, j in [(0, 1), (1, 2), (0, 3), (0, 4), (0, 5)]:
            bonds[i, j] = bonds[j, i] = 1
        ts.characterize(bond_mx=bonds)
        self.assertEqual(conformer_symmetry_numbers(reactant)['sigma_ext'], 12)
        self.assertEqual(conformer_symmetry_numbers(ts)['sigma_ext'], 3)
        self.assertEqual(12 / 3, 4)  # No extra raw site factor of four.

    def test_distorted_methane_keeps_legacy_symmetry_in_mc_and_non_mc_screening(self):
        atoms = molecule('CH4')
        xyz = atoms.positions + np.random.default_rng(40).normal(0, .04, atoms.positions.shape)
        point = StationaryPoint('distorted', 0, 1, atom=atoms.get_chemical_symbols(), geom=xyz)
        point.characterize()
        self.assertEqual(point.sigma_ext, -1)  # Symmetry finalization happens later.
        search = object.__new__(Conformers)
        search.species = point
        for strict in (False, True):
            search.strict_counting = strict
            result = search.find_unique([point.geom], [-1.], [[100.] * 9], [0], temp=300., boltz=.05)
            self.assertEqual(result[-1], [0])
            record = point.conformer_inventory[0]
            self.assertEqual(record.population_ratio, 1.)
            self.assertEqual(record.sigma_ext, 12)
            self.assertEqual(point.sigma_ext, -1)

    def test_repeated_writers_and_error_exports_preserve_excluded_member_sources(self):
        from kinbot.thermochemistry import thermochemistry_evidence
        for unsupported in (False, True):
            species = peroxide()
            species.energy, species.zpe = -1.01, .01
            species.freq = species.reduced_freqs = [100.] * 6
            species.conformer_geom = [species.geom.copy(), species.geom.copy()]
            species.conformer_index = [0, 1]
            species.conformer_zeroenergy = [-1., -1.]
            species.conformer_freq = [[100.] * 6] * 2
            species.conformer_representation = 'MC-RRHO'
            records = [replace(r, source_job=f'conf_{r.index}_high',
                electronic_energy_hartree=-1.01, zpe_hartree=.01)
                for r in self.records(species, species.conformer_geom)]
            retain(species, records, [0, 1])
            if unsupported:
                species.optical_reference = {'status': 'unsupported', 'reason': 'fixture'}
            for iteration in range(3):
                with self.subTest(unsupported=unsupported, iteration=iteration):
                    self.assertEqual(list(writer_members(species)), [0])
                    exported = thermochemistry_evidence(species)['conformer_inventory']
                    self.assertEqual([r['source_job'] for r in exported], ['conf_0_high', 'conf_1_high'])
                    self.assertTrue(all(r['electronic_energy_hartree'] == -1.01
                                        and r['zpe_hartree'] == .01 for r in exported))
                    self.assertFalse(exported[1]['retained'])


if __name__ == '__main__':
    unittest.main()
