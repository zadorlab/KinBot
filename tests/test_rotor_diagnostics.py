import copy
import json
from pathlib import Path
import unittest
import numpy as np
from kinbot.rotor_diagnostics import (hydroxymethyl_coordinates, continuation_diagnostics,
                                     cross_continuation_diagnostics)


class TestRotorDiagnostics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.evidence = json.loads((Path(__file__).parent / 'reference/ch2oh_continuations.json').read_text())

    def test_both_backends_and_both_seeds_close_but_are_not_half_turn_periodic(self):
        for backend, data in self.evidence['scans'].items():
            for seed in ('original', 'reflected'):
                for direction in (1, -1):
                    points = [point for point in data['points']
                              if point['branch_seed'] == seed and point['direction'] == direction]
                    with self.subTest(backend=backend, seed=seed, direction=direction):
                        result = continuation_diagnostics(points)
                        self.assertTrue(result['complete'])
                        self.assertFalse(result['sampled_half_turn_periodicity'])
                        self.assertGreater(result['maximum_half_turn_difference_kcal_mol'], 1.)
                        self.assertLess(result['closure_fixed_label_rmsd_angstrom'], .001)
                        self.assertEqual(result['branch_sign_changes'], 4)
                        json.dumps(result, allow_nan=False)

    def test_hydrogen_exchange_acts_on_both_coordinates(self):
        geometry = np.array(self.evidence['scans']['fc']['points'][0]['geometry'])
        original = geometry.copy()
        first = hydroxymethyl_coordinates(geometry)
        second = hydroxymethyl_coordinates(geometry[[0, 1, 3, 2, 4]])
        self.assertAlmostEqual((first['hh_frame_torsion_degrees']-
                               second['hh_frame_torsion_degrees']) % 360., 180.)
        self.assertAlmostEqual(first['signed_pyramidalization_angstrom'],
                               -second['signed_pyramidalization_angstrom'])
        self.assertGreater(abs((first['hcoh_dihedral_degrees']-
                                second['hcoh_dihedral_degrees']) % 360.-180.), 10.)
        np.testing.assert_array_equal(original, geometry)
        rotated = hydroxymethyl_coordinates(geometry @ np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]) + 7.)
        for key in first:
            self.assertAlmostEqual(first[key], rotated[key])

    def test_saved_seed_collapse_and_b3lyp_direction_dependence_are_explicit(self):
        for backend, data in self.evidence['scans'].items():
            result = cross_continuation_diagnostics(data['points'])
            self.assertFalse(result['independent_branch_coverage_established'])
            for start in result['seed_starts']:
                self.assertLess(start['accepted_wag_angstrom'], -.13)
                if start['branch_seed'] == 'reflected':
                    self.assertGreater(start['start_wag_angstrom'], .13)
            seed_pairs = [pair for comparison in result['comparisons']
                          if comparison['left'][1] == comparison['right'][1]
                          for pair in comparison['points']]
            self.assertLess(max(p['proper_fixed_label_rmsd_angstrom'] for p in seed_pairs), .0005)
            direction = next(c for c in result['comparisons']
                             if c['left'] == ['original', 1] and c['right'] == ['original', -1])
            if backend == 'orca':
                for index, expected in [(4, .069990), (6, -.379500), (7, .255086)]:
                    pair = next(p for p in direction['points'] if p['left_index'] == index)
                    self.assertEqual(pair['right_index'], (12-index) % 12)
                    self.assertAlmostEqual(pair['energy_left_minus_right_kcal_mol'], expected, places=5)
                    self.assertGreater(pair['proper_fixed_label_rmsd_angstrom'], .1)
                    self.assertLess(pair['left_wag_angstrom'] * pair['right_wag_angstrom'], 0.)
            else:
                self.assertLess(max(p['proper_fixed_label_rmsd_angstrom'] for p in direction['points']), .0004)
            expected = 1.3752415 if backend == 'fc' else 1.4727308
            self.assertAlmostEqual(result['lowest_observed_envelope_half_turn_difference_kcal_mol'], expected, places=5)
            json.dumps(result, allow_nan=False)

    def test_failed_or_missing_points_leave_periodicity_unknown(self):
        points = copy.deepcopy(self.evidence['scans']['fc']['points'][:13])
        points[5]['status'] = 'failed'
        self.assertIsNone(continuation_diagnostics(points)['sampled_half_turn_periodicity'])
        self.assertFalse(continuation_diagnostics(points[:-1])['complete'])

    def test_constraint_drift_is_not_a_valid_periodicity_test(self):
        points = copy.deepcopy(self.evidence['scans']['fc']['points'][:13])
        points[5]['angle_actual_degrees'] += 5.
        self.assertIsNone(continuation_diagnostics(points)['sampled_half_turn_periodicity'])

    def test_nonfinite_successful_values_leave_periodicity_unknown_and_json_finite(self):
        for key, value in [('energy_ev', float('nan')), ('angle_actual_degrees', float('inf')),
                           ('angle_target_degrees', float('nan')), ('geometry', [[float('nan')] * 3] * 5)]:
            points = copy.deepcopy(self.evidence['scans']['fc']['points'][:13])
            points[3][key] = value
            result = continuation_diagnostics(points)
            self.assertFalse(result['complete'])
            self.assertIsNone(result['sampled_half_turn_periodicity'])
            self.assertEqual(result['invalid_points'][0]['index'], 3)
            json.dumps(result, allow_nan=False)


if __name__ == '__main__':
    unittest.main()
