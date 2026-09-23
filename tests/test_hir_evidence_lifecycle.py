"""Preserve measured HIR data independently of acceptance and later selection."""
import copy
import json
import unittest
from unittest.mock import Mock, patch

import numpy as np

from kinbot import frequencies
from kinbot.thermochemistry import thermochemistry_evidence
from tests.test_thermochemistry_evidence import scanned_point


class TestHIREvidenceLifecycle(unittest.TestCase):
    def test_failed_getter_placeholders_are_not_measurements(self):
        for status in (-1, 1):
            species, _ = scanned_point()
            hir = species.hir
            hir.hir_status[0][1] = -1
            hir.qc = Mock(qc='test_backend')
            hir.qc.get_qc_geom.return_value = (0, species.geom.copy())
            hir.qc.get_qc_energy.return_value = (status, 0.)
            hir.test_hir()
            point = thermochemistry_evidence(species)['hir']['rotors'][0]['points'][1]
            self.assertIsNone(point['observation']['electronic_energy_hartree'])
            self.assertEqual(point['observation']['qc_energy_status'], 'failed' if status == -1 else 'running')
            self.assertIsNotNone(point['geometry_angstrom'])
        hir.hir_status[0][1] = -1
        hir.qc.get_qc_geom.return_value = (-1, np.zeros_like(species.geom))
        hir.test_hir()
        point = thermochemistry_evidence(species)['hir']['rotors'][0]['points'][1]
        self.assertIsNone(point['geometry_angstrom'])
        self.assertIsNone(point['observation']['geometry_angstrom'])
        hir.point_observations = []  # The legacy failure placeholder is also unavailable.
        self.assertIsNone(thermochemistry_evidence(species)['hir']['rotors'][0]['points'][1]['geometry_angstrom'])

    def test_bond_check_rejection_keeps_returned_geometry(self):
        species, _ = scanned_point()
        hir = species.hir
        hir.hir_status[0][1] = -1
        measured = species.geom * 2.
        hir.qc = Mock(qc='test_backend')
        hir.qc.get_qc_geom.return_value = (0, measured)
        hir.test_hir()
        hir.qc.get_qc_energy.assert_not_called()
        point = thermochemistry_evidence(species)['hir']['rotors'][0]['points'][1]
        self.assertIn('bond-length', point['observation']['rejection_reason'])
        np.testing.assert_array_equal(point['geometry_angstrom'], measured)
        self.assertIsNone(point['observation']['electronic_energy_hartree'])

    def test_rejected_high_point_keeps_measured_energy_and_geometry(self):
        species, _ = scanned_point()
        hir = species.hir
        hir.hir_status[0][1] = -1
        hir.qc = Mock(qc='test_backend')
        hir.qc.get_qc_geom.return_value = (0, species.geom.copy())
        hir.qc.get_qc_energy.return_value = (0, -99.)
        hir.test_hir()
        self.assertEqual(hir.hir_status[0][1], 1)
        self.assertEqual(hir.hir_energies[0][1], -1.)
        point = thermochemistry_evidence(species)['hir']['rotors'][0]['points'][1]
        self.assertIsNone(point['electronic_energy_hartree'])
        self.assertEqual(point['observation']['electronic_energy_hartree'], -99.)
        self.assertIn('20 kcal/mol', point['observation']['rejection_reason'])
        np.testing.assert_array_equal(point['observation']['geometry_angstrom'], species.geom)

    def test_global_demotion_keeps_prior_successful_observations(self):
        species, _ = scanned_point()
        species.energy -= .01
        measured = list(species.hir.hir_energies[0])
        with patch.object(species.hir, 'test_hir'):
            self.assertEqual(species.hir.check_hir(), 1)
        self.assertEqual(species.hir.hir_status[0], [1] * 12)
        data = thermochemistry_evidence(species)['hir']
        self.assertIn('rotor-zero', data['demotion_reason'])
        self.assertEqual([p['observation']['electronic_energy_hartree']
                          for p in data['rotors'][0]['points']], measured)
        self.assertFalse(data['rotors'][0]['usable'])

    def test_scan_and_projection_references_do_not_follow_mutable_species(self):
        species, _ = scanned_point()
        species.source_row_id = np.int64(7)
        species.hir.qc = Mock(qc='test_backend')
        species.hir.qc.qc_hir.return_value = 'actual_scan_job'
        original = species.geom.copy()
        species.hir.generate_hir_geoms(original, rigid=False)
        hessian = np.eye(species.natom * 3)
        frequencies.get_frequencies(species, hessian, original)
        scan = copy.deepcopy(species.hir.scan_reference)
        projection = copy.deepcopy(species.rotor_projection)
        species.geom = species.geom + 1.
        species.source_job = 'later_calculation'
        species.hir.qc.qc = 'later_backend'
        data = thermochemistry_evidence(species)
        self.assertEqual(species.hir.scan_reference, scan)
        self.assertEqual(species.rotor_projection, projection)
        self.assertEqual(data['hir']['scan_reference']['source_job'], 'selected_calculation')
        self.assertEqual(data['hir']['backend'], 'test_backend')
        self.assertEqual(data['frequency_projection']['reference']['source_row_id'], 7)
        self.assertFalse(data['frequency_projection']['reference']['hessian_massweighted'])
        json.dumps(data, allow_nan=False)

if __name__ == '__main__':
    unittest.main()
