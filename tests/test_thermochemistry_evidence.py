"""Evidence export must never turn missing HIR data into measured energies."""
import json
import unittest
from unittest.mock import patch

import numpy as np
from ase.build import molecule

from kinbot.hindered_rotors import HIR
from kinbot.stationary_pt import StationaryPoint
from kinbot import constants, frequencies, symmetry
from kinbot.thermochemistry import thermochemistry_evidence


def scanned_point():
    atoms = molecule('CH3OH')
    species = StationaryPoint('methanol', 0, 1,
                              atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    species.characterize()
    symmetry.calculate_symmetry(species)
    species.energy, species.zpe = -100., .05
    species.source_job = 'selected_calculation'
    species.freq, species.reduced_freqs = [100.] * 12, [100.] * 11
    hir = HIR(species, None, {'nrotation': 12, 'plot_hir_profiles': False, 'rotor_0_test': 1})
    species.hir = hir
    angles = np.arange(12) * 2 * np.pi / 12
    hir.hir_status = [[0] * 12]
    hir.hir_energies = [(-100. + (1 - np.cos(3 * angles)) / constants.AUtoKCAL).tolist()]
    hir.hir_geoms = [[species.geom.copy() for _ in angles]]
    return species, angles


class TestThermochemistryEvidence(unittest.TestCase):
    def test_raw_and_interpolated_points_remain_distinguishable(self):
        species, angles = scanned_point()
        species.hir.hir_status[0][2] = 1
        species.hir.hir_energies[0][2] = -1.
        species.hir.fourier_fit('test', angles, 0)
        species.optical_isomers = 2
        with patch('builtins.open', side_effect=AssertionError('No file reads allowed')):
            data = thermochemistry_evidence(species)
        json.dumps(data, allow_nan=False)
        rotor = data['hir']['rotors'][0]
        self.assertEqual(rotor['axis'], species.dihed[0][1:3])
        self.assertEqual(rotor['sigma_int'], 3)
        self.assertEqual(len(rotor['points']), 12)
        failed = rotor['points'][2]
        self.assertEqual(failed['status'], 'failed')
        self.assertIsNone(failed['electronic_energy_hartree'])
        self.assertIsNone(failed['relative_energy_kcal_mol'])
        # Master's six-term fit fills this missing point with 6/7 kcal/mol;
        # an interpolated value must remain distinct from a measured point.
        self.assertAlmostEqual(failed['fitted_relative_energy_kcal_mol'], 6. / 7.)
        self.assertIsNone(rotor['fourier']['diagnostics'])
        self.assertEqual(data['optical_counting']['total_optical_states'], 1)
        self.assertTrue(data['optical_counting']['reported_count_disagrees'])
        self.assertIsNone(data['optical_counting']['remaining_multiplier'])
        self.assertFalse(data['rate_model_ready'])
        self.assertEqual(data['source_job'], 'selected_calculation')
        self.assertEqual(len(data['raw_harmonic_frequencies_cm-1']), 12)
        self.assertEqual(len(data['thermochemical_frequencies_cm-1']), 11)

    def test_projection_and_demotion_are_exported_without_claiming_coverage(self):
        species, _ = scanned_point()
        species.hir.hir_status[0] = [1] * 12
        _, species.reduced_freqs = frequencies.get_frequencies(
            species, np.eye(3 * species.natom), species.geom)
        data = thermochemistry_evidence(species)
        self.assertFalse(data['hir']['rotors'][0]['usable'])
        self.assertEqual(data['frequency_projection']['internal_rank'], 0)
        self.assertIsNone(data['hir']['rotors'][0]['reference_electronic_energy_hartree'])
        self.assertEqual(len(data['thermochemical_frequencies_cm-1']), 12)

    def test_no_scan_does_not_manufacture_an_empty_successful_hir(self):
        species, _ = scanned_point()
        species.hir = None
        self.assertEqual(thermochemistry_evidence(species)['hir'],
                         {'status': 'not_recorded', 'rotors': []})


if __name__ == '__main__':
    unittest.main()
