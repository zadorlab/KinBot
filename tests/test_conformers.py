"""Conformer workflow regressions without electronic-structure jobs."""

import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from ase.build import molecule
from ase.thermochemistry import IdealGasThermo
import numpy as np

from kinbot.optimize import Optimize
from kinbot.conformers import Conformers
from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint


class TestConformerWeights(unittest.TestCase):
    def setUp(self):
        self.conformers = Conformers.__new__(Conformers)
        atoms = molecule('H2O')
        self.conformers.species = SimpleNamespace(atom=atoms.get_chemical_symbols(), mult=1)
        self.geometries = [atoms.positions.copy(), atoms.positions * 1.2]

    def test_real_saddle_populations_exclude_only_the_reaction_coordinate(self):
        self.conformers.species.wellorts = 1
        raw = [[-1000., 1500., 3000.]] * 2
        with self.assertWarnsRegex(UserWarning, '1 imag modes removed'):
            retained = self.conformers.find_unique(
                self.geometries, [-76., -75.99], raw, [0, 0], temp=300., boltz=1.e-4)
        self.assertEqual(retained[-1], [0])
        with self.assertWarnsRegex(UserWarning, '1 imag modes removed'):
            retained = self.conformers.find_unique(
                self.geometries, [-76., -75.99], raw, [0, 0], temp=300., boltz=1.e-6)
        self.assertEqual(retained[-1], [0, 1])
        self.assertEqual(retained[2], [-76., -75.99])

    def test_saddle_population_preserves_accepted_secondary_soft_mode(self):
        self.conformers.species.wellorts = 1
        raw = [[-1000., -20., 3000.]]
        with patch('kinbot.conformers.IdealGasThermo', wraps=IdealGasThermo) as thermo, \
                self.assertWarnsRegex(UserWarning, '1 imag modes removed'):
            retained = self.conformers.find_unique(
                [self.geometries[0]], [-76.], raw, [0], temp=300., boltz=.001)
        supplied = thermo.call_args.kwargs['vib_energies']
        self.assertEqual(sum(np.iscomplex(value) for value in supplied), 1)
        self.assertEqual(sum(np.real(value) > 0 for value in supplied), 2)
        self.assertEqual(raw, [[-1000., -20., 3000.]])
        self.assertEqual(retained[-1], [0])

    def test_equal_total_energies_do_not_acquire_a_second_zpe_bias(self):
        # Both supplied E + ZPE values are equal. Different high-frequency
        # modes must not add their ground-state energies a second time.
        result = self.conformers.find_unique(
            self.geometries, [0., 0.],
            [[1000., 1500., 3000.], [2000., 2500., 4000.]], [0, 0],
            temp=298.15, boltz=.1)
        self.assertEqual(result[-1], [0, 1])
        self.assertEqual(result[1], [0., 0.])

    def test_invalid_conformer_cannot_set_the_boltzmann_reference(self):
        result = self.conformers.find_unique(
            self.geometries, [0., -100.], [[1000., 1500., 3000.]] * 2,
            [0, 1], temp=298.15, boltz=.1)
        self.assertEqual(result[-1], [0])

    def test_final_l1_conformers_keep_their_zero_point_inclusive_energies(self):
        energies = [-76.0, -75.99]
        result = self.conformers.find_unique(
            self.geometries, energies, [[1000., 1500., 3000.]] * 2,
            [0, 0])
        self.assertEqual(result[-1], [0, 1])
        self.assertEqual(result[2], energies)
        # L2 refinement may later replace either array independently.
        self.assertIsNot(result[1], result[2])

    def test_failed_frequency_record_does_not_crash_population_filtering(self):
        for valid_index in (0, 1):
            with self.subTest(valid_index=valid_index):
                valid = [1, 1]
                valid[valid_index] = 0
                frequencies = [None, None]
                frequencies[valid_index] = [1000., 1500., 3000.]
                result = self.conformers.find_unique(
                    self.geometries, [-76., -76.], frequencies, valid,
                    temp=298.15, boltz=.1)
                self.assertEqual(result[-1], [valid_index])


class TestSemiEmpiricalConformers(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        Path('input.json').write_text(json.dumps({
            'barrier_threshold': 50., 'conformer_search': 1,
            'semi_emp_conformer_search': 1, 'semi_emp_confomer_threshold': 5.,
        }))
        self.par = Parameters('input.json', show_warnings=False).par
        atoms = molecule('CH4')
        self.species = StationaryPoint(
            'methane', 0, 1, atom=atoms.get_chemical_symbols(), geom=atoms.positions)
        self.species.characterize()
        self.optimization = Optimize(self.species, self.par, None)
        self.geom = self.species.geom.copy()
        self.regular = Mock(cyc_conf_geoms=[])
        self.regular.check_conformers.return_value = (0, '0000', self.geom, -10., [], [], [], [])
        self.semi = Mock()
        self.semi.check_conformers.side_effect = [
            (0, '0000', self.geom, -10., [], [], [], []),
            (1, '0000', self.geom, -10., [self.geom], [-9.95], [[100.]], [0]),
        ]
        factory = patch('kinbot.optimize.Conformers', side_effect=lambda *a, **kw:
                        self.semi if kw.get('semi_emp') else self.regular)
        factory.start()
        self.addCleanup(factory.stop)

    def test_running_search_is_polled_again_without_resubmission(self):
        self.optimization.do_optimization()
        self.assertEqual(self.optimization.ssemi_empconf, 0)
        self.regular.generate_conformers.assert_not_called()
        self.optimization.do_optimization()
        self.assertEqual(self.optimization.ssemi_empconf, 1)
        self.assertEqual(self.semi.generate_conformers.call_count, 1)
        self.assertEqual(self.semi.check_conformers.call_count, 2)
        self.assertEqual(self.regular.generate_conformers.call_count, 1)

    def test_waits_for_ring_conformers_before_starting_semi_empirical_jobs(self):
        self.species.cycle_chain = [[0, 1, 2, 3]]
        self.regular.check_ring_conformers.side_effect = [(0, []), (1, [self.geom])]
        self.optimization.do_optimization()
        self.semi.generate_conformers.assert_not_called()
        self.optimization.do_optimization()
        self.assertEqual(self.semi.generate_conformers.call_count, 1)

    def test_screening_compares_valid_total_energies(self):
        self.semi.check_conformers.side_effect = None
        self.semi.check_conformers.return_value = (
            1, '0000', self.geom, -10.,
            [self.geom, self.geom + .1, self.geom + .2],
            [-9.95, -9.94, -20.], [[100.]] * 3, [0, 0, 1])
        self.optimization.do_optimization()
        self.assertEqual(self.regular.generate_conformers.call_count, 1)
        np.testing.assert_array_equal(self.regular.generate_conformers.call_args.args[1], self.geom)
        # Surviving geometries are optimized as given, without a new scan.
        self.assertEqual(self.regular.generate_conformers.call_args.args[0], -999)

    def test_all_failed_search_falls_back_to_the_standard_search(self):
        self.semi.check_conformers.side_effect = None
        self.semi.check_conformers.return_value = (
            1, '0000', self.geom, -10., [], [], [], [1])
        self.optimization.do_optimization()
        self.assertEqual(self.regular.generate_conformers.call_count, 1)
        np.testing.assert_array_equal(self.regular.generate_conformers.call_args.args[1], self.geom)
        # Without valid seeds the dihedrals still have to be searched, so the
        # rotor index must start the recursion rather than skip it.
        self.assertEqual(self.regular.generate_conformers.call_args.args[0], 0)
        self.assertNotEqual(self.optimization.species.confs.nconfs, 1)




    def test_saddle_skips_l0_and_polls_l1_without_l0_seed_results(self):
        self.species.wellorts = 1
        self.regular.generate_conformers.return_value = 0
        for _ in range(2):
            self.optimization.do_optimization()
        self.semi.generate_conformers.assert_not_called()
        self.semi.check_conformers.assert_not_called()
        self.assertEqual(self.optimization.ssemi_empconf, 1)
        self.regular.generate_conformers.assert_called_once()
        self.assertEqual(self.regular.generate_conformers.call_args.args[0], 0)
        self.assertEqual(self.regular.check_conformers.call_count, 2)


    def test_ring_saddle_keeps_ring_sampling_then_proceeds_directly_to_l1(self):
        self.species.wellorts = 1
        self.species.cycle_chain = [[0, 1, 2, 3]]
        self.regular.generate_conformers.return_value = 0
        self.regular.check_ring_conformers.side_effect = [(0, []), (1, [self.geom])]
        self.optimization.do_optimization()
        self.regular.generate_conformers.assert_not_called()
        self.optimization.do_optimization()
        self.regular.generate_ring_conformers.assert_called_once()
        self.semi.generate_conformers.assert_not_called()
        self.regular.generate_conformers.assert_called_once()
        self.assertEqual(self.regular.generate_conformers.call_args.args[0], 0)



if __name__ == '__main__':
    unittest.main()
