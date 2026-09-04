"""Regression tests for hindered-rotor results, without QC calculations."""

import os
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, mock_open, patch

from ase.io import read
from ase.build import molecule
import numpy as np

from kinbot.hindered_rotors import HIR
from kinbot import constants, frequencies
from kinbot.mess import MESS
from kinbot.optimize import Optimize
from kinbot.stationary_pt import StationaryPoint


def completed_hir(status=None):
    species = SimpleNamespace(
        natom=4, atom=['C'] * 4, name='rotor_test', chemid=123,
        wellorts=0, energy=-100., dihed=[[0, 1, 2, 3], [1, 2, 3, 0]],
    )
    hir = HIR(species, None, {
        'nrotation': 12, 'plot_hir_profiles': False, 'rotor_0_test': True,
    })
    angles = np.arange(12) * 2 * np.pi / 12
    hir.hir_status = status or [[0] * 12, [0] * 12]
    hir.hir_energies = [
        (-100. + factor * (1 - np.cos(angles)) / constants.AUtoKCAL).tolist()
        for factor in (1., 3.)]
    hir.hir_geoms = [[np.zeros((4, 3)) for _ in angles] for _ in range(2)]
    species.hir = hir
    return hir


class TestHIRProfile(unittest.TestCase):
    def test_missing_failed_point_geometry_does_not_discard_other_frames(self):
        hir = HIR(SimpleNamespace(natom=1, atom=['H']), None, {
            'nrotation': 3, 'plot_hir_profiles': False, 'rotor_0_test': True,
        })
        hir.hir_energies = [[-100., -1., -99.99]]
        hir.hir_geoms = [[[[0., 0., 0.]], [], [[2., 0., 0.]]]]
        output = mock_open()
        with patch('builtins.open', output):
            hir.write_profile(0, 'partial')
        contents = ''.join(call.args[0] for call in output().write.call_args_list)
        frames = read(StringIO(contents), index=':', format='xyz')
        self.assertEqual(len(frames), 2)
        np.testing.assert_array_equal(frames[1].positions, [[2., 0., 0.]])

    def test_combined_xyz_contains_every_scan_point(self):
        with TemporaryDirectory() as directory:
            previous = Path.cwd()
            try:
                os.chdir(directory)
                Path('hir').mkdir()
                species = SimpleNamespace(natom=2, atom=['C', 'H'])
                hir = HIR(species, None, {
                    'nrotation': 3, 'plot_hir_profiles': False,
                    'rotor_0_test': True,
                })
                geometries = np.arange(18, dtype=float).reshape(3, 2, 3) / 10
                energies = [-40.0, -39.99, -39.98]
                hir.hir_geoms = [geometries]
                hir.hir_energies = [energies]

                hir.write_profile(0, 'profile')

                frames = read('hir/profile.xyz', index=':', format='xyz')
                self.assertEqual(len(frames), 3)
                for frame, positions in zip(frames, geometries):
                    self.assertEqual(frame.get_chemical_symbols(), ['C', 'H'])
                    np.testing.assert_allclose(frame.positions, positions)
                comments = Path('hir/profile.xyz').read_text().splitlines()[1::4]
                self.assertEqual(comments, [f'energy = {e}' for e in energies])
            finally:
                os.chdir(previous)


class TestHIRFits(unittest.TestCase):
    def test_completed_rotors_keep_distinct_coefficients(self):
        hir = completed_hir()
        with patch.object(hir, 'test_hir'), patch.object(hir, 'write_profile'), \
                patch.object(hir, 'fourier_fit', wraps=hir.fourier_fit) as fit:
            self.assertEqual(hir.check_hir(), 1)
        self.assertEqual(fit.call_count, 2)
        self.assertEqual(np.shape(hir.hir_fourier), (2, 12))
        self.assertAlmostEqual(hir.get_fit_value(np.pi, rotor=0), 2.)
        self.assertAlmostEqual(hir.get_fit_value(np.pi, rotor=1), 6.)

    def test_filling_failed_points_preserves_raw_results(self):
        hir = completed_hir()
        hir.n_terms = 2
        hir.hir_status[0][3] = 1
        hir.hir_energies[0][3] = -1.
        angles = np.arange(12) * 2 * np.pi / 12
        hir.fourier_fit('test', angles, 0)
        self.assertEqual(hir.hir_raw_energies[0][3], -1.)
        self.assertEqual(hir.hir_status[0][3], 1)
        self.assertAlmostEqual(
            hir.hir_energies[0][3], -100. + 1. / constants.AUtoKCAL)
        hir.fourier_fit('test', angles, 0)
        self.assertEqual(hir.hir_raw_energies[0][3], -1.)

    def test_mess_extra_point_uses_the_requested_rotor_fit(self):
        hir = completed_hir()
        angles = np.arange(12) * 2 * np.pi / 12
        for rotor in range(2):
            hir.fourier_fit('test', angles, rotor)
        writer = MESS.__new__(MESS)
        writer.par = {'free_rotor_thrs': 0.}
        writer.rotorsymm = lambda species, rotor: 6
        potential, kind = writer.make_rotorpot(hir.species, 0, hir.species.dihed[0], 1.)
        self.assertEqual(kind, 'hindered')
        self.assertEqual(float(potential.split()[1]), round(1 - np.cos(np.pi / 12), 2))


class TestHIRStatus(unittest.TestCase):
    def test_failed_first_rotor_does_not_trigger_a_conformer_restart(self):
        atoms = molecule('CH3OH')
        species = StationaryPoint('methanol', 0, 1,
                                  atom=atoms.get_chemical_symbols(), geom=atoms.positions)
        species.characterize()
        species.dihed = species.dihed * 2
        species.energy = -100.
        hir = HIR(species, None, {
            'nrotation': 12, 'plot_hir_profiles': False, 'rotor_0_test': True,
        })
        species.hir = hir
        hir.hir_status = [[1] * 12, [0] * 12]
        hir.hir_energies = [[-1.] * 12, [-100.] * 12]
        optimization = Optimize.__new__(Optimize)
        optimization.species = species
        optimization.name = species.name
        optimization.qc = SimpleNamespace(
            get_qc_geom=Mock(side_effect=AssertionError('Unexpected conformer restart')),
            read_qc_hess=lambda *args: np.eye(3 * species.natom),
            hessian_is_massweighted=lambda: False,
        )
        optimization.par = {
            'conformer_search': 0, 'rotation_restart': 3, 'high_level': 0,
            'rotor_scan': 1, 'multi_conf_tst': 0, 'L3_calc': 0,
        }
        optimization.shir = 0
        optimization.restart = 0
        optimization.just_high = False
        optimization.defer_hir = False
        optimization.wait = 0
        optimization.log_name = lambda *args, **kwargs: 'test'
        with patch.object(hir, 'test_hir'), patch.object(hir, 'write_profile'):
            optimization.do_optimization()
        self.assertEqual(optimization.restart, 0)
        self.assertEqual(optimization.shir, 1)

    def test_skipped_rotor_does_not_disable_successful_rotors(self):
        hir = completed_hir([[2] * 12, [0] * 12])
        hir.hir_energies[0] = [-1.] * 12
        with patch.object(hir, 'test_hir'), patch.object(hir, 'write_profile'):
            self.assertEqual(hir.check_hir(), 1)
        self.assertEqual(hir.hir_status, [[2] * 12, [0] * 12])
        self.assertFalse(hir.is_valid_rotor(0))
        self.assertTrue(hir.is_valid_rotor(1))

    def test_failure_count_warning_uses_each_rotors_status(self):
        hir = completed_hir()
        hir.hir_status[0][1:4] = [1, 1, 1]
        with patch.object(hir, 'test_hir'), patch.object(hir, 'write_profile'), \
                self.assertLogs('KinBot', level='WARNING') as logs:
            hir.check_hir()
        self.assertIn('More than 2 HIR calculations failed for 123_hir_0', '\n'.join(logs.output))

    def test_barrier_check_waits_for_reference_energy(self):
        atoms = molecule('CH3OH')
        species = StationaryPoint('methanol', 0, 1,
                                  atom=atoms.get_chemical_symbols(), geom=atoms.positions)
        species.characterize()
        qc = Mock()
        qc.get_qc_geom.side_effect = lambda job, natom: (
            1 if job.endswith('_00') else 0, species.geom)
        qc.get_qc_energy.side_effect = lambda job: (
            0, -100. if job.endswith('_00') else -99.95)
        hir = HIR(species, qc, {
            'nrotation': 12, 'plot_hir_profiles': False, 'rotor_0_test': True,
        })
        hir.hir_status = [[-1] * 12]
        hir.hir_energies = [[-1.] * 12]
        hir.hir_geoms = [[[] for _ in range(12)]]
        hir.test_hir()
        self.assertEqual(hir.hir_status, [[-1] * 12])
        qc.get_qc_geom.side_effect = lambda job, natom: (0, species.geom)
        hir.test_hir()
        self.assertEqual(hir.hir_status[0], [0] + [1] * 11)

    def test_failed_or_skipped_rotors_keep_their_harmonic_modes(self):
        atoms = molecule('CH3OH')
        for status in (-1, 0, 1, 2):
            with self.subTest(status=status):
                species = StationaryPoint('methanol', 0, 1,
                                          atom=atoms.get_chemical_symbols(), geom=atoms.positions)
                species.characterize()
                species.hir = HIR(species, None, {
                    'nrotation': 12, 'plot_hir_profiles': False, 'rotor_0_test': True,
                })
                species.hir.hir_status = [[status] * 12]
                raw, projected = frequencies.get_frequencies(
                    species, np.eye(3 * species.natom), species.geom)
                self.assertEqual(len(projected), len(raw) - (status == 0))
                if status != 0:
                    np.testing.assert_allclose(projected, raw, rtol=1.e-8)

                writer = MESS.__new__(MESS)
                writer.par = {'rotor_scan': 1}
                writer.freerotortpl = 'free rotor'
                writer.make_rotorpot = Mock(return_value=('0', 'free'))
                result = writer.make_rotors(species, 1.)
                self.assertEqual(result, 'free rotor' if status == 0 else '')
                self.assertEqual(writer.make_rotorpot.call_count, int(status == 0))


if __name__ == '__main__':
    unittest.main()
