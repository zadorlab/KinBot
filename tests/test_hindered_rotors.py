"""Regression tests for hindered-rotor results, without QC calculations."""

import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from ase.io import read
import numpy as np

from kinbot.hindered_rotors import HIR
from kinbot import constants
from kinbot.mess import MESS


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


if __name__ == '__main__':
    unittest.main()
