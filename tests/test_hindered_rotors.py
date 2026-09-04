"""Regression tests for hindered-rotor results, without QC calculations."""

import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from ase.io import read
import numpy as np

from kinbot.hindered_rotors import HIR


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


if __name__ == '__main__':
    unittest.main()
