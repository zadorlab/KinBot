"""Failure-path regressions for ASE vibration calculations."""

import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import numpy as np

from kinbot import frequencies


class TestVibrationWorkingDirectory(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        self.directory = Path.cwd()
        self.molecule = SimpleNamespace(calc=SimpleNamespace(parameters={}))

    def test_failure_restores_working_directory(self):
        for stage in ('construction', 'run', 'write_jmol'):
            with self.subTest(stage=stage):
                os.chdir(self.directory)
                vibrations = Mock()
                constructor = Mock(return_value=vibrations)
                target = constructor if stage == 'construction' else getattr(vibrations, stage)
                target.side_effect = RuntimeError('deliberate vibration failure')
                with patch.object(frequencies, 'Vibrations', constructor):
                    result = frequencies.calc_vibrations(self.molecule, stage)
                self.assertEqual(result, (None, None, None))
                self.assertEqual(Path.cwd(), self.directory)

    def test_success_restores_directory_and_preserves_results(self):
        vibrations = Mock(H=np.eye(9))
        vibrations.get_zero_point_energy.return_value = 0.5
        point = Mock(geom=np.zeros((3, 3)))
        with patch.object(frequencies, 'Vibrations', return_value=vibrations), \
                patch.object(frequencies.StationaryPoint, 'from_ase_atoms', return_value=point), \
                patch.object(frequencies, 'get_frequencies', return_value=([100., 200., 300.], [])):
            freq, zpe, hess = frequencies.calc_vibrations(self.molecule, 'success')
        self.assertEqual(Path.cwd(), self.directory)
        self.assertEqual(freq, [100., 200., 300.])
        self.assertEqual(zpe, 0.5 * frequencies.EVtoHARTREE)
        np.testing.assert_allclose(hess, np.eye(9) / 97.17370087)


if __name__ == '__main__':
    unittest.main()
