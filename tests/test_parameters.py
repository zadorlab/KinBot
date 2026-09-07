"""Backend-aware input validation regressions."""

from contextlib import redirect_stdout
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from kinbot.parameters import Parameters


class TestBackendWarnings(unittest.TestCase):
    def read_parameters(self, **options):
        with TemporaryDirectory() as directory:
            input_file = Path(directory) / 'input.json'
            input_file.write_text(json.dumps({
                'barrier_threshold': 50.,
                'sella_kwargs': {'internal': True, 'gamma': 0.},
                **options,
            }))
            return Parameters(str(input_file), show_warnings=True).par

    def test_fairchem_rotors_ignore_unused_method_and_basis(self):
        with patch('kinbot.parameters.logger.warning') as warning:
            par = self.read_parameters(qc='fc', fc_model_path='uma-m-1p1',
                                       rotor_scan=1)
        self.assertEqual(par['rotor_scan'], 1)
        warning.assert_not_called()

    def test_conventional_rotors_still_require_matching_levels(self):
        with self.assertRaises(SystemExit):
            self.read_parameters(qc='gauss', rotor_scan=1)
        with patch('kinbot.parameters.logger.warning') as warning:
            self.read_parameters(qc='gauss', rotor_scan=1,
                                 high_level_method='B3LYP', high_level_basis='6-31G')
        self.assertIn('B3LYP/6-31G', str(warning.call_args))

    def test_bimolecular_method_warning_uses_the_actual_backend(self):
        for backend in ['fc', 'gauss']:
            with self.subTest(backend=backend), redirect_stdout(io.StringIO()) as output:
                with patch('kinbot.parameters.logger.warning') as warning:
                    self.read_parameters(qc=backend, fc_model_path='uma-m-1p1',
                                         bimol=1, structure=[['H', 0., 0., 0.]] * 2)
                self.assertEqual('B3LYP is not recommended' in output.getvalue(),
                                 backend == 'gauss')
                self.assertEqual(warning.called, backend == 'gauss')


if __name__ == '__main__':
    unittest.main()
