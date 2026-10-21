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

    def test_rotdpy_mess_selection_is_validated(self):
        for options, message in (
                ({'rotdpy_mess_energy_index': True}, 'energy_index'),
                ({'rotdpy_mess_energy_index': -2}, 'energy_index'),
                ({'rotdpy_directory': ''}, 'directory')):
            with self.subTest(options=options), self.assertRaisesRegex(
                    ValueError, message):
                self.read_parameters(**options)

    def test_full_pes_accepts_deferred_required_rotd_models(self):
        for mode in ('required', 'production'):
            with self.subTest(mode=mode):
                self.assertEqual(self.read_parameters(pes=1, rotdpy_mess_mode=mode)['rotdpy_mess_mode'], mode)
        self.read_parameters(pes=1, vrc_tst_surface_model_file='reviewed.json')
        self.assertEqual(self.read_parameters(pes=1)['rotdpy_mess_mode'], 'auto')

    def test_spin_splitting_is_additive(self):
        direct = self.read_parameters()
        self.assertEqual(direct['vrc_tst_correction_model'], 'direct')
        spin = self.read_parameters(
            vrc_tst_correction_model='spin_splitting')
        self.assertEqual(spin['vrc_tst_correction_model'], 'spin_splitting')

    def test_vrc_long_walltime_default_preserves_explicit_input(self):
        self.assertEqual(self.read_parameters()['vrc_tst_walltime'], '7-00:00:00')
        self.assertEqual(self.read_parameters(vrc_tst_walltime='02:00:00')['vrc_tst_walltime'],
                         '02:00:00')

    def test_anl1_f12_vrc_preset_requires_explicit_orbital_policy(self):
        with self.assertRaisesRegex(ValueError, 'not a chemical active space'):
            self.read_parameters(vrc_tst_accuracy_preset='anl1-f12')
        parameters = self.read_parameters(
            vrc_tst_accuracy_preset='anl1-f12',
            vrc_tst_active_space='full_valence')
        self.assertEqual(parameters['vrc_tst_active_space'], 'full_valence')
        self.assertEqual(parameters['vrc_tst_scan_method'],
                         'caspt2-f12(2,2)')
        self.assertEqual(parameters['vrc_tst_scan_gradient'], 'numerical')
        self.assertEqual(parameters['vrc_tst_sample_basis'],
                         'avdz-f12')
        self.assertEqual(parameters['vrc_tst_scan_shift'], 0.1)
        self.assertEqual(parameters['vrc_tst_scan_ipea'], 0.25)
        gap, = parameters['vrc_tst_spin_splitting']['gap_levels']
        self.assertEqual(gap['basis'], ['avdz-f12', 'avtz-f12'])
        self.assertEqual(gap['extrapolation_power'], 3.7)
        self.assertEqual(parameters['vrc_tst_correction_model'],
                         'spin_splitting')

    def test_full_valence_rejects_nonaugmented_multireference_basis_early(self):
        with self.assertRaisesRegex(ValueError, 'augmented orbital bases'):
            self.read_parameters(
                vrc_tst_active_space='full_valence',
                vrc_tst_scan_method='caspt2-f12(2,2)',
                vrc_tst_scan_basis='vqz-f12')

    def test_invalid_spin_splitting_override_fails_during_input_read(self):
        with self.assertRaisesRegex(ValueError, 'adjacent'):
            self.read_parameters(
                vrc_tst_correction_model='spin_splitting',
                vrc_tst_spin_splitting={
                    'triplet_reference': {
                        'basis': ['avdz-f12', 'avqz-f12']}})

    def test_automatic_cbh_and_reaction_selection_parameters(self):
        parameters = self.read_parameters(
            composite_method='ANL0-F12', cbh_rung=1,
            fc_model_path='/cached/uma-s-1p2.pt',
            stop_product_reaction_search=1,
            selected_reactions={
                'initial': {'hom_sci': [[0, 1]],
                            'h2_elim': [[3, 0, 1, 6]]}})
        self.assertEqual(parameters['cbh_rung'], 1)
        self.assertEqual(parameters['stop_product_reaction_search'], 1)

    def test_automatic_cbh_defaults_are_hands_off_and_share_queue_cap(self):
        parameters = self.read_parameters(
            composite_method='ANL0-F12', theory_preset='uma-b2plyp-anl',
            queue_job_limit=8)
        self.assertEqual(parameters['cbh_rung'], 'auto')
        self.assertEqual(parameters['vrc_tst_max_nodes'], 8)
        self.assertEqual(parameters['rotdpy_max_jobs'], 8)

        with TemporaryDirectory() as directory:
            input_file = Path(directory) / 'input.json'
            input_file.write_text(json.dumps({
                'barrier_threshold': 50.,
                'sella_kwargs': {'internal': True, 'gamma': 0.},
                'composite_method': 'ANL0-F12',
                'theory_preset': 'uma-b2plyp-anl',
                'queue_job_limit': 8,
            }))
            parsed = Parameters(str(input_file), show_warnings=False)
        self.assertEqual(parsed.theory_profiles['l1'].model_path,
                         'uma-s-1p2')

    def test_automatic_queue_subsystem_overrides_are_preserved(self):
        parameters = self.read_parameters(
            composite_method='ANL0-F12', theory_preset='uma-b2plyp-anl',
            queue_job_limit=8, vrc_tst_max_nodes=3, rotdpy_max_jobs=5)
        self.assertEqual(parameters['vrc_tst_max_nodes'], 3)
        self.assertEqual(parameters['rotdpy_max_jobs'], 5)

    def test_reaction_selection_fails_closed_on_invalid_input(self):
        for selected in ([], {'initial': {}},
                         {'initial': {'hom_sci': []}},
                         {'initial': {'hom_sci': [[0, True]]}},
                         {'initial': {'hom_sci': [[0, 1], [0, 1]]}}):
            with self.subTest(selected=selected), self.assertRaises(ValueError):
                self.read_parameters(selected_reactions=selected)
        with self.assertRaisesRegex(ValueError, 'requires composite_method'):
            self.read_parameters(cbh_rung=1)


if __name__ == '__main__':
    unittest.main()
