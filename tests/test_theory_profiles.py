"""Opt-in theory configuration and calculator factory regressions.

These tests run without Gaussian, Molpro, CFOUR, MRCC, or FairChem installed.
"""

import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from ase import Atoms

from kinbot.ase_modules.calculators.factory import (
    build_calculator, calculator_spec, capabilities,
)
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.theory import TheoryProfile


class TestTheoryProfiles(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)

    def parameters(self, **options):
        Path('input.json').write_text(json.dumps({
            'barrier_threshold': 50., **options,
        }))
        return Parameters('input.json', show_warnings=False)

    def test_legacy_input_has_no_profile_and_keeps_gaussian_route(self):
        parameters = self.parameters()
        self.assertEqual(parameters.par['profiled_theory'], 0)
        self.assertEqual(parameters.theory_profiles, {})
        qc = QuantumChemistry(parameters.par)
        self.assertEqual(qc.qc, 'gauss')
        self.assertEqual(qc.get_qc_arguments('job', 1, 0, 2)['method'], 'b3lyp')
        self.assertEqual(qc.get_qc_arguments('job', 1, 0, 2,
                                             high_level=1)['method'], 'M062X')

    def test_preset_resolves_independent_l1_and_l2(self):
        parameters = self.parameters(theory_preset='uma-b2plyp-anl',
                                     fc_model_path='/site/uma.pt')
        self.assertEqual(parameters.par['profiled_theory'], 1)
        l1, l2 = (parameters.theory_profiles[level] for level in ('l1', 'l2'))
        self.assertEqual((l1.calculator, l1.model_path, l1.command),
                         ('fairchem', '/site/uma.pt', ''))
        self.assertEqual((l2.calculator, l2.method, l2.basis),
                         ('gaussian', 'B2PLYP', 'cc-pVTZ'))
        self.assertEqual(l2.calculator_kwargs['EmpiricalDispersion'], 'GD3BJ')
        self.assertEqual(l2.calculator_kwargs['integral'], 'UltraFine')
        self.assertEqual(l2.frequency_mode, 'native_hessian')

    def test_lower_and_higher_anl_tiers_select_matching_l2_surfaces(self):
        lower = self.parameters(composite_method='ANL0-F12',
                                fc_model_path='/site/uma.pt')
        self.assertEqual(lower.theory_profiles['l2'].method, 'B3LYP')
        self.assertNotIn('EmpiricalDispersion',
                         lower.theory_profiles['l2'].calculator_kwargs)
        higher = self.parameters(composite_method='ANL1',
                                 fc_model_path='/site/uma.pt')
        self.assertEqual(higher.theory_profiles['l2'].method, 'B2PLYP')
        self.assertEqual(higher.theory_profiles['l2'].calculator_kwargs[
            'EmpiricalDispersion'], 'GD3BJ')

    def test_explicit_profiles_override_preset_without_mutating_it(self):
        parameters = self.parameters(
            theory_preset='uma-b2plyp-anl', fc_model_path='/site/uma.pt',
            l1_profile={'calculator': 'gaussian', 'method': 'wB97XD',
                        'basis': '6-31+G(d,p)'},
            l2_profile={'basis': 'cc-pVQZ',
                        'calculator_kwargs': {'scf': 'tight'}},
        )
        self.assertEqual(parameters.theory_profiles['l1'].calculator, 'gaussian')
        self.assertEqual(parameters.theory_profiles['l1'].method, 'wB97XD')
        self.assertEqual(parameters.theory_profiles['l1'].calculator_kwargs, {})
        self.assertEqual(parameters.theory_profiles['l2'].basis, 'cc-pVQZ')
        self.assertEqual(parameters.theory_profiles['l2'].calculator_kwargs,
                         {'EmpiricalDispersion': 'GD3BJ', 'Symm': 'None',
                          'scf': 'tight', 'integral': 'UltraFine'})
        fresh = self.parameters(theory_preset='uma-b2plyp-anl',
                                fc_model_path='/site/uma.pt')
        self.assertEqual(fresh.theory_profiles['l2'].calculator_kwargs['scf'], 'xqc')

    def test_switching_backend_drops_stale_preset_method_and_frequency_policy(self):
        parameters = self.parameters(
            theory_preset='uma-b2plyp-anl', fc_model_path='/site/uma.pt',
            l1='gaussian', l2_profile={'calculator': 'qchem'})
        l1 = parameters.theory_profiles['l1']
        l2 = parameters.theory_profiles['l2']
        self.assertEqual((l1.method, l1.basis, l1.frequency_mode),
                         ('b3lyp', '6-31G', 'auto'))
        self.assertEqual((l2.method, l2.basis, l2.frequency_mode),
                         ('M062X', '6-311++G(d,p)', 'auto'))
        self.assertEqual(l2.calculator_kwargs, {})
        self.assertEqual(l2.command, '')

    def test_composite_method_selects_preset_and_requires_checkpoint(self):
        with self.assertRaisesRegex(ValueError, 'model_path'):
            self.parameters(composite_method='ANL0-F12')
        parameters = self.parameters(composite_method='ANL0-F12',
                                     fc_model_path='/site/uma.pt')
        self.assertEqual(parameters.theory_profiles['l1'].calculator, 'fairchem')
        self.assertEqual(parameters.theory_profiles['l2'].method, 'B3LYP')

    def test_invalid_profile_and_preset_fail_during_input_parsing(self):
        with self.assertRaisesRegex(ValueError, 'Unknown theory_preset'):
            self.parameters(theory_preset='unknown')
        with self.assertRaisesRegex(ValueError, 'frequency_mode'):
            self.parameters(profiled_theory=1,
                            l1_profile={'calculator': 'gaussian',
                                        'frequency_mode': 'imaginary'})
        with self.assertRaisesRegex(ValueError, 'requires composite_method'):
            self.parameters(l3_overrides={'dboc': {'basis': 'cc-pVTZ'}})

    def test_active_profiles_fail_closed_until_job_router_is_connected(self):
        parameters = self.parameters(theory_preset='uma-b2plyp-anl',
                                     fc_model_path='/site/uma.pt')
        with self.assertRaisesRegex(NotImplementedError, 'router'):
            QuantumChemistry(parameters.par)
        raw = dict(parameters.par, profiled_theory=0,
                   l1_profile={'calculator': 'gaussian'})
        with self.assertRaisesRegex(NotImplementedError, 'router'):
            QuantumChemistry(raw)
        self.assertFalse(Path('kinbot.db').exists())


class TestCalculatorFactory(unittest.TestCase):
    def test_capabilities_are_explicit_and_fairchem_is_lazy(self):
        self.assertEqual(calculator_spec('gauss'), calculator_spec('gaussian'))
        self.assertTrue(capabilities('fairchem').forces)
        self.assertFalse(capabilities('gaussian').native_hessian)
        self.assertTrue(capabilities('molpro').forces)
        self.assertTrue(capabilities('molpro').numerical_forces)

    def test_gaussian_builder_scopes_label_and_command_to_one_job(self):
        profile = TheoryProfile(name='l2', calculator='gaussian',
                                method='B2PLYP', basis='cc-pVTZ', command='g16',
                                calculator_kwargs={'EmpiricalDispersion': 'GD3BJ'})

        class StubCalculator:
            def __init__(self, **kwargs):
                self.kwargs = kwargs

        with TemporaryDirectory() as directory:
            with patch('kinbot.ase_modules.calculators.factory.calculator_class',
                       return_value=StubCalculator):
                calc = build_calculator(profile, directory,
                                        {'name': 'water', 'charge': -1, 'mult': 2})
            self.assertEqual(calc.kwargs['label'], str(Path(directory) / 'water'))
            self.assertEqual(calc.kwargs['charge'], -1)
            self.assertEqual(calc.kwargs['mult'], 2)
            self.assertEqual(calc.kwargs['EmpiricalDispersion'], 'GD3BJ')
            self.assertEqual(calc.command, 'g16 < PREFIX.com > PREFIX.log')

    def test_job_name_cannot_escape_directory(self):
        with TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, 'basename'):
                build_calculator({'calculator': 'gaussian'}, directory,
                                 {'name': '../outside'})

    def test_gaussian_preset_renders_method_basis_dispersion_and_spin(self):
        profile = TheoryProfile(
            name='l2', calculator='gaussian', method='B2PLYP',
            basis='cc-pVTZ', calculator_kwargs={'EmpiricalDispersion': 'GD3BJ'},
        )
        with TemporaryDirectory() as directory:
            calc = build_calculator(profile, directory,
                                    {'name': 'h2', 'charge': 0, 'mult': 1})
            atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
            calc.write_input(atoms, ['energy', 'forces'])
            generated = (Path(directory) / 'h2.com').read_text()
        self.assertIn('B2PLYP/cc-pVTZ', generated)
        self.assertIn('EmpiricalDispersion(GD3BJ)', generated)
        self.assertIn('force', generated.lower())
        self.assertIn('\n0 1\n', generated)


if __name__ == '__main__':
    unittest.main()
