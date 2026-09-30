"""Opt-in theory configuration and calculator factory regressions.

These tests run without Gaussian, Molpro, CFOUR, MRCC, or FairChem installed.
"""

import json
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from ase import Atoms

from kinbot.ase_modules.calculators.factory import (
    build_calculator, calculator_spec, capabilities,
)
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.fairchem_utils import load_predictor
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

    def test_active_profiles_route_l1_and_l2_without_loading_optional_codes(self):
        parameters = self.parameters(theory_preset='uma-b2plyp-anl',
                                     fc_model_path='/site/uma.pt')
        qc = QuantumChemistry(parameters.par)
        self.assertEqual(qc.qc, 'fc')
        self.assertEqual(qc.l1.qc, 'fc')
        self.assertEqual(qc.l2.qc, 'gauss')
        self.assertTrue(qc.l2.use_sella)
        self.assertEqual(qc.l2.qc_command, 'g16')
        self.assertEqual(qc.l2.get_qc_arguments(
            'job_high', 1, 0, 10, high_level=1)['method'], 'B2PLYP')
        self.assertEqual(qc.l2.par['calc_kwargs']['EmpiricalDispersion'],
                         'GD3BJ')
        self.assertEqual(qc.l2.frequency_mode, 'native_hessian')
        self.assertFalse(Path('kinbot.db').exists())

    def test_profiled_gaussian_sella_job_uses_one_native_frequency(self):
        parameters = self.parameters(theory_preset='uma-b2plyp-anl',
                                     fc_model_path='/site/uma.pt')
        qc = QuantumChemistry(parameters.par)
        species = SimpleNamespace(
            chemid=200000000000000000002, name='h2', mult=1, charge=0,
            nel=2, natom=2, wellorts=0, atom=['H', 'H'],
            geom=np.array([[0., 0., 0.], [0., 0., .74]]))
        with patch.object(qc.l2, 'submit_qc', return_value=1):
            qc.qc_opt(species, species.geom, high_level=1)
        generated = Path(f'{species.chemid}_well_high.py').read_text()
        compile(generated, f'{species.chemid}_well_high.py', 'exec')
        self.assertIn("frequency_mode = 'native_hessian'", generated)
        self.assertIn("frequency_kwargs['freq'] = ''", generated)
        self.assertIn("frequency_kwargs.pop('opt', None)", generated)

        Path('conf').mkdir()
        with patch.object(qc.l1, 'submit_qc', return_value=1):
            qc.qc_conf(species, species.geom, 0)
        generated_conf = Path(
            f'conf/{species.chemid}_0000.py').read_text()
        compile(generated_conf, 'profiled_conformer.py', 'exec')
        self.assertIn('calc_vibrations', generated_conf)

        with patch.object(qc.l2, 'submit_qc', return_value=1):
            qc.l2.qc_conf(species, species.geom, 0)
        generated_l2_conf = Path(
            f'conf/{species.chemid}_0000.py').read_text()
        compile(generated_l2_conf, 'profiled_l2_conformer.py', 'exec')
        self.assertIn("frequency_mode = 'native_hessian'",
                      generated_l2_conf)
        self.assertIn("'chk': '200000000000000000002_0000'",
                      generated_l2_conf)
        self.assertIn(
            "logfile='conf/200000000000000000002_0000_sella.log'",
            generated_l2_conf)
        self.assertIn(
            "frequency_kwargs['chk'] = "
            "os.path.basename('conf/200000000000000000002_0000')",
            generated_l2_conf)

        species.wellorts = 1
        species.name = 'h2_saddle'
        with patch.object(qc.l2, 'submit_qc', return_value=1):
            qc.qc_opt_ts(species, species.geom, high_level=1)
        generated_ts = Path('h2_saddle_high.py').read_text()
        compile(generated_ts, 'h2_saddle_high.py', 'exec')
        self.assertIn("frequency_mode = 'native_hessian'", generated_ts)
        self.assertIn("frequency_kwargs['freq'] = ''", generated_ts)
        self.assertIn("frequency_kwargs.pop('opt', None)", generated_ts)

        Path('hir').mkdir()
        rotor = SimpleNamespace(
            chemid=200000000000000000003, name='rotor', mult=1, charge=0,
            nel=18, natom=4, wellorts=0,
            atom=['C', 'C', 'H', 'H'],
            geom=np.array([[0., 0., 0.], [1.5, 0., 0.],
                           [-.5, 1., 0.], [2., 1., 1.]]))
        with patch.object(qc.l2, 'submit_qc', return_value=1):
            qc.qc_hir(rotor, rotor.geom, 0, 0, [[3, 1, 2, 4]], 0)
        generated_hir = Path(
            f'hir/{rotor.chemid}_hir_0_00.py').read_text()
        compile(generated_hir, 'profiled_l2_hir.py', 'exec')
        self.assertIn(
            f"'chk': '{rotor.chemid}_hir_0_00'", generated_hir)
        self.assertIn(
            f"logfile='hir/{rotor.chemid}_hir_0_00_sella.log'",
            generated_hir)
        self.assertIn('Constrained optimization failed:', generated_hir)

    def test_profiled_job_routes_and_scheduler_ids_survive_restart(self):
        parameters = self.parameters(theory_preset='uma-b2plyp-anl',
                                     fc_model_path='uma-s-1p2')
        qc = QuantumChemistry(parameters.par)
        qc._record('conf/a', 'l1', '123')
        qc._record('a_well_high', 'l2', '456')
        restarted = QuantumChemistry(parameters.par)
        self.assertIs(restarted._backend_for('conf/a'), restarted.l1)
        self.assertIs(restarted._backend_for('a_well_high'), restarted.l2)
        self.assertEqual(restarted.l1.job_ids['conf/a'], '123')
        self.assertEqual(restarted.l2.job_ids['a_well_high'], '456')
        self.assertEqual(restarted.job_ids['a_well_high'], '456')
        restarted._clear_job_id('a_well_high')
        final = QuantumChemistry(parameters.par)
        self.assertNotIn('a_well_high', final.job_ids)
        self.assertNotIn('a_well_high', final.l2.job_ids)

    def test_profiled_vrc_jobs_use_gaussian_context(self):
        parameters = self.parameters(theory_preset='uma-b2plyp-anl',
                                     fc_model_path='/site/uma.pt',
                                     queuing='local')
        qc = QuantumChemistry(parameters.par)
        fragment = SimpleNamespace(chemid=123)
        with patch.object(qc.l2, 'qc_vts_frag',
                          return_value='vrctst/123_vts') as submit:
            job = qc.qc_vts_frag(fragment)
        self.assertEqual(job, 'vrctst/123_vts')
        submit.assert_called_once_with(fragment)
        self.assertEqual(qc._routes[job]['level'], 'l2')

        # VRC templates always finish their native Gaussian ``.log`` file,
        # including when the L2 profile otherwise uses Sella.  Polling for a
        # ``_sella.log`` leaves a successfully completed VRC job waiting
        # forever.
        Path('vrctst').mkdir()
        qc.l2.db.write(Atoms('H'), name=job,
                       data={'energy': -1., 'status': 'normal'})
        Path(job + '.log').write_text('done\n')
        self.assertFalse(Path(job + '_sella.log').exists())
        self.assertEqual(qc.check_qc(job), 'normal')

    def test_anl1_f12_ladder_request_selects_higher_l2_surface(self):
        parameters = self.parameters(composite_method='ANL1-F12',
                                     fc_model_path='uma-s-1p2')
        self.assertEqual(parameters.theory_profiles['l2'].method, 'B2PLYP')


class TestCalculatorFactory(unittest.TestCase):
    def test_offline_fairchem_cache_ignores_invalid_generic_proxy(self):
        observed = {}
        expected = object()

        def load(*args, **kwargs):
            observed['ALL_PROXY'] = os.environ.get('ALL_PROXY')
            observed['all_proxy'] = os.environ.get('all_proxy')
            return expected

        fairchem = ModuleType('fairchem')
        core = ModuleType('fairchem.core')
        core.pretrained_mlip = SimpleNamespace(get_predict_unit=load)
        fairchem.core = core
        environment = {
            'HF_HUB_OFFLINE': '1',
            'ALL_PROXY': 'socks://proxy.example:80',
            'all_proxy': 'socks://proxy.example:80',
        }
        with patch.dict(sys.modules, {
                'fairchem': fairchem, 'fairchem.core': core}):
            with patch.dict(os.environ, environment, clear=True):
                self.assertIs(load_predictor('uma-s-1p2', 'cpu'), expected)
                self.assertEqual(os.environ['ALL_PROXY'],
                                 environment['ALL_PROXY'])
                self.assertEqual(os.environ['all_proxy'],
                                 environment['all_proxy'])
        self.assertEqual(observed, {'ALL_PROXY': None, 'all_proxy': None})

    def test_gated_fairchem_model_has_actionable_error(self):
        class GatedRepoError(Exception):
            pass

        fairchem = ModuleType('fairchem')
        core = ModuleType('fairchem.core')
        core.pretrained_mlip = SimpleNamespace(
            get_predict_unit=lambda *args, **kwargs: (_ for _ in ()).throw(
                GatedRepoError('403 Forbidden')))
        fairchem.core = core
        with patch.dict(sys.modules, {
                'fairchem': fairchem, 'fairchem.core': core}):
            with self.assertRaisesRegex(
                    RuntimeError, 'request or accept access.*facebook/UMA'):
                load_predictor('uma-s-1p2', 'cpu')

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
