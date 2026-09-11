"""Selected-source Hessians through real QC writers, databases and parsers."""
import ast
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock, patch

import numpy as np
from ase.build import molecule

from kinbot import constants
from kinbot.calculation import load_calculation_record
from kinbot.optimize import Optimize
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint


class TestSelectedHessian(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        Path('conf').mkdir()
        Path('input.json').write_text(json.dumps({'barrier_threshold': 100.,
            'method': 'HF', 'basis': 'sto-3g', 'high_level': 0,
            'high_level_method': 'MP2', 'high_level_basis': '6-31g',
            'rotor_scan': 0, 'conformer_search': 0, 'queuing': 'local'}))
        self.par = Parameters('input.json', show_warnings=False).par
        self.par['rotor_scan'] = 1
        self.qc = QuantumChemistry(self.par)
        self.qc.submit_qc = Mock(return_value=1)
        self.atoms = molecule('CH3OH')
        self.species = StationaryPoint('methanol', 0, 1,
            atom=self.atoms.get_chemical_symbols(), geom=self.atoms.positions)
        self.species.characterize()
        self.job = 'conf/selected_0001'
        self.freq = [100.] * (3 * len(self.atoms) - 6)
        self.hess = np.eye(3 * len(self.atoms)) * .02

    def record(self, job, hess=None):
        data = {'energy': -100. / constants.EVtoHARTREE, 'zpe': .01,
                'frequencies': self.freq, 'status': 'normal'}
        if hess is not None:
            data['hess'] = hess
        return self.qc.db.write(self.atoms, name=job, data=data)

    def write_hessian(self, job):
        if self.qc.qc == 'gauss':
            flat = self.hess[np.tril_indices(len(self.hess))]
            lines = ['Cartesian Force Constants R N= ' + str(len(flat))]
            lines += [' '.join(map(str, flat[i:i+5])) for i in range(0, len(flat), 5)]
            Path(job + '.fchk').write_text('\n'.join(lines) + '\n')
            Path(job + '.log').write_text('done\n')
        else:
            lines = [' Mass-Weighted Hessian Matrix']
            lines += [' '.join(map(str, row)) for row in self.hess]
            lines += [' Translations and Rotations']
            Path(job + '_freq.out').write_text('\n'.join(lines) + '\n')
            Path(job + '.out').write_text('done\n')

    def optimization(self):
        opt = Optimize(self.species, self.par, self.qc)
        opt.selected_job = self.job
        opt.shigh, opt.sconf = 1, 1
        return opt

    def test_conformers_drop_gaussian_checkpoints_but_keep_qchem_hessian_printing(self):
        for backend in ('gauss', 'qchem'):
            self.qc.qc = backend
            self.qc.qc_conf(self.species, self.species.geom, 1)
            job = f'conf/{self.species.chemid}_0001'
            tree = ast.parse(Path(job + '.py').read_text())
            kwargs = next(ast.literal_eval(node.value) for node in tree.body
                          if isinstance(node, ast.Assign) and any(
                              isinstance(t, ast.Name) and t.id == 'kwargs' for t in node.targets))
            if backend == 'gauss':
                self.assertNotIn('chk', kwargs)
            else:
                self.assertEqual(kwargs['vibman_print'], '4')
            self.assertEqual(kwargs['method'].lower(), 'hf')

    def test_loading_raw_thermochemistry_does_not_request_native_hessian(self):
        self.record(self.job)
        with patch.object(self.qc, 'read_qc_hess', side_effect=AssertionError('not needed')):
            load_calculation_record(self.species, self.qc, self.job)
        self.assertEqual(self.species.freq, self.freq)
        self.assertEqual(self.species.hess, [])

    def test_native_parsers_return_the_selected_hessian(self):
        for backend in ('gauss', 'qchem'):
            with self.subTest(backend=backend):
                self.qc.qc = backend
                self.record(self.job)
                self.write_hessian(self.job)
                np.testing.assert_allclose(self.qc.read_qc_hess(self.job, len(self.atoms)), self.hess)
                self.assertEqual(self.qc.hessian_is_massweighted(), backend == 'qchem')

    def test_legacy_missing_hessian_recovers_without_optimizing_or_using_parent(self):
        for backend in ('gauss', 'qchem'):
            with self.subTest(backend=backend):
                self.qc.qc = backend
                self.job = f'conf/{backend}_selected'
                source_id = self.record(self.job)
                Path(self.job + ('.log' if backend == 'gauss' else '.out')).write_text('done\n')
                # An available parent Hessian must not stand in for the selected source.
                self.record('parent', np.eye(len(self.hess)) * 99.)
                load_calculation_record(self.species, self.qc, self.job)
                opt = self.optimization()
                generated = Mock()
                hir = Mock(generate_hir_geoms=generated, check_hir=Mock(return_value=0))
                with patch('kinbot.optimize.HIR', return_value=hir):
                    opt.do_optimization()
                    generated.assert_not_called()
                    self.assertEqual(opt.shir, -1)
                    recovery = f'{self.job}_freq_recovery_{source_id}'
                    script = Path(recovery + '.py').read_text()
                    tree = ast.parse(script)
                    kwargs = (next(ast.literal_eval(node.value) for node in tree.body
                              if isinstance(node, ast.Assign) and any(
                                  isinstance(t, ast.Name) and t.id == 'kwargs' for t in node.targets))
                              if backend == 'gauss' else ast.literal_eval(next(
                                  node for node in ast.walk(tree) if isinstance(node, ast.Call)
                                  and isinstance(node.func, ast.Name) and node.func.id == 'QChem').keywords[0].value))
                    self.assertNotIn('opt', kwargs)
                    self.assertEqual(kwargs['method'].lower(), 'hf')
                    self.assertEqual(kwargs['basis'].lower(), 'sto-3g')
                    if backend == 'qchem':
                        self.assertEqual(kwargs['jobtype'], 'freq')
                    self.record(recovery)
                    self.write_hessian(recovery)
                    opt.do_optimization()
                generated.assert_called_once()
                self.assertEqual(opt.selected_job, recovery)
                self.assertEqual(self.species.source_job, recovery)
                np.testing.assert_allclose(self.species.hess, self.hess)
                np.testing.assert_allclose(self.species.geom, self.atoms.positions)
                self.assertEqual(self.species.freq, self.freq)
                self.assertEqual(self.qc.db.count(name=self.job), 1)
                self.assertEqual(self.qc.hessian_is_massweighted(), backend == 'qchem')

    def test_recovery_preserves_broken_symmetry_guess_options(self):
        self.qc.par['guessmix'] = 1
        for saddle in (0, 1):
            self.species.wellorts = saddle
            for extra, expected in (({}, 'Mix,Always'),
                                    ({'guess': 'Read,Mix,Always'}, 'Mix,Always'),
                                    ({'guess': '(read,Mix,Always)'}, 'Mix,Always'),
                                    ({'guess': 'Read'}, None)):
                with self.subTest(saddle=saddle, extra=extra):
                    self.record(self.job)
                    original = self.qc.get_qc_arguments
                    def arguments(*args, **kwargs):
                        return {**original(*args, **kwargs), **extra}
                    with patch.object(self.qc, 'get_qc_arguments', side_effect=arguments):
                        recovery = self.qc.qc_freq(self.species, self.job)
                    tree = ast.parse(Path(recovery + '.py').read_text())
                    kwargs = next(ast.literal_eval(node.value) for node in tree.body
                        if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name)
                            and target.id == 'kwargs' for target in node.targets))
                    self.assertEqual(kwargs.get('guess'), expected)
                    self.assertNotIn('opt', kwargs)

    def test_native_database_hessian_needs_no_checkpoint(self):
        self.record(self.job, self.hess)
        Path(self.job + '.log').write_text('done\n')
        np.testing.assert_allclose(self.qc.read_qc_hess(self.job, len(self.atoms)), self.hess)


if __name__ == '__main__':
    unittest.main()
