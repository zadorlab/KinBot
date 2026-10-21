"""IRC endpoint selection and L2 polling use the jobs that actually exist."""
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import numpy as np
import pytest
from ase import Atoms

from kinbot import constants
from kinbot.bimol import load_selected_species_properties
from kinbot.calculation import selected_calculation_job
from kinbot.optimize import Optimize
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.reaction_generator import ReactionGenerator
from kinbot.species_routing import routing_name
from kinbot.stationary_pt import StationaryPoint
from kinbot.thermochemistry import thermochemistry_evidence


class TestEndpointCalculationNames(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory()
        self.previous = Path.cwd()
        os.chdir(self.directory.name)
        Path('input.json').write_text('{"barrier_threshold": 50.0}')
        self.par = Parameters('input.json').par
        self.par.update(bimol=1, conformer_search=0, rotor_scan=0,
                        high_level=0, L3_calc=0, queuing='local')
        for method, result in (('submit_qc', 0), ('check_qc', 'normal')):
            patcher = patch.object(QuantumChemistry, method, return_value=result)
            setattr(self, method, patcher.start())
            self.addCleanup(patcher.stop)
        self.qc = QuantumChemistry(self.par)

    def tearDown(self):
        os.chdir(self.previous)
        self.directory.cleanup()

    def endpoint(self, direction):
        point = StationaryPoint(
            f'abstraction_final_IRC_{direction}_prod', 0, 1,
            atom=['H', 'H'], geom=np.array([[0., 0., 0.], [.74, 0., 0.]]))
        point.characterize()
        point.energy, point.zpe, point.freq = -1.02, .01, [4400.]
        self.record(point, point.name, point.energy)
        return point

    def record(self, point, job, energy):
        self.qc.db.write(Atoms(point.atom, positions=point.geom), name=job,
                         data={'status': 'normal',
                               'energy': energy / constants.EVtoHARTREE,
                               'zpe': .01, 'frequencies': [4400.]})

    def test_l1_reads_both_irc_endpoints_and_keeps_binding_depth_test(self):
        optimizations = []
        for direction in ('F', 'R'):
            endpoint = self.endpoint(direction)
            optimization = Optimize(endpoint, self.par, self.qc, just_high=True)
            optimization.do_optimization()
            self.assertEqual(optimization.shigh, 1)
            self.assertEqual(selected_calculation_job(optimization), endpoint.name)
            load_selected_species_properties(optimization, self.qc)
            self.assertEqual(endpoint.source_job, endpoint.name)
            self.assertEqual(self.qc.db.count(name=routing_name(endpoint)), 0)
            optimizations.append(optimization)
        self.submit_qc.assert_not_called()

        fragment = SimpleNamespace(energy=-.5, zpe=0.)
        generator = ReactionGenerator.__new__(ReactionGenerator)
        generator.par = self.par
        generator.species = SimpleNamespace(fragA=fragment, fragB=fragment)
        reaction = SimpleNamespace(
            instance_name='abstraction', do_reac_vdW=True, do_vdW=True,
            irc_reactant_opt=optimizations[0], irc_prod_opt=optimizations[1],
            prod_opt=[SimpleNamespace(species=fragment)] * 2)
        generator._finalize_bimol_complexes(reaction)
        self.assertTrue(reaction.do_reac_vdW)
        self.assertTrue(reaction.do_vdW)
        self.assertAlmostEqual(reaction.reac_vdW_depth, .01 * constants.AUtoKCAL)

    def test_l2_submits_polls_and_loads_the_same_endpoint_job(self):
        self.par['high_level'] = 1
        for direction in ('F', 'R'):
            endpoint = self.endpoint(direction)
            job = endpoint.name + '_high'
            optimization = Optimize(endpoint, self.par, self.qc, just_high=True)
            self.check_qc.return_value = 'running'
            optimization.do_optimization()
            self.assertEqual(optimization.shigh, 0)
            self.assertTrue(Path(job + '.py').is_file())
            self.assertEqual(self.submit_qc.call_args.args[0], job)
            self.assertEqual(self.check_qc.call_args.args[0], job)

            self.record(endpoint, job, -1.03)
            self.check_qc.return_value = 'normal'
            optimization.do_optimization()
            self.assertEqual(optimization.shigh, 1)
            self.assertEqual(selected_calculation_job(optimization), job)
            load_selected_species_properties(optimization, self.qc)
            self.assertEqual(endpoint.source_job, job)
            self.assertAlmostEqual(endpoint.energy, -1.03)
            forbidden = routing_name(endpoint) + '_high'
            self.assertNotIn(forbidden, [call.args[0] for call in
                                        self.check_qc.call_args_list])

    def test_complete_l2_irc_endpoints_do_not_repeat_optimization(self):
        self.par['high_level'] = 1
        endpoints = [self.endpoint(direction) for direction in ('F', 'R')]
        for endpoint in endpoints:
            self.record(endpoint, endpoint.name, -1.04)
            self.qc._record(endpoint.name, 'l2')
        reaction = SimpleNamespace(instance_name='abstraction',
            irc_reactant=endpoints[0], irc_prod=endpoints[1], prod_opt=[Mock(), Mock()])
        generator = ReactionGenerator(SimpleNamespace(), self.par, self.qc, 'input.json')
        with patch.object(self.qc, 'qc_opt') as optimize:
            generator._classify_bimol_complexes(reaction)
        optimize.assert_not_called()
        self.submit_qc.assert_not_called()
        for optimization, endpoint in zip(
                (reaction.irc_reactant_opt, reaction.irc_prod_opt), endpoints):
            load_selected_species_properties(optimization, self.qc)
            self.assertEqual(optimization.shigh, 1)
            self.assertEqual(endpoint.source_job, endpoint.name)
            self.assertEqual(self.qc.for_job(endpoint.source_job).theory_level, 'l2')
            self.assertAlmostEqual(endpoint.energy, -1.04)
            self.assertEqual(endpoint.freq, [4400.])
            self.assertEqual(self.qc.db.count(name=endpoint.name + '_high'), 0)
        self.assertEqual(self.par['high_level'], 1)

    def test_incomplete_or_l1_endpoint_still_requests_l2_properties(self):
        self.par['high_level'] = 1
        self.check_qc.return_value = 'running'
        endpoints = [self.endpoint(direction) for direction in ('F', 'R')]
        self.qc._record(endpoints[0].name, 'l1')
        self.qc._record(endpoints[1].name, 'l2')
        # Q-Chem's own IRC worker stores its endpoint energy only.
        self.qc.db.write(Atoms(endpoints[1].atom, positions=endpoints[1].geom),
            name=endpoints[1].name, data=dict(status='normal', energy=-1.04))
        reaction = SimpleNamespace(instance_name='abstraction',
            irc_reactant=endpoints[0], irc_prod=endpoints[1], prod_opt=[Mock(), Mock()])
        generator = ReactionGenerator(SimpleNamespace(), self.par, self.qc, 'input.json')
        with patch.object(self.qc, 'qc_opt', wraps=self.qc.qc_opt) as optimize:
            generator._classify_bimol_complexes(reaction)
        self.assertEqual(optimize.call_count, 2)
        self.assertTrue(all(call.kwargs['high_level'] == 1
                            for call in optimize.call_args_list))
        self.assertEqual(reaction.irc_reactant_opt.shigh, 0)
        self.assertEqual(reaction.irc_prod_opt.shigh, 0)

    def test_ordinary_well_and_ts_names_are_unchanged(self):
        point = self.endpoint('F')
        self.par['bimol'] = 0
        normal = Optimize(point, self.par, self.qc)
        self.assertEqual(normal.log_name(0), routing_name(point) + '_well')
        self.assertEqual(normal.log_name(1), routing_name(point) + '_well_high')
        complex_opt = Optimize(point, self.par, self.qc, just_high=True)
        self.assertEqual(complex_opt.log_name(1), routing_name(point) + '_high')
        point.wellorts = 1
        saddle = Optimize(point, self.par, self.qc)
        self.assertEqual(saddle.log_name(0), point.name)

    def test_ordinary_endpoint_keeps_its_conformer_representation(self):
        self.par.update(bimol=0, multi_conf_tst=1, conformer_search=1)
        point = self.endpoint('F')
        Optimize(point, self.par, self.qc, just_high=True)
        self.assertEqual(point.conformer_representation, 'MC-RRHO')

    def test_mc_run_exports_complexes_as_the_single_calculation_performed(self):
        self.par.update(multi_conf_tst=1, conformer_search=1)
        for high_level, direction in ((0, 'F'), (1, 'R')):
            self.par['high_level'] = high_level
            endpoint = self.endpoint(direction)
            if high_level:
                self.record(endpoint, endpoint.name + '_high', -1.03)
            optimization = Optimize(endpoint, self.par, self.qc, just_high=True)
            optimization.do_optimization()
            self.assertEqual(optimization.shigh, 1)
            load_selected_species_properties(optimization, self.qc)
            evidence = thermochemistry_evidence(endpoint)
            self.assertEqual(evidence['conformer_representation'], 'selected conformer')
            self.assertIsNone(evidence['conformer_counting_error'])
            self.assertEqual(evidence['optical_counting']['representation'], 'RRHO')
            self.assertEqual(evidence['optical_counting']['remaining_multiplier'], 1)


@pytest.mark.parametrize('converged', [None, False, True])
def test_fairchem_endpoint_reuse_requires_optimization_convergence(converged):
    from kinbot.bimol import has_complete_irc_endpoint
    data = dict(status='normal', energy=-1., zpe=.01, frequencies=[4400.])
    if converged is not None:
        data['endpoint_converged'] = converged
    qc = SimpleNamespace(for_job=lambda job: SimpleNamespace(theory_level='l2', qc='fc'),
        db=SimpleNamespace(select=lambda **kw: [SimpleNamespace(data=data)]))
    assert has_complete_irc_endpoint(qc, 'irc_prod', 'l2') is (converged is True)


@pytest.mark.parametrize('bimol', [False, True])
@pytest.mark.parametrize('converged', [False, True, 'error'])
def test_fairchem_worker_records_endpoint_convergence_without_rejecting_irc(
        tmp_path, monkeypatch, bimol, converged):
    import pickle
    import sys
    from types import ModuleType
    import ase.io
    import sella
    from kinbot import kb_path, fairchem_utils, frequencies

    monkeypatch.chdir(tmp_path)
    class Calculator:
        def __init__(self, *args, **kwargs):
            self.results = {}

        def get_potential_energy(self, atoms, **kwargs):
            self.results['forces'] = np.zeros((2, 3))
            return -1.

    class Optimizer:
        def __init__(self, *args, **kwargs):
            pass

        def run(self, **kwargs):
            if converged == 'error':
                raise RuntimeError('Endpoint optimization stopped')
            return converged

    module = ModuleType('fairchem.core')
    module.FAIRChemCalculator = Calculator
    monkeypatch.setitem(sys.modules, 'fairchem.core', module)
    monkeypatch.setattr(fairchem_utils, 'load_predictor', lambda *args, **kwargs: None)
    irc_run = Mock(return_value=True)
    monkeypatch.setattr(sella, 'IRC', lambda *args, **kwargs: SimpleNamespace(run=irc_run))
    monkeypatch.setattr(sella, 'Sella', Optimizer)
    monkeypatch.setattr(ase.io, 'read', lambda *args, **kwargs: [])
    monkeypatch.setattr(ase.io, 'write', lambda *args, **kwargs: None)
    monkeypatch.setattr(frequencies, 'calc_vibrations', lambda *args: ([4400.], .01, []))
    template = Path(kb_path, 'tpl/ase_fc_irc.tpl.py').read_text()
    script = template.format(label='irc_F', atom=['H', 'H'], geom=[[0., 0., 0.], [.74, 0., 0.]],
        kwargs=dict(charge=0, mult=1), prod_kwargs={}, fc_model_path='test-model',
        fc_device='cpu', fc_task_name='omol', sella_kwargs={}, fmax=.03, steps=10,
        bimol=bimol, fc_compile=False, irc_maxpoints=30 if bimol else 100)
    exec(compile(script, 'irc_F.py', 'exec'), {})
    irc_run.assert_called_once_with(fmax=.01, steps=30 if bimol else 100, direction='forward')
    with Path('irc_F_prod.pkl').open('rb') as handle:
        result = pickle.load(handle)['data']
    assert result['status'] == 'normal'
    if bimol:
        assert result['endpoint_converged'] is (converged is True)
    else:
        assert 'endpoint_converged' not in result


if __name__ == '__main__':
    unittest.main()
