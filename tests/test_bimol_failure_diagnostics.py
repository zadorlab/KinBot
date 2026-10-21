"""Keep rejected TS evidence and report its cause without changing retries."""
import pickle
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from kinbot import bimol, frequencies, kb_path, reaction_generator
from kinbot.reactions.reac_bimol import BimolecularReaction


@pytest.mark.parametrize('bimolecular', [False, True])
@pytest.mark.parametrize('convergence,spectrum', [
    ([True], [-600., -80., 200.]),
    ([True], [100., 200., 300.]),
    ([True], [-600., -20., 200.]),
    ([False, True], [-600., -80., 200.]),
    ([False, True], [-600., 100., 200.]),
    ([False, False], []),
])
def test_fairchem_worker_retains_rejected_spectrum_only_for_bimol(
        tmp_path, monkeypatch, bimolecular, convergence, spectrum):
    monkeypatch.chdir(tmp_path)

    class CalculatorStub(Calculator):
        implemented_properties = ['energy', 'forces']

        def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            self.results.update(energy=-10., forces=np.zeros((len(atoms), 3)))

    attempts = iter(convergence)
    calls = []

    class OptimizerStub:
        def __init__(self, atoms, **kwargs):
            self.atoms = atoms

        def run(self, **kwargs):
            calls.append(kwargs)
            self.atoms.calc.calculate(self.atoms)
            return next(attempts)

    core = ModuleType('fairchem.core')
    core.FAIRChemCalculator = lambda *args, **kwargs: CalculatorStub()
    package = ModuleType('fairchem')
    package.core = core
    monkeypatch.setitem(sys.modules, 'fairchem', package)
    monkeypatch.setitem(sys.modules, 'fairchem.core', core)
    monkeypatch.setattr('kinbot.fairchem_utils.load_predictor', lambda *args, **kwargs: None)
    monkeypatch.setattr('sella.Sella', OptimizerStub)
    monkeypatch.setattr('ase.io.read', lambda *args, **kwargs: [Atoms('H3')])
    vibration = Mock(return_value=(np.asarray(spectrum), .01, np.eye(9)))
    monkeypatch.setattr(frequencies, 'calc_vibrations', vibration)
    script = Path(kb_path, 'tpl/ase_fc_ts_end.tpl.py').read_text().format(
        label='candidate', atom=['H'] * 3, geom=[[-1., 0., 0.], [0., 0., 0.], [1., 0., 0.]],
        kwargs={'charge': 0, 'mult': 2}, fc_model_path='test-model', fc_device='cpu',
        fc_compile=False, fc_task_name='omol', sella_kwargs={'internal': False},
        fmax=.01, steps=10, bimol=bimolecular)
    exec(compile(script, 'candidate.py', 'exec'), {})
    with Path('candidate.pkl').open('rb') as handle:
        data = pickle.load(handle)['data']
    accepted = bool(spectrum and spectrum[0] < 0 and spectrum[1] >= -50.)
    assert data['status'] == ('normal' if accepted else 'error')
    assert len(calls) == len(convergence)
    assert vibration.call_count == int(convergence[-1])
    assert Path('candidate_sella.log').read_text().endswith('done\n')
    if bimolecular:
        assert data['converged'] is convergence[-1]
        np.testing.assert_array_equal(data['frequencies'], spectrum)
        if not accepted:
            assert data['error'] == ('TS frequency check failed.' if convergence[-1]
                                      else 'TS optimization did not converge.')
    elif not accepted:
        assert data == {'status': 'error'}
    else:
        assert 'converged' not in data and 'error' not in data


@pytest.mark.parametrize('restart', [False, True])
@pytest.mark.parametrize('bimolecular', [False, True])
@pytest.mark.parametrize('record', [
    {'status': 'error', 'error': 'TS frequency check failed.', 'frequencies': [-600., -80., 200.]},
    {'status': 'error', 'error': 'TS optimization did not converge.', 'frequencies': []},
    {'status': 'error', 'frequencies': None},
    {'status': 'error'},
])
def test_controller_reports_current_error_on_failure_and_restart(
        tmp_path, monkeypatch, caplog, restart, bimolecular, record):
    monkeypatch.chdir(tmp_path)
    qc = SimpleNamespace(check_qc=lambda name: 'error',
        db=SimpleNamespace(select=lambda **kwargs: iter([
            SimpleNamespace(data={'status': 'error', 'error': 'obsolete reason'}),
            SimpleNamespace(data=record)])))
    reaction = BimolecularReaction.__new__(BimolecularReaction)
    reaction.instance_name, reaction.qc = 'ethane_oh_candidate', qc
    species = SimpleNamespace(reac_obj=[reaction], reac_inst=[[0, 1, 2]],
        reac_step=[0 if restart else 1], reac_ts_done=[0], reac_type=['abstraction'])
    par = dict(bimol=bimolecular, skip_reactions=[], delete_intermediate_files=0,
               qc_command='', L3_calc=0)
    monkeypatch.setattr(bimol, 'prepare_bimolecular_reactions', lambda *args: None)
    monkeypatch.setattr(bimol, 'release_final_ts_clusters', lambda *args: None)
    monkeypatch.setattr(bimol, 'update_bimol_status', lambda *args, **kwargs: None)
    retry = Mock(return_value=False)
    monkeypatch.setattr(bimol, 'advance_pose_searches', retry)
    monkeypatch.setattr('kinbot.bimol_network.BimolecularNetwork', Mock())
    monkeypatch.setattr(reaction_generator.time, 'sleep', lambda *args: None)
    # Ordinary geometry failures still follow their existing message.
    if not bimolecular:
        monkeypatch.setattr(reaction_generator.reac_family, 'carry_out_reaction', lambda *args, **kwargs: -1)
    with caplog.at_level('INFO', logger='KinBot'):
        reaction_generator.ReactionGenerator(species, par, qc, 'input.json').generate()
    assert species.reac_ts_done == [-999]
    if bimolecular:
        retry.assert_called_once_with(species, par, qc)
        assert record.get('error', 'no reason was saved') in caplog.text
        assert 'Invalid geometry' not in caplog.text and 'obsolete reason' not in caplog.text
        if record.get('frequencies'):
            assert '-600.000, -80.000, 200.000' in caplog.text
    else:
        retry.assert_not_called()
        assert ('had error status on previous round' if restart else 'Invalid geometry') in caplog.text
