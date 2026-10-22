"""Nested Gaussian checkpoints, endpoint retry data and native IRC evidence."""
import ast
import json
from pathlib import Path
from unittest.mock import Mock

import numpy as np
import pytest
from ase import Atoms
from ase.io import read

from kinbot import reader_gauss
from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint


def test_generated_conformer_checkpoint_is_valid_from_gaussian_working_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('input.json').write_text(json.dumps({'barrier_threshold': 100.,
        'method': 'HF', 'basis': 'sto-3g', 'use_sella': 0}))
    qc = QuantumChemistry(Parameters('input.json', show_warnings=False).par)
    qc.par['bimol'] = 1
    qc.submit_qc = Mock()
    point = StationaryPoint('ts', 0, 2, atom=['H', 'H', 'H'],
                             geom=np.array([[0., 0., 0.], [.9, 0., 0.], [1.8, 0., 0.]]), wellorts=1)
    point.characterize()
    Path('conf').mkdir(exist_ok=True)
    qc.qc_conf(point, point.geom, 0)
    tree = ast.parse(Path('conf/ts_0000.py').read_text())
    kwargs = next(ast.literal_eval(item.value) for item in tree.body if isinstance(item, ast.Assign)
                  and isinstance(item.targets[0], ast.Name) and item.targets[0].id == 'kwargs')
    assert kwargs['chk'] == 'ts_0000'
    calc = Gaussian(**kwargs)
    calc.write_input(Atoms(point.atom, positions=point.geom), properties=['energy'])
    assert Path(calc.directory).resolve() == tmp_path / 'conf'
    assert '%chk=ts_0000' in Path('conf/ts_0000.com').read_text()
    assert (Path(calc.directory) / kwargs['chk']).resolve() == tmp_path / 'conf/ts_0000'


def test_failed_log_without_geometry_or_optimization_table_keeps_input(tmp_path):
    log = tmp_path / 'early_failure.log'
    log.write_text('Gaussian failed before the first SCF\n')
    atoms = Atoms('H2O', positions=[[1., 0., 0.], [0., 1., 0.], [0., 0., 0.]])
    np.testing.assert_array_equal(reader_gauss.read_geom(log, atoms, keep_input=True), atoms.positions)
    _, geometry = reader_gauss.read_bimol_retry_geom_energy(log, atoms)
    np.testing.assert_array_equal(geometry, atoms.positions)
    np.testing.assert_array_equal(reader_gauss.read_geom(log, atoms), np.zeros((3, 3)))
    assert reader_gauss.read_lowest_geom_energy(log, atoms) is None
    # Truncated optimization tables must have the same safe recovery.
    log.write_text('Item               Value     Threshold  Converged?\n')
    _, geometry = reader_gauss.read_bimol_retry_geom_energy(log, atoms)
    np.testing.assert_array_equal(geometry, atoms.positions)


def test_gaussian_retry_does_not_repeat_cartesian(tmp_path):
    log = tmp_path / 'failed.log'
    log.write_text('Error termination request processed by link 9999.\n')
    kwargs = dict(opt='CalcFC', freq='')
    for _ in range(3):
        kwargs = reader_gauss.correct_kwargs(str(log), kwargs, bimol=True)
    assert kwargs['opt'] == 'CalcAll,Cartesian'


def test_native_irc_export_uses_accepted_points_in_file_order(tmp_path):
    log = tmp_path / 'irc.log'
    def point(z):
        return ('CURRENT STRUCTURE\nheader\nheader\nheader\nheader\nheader\n'
                f' 1 1 0.0 0.0 {z}\n 2 1 0.0 0.0 1.5\n -------\n')
    log.write_text('Charge = 0 Multiplicity = 1\n' + point(0.) + point(.1))
    destination = tmp_path / 'irc.xyz'
    assert reader_gauss.write_irc_trajectory(log, destination)
    frames = read(destination, index=':')
    assert len(frames) == 2
    assert frames[0].positions[0, 2] == 0.
    assert frames[1].positions[0, 2] == .1


def test_incomplete_native_irc_table_does_not_interrupt_endpoint_work(tmp_path):
    log = tmp_path / 'irc.log'
    log.write_text('CURRENT STRUCTURE\ninterrupted\n')
    assert not reader_gauss.write_irc_trajectory(log, tmp_path / 'irc.xyz')
    assert not (tmp_path / 'irc.xyz').exists()


def _irc_point(z):
    return ('CURRENT STRUCTURE\nheader\nheader\nheader\nheader\nheader\n'
            f' 1 1 0.0 0.0 {z}\n 2 1 0.0 0.0 1.5\n -------\n')


def test_partial_irc_keeps_complete_points_before_truncated_table(tmp_path):
    log = tmp_path / 'irc.log'
    log.write_text(_irc_point(.1) + 'CURRENT STRUCTURE\ntruncated\n')
    _, frames, _, _ = reader_gauss.read_all_irc_geoms(log)
    assert frames.shape == (1, 2, 3)
    assert frames[-1, 0, 2] == .1


def _run_irc_worker(tmp_path, monkeypatch, outputs, bimol=True):
    from kinbot import kb_path
    import ase.db
    import kinbot.ase_modules.calculators.gaussian as module
    from kinbot import utils
    monkeypatch.chdir(tmp_path)
    calls, rows = [], []

    class FakeGaussian:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

        def get_potential_energy(self, atoms, **kwargs):
            label = self.kwargs['label']
            calls.append((self.kwargs.copy(), atoms.positions.copy()))
            if label == 'irc':
                output, failure = outputs.pop(0)
                Path('irc.log').write_text(output)
                Path('irc.com').write_text(str(self.kwargs))
                if failure:
                    raise RuntimeError('Gaussian corrector failed')
            else:
                Path(label + '.log').write_text('Normal termination\n')
            return -1.

    class Database:
        def write(self, mol, **kwargs):
            rows.append((mol.positions.copy(), kwargs))

    monkeypatch.setattr(module, 'Gaussian', FakeGaussian)
    monkeypatch.setattr(ase.db, 'connect', lambda *a, **k: Database())
    monkeypatch.setattr(utils, 'iowait', lambda *a: None)
    monkeypatch.setattr(reader_gauss, 'read_freq', lambda *a: [100.])
    monkeypatch.setattr(reader_gauss, 'read_zpe', lambda *a: .01)
    if not bimol:
        monkeypatch.setattr(reader_gauss, 'read_geom', lambda filename, atoms, **kwargs: atoms.positions.copy())
        monkeypatch.setattr(reader_gauss, 'read_all_irc_geoms', Mock(side_effect=AssertionError(
            'Ordinary Gaussian IRCs must retain their previous endpoint recovery.')))
    template = Path(kb_path, 'tpl', 'ase_gauss_irc.tpl.py').read_text()
    script = template.format(working_dir=tmp_path, label='irc', atom=['H', 'H'], bimol=bimol,
        geom=[[0., 0., 0.], [0., 0., 1.5]], qc_command='unused',
        kwargs=dict(label='irc', charge=0, mult=1, geom='AllCheck', guess='Read',
                    irc='RCFC,reverse,MaxPoints=60,StepSize=20,MaxCycle=50'),
        prod_kwargs=dict(label='irc_prod'))
    exec(compile(script, 'irc.py', 'exec'), {})
    return calls, rows


def test_no_accepted_irc_points_retry_from_ts_and_do_not_optimize_endpoint(tmp_path, monkeypatch):
    calls, rows = _run_irc_worker(tmp_path, monkeypatch,
        [('Input orientation: corrector trial\nError termination\n', True)] * 2)
    assert len(calls) == 2
    np.testing.assert_array_equal(calls[0][1], calls[1][1])
    assert 'CalcFC' in calls[1][0]['irc'] and 'RCFC' not in calls[1][0]['irc']
    assert 'StepSize=10' in calls[1][0]['irc'] and 'MaxCycle=100' in calls[1][0]['irc']
    assert 'geom' not in calls[1][0] and 'guess' not in calls[1][0]
    assert rows[-1][1]['data']['status'] == 'error'
    assert rows[-1][1]['data']['irc_accepted_points'] == 0
    assert Path('irc.irc_attempt1.log').exists()
    assert not Path('irc.xyz').exists()


def test_successful_retry_uses_accepted_point_for_endpoint(tmp_path, monkeypatch):
    calls, rows = _run_irc_worker(tmp_path, monkeypatch,
        [('Error termination\n', True), (_irc_point(.2), False)])
    assert len(calls) == 3
    assert calls[-1][0]['label'] == 'irc_prod'
    assert calls[-1][1][0, 2] == .2
    assert rows[0][1]['data']['irc_completed'] is True
    assert rows[-1][1]['data']['status'] == 'normal'


def test_partial_irc_uses_last_accepted_point_not_failed_trial(tmp_path, monkeypatch):
    calls, rows = _run_irc_worker(tmp_path, monkeypatch,
        [(_irc_point(.1) + _irc_point(.2) + 'CURRENT STRUCTURE\ntruncated\n', True)])
    assert len(calls) == 2
    assert calls[-1][1][0, 2] == .2
    assert rows[0][1]['data']['irc_completed'] is False
    assert rows[0][1]['data']['irc_accepted_points'] == 2
    assert 'energy' not in rows[0][1]['data']


def test_ordinary_irc_keeps_existing_endpoint_recovery(tmp_path, monkeypatch):
    calls, rows = _run_irc_worker(tmp_path, monkeypatch,
        [('Error termination\n', True)] * 2, bimol=False)
    assert len(calls) == 3
    assert calls[1][0]['irc'] == calls[0][0]['irc']
    assert calls[-1][0]['label'] == 'irc_prod'
    assert rows[0][1]['data'] == {'status': 'normal'}
    assert not Path('irc.irc_attempt1.log').exists()
    assert not Path('irc.xyz').exists()


@pytest.mark.parametrize('bimol', [False, True])
def test_saved_irc_accepted_path_check_is_bimolecular_only(tmp_path, monkeypatch, bimol):
    from types import SimpleNamespace
    from kinbot.irc import IRC
    monkeypatch.chdir(tmp_path)
    qc = SimpleNamespace(qc='gauss', use_sella=False, check_qc=lambda name: 'normal')
    qc.for_level = lambda level: qc
    rxn = SimpleNamespace(instance_name='ts', qc=qc)
    Path('ts_IRC_F.log').write_text(_irc_point(.1))
    Path('ts_IRC_R.log').write_text('Error termination after corrector trials\n')
    assert IRC(rxn, dict(bimol=bimol)).check_irc() == ['normal', 'error' if bimol else 'normal']
