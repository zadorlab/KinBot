"""Recorded levels take precedence over ordinary calculation name fallbacks."""
import unittest
from types import SimpleNamespace
from kinbot.profiled_qc import ProfiledQuantumChemistry

import json
from pathlib import Path
import pytest
from ase import Atoms
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry


@pytest.mark.parametrize('bimol, job, expected', [
    (False, 'methanol_well', 'error'),
    (True, 'methanol_well', 'error'),
    (True, 'bimol/methanol_search', 'normal'),
])
def test_only_bimolecular_searches_accept_separate_done_files(
        tmp_path, monkeypatch, bimol, job, expected):
    monkeypatch.chdir(tmp_path)
    Path('input.json').write_text(json.dumps(dict(
        barrier_threshold=100., queuing='local', error_missing_local=0)))
    qc = QuantumChemistry(Parameters('input.json', show_warnings=False).par)
    qc.par['bimol'] = bimol
    Path(job).parent.mkdir(exist_ok=True)
    qc.db.write(Atoms('H', positions=[[0., 0., 0.]]), name=job,
                data=dict(status='normal', energy=-1.))
    Path(job + '.done').write_text('done\n')
    assert qc.check_qc(job) == expected
    Path(job + '.log').write_text('done\n')
    assert qc.check_qc(job) == 'normal'


class TestBimolecularBackendRouting(unittest.TestCase):
    def test_only_ring_conformer_suffixes_select_l0(self):
        qc = ProfiledQuantumChemistry.__new__(ProfiledQuantumChemistry)
        qc.par = {'high_level': 0}
        qc._routes = {'ts_bimol_0001': {'level': 'l0'}}
        qc.backends = {'l0': SimpleNamespace(qc='gauss'),
                       'l1': SimpleNamespace(qc='qchem')}
        for job in ('conf/molecule_r0000_0001', 'conf/ts_r12_0345',
                    'conf/molecule_semi_emp_0001', 'ts_bimol_0001'):
            self.assertEqual(qc.for_job(job).qc, 'gauss', job)
        for job in ('conf/ts_0001_freq_recovery_10', 'conf/radical_reaction_0001',
                    'conf/ts_recombination_low', 'ts_hir_restart_1'):
            self.assertEqual(qc.for_job(job).qc, 'qchem', job)
        with self.assertRaisesRegex(ValueError, 'no recorded theory level'):
            qc.for_job('hir/rotor_01')
        qc._routes['conf/molecule_r0000_0001'] = {'level': 'l1'}
        self.assertEqual(qc.for_job('conf/molecule_r0000_0001').qc, 'qchem')


@pytest.mark.parametrize('bimol', [False, True])
@pytest.mark.parametrize('job', ['bimol/methanol_well', 'bimol/methanol_high'])
@pytest.mark.parametrize('recorded', [None, 'l0', 'l1', 'l2'])
def test_bimol_prefix_requires_bimol_and_preserves_saved_level(bimol, job, recorded):
    qc = ProfiledQuantumChemistry.__new__(ProfiledQuantumChemistry)
    qc.par = {'bimol': bimol}
    qc._routes = {job: {'level': recorded}} if recorded else {}
    expected = recorded or ('l0' if bimol else 'l2' if job.endswith('_high') else 'l1')
    assert qc._level_for(job) == expected
