"""Interrupted bimolecular workers must finish as failures without new jobs."""
import json
import pickle
from pathlib import Path
from types import SimpleNamespace as NS
from unittest.mock import Mock

import pytest
from ase import Atoms
from ase.db import connect

from kinbot import bimol
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry


SEARCH = bimol.search_job('ethyl_pair', sequential=True)
WORKERS = [SEARCH, bimol.search_job('ethyl_pair'),
           bimol.validation_job('bimol/unique_candidates/candidate_1'),
           'product_spin/ethylidene_m1']


def stopped_worker(tmp_path, monkeypatch, job=SEARCH, queue='slurm', bimolecular=True):
    monkeypatch.chdir(tmp_path)
    qc = QuantumChemistry.__new__(QuantumChemistry)
    qc.qc, qc.queuing, qc.theory_level = 'gauss', queue, None
    qc.par = {'bimol': bimolecular}
    qc.db = connect('kinbot.db')
    qc.job_ids = {job: '123'}
    qc.qc_bimol_search = Mock()
    process = Mock()
    process.communicate.return_value = (b'JOBID NAME\n', b'')
    monkeypatch.setattr('kinbot.qc.subprocess.Popen', Mock(return_value=process))
    monkeypatch.setattr('kinbot.qc.subprocess.call', Mock(return_value=1))
    clock = [0.]
    monkeypatch.setattr('kinbot.qc.time.monotonic', lambda: clock[0])
    return qc, clock


@pytest.mark.parametrize('queue', ['slurm', 'pbs'])
def test_interrupted_pose_search_does_not_restart(tmp_path, monkeypatch, caplog, queue):
    qc, clock = stopped_worker(tmp_path, monkeypatch, queue=queue)
    partial = Path(SEARCH).parent / 'search.json'
    partial.parent.mkdir(parents=True)
    partial.write_text('{"observations": [], "attempts": [{"pose": 0}]}\n')
    saved = partial.read_bytes()
    state = dict(reaction=NS(instance_name='ethyl_pair'), done=False,
                 active=None, attempt=0, cursor=None)
    species = NS(bimol_pose_searches=[state])
    assert bimol.advance_pose_searches(species, {}, qc)
    clock[0] = 59.
    assert bimol.advance_pose_searches(species, {}, qc)
    clock[0] = 61.
    assert not bimol.advance_pose_searches(species, {}, qc)
    assert state['done']
    clock[0] = 3600.
    assert not bimol.advance_pose_searches(species, {}, qc)
    qc.qc_bimol_search.assert_not_called()
    assert partial.read_bytes() == saved
    assert 'this run will not submit the job again' in caplog.text


@pytest.mark.parametrize('job', WORKERS)
def test_missing_bimol_result_stays_failed_after_grace(tmp_path, monkeypatch, job):
    qc, clock = stopped_worker(tmp_path, monkeypatch, job)
    assert qc.check_qc(job) == 'running'
    clock[0] = 61.
    assert qc.check_qc(job) == 'error'
    clock[0] = 3600.
    assert qc.check_qc(job) == 'error'
    assert list(qc.db.select()) == []


@pytest.mark.parametrize('arrival', [59., 61.])
def test_complete_result_wins_over_missing_result_deadline(tmp_path, monkeypatch, arrival):
    qc, clock = stopped_worker(tmp_path, monkeypatch)
    assert qc.check_qc(SEARCH) == 'running'
    Path(SEARCH).parent.mkdir(parents=True)
    atoms = Atoms('C2H5')
    with open(SEARCH + '.pkl', 'wb') as handle:
        pickle.dump(dict(sym=atoms.symbols, pos=atoms.positions,
                         name=SEARCH, data={'status': 'normal'}), handle)
    Path(SEARCH + '.done').write_text('done\n')
    clock[0] = arrival
    assert qc.check_qc(SEARCH) == 'normal'
    assert qc.db.get(name=SEARCH).data['status'] == 'normal'


@pytest.mark.parametrize('job,bimolecular', [(SEARCH, False),
    ('ethyl_pair_l0_ts_1', True), ('ethyl_well', False), ('hir/ethyl_hir_0_0', True)])
def test_existing_non_bimol_worker_status_is_unchanged(tmp_path, monkeypatch, job, bimolecular):
    qc, clock = stopped_worker(tmp_path, monkeypatch, job, bimolecular=bimolecular)
    assert qc.check_qc(job) == 'running'
    clock[0] = 61.
    assert qc.check_qc(job) == 0


def test_unsubmitted_search_can_still_start(tmp_path, monkeypatch):
    qc, clock = stopped_worker(tmp_path, monkeypatch)
    qc.job_ids = {}
    state = dict(reaction=NS(instance_name='ethyl_pair'), done=False,
                 active=None, attempt=0, cursor=None)
    assert bimol.advance_pose_searches(NS(bimol_pose_searches=[state]), {}, qc)
    qc.qc_bimol_search.assert_called_once_with(state['reaction'], attempt=0, start_after=None)


@pytest.mark.parametrize('job,level', [(SEARCH, 'l0'), (WORKERS[2], 'l0'),
                                     (WORKERS[3], 'l1')])
def test_restart_recognizes_the_interrupted_submission(tmp_path, monkeypatch, job, level):
    monkeypatch.chdir(tmp_path)
    Path('input.json').write_text(json.dumps(dict(
        barrier_threshold=50, profiled_theory=1, queuing='slurm')))
    par = Parameters('input.json').par
    par['bimol'] = 1
    original = QuantumChemistry(par)
    original._record(job, level, '123')
    restarted = QuantumChemistry(par)
    assert restarted.for_level(level).job_ids[job] == '123'
    # The queue and result files are absent after the interrupted submission.
    monkeypatch.setattr(QuantumChemistry, '_check_qc', lambda *args: 0)
    clock = [0.]
    monkeypatch.setattr('kinbot.qc.time.monotonic', lambda: clock[0])
    assert restarted.check_qc(job) == 'running'
    clock[0] = 61.
    assert restarted.check_qc(job) == 'error'
