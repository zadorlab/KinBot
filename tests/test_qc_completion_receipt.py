"""Interrupted admission must preserve completed generations across a refresh."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from ase import Atoms

from kinbot.execution import atomic_json
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.slurm_retry import completion_receipt


def receipt(root, name='path', job='42', **updates):
    directory = root/'.kinbot_slurm_io_retries'/hashlib.sha256(name.encode()).hexdigest()[:20]/job
    directory.mkdir(parents=True, exist_ok=True)
    (directory/'lock').touch()
    record = dict(schema=1, identity=name, slurm_job_id=job, slurm_restart_count=0,
                  status='complete', attempts=[dict(returncode=0, failed_result=False)])
    record.update(updates)
    atomic_json(directory/'execution.json', record)
    return directory


def test_admission_interruption_preserves_marker_and_result(tmp_path, monkeypatch):
    from kinbot import run_policy
    monkeypatch.chdir(tmp_path)
    Path('path.exitcode').write_text('42 0\n')
    Path('path.pkl').write_bytes(b'completed handoff')
    qc = SimpleNamespace(check_qc=lambda job: 'normal', queue_job_limit=1,
        par={'queue_template': '', 'queue_walltime': '01:00:00'}, queuing='slurm',
        qc='gauss', queue_name='cpu', slurm_feature='', job_ids={'path':'42'})
    monkeypatch.setattr(run_policy, 'find', lambda: tmp_path/'policy.json')
    def interrupted(*args, **kwargs):
        assert Path('path.exitcode').read_text() == '42 0\n'
        assert Path('path.pkl').read_bytes() == b'completed handoff'
        raise KeyboardInterrupt
    monkeypatch.setattr(run_policy, 'submit', interrupted)
    monkeypatch.setattr('kinbot.qc._ordinary_slurm_walltime', lambda *a: '#SBATCH --time=01:00:00')
    with pytest.raises(KeyboardInterrupt):
        QuantumChemistry._submit_qc(qc, 'path', 8, singlejob=0)
    assert qc.job_ids == {'path':'42'}
    assert Path('path.exitcode').read_text() == '42 0\n'
    assert Path('path.pkl').read_bytes() == b'completed handoff'
    assert '"$kinbot_exit_tmp" path.pkl' in Path('path.sbatch').read_text()


@pytest.mark.parametrize('active,expected', [(False, 'normal'), (True, 'running'), (None, 'running')])
def test_receipt_reuses_completed_database_only_after_scheduler_release(tmp_path, monkeypatch, active, expected):
    monkeypatch.chdir(tmp_path)
    atomic_json(tmp_path/'input.json', {'queuing':'slurm','barrier_threshold':100.})
    qc = QuantumChemistry(Parameters('input.json').par)
    qc.qc = 'fc'
    qc.job_ids['path'] = '42'
    qc.db.write(Atoms('H'), name='path', data={'status':'normal','energy':-1.})
    receipt(tmp_path)
    monkeypatch.setattr('kinbot.qc._slurm_job_active', lambda job: active)
    assert qc._check_qc('path') == expected
    assert not Path('path.exitcode').exists()  # Read-only adoption, not a forged marker.
    if active is False:
        Path('path.exitcode').write_text('43 0\n')
        assert qc._check_qc('path') == 0  # Never use an old receipt over a newer marker.


@pytest.mark.parametrize('updates', [
    {'identity':'different'}, {'slurm_job_id':'43'}, {'status':'running'},
    {'attempts':[]}, {'attempts':[{'returncode':1,'failed_result':False}]},
    {'slurm_restart_count':1}, {'attempts':[{'returncode':False,'failed_result':False}]},
])
def test_invalid_or_incomplete_receipt_does_not_prove_success(tmp_path, updates):
    receipt(tmp_path, **updates)
    assert completion_receipt('path', '42', tmp_path) is None


def test_latest_requeue_and_wrapper_lock_take_precedence(tmp_path):
    import fcntl
    folder = receipt(tmp_path)
    assert completion_receipt('path','42',tmp_path) == (0, True)
    with (folder/'lock').open('r+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert completion_receipt('path','42',tmp_path) is None
    (folder/'restart-1').mkdir()
    assert completion_receipt('path','42',tmp_path) is None
    record = json.loads((folder/'execution.json').read_text())
    record.update(slurm_restart_count=1, status='failed',
                  attempts=[dict(returncode=9,failed_result=True)])
    atomic_json(folder/'restart-1/execution.json', record)
    assert completion_receipt('path','42',tmp_path) == (9, True)
