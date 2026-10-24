"""Cleanup may affect only the named calculation, including on PID reuse."""
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import run_cleanup as cleanup


def test_job_inventory_scopes_directory_boundaries_and_owner(tmp_path, monkeypatch):
    root = tmp_path / 'run'
    root.mkdir()
    monkeypatch.setattr(cleanup.subprocess, 'run', lambda command, **kw: SimpleNamespace(
        stdout=f'11|RUNNING|node1|{root}\n12|COMPLETING|node2|{root}/reference/composite\n'
               f'13|RUNNING|node3|{root}_other\n14|RUNNING|node4|/unrelated\n'))
    jobs = cleanup._jobs(root)
    assert [job['job_id'] for job in jobs] == ['11', '12']
    assert jobs[1]['state'] == 'COMPLETING'


def test_stop_cancels_only_scoped_jobs_and_rechecks_pid_identity(tmp_path, monkeypatch):
    root = tmp_path / 'run'
    root.mkdir()
    record = {'pid': 500, 'starttime': 'original', 'directory': str(root)}
    inventory = {'directory': str(root), 'jobs': [{'job_id': '11'}], 'processes': [record]}
    monkeypatch.setattr(cleanup, 'inspect', lambda path, **kw: inventory if not (root / 'cleanup_inventory.json').exists()
                        else {'directory': str(root), 'jobs': [], 'processes': []})
    monkeypatch.setattr(cleanup, '_processes', lambda path, **kw: [{'pid': 500, 'starttime': 'reused'}])
    monkeypatch.setattr(cleanup.time, 'sleep', lambda delay: None)
    monkeypatch.setattr(cleanup.os, 'kill', lambda *args: pytest.fail('Reused PID must not be killed'))
    commands = []
    monkeypatch.setattr(cleanup.subprocess, 'run', lambda command, **kw: commands.append(command))
    result = cleanup.stop(root, archive=True)
    assert commands == [['scancel', '11']]
    assert result['status'] == 'stopped'
    assert not root.exists()
    assert (Path(result['archive']) / 'cleanup_inventory.json').is_file()
    assert json.loads((Path(result['archive']) / 'cleanup_status.json').read_text())['status'] == 'stopped'


def test_scheduler_error_cannot_authorize_cleanup(tmp_path, monkeypatch):
    root = tmp_path / 'run'
    root.mkdir()
    (root / '.kinbot_run.json').write_text('{}')
    def unavailable(*args):
        raise RuntimeError('scheduler unavailable')
    monkeypatch.setattr(cleanup, '_jobs', unavailable)
    monkeypatch.setattr(cleanup, '_processes', lambda *args: pytest.fail('Scheduler check must happen first'))
    with pytest.raises(RuntimeError, match='scheduler unavailable'):
        cleanup.stop(root, archive=True)
    assert root.is_dir()


def test_proc_inventory_preserves_invoking_shell_and_tracks_children_in_scratch(tmp_path):
    proc = tmp_path / 'proc'
    proc.mkdir()
    root = tmp_path / 'run'
    root.mkdir()
    scratch = tmp_path / 'scratch'
    scratch.mkdir()
    current = os.getpid()
    entries = [(current, 100, root), (100, 1, root), (200, 100, root), (201, 200, scratch),
               (300, 100, tmp_path / 'run_other')]
    for pid, parent, cwd in entries:
        path = proc / str(pid)
        path.mkdir()
        (path / 'cwd').symlink_to(cwd)
        # Fields after comm: state, ppid, ..., starttime (index 19).
        (path / 'stat').write_text(f'{pid} (worker (test)) S {parent} ' + '0 ' * 17 + f'{pid} 0')
        (path / 'cmdline').write_bytes(b'python\0worker.py')
    assert [item['pid'] for item in cleanup._processes(root, proc)] == [200, 201]
    known = cleanup._processes(root, proc)
    # The scratch worker outlives its parent. Tracking must survive reparenting,
    # but must not carry over to a different process reusing the same PID.
    (proc / '200/stat').write_text('200 (done) Z 100 ' + '0 ' * 17 + '200 0')
    (proc / '201/stat').write_text('201 (worker) S 1 ' + '0 ' * 17 + '201 0')
    assert cleanup._processes(root, proc) == []
    assert [item['pid'] for item in cleanup._processes(root, proc, known=known)] == [201]
    (proc / '201/stat').write_text('201 (other) S 1 ' + '0 ' * 17 + '999 0')
    assert cleanup._processes(root, proc, known=known) == []


def _fake_clock(monkeypatch):
    clock = [0.]
    monkeypatch.setattr(cleanup.time, 'monotonic', lambda: clock[0])
    def advance(delay):
        clock[0] += delay
    monkeypatch.setattr(cleanup.time, 'sleep', advance)
    return clock


def test_cancellation_waits_for_teardown_and_catches_late_submission(tmp_path, monkeypatch):
    clock = _fake_clock(monkeypatch)
    root = tmp_path / 'run'
    root.mkdir()
    def inventory(path, **kw):
        jobs = [{'job_id': '11', 'state': 'COMPLETING'}] if clock[0] < 42 else []
        if 2 <= clock[0] < 44:
            jobs.append({'job_id': '12', 'state': 'PENDING'})
        return {'directory': str(root), 'jobs': jobs, 'processes': []}
    monkeypatch.setattr(cleanup, 'inspect', inventory)
    commands = []
    monkeypatch.setattr(cleanup.subprocess, 'run', lambda command, **kw: commands.append(command))
    assert cleanup.stop(root)['status'] == 'stopped'
    assert clock[0] == 44
    assert commands == [['scancel', '11'], ['scancel', '12']]


def test_timeout_preserves_tree_and_reports_remaining_state(tmp_path, monkeypatch, capsys):
    _fake_clock(monkeypatch)
    root = tmp_path / 'run'
    root.mkdir()
    inventory = {'directory': str(root), 'jobs': [{'job_id': '11', 'state': 'COMPLETING'}],
                 'processes': []}
    monkeypatch.setattr(cleanup, 'inspect', lambda path, **kw: inventory)
    monkeypatch.setattr(cleanup.subprocess, 'run', lambda *args, **kw: None)
    assert cleanup.main([str(root), '--stop', '--archive', '--timeout', '4']) == 1
    output = json.loads(capsys.readouterr().out)
    assert output['status'] == 'stop_pending'
    assert output['jobs'][0]['state'] == 'COMPLETING'
    assert output['elapsed_seconds'] == 4
    assert root.is_dir()
    assert json.loads((root / 'cleanup_status.json').read_text()) == output


@pytest.mark.parametrize('still_queued', [False, True])
def test_scancel_exit_race_requires_verified_absence(tmp_path, monkeypatch, still_queued):
    def cancel(command, **kw):
        raise cleanup.subprocess.CalledProcessError(1, command, stderr='Invalid job id')
    monkeypatch.setattr(cleanup.subprocess, 'run', cancel)
    monkeypatch.setattr(cleanup, '_jobs', lambda path: [{'job_id': '11'}] if still_queued else [])
    if still_queued:
        with pytest.raises(RuntimeError, match='Slurm cancellation failed'):
            cleanup._cancel(tmp_path, [{'job_id': '11'}])
    else:
        cleanup._cancel(tmp_path, [{'job_id': '11'}])


@pytest.mark.parametrize('timeout', [0, -1, float('nan'), float('inf')])
def test_invalid_timeout_cannot_start_cleanup(monkeypatch, tmp_path, timeout):
    monkeypatch.setattr(cleanup, 'inspect', lambda *args: pytest.fail('No side effects allowed'))
    with pytest.raises(ValueError, match='finite positive'):
        cleanup.stop(tmp_path, timeout=timeout)


def test_cancel_allocations_rechecks_completion_and_late_job_without_process_signals(tmp_path,monkeypatch):
    rounds=iter([[{'job_id':'11'}],[{'job_id':'11'},{'job_id':'12'}],[{'job_id':'12'}],[]])
    monkeypatch.setattr(cleanup,'_jobs',lambda root:next(rounds))
    calls=[]
    monkeypatch.setattr(cleanup,'_cancel',lambda root,jobs:calls.append([j['job_id'] for j in jobs]))
    monkeypatch.setattr(cleanup.time,'sleep',lambda _:None)
    monkeypatch.setattr(cleanup.os,'kill',lambda *a:pytest.fail('No local process signals'))
    assert cleanup.cancel_allocations(tmp_path)['cancelled_job_ids']==['11','12']
    assert calls==[['11'],['12']]


def test_cancel_allocations_scheduler_failure_is_not_empty_queue(tmp_path,monkeypatch):
    def fail(root):raise RuntimeError('Slurm unavailable')
    monkeypatch.setattr(cleanup,'_jobs',fail)
    with pytest.raises(RuntimeError,match='Slurm unavailable'):
        cleanup.cancel_allocations(tmp_path)
