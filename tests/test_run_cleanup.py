"""Cleanup may affect only the named calculation, including on PID reuse."""
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import run_cleanup as cleanup


def test_job_inventory_scopes_directory_boundaries_and_owner(tmp_path, monkeypatch):
    root = tmp_path / 'run'
    root.mkdir()
    monkeypatch.setattr(cleanup.subprocess, 'run', lambda command, **kw: SimpleNamespace(
        stdout=f'11|{root}\n12|{root}/reference/composite\n13|{root}_other\n14|/unrelated\n'))
    assert [job['job_id'] for job in cleanup._jobs(root)] == ['11', '12']


def test_stop_cancels_only_scoped_jobs_and_rechecks_pid_identity(tmp_path, monkeypatch):
    root = tmp_path / 'run'
    root.mkdir()
    record = {'pid': 500, 'starttime': 'original', 'directory': str(root)}
    inventory = {'directory': str(root), 'jobs': [{'job_id': '11'}], 'processes': [record]}
    monkeypatch.setattr(cleanup, 'inspect', lambda path: inventory if not (root / 'cleanup_inventory.json').exists()
                        else {'directory': str(root), 'jobs': [], 'processes': []})
    monkeypatch.setattr(cleanup, '_processes', lambda path: [{'pid': 500, 'starttime': 'reused'}])
    monkeypatch.setattr(cleanup.os, 'kill', lambda *args: pytest.fail('Reused PID must not be killed'))
    commands = []
    monkeypatch.setattr(cleanup.subprocess, 'run', lambda command, **kw: commands.append(command))
    result = cleanup.stop(root, archive=True)
    assert commands == [['scancel', '11']]
    assert result['status'] == 'stopped'
    assert not root.exists()
    assert (Path(result['archive']) / 'cleanup_inventory.json').is_file()


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
