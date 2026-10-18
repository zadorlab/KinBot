"""Controller handoff preserves owned native work and rejects unsafe pauses."""
from copy import deepcopy
import ctypes
import errno
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
from types import SimpleNamespace

import pytest

from kinbot import run_cleanup, run_refresh, run_salvage
from kinbot.anl import dispatch
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task


@pytest.fixture
def live_tree(tmp_path, monkeypatch):
    root = tmp_path / 'run'; root.mkdir()
    pin = root / '.salvage_code.new'; pin.mkdir()
    inputs, processes, jobs, graphs = [], [], [], []
    for index, case in enumerate(('ch3s', 'ch3ss')):
        directory = root / case; directory.mkdir()
        path = directory / 'input.json'; path.write_text('{}')
        inputs.append(str(path))
        worker = directory / '123'; worker.mkdir()
        (worker / '.kinbot_run.json').write_text('{}')
        (worker / 'kinbot.db').write_bytes(b'accepted L1/L2 results')
        for folder, module, name in ((directory, 'pes', 'pes_status.json'),
                                      (worker, 'kb', 'kinbot_status.json')):
            pid = 100 + len(processes)
            processes.append(dict(pid=pid, ppid=1, starttime=str(pid * 10),
                                  cwd=str(folder), state='S',
                                  command=f'/python -m kinbot.{module} input.json'))
            atomic_json(folder / name, dict(pid=pid, host=socket.gethostname(), status='running'))
        atomic_json(directory / 'launch.json', dict(pid=processes[-2]['pid']))
        (directory / 'driver.log').write_text('original driver log')
        for branch in ('anl_automatic/states/test/composite', 'vrctst/molpro_scan/scan'):
            graph = worker / branch
            graph.parent.mkdir(parents=True)
            spec = dispatch_spec()
            spec['tasks'] = [molpro_task('energy', 'rhf\n', geometry_from='initial')]
            dispatch.prepare(spec, graph)
            state = dispatch._load(graph)[2]
            job = str(500 + len(jobs))
            state['tasks']['energy'].update(status='submitted', job_id=job)
            atomic_json(graph / 'state.json', state)
            jobs.append(dict(job_id=job, directory=str(graph / 'tasks/energy'), state='RUNNING'))
            graphs.append(graph)
    atomic_json(root / 'prepared.json', dict(inputs=inputs, sha256={
        str(Path(p).relative_to(root)): file_sha256(p) for p in inputs}))
    atomic_json(root / 'salvage.json', dict(status='launched', inputs=inputs,
                                         code=str(root / '.salvage_code.old')))
    monkeypatch.setattr(run_salvage, '_pin_code', lambda *a, **kw: pin)
    monkeypatch.setattr(run_refresh, '_preflight', lambda root: None)
    monkeypatch.setattr(run_salvage, '_preflight_worker_routes', lambda root: None)
    monkeypatch.setattr('kinbot.anl.automatic.validate_automatic_environment', lambda par: None)
    monkeypatch.setattr('kinbot.parameters.Parameters', lambda path: SimpleNamespace(par={}))
    monkeypatch.setattr(run_cleanup, '_jobs', lambda root: deepcopy(jobs))
    monkeypatch.setattr(run_cleanup, '_processes', lambda folder: [deepcopy(p) for p in processes
                        if Path(p['cwd']).is_relative_to(folder)])
    signals = []
    def send(fd, sig):
        if fd == os.getpid() and sig == 0:
            return
        record = next((p for p in processes if p['pid'] == fd), None)
        if record is None:
            raise ProcessLookupError
        if sig == 0:
            return
        signals.append((fd, sig))
        if sig == signal.SIGKILL:
            processes.remove(record)
        else:
            record['state'] = 'T' if sig == signal.SIGSTOP else 'S'
    monkeypatch.setattr(run_refresh.os, 'pidfd_open', lambda pid: pid, raising=False)
    monkeypatch.setattr(run_refresh.signal, 'pidfd_send_signal', send, raising=False)
    monkeypatch.setattr(run_refresh.os, 'close', lambda fd: None)
    monkeypatch.setattr(run_refresh, '_starttime', lambda pid: str(pid * 10))
    monkeypatch.setattr(run_refresh, '_stopped', lambda pid: all(p['state'] == 'T' for p in processes if p['pid'] == pid))
    launches = []
    def launch(command, **kwargs):
        launches.append((command, kwargs))
        return SimpleNamespace(pid=900 + len(launches))
    monkeypatch.setattr(run_salvage.subprocess, 'Popen', launch)
    return SimpleNamespace(root=root, processes=processes, jobs=jobs, graphs=graphs,
                           signals=signals, launches=launches, send=send)


def scientific_files(root):
    return {str(p.relative_to(root)): file_sha256(p) for p in root.rglob('*') if p.is_file()
            and ('composite' in p.parts or 'molpro_scan' in p.parts or p.name == 'kinbot.db')}


def test_refresh_keeps_live_anl_and_vrc_and_restarts_once(live_tree):
    t = live_tree
    before = scientific_files(t.root)
    original_jobs = deepcopy(t.jobs)
    review = run_refresh.refresh(t.root)
    assert review['phase'] == 'reviewed' and t.signals == []
    result = run_refresh.refresh(t.root, apply=True)
    assert result['phase'] == 'launched'
    assert t.jobs == original_jobs and scientific_files(t.root) == before
    assert len(t.launches) == 2
    assert {pid for pid, sig in t.signals if sig == signal.SIGKILL} == {100, 101, 102, 103}
    assert all(args['env']['PYTHONPATH'].startswith(result['code']) for _, args in t.launches)
    archive = Path(result['archive'])
    assert (archive / 'ch3ss/driver.log').read_text() == 'original driver log'
    assert run_refresh.refresh(t.root, apply=True) == result
    assert len(t.launches) == 2


@pytest.mark.parametrize('problem', ['child', 'unowned_job', 'submitting', 'changed_input',
                                     'scheduler', 'staging', 'partial_graph'])
def test_refresh_refuses_uncertain_ownership_without_stopping(live_tree, monkeypatch, problem):
    t = live_tree
    if problem == 'child':
        t.processes.append(dict(pid=105, cwd=str(t.root / 'ch3s/123'), command='/bin/sbatch job.slurm'))
    elif problem == 'unowned_job':
        t.jobs.append(dict(job_id='999', directory=str(t.root / 'ch3s/unowned')))
    elif problem == 'scheduler':
        monkeypatch.setattr(run_cleanup, '_jobs', lambda root: (_ for _ in ()).throw(RuntimeError('outage')))
    elif problem == 'submitting':
        path = t.graphs[0] / 'state.json'
        state = json.loads(path.read_text()); state['tasks']['energy']['status'] = 'submitting'
        atomic_json(path, state)
    elif problem == 'staging':
        (t.graphs[0] / 'tasks/uncommitted').mkdir()
    elif problem == 'partial_graph':
        (t.root / 'ch3s/123/vrctst/molpro_scan/new').mkdir()
    else:
        (t.graphs[0] / 'tasks/energy/energy.inp').write_text('wrong')
    with pytest.raises((RuntimeError, ValueError)):
        run_refresh.refresh(t.root, apply=True)
    assert t.signals == [] and not t.launches


def test_submission_race_after_stop_resumes_original_controllers(live_tree, monkeypatch):
    t = live_tree
    original = run_refresh._owned_jobs
    count = 0
    def inventory(root):
        nonlocal count
        count += 1
        if count == 2:
            raise RuntimeError('submission in flight')
        return original(root)
    monkeypatch.setattr(run_refresh, '_owned_jobs', inventory)
    with pytest.raises(RuntimeError, match='submission in flight'):
        run_refresh.refresh(t.root, apply=True)
    assert all(p['state'] == 'S' for p in t.processes)
    assert not any(sig == signal.SIGKILL for _, sig in t.signals)
    assert not t.launches


def test_pid_reuse_does_not_receive_any_signal(live_tree, monkeypatch):
    monkeypatch.setattr(run_refresh, '_starttime', lambda pid: 'reused')
    with pytest.raises(RuntimeError, match='identity changed'):
        run_refresh.refresh(live_tree.root, apply=True)
    assert not live_tree.signals


def test_interrupted_archive_resumes_without_touching_qc(live_tree, monkeypatch):
    t = live_tree
    before = scientific_files(t.root)
    original = run_salvage._save
    count = 0
    def archive(*args, **kwargs):
        nonlocal count
        count += 1
        if count == 3:
            raise RuntimeError('interrupted archiving')
        return original(*args, **kwargs)
    monkeypatch.setattr(run_salvage, '_save', archive)
    with pytest.raises(RuntimeError, match='interrupted archiving'):
        run_refresh.refresh(t.root, apply=True)
    assert not t.processes and not t.launches
    monkeypatch.setattr(run_salvage, '_save', original)
    assert run_refresh.refresh(t.root, apply=True)['phase'] == 'launched'
    assert scientific_files(t.root) == before
    assert len(t.launches) == 2


@pytest.mark.parametrize('missing', [('open',), ('send',), ('open', 'send')])
def test_missing_python_pidfd_apis_use_linux_syscalls_for_full_handoff(live_tree, monkeypatch, missing):
    t = live_tree
    before = scientific_files(t.root)
    if 'open' in missing:
        monkeypatch.delattr(run_refresh.os, 'pidfd_open')
    if 'send' in missing:
        monkeypatch.delattr(run_refresh.signal, 'pidfd_send_signal')
    monkeypatch.setattr(run_refresh.sys, 'platform', 'linux')
    monkeypatch.setattr(run_refresh.platform, 'machine', lambda: 'x86_64')
    calls = []
    def syscall(number, *args):
        values = [arg.value for arg in args]
        calls.append((number.value, values))
        if number.value == 434:
            assert values[1] == 0
            return values[0]
        assert number.value == 424 and values[2:] == [None, 0]
        try:
            t.send(values[0], values[1])
        except ProcessLookupError:
            ctypes.set_errno(errno.ESRCH)
            return -1
        return 0
    monkeypatch.setattr(run_refresh.ctypes, 'CDLL', lambda *a, **kw: SimpleNamespace(syscall=syscall))
    assert run_refresh.refresh(t.root, apply=True)['phase'] == 'launched'
    assert calls and len(t.launches) == 2
    assert scientific_files(t.root) == before


def test_kernel_probe_failure_happens_before_journal_or_stop(live_tree, monkeypatch):
    t = live_tree
    def unavailable(pid):
        raise OSError(errno.ENOSYS, 'kernel lacks pidfd_open')
    monkeypatch.setattr(run_refresh, '_pidfd_open', unavailable)
    with pytest.raises(RuntimeError, match='No controllers were paused'):
        run_refresh.refresh(t.root, apply=True)
    assert not t.signals and not t.launches
    assert not (t.root / 'controller_refresh.json').exists()


def test_all_handles_checked_before_first_stop(live_tree, monkeypatch):
    t = live_tree
    def unavailable(pid):
        if pid == 102:
            raise PermissionError('second case is inaccessible')
        return pid
    monkeypatch.setattr(run_refresh, '_pidfd_open', unavailable)
    with pytest.raises(PermissionError, match='inaccessible'):
        run_refresh.refresh(t.root, apply=True)
    assert not t.signals and not t.launches


def test_failed_old_preparation_accepts_fixed_code_without_manual_cleanup(live_tree):
    t = live_tree
    before = scientific_files(t.root)
    archive = t.root.parent / 'failed_refresh'; archive.mkdir()
    atomic_json(t.root / 'controller_refresh.json', {
        'phase': 'preparing', 'code': str(t.root / '.salvage_code.bad'),
        'previous_code': str(t.root / '.salvage_code.old'), 'archive': str(archive),
    })
    result = run_refresh.refresh(t.root, apply=True)
    assert result['phase'] == 'launched' and len(t.launches) == 2
    assert result['superseded_preparation_codes'] == [str(t.root / '.salvage_code.bad')]
    assert scientific_files(t.root) == before


@pytest.mark.parametrize('phase', ['stopping', 'stopped', 'prepared'])
def test_code_replacement_rejects_an_already_committed_handoff(live_tree, phase):
    t = live_tree
    archive = t.root.parent / 'failed_refresh'; archive.mkdir()
    atomic_json(t.root / 'controller_refresh.json', {
        'phase': phase, 'code': str(t.root / '.salvage_code.bad'),
        'previous_code': str(t.root / '.salvage_code.old'), 'archive': str(archive),
    })
    with pytest.raises(RuntimeError, match='already changed'):
        run_refresh.refresh(t.root, apply=True)
    assert not t.signals and not t.launches


@pytest.mark.parametrize('host,arch', [('darwin', 'arm64'), ('linux', 'unknown')])
def test_syscall_fallback_refuses_other_platforms_before_loading_libc(monkeypatch, host, arch):
    monkeypatch.setattr(run_refresh.sys, 'platform', host)
    monkeypatch.setattr(run_refresh.platform, 'machine', lambda: arch)
    def forbidden(*args, **kwargs):
        pytest.fail('must not issue a Linux syscall on an unrecognized ABI')
    monkeypatch.setattr(run_refresh.ctypes, 'CDLL', forbidden)
    with pytest.raises(RuntimeError, match='compatibility requires'):
        run_refresh._pidfd_syscall(434, ctypes.c_int(1), ctypes.c_uint(0))


@pytest.mark.skipif(sys.platform != 'linux', reason='requires the real Linux pidfd ABI')
def test_real_linux_process_handles_without_python_wrappers(monkeypatch):
    """Linux CI exercises the native ABI; no controller or scientific job is used."""
    monkeypatch.delattr(os, 'pidfd_open', raising=False)
    monkeypatch.delattr(signal, 'pidfd_send_signal', raising=False)
    run_refresh._probe_process_handles()
    child = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(30)'])
    fd = None
    try:
        fd = run_refresh._pidfd_open(child.pid)
        run_refresh._pidfd_send_signal(fd, signal.SIGSTOP)
        _, status = os.waitpid(child.pid, os.WUNTRACED)
        assert os.WIFSTOPPED(status)
        run_refresh._pidfd_send_signal(fd, signal.SIGCONT)
        run_refresh._pidfd_send_signal(fd, signal.SIGKILL)
        assert child.wait(timeout=5) == -signal.SIGKILL
        with pytest.raises(ProcessLookupError):
            run_refresh._pidfd_send_signal(fd, 0)
    finally:
        if child.poll() is None:
            child.kill(); child.wait(timeout=5)
        if fd is not None:
            os.close(fd)
