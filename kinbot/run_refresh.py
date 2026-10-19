"""Replace salvage controllers while retaining dispatcher-owned Slurm jobs.

Only idle PES/KB control processes are replaceable. Native calculations,
reference subprocesses and ambiguous submissions cause a refusal, not a kill.
No scancel, task reset, geometry migration or calculation deletion is used.
"""
import argparse
from contextlib import ExitStack, contextmanager
import ctypes
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import platform
import shlex
import signal
import socket
import sys
import time

from kinbot.execution import atomic_json
from kinbot import run_cleanup, run_salvage


def _pidfd_syscall(number, *args):
    """Use the Linux ABI when Python was built without pidfd wrappers.

    Linux arch/x86/entry/syscalls/syscall_64.tbl and asm-generic/unistd.h
    assign 434 to pidfd_open and 424 to pidfd_send_signal on these ABIs.
    Never issue these syscall numbers on macOS or a different word size/ABI.
    """
    if (sys.platform != 'linux' or platform.machine().lower() not in ('x86_64', 'aarch64')
            or ctypes.sizeof(ctypes.c_void_p) != 8 or ctypes.sizeof(ctypes.c_long) != 8):
        raise RuntimeError('PID-handle compatibility requires 64-bit x86_64/aarch64 Linux.')
    library = ctypes.CDLL(None, use_errno=True)
    call = library.syscall
    call.restype = ctypes.c_long
    call.argtypes = [ctypes.c_long]  # syscall's remaining arguments are variadic.
    result = call(ctypes.c_long(number), *args)
    if result == -1:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))
    return result


def _pidfd_open(pid):
    function = getattr(os, 'pidfd_open', None)
    if callable(function):
        return function(pid)
    return _pidfd_syscall(434, ctypes.c_int(pid), ctypes.c_uint(0))


def _pidfd_send_signal(fd, sig):
    function = getattr(signal, 'pidfd_send_signal', None)
    if callable(function):
        return function(fd, sig)
    return _pidfd_syscall(424, ctypes.c_int(fd), ctypes.c_int(sig),
                          ctypes.c_void_p(), ctypes.c_uint(0))


def _probe_process_handles():
    """Check kernel support and both bindings without signalling a controller."""
    try:
        fd = _pidfd_open(os.getpid())
        try:
            _pidfd_send_signal(fd, 0)  # Permission/existence check, no signal.
        finally:
            os.close(fd)
    except (OSError, RuntimeError, AttributeError) as error:
        raise RuntimeError('PID-safe process control is unavailable on this host: '
                           f'{error}. No controllers were paused or replaced.') from error


def _controllers(root):
    """Match owned local processes against recorded controller roles/PIDs."""
    _, inputs = run_salvage._cases(root)
    roles = {}
    for path in inputs:
        roles[path.parent] = ('kinbot.pes', 'pes_status.json')
        for worker in path.parent.glob('[0-9]*'):
            if worker.is_dir():
                roles[worker] = ('kinbot.kb', 'kinbot_status.json')
    for directory, (_, name) in roles.items():
        for marker in (directory / name, directory / 'launch.json'):
            if marker.exists():
                host = json.loads(marker.read_text()).get('host')
                if host and host.split('.')[0] != socket.gethostname().split('.')[0]:
                    raise RuntimeError(f'Run controller refresh on {host}, not this host.')
    records = run_cleanup._processes(root)
    for record in records:
        directory = Path(record['cwd']).resolve()
        role = roles.get(directory)
        args = shlex.split(record['command'])
        module = args[args.index('-m') + 1] if '-m' in args and args.index('-m') + 1 < len(args) else None
        if role is None or module != role[0]:
            raise RuntimeError(f'Local child/process {record["pid"]} is active; '
                               'leave it running and repeat refresh after it exits.')
        marker = directory / role[1]
        saved = json.loads(marker.read_text()) if marker.exists() else {}
        if saved.get('pid') != record['pid']:
            raise RuntimeError(f'{marker}: controller PID does not match; nothing stopped.')
    return records


def _graphs(root):
    for pattern in ('*/[0-9]*/anl_automatic/states/*/composite',
                    '*/[0-9]*/vrctst/molpro_scan/*',
                    '*/[0-9]*/vrctst/molpro_pivots/*',
                    '*/[0-9]*/vrctst/molpro/dispatch*'):
        for graph in sorted(root.glob(pattern)):
            if not graph.is_dir():
                continue
            if not all((graph / name).is_file() for name in ('workflow.json', 'state.json')):
                raise RuntimeError(f'{graph}: graph preparation is incomplete; repeat later.')
            yield graph / 'state.json'


def _owned_jobs(root):
    """Verify both composite and VRC job ownership without reconciling state."""
    from kinbot.anl import dispatch
    jobs = run_cleanup._jobs(root)  # A query failure must abort the handoff.
    owned = {}
    for path in _graphs(root):
        graph, spec, state = dispatch._load(path.parent)
        tasks = {task['id']: task for task in spec['tasks']}
        for directory in (graph / 'tasks').iterdir():
            if directory.is_dir() and directory.name not in state['tasks']:
                raise RuntimeError(f'{directory}: task staging is in flight; repeat later.')
        for ident, entry in state['tasks'].items():
            if entry['status'] == 'submitting':
                raise RuntimeError(f'{path}/{ident}: submission is in flight; repeat later.')
            if entry.get('task_sha256'):
                # Also detect interrupted staging/retry writes, before their
                # new in-memory entry reaches the atomic state ledger.
                dispatch._verify_stage_files(graph, tasks[ident], entry)
            if entry['status'] != 'submitted':
                continue
            if not entry.get('task_sha256'):
                raise RuntimeError(f'{graph}/{ident}: submitted job lacks staged provenance.')
            job = str(entry['job_id'])
            if job in owned:
                raise RuntimeError(f'Slurm job {job} has duplicate graph ownership.')
            owned[job] = (graph / 'tasks' / ident).resolve()
    for job in jobs:
        if owned.get(job['job_id']) != Path(job['directory']).resolve():
            raise RuntimeError(f'Live Slurm job {job["job_id"]} lacks matching '
                               'ANL/VRC dispatcher ownership; nothing reset.')
    return jobs


def _starttime(pid):
    return (Path('/proc') / str(pid) / 'stat').read_text().rsplit(')', 1)[1].split()[19]


def _stopped(pid):
    """SIGSTOP must have stopped every Python thread before reading graphs."""
    threads = list((Path('/proc') / str(pid) / 'task').glob('*/stat'))
    return bool(threads) and all(
        path.read_text().rsplit(')', 1)[1].split()[0] in ('T', 't') for path in threads)


@contextmanager
def _freeze(root, records):
    """Use Linux pidfds so neither STOP, KILL nor CONT can target a reused PID."""
    handles = []
    stopped = []
    try:
        for record in records:
            fd = _pidfd_open(record['pid'])
            handles.append((record, fd))
            if _starttime(record['pid']) != record['starttime']:
                raise RuntimeError('Controller identity changed; retry inventory.')
            _pidfd_send_signal(fd, 0)
        # Open/check every handle before pausing even the first controller.
        for _, fd in handles:
            _pidfd_send_signal(fd, signal.SIGSTOP)
            stopped.append(fd)
        deadline = time.monotonic() + 5
        while not all(_stopped(record['pid']) for record, _ in handles):
            if time.monotonic() >= deadline:
                raise RuntimeError('Controllers did not quiesce; original processes resumed.')
            time.sleep(.05)
        # A controller could have spawned a child between inventory and STOP.
        current = _controllers(root)
        if {(p['pid'], p['starttime']) for p in current} != {
                (p['pid'], p['starttime']) for p in records}:
            raise RuntimeError('Controller tree changed; original processes resumed.')
        yield handles
    finally:
        for _, fd in handles:
            try:
                if fd in stopped:
                    _pidfd_send_signal(fd, signal.SIGCONT)
            except ProcessLookupError:
                pass  # Already replaced. pidfd prevents signalling another PID.
            finally:
                os.close(fd)


def _preflight(root):
    from kinbot.parameters import Parameters
    from kinbot.anl.automatic import validate_automatic_environment
    run_salvage._preflight_worker_routes(root)
    _, paths = run_salvage._cases(root)
    for path in paths:
        print(f'Checking replacement controller environment: {path.parent.name}', flush=True)
        validate_automatic_environment(Parameters(str(path)).par)


def _archive_controllers(root, archive):
    _, paths = run_salvage._cases(root)
    for path in paths:
        for directory in (path.parent, *sorted(path.parent.glob('[0-9]*'))):
            if not directory.is_dir():
                continue
            for name in ('driver.log', 'pes.log', 'pes_summary.txt', 'launch.json',
                         'pes_status.json', 'kinbot.log', 'kinbot.err', 'kinbot.out',
                         'kinbot_status.json'):
                run_salvage._save(directory / name, root, archive, move=True)
            run_salvage._save(directory / '.kinbot_post_reaction_branches.json', root, archive)


def refresh(root, *, apply=False):
    """Journal a controller-only handoff; repeated same-code calls are no-ops."""
    root = Path(root).resolve(strict=True)
    marker = root / 'salvage.json'
    history = json.loads(marker.read_text())
    journal_path = root / 'controller_refresh.json'
    journal = json.loads(journal_path.read_text()) if journal_path.exists() else None
    pin = run_salvage._pin_code(root, refresh=True)
    if journal and journal['phase'] == 'launched' and journal['code'] == str(pin):
        return journal
    _probe_process_handles()  # Includes resumed transactions, before any STOP.
    print('Process-handle preflight passed.', flush=True)
    if journal and journal['phase'] != 'launched' and journal['code'] != str(pin):
        archive = Path(journal['archive'])
        # The observed missing-Python-API failure left only a preparing
        # journal. No controller/log migration has been committed in it.
        # Permit replacement code only at that boundary, after fresh checks.
        if (journal['phase'] != 'preparing' or history['status'] != 'launched'
                or history['code'] != journal['previous_code']
                or any(p.name != 'previous_refresh.json' for p in archive.iterdir())):
            raise RuntimeError('An interrupted refresh already changed controllers/files; '
                               'finish it with its original code before updating again.')
        _controllers(root)
        _owned_jobs(root)
        _preflight(root)
        if apply:
            journal.setdefault('superseded_preparation_codes', []).append(journal['code'])
            journal['code'] = str(pin)
            atomic_json(journal_path, journal)
        else:
            return {**journal, 'phase': 'reviewed', 'replacement_code': str(pin)}
    if not journal or journal['phase'] == 'launched':
        if history['status'] != 'launched':
            raise RuntimeError('Controller refresh requires a launched salvage run.')
        records = _controllers(root)
        jobs = _owned_jobs(root)
        _preflight(root)  # Do this while the original controllers keep working.
        if not apply:
            return {'phase': 'reviewed', 'controllers': records, 'jobs': jobs, 'code': str(pin)}
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        archive = root.with_name(root.name + '.controllers.' + stamp)
        archive.mkdir()
        if journal:
            atomic_json(archive / 'previous_refresh.json', journal)
        journal = {'phase': 'preparing', 'code': str(pin), 'previous_code': history['code'],
                   'archive': str(archive), 'created_utc': stamp}
        atomic_json(journal_path, journal)
    elif not apply:
        return journal
    archive = Path(journal['archive'])
    if journal['phase'] in ('preparing', 'stopping'):
        records = _controllers(root)
        print('Quiescing only PES/KB controllers; Slurm jobs remain running.', flush=True)
        with _freeze(root, records) as handles:
            jobs = _owned_jobs(root)  # No thread can submit/write while checked.
            journal.update(phase='stopping', controllers=records, jobs=jobs)
            atomic_json(journal_path, journal)
            for _, fd in handles:
                try:
                    _pidfd_send_signal(fd, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            deadline = time.monotonic() + 5
            while run_cleanup._processes(root):
                if time.monotonic() >= deadline:
                    raise RuntimeError('Controller exit pending; repeat the same refresh command.')
                time.sleep(.05)
        journal['phase'] = 'stopped'
        atomic_json(journal_path, journal)
    if journal['phase'] == 'stopped':
        if run_cleanup._processes(root):
            raise RuntimeError('Unexpected local processes after controller stop; refusing edits.')
        _owned_jobs(root)
        # Hold the original inode locks; do not remove or replace run markers.
        with ExitStack() as locks:
            for path in root.glob('*/[0-9]*/.kinbot_run.json'):
                handle = locks.enter_context(path.open())
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            _archive_controllers(root, archive)
            history.update(status='prepared', code=str(pin))
            atomic_json(marker, history)
            journal['phase'] = 'prepared'
            atomic_json(journal_path, journal)
    if journal['phase'] == 'prepared':
        # Launch uses the same input, DBs and task ledgers. In particular,
        # submitted jobs stay submitted and are adopted by their recorded IDs.
        run_salvage.launch(root)
        journal['phase'] = 'launched'
        atomic_json(journal_path, journal)
    return journal


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--apply', action='store_true', help='Replace controllers and launch with current code.')
    args = parser.parse_args(argv)
    root = args.root.absolute()
    # The same lock as run_salvage prevents simultaneous migration/launch.
    with root.with_suffix('.salvage.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        result = refresh(root, apply=args.apply)
    print(json.dumps({k: result[k] for k in ('phase', 'code', 'archive') if k in result}, indent=2))
    print('Retained Slurm IDs: ' + ', '.join(j['job_id'] for j in result.get('jobs', [])))


if __name__ == '__main__':
    main()
