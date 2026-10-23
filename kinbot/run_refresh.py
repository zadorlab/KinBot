"""Replace fresh-run or salvage controllers while retaining owned QC Slurm jobs.

Only verified PES/KB controllers (including CBH-reference controllers) are
replaceable. Native calculations and ambiguous submissions cause a refusal.
With --replace-pending, only jobs still PENDING are cancelled and restaged;
running allocations and completed calculations are retained.
"""
import argparse
from contextlib import ExitStack, contextmanager
from contextvars import ContextVar
import ctypes
from datetime import datetime, timezone
import errno
import fcntl
import json
import os
from pathlib import Path
import platform
import re
import shlex
import signal
import socket
import subprocess
import sys
import time

from kinbot.execution import atomic_json
from kinbot import run_cleanup, run_salvage

_POOL_GUARDED = ContextVar('refresh_pool_guarded', default=False)


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
    for inp in _ordinary_inputs(root):
        roles[inp.parent.resolve()] = ('kinbot.kb', 'kinbot_status.json')
    for directory, (_, name) in roles.items():
        for marker in (directory / name, directory / 'launch.json'):
            if marker.exists():
                host = json.loads(marker.read_text()).get('host')
                if host and host.split('.')[0] != socket.gethostname().split('.')[0]:
                    raise RuntimeError(f'Run controller refresh on {host}, not this host.')
    records = _replacement_processes(root)
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


def _replacement_processes(root):
    from kinbot.rotdpy_campaign import replacement_processes
    return replacement_processes(root, run_cleanup._processes(root))


@contextmanager
def _submission_guard(root):
    """Do not freeze a pool-lock holder, or race a retained sampler's sbatch."""
    from kinbot import run_policy
    policy = run_policy.find(root)
    if policy is None or _POOL_GUARDED.get():
        yield
        return
    # Bounded acquisition also works during recovery from an interrupted stop.
    with (policy.parent/'.kinbot_resource_policy.lock').open('a') as stream:
        deadline = time.monotonic() + 5
        while True:
            try:
                fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if time.monotonic() >= deadline:
                    raise RuntimeError('Pool submission is still in flight; repeat refresh after it finishes.')
                time.sleep(.05)
        token = _POOL_GUARDED.set(True)
        try:
            yield
        finally:
            _POOL_GUARDED.reset(token)


def _graphs(root):
    directories = set(root.glob('*/[0-9]*/anl_automatic/states/*/composite'))
    for inp in root.glob('*/[0-9]*/*.json'):
        if inp.name != inp.parent.name+'.json':
            continue
        parameters=json.loads(inp.read_text())
        base=inp.parent/parameters.get('anl_automatic_directory','anl_automatic')
        directories.update(base.glob('states/*/composite'))
    for pattern in ('*/[0-9]*/vrctst/molpro_scan/*',
                    '*/[0-9]*/vrctst/molpro_pivots/*',
                    '*/[0-9]*/vrctst/molpro/dispatch*'):
        directories.update(root.glob(pattern))
    for graph in sorted(directories):
        if not graph.is_dir():
            continue
        if not graph.resolve().is_relative_to(root.resolve()):
            raise RuntimeError(f'{graph}: dispatcher directory is outside this run; refusing ownership.')
        if not all((graph / name).is_file() for name in ('workflow.json', 'state.json')):
            raise RuntimeError(f'{graph}: graph preparation is incomplete; repeat later.')
        yield graph / 'state.json'


def _ordinary_inputs(root):
    prepared = json.loads((root/'prepared.json').read_text())
    inputs = [Path(p).resolve() for p in prepared['inputs']]
    if not inputs or any(not p.is_relative_to(root) for p in inputs):
        raise RuntimeError('Prepared inputs are outside the requested run.')
    result = set()
    for source in inputs:
        for inp in source.parent.glob('[0-9]*/*.json'):
            if inp.name != inp.parent.name+'.json':
                continue
            result.add(inp)
            settings = json.loads(inp.read_text())
            base = inp.parent/settings.get('anl_automatic_directory', 'anl_automatic')
            if not base.resolve().is_relative_to(root):
                raise RuntimeError('Automatic reference directory escapes the run.')
            result.update(base.glob('reference_l2/*/input.json'))
    return sorted(result)


def _ordinary_job_owners(root):
    """Read the persisted routes that replacement L1/L2 controllers adopt."""
    owners = {}
    # Ownership also runs while a journalled input migration is interrupted.
    # Its input hashes are validated by that migration, not against the old
    # prepared hash between its two atomic writes.
    for inp in _ordinary_inputs(root):
        settings = json.loads(inp.read_text())
        worker = inp.parent.resolve()
        route = Path(settings.get('profiled_route_file', '.kinbot_theory_jobs.json'))
        route = route if route.is_absolute() else worker/route
        if not route.resolve().is_relative_to(root):
            raise RuntimeError(f'{route}: ordinary job routes are outside this run.')
        if not route.is_file():
            continue
        payload = json.loads(route.read_text())
        if payload.get('schema') != 1 or not isinstance(payload.get('jobs'), dict):
            raise RuntimeError(f'{route}: invalid ordinary QC routes.')
        for name, entry in payload['jobs'].items():
            job = entry.get('job_id')
            if not job:
                continue
            target = (worker/name).resolve()
            if not target.is_relative_to(worker) or target == worker:
                raise RuntimeError(f'{route}: QC target is outside its worker.')
            job = str(job)
            if job in owners:
                raise RuntimeError(f'Slurm job {job} has duplicate ordinary QC ownership.')
            owners[job] = (worker, target)
    return owners


def _owned_jobs(root):
    with _submission_guard(root):
        return _owned_jobs_guarded(root)


def _owned_jobs_guarded(root):
    """Verify composite, VRC and profiled L1/L2 ownership without resetting it."""
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
    ordinary = _ordinary_job_owners(root)
    if owned.keys() & ordinary.keys():
        raise RuntimeError('An allocation has both dispatcher and ordinary QC ownership.')
    for job in jobs:
        ident = job['job_id']
        directory = Path(job['directory']).resolve()
        if ident in ordinary:
            worker, target = ordinary[ident]
            if (directory != worker or
                    not Path(str(target)+'.py').is_file() or
                    not Path(str(target)+'.sbatch').is_file()):
                raise RuntimeError(f'Job {ident}: ordinary QC directory/worker/script does not match its route.')
        elif owned.get(ident) != directory:
            from kinbot.rotdpy_campaign import owns_job
            if owns_job(root, job):
                continue
            raise RuntimeError(f'Live Slurm job {job["job_id"]} lacks matching '
                               'dispatcher, profiled QC or campaign ownership; nothing reset.')
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
        from kinbot.input_refresh import SNAPSHOT
        from kinbot.vrc_reference_review import validate_input_change
        snapshot = path.parent/SNAPSHOT
        if snapshot.exists():
            validate_input_change(json.loads(snapshot.read_text()), json.loads(path.read_text()), path.parent)


def _archive_controllers(root, archive):
    _, paths = run_salvage._cases(root)
    directories = {p.parent for p in _ordinary_inputs(root)}
    for path in paths:
        directories.add(path.parent)
        directories.update(p for p in path.parent.glob('[0-9]*') if p.is_dir())
    for directory in sorted(directories):
        for name in ('driver.log', 'driver.stdout', 'driver.stderr', 'driver.pid',
                     'pes.log', 'pes_summary.txt', 'launch.json', 'pes_status.json',
                     'kinbot.log', 'kinbot.err', 'kinbot.out', 'kinbot_status.json'):
            run_salvage._save(directory / name, root, archive, move=True)
        run_salvage._save(directory / '.kinbot_post_reaction_branches.json', root, archive)
        run_salvage._save(directory / 'chemids', root, archive)
        # This generation's owner record must not survive its controllers.
        # Startup checks it before creating the replacement generation, so a
        # truncated record after a filesystem outage otherwise prevents launch.
        # Preserve its exact bytes here, under the stopped-controller locks;
        # chemistry, discovery and task ledgers remain in their live locations.
        run_salvage._save(directory / 'pes_exploration.json', root, archive, move=True)


def _restore_discovery(root, archive):
    """Recover the discovery list lost by older PES startup code, once.

    The preceding handoff's summary predates that truncation. Current order
    remains a prefix so even an interrupted handoff cannot duplicate workers.
    Called with all controllers stopped, before launching the replacement.
    """
    previous = archive/'previous_refresh.json'
    if not previous.is_file():
        return {}
    record = json.loads(previous.read_text())
    source = Path(record['archive'])
    if (record.get('phase') != 'launched' or source.is_symlink()
            or source.parent != root.parent
            or not source.name.startswith(root.name+'.controllers.')):
        raise ValueError('Previous controller archive does not belong to this run.')
    from kinbot.pes_inventory import restore
    _, inputs = run_salvage._cases(root)
    reports = {}
    for path in inputs:
        summary = source/path.parent.relative_to(root)/'pes_summary.txt'
        if summary.is_file():
            reports[path.parent.name] = restore(path.parent, summary, apply=True)
    return reports


@contextmanager
def _worker_locks(root, *, timeout=30.):
    """Wait for controller-lock release without replacing the locked inode."""
    deadline = time.monotonic() + timeout
    with ExitStack() as locks:
        paths = set(root.glob('*/[0-9]*/.kinbot_run.json'))
        paths.update(p.parent/'.kinbot_run.json' for p in _ordinary_inputs(root) if (p.parent/'.kinbot_run.json').exists())
        for path in sorted(paths):
            handle = locks.enter_context(path.open('r+'))
            noticed = False
            while True:
                try:
                    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except OSError as error:
                    if error.errno not in (errno.EAGAIN, errno.EACCES):
                        raise
                    remaining = deadline - time.monotonic()
                    if remaining <= 0:
                        raise RuntimeError(
                            f'Worker lock still held: {path}. Controllers are stopped; '
                            'no allocations or inputs were changed at this step. '
                            'Repeat the same refresh command after the lock is released; '
                            'do not delete the run marker.') from error
                    if not noticed:
                        print(f'Waiting for worker lock release: {path.relative_to(root)}', flush=True)
                        noticed = True
                    time.sleep(min(.25, remaining))
        yield


def _unchanged_stopped_handoff(root, journal):
    """Allow fixed code at the lock boundary, before any scheduler/input edits."""
    if journal.get('phase') != 'stopped' or _replacement_processes(root):
        return False
    if any(journal.get(key) for key in (
            'pending_holds', 'pending_cancel_acknowledged', 'pending_plans',
            'pending_cancellation_evidence')):
        return False
    network = journal.get('network_revision') or {}
    if network.get('cancellations') or network.get('archive'):
        return False
    if network:
        if run_salvage.file_sha256(root/'prepared.json') != network['prepared_sha256']:
            return False
        for edit in network['inputs']:
            if run_salvage.file_sha256(root/edit['path']) != edit['before']:
                return False
    migration = journal.get('vrc_migration') or {}
    if migration.get('cancellations'):
        return False
    duplicate = journal.get('vrc_deduplication') or {}
    if duplicate.get('cancellations') or duplicate.get('phase') == 'archived':
        return False
    if migration:
        if run_salvage.file_sha256(root/'prepared.json') != migration['prepared_sha256']:
            return False
        for edit in migration['inputs']:
            if run_salvage.file_sha256(root/edit['path']) != edit['before']:
                return False
    return True


def _pending_preflight():
    result = subprocess.run(['scancel', '--help'], capture_output=True, text=True, timeout=10)
    if result.returncode:
        raise RuntimeError('Cannot inspect scancel support; no controllers or allocations were changed.')
    # Even some clients that advertise --ctld filter states client-side.
    # Use hold/check/cancel on all versions, not a help-based atomicity guess.
    result = subprocess.run(['scontrol', '--help'], capture_output=True, text=True, timeout=10)
    help_text = result.stdout + result.stderr
    if result.returncode or not all(re.search(r'\b' + word + r'\b', help_text, re.I)
                                   for word in ('hold', 'release')):
        raise RuntimeError('Pending replacement requires scontrol hold/release on this Slurm client; '
                           'no controllers or allocations were changed.')
    return 'hold_then_cancel'


def _pending_request(args):
    """Return a diagnostic; the scheduler record determines mutation outcome."""
    env = {k: v for k, v in os.environ.items() if not k.startswith('SCANCEL_')}
    try:
        result = subprocess.run(args, capture_output=True, text=True, timeout=30, env=env)
    except (OSError, subprocess.TimeoutExpired) as error:
        return str(error)
    return (result.stderr or result.stdout).strip() if result.returncode else ''


def _cancel_held_pending(candidates, journal, marker):
    """Older Slurm: hold, verify PENDING/priority zero, then cancel by ID.

    A hold does not stop a running job. If scheduling won that race, release
    only the hold we introduced and retain the allocation. Persist intent
    before each mutation so a timeout/interruption can be reconciled by the
    same command. No task files are changed here.
    """
    from kinbot.execution import slurm_job_record

    def record_for(job):
        record = slurm_job_record(job['job_id'])
        if (record.get('JobId') != job['job_id'] or not record.get('JobState')
                or Path(record.get('WorkDir', '/')).resolve() != Path(job['directory']).resolve()):
            raise RuntimeError(f'Job {job["job_id"]}: scheduler identity/state unavailable; '
                               'no task was reset. Repeat the same refresh command.')
        return record

    records = journal.setdefault('pending_holds', {})
    for job in candidates:
        ident = job['job_id']
        saved = records.get(ident)
        if saved and saved.get('phase') in ('retained', 'cancelled'):
            continue
        current = record_for(job)
        if saved is None:
            saved = records[ident] = {'before': current, 'phase': 'inspected'}
            atomic_json(marker, journal)
        if current['JobState'] == 'PENDING' and saved['phase'] != 'releasing':
            if not str(current.get('Priority', '')).isdecimal():
                raise RuntimeError(f'Job {ident}: pending priority unavailable; no cancellation attempted.')
            if current['Priority'] != '0':
                saved['phase'] = 'holding'
                atomic_json(marker, journal)
                detail = _pending_request(['scontrol', 'hold', ident])
                current = record_for(job)
                if current['JobState'] == 'PENDING' and current.get('Priority') != '0':
                    raise RuntimeError(f'Job {ident}: hold not confirmed ({detail}); '
                                       'no task was reset. Repeat the same refresh command.')
            if current['JobState'] == 'PENDING' and current.get('Priority') == '0':
                saved.update(phase='cancelling', held_record=current)
                atomic_json(marker, journal)
                detail = _pending_request(['scancel', ident])
                current = record_for(job)
                if current['JobState'] == 'PENDING':
                    raise RuntimeError(f'Job {ident}: cancellation not confirmed ({detail}); '
                                       'job remains held. Repeat the same refresh command.')
        if current['JobState'] == 'CANCELLED' and saved['phase'] == 'cancelling':
            saved.update(phase='cancelled', cancellation_record=current)
        else:
            # The job started before hold reached slurmctld. Holding RUNNING
            # changes its requeue priority, not its execution; undo that hold.
            # Never release a hold that existed before this refresh.
            if (saved['phase'] in ('holding', 'cancelling', 'releasing')
                    and saved['before'].get('Priority') != '0'
                    and current.get('Priority') == '0'
                    and current['JobState'] in ('PENDING', 'RUNNING', 'CONFIGURING', 'COMPLETING', 'SUSPENDED')):
                saved['phase'] = 'releasing'
                atomic_json(marker, journal)
                detail = _pending_request(['scontrol', 'release', ident])
                current = record_for(job)
                if current.get('Priority') == '0' and current['JobState'] in (
                        'PENDING', 'RUNNING', 'CONFIGURING', 'COMPLETING', 'SUSPENDED'):
                    raise RuntimeError(f'Job {ident}: releasing refresh hold not confirmed ({detail}); '
                                       'repeat the same refresh command. Allocation was not cancelled.')
            saved['phase'] = 'retained'
        atomic_json(marker, journal)


def _replace_pending(root, journal, marker):
    """Cancel pending jobs without a start/cancel race, then journal restaging.

    Never infer cancellation solely from queue absence. Require the retained
    Slurm job record to confirm CANCELLED, or preserve a completed execution.
    """
    from kinbot.execution import slurm_job_record
    from kinbot.run_recover import _apply_graph
    from kinbot.anl import dispatch
    if 'pending_plans' not in journal:
        if 'pending_candidates' not in journal:
            journal['pending_candidates'] = [j for j in _owned_jobs(root) if j['state'] == 'PENDING']
            atomic_json(marker, journal)
        candidates = journal['pending_candidates']
        _check_pending_scope(root, candidates)
        ids = [j['job_id'] for j in candidates]
        if any(not ident.isdecimal() for ident in ids):
            raise RuntimeError('Pending replacement requires individual numeric Slurm job IDs.')
        if ids and not journal.get('pending_cancel_acknowledged'):
            _cancel_held_pending(candidates, journal, marker)
            journal['pending_cancel_acknowledged'] = True
            atomic_json(marker, journal)
        live = {j['job_id']: j for j in _owned_jobs(root)}
        reset = {}
        evidence = {}
        for job in candidates:
            ident = job['job_id']
            if ident in live:
                if (live[ident]['state'] == 'PENDING' and journal.get(
                        'pending_holds', {}).get(ident, {}).get('phase') != 'retained'):
                    raise RuntimeError('Pending cancellation is not yet visible; repeat the same refresh command.')
                continue  # Includes RUNNING/CONFIGURING/COMPLETING: retain ownership.
            directory = Path(job['directory']).resolve()
            outcome = directory/'execution.json'
            if outcome.exists():
                # It started/finished before cancellation filtering. Let the
                # ordinary dispatcher validate it, not this cancellation path.
                continue
            record = journal.get('pending_holds', {}).get(ident, {}).get(
                'cancellation_record') or slurm_job_record(ident)
            if (record.get('JobState') != 'CANCELLED'
                    or record.get('JobId') != ident
                    or Path(record.get('WorkDir', '/')).resolve() != directory):
                raise RuntimeError(f'Job {ident}: cancellation outcome is not confirmed; '
                                   'no duplicate will be submitted. Repeat the same refresh command.')
            evidence[ident] = record
            reset.setdefault(directory.parent.parent, []).append(directory.name)
        plans = []
        for graph, selected in sorted(reset.items()):
            _, spec, state = dispatch._load(graph)
            plans.append({'graph': str(graph), 'spec': spec, 'state': state,
                'before': {name: run_salvage.file_sha256(graph/name) for name in ('state.json', 'workflow.json')},
                'reset': selected, 'reasons': {n: 'Pending allocation cancelled for current-code restart.' for n in selected},
                'complete': [n for n, e in state['tasks'].items() if e['status'] == 'complete'],
                'active': {n: e['job_id'] for n, e in state['tasks'].items()
                           if e['status'] == 'submitted' and n not in selected}})
        journal.update(pending_plans=plans, pending_cancellation_evidence=evidence)
        atomic_json(marker, journal)
    # Resume a partially committed archive before reading the changed ledgers.
    for plan in journal['pending_plans']:
        _apply_graph(plan, root, Path(journal['archive'])/'pending_allocations')
    journal['jobs'] = _owned_jobs(root)
    atomic_json(marker, journal)


def _check_pending_scope(root, jobs):
    ordinary = _ordinary_job_owners(root)
    if any(j['state'] == 'PENDING' and j['job_id'] in ordinary for j in jobs):
        raise RuntimeError('Ordinary L1/L2 allocations are retained by controller refresh. '
                           'Omit --replace-pending; use live QC pending-only requests '
                           'for an explicit ordinary allocation replacement.')


def refresh(root, *, apply=False, replace_pending=False, vrc_updates=None,
            cancel_superseded_vrc_scans=False, deduplicate_vrc=False, network_options=None):
    """Journal a controller-only handoff; repeated same-code calls are no-ops."""
    root = Path(root).resolve(strict=True)
    marker = root / 'salvage.json'
    history = run_salvage._controller_history(root)
    journal_path = root / 'controller_refresh.json'
    journal = json.loads(journal_path.read_text()) if journal_path.exists() else None
    vrc_updates = vrc_updates or {}
    if network_options and (vrc_updates or replace_pending or deduplicate_vrc):
        raise ValueError("Network revision must not be mixed with VRC/resource migration.")
    if (journal and journal["phase"] != "launched"
            and journal.get("network_options") != network_options):
        raise RuntimeError("Repeat the interrupted refresh with the same network options.")
    if deduplicate_vrc and (vrc_updates or replace_pending):
        raise ValueError('VRC deduplication cannot be combined with a scientific migration or pending replacement.')
    if (journal and journal['phase'] != 'launched'
            and bool(journal.get('deduplicate_vrc')) != deduplicate_vrc):
        raise RuntimeError('Repeat the interrupted refresh with the same VRC deduplication option.')
    if vrc_updates and replace_pending:
        raise ValueError('Targeted VRC migration cannot replace unrelated pending allocations.')
    if cancel_superseded_vrc_scans and not vrc_updates:
        raise ValueError('Cancelling superseded scans requires explicit VRC settings.')
    pin = run_salvage._pin_code(root, refresh=True)
    if (journal and journal['phase'] == 'launched' and journal['code'] == str(pin)
            and journal.get('network_options') == network_options
            and (not replace_pending or journal.get('replace_pending'))
            and (not deduplicate_vrc or journal.get('deduplicate_vrc'))
            and (not vrc_updates or journal.get('vrc_migration', {}).get('settings') == vrc_updates)):
        return journal
    if (journal and journal['phase'] != 'launched'
            and vrc_updates != journal.get('vrc_migration', {}).get('settings', {})):
        raise RuntimeError('Repeat the interrupted refresh with the same VRC settings.')
    pending_method = None
    if replace_pending or (journal and journal.get('replace_pending') and journal['phase'] != 'launched'):
        pending_method = _pending_preflight()
        print('Pending replacement: hold, verify, cancel; running allocations retained.', flush=True)
    _probe_process_handles()  # Includes resumed transactions, before any STOP.
    print('Process-handle preflight passed.', flush=True)
    if journal and journal['phase'] != 'launched' and journal['code'] != str(pin):
        archive = Path(journal['archive'])
        # Preparation failures and a stopped-controller lock failure can
        # accept fixed code only before any scheduler/input/log migration.
        unchanged = (journal['phase'] == 'preparing'
                     or _unchanged_stopped_handoff(root, journal))
        if (not unchanged or history['status'] != 'launched'
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
        if replace_pending:
            _check_pending_scope(root, jobs)
        _preflight(root)  # Do this while the original controllers keep working.
        migration = None
        duplicate = None
        network = None
        if network_options:
            from kinbot import network_revisit
            network = network_revisit.plan(root, jobs, network_options)
            print(f"Network revisit: {sum(len(w['rejected']) for w in network['workers'])} rejected channels; "
                  f"{len(network['selected_jobs'])} verified duplicate ANL allocations selected.", flush=True)
            local = sum(len(e['mapping']['aliases']) for e in network['aliases'])
            shared = sum(len(e['mapping']['aliases']) for e in network.get('external_ts', []))
            retired = {e['graph'] for group in [*network['aliases'], *network.get('external_ts', [])]
                       for e in group['remove']}
            print(f'ANL reuse mappings: local={local}; across wells={shared}; '
                  f'existing duplicate graphs to archive={len(retired)}. '
                  'Allocation cancellations are a separate count.', flush=True)
            recovered = sum(len(e['mapping'].get('frequency_migrations', {}))
                            for e in network.get('external_ts', []))
            if recovered:
                print(f'ANL reference frequency provenance: {recovered} verified native-row '
                      'recoveries recorded; original models and QC files retained.', flush=True)
        if deduplicate_vrc:
            from kinbot import vrc_deduplicate
            duplicate = vrc_deduplicate.plan(root, jobs)
            print(f'VRC deduplication: retain {len(duplicate["keep"])} original graphs; '
                  f'archive {len(duplicate["remove"])} duplicate/derived graphs; '
                  f'cancel {len(duplicate["selected_jobs"])} duplicate allocations.', flush=True)
        if vrc_updates:
            from kinbot import vrc_restart
            migration = vrc_restart.plan(root, vrc_updates)
            scans, selected = (vrc_restart.scan_inventory(root, jobs)
                               if migration['replace_geometry'] else ([], []))
            if selected and not cancel_superseded_vrc_scans:
                raise ValueError('Live old scans exist; use --cancel-superseded-vrc-scans '
                                 'to retain their files and free their allocations.')
            migration.update(scans=scans, selected_jobs=selected)
        if not apply:
            return {'phase': 'reviewed', 'controllers': records, 'jobs': jobs, 'code': str(pin),
                    'replace_pending': replace_pending, 'vrc_migration': migration,
                    'vrc_deduplication': duplicate, 'network_revision': network}
        if not marker.exists():
            # Bootstrap the normal handoff journal from verified fresh-launch
            # receipts before archiving any of them. A partial archive must be
            # resumable even if launch.json has already moved into the archive.
            atomic_json(marker, history)
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        archive = root.with_name(root.name + '.controllers.' + stamp)
        archive.mkdir()
        if journal:
            atomic_json(archive / 'previous_refresh.json', journal)
        journal = {'phase': 'preparing', 'code': str(pin), 'previous_code': history['code'],
                   'archive': str(archive), 'created_utc': stamp, 'replace_pending': replace_pending,
                   'pending_cancel_method': pending_method, 'deduplicate_vrc': deduplicate_vrc,
                   'network_options': network_options}
        if network is not None:
            journal['network_revision'] = network
        if duplicate is not None:
            journal['vrc_deduplication'] = duplicate
        if migration is not None:
            journal['vrc_migration'] = migration
            journal['cancel_superseded_vrc_scans'] = cancel_superseded_vrc_scans
        atomic_json(journal_path, journal)
    elif not apply:
        return journal
    archive = Path(journal['archive'])
    if journal['phase'] in ('preparing', 'stopping'):
        records = _controllers(root)
        print('Quiescing only PES/KB controllers; Slurm jobs remain running.', flush=True)
        with _submission_guard(root), _freeze(root, records) as handles:
            jobs = _owned_jobs(root)  # No thread can submit/write while checked.
            if journal.get('replace_pending'):
                _check_pending_scope(root, jobs)
            if journal.get('network_options'):
                from kinbot import network_revisit
                journal['network_revision'] = network_revisit.plan(root,jobs,journal['network_options'])
            if journal.get('deduplicate_vrc'):
                from kinbot import vrc_deduplicate
                journal['vrc_deduplication'] = vrc_deduplicate.plan(root, jobs)
            if journal.get('vrc_migration'):
                from kinbot import vrc_restart
                # Re-inventory after quiescing: a point may have finished and
                # the next one started during the read-only preflight.
                migration = vrc_restart.plan(root, journal['vrc_migration']['settings'])
                scans, selected = (vrc_restart.scan_inventory(root, jobs)
                                   if migration['replace_geometry'] else ([], []))
                if selected and not journal.get('cancel_superseded_vrc_scans'):
                    raise ValueError('New live scan detected; explicit targeted cancellation is required.')
                journal['vrc_migration'] = dict(migration, scans=scans, selected_jobs=selected)
            from kinbot import search_state
            search_checkpoints = []
            for inp in _ordinary_inputs(root):
                legacy = search_state.capture_legacy(inp, publish=False)
                if legacy is not None:
                    relative = (inp.parent/search_state.NAME).relative_to(root)
                    # Publish to the live worker only after its controller is
                    # gone. An aborted freeze must not install a stale record
                    # underneath an old controller that cannot maintain it.
                    atomic_json(archive/relative, legacy)
                    search_checkpoints.append(str(relative))
                for name in (search_state.NAME, 'kinbot_monitor.out'):
                    run_salvage._save(inp.parent/name, root, archive)
            journal.update(phase='stopping', controllers=records, jobs=jobs,
                           search_checkpoints=search_checkpoints)
            atomic_json(journal_path, journal)
            for _, fd in handles:
                try:
                    _pidfd_send_signal(fd, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            deadline = time.monotonic() + 5
            while _replacement_processes(root):
                if time.monotonic() >= deadline:
                    raise RuntimeError('Controller exit pending; repeat the same refresh command.')
                time.sleep(.05)
        journal['phase'] = 'stopped'
        atomic_json(journal_path, journal)
    if journal['phase'] == 'stopped':
        if _replacement_processes(root):
            raise RuntimeError('Unexpected local processes after controller stop; refusing edits.')
        # Disappearing from /proc need not coincide with filesystem lock
        # availability. Hold all original inodes before changing anything.
        with _worker_locks(root):
            for relative in journal.get('search_checkpoints', []):
                if not (root/relative).exists():
                    atomic_json(root/relative, json.loads((archive/relative).read_text()))
            if journal.get('replace_pending'):
                _replace_pending(root, journal, journal_path)
            else:
                _owned_jobs(root)
            if journal.get('deduplicate_vrc'):
                from kinbot import vrc_deduplicate
                vrc_deduplicate.apply(root, journal['vrc_deduplication'],
                    archive/'vrc_duplicates', journal, journal_path)
                journal['jobs'] = _owned_jobs(root)
                atomic_json(journal_path, journal)
            if journal.get('vrc_migration'):
                from kinbot import vrc_restart
                vrc_restart.cancel_scans(journal['vrc_migration'], journal, journal_path)
                vrc_restart.apply_inputs(root, journal['vrc_migration'], archive/'vrc_inputs')
                journal['jobs'] = _owned_jobs(root)
                atomic_json(journal_path, journal)
            if journal.get('network_revision'):
                from kinbot import network_revisit
                network_revisit.apply(root,journal['network_revision'],archive/'network_revision',journal,journal_path)
                journal['jobs'] = _owned_jobs(root)
                atomic_json(journal_path,journal)
            if journal.get('vrc_migration') or journal.get('network_revision'):
                from kinbot.search_state import rebase
                for inp in _ordinary_inputs(root):
                    rebase(inp, str(archive), reassess=bool(journal.get('network_revision')))
            _archive_controllers(root, archive)
            journal['discovery_restored'] = _restore_discovery(root, archive)
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
    parser.add_argument('--replace-pending', action='store_true',
                        help='Archive and replace pending allocations; preserve running allocations.')
    parser.add_argument('--deduplicate-vrc', action='store_true',
                        help='Verify original VRC results, cancel and archive resource-only duplicate chains.')
    parser.add_argument('--vdw-detection', type=float)
    parser.add_argument('--revisit-rejected', action='store_true')
    parser.add_argument('--deduplicate-ts', '--deduplicate-anl', action='store_true', help='Verify and share equivalent ANL stable/TS graphs within each PES.')
    parser.add_argument('--barrier-threshold-l2', type=float, help='L2 ceiling relative to the initial PES well, kcal/mol.')
    parser.add_argument('--vrc-scan-basis', help='Change only the VRC relaxed/asymptotic geometry basis.')
    parser.add_argument('--vrc-fragment-geometry', choices=['legacy', 'asymptotic_scan'])
    parser.add_argument('--vrc-pivot-level', choices=['L2', 'custom'])
    parser.add_argument('--vrc-direct-caspt2-cbs', action='store_true',
                        help='Use direct CASPT2-F12 aDZ/aTZ CBS for every channel; retain scans.')
    parser.add_argument('--vrc-reviewed-models', type=Path,
                        help='Install a hash-bound orbital model proposal in the same handoff.')
    parser.add_argument('--vrc-reference-basis', nargs='+',
                        help='One CC reference basis, or two adjacent bases for CBS extrapolation.')
    parser.add_argument('--vrc-reference-power', type=float,
                        help='CBS exponent for the CC reference (requires two --vrc-reference-basis values).')
    parser.add_argument('--cancel-superseded-vrc-scans', action='store_true',
                        help='Cancel only old relaxed-scan jobs; preserve their native files and all other jobs.')
    args = parser.parse_args(argv)
    root = args.root.absolute()
    # The same lock as run_salvage prevents simultaneous migration/launch.
    with root.with_suffix('.salvage.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        from kinbot.vrc_restart import settings
        from kinbot.network_revisit import settings as network_settings
        network_options = (network_settings(args.vdw_detection,args.revisit_rejected,args.deduplicate_ts,args.barrier_threshold_l2)
                           if args.vdw_detection is not None or args.revisit_rejected or args.deduplicate_ts or args.barrier_threshold_l2 is not None else None)
        result = refresh(root, apply=args.apply, replace_pending=args.replace_pending,
            vrc_updates=settings(args.vrc_scan_basis, args.vrc_fragment_geometry, args.vrc_pivot_level,
                                 args.vrc_reference_basis, args.vrc_reference_power,
                                 direct_caspt2_cbs=args.vrc_direct_caspt2_cbs,
                                 reviewed_models=args.vrc_reviewed_models),
            cancel_superseded_vrc_scans=args.cancel_superseded_vrc_scans,
            deduplicate_vrc=args.deduplicate_vrc, network_options=network_options)
    print(json.dumps({k: result[k] for k in ('phase', 'code', 'archive') if k in result}, indent=2))
    print('Retained Slurm IDs: ' + ', '.join(j['job_id'] for j in result.get('jobs', [])))
    if result.get('network_revision'):
        revision = result['network_revision']
        print('Network revisit: ' + json.dumps(revision['settings'],sort_keys=True))
        print('Rejected channels selected: '+str(sum(len(w['rejected']) for w in revision['workers'])))
        print('Verified duplicate ANL allocation IDs selected: '+
              ', '.join(j['job_id'] for j in revision['selected_jobs']))
        print('QC graphs retained except archived, verified ANL aliases; existing VRC paths retained.')
    if result.get('pending_cancellation_evidence'):
        print('Cancelled pending IDs (archived for replacement): ' +
              ', '.join(result['pending_cancellation_evidence']))
    if result.get('vrc_migration'):
        migration = result['vrc_migration']
        print('VRC changes: ' + json.dumps(migration['settings'], sort_keys=True))
        if migration['settings'].get('direct_caspt2_cbs') and result['phase']=='launched':
            from kinbot.vrc_restart import audit_direct_inputs
            print('Verified current VRC inputs: '+json.dumps(audit_direct_inputs(root),sort_keys=True))
            if migration['settings'].get('reviewed_models'):
                print('Reviewed orbital proposal installed in the same controller handoff.')
        print('Superseded scan IDs (native files retained): ' +
              ', '.join(j['job_id'] for j in migration['selected_jobs']))
    if result.get('vrc_deduplication'):
        repair = result['vrc_deduplication']
        print('VRC original graphs retained: ' + str(len(repair['keep'])))
        print('Duplicate/derived graphs archived: ' + str(len(repair['remove'])))
        print('Duplicate allocation IDs selected: ' + ', '.join(j['job_id'] for j in repair['selected_jobs']))
        if repair.get('archive'):
            print('VRC duplicate archive: ' + repair['archive'])


if __name__ == '__main__':
    main()
