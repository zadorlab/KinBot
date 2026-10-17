"""Inspect or stop one Linux/Slurm calculation tree without user-wide kills.

Inspection is the default. ``--stop --archive`` stops processes and Slurm
jobs scoped to the given working directory, then moves the stopped tree to
a dated sibling. No calculation files or shared scratch trees are deleted.
"""
import argparse
from datetime import datetime, timezone
import getpass
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from kinbot.execution import atomic_json


class CleanupTimeout(RuntimeError):
    """Cancellation is requested, but the tree is not yet safe to move."""

    def __init__(self, report):
        self.report = report
        super().__init__(report['error'])


def _inside(value, root):
    try:
        Path(value).resolve().relative_to(root.resolve())
        return True
    except (ValueError, OSError):
        return False


def _jobs(root):
    result = subprocess.run(
        ['squeue', '-h', '-u', getpass.getuser(), '-o', '%i|%T|%R|%Z'],
        capture_output=True, text=True, check=True, timeout=30)
    jobs = []
    for line in result.stdout.splitlines():
        fields = line.split('|', 3)
        if len(fields) != 4:
            if line.strip():
                raise ValueError(f'Unrecognized squeue row: {line!r}.')
            continue
        identifier, state, reason, directory = (field.strip() for field in fields)
        if directory and _inside(directory, root):
            if not identifier.strip() or any(char not in '0123456789_[]%,-+.' for char in identifier.strip()):
                raise ValueError(f'Unrecognized Slurm job identifier: {identifier!r}.')
            jobs.append({'job_id': identifier, 'directory': directory,
                         'state': state, 'reason': reason})
    return jobs


def _processes(root, proc=Path('/proc'), *, known=()):
    if not proc.is_dir():
        raise RuntimeError('Run cleanup on the Linux host where the controllers ran.')
    records = {}
    for path in proc.iterdir():
        if not path.name.isdigit():
            continue
        try:
            if path.stat().st_uid != os.getuid():
                continue
            # comm may contain spaces and parentheses; stat fields start
            # after its closing parenthesis. starttime prevents PID reuse.
            fields = (path / 'stat').read_text().rsplit(')', 1)[1].split()
            if fields[0] == 'Z':
                continue
            records[int(path.name)] = {
                'pid': int(path.name), 'ppid': int(fields[1]),
                'starttime': fields[19], 'state': fields[0],
                'cwd': str((path / 'cwd').resolve(strict=True)),
                'command': (path / 'cmdline').read_bytes().replace(b'\0', b' ').decode(errors='replace').strip(),
            }
        except (OSError, ValueError, IndexError):
            continue
    excluded = {os.getpid()}
    ancestor = os.getpid()
    while ancestor in records and records[ancestor]['ppid'] not in excluded:
        ancestor = records[ancestor]['ppid']
        excluded.add(ancestor)
    identities = {(record['pid'], record['starttime']) for record in known}
    selected = {pid for pid, record in records.items()
                if pid not in excluded and (_inside(record['cwd'], root)
                    or (pid, record['starttime']) in identities)}
    while True:
        descendants = {pid for pid, record in records.items()
                       if pid not in excluded and record['ppid'] in selected}
        if descendants <= selected:
            break
        selected |= descendants
    return [records[pid] for pid in sorted(selected)]


def inspect(root, *, known=()):
    root = Path(root).expanduser().resolve(strict=True)
    if (not root.is_dir() or root == Path.home().resolve() or root == Path('/')
            or not (root / '.kinbot_run.json').is_file()):
        raise ValueError('Expected one KinBot run directory with .kinbot_run.json.')
    # Query first: a scheduler outage must not be mistaken for an empty queue.
    jobs = _jobs(root)
    return {'directory': str(root), 'jobs': jobs, 'processes': _processes(root, known=known)}


def _cancel(root, jobs):
    """A job can leave squeue between inventory and scancel."""
    identifiers = {job['job_id'] for job in jobs}
    try:
        subprocess.run(['scancel', *sorted(identifiers)], check=True,
                       capture_output=True, text=True, timeout=30)
    except subprocess.CalledProcessError as error:
        # A failed scheduler query still raises; absence is never inferred
        # from an unavailable controller or an unsuccessful cancellation.
        if identifiers.intersection(job['job_id'] for job in _jobs(root)):
            raise RuntimeError('Slurm cancellation failed: '
                               + (error.stderr or str(error)).strip()) from error


def stop(root, *, archive=False, timeout=300., known=()):
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('Cleanup timeout must be a finite positive number of seconds.')
    # A scoped launcher may be working outside the calculation directory.
    # Callers can carry its verified PID/starttime identity into the same
    # descendant tracking and PID-reuse checks used for workers.
    report = inspect(root, known=known) if known else inspect(root)
    root = Path(report['directory'])
    atomic_json(root / 'cleanup_inventory.json', report)
    known = {}
    signalled = set()
    cancelled = set()
    started = time.monotonic()
    next_notice = 0.
    remaining = report
    while True:
        elapsed = time.monotonic() - started
        remaining = dict(remaining, elapsed_seconds=round(elapsed, 2))
        if not remaining['jobs'] and not remaining['processes']:
            remaining['status'] = 'stopped'
            atomic_json(root / 'cleanup_status.json', remaining)
            break
        if elapsed >= timeout:
            remaining.update(status='stop_pending', error=(
                f'Scoped processes/jobs remain after {timeout:g} seconds; tree preserved. '
                'See cleanup_status.json for remaining jobs and processes. '
                'Inspect or repeat --stop with a longer --timeout before restarting or archiving.'))
            atomic_json(root / 'cleanup_status.json', remaining)
            raise CleanupTimeout(remaining)
        remaining['status'] = 'stopping'
        atomic_json(root / 'cleanup_status.json', remaining)
        if elapsed >= next_notice:
            jobs = ', '.join(f'{job["job_id"]} ({job.get("state", "unknown")})'
                             for job in remaining['jobs']) or 'none'
            pids = ', '.join(str(record['pid']) for record in remaining['processes']) or 'none'
            print(f'Stopping {root.name}: {elapsed:.0f}s; Slurm jobs: {jobs}; '
                  f'local PIDs: {pids}', file=sys.stderr, flush=True)
            next_notice = elapsed + 10.
        # Keep identities of descendants even if the parent exits and the
        # descendants are reparented while working in a separate scratch tree.
        for record in remaining['processes']:
            known[(record['pid'], record['starttime'])] = record
        for sig in (signal.SIGTERM, signal.SIGKILL):
            if sig == signal.SIGKILL and elapsed < min(5., timeout / 2):
                continue
            for record in remaining['processes']:
                identity = (record['pid'], record['starttime'], sig)
                if identity in signalled:
                    continue
                current = {item['pid']: item for item in _processes(root, known=known.values())}
                if current.get(record['pid'], {}).get('starttime') == record['starttime']:
                    try:
                        os.kill(record['pid'], sig)
                    except ProcessLookupError:
                        pass
                signalled.add(identity)
        # Catch work submitted between the initial inventory and controller
        # termination. Cancel only jobs whose current WorkDir is in this tree.
        pending = [job for job in remaining['jobs'] if job['job_id'] not in cancelled]
        if pending:
            _cancel(root, pending)
            cancelled.update(job['job_id'] for job in pending)
        # Accepted cancellation can remain COMPLETING during Slurm teardown.
        # Do not hammer scancel every poll while waiting for it to finish.
        time.sleep(min(2., max(0., timeout - (time.monotonic() - started))))
        remaining = inspect(root, known=known.values())
    report['status'] = 'stopped'
    if archive:
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        destination = root.with_name(root.name + '.stopped.' + stamp)
        root.rename(destination)
        report['archive'] = str(destination)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory')
    parser.add_argument('--stop', action='store_true')
    parser.add_argument('--archive', action='store_true')
    parser.add_argument('--timeout', type=float, default=300., metavar='SECONDS',
                        help='Wait for scoped jobs/processes to disappear (default: 300).')
    args = parser.parse_args(argv)
    if args.archive and not args.stop:
        parser.error('--archive requires --stop')
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error('--timeout must be a finite positive number of seconds')
    try:
        result = (stop(args.run_directory, archive=args.archive, timeout=args.timeout)
                  if args.stop else inspect(args.run_directory))
    except CleanupTimeout as error:
        print(json.dumps(error.report, indent=2))
        return 1
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
