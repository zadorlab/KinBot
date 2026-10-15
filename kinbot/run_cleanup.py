"""Inspect or stop one Linux/Slurm calculation tree without user-wide kills.

Inspection is the default. ``--stop --archive`` stops processes and Slurm
jobs scoped to the given working directory, then moves the stopped tree to
a dated sibling. No calculation files or shared scratch trees are deleted.
"""
import argparse
from datetime import datetime, timezone
import getpass
import json
import os
from pathlib import Path
import signal
import subprocess
import time


def _inside(value, root):
    try:
        Path(value).resolve().relative_to(root.resolve())
        return True
    except (ValueError, OSError):
        return False


def _jobs(root):
    result = subprocess.run(
        ['squeue', '-h', '-u', getpass.getuser(), '-o', '%i|%Z'],
        capture_output=True, text=True, check=True, timeout=30)
    jobs = []
    for line in result.stdout.splitlines():
        identifier, separator, directory = line.partition('|')
        if separator and directory.strip() and _inside(directory.strip(), root):
            if not identifier.strip() or any(char not in '0123456789_[]%,-+.' for char in identifier.strip()):
                raise ValueError(f'Unrecognized Slurm job identifier: {identifier!r}.')
            jobs.append({'job_id': identifier.strip(), 'directory': directory.strip()})
    return jobs


def _processes(root, proc=Path('/proc')):
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
                'starttime': fields[19],
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
    selected = {pid for pid, record in records.items()
                if pid not in excluded and _inside(record['cwd'], root)}
    while True:
        descendants = {pid for pid, record in records.items()
                       if pid not in excluded and record['ppid'] in selected}
        if descendants <= selected:
            break
        selected |= descendants
    return [records[pid] for pid in sorted(selected)]


def inspect(root):
    root = Path(root).expanduser().resolve(strict=True)
    if (not root.is_dir() or root == Path.home().resolve() or root == Path('/')
            or not (root / '.kinbot_run.json').is_file()):
        raise ValueError('Expected one KinBot run directory with .kinbot_run.json.')
    # Query first: a scheduler outage must not be mistaken for an empty queue.
    jobs = _jobs(root)
    return {'directory': str(root), 'jobs': jobs, 'processes': _processes(root)}


def stop(root, *, archive=False, timeout=30.):
    report = inspect(root)
    root = Path(report['directory'])
    (root / 'cleanup_inventory.json').write_text(json.dumps(report, indent=2) + '\n')
    for record in report['processes']:
        current = {item['pid']: item for item in _processes(root)}
        if current.get(record['pid'], {}).get('starttime') == record['starttime']:
            try:
                os.kill(record['pid'], signal.SIGTERM)
            except ProcessLookupError:
                pass
    if report['jobs']:
        subprocess.run(['scancel', *[job['job_id'] for job in report['jobs']]], check=True, timeout=30)
    started = time.monotonic()
    while True:
        remaining = inspect(root)
        if not remaining['jobs'] and not remaining['processes']:
            break
        elapsed = time.monotonic() - started
        if elapsed >= timeout:
            raise RuntimeError('Scoped processes/jobs remain; tree preserved. Inspect again before archiving.')
        # Catch work submitted between the initial inventory and controller
        # termination. Cancel only jobs whose current WorkDir is in this tree.
        if remaining['jobs']:
            subprocess.run(['scancel', *[job['job_id'] for job in remaining['jobs']]],
                           check=True, timeout=30)
        if elapsed >= min(5., timeout / 2):
            for record in remaining['processes']:
                current = {item['pid']: item for item in _processes(root)}
                if current.get(record['pid'], {}).get('starttime') == record['starttime']:
                    try:
                        os.kill(record['pid'], signal.SIGKILL)
                    except ProcessLookupError:
                        pass
        time.sleep(.5)
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
    args = parser.parse_args(argv)
    if args.archive and not args.stop:
        parser.error('--archive requires --stop')
    result = stop(args.run_directory, archive=args.archive) if args.stop else inspect(args.run_directory)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
