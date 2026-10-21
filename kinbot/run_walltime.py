"""Inspect or extend Slurm time limits without touching QC files or controllers.

The default is read-only. Running-job increases can require a Slurm admin;
always verify the scheduler's response rather than treating a request as success.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import pwd
import re
import shlex
import subprocess
import time

_LIMITS = (0., {})
_DEFAULTS = (0., {})


def seconds(value):
    if value.upper() in ('UNLIMITED', 'INFINITE'):
        return math.inf
    if not re.fullmatch(r'(?:\d+-)?\d+(?::\d{1,2}){0,2}', value):
        raise ValueError(f'Unrecognized Slurm time: {value!r}')
    days, clock = value.split('-', 1) if '-' in value else ('0', value)
    parts = [int(x) for x in clock.split(':')]
    if '-' in value:
        hours, minutes, secs = (parts + [0, 0])[:3]
    elif len(parts) == 3:
        hours, minutes, secs = parts
    elif len(parts) == 2:
        hours, minutes, secs = 0, *parts
    else:
        hours, minutes, secs = 0, parts[0], 0
    return ((int(days)*24 + hours)*60 + minutes)*60 + secs


def duration(value):
    if value is None:
        return 'unknown'
    if not math.isfinite(value):
        return 'UNLIMITED'
    minutes = max(0, math.ceil(value/60))
    days, minutes = divmod(minutes, 1440)
    hours, minutes = divmod(minutes, 60)
    return f'{days}-{hours:02d}:{minutes:02d}:00'


def validate_minimum(value):
    if not isinstance(value, str) or not re.fullmatch(r'(?:\d+-)?\d+:\d{2}:\d{2}', value):
        raise ValueError('minimum_walltime must be HH:MM:SS or D-HH:MM:SS.')
    if seconds(value) <= 0 or any(int(x) > 59 for x in value.split(':')[-2:]):
        raise ValueError('minimum_walltime must be a positive valid duration.')
    return value


def command(args):
    return subprocess.run(args, capture_output=True, text=True, check=True, timeout=20)


def partition_limits(*, cached=False):
    global _LIMITS
    if cached and time.monotonic()-_LIMITS[0] < 60 and _LIMITS[1]:
        return _LIMITS[1]
    limits = {}
    for line in command(['sinfo', '-h', '-o', '%P|%l']).stdout.splitlines():
        name, value = line.strip().split('|')
        limit = seconds(value)
        key = name.rstrip('*')
        limits[key] = min(limits.get(key, math.inf), limit)
        if name.endswith('*'):
            limits[''] = limits[key]
    if not limits:
        raise RuntimeError('Slurm returned no partition time limits.')
    _LIMITS = time.monotonic(), limits
    return limits


def partition_walltime(partition=''):
    """Resolve maximum time for every eligible partition; zero means unlimited."""
    limits = partition_limits(cached=True)
    try:
        cap = min(limits[name] for name in partition.split(','))
    except KeyError as error:
        raise RuntimeError(f'Slurm did not report the selected partition: {partition!r}') from error
    if cap <= 0:
        raise RuntimeError(f'Partition {partition!r} has no positive time limit.')
    if math.isinf(cap):
        return '00:00:00'
    days, remainder = divmod(int(cap), 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, secs = divmod(remainder, 60)
    clock = f'{hours:02d}:{minutes:02d}:{secs:02d}'
    return f'{days}-{clock}' if days else clock


def partition_defaults():
    """Used only when the script/environment leaves the time to Slurm."""
    global _DEFAULTS
    if time.monotonic()-_DEFAULTS[0] < 60 and _DEFAULTS[1]:
        return _DEFAULTS[1]
    result = {}
    for line in command(['sinfo', '-h', '-o', '%P|%L']).stdout.splitlines():
        name, value = line.strip().split('|')
        result[name.rstrip('*')] = seconds(value)
        if name.endswith('*'):
            result[''] = seconds(value)
    if not result:
        raise RuntimeError('Slurm returned no partition default times.')
    _DEFAULTS = time.monotonic(), result
    return result


def submission_override(script, minimum, *, env=None):
    """Return a bounded CLI override; do not change the hash-bound script."""
    validate_minimum(minimum)
    options = {}
    for line in Path(script).read_text().splitlines():
        if line.strip() and not line.lstrip().startswith('#'):
            break  # Slurm ignores directives after the first executable line.
        if not line.lstrip().startswith('#SBATCH'):
            continue
        words = shlex.split(line.lstrip()[len('#SBATCH'):])
        for index, word in enumerate(words):
            key, sep, value = word.partition('=')
            if key in ('--time', '-t', '--partition', '-p'):
                if not sep:
                    if index+1 == len(words):
                        raise ValueError(f'Missing value for {key} in {script}')
                    value = words[index+1]
                options['time' if key in ('--time', '-t') else 'partition'] = value
    # SBATCH environment options take precedence over the script itself.
    env = os.environ if env is None else env
    requested = env.get('SBATCH_TIMELIMIT', options.get('time'))
    partition = env.get('SBATCH_PARTITION', options.get('partition', ''))
    if requested is None:
        try:
            defaults = partition_defaults()
            original = max(defaults[name] for name in partition.split(','))
        except KeyError as error:
            raise RuntimeError(f'Slurm omitted a default time for partition {partition!r}') from error
    else:
        original = seconds(requested)
    if requested is not None and original == 0:  # A requested zero means unlimited.
        return {}
    if original >= seconds(minimum):
        return {}
    limits = partition_limits(cached=True)
    try:
        cap = min(limits[name] for name in partition.split(','))
    except KeyError as error:
        raise RuntimeError(f'Slurm did not report the selected partition: {partition!r}') from error
    selected = max(original, min(seconds(minimum), cap))
    if selected <= original:
        return {}
    return {'requested_minimum': minimum, 'original': requested,
            'effective': duration(selected), 'partition': partition,
            'capped_to_partition': cap < seconds(minimum)}


def jobs(root):
    root = Path(root).resolve(strict=True)
    owner = pwd.getpwuid(os.getuid()).pw_name
    output = command(['squeue', '-h', '-u', owner, '-o', '%i|%T|%M|%l|%P|%Z']).stdout
    result = []
    for line in output.splitlines():
        if not line.strip():
            continue
        ident, state, elapsed, limit, partition, directory = [x.strip() for x in line.split('|', 5)]
        if directory in ('', '(null)', 'N/A'):
            raise ValueError('Slurm omitted a work directory; cannot safely scope jobs to this run.')
        if not Path(directory).is_absolute() or not Path(directory).resolve().is_relative_to(root):
            continue
        # Never let an array master update affect tasks outside this scope.
        if not re.fullmatch(r'\d+(?:_\d+)?', ident):
            raise ValueError(f'Unsupported aggregate/heterogeneous job ID: {ident}')
        result.append(dict(id=ident, state=state, elapsed=elapsed, limit=limit,
                           partition=partition, directory=directory,
                           remaining_seconds=(None if limit == 'NOT_SET' or elapsed == 'INVALID'
                                              else seconds(limit)-seconds(elapsed))))
    return result


def extend(root, minimum, *, apply=False):
    validate_minimum(minimum)
    limits = partition_limits()
    before = jobs(root)
    selected, skipped = [], []
    for job in before:
        if job['state'] not in ('RUNNING', 'PENDING', 'CONFIGURING', 'SUSPENDED'):
            continue
        if job['limit'] == 'NOT_SET':
            skipped.append({**job, 'reason': 'Slurm has not established a limit yet'})
            continue
        if seconds(job['limit']) >= seconds(minimum):
            continue
        cap = min(limits[p] for p in job['partition'].split(','))
        if seconds(minimum) > cap:
            skipped.append({**job, 'reason': 'requested limit exceeds partition maximum'})
        else:
            selected.append(job)
    # One bulk update; no cancel/requeue, no repeat after an uncertain request.
    args = ['scontrol', 'update', 'JobId='+','.join(j['id'] for j in selected),
            'TimeLimit='+minimum] if selected else []
    report = dict(before=before, selected=selected, skipped=skipped, command=args,
                  applied=False, verified=[], unresolved=[], request_error=None)
    if not apply or not args:
        return report
    # Recheck immediately before the request: do not lower a limit that another
    # operator already raised, or act on a job that left this run/the queue.
    latest = {j['id']: j for j in jobs(root)}
    selected = [j for j in selected if j['id'] in latest
                and latest[j['id']]['state'] in ('RUNNING', 'PENDING', 'CONFIGURING', 'SUSPENDED')
                and latest[j['id']]['limit'] != 'NOT_SET'
                and seconds(latest[j['id']]['limit']) < seconds(minimum)
                and latest[j['id']]['partition'] == j['partition']]
    report['selected'] = selected
    if not selected:
        report['command'] = []
        return report
    args[2] = 'JobId='+','.join(j['id'] for j in selected)
    try:
        command(args)
    except (OSError, subprocess.SubprocessError) as error:
        report['request_error'] = str(getattr(error, 'stderr', '') or error)
    report['applied'] = True
    try:
        after = {j['id']: j for j in jobs(root)}
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        report['verification_error'] = str(error)
        report['unresolved'] = [{**j, 'reason': 'scheduler verification unavailable'} for j in selected]
        return report
    for job in selected:
        current = after.get(job['id'])
        if current and current['limit'] != 'NOT_SET' and seconds(current['limit']) >= seconds(minimum):
            report['verified'].append(current)
        else:
            report['unresolved'].append({**job, 'after': current,
                'reason': 'not increased' if current else 'left queue; extension unverified'})
    return report


def tail(path, limit=32768):
    try:
        with path.open('rb') as stream:
            stream.seek(max(0, path.stat().st_size-limit))
            return stream.read(limit).decode(errors='replace')
    except FileNotFoundError:
        return ''


def evidence(job):
    directory = Path(job['directory'])
    files = [p for p in directory.glob('*') if p.suffix in ('.log', '.out')
             and p.name != 'optimization.log']
    native = max(files, key=lambda p: p.stat().st_mtime) if files else None
    text = tail(native) if native else ''
    patterns = r'Iteration\s+\d+\s+CC energy|Norm of residual|Wall time \[min\]|Normal termination|Convergence criterion'
    native_rows = [s.strip() for s in text.splitlines() if re.search(patterns, s)]
    optimizer = [s.strip() for s in tail(directory/'optimization.log').splitlines()
                 if re.match(r'\s*(Sella|BFGS)\s+\d+', s)]
    curvature = sorted((p for p in directory.glob('minimum_curvature_*.json')
                        if p.stem.rsplit('_', 1)[1].isdigit()),
                       key=lambda p: int(p.stem.rsplit('_', 1)[1]))
    probes = None
    if curvature:
        data = json.loads(curvature[-1].read_text())
        probes = dict(status=data.get('status'), completed=len(data.get('probes', [])),
                      total=6*len(data.get('geometry_angstrom', [])),
                      attempt=curvature[-1].name)
    return dict(native_file=native.name if native else None,
                native_hash=hashlib.sha256(text.encode()).hexdigest(),
                native_rows=native_rows[-9:] or text.splitlines()[-3:],
                optimizer_rows=optimizer[-3:], curvature=probes)


def counters(identifiers):
    if not identifiers:
        return {}, None
    try:
        output = command(['sstat', '-a', '-j', ','.join(identifiers), '-n', '-P',
                          '--format=JobID,AveCPU,MaxRSS,AveDiskRead,AveDiskWrite']).stdout
        rows = {}
        for line in output.splitlines():
            values = line.rstrip('|').split('|')
            if len(values) != 5:
                raise ValueError('Unrecognized sstat output')
            rows[values[0]] = values[1:]
        return rows, None
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        return {}, str(getattr(error, 'stderr', '') or error)


def progress(root, sample_seconds=20):
    initial = jobs(root)
    running = [j['id'] for j in initial if j['state'] == 'RUNNING']
    errors, first = [], {}
    for job in initial:
        try:
            first[job['id']] = evidence(job)
        except (OSError, ValueError) as error:
            errors.append(f"{job['id']}: {error}")
    a, error_a = counters(running)
    time.sleep(sample_seconds)
    b, error_b = counters(running)
    current = jobs(root)
    for job in current:
        try:
            job['evidence'] = evidence(job)
            old = first.get(job['id'])
            job['native_changed'] = (old['native_hash'] != job['evidence']['native_hash']) if old else None
        except (OSError, ValueError) as error:
            errors.append(f"{job['id']}: {error}")
        steps = sorted(k for k in a if k in b and k.split('.', 1)[0] == job['id'])
        job['cpu_changed'] = any(a[k][0] != b[k][0] for k in steps) if steps else None
        job['io_changed'] = any(a[k][2:] != b[k][2:] for k in steps) if steps else None
    return dict(captured_utc=datetime.now(timezone.utc).isoformat(), jobs=current,
                sample_seconds=sample_seconds, counters_before=a, counters_after=b,
                scheduler_errors=[e for e in (error_a, error_b) if e], read_errors=errors)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--minimum', default='7-00:00:00')
    parser.add_argument('--apply', action='store_true')
    parser.add_argument('--progress', action='store_true', help='Read-only 20-second activity sample.')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    root = args.root.resolve(strict=True)
    if args.progress and args.apply:
        parser.error('--progress is read-only; do not combine it with --apply.')
    if args.output and (args.output.resolve().is_relative_to(root) or args.output.exists()):
        parser.error('--output must be a new file outside the run.')
    if args.progress:
        print('Sampling actual output content and CPU/I/O counters for 20 seconds...', flush=True)
        report = progress(root)
        print(report['captured_utc'], dict(Counter(j['state'] for j in report['jobs'])))
        for job in report['jobs']:
            directory = Path(job['directory']).resolve()
            case = directory.relative_to(root).parts[0] if directory != root else '.'
            graph = directory.parent.parent if directory.parent.name == 'tasks' else directory.parent
            if graph.name == 'composite':
                graph = graph.parent
            label = f'{case}/{graph.name}/{directory.name}'
            print(f"J{job['id']} {job['state']} elapsed={job['elapsed']} "
                  f"remaining={duration(job['remaining_seconds'])} "
                  f"CPU_changed={job.get('cpu_changed')} IO_changed={job.get('io_changed')} "
                  f"output_changed={job.get('native_changed')} {label}")
            ev = job.get('evidence', {})
            detail = (str(ev['curvature']) if ev.get('curvature') else
                      ' | '.join(ev.get('optimizer_rows', [])[-1:] or ev.get('native_rows', [])[-3:]))
            if detail:
                print('  '+detail[:330])
        print('Read errors:', report['read_errors'])
        print('sstat errors:', report['scheduler_errors'])
        print('Activity is not an ETA; unchanged output during 20 seconds does not prove a stall.')
    else:
        report = extend(root, args.minimum, apply=args.apply)
        print('Selected:', len(report['selected']), 'Partition-limited:', len(report['skipped']))
        if report['command']:
            print('Command:', shlex.join(report['command']))
        if report['applied']:
            print('Verified increases:', ', '.join(j['id'] for j in report['verified']) or 'none')
            print('Unverified:', ', '.join(j['id'] for j in report['unresolved']) or 'none')
        if report['request_error']:
            print('Slurm response:', report['request_error'].strip())
            print('Running-job increases normally require a Slurm administrator. No jobs were cancelled or requeued.')
        if report.get('verification_error'):
            print('Verification failed:', report['verification_error'])
        for job in report['skipped']:
            print(f"J{job['id']}: {job['reason']} ({job['partition']}).")
    if args.output:
        with args.output.open('x') as stream:
            json.dump(report, stream, indent=2)
            stream.write('\n')
        print('Full report:', args.output)
    return int(bool(report.get('unresolved') or report.get('skipped') or
                    report.get('read_errors') or report.get('scheduler_errors')))


if __name__ == '__main__':
    raise SystemExit(main())
