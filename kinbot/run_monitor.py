"""Read-only PES dashboard and log follower, independent of Slurm accounting.

Usage: python -m kinbot.run_monitor RUN [--all] [--follow] [--json]
RUN can be one PES or a handoff directory containing prepared.json.
No controller, dispatcher, database, or restart code is executed by this view.
"""
import argparse
from collections import Counter
from datetime import datetime
import getpass
import json
import os
from pathlib import Path
import re
import subprocess
import time

TAIL_BYTES = 65536
EVENT = re.compile(r'^\d{2}-\w{3}-\d{2} \d\d:\d\d:\d\d-(?:INFO|DEBUG|WARNING|ERROR|CRITICAL):')
ERROR = re.compile(r'-ERROR:|-CRITICAL:|Traceback \(most recent|QPErr|'
                   r'Error termination|segmentation (?:violation|fault)|'
                   r'GLOBAL ERROR|\?\s*Error|\b\w*(?:Error|Exception):|'
                   r'Permission denied|No space left|oom.kill|OUT_OF_MEMORY|'
                   r'Worker failure:|Frequency validation failed|'
                   r'Converged but frequency check failed', re.I)
PROGRESS = re.compile(r'PROGRAM\s*\*|SCF Done:|Step number|Optimization completed|'
                      r'Normal termination|Maximum Force|RMS Force|Maximum Displacement|'
                      r'RMS Displacement|^\s*(?:Sella|BFGS|LBFGS|FIRE)(?::|\s+\d+\s)|'
                      r'Frequencies\s+--|Zero-point correction|Minimum recovery \d|'
                      r'^\s*\d+\s+\d+\s+\d+\s+[-+]?\d+\.\d+|'
                      r'^\s*\?\s*\S|convergence|\bKBGAP\w*\s*=', re.I)
LOG_NAMES = {'driver.log', 'pes.log', 'kinbot.log', 'kinbot.err', 'kinbot.out'}
STATE_NAMES = {'state.json', 'final_report.json', 'execution.json',
               'mess_execution.json', '.kinbot_post_reaction_branches.json',
               'kinbot_status.json', 'pes_status.json', 'recovery.json',
               'minimum_recovery.json'}
REACTION_STAGES = {-999: 'failed/rejected', -1: 'accepted', 0: 'TS search',
                   1: 'barrier/IRC', 2: 'product characterization',
                   3: 'start product/TS optimization', 4: 'product/TS optimization',
                   5: 'final validation'}


def inside(path, root):
    return Path(path).resolve().is_relative_to(root)


def tail(path, limit=TAIL_BYTES):
    with path.open('rb') as handle:
        size = os.fstat(handle.fileno()).st_size
        handle.seek(max(0, size-limit))
        data = handle.read(limit)
    if size > limit:
        data = data.partition(b'\n')[2]
    return data.decode(errors='replace')


def read_json(path, warnings):
    try:
        if path.stat().st_size > 8 * 1024 * 1024:
            raise ValueError('JSON exceeds monitor read limit (8 MiB)')
        value = json.loads(path.read_text())
        if not isinstance(value, dict):
            raise ValueError('expected a JSON object')
        return value
    except FileNotFoundError:
        return {}
    except (OSError, ValueError) as error:
        warnings.append(f'{path}: {error}')
        return {}


def queue():
    fields = ('job_id', 'state', 'nodes', 'elapsed', 'limit', 'cpus',
              'memory', 'excluded_nodes', 'name', 'reason', 'directory')
    try:
        result = subprocess.run(['squeue', '-h', '-u', getpass.getuser(), '-o',
            '%i|%T|%N|%M|%l|%C|%m|%x|%j|%R|%Z'],
            capture_output=True, text=True, timeout=15, check=False)
        if result.returncode:
            raise RuntimeError(result.stderr.strip() or f'squeue exit {result.returncode}')
        rows = []
        for line in result.stdout.splitlines():
            if not line.strip():
                continue
            parts = [part.strip() for part in line.split('|')]
            if len(parts) != len(fields):
                raise ValueError(f'Unrecognized squeue row: {line!r}')
            rows.append(dict(zip(fields, parts)))
        return rows, None
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
        return [], str(error)


def inventory(root, warnings):
    """Bound metadata traversal; exclude archives, scratch and sampling payloads."""
    paths, visited = [], 0
    def walk_error(error):
        warnings.append(str(error))
    for directory, dirs, files in os.walk(root, followlinks=False, onerror=walk_error):
        if Path(directory).name.endswith('_minimum_recovery'):
            # Retained rejected Hessians and trial logs are historical evidence,
            # not new failures of the current optimization. Show their reports.
            paths.extend(sorted(Path(directory).glob('attempt_*/recovery.json')))
            dirs[:] = []
            continue
        dirs[:] = sorted(d for d in dirs if not (
            d.startswith('.') or '.stopped.' in d or d.startswith('restart_')
            or d in {'assets', 'orbitals', 'scratch', '__pycache__', 'archive', 'archives'}
            or Path(directory, d).is_symlink()))
        visited += len(files) + 1
        if visited > 30000:
            warnings.append('Monitor inventory capped at 30000 entries; some files omitted.')
            break
        for name in files:
            path = Path(directory, name)
            if (name.endswith(('.execution.json', '.rotdpy.json', '.result.json')) or
                    name in LOG_NAMES | STATE_NAMES | {'kinbot_monitor.out', 'pes_summary.txt',
                'kinbot_optimization_monitor.json'} or
                    path.suffix in {'.log', '.err', '.stdout', '.stderr', '.exitcode', '.out'}):
                if not path.is_symlink():
                    paths.append(path)
        # RotdPy sample trees are very large. Report surface output counts
        # separately, without traversing individual sample/checkpoint files.
        if Path(directory).name.startswith('kb_'):
            paths.append(Path(directory))
            dirs[:] = []
    return paths


def log_record(path, root, now):
    raw = tail(path)
    lines = raw.splitlines()
    events = [line.strip() for line in lines if EVENT.match(line) and
              line.split(': ', 1)[-1].strip() and
              not re.search(r'INFO: (?:KinBot \d|RDKit |Input parameters)', line)]
    errors = []
    for index, line in enumerate(lines):
        if ERROR.search(line):
            errors.append(line.strip())
            if 'Converged but frequency check failed' in line and index + 1 < len(lines):
                errors.append(lines[index + 1].strip())
    progress = [line.strip() for line in lines if PROGRESS.search(line)]
    return {'path': str(path.relative_to(root)), 'mtime': path.stat().st_mtime,
            'age_seconds': max(0, now-path.stat().st_mtime),
            'bytes': path.stat().st_size, 'events': events[-8:],
            'errors': errors[-6:], 'progress': progress[-6:]}


def pes_summary(path, warnings):
    result = {'counts': {}, 'workers': {}}
    try:
        section = None
        for line in tail(path).splitlines():
            match = re.fullmatch(r'(Total|Running|Finished)\s+(\d+)', line.strip())
            if match:
                result['counts'][match[1]] = int(match[2])
            elif line.strip() in {'Running:', 'Finished:', 'Waiting:', 'Failed:'}:
                section = line.strip()[:-1]
            elif section and line.strip():
                result['workers'][line.strip()] = section
    except FileNotFoundError:
        pass
    except OSError as error:
        warnings.append(str(error))
    return result


def reaction_summary(path, warnings):
    rows = []
    try:
        for line in tail(path).splitlines():
            parts = line.split()
            if len(parts) < 3:
                continue
            state = int(parts[0])
            rows.append({'name': parts[2], 'step': parts[1], 'code': state,
                         'stage': REACTION_STAGES.get(state, f'state {state}')})
    except (OSError, ValueError) as error:
        warnings.append(f'{path}: {error}')
    return {'path': str(path), 'counts': dict(Counter(r['stage'] for r in rows)), 'rows': rows}


def artifacts(paths, root, warnings):
    results = []
    for path in paths:
        if path.name not in STATE_NAMES and not path.name.endswith(('.execution.json', '.rotdpy.json', '.result.json')):
            continue
        data = read_json(path, warnings)
        if not data:
            continue
        item = {'path': str(path.relative_to(root)), 'status': data.get('status', 'recorded')}
        if path.name == 'state.json' and isinstance(data.get('tasks'), dict):
            spec = read_json(path.with_name('workflow.json'), warnings)
            tasks = {t['id']: t for t in spec.get('tasks', [])}
            tasks.update({k: tasks.get(k, {}) for k in data['tasks']})
            rows = []
            for name, task in tasks.items():
                entry = data['tasks'].get(name, {})
                deps = entry.get('blocked_by', task.get('depends_on', []))
                blockers = {dep: data['tasks'].get(dep, {}).get('status', 'waiting')
                            for dep in deps if data['tasks'].get(dep, {}).get('status') != 'complete'}
                rows.append({'name': name, 'status': entry.get('status', 'waiting'),
                             'job_id': entry.get('job_id'), 'blocked_by': blockers,
                             'error': entry.get('error') or entry.get('skip_reason')})
            item.update(counts=dict(Counter(row['status'] for row in rows)), tasks=rows)
        elif path.name == '.kinbot_post_reaction_branches.json':
            item.update(branches=data.get('branches'), mess_gate=data.get('mess_gate'),
                        final_solver=data.get('final_solver'))
        elif path.name == 'recovery.json':
            item['recovery'] = {key: data.get(key) for key in
                ('attempt', 'frequency_cm-1', 'displacement_angstrom',
                 'selected', 'energy_change_eV')}
        elif path.name == 'minimum_recovery.json':
            item['recovery'] = {'controls': data.get('controls'),
                                'curvature_checks': data.get('attempts', [])}
        elif path.name.endswith('.result.json'):
            item['recovery'] = {'completed_attempts': data.get('minimum_recovery_attempts', 0)}
        elif path.name == 'execution.json' and data.get('status') not in ('failed', 'unavailable'):
            continue  # Per-task successful execution is covered by state.json.
        if data.get('error'):
            item['error'] = str(data['error'])[-1200:]
        results.append(item)
    return results


def case_snapshot(root, jobs, processes, now, *, log_limit=40):
    warnings = []
    files = inventory(root, warnings)
    exit_records = []
    failed_log_names = set()
    for path in files:
        if path.suffix != '.exitcode':
            continue
        try:
            values = tail(path, 1024).split()
            if len(values) not in (1, 2):
                raise ValueError('expected exit code or job-id/exit-code')
            code = int(values[-1])
            exit_records.append({'path': str(path.relative_to(root)), 'returncode': code,
                                 'job_id': values[0] if len(values) == 2 else None})
            if code:
                failed_log_names.update(path.stem + suffix for suffix in
                    ('.log', '_sella.log', '.err', '.stdout', '.stderr', '.out'))
        except (ValueError, OSError) as error:
            warnings.append(f'{path}: {error}')
    par = read_json(root/'input.json', warnings)
    summary = pes_summary(root/'pes_summary.txt', warnings)
    logs = [p for p in files if p.name in LOG_NAMES or p.suffix in {'.log', '.err', '.stderr', '.stdout'}
            or (p.suffix == '.out' and p.name != 'kinbot_monitor.out')]
    mtimes = {}
    for path in logs:
        try:
            mtimes[path] = path.stat().st_mtime
        except OSError as error:
            warnings.append(str(error))
    # Failed jobs' native/worker diagnostics must not be crowded out by many
    # newer successful calculations. Controller logs retain first priority.
    chosen = sorted(mtimes, key=lambda p: (p.name in LOG_NAMES,
                    p.name in failed_log_names, mtimes[p]), reverse=True)
    records = []
    for path in chosen[:log_limit]:
        try:
            records.append(log_record(path, root, now))
        except OSError as error:
            warnings.append(f'{path}: {error}')
    if len(chosen) > log_limit:
        warnings.append(f'Reading {log_limit}/{len(chosen)} log tails; use --all for up to 200.')
    record_by_path = {r['path']: r for r in records}
    workers = []
    worker_dirs = {p.parent for p in files if p.name in
                   {'kinbot.log', 'kinbot_monitor.out', 'kinbot_status.json'}}
    for directory in sorted(worker_dirs):
        name = str(directory.relative_to(root))
        status = read_json(directory/'kinbot_status.json', warnings)
        log = record_by_path.get(str((directory/'kinbot.log').relative_to(root)), {})
        rows = reaction_summary(directory/'kinbot_monitor.out', warnings) if (directory/'kinbot_monitor.out').exists() else {}
        optimization = read_json(directory/'kinbot_optimization_monitor.json', warnings)
        workers.append({'directory': name, 'reported_pes_state': summary['workers'].get(name),
                        'status': status, 'events': log.get('events', []),
                        'reactions': rows, 'optimizations': optimization.get('reactions', [])})
    sampling = []
    for path in files:
        if path.name.startswith('kb_') and path.is_dir():
            sampling.append({'directory': str(path.relative_to(root)),
                'surface_directories': len(list(path.glob('Surface_*'))),
                'surface_output_files': len(list(path.glob('output/surface_*.dat')))})
    return {'directory': str(root), 'name': root.name, 'sampling': sampling,
            'allocation_cap': par.get('queue_job_limit'), 'configured_exclusions': par.get('slurm_exclude_nodes', []),
            'jobs': [job for job in jobs if inside(job['directory'], root)],
            'processes': [p for p in processes if inside(p['cwd'], root)],
            'pes': summary, 'workers': workers, 'logs': records,
            'exit_records': exit_records, 'artifacts': artifacts(files, root, warnings),
            'warnings': warnings}


def snapshot(root, *, include_all=False):
    root = Path(root).expanduser().resolve(strict=True)
    warnings = []
    prepared = read_json(root/'prepared.json', warnings)
    directories = [Path(p).parent.resolve() for p in prepared.get('inputs', [])] or [root]
    if any(not inside(d, root) for d in directories):
        raise ValueError('Prepared input directory is outside the monitored tree.')
    jobs, queue_error = queue()
    process_error = None
    try:
        from kinbot.run_cleanup import _processes
        processes = _processes(root)
    except (OSError, RuntimeError) as error:
        processes, process_error = [], str(error)
    now = time.time()
    return {'time': datetime.fromtimestamp(now).astimezone().isoformat(timespec='seconds'),
            'directory': str(root), 'scheduler_error': queue_error, 'process_error': process_error,
            'warnings': warnings, 'cases': [case_snapshot(d, jobs, processes, now,
                log_limit=200 if include_all else 40) for d in directories]}


def age(seconds):
    if seconds < 60:
        return f'{seconds:.0f}s'
    if seconds < 3600:
        return f'{seconds/60:.1f}m'
    return f'{seconds/3600:.1f}h'


def render(report, *, include_all=False, lines=4):
    out = [report['time'], report['directory']]
    for key in ('scheduler_error', 'process_error'):
        if report[key]:
            out.append(f'{key.upper()}: UNKNOWN — {report[key]}')
    for case in report['cases']:
        out += ['', f'=== {case["name"]} ===']
        out.append(f'Slurm: {"UNKNOWN" if report["scheduler_error"] else len(case["jobs"])} allocations'
                   f' / cap {case["allocation_cap"] or "unspecified"}; '
                   f'local processes: {"UNKNOWN" if report["process_error"] else len(case["processes"])}')
        out.append('Configured exclusions: '+(','.join(case['configured_exclusions']) or '(none)'))
        out.append('JOB       STATE       NODES       ELAPSED / LIMIT   CPU  MEMORY  NAME')
        for job in case['jobs']:
            out.append(f'{job["job_id"]:<9} {job["state"]:<11} {job["nodes"]:<11} '
                       f'{job["elapsed"]} / {job["limit"]}  {job["cpus"]}  {job["memory"]}  {job["name"]}')
            out.append(f'  WorkDir: {job["directory"]}; Slurm exclusions: {job["excluded_nodes"] or "(none)"}')
            if job['state'] != 'RUNNING':
                out.append('  Reason: '+job['reason'])
        for process in case['processes']:
            out.append(f'PID {process["pid"]} [{process["state"]}] {process["command"]}')
        if case['pes']['counts']:
            out.append('PES worker counts: '+str(case['pes']['counts'])+' (legacy Finished means exited, not successful)')
        failures = [r for r in case['exit_records'] if r['returncode']]
        for failure in failures:
            out.append(f'FAILED QC EXIT: {failure["path"]} job={failure["job_id"] or "unknown"} exit={failure["returncode"]}')
        for worker in case['workers']:
            status = worker['status'].get('status', 'no completion record')
            out.append(f'WORKER {worker["directory"]}: {status}; PES says {worker["reported_pes_state"] or "unlisted"}')
            if worker['reactions']:
                out.append('  Reactions: '+str(worker['reactions']['counts']))
                rows = worker['reactions']['rows']
                for row in rows if include_all else [r for r in rows if r['code'] != -1][:12]:
                    out.append(f'    {row["name"]}: {row["stage"]}; step={row["step"]}')
            for row in worker['optimizations']:
                out.append(f'  Optimizing {row["reaction"]}: pending={row.get("pending")} failures={row.get("failures")}')
                if include_all:
                    out.append('    '+json.dumps(row.get('optimizers', [])))
        for item in case['artifacts']:
            out.append(f'RECORD {item["path"]}: {item["status"]}')
            for key in ('counts', 'branches', 'mess_gate', 'final_solver', 'recovery', 'error'):
                if item.get(key) is not None:
                    out.append(f'  {key}: {item[key]}')
            tasks = item.get('tasks', [])
            for task in tasks if include_all else [t for t in tasks if t['status'] in ('failed', 'blocked', 'submitted', 'unavailable')][:12]:
                out.append(f'  {task["name"]}: {task["status"]} job={task["job_id"]}; waiting for={task["blocked_by"]}; {task["error"] or ""}')
        for sampling in case['sampling']:
            out.append(f'ROTD {sampling["directory"]}: {sampling["surface_output_files"]}/'
                       f'{sampling["surface_directories"]} surface output files; '
                       'file counts only, final flux validation is separate.')
        out.append('RECENT LOG ACTIVITY (bounded tails; warnings/errors may belong to earlier attempts):')
        displayed = 0
        for log in sorted(case['logs'], key=lambda r: (bool(r['errors']), r['mtime']), reverse=True):
            content = log['events'][-lines:] + log['progress'][-lines:] + log['errors']
            content = list(dict.fromkeys(content))
            if not content:
                continue
            if displayed >= (200 if include_all else 12):
                out.append('  More logs omitted; use --all.')
                break
            displayed += 1
            out.append(f'  {log["path"]} [modified {age(log["age_seconds"])} ago; {log["bytes"]} bytes]')
            out.extend('    '+line for line in content)
        out.extend('READ WARNING: '+warning for warning in case['warnings'])
        if not case['jobs'] and not case['processes'] and not report['scheduler_error'] and not report['process_error']:
            out.append('INACTIVE: no live jobs/controllers observed. This is not a completion verdict.')
    out.extend('READ WARNING: '+warning for warning in report['warnings'])
    out.append('\nRead-only. Record statuses are not a new scientific audit. Quiet output is not proof of a stall.')
    return '\n'.join(out)


class LogFollower:
    """Track append offsets; rediscover new logs and handle rotation/truncation."""
    def __init__(self):
        self.positions = {}

    def poll(self, root, initial_lines=4):
        warnings, output = [], []
        paths = inventory(root, warnings)
        for path in sorted(paths):
            if (path.name == 'kinbot_monitor.out' or
                    (path.name not in LOG_NAMES and path.suffix not in {'.log', '.err', '.stderr', '.stdout', '.out'})):
                continue
            try:
                with path.open('rb') as stream:
                    stat = os.fstat(stream.fileno())
                    identity = (stat.st_dev, stat.st_ino)
                    previous = self.positions.get(path)
                    fresh = previous is None or previous[0] != identity or stat.st_size < previous[1]
                    offset = max(0, stat.st_size-TAIL_BYTES) if fresh else previous[1]
                    stream.seek(offset)
                    raw = stream.read(TAIL_BYTES)
                    if fresh and offset:
                        raw = raw.partition(b'\n')[2]
                    complete, _, pending = raw.rpartition(b'\n')
                    if len(pending) == TAIL_BYTES:
                        # A single very long line must not pin the cursor forever.
                        warnings.append(f'{path}: line exceeds {TAIL_BYTES} bytes; skipping this chunk.')
                        pending = b''
                    self.positions[path] = (identity, stream.tell()-len(pending))
                    text = complete.decode(errors='replace').splitlines()
                    if fresh:
                        text = text[-initial_lines:] if initial_lines else []
                    if text:
                        output.append((str(path.relative_to(root)), text))
            except OSError as error:
                warnings.append(f'{path}: {error}')
        return output, warnings


def add_arguments(parser):
    parser.add_argument('--all', action='store_true', help='Show completed tasks, all reaction rows and more log tails.')
    parser.add_argument('--lines', type=int, default=4, help='Recent event/native progress lines per log (default 4).')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--json', action='store_true', help='Machine-readable snapshot.')
    mode.add_argument('--follow', action='store_true', help='Snapshot then follow new controller and QC log text.')
    parser.add_argument('--interval', type=float, default=3, help='Log-follow poll interval, seconds (minimum 1).')


def run(root, args):
    root = Path(root).expanduser().resolve(strict=True)
    if args.lines < 1 or not 1 <= args.interval <= 60:
        raise ValueError('--lines must be positive; --interval must be between 1 and 60 seconds.')
    report = snapshot(root, include_all=args.all)
    print(json.dumps(report, indent=2) if args.json else render(report, include_all=args.all, lines=args.lines), flush=True)
    if args.follow:
        follower = LogFollower()
        # Snapshot already showed existing activity; follow only newly written text.
        follower.poll(root, 0)
        print('\nFollowing logs; Ctrl-C stops only this viewer.', flush=True)
        try:
            while True:
                entries, warnings = follower.poll(root, args.lines)
                for name, text in entries:
                    print('\n==> '+name+' <==\n'+'\n'.join(text), flush=True)
                for warning in warnings:
                    print('READ WARNING: '+warning, flush=True)
                time.sleep(args.interval)
        except KeyboardInterrupt:
            pass


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    add_arguments(parser)
    args = parser.parse_args()
    run(args.root, args)


if __name__ == '__main__':
    main()
