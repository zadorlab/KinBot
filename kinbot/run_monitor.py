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
import shlex
import shutil
import subprocess
import textwrap
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
                      r'^\s*\?\s*\S|convergence|\b(?:KBGAP|RECONC)\w*\s*=|'
                      r'E\((?:U?B2PLYP|R?MP2)\)\s*=|!.*Energy|'
                      r'SCF failed|Convergence failure', re.I)
# Require a decimal/exponent after the integer iteration/state columns.
# Molpro's all-integer dataset record lists are not solver iterations.
NUMBER = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[DdEe][-+]?\d+)?'
REAL = r'[-+]?(?:\d+\.\d*|\.\d+|\d+[DdEe][-+]?\d+)(?:[DdEe][-+]?\d+)?'
ITERATION = re.compile(r'^\s*\d+\s+(?:\d+\s+\d+\s+)?' +
                       rf'{REAL}\s+(?:{NUMBER}\s+){{2}}{NUMBER}(?:\s|$)')
TABLE_HEADER = re.compile(r'^\s*(?:Step\s+Time\s+Energy|ITER\.?\s+.*(?:ENERGY|STATE)|'
                          r'ITERATION\s+.*ENERGY)', re.I)
LOG_NAMES = {'driver.log', 'pes.log', 'kinbot.log', 'kinbot.err', 'kinbot.out'}
STATE_NAMES = {'state.json', 'final_report.json', 'execution.json',
               'mess_execution.json', '.kinbot_post_reaction_branches.json',
               'kinbot_status.json', 'pes_status.json', 'recovery.json',
               'minimum_recovery.json', 'composite_preflight.json', 'vrc_model_review.json',
               'scan_review.json'}
INVENTORY_DIRECTORY_LIMIT = 30000
INVENTORY_LOG_LIMIT = 30000
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


def inventory(root, warnings, archives=None):
    """Bound metadata traversal; exclude archives, scratch and sampling payloads."""
    paths, visited, log_count = [], 0, 0
    def walk_error(error):
        warnings.append(str(error))
    for directory, dirs, files in os.walk(root, followlinks=False, onerror=walk_error):
        if 'attempts' in dirs and 'workflow.json' in files:
            # Dispatcher retries move obsolete inputs/results/logs here. Never
            # interpret them as current execution records, even in --details.
            dirs.remove('attempts')
            if archives is not None:
                archives.append(str((Path(directory)/'attempts').relative_to(root)))
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
        # Thousands of finite-difference step files in one geometry job must
        # not hide the state/review records of every subsequent well. Bound
        # traversal and optional log inventory separately from metadata.
        visited += 1
        if visited > INVENTORY_DIRECTORY_LIMIT:
            warnings.append('Monitor directory inventory capped; some records omitted.')
            break
        for name in files:
            path = Path(directory, name)
            metadata = (name.endswith(('.execution.json', '.rotdpy.json', '.result.json', '.scf_recovery.json', '.exitcode')) or
                    name in LOG_NAMES | STATE_NAMES | {'kinbot_monitor.out', 'pes_summary.txt',
                'kinbot_optimization_monitor.json'})
            log = path.suffix in {'.log', '.err', '.stdout', '.stderr', '.out'}
            if metadata or (log and log_count < INVENTORY_LOG_LIMIT):
                if not path.is_symlink():
                    paths.append(path)
            if log and not metadata:
                log_count += 1
                if log_count == INVENTORY_LOG_LIMIT + 1:
                    warnings.append('QC log inventory capped; state/review inventory continues.')
        # RotdPy sample trees are very large. Report surface output counts
        # separately, without traversing individual sample/checkpoint files.
        if Path(directory).name.startswith('kb_'):
            paths.append(Path(directory))
            dirs[:] = []
    return paths


def log_record(path, root, now, *, tail_lines=4):
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
    # KinBot parameter dumps are not QC progress. In particular "convergence"
    # can occur inside a many-kilobyte inline state model in kinbot.log.
    progress = [] if path.name in LOG_NAMES else [line.strip() for line in lines
                if PROGRESS.search(line) or ITERATION.match(line)]
    headers = [line.strip() for line in lines if TABLE_HEADER.match(line)]
    stat = path.stat()
    return {'path': str(path.relative_to(root)), 'mtime': stat.st_mtime,
            'age_seconds': max(0, now-stat.st_mtime),
            'bytes': stat.st_size, 'events': events[-max(8, tail_lines):],
            'errors': errors[-6:], 'progress': progress[-max(6, tail_lines):],
            'table_header': headers[-1] if headers else None,
            'raw_tail': lines[-tail_lines:]}


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
                         'stage': REACTION_STAGES.get(state, f'state {state}'),
                         'products': parts[3:]})
    except (OSError, ValueError) as error:
        warnings.append(f'{path}: {error}')
    return {'path': str(path), 'counts': dict(Counter(r['stage'] for r in rows)), 'rows': rows,
            'accepted_product_groups': len({tuple(sorted(r['products'])) for r in rows
                                            if r['code'] == -1 and r['products']})}


def previous_controller_record(path, root, warnings):
    """Identify retained branch/preflight reports from before a worker restart.

    Refresh preserves these reports until the replacement reaches the branch.
    Compare against the worker's start, never its status-file modification time
    (which also changes on exit). Unknown provenance keeps diagnostics visible.
    """
    starts = []
    for directory in path.parents:
        if not directory.is_relative_to(root):
            break
        # A replacement PES can be exploring other wells before it launches
        # this child. In that interval there is no new child status file.
        for name in ('kinbot_status.json', 'pes_status.json'):
            marker = directory/name
            try:
                if not marker.exists():
                    continue
                status = read_json(marker, warnings)
                started = datetime.fromisoformat(status['started_utc'])
                if started.tzinfo is not None:
                    starts.append(started.timestamp())
            except (KeyError, TypeError, ValueError):
                continue
            except OSError as error:
                warnings.append(f'{marker}: {error}')
    try:
        return bool(starts and path.stat().st_mtime < max(starts))
    except OSError as error:
        warnings.append(f'{path}: {error}')
        return False


def vrc_review_status(path, data, warnings):
    """Reconcile a retained well-level hold with its newer channel receipts.

    Older controllers clear the aggregate only after *all calculations* finish,
    although each channel request is marked model_accepted before dispatch.
    Never infer approval from a running job or an input refresh alone. Missing,
    unreadable, mismatched or older receipts leave the hold visible.
    """
    if data.get('status') != 'needs_review':
        return {}
    channels = data.get('channels')
    if not isinstance(channels, list) or not channels:
        return {}
    if any(isinstance(c, dict) and c.get('kind') == 'scan_connectivity' for c in channels):
        # A model-accepted orbital receipt never resolves a scan-path hold.
        return {}
    try:
        held_at = path.stat().st_mtime_ns
    except OSError as error:
        warnings.append(f'{path}: {error}')
        return {}
    receipts, outstanding = [], []
    for channel in channels:
        name = channel.get('reaction') if isinstance(channel, dict) else None
        if not isinstance(name, str):
            name = None
        record = {'reaction': name, 'status': 'unknown'}
        if isinstance(name, str) and name not in ('', '.', '..') and Path(name).name == name:
            request = path.parent/'vrctst/orbital_review'/name/'request.json'
            receipt = read_json(request, warnings)
            try:
                current = (receipt.get('reaction') == name
                           and request.stat().st_mtime_ns >= held_at)
            except FileNotFoundError:
                current = False
            except OSError as error:
                warnings.append(f'{request}: {error}')
                current = False
            record.update(status=receipt.get('status', 'unknown'),
                          current=current)
        if record['status'] != 'model_accepted' or not record.get('current'):
            outstanding.append(name or '(unidentified channel)')
        receipts.append(record)
    result = {'review_channels': receipts, 'recorded_status': data['status'],
              'recorded_error': data.get('error')}
    if outstanding:
        result['error'] = 'Orbital/state review still required or unconfirmed for: ' + ', '.join(outstanding)
    else:
        result.update(status='model_accepted', error=None,
                      stale_review_hold=True)
    return result


def artifacts(paths, root, warnings):
    results = []
    for path in paths:
        if path.name not in STATE_NAMES and not path.name.endswith(('.execution.json', '.rotdpy.json', '.result.json', '.scf_recovery.json')):
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
                error = entry.get('error') or entry.get('unavailable_reason') or entry.get('skip_reason')
                if not error and entry.get('status') in ('failed', 'unavailable'):
                    outcome = path.parent/'tasks'/name/'execution.json'
                    if outcome.is_file():
                        error = read_json(outcome, warnings).get('error')
                deps = entry.get('blocked_by', task.get('depends_on', []))
                blockers = {dep: data['tasks'].get(dep, {}).get('status', 'waiting')
                            for dep in deps if data['tasks'].get(dep, {}).get('status') != 'complete'}
                rows.append({'name': name, 'status': entry.get('status', 'waiting'),
                             'job_id': entry.get('job_id'), 'blocked_by': blockers,
                             'error': error,
                             'recovery': task.get('recovery')})
            item.update(counts=dict(Counter(row['status'] for row in rows)), tasks=rows)
        elif path.name == '.kinbot_post_reaction_branches.json':
            item.update(branches=data.get('branches'), mess_gate=data.get('mess_gate'),
                        final_solver=data.get('final_solver'))
            item['previous_controller'] = previous_controller_record(path, root, warnings)
        elif path.name == 'composite_preflight.json':
            item['previous_controller'] = previous_controller_record(path, root, warnings)
            item['preflight'] = {key: data.get(key) for key in (
                'controller_pid', 'phase', 'graph', 'task', 'graphs_complete',
                'graphs_total', 'environment_checks', 'environment_checks_reused',
                'slurm_scripts_tested', 'elapsed_seconds', 'updated_at', 'detail')}
        elif path.name == 'recovery.json':
            item['recovery'] = {key: data.get(key) for key in
                ('attempt', 'frequency_cm-1', 'displacement_angstrom',
                 'selected', 'energy_change_eV')}
        elif path.name == 'minimum_recovery.json':
            item['recovery'] = {'controls': data.get('controls'),
                                'curvature_checks': data.get('attempts', [])}
        elif path.name.endswith('.result.json'):
            item['recovery'] = {'completed_attempts': data.get('minimum_recovery_attempts', 0)}
        elif path.name.endswith('.scf_recovery.json'):
            item['recovery'] = {k: data.get(k) for k in ('reason', 'original_scf', 'retry_scf')}
        elif path.name == 'execution.json' and data.get('status') not in ('failed', 'unavailable'):
            continue  # Per-task successful execution is covered by state.json.
        if data.get('error'):
            item['error'] = str(data['error'])
        if path.name == 'vrc_model_review.json':
            item.update(vrc_review_status(path, data, warnings))
        results.append(item)
    return results


def case_snapshot(root, jobs, processes, now, *, log_limit=40, log_lines=4):
    warnings = []
    archives = []
    files = inventory(root, warnings, archives)
    case_jobs = [job for job in jobs if inside(job['directory'], root)]
    live_dirs = {Path(job['directory']).resolve() for job in case_jobs}
    # A targeted basis migration retains old scans in place for comparison.
    # Exclude only the journaled, unchanged, non-live graphs from the current
    # dashboard; their native files remain available for direct inspection.
    retired = []
    from kinbot.execution import file_sha256
    for marker in (root/'vrc_restart.json', root.parent/'vrc_restart.json'):
        for entry in read_json(marker, warnings).get('retired_scans', []):
            graph = Path(entry.get('graph', '/')).resolve()
            if (not graph.is_relative_to(root.resolve())
                    or 'vrctst/molpro_scan' not in graph.as_posix()
                    or any(d.is_relative_to(graph) for d in live_dirs)):
                continue
            try:
                if file_sha256(graph/'workflow.json') == entry.get('workflow_sha256'):
                    retired.append(graph)
            except OSError as error:
                warnings.append(str(error))
    if retired:
        files = [p for p in files if not any(p.is_relative_to(g) for g in retired)]
    task_records = artifacts(files, root, warnings)
    failed_dirs = {root/Path(item['path']).parent/'tasks'/task['name']
                   for item in task_records for task in item.get('tasks', [])
                   if task['status'] in ('failed', 'unavailable', 'blocked')}
    completed_workers = {Path(item['path']).parent for item in task_records
                         if Path(item['path']).name == 'kinbot_status.json'
                         and item['status'] == 'complete'
                         and len(Path(item['path']).parts) > 2}
    exit_records = []
    failed_log_names = set()
    archived_exits = read_json(root/'resume_provenance.json', warnings).get(
        'archived_failure_exits', {})
    for path in files:
        if path.suffix != '.exitcode':
            continue
        try:
            values = tail(path, 1024).split()
            if len(values) not in (1, 2):
                raise ValueError('expected exit code or job-id/exit-code')
            code = int(values[-1])
            from kinbot.execution import file_sha256
            relative = str(path.relative_to(root))
            historical = (relative in archived_exits
                          and archived_exits[relative] == file_sha256(path))
            exit_records.append({'path': str(path.relative_to(root)), 'returncode': code,
                                 'job_id': values[0] if len(values) == 2 else None,
                                 'archived_failure': historical})
            if code and not historical:
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
    chosen = sorted(mtimes, key=lambda p: (
        p.name in LOG_NAMES and p.parent.relative_to(root) not in completed_workers,
        p.name in failed_log_names or p.parent in failed_dirs,
        p.parent in live_dirs, mtimes[p]), reverse=True)
    # Read one tail per task before filling the remaining budget with older
    # force steps. A busy geometry must not crowd out a slower live calculation.
    first, rest, seen = [], [], set()
    for path in chosen:
        group = log_group(str(path.relative_to(root)))
        (rest if group in seen else first).append(path)
        seen.add(group)
    chosen = first + rest
    records = []
    for path in chosen[:log_limit]:
        try:
            record = log_record(path, root, now, tail_lines=log_lines)
            record['completed_worker'] = (path.name in LOG_NAMES and
                                          path.parent.relative_to(root) in completed_workers)
            record['live_task'] = path.parent in live_dirs
            record['historical_events'] = []
            if path.name in LOG_NAMES:
                for line in set(record['errors'] + record['events']):
                    match = re.search(r'-ERROR: (\S+) exited with worker status (\d+)\.$', line)
                    if match:
                        marker = str((path.parent/(match[1]+'.exitcode')).relative_to(root))
                        if any(e['path'] == marker and e['archived_failure'] and
                               e['returncode'] == int(match[2]) for e in exit_records):
                            record['historical_events'].append(line)
            records.append(record)
        except OSError as error:
            warnings.append(f'{path}: {error}')
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
            'exploration': read_json(root/'pes_exploration.json', warnings),
            'anl_selection': read_json(root/'pes_l2_selection.json', warnings),
            'allocation_cap': par.get('queue_job_limit'), 'configured_exclusions': par.get('slurm_exclude_nodes', []),
            'jobs': case_jobs,
            'processes': [p for p in processes if inside(p['cwd'], root)],
            'pes': summary, 'workers': workers, 'logs': records,
            'exit_records': exit_records, 'artifacts': task_records,
            'archived_attempt_directories': archives,
            'retired_scan_directories': [str(p) for p in retired],
            'log_coverage': {'read': len(records), 'available': len(chosen)},
            'warnings': warnings}


def snapshot(root, *, include_all=False, log_lines=4):
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
    from kinbot import run_policy
    policy = run_policy.read(root)
    policy_root = run_policy.find(root)
    from kinbot.live_control import summary as control_summary
    try:
        controls = control_summary(root)
    except (OSError, ValueError) as error:
        controls = {'error': str(error)}
    return {'time': datetime.fromtimestamp(now).astimezone().isoformat(timespec='seconds'),
            'timestamp': now,
            'shared_policy': policy,
            'live_controls': controls,
            'shared_policy_root': str(policy_root.parent) if policy_root else None,
            'directory': str(root), 'scheduler_error': queue_error, 'process_error': process_error,
            'warnings': warnings, 'cases': [case_snapshot(d, jobs, processes, now,
                log_limit=200 if include_all else 40, log_lines=log_lines) for d in directories]}


def age(seconds):
    if seconds < 60:
        return f'{seconds:.0f}s'
    if seconds < 3600:
        return f'{seconds/60:.1f}m'
    return f'{seconds/3600:.1f}h'


def relative_path(value, root):
    """Shorten display paths only; snapshot/JSON retain original metadata."""
    try:
        return str(Path(value).relative_to(root))
    except ValueError:
        try:
            # Adopted allocations can still report the old run's symlink.
            return str(Path(value).resolve().relative_to(Path(root).resolve()))
        except ValueError:
            return str(value)


def display_path(value, root):
    """Compact task paths; keep state/channel identity and the original in JSON."""
    path = relative_path(value, root)
    for marker, label in (('/anl_automatic/states/', 'ANL/'),
                          ('/vrctst/molpro_scan/', 'VRC/')):
        if marker in path:
            path = label + path.split(marker, 1)[1]
            path = path.replace('/composite/', '/').replace('/tasks/', '/')
            if label == 'VRC/':
                path = re.sub(r'^VRC/\d+_', 'VRC/', path)
            break
    return path


def log_group(path):
    """Native force steps of one dispatch task share a display slot."""
    parts = Path(path).parts
    if 'tasks' in parts:
        index = parts.index('tasks')
        if len(parts) > index + 2:
            return str(Path(*parts[:index+2]))
    return path


def process_label(command):
    try:
        words = shlex.split(command)
    except ValueError:
        return command
    if '-m' in words and words.index('-m') + 1 < len(words):
        return words[words.index('-m') + 1]
    return Path(words[0]).name if words else 'process'


def log_content(log, *, lines=4, raw=False, details=False):
    if raw:
        return log.get('raw_tail', [])[-lines:]
    historical = set() if details else set(log.get('historical_events', []))
    events = [s for s in log['events'] if s not in historical]
    errors = [s for s in log['errors'] if s not in historical]
    content = events[-lines:] + log['progress'][-lines:] + errors
    if log.get('table_header') and any(
            ITERATION.match(line) or re.match(r'^(?:Sella|BFGS|LBFGS|FIRE)(?::|\s)', line)
            for line in content):
        content.insert(0, log['table_header'])
    return list(dict.fromkeys(content))


def render_logs(case, *, include_all=False, lines=4, raw=False, details=False):
    """Keep controller tails above QC, irrespective of native output volume."""
    out = []
    for controller, title in ((True, 'KINBOT LOGS'), (False, 'QC OUTPUT')):
        logs = [log for log in case['logs'] if (Path(log['path']).name in LOG_NAMES) == controller]
        if not details:
            logs = [log for log in logs if not log.get('completed_worker')]
        if controller:
            logs.sort(key=lambda r: (Path(r['path']).name != 'pes.log',
                                    len(Path(r['path']).parts), r['path']))
        else:
            logs.sort(key=lambda r: (bool(r['errors']), r.get('live_task', False), r['mtime']), reverse=True)
        shown = 0
        groups = set()
        for log in logs:
            group = log_group(log['path'])
            if not details and group in groups and not log['errors']:
                continue
            content = log_content(log, lines=lines, raw=raw, details=details)
            if not content:
                continue
            if not shown:
                out.append('  '+title+(' (raw tails)' if raw else ''))
            if shown >= (200 if details else (6 if include_all else 4)):
                out.append('    More log tails available with --all --details (or --json).')
                break
            shown += 1
            groups.add(group)
            path = log['path'] if details else display_path(log['path'], case['directory'])
            out.append(f'    {path} | {age(log["age_seconds"])} ago' +
                       (f' | {log["bytes"]:,} B' if details else ''))
            for line in content:
                # Show the event's clock and severity; date/root are already
                # in the dashboard header. Raw tails and JSON stay untouched.
                if not raw and EVENT.match(line):
                    line = line.split(' ', 1)[1].replace(case['directory']+'/', '')
                out.append('      '+line)
    return out


def render_status_summary(case):
    """Bound repeated failures and dependency cascades; full records stay in JSON."""
    out, counts = [], Counter()
    graphs = [item for item in case['artifacts'] if 'tasks' in item]
    owners, failures, retries, intermediate = {}, {}, [], set()
    live = {str(j['job_id']) for j in case['jobs']}
    blocked = 0
    def failure(label, error):
        error = ' '.join(str(error or 'No diagnostic recorded').split())
        failures.setdefault(error, set()).add(label)
    for graph in graphs:
        counts.update(graph['counts'])
        parent = Path(graph['path']).parent
        for task in graph['tasks']:
            owners[parent/'tasks'/task['name']] = task
            label = display_path(str(parent/'tasks'/task['name']), case['directory'])
            if task['status'] in ('failed', 'unavailable'):
                failure(label, task.get('error'))
            if task['status'] == 'blocked' or any(
                    s in ('failed', 'unavailable', 'blocked') for s in task['blocked_by'].values()):
                blocked += 1
            if task.get('recovery') and task['status'] not in ('complete', 'skipped', 'blocked'):
                retries.append(f'  RETRY {label}: {task["status"]}' +
                               (f'; job={task["job_id"]}' if task['job_id'] else ''))
    if graphs:
        out.append(f'Dispatch tasks ({len(graphs)} graphs): ' +
                   '  '.join(f'{k}={v}' for k, v in sorted(counts.items())))
    for item in case['artifacts']:
        label = display_path(item['path'], case['directory'])
        if item.get('stale_review_hold'):
            out.append(f'  VRC REVIEW {label}: listed channel models accepted; '
                       'old aggregate hold superseded (not a QC completion verdict).')
            continue
        if item['status']=='needs_review':
            out.append(f'  REVIEW HOLD {label}: {item.get("error", "reviewed input required")}')
            continue
        if item.get('previous_controller'):
            if item.get('branches'):
                worker = str(Path(item['path']).parent)
                stage = case.get('exploration', {}).get('workers', {}).get(worker, {}).get('status')
                out.append(f'Branch report for {worker} belongs to the previous controller; '
                           + (f'current PES stage={stage}; ' if stage else '') +
                           'current outcomes not reported yet (--details retains history).')
            continue
        if item.get('branches'):
            for name, branch in item['branches'].items():
                if branch.get('error'):
                    out.append(f'  {name} [{branch.get("status", "unknown")}]: {branch["error"]}')
            if item.get('mess_gate') is not None:
                out.append(f'MESS gate: {item["mess_gate"]}; solver: {item.get("final_solver")}')
        owner = owners.get(Path(item['path']).parent)
        if owner and (item.get('error') or item['status'] in ('failed', 'unavailable')):
            if owner['status'] in ('failed', 'unavailable'):
                continue  # Native and minimum-recovery records describe this same task.
            if owner['status'] == 'submitted' and str(owner['job_id']) in live:
                intermediate.add(str(Path(item['path']).parent))
                continue  # A live retry can still have its earlier failed execution record.
        if not item.get('branches') and not item.get('preflight') and not item.get('tasks') and (
                item.get('error') or item['status'] in ('failed', 'unavailable', 'blocked')):
            failure(label, item.get('error'))
        elif item.get('recovery') and item['status'] in ('retrying', 'running'):
            retries.append(f'  RECOVERY {label}: {item["status"]}')
        elif item.get('recovery', {}).get('selected'):
            retries.append(f'  RECOVERY {label}: {item["status"]}; selected={item["recovery"]["selected"]}')
    if failures:
        out.append(f'Failures: {sum(len(v) for v in failures.values())} tasks/records; '
                   f'{len(failures)} distinct diagnostics')
        for error, labels in sorted(failures.items(), key=lambda pair: (-len(pair[1]), pair[0]))[:5]:
            excerpt = error if len(error) <= 240 else error[:237] + '...'
            out.append(f'  {len(labels)} affected: {excerpt}')
            out.append('    '+', '.join(sorted(labels)[:2]) +
                       (f'; +{len(labels)-2} more' if len(labels) > 2 else ''))
        if len(failures) > 5:
            out.append(f'  {len(failures)-5} additional diagnostics in --details / --json.')
    if blocked:
        out.append(f'Dependency blocks: {blocked} tasks waiting on failed prerequisites (not independent failures).')
    if intermediate:
        out.append(f'Live tasks with earlier/intermediate failure records: {len(intermediate)}; final outcomes pending.')
    out.extend(retries[:3])
    if len(retries) > 3:
        out.append(f'  {len(retries)-3} additional retry/recovery records in --details / --json.')
    return out


def wrap_display(out, width):
    """Wrap, never truncate, long diagnostics; keep continuation indentation."""
    wrapped = []
    for line in out:
        indent = len(line) - len(line.lstrip())
        wrapped.extend(textwrap.wrap(line, width=max(60, width),
                       subsequent_indent=' '*(indent+2),
                       replace_whitespace=False, drop_whitespace=True) or [''])
    return '\n'.join(wrapped)


def render(report, *, include_all=False, lines=4, details=False, raw_logs=False, width=None):
    out = [f'KINBOT  {report["time"]}', report['directory']]
    for key in ('scheduler_error', 'process_error'):
        if report[key]:
            out.append(f'{key.upper()}: UNKNOWN — {report[key]}')
    common_exclusions = report['cases'][0]['configured_exclusions'] if report['cases'] else []
    same_exclusions = all(c['configured_exclusions'] == common_exclusions for c in report['cases'])
    shared = report.get('shared_policy')
    controls = report.get('live_controls')
    if controls:
        if controls.get('error'):
            out.append('Live controls: UNKNOWN — '+controls['error'])
        else:
            out.append(f'Live controls: revision {controls["revision"]}; '
                       f'hold failures={controls["hold_failures"]}; requests={controls["requests"]}')
            if controls.get('drain',{}).get('status') == 'blocked':
                out.append('Pending drain unconfirmed: '+str(controls['drain']['error']))
    if shared:
        out.append(f'Shared pool: cap {shared["max_jobs"]} across all runs/branches under '
                   f'{report["shared_policy_root"]} (revision {shared["revision"]})')
        out.append('New-submission exclusions: '+(','.join(shared['excluded_nodes']) or '(none)'))
    elif same_exclusions:
        out.append('Configured exclusions: '+(','.join(common_exclusions) or '(none)'))
    # Show every PES before expanding either one; a long failed QC log must
    # never make the other run invisible on a terminal/watch display.
    for case in report['cases']:
        branches = [item for item in case['artifacts'] if item.get('branches')]
        out.append(f'OVERVIEW {case["name"]}: '
                   f'jobs={len(case["jobs"]) if not report["scheduler_error"] else "UNKNOWN"}; '
                   f'controllers={len(case["processes"]) if not report["process_error"] else "UNKNOWN"}'
                   + ('' if branches else '; no post-reaction branch record yet'))
        for item in branches:
            worker = str(Path(item['path']).parent)
            states = '; '.join(f'{name}=' + ('not_reported_yet' if item.get('previous_controller')
                               else value.get('status', 'unknown'))
                               for name, value in item['branches'].items())
            stage = case.get('exploration', {}).get('workers', {}).get(worker, {}).get('status')
            out.append(f'  {worker}: {states}' + (f'; PES stage={stage}' if stage else ''))
    for case in report['cases']:
        out += ['', f'=== {case["name"]} ===']
        exploration = case.get('exploration', {})
        if exploration:
            counts = dict(Counter(w.get('status', 'unknown') for w in exploration.get('workers', {}).values()))
            owner_alive = any(p['pid']==exploration.get('owner_pid') for p in case['processes'])
            out.append(f'PES exploration: {exploration.get("phase", "unknown")}; '
                       f'epoch={exploration.get("epoch")}; wells={counts}; '
                       f'owner={"present" if owner_alive else "not observed"}')
            if exploration.get('quiescing'):
                owners = sorted(n for n, w in exploration.get('workers', {}).items()
                                if w.get('status') in ('vrc', 'postprocessing'))
                out.append('PES revision waiting for branch controllers to yield: '
                           + (', '.join(owners) if owners else 'owner acknowledgment')
                           + '. Existing Slurm jobs are retained.')
            if 'max_searches' in exploration:
                exploring = exploration.get('exploring_workers', [])
                waiting = exploration.get('waiting_workers', [])
                out.append(f'Exploration slots: {len(exploring)}/{exploration["max_searches"]}; '
                           f'unstarted wells={len(waiting)} (separate from the QC node pool).')
                if waiting and len(exploring) >= exploration['max_searches']:
                    out.append('Next well waits for an exploring well to finish: '+', '.join(exploring))
            selected = case.get('anl_selection', {})
            current = (selected.get('token'), selected.get('epoch')) == (
                exploration.get('token'), exploration.get('epoch'))
            if current and selected.get('status') == 'selected':
                chosen = sum(len(w['selected']) for w in selected['workers'].values())
                omitted = sum(len(w['omitted']) for w in selected['workers'].values())
                out.append(f'Outer PES L2 selection: {chosen} retained; {omitted} omitted by pathway policy.')
            elif current and selected.get('status') == 'blocked':
                out.append('ANL selection blocked: '+selected.get('error', 'incomplete L2 evidence'))
            elif (exploration.get('phase') == 'exploring'
                  and exploration.get('anl_selection_policy') == 'lowest_l2'):
                out.append('New ANL work waits for outer PES L2 pathway selection; retained QC may still be running.')
        out.append(f'Slurm: {"UNKNOWN" if report["scheduler_error"] else len(case["jobs"])} allocations'
                   + (' / shared pool; ' if shared else f' / cap {case["allocation_cap"] or "unspecified"}; ')
                   + f'local processes: {"UNKNOWN" if report["process_error"] else len(case["processes"])}')
        if not same_exclusions and not shared:
            out.append('Configured exclusions: '+(','.join(case['configured_exclusions']) or '(none)'))
        if case['jobs']:
            out.append('JOB       STATE       NODES       ELAPSED / LIMIT   CPU  MEMORY  NAME')
        for job in case['jobs']:
            out.append(f'{job["job_id"]:<9} {job["state"]:<11} {job["nodes"]:<11} '
                       f'{job["elapsed"]} / {job["limit"]}  {job["cpus"]}  {job["memory"]}  {job["name"]}')
            out.append('  WorkDir: '+(job['directory'] if details else display_path(job['directory'], case['directory'])))
            if job['state'] != 'RUNNING':
                out.append('  Reason: '+job['reason'])
        # squeue uses a compressed host list; display that exact value once
        # per group instead of repeating it for every allocation.
        for exclusion in dict.fromkeys(j['excluded_nodes'] for j in case['jobs']):
            ids = ','.join(j['job_id'] for j in case['jobs'] if j['excluded_nodes'] == exclusion)
            out.append(f'  Slurm exclusions: {exclusion or "(none)"}' +
                       (f' [jobs {ids}]' if details else ''))
        if case['processes']:
            if details:
                for process in case['processes']:
                    out.append(f'PID {process["pid"]} [{process["state"]}] {process["command"]}')
            else:
                out.append('Local: '+'; '.join(f'{process_label(p["command"])} PID {p["pid"]} [{p["state"]}]'
                                              for p in case['processes']))
        preflights = [item for item in case['artifacts'] if item.get('preflight')]
        for item in preflights:
            p = item['preflight']
            alive = any(proc['pid'] == p['controller_pid'] for proc in case['processes'])
            status = item['status']
            if item.get('previous_controller'):
                out.append('Composite preflight: previous controller report; '
                           'current preflight not reported yet')
                continue
            if status == 'running' and not alive:
                status = 'controller not observed; recorded state only'
            out.append(f'Composite preflight: {status}; graphs checked '
                       f'{p["graphs_complete"]}/{p["graphs_total"]}; {p["phase"]}')
            if item['status'] != 'complete':
                updated = (age(max(0, report.get('timestamp', time.time())-p['updated_at']))
                           if isinstance(p.get('updated_at'), (int, float)) else 'unknown')
                out.append(f'  {p["graph"]}; task={p["task"] or "—"}; '
                           f'updated {updated} ago')
                if p.get('detail'):
                    out.append('  '+p['detail'])
            if item.get('error'):
                out.append('  '+item['error'])
        # Older pinned controllers have no progress JSON. Recognize their
        # exact probe command without executing any dispatcher code.
        if not preflights:
            for process in case['processes']:
                _, separator, body = process['command'].partition(' -c ')
                if separator and body == 'import ase, sella; import kinbot.anl.dispatch':
                    out.append('Composite preflight: checking Python environment in '
                               +relative_path(process['cwd'], case['directory']))
        if case['pes']['counts']:
            out.append('PES workers: '+'  '.join(f'{k.lower()}={v}' for k, v in case['pes']['counts'].items()))
        failures = [r for r in case['exit_records'] if r['returncode']
                    and not r.get('archived_failure')]
        historical = [r for r in case['exit_records'] if r['returncode']
                      and r.get('archived_failure')]
        if historical:
            out.append(f'Cached rejected QC outcomes: {len(historical)}; native attempts archived before restart.')
        if case.get('archived_attempt_directories'):
            out.append(f'Archived task attempts excluded: {len(case["archived_attempt_directories"])} graph(s).')
        if case.get('retired_scan_directories'):
            out.append(f'Old VRC scan graphs retained for comparison: {len(case["retired_scan_directories"])} (excluded from current progress).')
        for failure in failures:
            out.append(f'FAILED QC EXIT: {failure["path"]} job={failure["job_id"] or "unknown"} exit={failure["returncode"]}')
        workers = case['workers']
        if not details:
            completed = [w for w in workers if '/' in w['directory']
                         and w['status'].get('status') == 'complete']
            workers = [w for w in workers if w not in completed]
            if completed:
                out.append(f'Reference/subtask workers: {len(completed)} complete (--details lists them).')
        for worker in workers:
            status = worker['status'].get('status', 'no completion record')
            out.append(f'WORKER {worker["directory"]}: {status}; PES says {worker["reported_pes_state"] or "unlisted"}')
            if worker['reactions']:
                out.append('  Reaction search outcomes: '+str(worker['reactions']['counts']))
                groups = worker['reactions'].get('accepted_product_groups')
                if groups:
                    out.append(f'  Accepted searches reach {groups} product group(s); counts are not unique saddles.')
                rows = worker['reactions']['rows']
                for row in rows if include_all or details else [r for r in rows if r['code'] >= 0][:12]:
                    out.append(f'    {row["name"]}: {row["stage"]}; step={row["step"]}')
            for row in worker['optimizations']:
                out.append(f'  Optimizing {row["reaction"]}: pending={row.get("pending")} failures={row.get("failures")}')
                if include_all:
                    out.append('    '+json.dumps(row.get('optimizers', [])))
        if not details and not include_all:
            out.extend(render_status_summary(case))
        items = [item for item in case['artifacts']
                 if details or not item.get('preflight')]
        if not details and not include_all:
            items = []
        if not details:
            healthy = [item for item in items if not item.get('error')
                       and not item.get('branches') and not item.get('tasks')
                       and item['status'] in ('normal', 'complete', 'running', 'recorded')
                       and not item.get('recovery', {}).get('completed_attempts')
                       and not item.get('recovery', {}).get('reason')]
            items = [item for item in items if item not in healthy]
            if healthy:
                out.append(f'Other status/result records: {dict(Counter(i["status"] for i in healthy))}')
            items = sorted(items, key=lambda i: (not bool(i.get('branches')),
                           not bool(i.get('error') or any(t['status'] in ('failed', 'unavailable')
                                                        for t in i.get('tasks', [])))))
            if len(items) > 8:
                out.append(f'Displaying 8/{len(items)} detailed records; all retained in --details / --json.')
            items = items[:8]
        for item in items:
            if not details and not item.get('error') and Path(item['path']).name == 'kinbot_status.json':
                continue  # Already in the WORKER row above.
            historical = ' [previous controller]' if item.get('previous_controller') else ''
            out.append(f'RECORD {item["path"]}{historical}: {item["status"]}')
            for key in ('counts', 'branches', 'mess_gate', 'final_solver', 'recovery', 'error',
                        'review_channels', 'recorded_status', 'recorded_error'):
                if item.get(key) is not None:
                    value = item[key]
                    if key == 'recovery' and not details and 'curvature_checks' in value:
                        out.append(f'  curvature_checks: {len(value["curvature_checks"])}; controls: {value.get("controls")}')
                        for check in value['curvature_checks']:
                            out.append('    '+str({k: v for k, v in check.items() if k in
                                ('attempt', 'status', 'lowest_frequency_cm_inverse', 'selected_energy_eV', 'error')}))
                    else:
                        out.append(f'  {key}: {value}')
            tasks = item.get('tasks', [])
            visible = tasks if details else [t for t in tasks if t['status'] in ('failed', 'unavailable', 'submitted')][:4]
            for task in visible:
                out.append(f'  {task["name"]}: {task["status"]} job={task["job_id"]}; waiting for={task["blocked_by"]}; {task["error"] or ""}')
            if not details:
                blocked = [t for t in tasks if t['blocked_by']]
                if blocked:
                    blockers = {k: v for t in blocked for k, v in t['blocked_by'].items()}
                    out.append(f'  {len(blocked)} dependent tasks waiting for={blockers}')
        for sampling in case['sampling']:
            out.append(f'ROTD {sampling["directory"]}: {sampling["surface_output_files"]}/'
                       f'{sampling["surface_directories"]} surface output files; '
                       'file counts only, final flux validation is separate.')
        out.extend(render_logs(case, include_all=include_all, lines=lines, raw=raw_logs, details=details))
        coverage = case.get('log_coverage', {})
        if details and coverage:
            out.append(f'Log coverage: {coverage["read"]}/{coverage["available"]} bounded tails read.')
        out.extend('READ WARNING: '+warning for warning in case['warnings'])
        if not case['jobs'] and not case['processes'] and not report['scheduler_error'] and not report['process_error']:
            out.append('INACTIVE: no live jobs/controllers observed. This is not a completion verdict.')
    out.extend('READ WARNING: '+warning for warning in report['warnings'])
    out += ['', 'Read-only | --all: reaction history / more logs | --details: full current records | --raw-logs: unfiltered tails | --json: snapshot',
            'Archived task attempts are excluded; retained files are unchanged. Log reads are bounded (40; --all: 200).',
            'Tails may include earlier errors; quiet output is not proof of a stall. Statuses are not a scientific audit.']
    if any(c['pes']['counts'] for c in report['cases']):
        out.append('Legacy PES Finished means exited, not successful.')
    return wrap_display(out, width or min(160, shutil.get_terminal_size((120, 24)).columns))


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
    parser.add_argument('--all', action='store_true', help='Show all reaction rows and more log tails; use --details for every task/record.')
    parser.add_argument('--lines', type=int, default=4, help='Recent event/native progress lines per log (default 4).')
    parser.add_argument('--details', action='store_true', help='Expand process commands, paths and recovery metadata.')
    parser.add_argument('--raw-logs', action='store_true', help='Show unfiltered last --lines lines instead of selected log progress.')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--json', action='store_true', help='Machine-readable snapshot.')
    mode.add_argument('--follow', action='store_true', help='Snapshot then follow new controller and QC log text.')
    parser.add_argument('--interval', type=float, default=3, help='Log-follow poll interval, seconds (minimum 1).')


def run(root, args):
    root = Path(root).expanduser().resolve(strict=True)
    if args.lines < 1 or not 1 <= args.interval <= 60:
        raise ValueError('--lines must be positive; --interval must be between 1 and 60 seconds.')
    report = snapshot(root, include_all=args.all or args.details, log_lines=args.lines)
    print(json.dumps(report, indent=2) if args.json else render(report, include_all=args.all,
          lines=args.lines, details=args.details, raw_logs=args.raw_logs), flush=True)
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
