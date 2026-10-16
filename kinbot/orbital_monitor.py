"""Read-only live campaign view: scheduler facts and bounded native log tails.

No audits, hashes, locks, retries or manifest writes. Log matches are evidence
to inspect, not a replacement for the controller's convergence validation.
"""
from collections import Counter
from datetime import datetime
import json
from pathlib import Path
import re
import subprocess
import time


_ERROR = re.compile(
    r'no convergence|not converged|error exit|\?\s*error|global error|'
    r'input/output error|stale file handle|no space left|disk quota exceeded|'
    r'error (?:writing|reading)\s+\d+\s+words|erroneous (?:write|read)|'
    r'\?\s*i/o error|traceback \(most recent|mpi_abort|segmentation fault|'
    r'out.of.memory|oom.kill|killed by signal|bad termination|permission denied', re.I)
_PROGRAM = re.compile(r'^\s*PROGRAM\s*\*\s*(.+)', re.I)


def _read_json(path):
    try:
        value = json.loads(path.read_text())
        if not isinstance(value, dict):
            raise ValueError('Expected a JSON object')
        return value, None
    except FileNotFoundError:
        return {}, None
    except (OSError, ValueError) as error:
        return {}, f'{path}: {error}'


def _tail(path, now, limit=131072):
    try:
        with path.open('rb') as stream:
            stream.seek(0, 2)
            size = stream.tell()
            stream.seek(max(0, size - limit))
            raw = stream.read(limit)
        # Do not match the first fragment of a truncated line.
        if size > limit:
            raw = raw.partition(b'\n')[2]
        return {'text': raw.decode(errors='replace'), 'bytes': size,
                'age_seconds': max(0., now - path.stat().st_mtime),
                'truncated': size > limit}
    except FileNotFoundError:
        return {'text': '', 'bytes': None, 'age_seconds': None}
    except OSError as error:
        return {'text': '', 'bytes': None, 'age_seconds': None, 'read_error': str(error)}


def _queue():
    try:
        result = subprocess.run(
            ['squeue', '-h', '-o', '%i|%T|%N|%M|%R|%m|%x|%Z'],
            capture_output=True, text=True, timeout=15, check=False)
        if result.returncode:
            return {}, result.stderr.strip() or f'squeue exit {result.returncode}'
    except (OSError, subprocess.TimeoutExpired) as error:
        return {}, str(error)
    rows = {}
    for line in result.stdout.splitlines():
        fields = [v.strip() for v in line.split('|', 7)]
        if len(fields) != 8:
            return {}, 'Unrecognized squeue output; scheduler state is unknown.'
        rows[fields[0]] = dict(zip(('state', 'nodes', 'elapsed', 'reason',
                                   'memory', 'excluded_nodes', 'directory'), fields[1:]))
    return rows, None


def snapshot(root, *, now=None):
    root = Path(root).expanduser().resolve()
    now = time.time() if now is None else now
    campaign, error = _read_json(root / 'campaign.json')
    if error or not campaign:
        raise ValueError(error or f'Campaign missing: {root}')
    queue, queue_error = _queue()
    rows, counts, warnings, excluded, retries = [], {}, [], set(), Counter()
    if queue_error:
        warnings.append('SCHEDULER UNAVAILABLE: ' + queue_error)
    for system in campaign['systems']:
        config, error = _read_json(Path(system['configuration']))
        if error or not config:
            warnings.append(error or f'Configuration missing: {system["configuration"]}')
        selected = campaign.get('selected_probe_method')
        if config and selected and config.get('probe_method', 'casci') != selected:
            continue
        for phase in ('pilot_manifest', 'manifest'):
            manifest, error = _read_json(Path(system[phase]))
            if error or not manifest:
                warnings.append(error or f'Manifest missing: {system[phase]}')
                continue
            excluded.update(manifest.get('slurm_excluded_nodes', []))
            count = counts.setdefault(phase, Counter())
            for job in manifest['jobs']:
                identifier = str(job.get('job_id', ''))
                q = queue.get(identifier, {})
                directory = Path(job['directory'])
                alerts = []
                if q and Path(q['directory']).resolve() != directory.resolve():
                    alerts.append('Slurm WorkDir differs from this probe; scheduler identity mismatch.')
                    q = {}
                state = q.get('state') or (
                    'FAILED' if job.get('terminal_failure') else
                    'VALIDATED' if job.get('status') == 'complete' else
                    'QUEUE-UNKNOWN' if identifier and queue_error else
                    'LEFT-QUEUE' if identifier else 'WAITING')
                count['total'] += 1
                count[state] += 1
                retries['startup'] += len(job.get('startup_retry_history', []))
                retries['manual_startup'] += len(job.get('manual_startup_retry_history', []))
                retries['allocation'] += len(job.get('allocation_retry_history', []))
                # Unsubmitted duplicate panels stay in the totals, not the
                # live table. Completed files were validated by the controller.
                if state in ('WAITING', 'VALIDATED'):
                    continue
                execution, error = _read_json(directory / 'execution.json')
                if error:
                    alerts.append(error)
                if job.get('terminal_failure'):
                    alerts.append(str(job['terminal_failure']).splitlines()[0][:300])
                if execution.get('error'):
                    alerts.append(str(execution['error']).splitlines()[0][:300])
                name = job['name']
                output = _tail(directory / (name + '.out'), now)
                lines = output['text'].splitlines()
                stage, iteration = '-', ''
                solver = {}
                for line in lines:
                    match = _PROGRAM.match(line)
                    if match:
                        stage, iteration = match[1].strip(), ''
                    elif re.match(r'^\s*\d+\s+\d+\s+\d+\s+[-+\d]+\.\d+', line):
                        iteration = line.strip()
                    for key, value in re.findall(r'\b(MAXIT|MAXITI|MAXVI|THRVAR|THRDEN)\s*=\s*(\S+)', line):
                        solver[key] = value
                guard_status, guard_retries, slurm_restarts = '-', 0, 0
                logs = [(name + '.out', output)]
                for suffix in ('.stderr', '.stdout', '.launcher.stderr', '.launcher.stdout'):
                    logs.append((name + suffix, _tail(directory / (name + suffix), now)))
                # Only the current allocation's retry journal. Old failures
                # are counted separately, never presented as a new log error.
                if identifier.isdigit():
                    for allocation in sorted((directory / '.kinbot_slurm_io_retries').glob('*/' + identifier)):
                        folders = [allocation, *[p for p in allocation.glob('restart-*')
                            if p.is_dir() and re.fullmatch(r'restart-[0-9]+', p.name)]]
                        latest = max(folders, key=lambda p: int(p.name.split('-')[1]) if p != allocation else 0)
                        path = latest / 'execution.json'
                        guard, error = _read_json(path)
                        if error:
                            alerts.append(error)
                        guard_status = guard.get('status', '?')
                        slurm_restarts = int(guard.get('slurm_restart_count', 0))
                        attempts = guard.get('attempts', [])
                        guard_retries = max(0, len(attempts) - (guard_status != 'running'))
                        current = len(attempts) + 1 if guard_status == 'running' else len(attempts)
                        for stream in ('stdout.txt', 'stderr.txt'):
                            p = path.parent / f'attempt-{current}' / stream
                            logs.append((str(p.relative_to(directory)), _tail(p, now)))
                matches, seen = [], set()
                for filename, log in logs:
                    if log.get('read_error'):
                        alerts.append(f'{filename}: {log["read_error"]}')
                    for line in log['text'].splitlines():
                        value = line.strip()[:400]
                        if _ERROR.search(line) and value not in seen:
                            seen.add(value)
                            matches.append({'file': filename, 'line': value})
                # Preserve the first native cause ahead of repetitive MPI
                # teardown messages; signal 9 during abort is not proof of OOM.
                matches.sort(key=lambda m: bool(re.search(
                    r'mpi_abort|global error|killed by signal|bad termination|traceback', m['line'], re.I)))
                host = q.get('nodes') or execution.get('host', '-')
                for node in manifest.get('slurm_excluded_nodes', []):
                    if q.get('state') == 'RUNNING' and node == host:
                        alerts.append(f'ASSIGNED TO EXCLUDED NODE: {node}')
                if q and manifest.get('slurm_excluded_nodes') and q.get('excluded_nodes') in ('', '(null)', 'N/A'):
                    alerts.append('Slurm has no node exclusion despite campaign exclusions.')
                rows.append({'system': system['title'], 'phase': phase, 'name': name,
                    'directory': str(directory), 'relative_directory': str(directory.relative_to(root)),
                    'job_id': identifier, 'state': state, 'scheduler': q,
                    'host': host, 'execution_status': execution.get('status', '-'),
                    'guard_status': guard_status, 'guard_retries': guard_retries,
                    'slurm_restarts': slurm_restarts,
                    'allocation_retries': len(job.get('allocation_retry_history', [])),
                    'startup_retries': len(job.get('startup_retry_history', [])),
                    'stage': stage, 'iteration': iteration,
                    'solver_parameters': solver,
                    'last_output': next((s.strip() for s in reversed(lines) if s.strip()), ''),
                    'output_bytes': output['bytes'], 'output_age_seconds': output['age_seconds'],
                    'tail_truncated': output.get('truncated', False),
                    'alerts': alerts, 'log_matches': matches[:8]})
    rows.sort(key=lambda r: (not bool(r['alerts'] or r['log_matches']), r['state'] != 'RUNNING',
                             int(r['job_id']) if r['job_id'].isdigit() else 0))
    controller = _tail(root / 'controller.log', now)
    if controller.get('read_error'):
        warnings.append('controller.log: ' + controller['read_error'])
    return {'root': str(root), 'timestamp': now, 'status': campaign['status'],
            'max_jobs': campaign.get('max_jobs'), 'counts': counts, 'jobs': rows,
            'retries': retries,
            'excluded_nodes': sorted(excluded), 'warnings': warnings,
            'controller_tail': [line[:300] for line in controller['text'].splitlines()[-6:]]}


def render(data, *, details=False):
    lines = [f'{datetime.fromtimestamp(data["timestamp"]).astimezone().isoformat(timespec="seconds")}  '
             f'campaign={data["status"]}  cap={data["max_jobs"]}', data['root']]
    for phase, count in data['counts'].items():
        lines.append(f'{phase}: ' + '  '.join(f'{k}={v}' for k, v in count.items()))
    lines.append('Campaign excluded nodes: ' + (','.join(data['excluded_nodes']) or 'none'))
    lines.append('Campaign retries (cumulative): ' + '  '.join(f'{key}={value}' for key, value in data['retries'].items()))
    lines.extend('WARNING: ' + warning for warning in data['warnings'])
    lines.append('JOB       STATE          HOST       ELAPSED    OUT-AGE  MEM-REQ  EXCL(SLURM)  NATIVE/GUARD(RETRIES)')
    for row in data['jobs']:
        q = row['scheduler']
        age = row['output_age_seconds']
        age = '-' if age is None else f'{int(age)}s'
        lines.append(f'{row["job_id"]:<9} {row["state"]:<14} {row["host"]:<10} '
                     f'{q.get("elapsed", "-"):<10} {age:>7} {q.get("memory", "-"):>8} '
                     f'{q.get("excluded_nodes", "-")}  {row["execution_status"]}/'
                     f'{row["guard_status"]}({row["guard_retries"]}) '
                     f'requeues={row["slurm_restarts"]} allocation-retries={row["allocation_retries"]}')
        lines.append('  ' + row['relative_directory'])
        if q.get('state') == 'PENDING':
            lines.append('  Waiting: ' + q.get('reason', '?'))
            if row['log_matches']:
                lines.append('  On-disk diagnostics below may precede a pending requeue.')
        elif details:
            lines.append('  Stage: ' + row['stage'])
            if row['solver_parameters']:
                lines.append('  Native solver: ' + '  '.join(f'{k}={v}' for k, v in row['solver_parameters'].items()))
            lines.append('  ' + ('Iteration: ' + row['iteration'] if row['iteration']
                                else 'Last output: ' + (row['last_output'] or '(none yet)')))
        for alert in row['alerts']:
            lines.append('  ALERT: ' + alert)
        for match in row['log_matches'][:4]:
            lines.append(f'  LOG [{match["file"]}]: {match["line"]}')
    lines.append('OUT-AGE = seconds since output write; quiet output alone does not prove a stall.')
    lines.append('LOG = diagnostic in current file tails (128 KiB/file), not a terminal-failure verdict.')
    lines.append('Controller tail:')
    lines.extend(data['controller_tail'])
    return '\n'.join(lines)
