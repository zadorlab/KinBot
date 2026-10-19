"""Bounded, provenance-checked allocation retries for unusable local scratch."""
from copy import deepcopy
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from kinbot.execution import atomic_json, file_sha256, slurm_job_record
from kinbot.qc_scratch import scratch_retry_limit


def resume(run_dir):
    marker = Path(run_dir) / 'scratch_recovery.json'
    if not marker.exists():
        return
    record = json.loads(marker.read_text())
    if record['phase'] == 'preparing':
        from kinbot.run_recover import _apply_graph
        _apply_graph(record['plan'], Path(run_dir).resolve(), Path(record['archive']))
        record['phase'] = 'complete'
        atomic_json(marker, record)


def eligible(execution):
    if execution.get('status') != 'failed':
        return False
    if 'retryable_scratch_node' in execution:
        return execution.get('failure_phase') == 'scratch_setup' and execution['retryable_scratch_node'] is True
    # Recognize old pinned workers by the actual selector traceback, not
    # arbitrary PermissionError in an executable, input or result directory.
    error, trace = execution.get('error', ''), execution.get('traceback', '')
    return ('No verified local QC scratch is usable:' in error
            and 'qc_scratch.py' in trace and 'allocate' in trace
            and bool(re.search(r'Permission denied|Operation not permitted|Input/output error|'
                               r'Stale file handle|Read-only file system|No space left on device|'
                               r'Disk quota exceeded|requires \d+ MB; \d+ MB free', error, re.I)))


def _node(directory, ident, entry, execution):
    job = str(entry.get('job_id') or execution.get('slurm_job_id') or '')
    if not job.isdecimal():
        return None
    node = None
    if str(execution.get('slurm_job_id')) == job:
        node = execution.get('execution_node')
    if not node:
        # The guard records hostname/job ID even if scratch allocation failed
        # before native runtime provenance was available in older workers.
        records = []
        for path in (directory / '.kinbot_slurm_io_retries').glob(f'*/{job}/execution.json'):
            record = json.loads(path.read_text())
            if (record.get('identity') == ident and str(record.get('slurm_job_id')) == job
                    and record.get('status') == 'failed'):
                records.append(record)
        if len(records) == 1:
            node = records[0].get('host')
    if not node:
        record = slurm_job_record(job)
        if (record.get('JobId') == job
                and Path(record.get('WorkDir', '/')).resolve() == directory.resolve()):
            node = record.get('NodeList')
    return node if isinstance(node, str) and re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', node) else None


def apply_exclusions(script, state):
    """Refresh unsubmitted scripts; preserve every original exclusion/resource."""
    from kinbot.qc import _load_slurm_node_state, _slurm_node_names
    nodes = _slurm_node_names(state.get('scratch_excluded_nodes', []))
    nodes.update(_load_slurm_node_state())
    if not nodes:
        return script
    pattern = r'^#SBATCH\s+--exclude(?:=|\s+)([^\n]+)\n?'
    existing = re.findall(pattern, script, flags=re.M)
    # Preserve bracketed Slurm hostlists while deduplicating explicit names.
    # This must be idempotent: a staged script can be inspected many times.
    values = sorted({token.strip() for value in existing
                     for token in re.split(r',(?=[^\]]*(?:\[|$))', value)
                     if token.strip()} | nodes)
    directive = '#SBATCH --exclude=' + ','.join(values) + '\n'
    stripped = re.sub(pattern, '', script, flags=re.M)
    first, separator, rest = stripped.partition('\n')
    return first + separator + directive + rest


def recover(run_dir, spec, state):
    from kinbot.anl import dispatch
    from kinbot.anl.salvage import descendants
    from kinbot.qc import _export_slurm_excluded_nodes

    limit = scratch_retry_limit()
    _export_slurm_excluded_nodes(state.get('scratch_excluded_nodes', []))
    selected, reasons, failures = [], {}, {}
    tasks = {t['id']: t for t in spec['tasks']}
    for ident, entry in state['tasks'].items():
        history = state.get('scratch_node_retries', {}).get(ident, [])
        if entry['status'] not in ('failed', 'unavailable'):
            continue
        directory = run_dir / 'tasks' / ident
        path = directory / 'execution.json'
        if not path.exists():
            continue
        execution = json.loads(path.read_text())
        if not eligible(execution):
            continue
        dispatch._verify_stage_files(run_dir, tasks[ident], entry)
        dispatch._verify_execution(run_dir, tasks[ident], entry, execution)
        if entry.get('job_id') and dispatch._job_active(entry['job_id']):
            continue
        # A site scratch fault is not evidence that a scientific component is
        # unavailable. Keep it a hard failure after this separate budget ends,
        # rather than triggering L3 retries or a lower composite fallback.
        entry.update(status='failed', failure_phase='scratch_setup')
        if len(history) >= limit:
            entry['scratch_retry_blocked'] = f'Scratch allocation retry limit reached ({limit}).'
            continue
        node = _node(directory, ident, entry, execution)
        if node is None:
            entry['scratch_retry_blocked'] = 'Execution node cannot be verified; automatic retry withheld.'
            continue
        if any(state['tasks'].get(child, {}).get('status') in
               ('complete', 'submitted', 'submitting', 'staged') for child in descendants(spec, ident)):
            raise RuntimeError(f'{ident}: scratch retry has existing downstream work; review provenance.')
        selected.append(ident)
        failures[ident] = {'node': node, 'job_id': entry.get('job_id'), 'error': execution.get('error', '')}
        reasons[ident] = f'Local scratch unavailable on {node}; retry on another node ({len(history)+1}/{limit}).'
    if not selected:
        return state
    updated = deepcopy(state)
    history = updated.setdefault('scratch_node_retries', {})
    for ident, evidence in failures.items():
        history.setdefault(ident, []).append(evidence)
    updated['scratch_excluded_nodes'] = sorted(set(updated.get('scratch_excluded_nodes', [])) |
                                              {e['node'] for e in failures.values()})
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    archive = run_dir / 'attempts' / ('scratch_node_' + stamp)
    archive.mkdir(parents=True)
    plan = {'graph': str(run_dir), 'before': {name: file_sha256(run_dir/name)
            for name in ('workflow.json', 'state.json')}, 'spec': spec, 'state': updated,
            'reset': selected, 'reasons': reasons,
            'complete': [k for k, v in state['tasks'].items() if v['status'] == 'complete'],
            'active': {k: v['job_id'] for k, v in state['tasks'].items() if v['status'] == 'submitted'}}
    marker = run_dir / 'scratch_recovery.json'
    if marker.exists():
        atomic_json(archive / 'previous_recovery.json', json.loads(marker.read_text()))
    atomic_json(marker, {'phase': 'preparing', 'archive': str(archive), 'plan': plan})
    resume(run_dir)
    _export_slurm_excluded_nodes(updated['scratch_excluded_nodes'])
    return dispatch._load(run_dir)[2]
