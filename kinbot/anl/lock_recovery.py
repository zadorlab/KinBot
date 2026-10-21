"""Recover the write-only NFS shared-lock startup bug once per task.

This is not a general I/O retry: only the exact pre-QC worker_load traceback
qualifies. Native errors and all completed/live allocations remain untouched.
"""
from copy import deepcopy
import json
from pathlib import Path
import uuid

from kinbot.execution import atomic_json, file_sha256


def eligible(outcome):
    trace = outcome.get('traceback', '')
    return (outcome.get('status') == 'failed'
            and outcome.get('error') == '[Errno 9] Bad file descriptor'
            and 'live_control.py' in trace and 'in worker_load' in trace
            and 'fcntl.flock(handle, fcntl.LOCK_SH)' in trace)


def resume(graph):
    marker = Path(graph)/'manifest_lock_recovery.json'
    if marker.exists():
        record = json.loads(marker.read_text())
        if record['phase'] == 'preparing':
            from kinbot.run_recover import _apply_graph
            _apply_graph(record['plan'], Path(graph).resolve(), Path(record['archive']))
            record['phase'] = 'complete'
            atomic_json(marker, record)


def recover(graph, spec, state):
    from kinbot.anl import dispatch
    from kinbot.anl.salvage import descendants
    from kinbot.anl.timeout_recovery import _job_id

    selected = []
    tasks = {t['id']: t for t in spec['tasks']}
    for ident, entry in state['tasks'].items():
        if (entry['status'] not in ('failed', 'unavailable')
                or ident in state.get('manifest_lock_retries', {})):
            continue
        directory = graph/'tasks'/ident
        path = directory/'execution.json'
        if not path.is_file():
            continue
        outcome = json.loads(path.read_text())
        if not eligible(outcome):
            continue
        job = _job_id(directory, ident, entry)
        if not job or dispatch._job_active(job):
            continue
        dispatch._verify_stage_files(graph, tasks[ident], entry)
        dispatch._verify_execution(graph, tasks[ident], entry, outcome)
        if any(state['tasks'].get(n, {}).get('status') in
               ('complete', 'submitted', 'submitting', 'staged')
               for n in descendants(spec, ident)):
            raise RuntimeError(f'{ident}: lock-startup recovery has downstream work; review provenance.')
        selected.append(ident)
    if not selected:
        return state
    archive = graph/'attempts'/('manifest_lock_'+uuid.uuid4().hex)
    updated = deepcopy(state)
    updated.setdefault('manifest_lock_retries', {}).update(
        {n: str(archive) for n in selected})
    plan = dict(graph=str(graph), spec=spec, state=updated, reset=selected,
                before={n: file_sha256(graph/n) for n in ('workflow.json', 'state.json')},
                reasons={n: 'Corrected NFS shared-lock descriptor; no native QC started.' for n in selected},
                complete=[n for n,e in state['tasks'].items() if e['status']=='complete'],
                active={n:e['job_id'] for n,e in state['tasks'].items()
                        if e['status']=='submitted' and n not in selected})
    marker = graph/'manifest_lock_recovery.json'
    if marker.exists():
        atomic_json(archive/'previous_recovery.json', json.loads(marker.read_text()))
    atomic_json(marker, dict(phase='preparing', archive=str(archive), plan=plan))
    resume(graph)
    return dispatch._load(graph)[2]
