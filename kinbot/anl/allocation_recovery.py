"""Persistent, bounded retries of positively identified infrastructure failures.

Unknown disappearance, user cancellation and native convergence failures do not
authorize retries. Scheduler outages leave the durable wait and attempt budget
unchanged. This module never changes the scientific task or resource request.
"""
from copy import deepcopy
from datetime import datetime, timezone
import json
from pathlib import Path
import time
import uuid

from kinbot.execution import atomic_json, file_sha256, slurm_job_record
from kinbot.run_policy import read as read_policy

STATES = {'BOOT_FAIL', 'NODE_FAIL', 'PREEMPTED'}


def evidence(directory, job):
    try:
        record = slurm_job_record(job, strict=True)
    except OSError as error:
        from kinbot.anl.dispatch import SchedulerUnavailableError
        raise SchedulerUnavailableError(str(error)) from error
    if (record.get('JobId') != str(job) or not record.get('WorkDir')
            or Path(record['WorkDir']).resolve() != directory.resolve()
            or record.get('JobState', '').split('+')[0] not in STATES):
        return None
    return dict(job_id=str(job), directory=str(directory.resolve()),
                scheduler_record=record, observed_at=time.time(), source='scontrol')


def resume(graph):
    marker = Path(graph)/'allocation_recovery.json'
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

    settings = read_policy(graph) or {}
    limit = settings.get('allocation_retries', 0)
    if not limit and not any(e.get('status') == 'retry_wait' for e in state['tasks'].values()):
        return state
    # A deliberately bounded diagnostic has a separate explicit retry policy.
    if spec.get('intent', {}).get('vrc_fragment_retry'):
        return state
    updated = deepcopy(state)
    tasks = {t['id']: t for t in spec['tasks']}
    selected, reasons = [], {}
    for ident, entry in state['tasks'].items():
        if entry.get('status') not in ('failed', 'unavailable', 'retry_wait'):
            continue
        directory = graph/'tasks'/ident
        job = _job_id(directory, ident, entry)
        if not job:
            continue
        current = updated['tasks'][ident]
        if dispatch._job_active(job):
            # An administrator may have requeued the same allocation while
            # KinBot was waiting. Reattach instead of creating another job.
            if entry.get('allocation_retry_wait', {}).get('job_id') == job:
                current.update(status='submitted', job_id=job)
                current.pop('allocation_retry_wait', None)
            continue
        wait = entry.get('allocation_retry_wait')
        proof = wait['evidence'] if wait and wait['job_id'] == job else evidence(directory, job)
        if proof is None:
            continue
        dispatch._verify_stage_files(graph, tasks[ident], entry)
        outcome = directory/'execution.json'
        if outcome.exists():
            execution = json.loads(outcome.read_text())
            dispatch._verify_execution(graph, tasks[ident], entry, execution)
            if execution.get('status') == 'executed':
                raise RuntimeError(f'{ident}: reconcile completed native work before allocation recovery.')
        current.update(status='failed', failure_phase='allocation', job_id=job,
                       allocation_evidence=proof)
        history = updated.setdefault('allocation_retries', {}).setdefault(ident, [])
        reason = None
        if len(history) >= limit:
            reason = f'Allocation retry limit reached ({limit}).'
        elif tasks[ident].get('restart'):
            reason = 'Explicit native restart needs its checkpoint provenance reviewed.'
        elif any(state['tasks'].get(child, {}).get('status') in
                 ('complete', 'submitted', 'submitting', 'staged', 'retry_wait')
                 for child in descendants(spec, ident)):
            reason = 'Existing downstream work requires provenance review.'
        if reason:
            current.pop('allocation_retry_wait', None)
            current['allocation_retry_blocked'] = reason
            continue
        if wait is None or wait['job_id'] != job:
            delay = min(settings.get('allocation_retry_delay', 300) * 2**len(history),
                        settings.get('allocation_retry_max_delay', 3600))
            wait = dict(job_id=job, evidence=proof, retry_after=time.time()+delay)
        current.update(status='retry_wait', allocation_retry_wait=wait)
        current.pop('allocation_retry_blocked', None)
        if time.time() < wait['retry_after']:
            continue
        history.append(dict(**wait, retired_at=datetime.now(timezone.utc).isoformat(),
                            mode='same_task_preserve_checkpoint'))
        selected.append(ident)
        reasons[ident] = f'Confirmed {proof["scheduler_record"]["JobState"]}; allocation retry {len(history)}/{limit}.'
    if not selected:
        return updated
    archive = graph/'attempts'/('allocation_'+uuid.uuid4().hex)
    archive.mkdir(parents=True)
    plan = dict(graph=str(graph), before={name:file_sha256(graph/name)
                for name in ('workflow.json', 'state.json')}, spec=spec, state=updated,
                reset=selected, reasons=reasons,
                complete=[k for k,v in state['tasks'].items() if v['status']=='complete'],
                active={k:v['job_id'] for k,v in state['tasks'].items() if v['status']=='submitted'})
    marker = graph/'allocation_recovery.json'
    if marker.exists():
        atomic_json(archive/'previous_recovery.json', json.loads(marker.read_text()))
    atomic_json(marker, dict(phase='preparing', archive=str(archive), plan=plan))
    resume(graph)
    return dispatch._load(graph)[2]
