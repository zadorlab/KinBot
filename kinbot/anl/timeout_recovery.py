"""Transactional clean retries of timed-out dispatcher tasks, for every backend."""
from copy import deepcopy
from datetime import datetime, timezone
import json
from pathlib import Path

from kinbot import timeout_retry
from kinbot.execution import atomic_json, file_sha256


def _job_id(directory, ident, entry):
    if entry.get('job_id'):
        return str(entry['job_id'])
    # Older optional-component handling dropped job_id when marking a task
    # unavailable. Recover only one unambiguous generation from native/guard
    # provenance, never the newest guessed filename or another task's ID.
    jobs = set()
    paths = [directory/'execution.json',
             *(directory/'.kinbot_slurm_io_retries').glob('*/*/execution.json')]
    for path in paths:
        if not path.is_file():
            continue
        record = json.loads(path.read_text())
        if record.get('task_id', record.get('identity')) != ident:
            continue
        value = str(record.get('slurm_job_id', ''))
        if value.isdecimal():
            jobs.add(value)
    return next(iter(jobs)) if len(jobs) == 1 else None


def resume(run_dir):
    marker = Path(run_dir) / 'timeout_recovery.json'
    if marker.exists():
        record = json.loads(marker.read_text())
        if record['phase'] == 'preparing':
            from kinbot.run_recover import _apply_graph
            _apply_graph(record['plan'], Path(run_dir).resolve(), Path(record['archive']))
            record['phase'] = 'complete'
            atomic_json(marker, record)


def recover(run_dir, spec, state):
    from kinbot.anl import dispatch
    from kinbot.anl.salvage import descendants

    if spec.get('intent', {}).get('vrc_fragment_retry'):
        return spec, state  # One explicitly bounded diagnostic allocation.
    settings = timeout_retry.policy(run_dir)
    limit = settings.get('timeout_retries', 0)
    if not limit:
        return spec, state
    selected, reasons = [], {}
    updated, replacement = deepcopy(state), deepcopy(spec)
    tasks = {t['id']: t for t in replacement['tasks']}
    for ident, entry in state['tasks'].items():
        if entry.get('status') not in ('failed', 'unavailable'):
            continue
        directory = run_dir / 'tasks' / ident
        job = _job_id(directory, ident, entry)
        if not job or dispatch._job_active(job):
            continue
        proof = timeout_retry.evidence(directory, job, [directory/'slurm.stderr'])
        if proof is None:
            continue
        dispatch._verify_stage_files(run_dir, tasks[ident], entry)
        execution_file = directory/'execution.json'
        if execution_file.exists():
            execution = json.loads(execution_file.read_text())
            dispatch._verify_execution(run_dir, tasks[ident], entry, execution)
            if execution.get('status') == 'executed':
                raise RuntimeError(f'{ident}: completed timeout result needs reconciliation before retry.')
        # Walltime exhaustion is not scientific unavailability. Prevent both
        # optional-component fallbacks and independent geometry retry budgets.
        current = updated['tasks'][ident]
        current.update(status='failed', failure_phase='walltime', timeout_evidence=proof, job_id=job)
        history = updated.setdefault('timeout_retries', {}).setdefault(ident, [])
        if len(history) >= limit:
            current['timeout_retry_blocked'] = f'Timeout retry limit reached ({limit}).'
            continue
        if any(state['tasks'].get(child, {}).get('status') in
               ('complete', 'submitted', 'submitting', 'staged') for child in descendants(spec, ident)):
            current['timeout_retry_blocked'] = 'Existing downstream work requires provenance review.'
            continue
        walltime = timeout_retry.next_walltime(settings, proof)
        if walltime is None:
            current['timeout_retry_blocked'] = 'Partition maximum prevents a longer allocation.'
            continue
        # A clean retry must not accidentally inherit a manual amplitude
        # restart requiring a checkpoint from the archived attempt.
        if tasks[ident].get('restart'):
            current['timeout_retry_blocked'] = 'Checkpoint restart task requires manual review.'
            continue
        tasks[ident]['resources']['walltime'] = walltime
        tasks[ident]['resources']['partition'] = proof['partition']
        history.append({**proof, 'walltime': walltime, 'mode': 'clean_retry'})
        selected.append(ident)
        reasons[ident] = (f'Confirmed Slurm timeout; clean retry {len(history)}/{limit} '
                          f'with walltime {walltime}.')
    if not selected:
        return spec, updated
    dispatch.validate_spec(replacement)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    archive = run_dir/'attempts'/('timeout_' + stamp)
    archive.mkdir(parents=True)
    plan = {'graph': str(run_dir), 'before': {name: file_sha256(run_dir/name)
            for name in ('workflow.json', 'state.json')}, 'spec': replacement, 'state': updated,
            'reset': selected, 'reasons': reasons,
            'complete': [k for k, v in state['tasks'].items() if v['status'] == 'complete'],
            'active': {k: v['job_id'] for k, v in state['tasks'].items() if v['status'] == 'submitted'}}
    marker = run_dir/'timeout_recovery.json'
    if marker.exists():
        atomic_json(archive/'previous_recovery.json', json.loads(marker.read_text()))
    atomic_json(marker, {'phase': 'preparing', 'archive': str(archive), 'plan': plan})
    resume(run_dir)
    _, spec, state = dispatch._load(run_dir)
    return spec, state
