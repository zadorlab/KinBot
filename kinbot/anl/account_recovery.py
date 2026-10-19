"""One-time migration of the legacy compute-node passwd lookup failure.

This repairs a KinBot startup exception, not a native convergence failure.
Completed results and live allocations are never reset by this migration.
"""
from copy import deepcopy
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from kinbot.execution import atomic_json, file_sha256


def resume(run_dir):
    """Finish a journaled archive before the dispatcher reads its ledger."""
    marker = Path(run_dir) / 'account_lookup_recovery.json'
    if not marker.exists():
        return
    record = json.loads(marker.read_text())
    if record['phase'] != 'preparing':
        return
    from kinbot.run_recover import _apply_graph
    _apply_graph(record['plan'], Path(run_dir).resolve(), Path(record['archive']))
    record['phase'] = 'complete'
    atomic_json(marker, record)


def recover(run_dir, spec, state):
    from kinbot.anl import dispatch
    from kinbot.anl.salvage import descendants
    from kinbot.qc_scratch import account_identity

    selected, reasons = [], {}
    tasks = {t['id']: t for t in spec['tasks']}
    for ident, entry in state['tasks'].items():
        if (entry['status'] not in ('failed', 'unavailable')
                or ident in state.get('account_lookup_retries', {})):
            continue
        path = run_dir / 'tasks' / ident / 'execution.json'
        if not path.exists():
            continue
        outcome = json.loads(path.read_text())
        error = outcome.get('error', '')
        trace = outcome.get('traceback', '')
        # This KeyError was incorrectly labelled "validation" by old workers.
        # Require both its exact message and its actual scratch-code traceback.
        if (outcome.get('status') != 'failed'
                or not re.fullmatch(r"['\"]?getpwuid\(\): uid not found: [0-9]+['\"]?", error)
                or 'qc_scratch.py' not in trace or 'getpwuid' not in trace):
            continue
        dispatch._verify_stage_files(run_dir, tasks[ident], entry)
        dispatch._verify_execution(run_dir, tasks[ident], entry, outcome)
        if entry.get('job_id') and dispatch._job_active(entry['job_id']):
            continue  # A live retry's intermediate failure is not terminal.
        consumers = descendants(spec, ident)
        if any(state['tasks'].get(child, {}).get('status') in
               ('complete', 'submitted', 'submitting', 'staged') for child in consumers):
            raise RuntimeError(f'{ident}: account-lookup retry has existing downstream '
                               'work; review geometry provenance before resetting it.')
        selected.append(ident)
        reasons[ident] = 'Retry legacy passwd lookup failure with exported account-name support.'
    if not selected:
        return state
    account_identity()  # Refuse migration if the replacement cannot resolve the account.
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    archive = run_dir / 'attempts' / ('account_lookup_' + stamp)
    archive.mkdir(parents=True)
    updated = deepcopy(state)
    retries = updated.setdefault('account_lookup_retries', {})
    retries.update({ident: str(archive) for ident in selected})
    plan = {'graph': str(run_dir), 'before': {name: file_sha256(run_dir/name)
            for name in ('workflow.json', 'state.json')}, 'spec': spec, 'state': updated,
            'reset': selected, 'reasons': reasons,
            'complete': [k for k, v in state['tasks'].items() if v['status'] == 'complete'],
            'active': {k: v['job_id'] for k, v in state['tasks'].items() if v['status'] == 'submitted'}}
    marker = run_dir / 'account_lookup_recovery.json'
    if marker.exists():
        atomic_json(archive/'previous_recovery.json', json.loads(marker.read_text()))
    atomic_json(marker, {'phase': 'preparing', 'archive': str(archive), 'plan': plan})
    resume(run_dir)
    return dispatch._load(run_dir)[2]
