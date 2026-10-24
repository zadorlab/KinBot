"""Reconcile pool reservations using verified dispatcher attempt evidence.

An absent job is not a failed submission. Unknown sbatch outcomes stay fenced.
Only a known allocation with a matching terminal execution can retire a purged
reservation. Older pre-admission state writes have a separate, stricter repair.
"""
from contextlib import nullcontext
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess

from kinbot.execution import atomic_json, file_sha256


def _read(path):
    return json.loads(Path(path).read_text())


def legacy_pending(claim):
    """Do not let another pool poll discard a legacy owner's blocking proof."""
    directory = Path(claim.get('directory', '/'))
    if directory.parent.name != 'tasks':
        return False
    try:
        entry = _read(directory.parent.parent/'state.json')['tasks'][directory.name]
        return (entry.get('status') == 'submitting' and not entry.get('job_id')
                and bool(entry.get('previous_attempt')))
    except (OSError, ValueError, KeyError, TypeError):
        return False


def _attempt(root, claim):
    """Locate this allocation in the current generation or its linked archive."""
    from kinbot.anl import dispatch
    job = str(claim.get('job_id', ''))
    directory = Path(claim.get('directory', '/')).resolve()
    if (not job.isdecimal() or claim.get('status') != 'submitted'
            or not directory.is_relative_to(root) or directory.parent.name != 'tasks'):
        return None
    receipt = root/'.kinbot_policy_submissions'/(job+'.json')
    if not receipt.is_file() or _read(receipt) != claim:
        return None
    graph = directory.parent.parent
    current_graph, spec, state = dispatch._load(graph)
    current = state['tasks'].get(directory.name, {})
    entry = current
    visited = set()
    while str(entry.get('job_id')) != job:
        previous = entry.get('previous_attempt')
        if not previous:
            return None
        archived = Path(previous).resolve()
        if (archived in visited or not archived.is_relative_to(graph/'attempts')
                or archived.name != directory.name or archived.parent.name != 'tasks'):
            return None
        visited.add(archived)
        current_graph = archived.parent.parent
        if directory.name not in _read(current_graph/'committed.json').get('reset', []):
            return None
        _, spec, state = dispatch._load(current_graph)
        entry = state['tasks'].get(directory.name, {})
    task = next(t for t in spec['tasks'] if t['id'] == directory.name)
    execution_path = current_graph/'tasks'/directory.name/'execution.json'
    execution = _read(execution_path)
    if (str(execution.get('slurm_job_id')) != job
            or execution.get('status') not in ('failed', 'executed')):
        return None
    dispatch._verify_stage_files(current_graph, task, entry)
    dispatch._verify_execution(current_graph, task, entry, execution)
    return dict(graph=graph, directory=directory, current=current, entry=entry,
                attempt=current_graph, receipt=receipt, execution=execution_path,
                execution_status=execution['status'])


def _purged(job):
    """Require an explicit invalid-ID reply, never an empty best-effort lookup."""
    try:
        result = subprocess.run(['scontrol', 'show', 'job', '-o', job],
                                capture_output=True, text=True, timeout=15)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if (result.returncode and not result.stdout.strip() and re.fullmatch(
            r'(?:slurm_load_jobs error:\s*)?Invalid job id specified',
            result.stderr.strip(), re.I)):
        return result.stderr.strip()
    return None


def terminal_evidence(root, claim, *, repair_pre_admission=False):
    """Read-only evidence; caller must first obtain a healthy queue snapshot."""
    try:
        attempt = _attempt(root, claim)
        if attempt is None:
            return None
        current = attempt['current']
        # Leave the blocking claim intact until its owner repairs an old
        # pre-admission state. A generic pool poll must not erase that proof.
        if current.get('status') == 'submitting' and not repair_pre_admission:
            return None
        diagnostic = _purged(str(claim['job_id']))
        if not diagnostic:
            return None
        return dict(kind='purged_allocation_with_verified_execution',
                    job_id=str(claim['job_id']), queue_absent=True,
                    checked_at=datetime.now(timezone.utc).isoformat(),
                    scheduler_diagnostic=diagnostic,
                    execution=str(attempt['execution']),
                    execution_sha256=file_sha256(attempt['execution']),
                    execution_status=attempt['execution_status'],
                    receipt=str(attempt['receipt']), receipt_sha256=file_sha256(attempt['receipt']))
    except (OSError, ValueError, KeyError, TypeError, RuntimeError, StopIteration):
        # A concurrent archive or any provenance mismatch retains the fence.
        return None


def reconcile_pre_admission(graph, *, pool_locked=False):
    """Repair only the legacy write-before-admission case, under owner locks.

The caller owns the graph (or has frozen its owner at a pool-locked boundary).
The old claim must still block this exact directory/command; its job must be a
verified archived terminal attempt. No current worker artifacts, competing
claim/receipt or live allocation may exist. Real sbatch uncertainty is retained.
"""
    from kinbot import run_policy, run_cleanup
    from kinbot.anl import dispatch
    graph, spec, state = dispatch._load(graph)
    pending = [n for n, e in state['tasks'].items() if e.get('status') == 'submitting']
    if not pending:
        return state
    policy = run_policy.find(graph)
    if policy is None:
        return state
    root = policy.parent
    with (nullcontext() if pool_locked else run_policy._locked(policy)):
        ledger = root/'.kinbot_policy_admissions.json'
        if not ledger.exists():
            return state
        claims = _read(ledger)
        jobs = run_cleanup._jobs(root)  # Healthy query, including work directories.
        tasks = {t['id']: t for t in spec['tasks']}
        for ident in pending:
            entry = state['tasks'][ident]
            directory = graph/'tasks'/ident
            if entry.get('job_id') or entry.get('submission_started_utc'):
                continue
            matching = [(token, c) for token, c in claims.items()
                        if c.get('directory') == str(directory)]
            if len(matching) != 1:
                continue
            token, claim = matching[0]
            if claim.get('command') != ['sbatch', '--parsable', 'job.slurm']:
                continue
            if any(j['job_id'] == str(claim.get('job_id'))
                   or Path(j['directory']).resolve() == directory for j in jobs):
                continue
            try:
                attempt = _attempt(root, claim)
                if attempt is None or attempt['attempt'] == graph:
                    continue
                previous = graph/'tasks'/ident
                allowed = {'task.json', 'geometry.xyz', 'job.slurm', tasks[ident].get('input_name')}
                if any(p.name not in allowed or not p.is_file() for p in previous.iterdir()):
                    continue
                dispatch._verify_stage_files(graph, tasks[ident], entry)
                submitted = datetime.fromisoformat(claim['created_at'])
                receipts = [_read(p) for p in (root/'.kinbot_policy_submissions').glob('*.json')]
                if any(c.get('directory') == str(directory) and c != claim
                       and datetime.fromisoformat(c['created_at']) >= submitted for c in receipts):
                    continue
                proof = terminal_evidence(root, claim, repair_pre_admission=True)
                if proof is None:
                    continue
            except (OSError, ValueError, KeyError, TypeError, RuntimeError, StopIteration):
                continue
            # Journal and state precede retirement. On interruption the same
            # old reservation either remains repairable or can retire normally.
            proof.update(kind='legacy_pre_admission_state', claim=claim, before=dict(entry))
            journal = root/'.kinbot_policy_reconciled'/(token+'.json')
            atomic_json(journal, proof)
            entry.update(status='staged', submission_reconciliation=str(journal))
            atomic_json(graph/'state.json', state)
            del claims[token]
            atomic_json(ledger, claims)
    return state
