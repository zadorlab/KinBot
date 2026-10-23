"""Evidence for explicit owner retries of jobs lost before any worker output.

Slurm may purge terminal records while accounting is disabled. A best-effort
empty lookup is not proof: require an explicit invalid-ID response, a healthy
queue query, the saved submission receipt, unchanged staged files, and no
evidence that a worker started. This does not classify the original exit.
"""
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess
import time

from kinbot.execution import file_sha256


def evidence(directory, job, task, entry):
    """Return auditable evidence only for an explicit, generation-fenced retry."""
    from kinbot.anl import dispatch
    from kinbot.run_policy import find

    directory, job = Path(directory).resolve(), str(job)
    if (not job.isdecimal() or str(entry.get('job_id')) != job
            or entry.get('status') not in ('failed', 'unavailable')
            or entry.get('failure_phase') != 'scheduler'
            or entry.get('error') != 'Slurm job disappeared without execution.json.'):
        return None
    missing = entry.get('missing_since')
    if (not isinstance(missing, (int, float))
            or not dispatch._RESULT_PUBLICATION_GRACE < time.time() - missing):
        return None
    # Native artifacts, even empty scheduler logs or a guard journal, require
    # their own diagnosis. Do not erase a partially executed calculation.
    staged = {'task.json', 'geometry.xyz', 'job.slurm', task.get('input_name')}
    if any(p.name not in staged or not p.is_file() for p in directory.iterdir()):
        return None
    policy = find(directory)
    if policy is None:
        return None
    receipt = policy.parent / '.kinbot_policy_submissions' / (job + '.json')
    if not receipt.is_file():
        return None
    saved = json.loads(receipt.read_text())
    if (str(saved.get('job_id')) != job or saved.get('status') != 'submitted'
            or not saved.get('directory')
            or Path(saved['directory']).resolve() != directory
            or not saved.get('command') or saved['command'][0] != 'sbatch'):
        return None
    try:
        submitted = datetime.fromisoformat(saved['created_at']).timestamp()
    except (KeyError, TypeError, ValueError):
        return None
    if submitted > missing:
        return None
    dispatch._verify_stage_files(directory.parent.parent, task, entry)
    try:
        query = subprocess.run(['scontrol', 'show', 'job', '-o', job],
                               capture_output=True, text=True, timeout=15)
        if (not query.returncode or query.stdout.strip() or not re.fullmatch(
                r'(?:slurm_load_jobs error:\s*)?Invalid job id specified',
                query.stderr.strip(), re.I)):
            return None
        queue = subprocess.run(['squeue', '-h', '-o', '%i'],
                               capture_output=True, text=True, timeout=15)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if queue.returncode or queue.stderr.strip() or job in queue.stdout.split():
        return None
    return dict(kind='purged_submission_without_worker_output', job_id=job,
                checked_at=datetime.now(timezone.utc).isoformat(),
                admission=str(receipt), admission_sha256=file_sha256(receipt),
                task_sha256=entry['task_sha256'], missing_since=missing,
                scheduler_diagnostic=query.stderr.strip(),
                original_exit='unknown', queue_absent=True)
