"""Shared, dependency-free helpers for restartable execution records."""

import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile


SLURM_TERMINAL_STATES = {
    'BOOT_FAIL', 'CANCELLED', 'COMPLETED', 'DEADLINE', 'FAILED',
    'NODE_FAIL', 'OUT_OF_MEMORY', 'PREEMPTED', 'TIMEOUT',
}


def slurm_job_record(job_id, *, timeout=15, strict=False):
    """Retained scheduler evidence; strict mode distinguishes outage from purge."""
    try:
        result = subprocess.run(
            ['scontrol', 'show', 'job', '-o', str(job_id)],
            capture_output=True, text=True, check=False, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as error:
        if strict:
            raise OSError(f'Cannot query Slurm job {job_id}: {error}') from error
        return {}
    if result.returncode or not result.stdout.strip():
        if strict and not (result.returncode and not result.stdout.strip() and re.fullmatch(
                r'(?:slurm_load_jobs error:\s*)?Invalid job id specified',
                result.stderr.strip(), re.I)):
            raise OSError(f'Slurm job query was not observed: {result.stderr.strip()}')
        return {}
    return {match.group(1): match.group(2)
            for match in re.finditer(r'(\w+)=([^\s]+)', result.stdout)}


_SCHEDULER_UNAVAILABLE = re.compile(
    r'(?:unable to contact (?:slurm )?controller|'
    r'unable to establish controller connection|'
    r'(?:scheduler )?controller unavailable|'
    r'slurm_load_jobs error|socket timed out|connection refused|'
    r'resource temporarily unavailable|network is unreachable|'
    r'host is down|no route to host|'
    r'communication connection failure|controller not responding|'
    r'scheduler query timed out)', re.IGNORECASE)


def scheduler_unavailable(detail):
    """Return whether a diagnostic means the scheduler did not answer."""
    return bool(_SCHEDULER_UNAVAILABLE.search(str(detail or '')))


def file_sha256(path):
    """Hash a file without loading large native outputs into memory."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, payload):
    """Publish canonical JSON through an atomic same-directory rename."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(
        payload, indent=2, sort_keys=True) + '\n')
    os.replace(temporary, path)
    return path


def atomic_text(path, text):
    """Publish shared text without truncating files read by other owners."""
    path = Path(path)
    if path.is_file() and path.read_text() == text:
        return path
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as f:
        temporary = Path(f.name)
        f.write(text)
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return path


def resolve_path(base, value, label):
    """Resolve a required path relative to a manifest's directory."""
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{label} must be a nonempty path.')
    path = Path(value).expanduser()
    return (path if path.is_absolute() else Path(base) / path).resolve()
