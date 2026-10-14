"""Shared, dependency-free helpers for restartable execution records."""

import hashlib
import json
import os
from pathlib import Path
import re


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


def resolve_path(base, value, label):
    """Resolve a required path relative to a manifest's directory."""
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{label} must be a nonempty path.')
    path = Path(value).expanduser()
    return (path if path.is_absolute() else Path(base) / path).resolve()
