"""Run generated rotdPy inputs and retain a restartable execution record."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def result_path(input_file: str | Path) -> Path:
    input_file = Path(input_file)
    return input_file.with_name(f'{input_file.stem}.rotdpy.json')


def execution_path(input_file: str | Path) -> Path:
    input_file = Path(input_file)
    return input_file.with_name(f'{input_file.stem}.execution.json')


def ensure_available() -> None:
    """Fail before QC work when an explicitly requested rotdPy is absent."""
    if importlib.util.find_spec('rotd_py') is None:
        raise RuntimeError(
            "rotdPy is required because 'rotdpy_run' is enabled. Install the "
            "rotdPy source distribution into the same environment as KinBot "
            "and verify that 'import rotd_py' succeeds.")


def read_result(input_file: str | Path) -> dict:
    """Read and validate the manifest written by a completed rotdPy input."""
    path = result_path(input_file)
    if not path.is_file():
        raise RuntimeError(f'rotdPy did not write {path.name}.')
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise RuntimeError(f'Cannot read rotdPy result {path}: {error}') from error
    if payload.get('status') != 'complete':
        raise RuntimeError(f'rotdPy result {path} is not complete.')
    if (not isinstance(payload.get('surface_count'), int)
            or payload['surface_count'] < 1):
        raise RuntimeError(f'rotdPy result {path} has no sampled surfaces.')
    if not isinstance(payload.get('result_files'), list):
        raise RuntimeError(f'rotdPy result {path} has no result file list.')
    return payload


def run(input_file: str | Path, python: str | Path | None = None) -> dict:
    """Execute one generated input with the active KinBot interpreter.

    A successful record is reused only while the generated input hash and
    result manifest both remain valid.  rotdPy's own database therefore stays
    available for its normal restart behavior after interrupted sampling.
    """
    input_file = Path(input_file).resolve()
    if not input_file.is_file():
        raise FileNotFoundError(input_file)
    ensure_available()

    input_hash = _sha256(input_file)
    record_file = execution_path(input_file)
    if record_file.is_file():
        try:
            previous = json.loads(record_file.read_text())
            if (previous.get('status') == 'complete'
                    and previous.get('input_sha256') == input_hash):
                read_result(input_file)
                return previous
        except (OSError, json.JSONDecodeError, RuntimeError):
            pass

    result_file = result_path(input_file)
    result_file.unlink(missing_ok=True)
    stdout_file = input_file.with_name(f'{input_file.stem}.rotdpy.stdout')
    stderr_file = input_file.with_name(f'{input_file.stem}.rotdpy.stderr')
    command = [str(python or sys.executable), input_file.name]
    completed = subprocess.run(
        command, cwd=input_file.parent, text=True, capture_output=True,
        check=False)
    stdout_file.write_text(completed.stdout)
    stderr_file.write_text(completed.stderr)

    record = {
        'schema': 1,
        'status': 'failed',
        'input': input_file.name,
        'input_sha256': input_hash,
        'command': command,
        'returncode': completed.returncode,
        'stdout': stdout_file.name,
        'stderr': stderr_file.name,
    }
    error = None
    if completed.returncode:
        error = f'rotdPy exited with status {completed.returncode}.'
    else:
        try:
            result = read_result(input_file)
            record['result'] = result_file.name
            record['surface_count'] = result['surface_count']
            record['status'] = 'complete'
        except RuntimeError as caught:
            error = str(caught)
    if error:
        record['error'] = error
    record_file.write_text(json.dumps(record, indent=2) + '\n')
    if error:
        raise RuntimeError(f'{input_file.stem}: {error}')
    return record
