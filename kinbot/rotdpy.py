"""Run generated rotdPy inputs and retain a restartable execution record."""

from __future__ import annotations

import hashlib
import importlib.util
from importlib import metadata
import json
from pathlib import Path
import subprocess
import sys


SUPPORTED_REVISION = '245317ca46c3324da4f71a80ff7263cbe7ceeac1'


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


def dependency_provenance() -> dict:
    """Return the installed ROTD_py version, location, and Git revision."""
    spec = importlib.util.find_spec('rotd_py')
    if spec is None or not spec.submodule_search_locations:
        raise RuntimeError(
            "rotdPy is required because 'rotdpy_run' is enabled. Install the "
            "pinned external/ROTD_py submodule into the same environment as "
            "KinBot and verify that 'import rotd_py' succeeds.")
    location = Path(next(iter(spec.submodule_search_locations))).resolve()
    try:
        version = metadata.version('rotd-py')
    except metadata.PackageNotFoundError:
        version = None
    checkout = location.parent
    revision = None
    try:
        result = subprocess.run(
            ['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True,
            capture_output=True, check=False)
        if result.returncode == 0:
            revision = result.stdout.strip()
    except OSError:
        pass
    return {'version': version, 'location': str(location),
            'git_revision': revision}


def ensure_available() -> dict:
    """Fail before QC work when the pinned rotdPy runtime is incomplete."""
    provenance = dependency_provenance()
    if provenance['git_revision'] != SUPPORTED_REVISION:
        raise RuntimeError(
            'Installed rotdPy revision does not match KinBot\'s pinned '
            f'external/ROTD_py revision {SUPPORTED_REVISION}.')
    try:
        from rotd_py.flux.fluxbase import FluxBase  # noqa: F401
        from rotd_py.new_multi import Multi  # noqa: F401
        from rotd_py.sample.multi_sample import MultiSample  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(f'rotdPy dependency import failed: {exc}') from exc
    return provenance


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
    hashes = payload.get('result_sha256')
    if not isinstance(hashes, dict):
        raise RuntimeError(f'rotdPy result {path} has no output hashes.')
    root = path.parent.resolve()
    outputs = []
    for relative in payload['result_files']:
        if not isinstance(relative, str):
            raise RuntimeError(f'rotdPy result {path} has an invalid file name.')
        output = (root / relative).resolve()
        if root not in output.parents or not output.is_file():
            raise RuntimeError(f'rotdPy result file is missing or unsafe: {relative}.')
        if hashes.get(relative) != _sha256(output):
            raise RuntimeError(f'rotdPy result file hash changed: {relative}.')
        outputs.append(output)
    if not any(output.name.startswith('Ne_') and output.suffix == '.out'
               for output in outputs):
        raise RuntimeError(f'rotdPy result {path} has no number-of-states output.')
    if not any(output.parent.name == 'output'
               and output.name.startswith('surface_')
               and output.suffix == '.dat' for output in outputs):
        raise RuntimeError(f'rotdPy result {path} has no surface flux output.')
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
    provenance = ensure_available()
    if not input_file.with_name('qu.tpl').is_file():
        raise RuntimeError(f'{input_file.stem}: missing rotdPy Slurm qu.tpl.')

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
        'rotdpy': provenance,
    }
    error = None
    if completed.returncode:
        detail = (completed.stderr.strip() or completed.stdout.strip())[-2000:]
        suffix = f' Last output: {detail}' if detail else ''
        error = f'rotdPy exited with status {completed.returncode}.{suffix}'
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
