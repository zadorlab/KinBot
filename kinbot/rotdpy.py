"""Run generated rotdPy inputs and retain a restartable execution record."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
from importlib import metadata
import json
import pickle
from pathlib import Path
import re
import sqlite3
import subprocess
import sys
import time


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


def sampling_progress(input_file: str | Path) -> dict:
    """Report live ROTD_py surface progress without opening its pickle DB.

    A surface file is written atomically by ROTD_py when that dividing surface
    converges.  Counting those files is safer during a live run than
    unpickling the restart database while its driver is updating it.  The ETA
    is deliberately labelled as linear because adaptive Monte Carlo surfaces
    can require very different numbers of samples.
    """
    input_file = Path(input_file).resolve()
    if not input_file.is_file():
        raise FileNotFoundError(input_file)
    root = input_file.with_name(f'kb_{input_file.stem}')
    surfaces = sorted(
        (path for path in root.glob('Surface_*') if path.is_dir()),
        key=lambda path: int(path.name.removeprefix('Surface_')))
    completed = sorted(
        root.glob('output/surface_*.dat'),
        key=lambda path: int(path.stem.removeprefix('surface_')))
    active = []
    pattern = re.compile(r'surf(\d+)_face(\d+)_samp(\d+)\.pkl')
    for path in root.glob('Surface_*/jobs/*.pkl'):
        match = pattern.fullmatch(path.name)
        if match:
            active.append({
                'surface': int(match.group(1)),
                'face': int(match.group(2)),
                'sample': int(match.group(3)),
            })
    active.sort(key=lambda item: (
        item['surface'], item['face'], item['sample']))
    total = len(surfaces)
    done = len(completed)
    elapsed = max(0., time.time() - input_file.stat().st_mtime)
    eta = (elapsed * (total - done) / done
           if 0 < done < total else 0. if done == total and total else None)
    database = root / 'rotdPy_restart.db'
    snapshot = None
    if database.is_file():
        try:
            with sqlite3.connect(
                    f'file:{database.resolve()}?mode=ro', uri=True,
                    timeout=10) as connection:
                rows = connection.execute(
                    'SELECT surf_id, multi_flux, run_index '
                    'FROM rotdpy_saved_runs').fetchall()
            latest = {}
            for surface, blob, run_index in rows:
                surface = int(surface)
                if (surface not in latest
                        or run_index >= latest[surface][1]):
                    latest[surface] = (blob, run_index)
            accepted_by_surface = {}
            failed = 0
            space = 0
            ceiling = 0
            database_converged = 0
            for surface, (blob, _) in latest.items():
                flux = pickle.loads(blob)
                faces = list(flux.flux_array)
                accepted = sum(int(getattr(face, '_acct_num', 0))
                               for face in faces)
                accepted_by_surface[str(surface)] = accepted
                failed += sum(int(getattr(face, '_fail_num', 0))
                              for face in faces)
                space += sum(int(getattr(face, '_close_num', 0))
                             + int(getattr(face, '_face_num', 0))
                             for face in faces)
                selected = list(getattr(
                    flux, 'selected_faces', range(len(faces))))
                ceiling += int(getattr(flux, 'pot_max', 0)) * len(selected)
                database_converged += bool(getattr(flux, 'converged', False))
            accepted_values = list(accepted_by_surface.values())
            accepted_total = sum(accepted_values)
            snapshot = {
                'saved_at_epoch': round(database.stat().st_mtime),
                'surfaces_in_snapshot': len(latest),
                'database_converged_surfaces': database_converged,
                'accepted_samples': accepted_total,
                'accepted_samples_per_surface_min': (
                    min(accepted_values) if accepted_values else 0),
                'accepted_samples_per_surface_max': (
                    max(accepted_values) if accepted_values else 0),
                'failed_samples': failed,
                'space_rejections': space,
                'potential_sample_ceiling': ceiling,
                'accepted_samples_to_ceiling': max(
                    0, ceiling - accepted_total),
                'note': ('The restart snapshot is periodic and may lag '
                         'currently running sample jobs.'),
            }
        except (OSError, sqlite3.Error, pickle.UnpicklingError,
                AttributeError, EOFError, ImportError) as error:
            snapshot = {'error': str(error)}
    manifest = result_path(input_file)
    result = {
        'status': ('complete' if manifest.is_file()
                   else 'running' if root.is_dir() else 'not_started'),
        'input': str(input_file),
        'sample_root': str(root),
        'surface_count': total,
        'converged_surfaces': done,
        'remaining_surfaces': max(0, total - done),
        'convergence_fraction': done / total if total else 0.,
        'active_samples': active,
        'elapsed_seconds': round(elapsed),
        'linear_eta_seconds': round(eta) if eta is not None else None,
        'eta_note': ('Linear estimate from completed surfaces; adaptive Monte '
                     'Carlo convergence can vary substantially by surface.'),
    }
    if snapshot is not None:
        result['restart_snapshot'] = snapshot
    return result


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


def _has_finite_numeric_data(path: Path) -> bool:
    """Return whether a text output contains at least one finite number."""
    try:
        text = path.read_text(errors='replace')
    except OSError:
        return False
    for token in re.findall(
            r'(?<![A-Za-z])[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?',
            text):
        try:
            value = float(token)
        except ValueError:
            continue
        if value == value and abs(value) != float('inf'):
            return True
    return False


def _validate_production_result(input_file: Path, payload: dict,
                                outputs: list[Path]) -> None:
    """Reject reduced interface samples from production MESS calculations."""
    validation = payload.get('validation')
    if not isinstance(validation, dict) or validation.get('name') != 'production':
        raise RuntimeError('rotdPy result is not production validated.')
    correction = validation.get('correction')
    if (not isinstance(correction, dict)
            or correction.get('kind') != 'one_dimensional'
            or not isinstance(correction.get('point_count'), int)
            or correction['point_count'] < 4
            or not re.fullmatch(r'[0-9a-f]{64}', str(
                correction.get('source_sha256', '')))):
        raise RuntimeError('Production rotdPy result lacks a multipoint, '
                           'hash-identified correction potential.')
    surfaces = validation.get('dividing_surfaces')
    if (not isinstance(surfaces, dict)
            or not isinstance(surfaces.get('requested_distances_angstrom'), list)
            or len(surfaces['requested_distances_angstrom']) < 3
            or surfaces.get('generated_count') != payload['surface_count']
            or payload['surface_count'] < 3):
        raise RuntimeError('Production rotdPy result lacks multiple verified '
                           'dividing surfaces.')
    levels = (validation.get('sampling_level'),
              validation.get('trusted_correction_level'))
    if any(not isinstance(level, dict)
           or not str(level.get('method', '')).strip()
           or not str(level.get('basis', '')).strip() for level in levels):
        raise RuntimeError('Production rotdPy method provenance is incomplete.')
    statistics = payload.get('surface_statistics')
    if not isinstance(statistics, dict) or len(statistics) != payload['surface_count']:
        raise RuntimeError('Production rotdPy surface statistics are incomplete.')
    for statistic in statistics.values():
        if (not isinstance(statistic, dict)
                or statistic.get('converged') is not True
                or not isinstance(statistic.get('accepted_samples'), int)
                or statistic['accepted_samples'] < 1):
            raise RuntimeError('Production rotdPy Monte Carlo sampling is '
                               'not converged on every surface.')
    if payload.get('input_sha256') != _sha256(input_file):
        raise RuntimeError('Production rotdPy input hash changed.')
    if any(path.stat().st_size == 0 or not _has_finite_numeric_data(path)
           for path in outputs):
        raise RuntimeError('Production rotdPy output lacks numeric data.')


def read_result(input_file: str | Path,
                required_profile: str | None = None) -> dict:
    """Read and validate the manifest written by a completed rotdPy input."""
    input_file = Path(input_file).resolve()
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
    manifest_profile = (payload.get('validation') or {}).get('name')
    if manifest_profile == 'production' or required_profile == 'production':
        _validate_production_result(input_file, payload, outputs)
    elif required_profile not in (None, 'interface'):
        raise ValueError('required ROTD_py profile must be interface or '
                         'production.')
    return payload


def number_of_states_file(input_file: str | Path, energy_index: int = -1,
                          required_profile: str | None = None
                          ) -> tuple[Path, dict]:
    """Return one hash-verified MESS ``Rotd`` number-of-states file.

    ROTD_py numbers its electronic/correction surfaces as ``Ne_0.out``,
    ``Ne_1.out``, and so on.  Index ``-1`` selects the highest available
    index, which is the final corrected surface written by ROTD_py.  An
    explicit nonnegative index can be used when a different surface is
    scientifically intended.  The returned metadata is suitable for a MESS
    provenance sidecar.
    """
    if isinstance(energy_index, bool) or not isinstance(energy_index, int) \
            or energy_index < -1:
        raise ValueError(
            'ROTD_py MESS energy index must be -1 or a nonnegative integer.')
    input_file = Path(input_file).resolve()
    payload = read_result(input_file, required_profile=required_profile)
    candidates = {}
    for relative in payload['result_files']:
        name = Path(relative).name
        match = re.fullmatch(r'Ne_(\d+)\.out', name)
        if match:
            candidates[int(match.group(1))] = relative
    if not candidates:
        raise RuntimeError(
            f'rotdPy result for {input_file.stem} has no MESS '
            'number-of-states file.')
    selected = max(candidates) if energy_index == -1 else energy_index
    if selected not in candidates:
        raise RuntimeError(
            f'rotdPy result for {input_file.stem} has no Ne_{selected}.out; '
            f'available indices are {sorted(candidates)}.')
    relative = candidates[selected]
    source = input_file.parent / relative
    return source.resolve(), {
        'schema': 1,
        'reaction': payload.get('reaction', input_file.stem),
        'surface_count': payload['surface_count'],
        'validation_profile': (payload.get('validation') or {}).get(
            'name', 'interface'),
        'energy_index': selected,
        'selection': ('highest_available_correction'
                      if energy_index == -1 else 'explicit'),
        'source': str(source.resolve()),
        'source_sha256': payload['result_sha256'][relative],
        'result_manifest': str(result_path(input_file).resolve()),
        'result_manifest_sha256': _sha256(result_path(input_file)),
    }


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


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Run and verify generated KinBot ROTD_py inputs')
    commands = parser.add_subparsers(dest='action', required=True)
    run_parser = commands.add_parser('run')
    run_parser.add_argument('input', type=Path)
    check_parser = commands.add_parser('check')
    check_parser.add_argument('input', type=Path)
    check_parser.add_argument('--profile', choices=('interface', 'production'))
    select_parser = commands.add_parser('select')
    select_parser.add_argument('input', type=Path)
    select_parser.add_argument('--energy-index', type=int, default=-1)
    select_parser.add_argument('--profile', choices=('interface', 'production'))
    progress_parser = commands.add_parser('progress')
    progress_parser.add_argument('input', type=Path)
    args = parser.parse_args(argv)
    if args.action == 'run':
        result = run(args.input)
    elif args.action == 'check':
        result = read_result(args.input, required_profile=args.profile)
    elif args.action == 'select':
        source, provenance = number_of_states_file(
            args.input, args.energy_index, required_profile=args.profile)
        result = {'source': str(source), **provenance}
    else:
        result = sampling_progress(args.input)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
