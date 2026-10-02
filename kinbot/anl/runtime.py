"""Resolve native QC runtime settings without changing the Python process.

Some site installations omit a compiler runtime from their module.  Only the
QC child receives a selected runtime library; adding an entire Conda lib
directory to LD_LIBRARY_PATH can replace unrelated system libraries.

Single-node Molpro jobs do not need a network fabric.  Intel MPI can otherwise
select a broken PSM3/OFI device and abort during MPI initialization before
Molpro reads its input.  Default those jobs to shared-memory transport while
preserving an explicit site or user selection.
"""

from __future__ import annotations

import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


_MISSING = re.compile(r'^\s*(\S+)\s+=>\s+not found\s*$', re.MULTILINE)
_FORTRAN = re.compile(r'libgfortran\.so\.\d+\Z')
_DEFAULT_MOLPRO_SCRATCH_MB = 4096


def _positive_scratch_mb(value):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError('Molpro scratch minimum must be a positive integer.')
    return value


def _molpro_scratch(child, work_directory, minimum_mb):
    """Create an isolated Molpro repository on a sufficiently large filesystem."""
    minimum_mb = _positive_scratch_mb(
        _DEFAULT_MOLPRO_SCRATCH_MB if minimum_mb is None else minimum_mb)
    override = child.get('KINBOT_MOLPRO_SCRATCH')
    if override:
        candidates = [('KINBOT_MOLPRO_SCRATCH', Path(override).expanduser())]
    else:
        candidates = []
        for name in ('SLURM_TMPDIR', 'SCRATCH', 'TMPDIR'):
            value = child.get(name)
            if value:
                candidates.append((name, Path(value).expanduser()))
        home = child.get('HOME')
        if home:
            candidates.append(('HOME', Path(home).expanduser()
                               / '.cache' / 'kinbot' / 'molpro'))
        if work_directory is not None:
            candidates.append(('work_directory', Path(work_directory).resolve()
                               / '.molpro_scratch'))
        candidates.append(('/tmp', Path('/tmp')))

    usable = []
    seen = set()
    errors = []
    for source, root in candidates:
        try:
            root = root.resolve()
        except OSError:
            root = root.absolute()
        if root in seen:
            continue
        seen.add(root)
        try:
            root.mkdir(parents=True, exist_ok=True)
            if not root.is_dir() or not os.access(root, os.W_OK | os.X_OK):
                raise OSError('directory is not writable')
            free_mb = shutil.disk_usage(root).free // (1024 * 1024)
        except OSError as exc:
            errors.append(f'{source}={root} ({exc})')
            continue
        usable.append((source, root, free_mb))
        if free_mb >= minimum_mb:
            scratch = Path(tempfile.mkdtemp(prefix='kinbot-molpro-', dir=root))
            return scratch, {
                'scratch_directory': str(scratch),
                'scratch_root': str(root),
                'scratch_source': source,
                'scratch_available_mb': free_mb,
                'scratch_minimum_mb': minimum_mb,
            }

    capacities = ', '.join(
        f'{source}={root} ({free_mb} MB free)'
        for source, root, free_mb in usable)
    detail = '; '.join(filter(None, (capacities, ', '.join(errors))))
    if override:
        raise RuntimeError('KINBOT_MOLPRO_SCRATCH cannot provide the required '
                           f'{minimum_mb} MB: {detail or override}')
    raise RuntimeError(f'No Molpro scratch filesystem has the required '
                       f'{minimum_mb} MB free: {detail or "no usable candidates"}. '
                       'Set KINBOT_MOLPRO_SCRATCH to a suitable directory.')


def cleanup_qc_runtime(provenance):
    """Remove scratch created by :func:`qc_runtime_environment`."""
    scratch = provenance.get('scratch_directory')
    if scratch:
        shutil.rmtree(scratch, ignore_errors=True)


def _missing_libraries(executable, env):
    """Read ELF dependencies under the same environment as the QC child."""
    with Path(executable).open('rb') as stream:
        if stream.read(4) != b'\x7fELF':
            return []  # Shell launchers and test fixtures are checked by execution.
    if not shutil.which('ldd', path=env.get('PATH')):
        raise RuntimeError('ldd is required to check native QC dependencies.')
    result = subprocess.run(['ldd', str(executable)], env=env,
                            capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError(f'ldd failed for {executable}: '
                           + (result.stderr.strip() or result.stdout.strip()))
    return sorted(set(_MISSING.findall(result.stdout + '\n' + result.stderr)))


def _runtime_candidates(soname, executable, env):
    """Search shallow installation roots; never add their whole lib dirs."""
    override = env.get('CFOUR_LIBGFORTRAN')
    if override:
        path = Path(override).expanduser()
        if path.name != soname or not path.is_file():
            raise RuntimeError(f'CFOUR_LIBGFORTRAN must name an existing {soname}.')
        yield path
        return
    prefixes = [Path(executable).resolve().parent.parent, Path(sys.prefix),
                Path(sys.base_prefix)]
    for key in ('CONDA_PREFIX', 'CONDA_EXE'):
        value = env.get(key)
        if value:
            path = Path(value).expanduser()
            prefixes.append(path.parent.parent if key == 'CONDA_EXE' else path)
    conda = shutil.which('conda', path=env.get('PATH'))
    if conda:
        prefixes.append(Path(conda).resolve().parent.parent)
    seen = set()
    for prefix in prefixes:
        candidate = prefix / 'lib' / soname
        if candidate not in seen:
            seen.add(candidate)
            yield candidate
    for root in (Path('/opt'), Path('/usr/local'), Path.home()):
        for candidate in sorted(root.glob(f'*/lib/{soname}')):
            if candidate not in seen:
                seen.add(candidate)
                yield candidate


def qc_runtime_environment(program, backend, env=None, *, work_directory=None,
                           scratch_min_mb=None):
    """Return a child-only environment and selected runtime provenance.

    A loaded module or normal system runtime wins.  For CFOUR installations
    missing libgfortran, test one matching SONAME with ldd before launching.
    """
    child = dict(os.environ if env is None else env)
    backend = backend.lower()
    if backend == 'molpro':
        nodes = child.get('SLURM_NNODES', child.get('SLURM_JOB_NUM_NODES'))
        if nodes == '1' and not child.get('I_MPI_FABRICS'):
            child['I_MPI_FABRICS'] = 'shm'
        provenance = ({'mpi_fabrics': child['I_MPI_FABRICS']}
                      if child.get('I_MPI_FABRICS') else {})
        if work_directory is not None:
            scratch, scratch_provenance = _molpro_scratch(
                child, work_directory, scratch_min_mb)
            child['TMPDIR'] = str(scratch)
            provenance.update(scratch_provenance)
        return child, provenance
    if backend != 'cfour':
        return child, {}
    executable = shutil.which(program, path=child.get('PATH'))
    if not executable:
        raise RuntimeError(f'CFOUR executable unavailable: {program}.')
    missing = _missing_libraries(executable, child)
    if not missing:
        return child, {}
    if len(missing) != 1 or not _FORTRAN.fullmatch(missing[0]):
        raise RuntimeError('CFOUR executable has unresolved shared libraries: '
                           + ', '.join(missing))
    soname = missing[0]
    for candidate in _runtime_candidates(soname, executable, child):
        if not candidate.is_file():
            continue
        trial = dict(child)
        trial['LD_PRELOAD'] = ' '.join(filter(None, (
            str(candidate), child.get('LD_PRELOAD', ''))))
        if not _missing_libraries(executable, trial):
            return trial, {'runtime_library': str(candidate.resolve())}
    raise RuntimeError(f'CFOUR needs {soname}, but no usable copy was found. '
                       'Load its compiler runtime module or set '
                       f'CFOUR_LIBGFORTRAN to an absolute {soname} path.')
