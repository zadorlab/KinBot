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

import errno
import json
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import sys
import tempfile

from kinbot.execution import atomic_json
from kinbot.slurm_retry import io_failure as molpro_io_failure


_MISSING = re.compile(r'^\s*(\S+)\s+=>\s+not found\s*$', re.MULTILINE)
_FORTRAN = re.compile(r'libgfortran\.so\.\d+\Z')
_DEFAULT_MOLPRO_SCRATCH_MB = 4096


class MolproScratchError(RuntimeError):
    """Scratch setup failed, with separately classified transient I/O evidence."""

    def __init__(self, message, *, retryable_io=False):
        super().__init__(message)
        self.retryable_io = retryable_io


def _positive_scratch_mb(value):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError('Molpro scratch minimum must be a positive integer.')
    return value


def _check_scratch_io(directory):
    """Exercise the nested directories and file I/O needed by a QC launcher."""
    probe = Path(directory) / '.kinbot-io-check'
    probe.mkdir()
    payload = b'KinBot Molpro scratch check\n' * 128
    original = probe / 'write'
    renamed = probe / 'read'
    with original.open('xb') as stream:
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    original.rename(renamed)
    if renamed.read_bytes() != payload:
        raise OSError(errno.EIO, 'scratch read-back differs from written data')
    renamed.unlink()
    probe.rmdir()


def _molpro_scratch(child, work_directory, minimum_mb, *, avoid_roots=()):
    """Create and probe an isolated repository on a sufficiently large filesystem.

    Prefer the executing account's /scratch directory. Every alternative
    must also be verified local storage; home/work-directory fallbacks are
    forbidden. No filesystem probe guarantees subsequent native I/O.
    """
    minimum_mb = _positive_scratch_mb(
        _DEFAULT_MOLPRO_SCRATCH_MB if minimum_mb is None else minimum_mb)
    configured_minimum = child.get('KINBOT_MOLPRO_SCRATCH_MIN_MB')
    if configured_minimum is not None:
        if not re.fullmatch(r'[0-9]+', str(configured_minimum)) or int(configured_minimum) < 1:
            raise ValueError('KINBOT_MOLPRO_SCRATCH_MIN_MB must be a positive integer.')
        # The environment may raise a task's capacity requirement, never
        # weaken an explicit resource request recorded with that task.
        minimum_mb = max(minimum_mb, int(configured_minimum))
    from kinbot.qc_scratch import allocate
    try:
        return allocate(child, 'molpro', minimum_mb, probe=_check_scratch_io,
                        work_directory=work_directory, avoid_roots=avoid_roots)
    except OSError as error:
        raise MolproScratchError(str(error), retryable_io=error.errno in
                                 (errno.EIO, errno.ESTALE, errno.ETIMEDOUT)) from error


def cleanup_qc_runtime(provenance):
    """Remove scratch created by :func:`qc_runtime_environment`."""
    scratch = provenance.get('scratch_directory')
    if scratch and not provenance.get('retain_scratch'):
        shutil.rmtree(scratch, ignore_errors=True)


def preserve_qc_failure(provenance, directory, label, returncode):
    """Keep bounded native/rank diagnostics before removing private scratch.

    Large integrals/checkpoints are not copied. MPI abort alone is never
    labelled as an out-of-memory, disk, or electronic-convergence failure.
    Failure of this best-effort capture must not hide the native exit status.
    """
    try:
        return _preserve_qc_failure(provenance, directory, label, returncode)
    except OSError:
        return None


def _preserve_qc_failure(provenance, directory, label, returncode):
    directory = Path(directory)
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', label):
        raise ValueError('Failure report label must be a basename.')
    report = {'returncode': returncode, 'runtime': dict(provenance),
              'host': socket.gethostname(), 'files': {}, 'capture_errors': []}
    root = Path(provenance['scratch_directory']) if provenance.get('scratch_directory') else None
    native_output = re.compile(r'\.(?:out|log|err|stderr|stdout)(?:[._]\d+)?$', re.I)
    candidates = [(f'work/{p.name}', p) for p in directory.glob(label + '.*')
                  if native_output.search(p.name)]
    if root and root.is_dir():
        try:
            report['scratch_free_bytes_at_exit'] = shutil.disk_usage(root).free
            for parent, folders, files in os.walk(root, followlinks=False):
                folders[:] = [n for n in folders if not Path(parent, n).is_symlink()]
                for name in files:
                    if native_output.search(name):
                        p = Path(parent, name)
                        candidates.append(('scratch/' + str(p.relative_to(root)), p))
                        if len(candidates) >= 64:
                            break
                if len(candidates) >= 64:
                    report['file_limit_reached'] = True
                    break
        except OSError as error:
            report['capture_errors'].append(str(error))
    context_budget = 64*1024*1024
    for name, path in candidates[:64]:
        if path.is_symlink():
            continue
        try:
            from kinbot.native_failure_report import excerpt
            # Bound aggregate reads too; every remaining file still keeps its
            # ordinary tail once the larger error-context budget is used.
            limit = max(65536, min(16*1024*1024, context_budget))
            evidence = excerpt(path, limit=limit)
            report['files'][name] = evidence
            context_budget -= min(evidence['size_bytes'], limit)
        except OSError as error:
            report['capture_errors'].append(f'{name}: {error}')
    report['explicit_io_failure'] = any(
        molpro_io_failure(e['tail']+'\n'+e.get('error_context',''))
        for e in report['files'].values())
    destination = directory / f'{label}.failure.json'
    try:
        atomic_json(destination, report)
    except OSError:
        return None
    return destination


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
                           scratch_min_mb=None, avoid_scratch_roots=()):
    """Return a child-only environment and selected runtime provenance.

    A loaded module or normal system runtime wins.  For CFOUR installations
    missing libgfortran, test one matching SONAME with ldd before launching.
    """
    child = dict(os.environ if env is None else env)
    backend = backend.lower()
    # Resolve libraries before allocating a workspace so a failed loader
    # preflight cannot leak a scratch directory.
    if backend == 'cfour' and work_directory is not None:
        child, libraries = qc_runtime_environment(program, backend, child)
        from kinbot.qc_scratch import allocate
        scratch, provenance = allocate(child, backend, scratch_min_mb or 4096,
                                       probe=_check_scratch_io,
                                       work_directory=work_directory)
        child['TMPDIR'] = str(scratch)
        return child, {**libraries, **provenance}
    if backend == 'molpro':
        nodes = child.get('SLURM_NNODES', child.get('SLURM_JOB_NUM_NODES'))
        if nodes == '1' and not child.get('I_MPI_FABRICS'):
            child['I_MPI_FABRICS'] = 'shm'
        provenance = ({'mpi_fabrics': child['I_MPI_FABRICS']}
                      if child.get('I_MPI_FABRICS') else {})
        if work_directory is not None:
            scratch, scratch_provenance = _molpro_scratch(
                child, work_directory, scratch_min_mb, avoid_roots=avoid_scratch_roots)
            child['TMPDIR'] = str(scratch)
            provenance.update(scratch_provenance)
            attempt = child.get('KINBOT_SLURM_IO_ATTEMPT_DIRECTORY')
            if attempt:
                atomic_json(Path(attempt) / f'scratch-{scratch.name}.json', scratch_provenance)
        return child, provenance
    provenance = {}
    if work_directory is not None:
        from kinbot.qc_scratch import allocate
        scratch, provenance = allocate(child, backend, scratch_min_mb or 4096,
                                       probe=_check_scratch_io,
                                       work_directory=work_directory)
        child['TMPDIR'] = str(scratch)
        if backend in ('gaussian', 'gauss'):
            child['GAUSS_SCRDIR'] = str(scratch)
    if backend != 'cfour':
        return child, provenance
    executable = shutil.which(program, path=child.get('PATH'))
    if not executable:
        raise RuntimeError(f'CFOUR executable unavailable: {program}.')
    missing = _missing_libraries(executable, child)
    if not missing:
        return child, provenance
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
            return trial, {**provenance, 'runtime_library': str(candidate.resolve())}
    raise RuntimeError(f'CFOUR needs {soname}, but no usable copy was found. '
                       'Load its compiler runtime module or set '
                       f'CFOUR_LIBGFORTRAN to an absolute {soname} path.')
