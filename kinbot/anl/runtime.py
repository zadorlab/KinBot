"""Resolve native QC program dependencies without changing the Python process.

Some site installations omit a compiler runtime from their module.  Only the
QC child receives a selected runtime library; adding an entire Conda lib
directory to LD_LIBRARY_PATH can replace unrelated system libraries.
"""

from __future__ import annotations

import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


_MISSING = re.compile(r'^\s*(\S+)\s+=>\s+not found\s*$', re.MULTILINE)
_FORTRAN = re.compile(r'libgfortran\.so\.\d+\Z')


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


def qc_runtime_environment(program, backend, env=None):
    """Return a child-only environment and selected runtime provenance.

    A loaded module or normal system runtime wins.  For CFOUR installations
    missing libgfortran, test one matching SONAME with ldd before launching.
    """
    child = dict(os.environ if env is None else env)
    if backend.lower() != 'cfour':
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
