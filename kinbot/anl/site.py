"""Discover portable site setup for prepared ANL batch jobs.

Discovery runs on the login node during ``prepare``.  The resulting shell
script pins the installations that were visible there, so a batch shell does
not need to inherit the login shell's PATH or loaded modules.
"""

from __future__ import annotations

import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess


_MODULE_NAME = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.+/-]*\Z')
_PARTITION_NAME = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*\Z')


def _duration_seconds(value):
    """Parse Slurm's day-hour:minute:second or minute:second display."""
    if value.lower() in ('infinite', 'unlimited'):
        return float('inf')
    days = 0
    if '-' in value:
        day_text, value = value.split('-', 1)
        days = int(day_text)
    parts = [int(part) for part in value.split(':')]
    if len(parts) == 3:
        hours, minutes, seconds = parts
    elif len(parts) == 2:
        hours, minutes, seconds = 0, *parts
    elif len(parts) == 1:
        hours, minutes, seconds = 0, parts[0], 0
    else:
        raise ValueError(f'Invalid Slurm duration: {value!r}')
    return ((days * 24 + hours) * 60 + minutes) * 60 + seconds


def _partitions():
    """Return available partition/node groups; no Slurm means no discovery."""
    if not shutil.which('sinfo'):
        return []
    result = subprocess.run(
        ['sinfo', '-h', '-o', '%P|%c|%m|%l|%a'],
        capture_output=True, text=True, check=False, timeout=15)
    if result.returncode:
        raise RuntimeError('sinfo failed during partition discovery: '
                           + (result.stderr.strip() or result.stdout.strip()))
    entries = []
    for line in result.stdout.splitlines():
        fields = [field.strip() for field in line.split('|')]
        if len(fields) != 5 or fields[4].lower() != 'up':
            continue
        name = fields[0].removesuffix('*')
        if not _PARTITION_NAME.fullmatch(name):
            continue
        try:
            entries.append({
                'name': name,
                'default': fields[0].endswith('*'),
                'cores': int(fields[1].removesuffix('+')),
                'memory_mb': int(fields[2].removesuffix('+')),
                'seconds': _duration_seconds(fields[3]),
            })
        except ValueError:
            continue
    return entries


def assign_partitions(spec):
    """Fill missing Slurm partitions with the shortest fitting available one.

    Explicit task partitions always win.  When Slurm is absent (for example
    on a development laptop), leave them unset for the target site to choose.
    """
    pending = [task for task in spec['tasks']
               if not task['resources'].get('partition')]
    if not pending:
        return
    entries = _partitions()
    if not entries:
        if shutil.which('sinfo'):
            raise RuntimeError('sinfo returned no usable up partitions; set '
                               'resources.partition explicitly for each task.')
        return
    for task in pending:
        resources = task['resources']
        needed = _duration_seconds(resources['walltime'])
        fits = [entry for entry in entries
                if entry['cores'] >= resources['cores']
                and entry['memory_mb'] >= resources['memory_mb']
                and entry['seconds'] >= needed]
        if not fits:
            raise RuntimeError(f"{task['id']}: no available Slurm partition "
                               'fits requested cores, memory, and walltime; '
                               'adjust resources or set a partition explicitly.')
        chosen = min(fits, key=lambda item: (
            item['seconds'], not item['default'], item['name']))
        resources['partition'] = chosen['name']


def _loaded_module(backend):
    aliases = {'gauss': 'gaussian'}
    target = aliases.get(backend, backend)
    matches = [name for name in os.environ.get('LOADEDMODULES', '').split(':')
               if name and name.split('/', 1)[0].lower() == target]
    if len(matches) > 1:
        raise RuntimeError(f'{backend}: multiple loaded modules match; '
                           'unload the unwanted version before prepare.')
    if matches and not _MODULE_NAME.fullmatch(matches[0]):
        raise RuntimeError(f'{backend}: unsafe loaded module name.')
    return matches[0] if matches else None


def _executable(program):
    found = shutil.which(program)
    if not found:
        return None
    return Path(found).absolute()


def _installation_paths(path):
    """Check the PATH entry before a symlink target for adjacent support files."""
    return dict.fromkeys((path, path.resolve()))


def render_site_setup(programs_by_backend):
    """Build a user-editable setup script from visible programs and modules."""
    lines = [
        '#!/usr/bin/env bash',
        '# Generated from the environment visible during prepare.',
        '# Review this file and preflight before submitting licensed QC jobs.',
        'case "${KINBOT_BACKEND:-}" in',
    ]
    for backend, programs in sorted(programs_by_backend.items()):
        lines.append(f'  {backend})')
        module = _loaded_module(backend)
        if module:
            lines += [
                '    if ! type module >/dev/null 2>&1; then',
                '      if [ -r "${MODULESHOME:-}/init/bash" ]; then',
                '        source "$MODULESHOME/init/bash"',
                '      else',
                '        echo "Module initialization is unavailable" >&2; return 1',
                '      fi',
                '    fi',
                f'    module load {shlex.quote(module)}',
            ]
        found = [_executable(program) for program in sorted(programs)]
        dirs = []
        for path in found:
            if path and str(path.parent) not in dirs:
                dirs.append(str(path.parent))
        for directory in reversed(dirs):
            lines.append(f'    export PATH={shlex.quote(directory)}:"$PATH"')
        if backend in ('gaussian', 'gauss'):
            profiles = (candidate.parent / 'bsd' / 'g16.profile'
                        for path in found if path and path.name == 'g16'
                        for candidate in _installation_paths(path))
            profile = next((candidate for candidate in profiles
                            if candidate.is_file()), None)
            if profile:
                lines += [
                    f'    export g16root={shlex.quote(str(profile.parent.parent.parent))}',
                    '    export GAUSS_SCRDIR="${GAUSS_SCRDIR:-${SLURM_TMPDIR:-${TMPDIR:-$PWD}}}"',
                    '    mkdir -p "$GAUSS_SCRDIR"',
                    '    set +u',
                    f'    source {shlex.quote(str(profile))}',
                    '    set -u',
                ]
        if backend == 'cfour':
            candidates = (candidate.parent.parent / 'basis' / 'GENBAS'
                          for path in found if path and path.name == 'xcfour'
                          for candidate in _installation_paths(path))
            configured = os.environ.get('CFOUR_GENBAS')
            derived = next((candidate for candidate in candidates
                            if candidate.is_file()), None)
            genbas = (str(Path(configured).expanduser().absolute()) if configured
                      else str(derived) if derived else None)
            if genbas:
                lines += [
                    '    if [ -z "${CFOUR_GENBAS:-}" ]; then',
                    f'      export CFOUR_GENBAS={shlex.quote(genbas)}',
                    '    fi',
                ]
        lines.append('    ;;')
    lines += ['  *) ;;', 'esac', '']
    return '\n'.join(lines)
