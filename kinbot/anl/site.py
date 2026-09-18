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
_MOLPRO_OVERHEAD_MW = 200
_DEFAULT_MOLPRO_STACK_MW = 1024
_DEFAULT_MEMORY_PER_CORE_MB = 4096
_DEFAULT_MAX_CORES = 16
_EFFICIENT_CORE_COUNTS = (1, 2, 4, 8, 12, 16)


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


def _automatic_core_count(task, memory_mb, node_cores, limits):
    """Return a memory-safe count and the policy used to obtain it."""
    resources = task['resources']
    core_cap = min(node_cores, limits.get('max_cores_per_node', node_cores))
    performance_cap = resources.get('max_cores', _DEFAULT_MAX_CORES)
    if isinstance(performance_cap, bool) or not isinstance(performance_cap, int) \
            or performance_cap < 1:
        raise ValueError(f"{task['id']}: max_cores must be positive.")
    core_cap = min(core_cap, performance_cap)
    backend = (task.get('backend') if task.get('kind') == 'external'
               else task.get('profile', {}).get('calculator', ''))
    if isinstance(backend, str) and backend.lower() == 'molpro':
        minimum = resources.get('min_stack_mw', _DEFAULT_MOLPRO_STACK_MW)
        if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum < 32:
            raise ValueError(f"{task['id']}: min_stack_mw must be at least 32.")
        affordable = int(memory_mb * 0.85 / 8) // (_MOLPRO_OVERHEAD_MW + minimum)
        policy = ('min_stack_mw', minimum)
    else:
        minimum = resources.get('min_memory_mb_per_core',
                                _DEFAULT_MEMORY_PER_CORE_MB)
        if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum < 1:
            raise ValueError(f"{task['id']}: min_memory_mb_per_core must be positive.")
        affordable = int(memory_mb * 0.85) // minimum
        policy = ('min_memory_mb_per_core', minimum)
    affordable = min(core_cap, affordable)
    if affordable < 1:
        return None, policy
    if affordable > _EFFICIENT_CORE_COUNTS[-1]:
        return affordable, policy
    candidates = [count for count in _EFFICIENT_CORE_COUNTS if count <= affordable]
    return (max(candidates) if candidates else affordable), policy


def assign_partitions(spec):
    """Select partitions and resolve optional node-memory-based core counts.

    Explicit task partitions always win. A task may omit ``cores`` and/or
    ``memory_mb``; preparation then sizes it from the smallest node reported
    for the chosen partition, and requests all node memory when uncapped.
    Explicit resources remain useful for small smoke tests and local mocks.
    """
    tasks = spec.get('tasks', [])
    if not isinstance(tasks, list):
        raise ValueError('tasks must be a list.')
    pending = []
    for task in tasks:
        if not isinstance(task, dict) or not isinstance(task.get('resources'), dict):
            raise ValueError('Every task needs a resources object.')
        resources = task['resources']
        auto_cores = resources.get('cores') in (None, 'auto')
        auto_memory = resources.get('memory_mb') in (None, 'auto', 'node')
        if not resources.get('partition') or auto_cores or auto_memory:
            pending.append((task, auto_cores, auto_memory))
    if not pending:
        return
    entries = _partitions()
    if not entries:
        if any(auto_cores or auto_memory for _, auto_cores, auto_memory in pending):
            raise RuntimeError('Automatic cores or node memory require sinfo on '
                               'the target Slurm site during prepare.')
        if shutil.which('sinfo'):
            raise RuntimeError('sinfo returned no usable up partitions; set '
                               'resources.partition explicitly for each task.')
        return
    limits = spec.get('limits', {})
    if not isinstance(limits, dict):
        raise ValueError('limits must be an object.')
    for name in ('max_cores_per_node', 'max_memory_mb_per_node'):
        if name in limits and (isinstance(limits[name], bool)
                               or not isinstance(limits[name], int)
                               or limits[name] < 1):
            raise ValueError(f'{name} must be a positive integer.')
    for task, auto_cores, auto_memory in pending:
        resources = task['resources']
        needed = _duration_seconds(resources['walltime'])
        selected = resources.get('partition')
        fits = [entry for entry in entries
                if (not selected or entry['name'] == selected)
                and (auto_cores or entry['cores'] >= resources['cores'])
                and (auto_memory or entry['memory_mb'] >= resources['memory_mb'])
                and entry['seconds'] >= needed]
        if not fits:
            raise RuntimeError(f"{task['id']}: no available Slurm partition "
                               'fits requested cores, memory, and walltime; '
                               'adjust resources or set a partition explicitly.')
        groups = {}
        for entry in fits:
            groups.setdefault(entry['name'], []).append(entry)
        choices = []
        for group in groups.values():
            # A partition can contain several node types. The job can land
            # on any eligible node, so every such type must meet the floor.
            node_cores = min(entry['cores'] for entry in group)
            node_memory = min(entry['memory_mb'] for entry in group)
            memory = node_memory if auto_memory else resources['memory_mb']
            if auto_memory and 'max_memory_mb_per_node' in limits:
                memory = min(memory, limits['max_memory_mb_per_node'])
            cores, policy = (_automatic_core_count(task, memory, node_cores, limits)
                             if auto_cores else (resources['cores'], None))
            if cores is None:
                continue
            representative = min(group, key=lambda item: (
                item['seconds'], not item['default'], item['name']))
            choices.append((representative, memory, cores, policy))
        if not choices:
            raise RuntimeError(f"{task['id']}: node memory cannot support one "
                               'core at the configured per-core minimum.')
        chosen, memory, cores, policy = min(choices, key=lambda item: (
            item[0]['seconds'], not item[0]['default'], item[0]['name']))
        resources['partition'] = chosen['name']
        if auto_memory:
            resources['memory_mb'] = memory
            resources['use_all_node_memory'] = 'max_memory_mb_per_node' not in limits
        if auto_cores:
            resources['cores'] = cores
            resources[policy[0]] = policy[1]


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
                    '    for kinbot_gaussian_scratch in "${GAUSS_SCRDIR:-}" '
                    '"${SLURM_TMPDIR:-}" "${SCRATCH:-}" "${TMPDIR:-}" "$PWD"; do',
                    '      if [ -n "$kinbot_gaussian_scratch" ] && '
                    'mkdir -p "$kinbot_gaussian_scratch" 2>/dev/null && '
                    '[ -w "$kinbot_gaussian_scratch" ]; then',
                    '        export GAUSS_SCRDIR="$kinbot_gaussian_scratch"',
                    '        break',
                    '      fi',
                    '    done',
                    '    unset kinbot_gaussian_scratch',
                    '    if [ ! -w "${GAUSS_SCRDIR:-}" ]; then',
                    '      echo "No writable Gaussian scratch directory" >&2; return 1',
                    '    fi',
                    '    kinbot_gaussian_scratch_selected="$GAUSS_SCRDIR"',
                    '    set +u',
                    f'    if source {shlex.quote(str(profile))}; then',
                    '      :',
                    '    else',
                    '      kinbot_gaussian_profile_status=$?',
                    '      set -euo pipefail',
                    '      echo "Gaussian profile failed with exit status '
                    '$kinbot_gaussian_profile_status" >&2',
                    '      return "$kinbot_gaussian_profile_status"',
                    '    fi',
                    '    set -euo pipefail',
                    '    export GAUSS_SCRDIR="$kinbot_gaussian_scratch_selected"',
                    '    unset kinbot_gaussian_scratch_selected',
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
