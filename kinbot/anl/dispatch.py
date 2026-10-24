"""Molecule-independent, restartable QC task dispatch for exclusive Slurm jobs.

This is an execution layer. An external task is *executed* when its process and
declared artifacts pass checks; that does not validate a scientific energy.
Program-specific result parsers and MESS handoff are separate gates.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import socket
import subprocess
import sys
import time
import traceback

from ase import Atoms
from ase.io import read, write
from ase.units import Hartree
import numpy as np

from kinbot.ase_modules.calculators.factory import capabilities
from kinbot.anl.cfour import normalize_cfour_zmat
from kinbot.anl.runtime import cleanup_qc_runtime, qc_runtime_environment
from kinbot.anl.site import assign_partitions, render_site_setup
from kinbot.anl.tasks import mrcc_task
from kinbot.execution import (atomic_json as _atomic_json,
                              file_sha256 as _file_hash,
                              scheduler_unavailable as _scheduler_unavailable)
from kinbot.theory import TheoryProfile
from kinbot.gaussian_input import has_anharmonic_frequency
from kinbot.anl.vpt2_recovery import (
    VPT2RecoveryUnavailable,
    framework_history as _vpt2_framework_history,
    plan as _vpt2_framework_recovery,
)


_SAFE_NAME = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*\Z')
_TOKEN = re.compile(r'\{\{([A-Z_]+)\}\}')
_TIME = re.compile(r'(?:\d+-)?\d{1,3}:[0-5]\d:[0-5]\d\Z')
_INPUT_TOKENS = {'CARTESIAN', 'XYZ', 'MRCC_XYZ', 'NATOMS', 'CHARGE',
                 'MULT', 'NELECTRONS', 'SPIN', 'CORES', 'MEMORY_MB',
                 'WORK_MEMORY_MB'}
_RESERVED_FILES = {'geometry.xyz', 'task.json', 'job.slurm', 'execution.json',
                   'slurm.stdout', 'slurm.stderr', 'drive.lock'}
_SCHEDULER_TIMEOUT = 30
_SCHEDULER_QUERY_ATTEMPTS = 3
_RESULT_PUBLICATION_GRACE = 60.


class CalculationUnavailableError(RuntimeError):
    """A native calculation could not be completed with its resources.

    This is deliberately narrower than a result-parser or provenance error.
    Optional ANL component work may fall back after this exception, whereas
    malformed or scientifically inconsistent output must remain a hard
    failure that is fixed rather than silently omitted.
    """


class SchedulerUnavailableError(RuntimeError):
    """A scheduler status query failed without proving a job has ended.

    Automatic drivers must retain submitted state and retry this condition.
    Treating a Slurm controller outage as an absent job can otherwise turn a
    healthy calculation into a false scientific failure.
    """


def preflight_until_scheduler_ready(run_dir, *, skip_optional=False,
                                    interval=20., progress=None):
    """Retry only scheduler availability failures during preflight.

    Executable, input, resource, and scientific validation failures remain
    immediate hard errors. This wrapper is for long-lived automatic drivers
    whose prepared work must survive a temporary slurmctld outage.
    """
    previous = None
    while True:
        try:
            return preflight(run_dir, skip_optional=skip_optional,
                             progress=progress)
        except SchedulerUnavailableError as error:
            message = str(error)
            if progress is not None:
                progress({'phase': 'waiting_for_scheduler', 'detail': message})
            if message != previous:
                print(
                    'Scheduler preflight is unavailable; retaining staged '
                    f'tasks and retrying: {message}',
                    file=sys.stderr, flush=True)
                previous = message
            time.sleep(float(interval))


def _basename(value, label):
    if not isinstance(value, str) or not _SAFE_NAME.fullmatch(value):
        raise ValueError(f'{label} must be a safe basename.')
    return value


def _positive_integer(value, label):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f'{label} must be a positive integer.')
    return value


def _molpro_total_mw(resources):
    """Budget aggregate Molpro stack after documented per-process overhead."""
    total_mw = int(resources['memory_mb'] * 0.85 / 8) - 200 * resources['cores']
    if total_mw // resources['cores'] < 32:
        raise ValueError('Molpro memory allocation leaves no usable stack '
                         'after per-process overhead; request more node memory '
                         'or fewer cores.')
    return total_mw


def _molpro_stack_mw(resources):
    """Per-process -m allocation for Molpro's single-node disk mode."""
    return _molpro_total_mw(resources) // resources['cores']


def _atoms(molecule):
    if not isinstance(molecule, dict):
        raise ValueError('molecule must be an object.')
    symbols = molecule.get('symbols')
    positions = molecule.get('positions')
    if not isinstance(symbols, list) or not symbols or not isinstance(positions, list):
        raise ValueError('molecule needs symbols and positions lists.')
    try:
        atoms = Atoms(symbols=symbols, positions=positions)
    except (TypeError, ValueError, KeyError) as exc:
        raise ValueError('Invalid molecule symbols or positions.') from exc
    if any(number < 1 for number in atoms.numbers) or not np.isfinite(atoms.positions).all():
        raise ValueError('Molecule has an invalid element or coordinate.')
    charge = molecule.get('charge', 0)
    mult = molecule.get('multiplicity', 1)
    if isinstance(charge, bool) or not isinstance(charge, int):
        raise ValueError('Molecular charge must be an integer.')
    _positive_integer(mult, 'Molecular multiplicity')
    electrons = int(sum(atoms.numbers)) - charge
    if electrons < mult - 1 or (electrons - mult + 1) % 2:
        raise ValueError('Charge and multiplicity conflict with electron count.')
    return atoms


def _geometry_hash(atoms):
    payload = json.dumps({'symbols': atoms.get_chemical_symbols(),
                          'positions': atoms.positions.tolist()},
                         sort_keys=True, separators=(',', ':'))
    return hashlib.sha256(payload.encode()).hexdigest()


def _validate_task(task, ids, limits, *, allow_obsolete=False):
    if not isinstance(task, dict):
        raise ValueError('Every task must be an object.')
    ident = _basename(task.get('id'), 'Task id')
    kind = task.get('kind')
    if kind not in ('external', 'ase_optimize'):
        raise ValueError(f'{ident}: unsupported task kind {kind!r}.')
    source = task.get('geometry_from', 'initial')
    if source != 'initial' and source not in ids:
        raise ValueError(f'{ident}: unknown geometry source {source!r}.')
    fallback = task.get('geometry_fallback')
    if fallback is not None:
        if (source == 'initial' or fallback == source
                or (fallback != 'initial' and fallback not in ids)):
            raise ValueError(f'{ident}: invalid geometry_fallback {fallback!r}.')
    dependencies = task.get('depends_on', [])
    if not isinstance(dependencies, list) or any(dep not in ids for dep in dependencies):
        raise ValueError(f'{ident}: depends_on must contain task ids.')
    fallback_on = task.get('fallback_on', [])
    if ('fallback_on' in task
            and (not isinstance(fallback_on, list) or not fallback_on)):
        raise ValueError(f'{ident}: fallback_on must contain task ids.')
    if (any(dep not in ids or dep == ident for dep in fallback_on)
            or len(fallback_on) != len(set(fallback_on))
            or set(fallback_on) & set(dependencies)):
        raise ValueError(
            f'{ident}: fallback_on must contain distinct alternative task ids.')
    failure_policy = task.get('failure_policy')
    if failure_policy not in (None, 'skip_component'):
        raise ValueError(f'{ident}: unsupported failure_policy.')
    component_families = task.get('component_families', [])
    if (not isinstance(component_families, list)
            or any(not isinstance(value, str)
                   or not _SAFE_NAME.fullmatch(value)
                   for value in component_families)
            or len(component_families) != len(set(component_families))):
        raise ValueError(
            f'{ident}: component_families must contain unique safe names.')
    failure_skips = task.get('failure_skips', [])
    if (not isinstance(failure_skips, list)
            or any(not isinstance(value, str)
                   or not _SAFE_NAME.fullmatch(value) or value == ident
                   for value in failure_skips)
            or len(failure_skips) != len(set(failure_skips))):
        raise ValueError(f'{ident}: failure_skips must contain other task ids.')
    if failure_policy == 'skip_component' and not component_families:
        raise ValueError(
            f'{ident}: skip_component requires a component family.')
    if failure_policy is None and (component_families or failure_skips):
        raise ValueError(
            f'{ident}: component fallback metadata requires failure_policy.')
    resources = task.get('resources')
    if not isinstance(resources, dict):
        raise ValueError(f'{ident}: resources must be an object.')
    cores = _positive_integer(resources.get('cores'), f'{ident} cores')
    memory = _positive_integer(resources.get('memory_mb'), f'{ident} memory_mb')
    if (cores > limits.get('max_cores_per_node', cores)
            or memory > limits.get('max_memory_mb_per_node', memory)):
        raise ValueError(f'{ident}: task exceeds per-node resource limits.')
    if 'use_all_node_memory' in resources and not isinstance(resources['use_all_node_memory'], bool):
        raise ValueError(f'{ident}: use_all_node_memory must be boolean.')
    for key in ('max_cores', 'min_memory_mb_per_core', 'min_scratch_mb'):
        if key in resources:
            _positive_integer(resources[key], f'{ident} {key}')
    if 'min_stack_mw' in resources:
        if (isinstance(resources['min_stack_mw'], bool)
                or not isinstance(resources['min_stack_mw'], int)
                or resources['min_stack_mw'] < 32):
            raise ValueError(f'{ident}: min_stack_mw must be at least 32.')
    if not isinstance(resources.get('walltime'), str) or not _TIME.fullmatch(resources['walltime']):
        raise ValueError(f'{ident}: walltime must be HH:MM:SS or D-HH:MM:SS.')
    if 'partition' in resources:
        _basename(resources['partition'], f'{ident} partition')
    setup = task.get('setup', [])
    if (not isinstance(setup, list) or
            any(not isinstance(line, str) or '\n' in line or '\r' in line
                for line in setup)):
        raise ValueError(f'{ident}: setup must be a list of single shell lines.')
    runtime_setup = task.get('runtime_setup', [])
    if (not isinstance(runtime_setup, list) or any(not isinstance(line,str)
            or any(c in line for c in '\n\r\x00') for line in runtime_setup)):
        raise ValueError(f'{ident}: runtime_setup must contain single shell lines.')
    geometry_output = task.get('geometry_output')
    if kind == 'ase_optimize':
        if 'result_parser' in task:
            raise ValueError(f'{ident}: result_parser requires an external task.')
        if geometry_output != 'final.xyz':
            raise ValueError(f'{ident}: ase_optimize geometry_output must be final.xyz.')
        profile = _runtime_profile(task)
        if not capabilities(profile.calculator).forces:
            raise ValueError(f'{ident}: selected ASE calculator has no force capability.')
        options = task.get('optimizer', {})
        if (not isinstance(options, dict) or
                set(options) - {'fmax', 'steps', 'sella_kwargs', 'order',
                                'constraints', 'minimum_recovery', 'step_retries'}):
            raise ValueError(f'{ident}: invalid optimizer options.')
        if ('step_retries' in options and (type(options['step_retries']) is not int
                or not 0 <= options['step_retries'] <= 10)):
            raise ValueError(f'{ident}: step_retries must be an integer from 0 to 10.')
        if ('sella_kwargs' in options and not isinstance(options['sella_kwargs'], dict)):
            raise ValueError(f'{ident}: sella_kwargs must be an object.')
        if ('fmax' in options and
                (isinstance(options['fmax'], bool) or
                 not isinstance(options['fmax'], (int, float)) or
                 not np.isfinite(options['fmax']) or options['fmax'] <= 0)):
            raise ValueError(f'{ident}: fmax must be a positive finite number.')
        if 'steps' in options:
            _positive_integer(options['steps'], f'{ident} steps')
        if (isinstance(options.get('order', 0), bool)
                or options.get('order', 0) not in (0, 1)):
            raise ValueError(f'{ident}: optimizer order must be 0 or 1.')
        if options.get('order', 0) and profile.optimizer != 'sella':
            raise ValueError(f'{ident}: a saddle optimization requires Sella.')
        if 'minimum_recovery' in options:
            from kinbot.minimum_recovery import validate_anl_recovery
            validate_anl_recovery(options['minimum_recovery'])
            if options.get('order', 0) or options.get('constraints'):
                raise ValueError(f'{ident}: minimum recovery requires an unconstrained minimum.')
        constraints = options.get('constraints', {})
        if not isinstance(constraints, dict) or set(constraints) - {
                'fixed_bonds', 'minimum_bonds'}:
            raise ValueError(f'{ident}: invalid optimizer constraints.')
        fixed = constraints.get('fixed_bonds', [])
        minimum = constraints.get('minimum_bonds', [])
        if (not isinstance(fixed, list)
                or any(not isinstance(pair, list) or len(pair) != 2
                       or any(isinstance(index, bool) or not isinstance(index, int)
                              or index < 0 for index in pair)
                       or pair[0] == pair[1] for pair in fixed)):
            raise ValueError(f'{ident}: fixed_bonds must contain atom-index pairs.')
        if (not isinstance(minimum, list)
                or any(not isinstance(item, list) or len(item) != 3
                       or any(isinstance(index, bool) or not isinstance(index, int)
                              or index < 0 for index in item[:2])
                       or item[0] == item[1]
                       or isinstance(item[2], bool)
                       or not isinstance(item[2], (int, float))
                       or not np.isfinite(item[2]) or item[2] <= 0
                       for item in minimum)):
            raise ValueError(
                f'{ident}: minimum_bonds must contain [atom, atom, angstrom].')
    else:
        _basename(task.get('backend'), f'{ident} backend')
        _basename(task.get('input_name'), f'{ident} input_name')
        if not isinstance(task.get('input_template'), str) or not task['input_template'].strip():
            raise ValueError(f'{ident}: input_template is required.')
        unknown_tokens = set(_TOKEN.findall(task['input_template'])) - _INPUT_TOKENS
        if (unknown_tokens or
                _TOKEN.sub('', task['input_template']).find('{{') >= 0 or
                _TOKEN.sub('', task['input_template']).find('}}') >= 0):
            raise ValueError(f'{ident}: unknown input tokens {sorted(unknown_tokens)}.')
        command = task.get('command')
        if not isinstance(command, list) or not command or any(
                not isinstance(arg, str) or not arg for arg in command):
            raise ValueError(f'{ident}: command must be an argument list.')
        command_tokens = {token for arg in command
                          for token in re.findall(r'\{([^{}]+)\}', arg)}
        if command_tokens - {'cores', 'input', 'molpro_stack_mw'}:
            raise ValueError(f'{ident}: unsupported command placeholders.')
        required_executables = task.get('required_executables', [])
        if (not isinstance(required_executables, list)
                or any(not isinstance(program, str)
                       or not _SAFE_NAME.fullmatch(program)
                       for program in required_executables)
                or len(required_executables) != len(set(required_executables))):
            raise ValueError(f'{ident}: required_executables must contain '
                             'unique executable basenames.')
        if task.get('stdin'):
            if task['stdin'] != task['input_name']:
                raise ValueError(f'{ident}: stdin must equal input_name.')
        for key in ('stdout', 'stderr'):
            _basename(task.get(key, f'{key}.txt'), f'{ident} {key}')
        outputs = task.get('required_outputs')
        if not isinstance(outputs, list) or not outputs or any(
                not isinstance(item, str) or not _SAFE_NAME.fullmatch(item)
                for item in outputs):
            raise ValueError(f'{ident}: required_outputs must be basenames.')
        aliases = task.get('output_aliases', {})
        if not isinstance(aliases, dict):
            raise ValueError(f'{ident}: output_aliases must be an object.')
        for target, candidates in aliases.items():
            if (target not in outputs or not isinstance(candidates, list) or not candidates
                    or any(not isinstance(c, str) for c in candidates)
                    or len(set(candidates)) != len(candidates)):
                raise ValueError(f'{ident}: output aliases require one declared target and unique candidates.')
            for candidate in candidates:
                _basename(candidate, 'output alias')
                if candidate in {*outputs, *_RESERVED_FILES, task['input_name'],
                                  task.get('stdout', 'stdout.txt'), task.get('stderr', 'stderr.txt')}:
                    raise ValueError(f'{ident}: output alias collides with a task file.')
        stdout_name = task.get('stdout', 'stdout.txt')
        stderr_name = task.get('stderr', 'stderr.txt')
        if ({task['input_name'], stdout_name, stderr_name, *outputs} & _RESERVED_FILES):
            raise ValueError(f'{ident}: task filenames use reserved dispatch files.')
        if (stdout_name == stderr_name or task['input_name'] in (stdout_name, stderr_name)
                or task['input_name'] in outputs
                or (task['backend'].lower() == 'molpro' and stdout_name in outputs)):
            raise ValueError(f'{ident}: input, output, stdout, and stderr filenames must not collide.')
        backend = task['backend'].lower()
        if backend == 'molpro':
            _molpro_total_mw(resources)
            if (_molpro_stack_mw(resources)
                    < resources.get('min_stack_mw', 32)):
                raise ValueError(f'{ident}: Molpro stack per rank is below min_stack_mw.')
            native_output = Path(task['input_name']).with_suffix('.out').name
            if (not task['input_name'].endswith('.inp')
                    or native_output not in outputs or stdout_name == native_output):
                raise ValueError(f'{ident}: Molpro needs .inp, its native .out, '
                                 'and separate launcher stdout.')
        if backend == 'cfour' and task['input_name'] != 'ZMAT':
            raise ValueError(f'{ident}: CFOUR input must be ZMAT.')
        if backend == 'mrcc' and task['input_name'] != 'MINP':
            raise ValueError(f'{ident}: direct MRCC input must be MINP.')
        if backend == 'gaussian' and (
                not task['input_name'].endswith('.com')
                or task.get('stdin') != task['input_name']):
            raise ValueError(f'{ident}: Gaussian needs a .com file on stdin.')
        export = task.get('gaussian_orbital_export')
        if export is not None:
            if (backend != 'gaussian' or not isinstance(export, dict)
                    or set(export) != {'checkpoint', 'formatted', 'cube', 'log', 'orbital'}
                    or export.get('orbital') not in ('MO=HOMO', 'AMO=HOMO')):
                raise ValueError(f'{ident}: invalid Gaussian orbital export.')
            names = [export[key] for key in ('checkpoint', 'formatted', 'cube', 'log')]
            if len(set(names)) != 4 or any(name not in outputs for name in names):
                raise ValueError(f'{ident}: orbital export needs distinct declared output files.')
            if not {'formchk', 'cubegen'} <= set(required_executables):
                raise ValueError(f'{ident}: orbital export needs formchk/cubegen preflight.')
        staged_files = task.get('files_from_env', {})
        if (not isinstance(staged_files, dict) or any(
                not isinstance(env, str) or not re.fullmatch(r'[A-Z_][A-Z_0-9]*', env)
                for env in staged_files.values())):
            raise ValueError(f'{ident}: files_from_env must map filenames to env vars.')
        for filename in staged_files:
            _basename(filename, f'{ident} staged filename')
            if filename in (task['input_name'], stdout_name, stderr_name) \
                    or filename in _RESERVED_FILES or filename in outputs:
                raise ValueError(f'{ident}: staged filename collides with task files.')
        if geometry_output:
            _basename(geometry_output, f'{ident} geometry_output')
            if geometry_output not in outputs:
                raise ValueError(f'{ident}: geometry_output must be required.')
            if not task.get('success_marker'):
                raise ValueError(f'{ident}: geometry needs a success_marker.')
        if task.get('success_marker'):
            marker = task['success_marker']
            if (not isinstance(marker, dict) or set(marker) != {'file', 'contains'}
                    or not isinstance(marker['contains'], str)
                    or not marker['contains']):
                raise ValueError(f'{ident}: invalid success_marker.')
            _basename(marker['file'], f'{ident} success_marker file')
            if marker['file'] not in set(outputs) | {stdout_name, stderr_name}:
                raise ValueError(f'{ident}: success_marker file must be a required output or stream.')
        result_parser = task.get('result_parser')
        if result_parser is not None:
            from kinbot.anl.results import validate_result_parser
            try:
                validate_result_parser(result_parser, backend=backend,
                                       template=task['input_template'],
                                       outputs=outputs,
                                       allow_obsolete=allow_obsolete)
            except ValueError as exc:
                raise ValueError(f'{ident}: {exc}') from exc
        failure_markers = task.get('failure_markers', [])
        if (not isinstance(failure_markers, list) or any(
                not isinstance(marker, dict) or set(marker) != {'file', 'contains'}
                or marker.get('file') not in set(outputs) | {stdout_name, stderr_name}
                or not isinstance(marker.get('contains'), str)
                or not marker['contains'] for marker in failure_markers)):
            raise ValueError(f'{ident}: failure_markers must refer to output files.')


def validate_spec(spec, *, allow_obsolete=False):
    """Validate a declarative graph before creating or submitting any job."""
    if not isinstance(spec, dict) or spec.get('schema') != 1:
        raise ValueError('Workflow must be a schema 1 object.')
    _basename(spec.get('name'), 'Workflow name')
    atoms = _atoms(spec.get('molecule'))
    limits = spec.get('limits')
    if not isinstance(limits, dict):
        raise ValueError('limits must be an object.')
    _positive_integer(limits.get('max_nodes'), 'max_nodes')
    for key in ('max_cores_per_node', 'max_memory_mb_per_node'):
        if key in limits:
            _positive_integer(limits[key], key)
    tasks = spec.get('tasks')
    if not isinstance(tasks, list):
        raise ValueError('tasks must be a list.')
    if not tasks:
        rank_exact = spec.get('intent', {}).get('rank_exact_higher_order')
        electrons = int(sum(atoms.numbers)) - spec['molecule'].get('charge', 0)
        if (not isinstance(rank_exact, dict)
                or rank_exact.get('electron_count') not in (1, 2)
                or rank_exact.get('electron_count') != electrons
                or rank_exact.get('correction_hartree') != 0.):
            raise ValueError('tasks may be empty only for a rank-exact '
                             'higher-order workflow.')
    ids = [_basename(task.get('id') if isinstance(task, dict) else None, 'Task id')
           for task in tasks]
    if len(ids) != len(set(ids)):
        raise ValueError('Task ids must be unique.')
    by_id = dict(zip(ids, tasks))
    for task in tasks:
        _validate_task(task, by_id, limits, allow_obsolete=allow_obsolete)
        source = task.get('geometry_from', 'initial')
        if source != 'initial' and not by_id[source].get('geometry_output'):
            raise ValueError(f"{task['id']}: geometry source has no geometry output.")
        parser = task.get('result_parser', {})
        if 'rank_exact_electrons' in parser:
            atoms = _atoms(spec['molecule'])
            electrons = int(sum(atoms.numbers)) - spec['molecule'].get(
                'charge', 0)
            if parser['rank_exact_electrons'] != electrons:
                raise ValueError(
                    f"{task['id']}: rank-exact electron count disagrees "
                    'with the molecular state.')
        if parser.get('kind') in ('mrcc_energy',
                                  'molpro_energy') and 'reference' in parser:
            multiplicity = spec['molecule'].get('multiplicity', 1)
            reference = parser['reference']
            # Explicit historical-reference reproductions are diagnostics.
            # Validate their declared provenance, never a particular molecule
            # or benchmark name. Automatic ANL still emits RHF/ROHF throughout.
            intent = spec.get('intent', {})
            benchmark = intent.get('literature_benchmark', {})
            source_validation = task.get('source_reference_validation')
            declared_source_uhf = (
                parser.get('kind') == 'mrcc_energy'
                and reference == 'UHF'
                and isinstance(source_validation, str) and bool(source_validation.strip())
                and intent.get('source_reference_validation') == source_validation
                and intent.get('claim') == 'pinned-literature-source-reference-validation-only'
                and isinstance(benchmark, dict) and bool(benchmark.get('name'))
                and bool(benchmark.get('source'))
                and multiplicity > 1)
            if (not declared_source_uhf
                    and (reference not in ('RHF', 'ROHF')
                    or (reference == 'RHF' and multiplicity != 1)
                    or (reference == 'ROHF' and multiplicity == 1))):
                raise ValueError(f"{task['id']}: declared reference conflicts "
                                 'with molecular multiplicity.')
        if (parser.get('kind') == 'gaussian_vpt2'
                and re.search(r'\bOpt\s*(?:=|\()', task['input_template'],
                              re.IGNORECASE) is None):
            source_task = by_id.get(source)
            profile = (source_task or {}).get('profile', {})
            keywords = profile.get('calculator_kwargs', {}) if isinstance(profile, dict) else {}
            task_source_matches = (
                source_task is not None
                and source_task.get('kind') == 'ase_optimize'
                and isinstance(profile, dict)
                and isinstance(keywords, dict)
                and profile.get('calculator', '').lower()
                in ('gaussian', 'gauss')
                and profile.get('method', '').casefold()
                == parser['method'].casefold()
                and profile.get('basis', '').casefold()
                == parser['basis'].casefold()
                and keywords.get('EmpiricalDispersion', '').casefold()
                == parser.get('dispersion', '').casefold())
            imported_source = spec.get('molecule', {}).get('source', {})
            imported_profile = (imported_source.get('l2_theory', {})
                                if isinstance(imported_source, dict) else {})
            imported_source_matches = (
                source == 'initial'
                and isinstance(imported_profile, dict)
                and re.fullmatch(r'[0-9a-f]{64}', str(
                    imported_source.get('geometry_sha256', '')))
                and re.fullmatch(r'[0-9a-f]{64}', str(
                    imported_source.get('record_sha256', '')))
                and imported_source.get('geometry_sha256')
                == _geometry_hash(atoms)
                and imported_profile.get('calculator', '').lower()
                in ('gaussian', 'gauss')
                and imported_profile.get('method', '').casefold()
                == parser['method'].casefold()
                and imported_profile.get('basis', '').casefold()
                == parser['basis'].casefold()
                and imported_profile.get('dispersion', '').casefold()
                == parser.get('dispersion', '').casefold())
            if not (task_source_matches or imported_source_matches):
                raise ValueError(f"{task['id']}: frequency-only VPT2 must use "
                                 'a matching Gaussian geometry source.')
    visiting, visited = set(), set()

    def visit(ident):
        if ident in visiting:
            raise ValueError('Task dependency cycle detected.')
        if ident in visited:
            return
        visiting.add(ident)
        task = by_id[ident]
        predecessors = set(task.get('depends_on', []))
        if task.get('geometry_from', 'initial') != 'initial':
            predecessors.add(task['geometry_from'])
        for predecessor in predecessors:
            visit(predecessor)
        visiting.remove(ident)
        visited.add(ident)

    for ident in ids:
        visit(ident)
    return spec


def _load(run_dir):
    run_dir = Path(run_dir).resolve()
    marker = run_dir/'.manifest.lock'
    if marker.exists():
        with marker.open('r') as handle:
            fcntl.flock(handle,fcntl.LOCK_SH)
            return _load_unlocked(run_dir)
    return _load_unlocked(run_dir)


def _load_unlocked(run_dir):
    spec = json.loads((run_dir / 'workflow.json').read_text())
    validate_spec(spec, allow_obsolete=True)
    state = json.loads((run_dir / 'state.json').read_text())
    checksum = hashlib.sha256((run_dir / 'workflow.json').read_bytes()).hexdigest()
    if state['spec_sha256'] != checksum:
        raise RuntimeError('Workflow changed after preparation; start a new run.')
    return run_dir, spec, state


def _geometry_source(task, state=None):
    """Return the primary or declared fallback geometry source.

    The choice becomes deterministic only after the primary geometry task is
    terminal. It is recorded in the staged state, so retries and provenance
    checks never infer a different geometry later.
    """
    source = task.get('geometry_from', 'initial')
    fallback = task.get('geometry_fallback')
    if (state is not None and fallback is not None and source != 'initial'
            and state['tasks'].get(source, {}).get('status') in
            ('unavailable', 'skipped')):
        return fallback
    return source


def _task_dependencies(task, state=None):
    dependencies = set(task.get('depends_on', []))
    source = _geometry_source(task, state)
    if source != 'initial':
        dependencies.add(source)
    return dependencies


def _render_input(template, atoms, molecule, resources, backend=None):
    symbols = atoms.get_chemical_symbols()
    coordinates = '\n'.join(
        f'{symbol:<2} {x: .12f} {y: .12f} {z: .12f}'
        for symbol, (x, y, z) in zip(symbols, atoms.positions))
    replacements = {
        'CARTESIAN': coordinates,
        'XYZ': f'{len(atoms)}\nKinBot accepted geometry\n{coordinates}',
        'MRCC_XYZ': f'{len(atoms)}\n\n{coordinates}',
        'NATOMS': str(len(atoms)),
        'CHARGE': str(molecule.get('charge', 0)),
        'MULT': str(molecule.get('multiplicity', 1)),
        'NELECTRONS': str(int(sum(atoms.numbers)) - molecule.get('charge', 0)),
        'SPIN': str(molecule.get('multiplicity', 1) - 1),
        'CORES': str(resources['cores']),
        'MEMORY_MB': str(resources['memory_mb']),
        'WORK_MEMORY_MB': str(max(1, int(resources['memory_mb'] * 0.7))),
    }

    def replace(match):
        try:
            return replacements[match.group(1)]
        except KeyError as exc:
            raise ValueError(f'Unknown input token {match.group(0)}.') from exc

    rendered = _TOKEN.sub(replace, template)
    if (backend and backend.lower() in ('gaussian', 'gauss')
            and has_anharmonic_frequency(rendered)
            and len(atoms) > 1):
        # Symmetry=None fixes nonlinear framework-group mismatches, but also
        # suppresses the linear classification needed by Gaussian Link 716.
        # Keep the established nonlinear workaround; only linear geometries
        # retain point-group detection through NoSymm.
        singular = np.linalg.svd(atoms.positions - atoms.positions.mean(axis=0),
                                 compute_uv=False)
        if singular[1] < 1.e-7:
            rendered = re.sub(r'(?i)\bSymm(?:etry)?\s*=\s*(?:\(\s*None\s*\)|None\b)',
                              'NoSymm', rendered)
    if backend and backend.lower() == 'molpro':
        from kinbot.molpro_diagnostics import ensure_native_printing
        rendered = ensure_native_printing(rendered)
    return rendered


def _runtime_profile(task):
    """Bind calculator process/memory settings to this task's resources."""
    values = deepcopy(task.get('profile'))
    if not isinstance(values, dict):
        raise ValueError(f"{task['id']}: profile must be an object.")
    backend = values.get('calculator', '').lower()
    if backend == 'molpro':
        resources = task['resources']
        kwargs = values.setdefault('calculator_kwargs', {})
        if not isinstance(kwargs, dict):
            raise ValueError(f"{task['id']}: calculator_kwargs must be an object.")
        cores = resources['cores']
        stack_mw = _molpro_stack_mw(resources)
        if stack_mw < resources.get('min_stack_mw', 32):
            raise ValueError(f"{task['id']}: Molpro stack per rank is below min_stack_mw.")
        if resources.get('live_override'):
            kwargs.update(nproc=cores, stack_mw=stack_mw)
        if 'nproc' in kwargs and kwargs['nproc'] != cores:
            raise ValueError(f"{task['id']}: Molpro nproc must match Slurm cores.")
        if 'stack_mw' in kwargs and kwargs['stack_mw'] != stack_mw:
            raise ValueError(f"{task['id']}: Molpro stack_mw must match the node budget.")
        kwargs.update(nproc=cores, stack_mw=stack_mw)
        kwargs.setdefault(
            'scratch_min_mb',
            resources.get('min_scratch_mb', max(4096, resources['memory_mb'] // 2)))
    if backend in ('gaussian', 'gauss'):
        kwargs = values.setdefault('calculator_kwargs', {})
        if not isinstance(kwargs, dict):
            raise ValueError(f"{task['id']}: calculator_kwargs must be an object.")
        cores = task['resources']['cores']
        if task['resources'].get('live_override'):
            kwargs.update(nprocshared=cores,
                          mem=f"{max(1, int(task['resources']['memory_mb'] * 0.7))}MB")
        if 'nprocshared' in kwargs and kwargs['nprocshared'] != cores:
            raise ValueError(f"{task['id']}: Gaussian nprocshared must match Slurm cores.")
        kwargs.setdefault('nprocshared', cores)
        kwargs.setdefault('mem',
                          f"{max(1, int(task['resources']['memory_mb'] * 0.7))}MB")
        match = re.fullmatch(r'(\d+)(MB|GB)', str(kwargs['mem']), re.I)
        if not match:
            raise ValueError(f"{task['id']}: Gaussian mem must use MB or GB.")
        memory_mb = int(match.group(1)) * (1024 if match.group(2).upper() == 'GB' else 1)
        if memory_mb > task['resources']['memory_mb'] * 0.8:
            raise ValueError(f"{task['id']}: Gaussian mem leaves too little node headroom.")
    return TheoryProfile.from_dict(task['id'], values)


def _backend(task):
    return (task['backend'] if task['kind'] == 'external'
            else _runtime_profile(task).calculator).lower()


def _omp_threads(task):
    # Molpro -n starts one MPI process per requested core. Giving every MPI
    # rank all cores again through OpenMP oversubscribes the exclusive node.
    return 1 if _backend(task) == 'molpro' else task['resources']['cores']


def _slurm_script(task, directory, python):
    from kinbot.qc import _slurm_exclude_directive

    resources = task['resources']
    mpi_ranks = resources['cores'] if _backend(task) == 'molpro' else 1
    cpus_per_rank = 1 if _backend(task) == 'molpro' else resources['cores']
    lines = [
        '#!/usr/bin/env bash',
        f"#SBATCH --job-name=kb-{task['id']}",
        '#SBATCH --nodes=1',
        f'#SBATCH --ntasks={mpi_ranks}',
        f'#SBATCH --cpus-per-task={cpus_per_rank}',
        ("#SBATCH --mem=0" if resources.get('use_all_node_memory')
         else f"#SBATCH --mem={resources['memory_mb']}M"),
        f"#SBATCH --time={resources['walltime']}",
        '#SBATCH --exclusive',
        '#SBATCH --output=slurm.stdout',
        '#SBATCH --error=slurm.stderr',
    ]
    if resources.get('partition'):
        lines.append(f"#SBATCH --partition={resources['partition']}")
    exclusion = _slurm_exclude_directive()
    if exclusion:
        lines.append(exclusion)
    lines += ['', 'set -euo pipefail', f'cd {shlex.quote(str(directory))}',
              f"export OMP_NUM_THREADS={_omp_threads(task)}",
              f"export KINBOT_BACKEND={shlex.quote(_backend(task))}",
              'source ../../site_setup.sh']
    lines += task.get('setup', [])
    lines.append(f'export OMP_NUM_THREADS={_omp_threads(task)}')
    if _backend(task) == 'molpro':
        lines.append('export MKL_NUM_THREADS=1')
    elif _backend(task) == 'mrcc':
        lines.append(f'export MKL_NUM_THREADS={_omp_threads(task)}')
    lines += task.get('runtime_setup', [])
    lines.append(f'{shlex.quote(python)} -m kinbot.anl.dispatch run-task task.json')
    from kinbot.slurm_retry import wrap_script
    return wrap_script('\n'.join(lines) + '\n', identity=task['id'],
                       directories=[str(directory)], python=python)


def _stage_task(run_dir, spec, state, task):
    (run_dir / '.manifest.lock').touch(exist_ok=True)
    ident = task['id']
    directory = run_dir / 'tasks' / ident
    if directory.exists():
        raise RuntimeError(f'{ident}: existing task directory needs manual review.')
    source = _geometry_source(task, state)
    if source == 'initial':
        atoms = _atoms(spec['molecule'])
    else:
        source_entry = state['tasks'][source]
        source_file = run_dir / 'tasks' / source / next(
            item['geometry_output'] for item in spec['tasks'] if item['id'] == source)
        atoms = read(source_file)
        if (_geometry_hash(atoms) != source_entry.get('final_geometry_sha256')
                or source_entry['status'] != 'complete'):
            raise RuntimeError(f'{ident}: accepted source geometry changed.')
    initial = _atoms(spec['molecule'])
    if atoms.get_chemical_symbols() != initial.get_chemical_symbols():
        raise RuntimeError(f'{ident}: geometry changed atom identity or order.')
    if not np.isfinite(atoms.positions).all():
        raise RuntimeError(f'{ident}: geometry contains nonfinite coordinates.')
    directory.mkdir(parents=True)
    write(directory / 'geometry.xyz', atoms)
    atoms = read(directory / 'geometry.xyz')
    if task['kind'] == 'external':
        rendered = _render_input(task['input_template'], atoms,
                                 spec['molecule'], task['resources'], _backend(task))
        if _backend(task) == 'cfour':
            rendered = normalize_cfour_zmat(rendered)
        (directory / task['input_name']).write_text(rendered)
    from kinbot.anl.optimization_recovery import stage_seed
    continuation = stage_seed(run_dir, directory, task, state, atoms)
    _atomic_json(directory / 'task.json', {
        'schema': 1, 'live_control_schema': 1, 'task': task,
        **({'optimizer_continuation': continuation} if continuation else {}),
        'molecule': {'symbols': initial.get_chemical_symbols(),
                     'charge': spec['molecule'].get('charge', 0),
                     'multiplicity': spec['molecule'].get('multiplicity', 1)},
        'geometry_sha256': _geometry_hash(read(directory / 'geometry.xyz')),
    })
    (directory / 'job.slurm').write_text(
        _slurm_script(task, directory, state['python']))
    from kinbot.anl.scratch_recovery import apply_exclusions
    script = directory / 'job.slurm'
    script.write_text(apply_exclusions(script.read_text(), state))
    state['tasks'][ident] = {
        'status': 'staged', 'geometry_from': source,
        'geometry_sha256': _geometry_hash(atoms),
        'task_sha256': _file_hash(directory / 'task.json'),
        'job_sha256': _file_hash(directory / 'job.slurm'),
    }
    if task['kind'] == 'external':
        state['tasks'][ident]['input_sha256'] = _file_hash(
            directory / task['input_name'])
    # Explicit outage recovery may defer a retry until its L3 source exists.
    # Keep the attempt/archive lineage when that dependency later completes.
    pending = state.get('pending_retries', {})
    if ident in pending:
        lineage = pending.pop(ident)
        state['tasks'][ident].update({key: lineage[key] for key in
            ('attempt', 'previous_attempt', 'outage_recovery_reason') if key in lineage})


def _programs_by_backend(spec):
    """Collect every executable that a prepared workflow must expose."""
    programs_by_backend = {}
    for task in spec['tasks']:
        backend = _backend(task)
        if task['kind'] == 'external':
            program = task['command'][0]
        else:
            parts = shlex.split(_runtime_profile(task).command or
                                ('molpro' if backend == 'molpro' else 'g16'))
            if not parts:
                raise ValueError(f"{task['id']}: calculator command is empty.")
            program = parts[0]
        programs_by_backend.setdefault(backend, set()).add(program)
        programs_by_backend[backend].update(task.get('required_executables', []))
    return programs_by_backend


def prepare(spec_path, run_dir):
    if isinstance(spec_path, dict):
        spec = deepcopy(spec_path)
    else:
        spec = json.loads(Path(spec_path).read_text())
    assign_partitions(spec)
    validate_spec(spec)
    programs_by_backend = _programs_by_backend(spec)
    site_setup = render_site_setup(programs_by_backend)
    run_dir = Path(run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=False)
    (run_dir / 'tasks').mkdir()
    (run_dir / 'site_setup.sh').write_text(site_setup)
    _atomic_json(run_dir / 'workflow.json', spec)
    state = {'schema': 1,
             'python': sys.executable,
             'spec_sha256': hashlib.sha256((run_dir / 'workflow.json').read_bytes()).hexdigest(),
             'tasks': {}}
    for task in spec['tasks']:
        if not _task_dependencies(task) and not task.get('fallback_on'):
            _stage_task(run_dir, spec, state, task)
    _atomic_json(run_dir / 'state.json', state)
    return run_dir


def _verify_stage_files(run_dir, task, entry):
    """Stop submission or acceptance if prepared inputs have been edited."""
    directory = run_dir / 'tasks' / task['id']
    expected = {'task.json': 'task_sha256', 'job.slurm': 'job_sha256'}
    if task['kind'] == 'external':
        expected[task['input_name']] = 'input_sha256'
    for filename, key in expected.items():
        path = directory / filename
        if not path.is_file() or _file_hash(path) != entry.get(key):
            raise RuntimeError(f"{task['id']}: staged {filename} changed or is missing.")
    record = json.loads((directory/'task.json').read_text())
    if record.get('optimizer_continuation'):
        from kinbot.anl.optimization_recovery import load_seed
        load_seed(directory, record, read(directory/'geometry.xyz'))
    geometry = directory / 'geometry.xyz'
    if (not geometry.is_file() or
            _geometry_hash(read(geometry)) != entry['geometry_sha256']):
        raise RuntimeError(f"{task['id']}: staged geometry changed or is missing.")


def _verify_execution(run_dir, task, entry, execution):
    ident = task['id']
    if (execution.get('schema') != 1 or execution.get('task_id') != ident
            or execution.get('geometry_sha256') != entry['geometry_sha256']
            or execution.get('status') not in ('executed', 'failed')):
        raise RuntimeError(f'{ident}: execution record does not match the task.')
    if (entry.get('job_id') and execution.get('slurm_job_id')
            and str(entry['job_id']) != str(execution['slurm_job_id'])):
        raise RuntimeError(f'{ident}: execution record belongs to a different allocation.')
    if execution['status'] == 'failed':
        return None
    artifacts = execution.get('artifacts')
    if not isinstance(artifacts, dict) or not artifacts:
        raise RuntimeError(f'{ident}: execution has no artifact hashes.')
    directory = run_dir / 'tasks' / ident
    for name, expected in artifacts.items():
        if (not isinstance(name, str) or not _SAFE_NAME.fullmatch(name)
                or not isinstance(expected, str) or
                not re.fullmatch(r'[0-9a-f]{64}', expected)):
            raise RuntimeError(f'{ident}: invalid artifact record.')
        path = directory / name
        if not path.is_file() or _file_hash(path) != expected:
            raise RuntimeError(f'{ident}: artifact {name} changed or is missing.')
    needed = {'task.json', 'geometry.xyz'}
    if task['kind'] == 'external':
        needed.update({task['input_name'], task.get('stdout', 'stdout.txt'),
                       task.get('stderr', 'stderr.txt')})
        needed.update(task['required_outputs'])
        needed.update(task.get('files_from_env', {}))
    else:
        needed.update({'final.xyz', 'optimization.log'})
        details = execution.get('details', {})
        monatomic_identity = details.get('optimizer') == 'identity_monatomic'
        if monatomic_identity:
            if (details.get('calculation_performed') is not False
                    or details.get('generated_files') != []
                    or len(read(directory / 'geometry.xyz')) != 1):
                raise RuntimeError(
                    f'{ident}: invalid monatomic identity geometry record.')
        elif _backend(task) == 'molpro':
            generated = execution.get('details', {}).get('generated_files')
            if (not isinstance(generated, list) or
                    not all(isinstance(name, str) and _SAFE_NAME.fullmatch(name)
                            for name in generated)):
                raise RuntimeError(f'{ident}: Molpro ASE artifact list is missing or invalid.')
            first = f'{ident}_step_0001'
            if not {first + suffix for suffix in ('.inp', '.out', '.log', '.xyz')} <= set(generated):
                raise RuntimeError(f'{ident}: first Molpro force step is missing.')
            needed.update(generated)
    if needed - artifacts.keys():
        raise RuntimeError(f'{ident}: execution is missing required artifact hashes.')
    geometry_output = task.get('geometry_output')
    if geometry_output:
        final = directory / geometry_output
        actual_hash = _geometry_hash(read(final))
        if actual_hash != execution.get('details', {}).get('final_geometry_sha256'):
            raise RuntimeError(f'{ident}: accepted geometry hash mismatch.')
        return actual_hash
    return None


def _run_external(directory, record):
    task = record['task']
    for filename, env_name in task.get('files_from_env', {}).items():
        source = os.environ.get(env_name)
        if not source or not Path(source).is_file():
            raise RuntimeError(f'{filename} requires existing file in ${env_name}.')
        shutil.copyfile(source, directory / filename)
    # Molpro 2024 defaults to disk mode on one node. -m is per process, so
    # reserve 200 MW per MPI process and additional node headroom first.
    molpro_stack_mw = (_molpro_stack_mw(task['resources'])
                       if task['backend'].lower() == 'molpro' else 0)
    command = _external_command(task, molpro_stack_mw)
    scratch_min_mb = task['resources'].get(
        'min_scratch_mb', max(4096, task['resources']['memory_mb'] // 2))
    child_env, runtime = qc_runtime_environment(
        command[0], task['backend'], work_directory=directory,
        scratch_min_mb=scratch_min_mb)
    child_env['OMP_NUM_THREADS'] = str(_omp_threads(task))
    from kinbot.qc_scratch import workspace, molpro_command
    result = None
    native_completed = False
    export_files = ()
    try:
        inputs = [task['input_name'], 'geometry.xyz', *task.get('files_from_env', {})]
        outputs = [*task['required_outputs'], task.get('stdout', 'stdout.txt'),
                   task.get('stderr', 'stderr.txt')]
        outputs.extend(name for alternatives in task.get('output_aliases', {}).values()
                       for name in alternatives)
        if task['backend'].lower() == 'mrcc':
            restart = task.get('restart', {})
            if restart.get('kind') == 'mrcc_amplitudes':
                if (restart.get('checkpoint') != 'fort.16'
                        or _file_hash(directory/'fort.16') != restart.get('sha256')):
                    raise RuntimeError('MRCC restart checkpoint changed or is not hash-bound.')
                inputs.append('fort.16')
                if restart.get('fort17_sha256'):
                    if _file_hash(directory/'fort.17') != restart['fort17_sha256']:
                        raise RuntimeError('MRCC auxiliary restart checkpoint changed.')
                    inputs.append('fort.17')
        if task.get('gaussian_orbital_export'):
            from kinbot.vrc_orbitals import GAUSSIAN_EXPORT_FILES
            export_files = GAUSSIAN_EXPORT_FILES
            outputs.extend(export_files)
        with workspace(directory, runtime, backend=task['backend'], inputs=inputs,
                       outputs=outputs, label=task['id']) as local:
            if task['backend'].lower() == 'molpro':
                command = molpro_command(command, runtime, local)
            stdin = (local / task['input_name']).open('rb') if task.get('stdin') else None
            try:
                with (local / task.get('stdout', 'stdout.txt')).open('wb') as stdout, \
                        (local / task.get('stderr', 'stderr.txt')).open('wb') as stderr:
                    result = subprocess.run(command, cwd=local, stdin=stdin, stdout=stdout,
                                            stderr=stderr, check=False, env=child_env)
                    if result.returncode and task['backend'].lower() == 'mrcc':
                        # Retain once after a handled native failure for manual
                        # review. Never mirror large amplitudes during a run or
                        # add them to every successful job's shared outputs.
                        outputs.extend(['fort.16', 'fort.17'])
            except OSError as exc:
                raise CalculationUnavailableError(f'Could not start native program: {exc}') from exc
            finally:
                if stdin:
                    stdin.close()
            if not result.returncode and task.get('gaussian_orbital_export'):
                from kinbot.vrc_orbitals import export_gaussian_cube
                export_gaussian_cube(local, task['gaussian_orbital_export'], child_env)
        native_completed = not result.returncode
    finally:
        if not native_completed:
            from kinbot.anl.runtime import preserve_qc_failure
            preserve_qc_failure(runtime, directory, Path(task['input_name']).stem,
                                None if result is None else result.returncode)
        cleanup_qc_runtime(runtime)
    orbital_report = None
    if task['backend'].lower() == 'molpro':
        from kinbot.molpro_diagnostics import write_native_report
        orbital_report = write_native_report(
            (directory / task['input_name']).with_suffix('.out'),
            directory / task['input_name'])
    if result.returncode:
        error_file = directory / task.get('stderr', 'stderr.txt')
        detail = error_file.read_text(errors='replace').strip()[-1200:]
        suffix = f' Last stderr: {detail}' if detail else ''
        raise CalculationUnavailableError(
            f'Program exited with status {result.returncode}.{suffix}')
    aliases = _normalize_output_aliases(directory, task)
    _validate_external_outputs(directory, task)
    details = {'command': command, 'returncode': result.returncode, **runtime}
    if export_files:
        details['generated_files'] = [name for name in export_files
                                      if (directory / name).is_file()]
    if orbital_report is not None:
        details['orbital_audit'] = orbital_report.name
        details['generated_files'] = [orbital_report.name]
    if aliases:
        details['output_aliases'] = aliases
        details.setdefault('generated_files', []).extend(aliases.values())
    if task.get('result_parser'):
        from kinbot.anl.results import parse_result
        requested = task['result_parser']
        details['parsed_result'] = parse_result(
            (directory / requested['file']).read_text(errors='replace'), requested)
    return details


def _external_command(task, molpro_stack_mw=None):
    if molpro_stack_mw is None:
        molpro_stack_mw = (_molpro_stack_mw(task['resources'])
                           if task['backend'].lower() == 'molpro' else 0)
    return [arg.replace('{cores}', str(task['resources']['cores']))
            .replace('{input}', task['input_name'])
            .replace('{molpro_stack_mw}', str(molpro_stack_mw))
            for arg in task['command']]


def _normalize_output_aliases(directory, task):
    """Retain a native version's filename and hash-bound canonical copy."""
    selected = {}
    for target, alternatives in task.get('output_aliases', {}).items():
        present = [name for name in (target, *alternatives)
                   if (directory / name).is_file() and (directory / name).stat().st_size]
        if len(present) > 1:
            raise RuntimeError(f'Ambiguous native output files for {target}: {present}')
        if present and present[0] != target:
            shutil.copy2(directory / present[0], directory / target)
            selected[target] = present[0]
    return selected


def _validate_external_outputs(directory, task):
    missing = [name for name in task['required_outputs']
               if not (directory / name).is_file()
               or not (directory / name).stat().st_size]
    if missing:
        raise RuntimeError(f'Missing or empty required output: {missing}')
    if task.get('gaussian_orbital_export'):
        from kinbot.vrc_orbitals import validate_cube_geometry
        from kinbot.orbital_cube import read_orbital_cube
        settings = task['gaussian_orbital_export']
        text = (directory/settings['log']).read_text(errors='replace')
        if re.search(r'Convergence criterion not met|Convergence failure|Error termination', text, re.I):
            raise RuntimeError('Gaussian orbital calculation did not converge normally.')
        validate_cube_geometry(read_orbital_cube(directory/settings['cube']), read(directory/'geometry.xyz'))
    if task['backend'].lower() == 'molpro':
        from kinbot.ase_modules.calculators.molpro import check_process_count
        check_process_count(directory / Path(task['input_name']).with_suffix('.out').name,
                            task['resources']['cores'])
    if task.get('success_marker'):
        marker = task['success_marker']
        output = directory / marker['file']
        if not output.is_file() or marker['contains'] not in output.read_text(errors='replace'):
            raise RuntimeError(f"Success marker absent from {marker['file']}.")
    for marker in task.get('failure_markers', []):
        output = directory / marker['file']
        if output.is_file() and marker['contains'] in output.read_text(errors='replace'):
            raise RuntimeError(f"Failure marker present in {marker['file']}: "
                               f"{marker['contains']}")


def _native_success_recovery_candidate(task, execution):
    """Whether a failed task can be checked for native-output recovery.

    Eligibility deliberately does not depend on where the old Python stack
    classified the failure.  Some site launchers return a nonzero status even
    though the native program has written its unambiguous normal-termination
    marker and complete result.  The recovery itself still verifies the
    staged input and geometry, required output files, native success and
    failure markers, and the current result parser before accepting anything.
    """
    return (task.get('kind') == 'external'
            and bool(task.get('result_parser'))
            and bool(task.get('success_marker'))
            and execution.get('status') == 'failed')


def _native_success_present(directory, task):
    """Cheaply identify a native normal termination before full validation."""
    marker = task.get('success_marker')
    if not marker:
        return False
    output = directory / marker['file']
    if not output.is_file() or not output.stat().st_size:
        return False
    text_by_file = {marker['file']: output.read_text(errors='replace')}
    if marker['contains'] not in text_by_file[marker['file']]:
        return False
    for failure in task.get('failure_markers', []):
        name = failure['file']
        path = directory / name
        if name not in text_by_file and path.is_file():
            text_by_file[name] = path.read_text(errors='replace')
        if failure['contains'] in text_by_file.get(name, ''):
            return False
    return True


def _recorded_failed_returncode(execution):
    """Recover a subprocess status from an old failed execution, if present."""
    match = re.search(r'Program exited with status\s+(\d+)',
                      execution.get('error', ''))
    return int(match.group(1)) if match else None


def _recover_one_electron_f12_nan(run_dir, spec, state, task, entry,
                                  previous):
    """Accept Molpro's completed one-electron F12 value before its 0/0 abort."""
    parser = task.get('result_parser', {})
    molecule = spec.get('molecule', {})
    atoms = _atoms(molecule)
    electrons = int(sum(atoms.numbers)) - molecule.get('charge', 0)
    if (electrons != 1 or task.get('kind') != 'external'
            or task.get('backend', '').lower() != 'molpro'
            or parser.get('kind') != 'molpro_energy'
            or parser.get('method') != 'CCSD(T)-F12b'
            or entry.get('status') != 'failed'
            or previous.get('status') != 'failed'):
        raise ValueError('Task is not a failed one-electron Molpro F12 job.')
    if entry.get('job_id') and _job_active(entry['job_id']):
        raise RuntimeError(f"{task['id']}: Slurm job is still active.")
    _verify_stage_files(run_dir, task, entry)
    directory = run_dir / 'tasks' / task['id']
    native = directory / parser['file']
    if not native.is_file():
        raise ValueError('One-electron Molpro F12 output is missing.')
    replacement_parser = {**parser, 'rank_exact_electrons': 1}
    from kinbot.anl.results import parse_result
    parsed = parse_result(native.read_text(errors='replace'),
                          replacement_parser)
    if not parsed.get('rank_exact', {}).get(
            'recovered_molpro_scale_trip_nan'):
        raise ValueError('Molpro output is not the exact scale-trip NaN case.')

    replacement = deepcopy(task)
    replacement['result_parser'] = replacement_parser
    replacement['recovery'] = {
        'kind': 'one_electron_f12_scale_trip_nan',
        'calculation_changed': False,
        'reason': ('Molpro evaluated the CABS reference and exact-zero F12 '
                   'correlation before an undefined zero-over-zero triples '
                   'scale factor aborted the process'),
    }
    candidate = deepcopy(spec)
    candidate['tasks'][candidate['tasks'].index(task)] = replacement
    validate_spec(candidate)

    previous_name = 'execution.failed.json'
    previous_path = directory / previous_name
    if previous_path.exists():
        raise RuntimeError(
            f"{task['id']}: previous failed execution archive already exists.")
    names = {'geometry.xyz', 'task.json', task['input_name'], parser['file'],
             task.get('stdout', 'stdout.txt'), task.get('stderr', 'stderr.txt'),
             previous_name}
    missing = sorted(name for name in names - {previous_name}
                     if not (directory / name).is_file())
    if missing:
        raise RuntimeError(
            f"{task['id']}: recovery artifacts are missing: {missing}")
    _atomic_json(previous_path, previous)
    task.clear()
    task.update(replacement)
    _atomic_json(run_dir / 'workflow.json', spec)
    state['spec_sha256'] = _file_hash(run_dir / 'workflow.json')
    result = {
        'schema': 1, 'task_id': task['id'],
        'geometry_sha256': entry['geometry_sha256'], 'status': 'executed',
        'details': {
            'command': _external_command(task),
            'returncode': 255,
            'returncode_source': 'preserved_molpro_scale_trip_nan',
            'parsed_result': parsed,
            'recovered_without_execution': True,
            'previous_error': previous.get('error', ''),
        },
        'artifacts': {name: _file_hash(directory / name)
                      for name in sorted(names)},
    }
    _verify_execution(run_dir, task, entry, result)
    _atomic_json(directory / 'execution.json', result)
    entry.update(status='complete', execution='executed',
                 automatically_recovered_rank_exact=True)
    entry.pop('error', None)
    entry.pop('job_id', None)
    return result


def _reparse_execution(run_dir, task, entry, previous):
    """Accept one fully checked native success without executing QC again."""
    ident = task['id']
    directory = run_dir / 'tasks' / ident
    outcome = directory / 'execution.json'
    if not _native_success_recovery_candidate(task, previous):
        raise ValueError(
            f'{ident}: task is not eligible for native-success recovery.')
    _verify_execution(run_dir, task, entry, previous)
    _verify_stage_files(run_dir, task, entry)
    record = json.loads((directory / 'task.json').read_text())
    if record.get('geometry_sha256') != entry.get('geometry_sha256'):
        raise RuntimeError(f'{ident}: staged geometry provenance changed.')
    previous_name = 'execution.failed.json'
    previous_path = directory / previous_name
    if previous_path.exists():
        raise RuntimeError(f'{ident}: previous failed execution archive already exists.')
    names = {'geometry.xyz', 'task.json', task['input_name'],
             task.get('stdout', 'stdout.txt'), task.get('stderr', 'stderr.txt')}
    names.update(task['required_outputs'])
    names.update(task.get('files_from_env', {}))
    missing = sorted(name for name in names if not (directory / name).is_file())
    if missing:
        raise RuntimeError(f'{ident}: execution artifacts are missing: {missing}')
    _validate_external_outputs(directory, task)
    requested = task['result_parser']
    from kinbot.anl.results import parse_result
    parsed = parse_result(
        (directory / requested['file']).read_text(errors='replace'), requested)
    _atomic_json(previous_path, previous)
    failed_returncode = _recorded_failed_returncode(previous)
    parser_failure = 'parse_result(' in previous.get('traceback', '')
    details = {
        'command': _external_command(task),
        'returncode': (0 if parser_failure else failed_returncode),
        'returncode_source': (
            'original_parser_path_and_native_success_marker'
            if parser_failure else
            'preserved_dispatch_failure_and_native_success_marker'),
        'parsed_result': parsed,
        'reparsed_without_execution': True,
        'native_success_recovered': True,
        'previous_error': previous.get('error', ''),
    }
    names.add(previous_name)
    result = {
        'schema': 1,
        'task_id': ident,
        'geometry_sha256': record['geometry_sha256'],
        'status': 'executed',
        'details': details,
        'artifacts': {name: _file_hash(directory / name)
                      for name in sorted(names)},
    }
    actual_hash = _verify_execution(run_dir, task, entry, result)
    _atomic_json(outcome, result)
    return result, actual_hash


def reparse_failed(run_dir, ident):
    """Recover a verified native success without executing QC again."""
    run_dir, spec, state = _load(run_dir)
    by_id = {task['id']: task for task in spec['tasks']}
    task = by_id.get(ident)
    entry = state['tasks'].get(ident)
    if task is None or entry is None:
        raise ValueError(f'{ident}: task is unavailable.')
    if task['kind'] != 'external' or not task.get('result_parser'):
        raise ValueError(f'{ident}: reparse requires an external parsed task.')
    if not task.get('success_marker'):
        raise ValueError(f'{ident}: reparse requires a native success marker.')
    directory = run_dir / 'tasks' / ident
    outcome = directory / 'execution.json'
    if not outcome.is_file():
        raise RuntimeError(f'{ident}: failed execution record is missing.')
    previous = json.loads(outcome.read_text())
    if entry.get('status') == 'complete':
        _verify_stage_files(run_dir, task, entry)
        _verify_execution(run_dir, task, entry, previous)
        if previous.get('status') != 'executed':
            raise RuntimeError(f'{ident}: complete task has no executed result.')
        requested = task['result_parser']
        from kinbot.anl.results import parse_result, parsed_result_extends
        parsed = parse_result(
            (directory / requested['file']).read_text(errors='replace'),
            requested)
        if not parsed_result_extends(
                previous.get('details', {}).get('parsed_result'), parsed):
            raise RuntimeError(f'{ident}: saved parsed result differs from native output.')
        return previous
    if entry.get('status') != 'failed':
        raise ValueError(f'{ident}: only a failed or complete task can be reparsed.')
    result, actual_hash = _reparse_execution(
        run_dir, task, entry, previous)
    entry.update(status='complete', execution='executed')
    entry.pop('error', None)
    if actual_hash is not None:
        entry['final_geometry_sha256'] = actual_hash
    _atomic_json(run_dir / 'state.json', state)
    return result


def _run_ase_optimize(directory, record):
    task = record['task']
    atoms = read(directory / 'geometry.xyz')
    if len(atoms) == 1:
        # A monatomic state has no internal or Cartesian geometry degree of
        # freedom after removal of translation.  Calling an electronic
        # structure calculator here adds no geometry information and several
        # otherwise valid gradient implementations reject the empty problem.
        # Preserve the exact accepted position and record the identity step;
        # the independent energy/correction tasks still run normally.
        (directory / 'optimization.log').write_text(
            'Monatomic state: geometry optimization is an exact identity; '
            'no calculator was invoked.\n')
        write(directory / 'final.xyz', atoms)
        return {
            'optimizer': 'identity_monatomic',
            'calculation_performed': False,
            'generated_files': [],
        }
    from ase.optimize import BFGS
    from sella import Sella, Constraints
    from kinbot.ase_modules.calculators.factory import build_calculator

    initial_positions = atoms.positions.copy()
    from kinbot.anl.optimization_recovery import load_seed
    used_chunks = load_seed(directory, record, atoms)
    # UMA/OMol reads these from Atoms.info; XYZ does not preserve them.
    atoms.info['charge'] = record['molecule']['charge']
    atoms.info['spin'] = record['molecule']['multiplicity']
    profile = _runtime_profile(task)
    atoms.calc = build_calculator(profile, directory, {
        'name': task['id'], 'charge': record['molecule']['charge'],
        'mult': record['molecule']['multiplicity']})
    options = task.get('optimizer', {})
    common = {'trajectory': str(directory / 'optimization.traj'),
              'logfile': str(directory / 'optimization.log')}

    def optimize():
        if len(atoms) > 2 and profile.optimizer == 'sella':
            specification = options.get('constraints', {})
            constraints = None
            if specification:
                constraints = Constraints(atoms)
                for pair in specification.get('fixed_bonds', ()):
                    if max(pair) >= len(atoms):
                        raise ValueError('Fixed-bond atom index is out of range.')
                    constraints.fix_bond(tuple(pair), target=float(np.linalg.norm(
                        initial_positions[pair[0]] - initial_positions[pair[1]])))
                for first, second, distance in specification.get(
                        'minimum_bonds', ()):
                    if max(first, second) >= len(atoms):
                        raise ValueError('Minimum-bond atom index is out of range.')
                    constraints.fix_bond(
                        (first, second), target=float(distance), comparator='gt')
            optimizer = Sella(
                atoms, order=options.get('order', 0), constraints=constraints,
                **common, **options.get('sella_kwargs', {}))
        else:
            optimizer = BFGS(atoms, **common)
        from kinbot.optimization_budget import run_optimizer, step_retries
        return run_optimizer(optimizer, fmax=options.get('fmax', 0.03),
                             steps=options.get('steps', 100),
                             chunks=max(1, 1 + options.get('step_retries',
                                 step_retries(directory)) - used_chunks),
                             used=used_chunks, report=directory/'optimization_budget.json')

    recovery_files = []
    if 'minimum_recovery' in options:
        from kinbot.minimum_recovery import optimize_anl_minimum
        recovery_files = optimize_anl_minimum(
            atoms, optimize, directory, options['minimum_recovery'])
    elif not optimize():
        raise RuntimeError('ASE geometry optimization did not converge.')
    energy_ev = float(atoms.get_potential_energy())
    write(directory / 'final.xyz', atoms)
    details = {
        'energy_ev': energy_ev, 'energy_hartree': energy_ev / Hartree,
        'optimizer': profile.optimizer if len(atoms) > 2 else 'bfgs_or_atom',
        'order': options.get('order', 0),
        'generated_files': [*getattr(atoms.calc, 'generated_files', []),
                            *recovery_files],
    }
    checkpoint = getattr(atoms.calc, 'orbital_checkpoint_model', None)
    if checkpoint is not None:
        details['orbital_checkpoint_model'] = checkpoint
    return details


def run_task(task_file):
    task_file = Path(task_file).resolve()
    directory = task_file.parent
    record = json.loads(task_file.read_text())
    task = record['task']
    outcome = directory / 'execution.json'
    if outcome.exists():
        raise RuntimeError('Task already has execution.json; review before rerunning.')
    result = {'schema': 1, 'task_id': task['id'],
              'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
              'execution_node': os.environ.get('SLURMD_NODENAME', socket.gethostname()),
              'geometry_sha256': record['geometry_sha256']}
    try:
        from kinbot.anl.live_control import worker_load
        run_dir, spec, state = worker_load(directory.parents[1], task['id'])
        task_spec = next(item for item in spec['tasks'] if item['id'] == task['id'])
        if task != task_spec or task['id'] != directory.name:
            raise RuntimeError('Task record does not match prepared workflow.')
        _verify_stage_files(run_dir, task, state['tasks'][task['id']])
        atoms = read(directory / 'geometry.xyz')
        if _geometry_hash(atoms) != record['geometry_sha256']:
            raise RuntimeError('Staged geometry changed after preparation.')
        if task['kind'] == 'external':
            details = _run_external(directory, record)
        else:
            try:
                details = _run_ase_optimize(directory, record)
            except Exception as exc:
                raise CalculationUnavailableError(str(exc)) from exc
        if task.get('geometry_output'):
            final = read(directory / task['geometry_output'])
            if (final.get_chemical_symbols() != record['molecule']['symbols']
                    or not np.isfinite(final.positions).all()):
                raise RuntimeError('Final geometry has changed atom identity/order or is nonfinite.')
            details['final_geometry_sha256'] = _geometry_hash(final)
        names = {'geometry.xyz', 'task.json'}
        if task['kind'] == 'external':
            names |= {task['input_name'], task.get('stdout', 'stdout.txt'),
                      task.get('stderr', 'stderr.txt')}
            names.update(task['required_outputs'])
            names.update(task.get('files_from_env', {}))
            names.update(details.get('generated_files', []))
        else:
            names |= {'final.xyz', 'optimization.log', 'continuation.xyz', 'optimization_budget.json'}
            if (directory / 'optimization.traj').is_file():
                names.add('optimization.traj')
            names.update(_basename(name, 'calculator artifact')
                         for name in details.get('generated_files', []))
        artifacts = {name: _file_hash(directory / name)
                     for name in sorted(names) if (directory / name).is_file()}
        result.update(status='executed', details=details, artifacts=artifacts)
    except Exception as exc:
        result.update(status='failed', error=str(exc),
                      traceback=traceback.format_exc())
        result['failure_phase'] = (
            'execution' if isinstance(exc, CalculationUnavailableError)
            else 'validation')
        from kinbot.qc_scratch import ScratchUnavailableError
        cause = exc
        while cause is not None:
            if isinstance(cause, ScratchUnavailableError):
                result.update(failure_phase='scratch_setup',
                              retryable_scratch_node=cause.retryable_node)
                break
            cause = cause.__cause__
    _atomic_json(outcome, result)
    if result['status'] != 'executed':
        raise RuntimeError(f"{task['id']}: {result['error']}")
    return result


def _monatomic_geometry_task(spec, task):
    return (task.get('kind') == 'ase_optimize'
            and len(spec.get('molecule', {}).get('symbols', ())) == 1)


def _complete_monatomic_geometry(run_dir, spec, task, entry, previous=None):
    """Complete an atomic geometry node as an exact identity operation."""
    ident = task['id']
    if not _monatomic_geometry_task(spec, task):
        raise ValueError(f'{ident}: task is not a monatomic geometry operation.')
    _verify_stage_files(run_dir, task, entry)
    directory = run_dir / 'tasks' / ident
    record = json.loads((directory / 'task.json').read_text())
    if record.get('geometry_sha256') != entry.get('geometry_sha256'):
        raise RuntimeError(f'{ident}: staged geometry provenance changed.')
    if previous is not None and previous.get('status') != 'failed':
        raise RuntimeError(f'{ident}: only a failed atomic attempt can be recovered.')

    details = _run_ase_optimize(directory, record)
    final_hash = _geometry_hash(read(directory / 'final.xyz'))
    if final_hash != entry['geometry_sha256']:
        raise RuntimeError(f'{ident}: monatomic identity changed the geometry.')
    details['final_geometry_sha256'] = final_hash
    names = {'geometry.xyz', 'task.json', 'final.xyz', 'optimization.log'}
    if previous is not None:
        previous_name = 'execution.failed.json'
        previous_path = directory / previous_name
        if previous_path.exists():
            raise RuntimeError(
                f'{ident}: previous failed execution archive already exists.')
        _atomic_json(previous_path, previous)
        names.add(previous_name)
        details.update(
            recovered_without_qc_execution=True,
            previous_error=previous.get('error', ''))
    result = {
        'schema': 1,
        'task_id': ident,
        'geometry_sha256': record['geometry_sha256'],
        'status': 'executed',
        'details': details,
        'artifacts': {name: _file_hash(directory / name)
                      for name in sorted(names)},
    }
    _verify_execution(run_dir, task, entry, result)
    _atomic_json(directory / 'execution.json', result)
    entry.update(status='complete', execution='executed',
                 final_geometry_sha256=final_hash)
    entry.pop('error', None)
    entry.pop('job_id', None)
    return result


def _job_active(job_id):
    for attempt in range(_SCHEDULER_QUERY_ATTEMPTS):
        try:
            result = subprocess.run(
                ['squeue', '-h', '-j', str(job_id), '-o', '%i'],
                capture_output=True, text=True, check=False,
                timeout=_SCHEDULER_TIMEOUT)
            break
        except subprocess.TimeoutExpired as error:
            if attempt + 1 == _SCHEDULER_QUERY_ATTEMPTS:
                raise SchedulerUnavailableError(
                    f'squeue timed out {_SCHEDULER_QUERY_ATTEMPTS} times '
                    f'while checking job {job_id}.') from error
            time.sleep(5)
        except FileNotFoundError as error:
            raise RuntimeError('squeue is unavailable on PATH.') from error
        except OSError as error:
            if attempt + 1 == _SCHEDULER_QUERY_ATTEMPTS:
                raise SchedulerUnavailableError(
                    f'squeue could not start {_SCHEDULER_QUERY_ATTEMPTS} '
                    f'times while checking job {job_id}: {error}') from error
            time.sleep(5)
    if result.returncode:
        # A completed job can disappear from slurmctld before its execution
        # record is reconciled. Some Slurm versions return an error rather
        # than an empty queue for this specific case.
        if re.search(r'\binvalid job id specified\b', result.stderr, re.I):
            return False
        raise SchedulerUnavailableError(
            f'squeue failed: {result.stderr.strip()}')
    return str(job_id) in result.stdout.split()


def _refresh_generated_site_setup(run_dir, backend, programs):
    """Append a verified setup overlay when site configuration appeared late.

    Prepared graphs remain reusable when a module, PATH entry, or advertised
    backend root becomes available after ``prepare``.  Existing setup is
    preserved verbatim; the overlay is appended only after a shell probe proves
    that it resolves every missing executable for the affected backend.
    """
    path = run_dir / 'site_setup.sh'
    current = path.read_text()
    generated_marker = '# Generated from the environment visible during prepare.'
    if generated_marker not in current:
        return False
    overlay = render_site_setup({backend: set(programs)})
    digest = hashlib.sha256(overlay.encode()).hexdigest()
    marker = f'# KinBot preflight site refresh {digest}'
    if marker in current:
        return False
    temporary = run_dir / '.site_setup.refresh.sh'
    temporary.write_text(overlay)
    lines = [
        'set -euo pipefail',
        f'export KINBOT_BACKEND={shlex.quote(backend)}',
        f'source {shlex.quote(str(path))}',
        f'source {shlex.quote(str(temporary))}',
    ]
    for program in sorted(programs):
        lines.append(f'command -v {shlex.quote(program)} >/dev/null')
    try:
        result = subprocess.run(
            ['bash', '-c', '\n'.join(lines)], cwd=run_dir,
            capture_output=True, text=True, check=False)
    finally:
        temporary.unlink(missing_ok=True)
    if result.returncode:
        return False
    with path.open('a') as stream:
        stream.write(f'\n{marker}\n{overlay}')
    return True


def _mark_preflight_unavailable(spec, state, ident, reason):
    """Record one optional environment failure and its unsent dependants."""
    by_id = {task['id']: task for task in spec['tasks']}
    task = by_id[ident]
    if task.get('failure_policy') != 'skip_component':
        raise ValueError(f'{ident}: preflight failure is not optional.')
    entry = dict(state['tasks'].get(ident, {}))
    if entry.get('status') in ('submitted', 'complete'):
        raise RuntimeError(
            f'{ident}: cannot make a submitted or complete task unavailable.')
    entry.update({
        'status': 'unavailable', 'failure_phase': 'preflight',
        'unavailable_reason': reason,
        'component_families': list(task.get('component_families', ())),
    })
    entry.pop('job_id', None)
    state['tasks'][ident] = entry

    queue = [(target, ident) for target in task.get('failure_skips', ())]
    while queue:
        target, source = queue.pop(0)
        if target not in by_id:
            continue
        current = state['tasks'].get(target)
        if current is not None and current.get('status') in (
                'submitted', 'complete', 'unavailable', 'failed'):
            continue
        target_task = by_id[target]
        state['tasks'][target] = {
            'status': 'skipped',
            'geometry_from': target_task.get('geometry_from', 'initial'),
            'skipped_by': [source],
            'skip_reason': 'requested component is unavailable',
            'component_families': list(
                target_task.get('component_families', ())),
        }
        queue.extend((child, target)
                     for child in target_task.get('failure_skips', ()))

    changed = True
    while changed:
        changed = False
        for dependent in spec['tasks']:
            target = dependent['id']
            current = state['tasks'].get(target)
            if current is not None and current.get('status') in (
                    'submitted', 'complete', 'unavailable', 'failed', 'skipped'):
                continue
            unavailable = [
                dependency for dependency in _task_dependencies(
                    dependent, state)
                if state['tasks'].get(dependency, {}).get('status') in
                ('unavailable', 'skipped')]
            if not unavailable:
                continue
            state['tasks'][target] = {
                'status': 'skipped',
                'geometry_from': dependent.get('geometry_from', 'initial'),
                'skipped_by': unavailable,
                'skip_reason': 'dependency is unavailable',
            }
            changed = True


def preflight(run_dir, _site_refresh_attempts=None, *, skip_optional=False,
              progress=None):
    """Check each task, reusing identical successful environment probes.

    Reuse lasts only for this call and this graph's working directory. Input
    hashes, geometry, resource layout and Slurm tests are always task-specific.
    No success is reused across graphs, scheduler retries or later launches.
    """
    from kinbot.pes_phase import checkpoint
    checkpoint()
    from kinbot.anl import live_control
    live_control.resume(run_dir)
    from kinbot.anl import account_recovery, scratch_recovery, timeout_recovery, lock_recovery, optimization_recovery
    account_recovery.resume(run_dir)
    scratch_recovery.resume(run_dir)
    timeout_recovery.resume(run_dir)
    lock_recovery.resume(run_dir)
    optimization_recovery.resume(run_dir)
    run_dir, spec, state = _load(run_dir)
    refresh_attempts = (set() if _site_refresh_attempts is None
                        else _site_refresh_attempts)
    for name in ('sbatch', 'squeue', 'bash'):
        if not shutil.which(name):
            raise RuntimeError(f'Preflight: {name} is unavailable on PATH.')
    programs = set()
    env_files = set()
    environments = set()
    environment_checks = 0
    environment_reuses = 0
    runtime_environment_deferred = []
    for task in spec['tasks']:
        entry = state['tasks'].get(task['id'])
        if entry and entry.get('status') in ('unavailable', 'skipped'):
            continue
        if entry and entry['status'] in ('staged', 'submitted', 'complete'):
            _verify_stage_files(run_dir, task, entry)
            script = (run_dir / 'tasks' / task['id'] / 'job.slurm').read_text()
            if '#SBATCH --exclusive\n' not in script:
                raise RuntimeError(f"{task['id']}: job does not request an exclusive node.")
            if entry['status'] == 'staged':
                ranks = task['resources']['cores'] if _backend(task) == 'molpro' else 1
                cpus = 1 if _backend(task) == 'molpro' else task['resources']['cores']
                if (f'#SBATCH --ntasks={ranks}\n' not in script
                        or f'#SBATCH --cpus-per-task={cpus}\n' not in script):
                    raise RuntimeError(f"{task['id']}: Slurm task/CPU layout does not match the backend.")
        task_programs = set()
        task_files = set()
        if task['kind'] == 'external':
            task_programs.add(task['command'][0])
            task_programs.update(task.get('required_executables', []))
            task_files.update(task.get('files_from_env', {}).values())
        elif _backend(task) in ('gaussian', 'gauss', 'molpro'):
            backend = _backend(task)
            parts = shlex.split(_runtime_profile(task).command or
                                ('molpro' if backend == 'molpro' else 'g16'))
            if not parts:
                raise ValueError(f"{task['id']}: calculator command is empty.")
            task_programs.add(parts[0])
        programs.update(task_programs)
        env_files.update(task_files)
        if task.get('runtime_setup'):
            # These are explicit operator commands for the allocation only.
            # Login-host probes cannot validate that changed environment.
            runtime_environment_deferred.append(task['id'])
            continue
        lines = [
            'set -euo pipefail',
            f"export OMP_NUM_THREADS={_omp_threads(task)}",
            f"export KINBOT_BACKEND={shlex.quote(_backend(task))}",
            f'source {shlex.quote(str(run_dir / "site_setup.sh"))}',
            *task.get('setup', []),
            f"export OMP_NUM_THREADS={_omp_threads(task)}",
        ]
        if _backend(task) == 'molpro':
            lines.append('export MKL_NUM_THREADS=1')
        elif _backend(task) == 'mrcc':
            lines.append(f'export MKL_NUM_THREADS={_omp_threads(task)}')
        for program in sorted(task_programs):
            lines.append(f'command -v {shlex.quote(program)} >/dev/null || '
                         f'{{ echo {shlex.quote("Missing executable: " + program)} >&2; exit 1; }}')
            if _backend(task) == 'cfour':
                lines.append(f'{shlex.quote(state["python"])} -c '
                             + shlex.quote('from kinbot.anl.runtime import '
                                           'qc_runtime_environment; '
                                           f'qc_runtime_environment({program!r}, "cfour")'))
        for env_name in sorted(task_files):
            lines.append(f'test -f "${{{env_name}:-}}" || '
                         f'{{ echo {shlex.quote("Missing file in $" + env_name)} >&2; exit 1; }}')
        python = shlex.quote(state['python'])
        lines.append(f'{python} -c '
                     + shlex.quote('import ase, sella; import kinbot.anl.dispatch'))
        shell = '\n'.join(lines)
        # The shell captures backend, executable/file requirements, task setup,
        # interpreter and thread counts. Include the site file and inherited
        # environment as well; a changed setup must never inherit an old pass.
        key = (shell, _file_hash(run_dir / 'site_setup.sh'),
               tuple(sorted(os.environ.items())))
        if key in environments:
            environment_reuses += 1
            continue
        environment_checks += 1
        if progress is not None:
            progress({'phase': 'checking_environment', 'task': task['id'],
                      'backend': _backend(task),
                      'environment_checks': environment_checks,
                      'environment_checks_reused': environment_reuses})
        result = subprocess.run(['bash', '-c', shell], cwd=run_dir,
                                capture_output=True, text=True, check=False)
        if result.returncode:
            detail = (result.stderr.strip() or result.stdout.strip()
                      or f'exit status {result.returncode} without diagnostics')
            refresh_key = (_backend(task), tuple(sorted(task_programs)))
            if ('Missing executable:' in detail
                    and refresh_key not in refresh_attempts):
                refresh_attempts.add(refresh_key)
                if _refresh_generated_site_setup(
                        run_dir, refresh_key[0], refresh_key[1]):
                    refreshed = preflight(
                        run_dir, _site_refresh_attempts=refresh_attempts,
                        skip_optional=skip_optional, progress=progress)
                    refreshed['site_setup_refreshed'] = True
                    return refreshed
            if (skip_optional
                    and task.get('failure_policy') == 'skip_component'):
                _mark_preflight_unavailable(
                    spec, state, task['id'], detail)
                _atomic_json(run_dir / 'state.json', state)
                continue
            raise RuntimeError(f"{task['id']}: preflight failed after site setup: "
                               + detail)
        environments.add(key)
    checked = 0
    for task in spec['tasks']:
        entry = state['tasks'].get(task['id'])
        if not entry or entry['status'] != 'staged':
            continue
        if progress is not None:
            progress({'phase': 'checking_slurm_script', 'task': task['id'],
                      'slurm_scripts_tested': checked,
                      'environment_checks': environment_checks,
                      'environment_checks_reused': environment_reuses})
        try:
            result = subprocess.run(
                ['sbatch', '--test-only', 'job.slurm'],
                cwd=run_dir / 'tasks' / task['id'], capture_output=True,
                text=True, check=False, timeout=_SCHEDULER_TIMEOUT)
        except subprocess.TimeoutExpired as error:
            raise SchedulerUnavailableError(
                f'{task["id"]}: Slurm test-only submission timed out; no '
                'job was created.') from error
        except FileNotFoundError as error:
            raise RuntimeError(
                'Preflight: sbatch is unavailable on PATH.') from error
        except OSError as error:
            raise SchedulerUnavailableError(
                f'{task["id"]}: Slurm test-only submission could not start: '
                f'{error}') from error
        if result.returncode:
            detail = result.stderr.strip() or result.stdout.strip()
            if _scheduler_unavailable(detail):
                raise SchedulerUnavailableError(
                    f'{task["id"]}: Slurm test-only submission is '
                    f'unavailable: {detail}')
            if (skip_optional
                    and task.get('failure_policy') == 'skip_component'):
                _mark_preflight_unavailable(
                    spec, state, task['id'], detail)
                _atomic_json(run_dir / 'state.json', state)
                continue
            raise RuntimeError(f"{task['id']}: Slurm rejected test-only batch script: "
                               + detail)
        checked += 1
    return {'programs': sorted(programs), 'file_variables': sorted(env_files),
            'exclusive_jobs_checked': len(state['tasks']),
            'environment_checks': environment_checks,
            'environment_checks_reused': environment_reuses,
            'runtime_environment_deferred': runtime_environment_deferred,
            'slurm_scripts_tested': checked}


def _clear_blocked_tasks(state):
    """Discard derived blocks/skips so a retried graph can be reevaluated."""
    for task_id, entry in list(state['tasks'].items()):
        if entry.get('status') == 'blocked':
            state['tasks'].pop(task_id)
        elif entry.get('status') == 'skipped':
            if entry.get('resume_status') == 'staged':
                entry['status'] = 'staged'
                for key in ('resume_status', 'skipped_by', 'skip_reason',
                            'component_families'):
                    entry.pop(key, None)
            else:
                state['tasks'].pop(task_id)


def retry_failed(run_dir, ident):
    """Archive one failed attempt and restage it without rerunning siblings."""
    run_dir, spec, state = _load(run_dir)
    by_id = {task['id']: task for task in spec['tasks']}
    if (ident not in by_id or state['tasks'].get(ident, {}).get('status')
            not in ('failed', 'unavailable')):
        raise ValueError(
            f'{ident}: only a failed or unavailable task can be retried.')
    entry = state['tasks'][ident]
    if entry.get('job_id') and _job_active(entry['job_id']):
        raise RuntimeError(f'{ident}: Slurm job is still active.')
    if any(other != ident and other in state['tasks']
           and state['tasks'][other].get('status') not in ('blocked', 'skipped')
           and ident in _task_dependencies(task)
           for other, task in by_id.items()):
        raise RuntimeError(f'{ident}: downstream task already exists; review run.')
    old = run_dir / 'tasks' / ident
    attempt = entry.get('attempt', 1)
    archive = run_dir / 'attempts' / ident / str(attempt)
    if archive.exists():
        raise RuntimeError(f'{ident}: retry archive already exists.')
    archive.parent.mkdir(parents=True, exist_ok=True)
    old.replace(archive)
    try:
        _stage_task(run_dir, spec, state, by_id[ident])
    except Exception:
        state['tasks'][ident] = {**entry, 'archive': str(archive)}
        _atomic_json(run_dir / 'state.json', state)
        raise
    state['tasks'][ident]['attempt'] = attempt + 1
    state['tasks'][ident]['previous_attempt'] = str(archive)
    _clear_blocked_tasks(state)
    _atomic_json(run_dir / 'state.json', state)
    return archive


def can_recover_gaussian_vpt2_symmetry(run_dir, task, entry):
    """Read-only eligibility check; recovery still verifies hashes/ownership."""
    if (entry.get('status') not in ('failed', 'unavailable')
            or task.get('kind') != 'external'
            or task.get('backend', '').lower() != 'gaussian'
            or task.get('result_parser', {}).get('kind') != 'gaussian_vpt2'):
        return False
    native = Path(run_dir) / 'tasks' / task['id'] / task['result_parser']['file']
    if not native.is_file():
        return False
    try:
        from kinbot.anl.vpt2_recovery import context
        history, linear = context(run_dir, task, entry)
        _vpt2_framework_recovery(task, native.read_text(errors='replace'),
                                 history=history, linear=linear)
    except VPT2RecoveryUnavailable:
        return False
    return True


def recover_gaussian_vpt2_symmetry(run_dir, ident):
    """Apply the common, history-aware symmetry plan under the graph owner.

    Preserve geometry, theory and validators; archive only the failed task
    and affected blocked descendants. Native successes are reconciled first.
    """
    from kinbot.anl import live_control
    from kinbot.anl.salvage import descendants
    from kinbot.anl.timeout_recovery import _job_id
    from kinbot.anl.vpt2_recovery import context, replace_symmetry, symmetry
    live_control.resume(run_dir)
    run_dir, spec, state = _load(run_dir)
    by_id = {task['id']: task for task in spec['tasks']}
    task = by_id.get(ident)
    entry = state['tasks'].get(ident)
    if task is None or entry is None:
        raise ValueError(f'{ident}: task is unavailable.')
    parser = task.get('result_parser', {})
    if (task.get('kind') != 'external'
            or task.get('backend', '').lower() != 'gaussian'
            or parser.get('kind') != 'gaussian_vpt2'):
        raise ValueError(f'{ident}: recovery requires a Gaussian VPT2 task.')
    if entry.get('status') not in ('failed', 'unavailable'):
        raise ValueError(
            f'{ident}: only a failed or unavailable Gaussian VPT2 task can recover.')
    directory = run_dir / 'tasks' / ident
    job = _job_id(directory, ident, entry)
    if job and _job_active(job):
        raise RuntimeError(f'{ident}: Slurm job is still active.')
    children = descendants(spec, ident)
    if any(state['tasks'].get(other, {}).get('status') not in
           (None, 'waiting', 'blocked', 'skipped') for other in children):
        raise RuntimeError(f'{ident}: downstream task already exists; review run.')
    _verify_stage_files(run_dir, task, entry)
    outcome = directory/'execution.json'
    if outcome.is_file():
        execution = json.loads(outcome.read_text())
        _verify_execution(run_dir, task, entry, execution)
        if execution.get('slurm_job_id') and not job:
            raise VPT2RecoveryUnavailable(f'{ident}: allocation history is ambiguous; reconcile it first.')
        if execution['status'] != 'failed':
            raise RuntimeError(f'{ident}: a native success must be reconciled, not retried.')
    native = directory / parser['file']
    if not native.is_file():
        raise RuntimeError(f'{ident}: native Gaussian output is missing.')
    output = native.read_text(errors='replace')
    history, linear = context(run_dir, task, entry)
    keyword, assigned_native, framework_change = _vpt2_framework_recovery(
        task, output, history=history, linear=linear)
    history.append(symmetry('#p ' + keyword)[0])
    groups = list(dict.fromkeys(s.split('=', 1)[1].split(',')[0]
                               for s in history if s.startswith('PG=')))
    replacement = deepcopy(task)
    replacement['input_template'] = replace_symmetry(task['input_template'], keyword)
    replacement['recovery'] = {
        'kind': 'gaussian_vpt2_framework_group',
        'change': keyword,
        'assigned_native_framework': assigned_native,
        'attempted_framework_groups': groups,
        'attempted_symmetry_settings': list(dict.fromkeys(history)),
        'observed_framework_change': framework_change,
        'symmetry_detection': 'loose' if ',Loose' in keyword else 'default',
        'orientation': ('center_of_mass' if 'CenterOfMass' in keyword
                        else 'input' if keyword == 'NoSymm' else 'default'),
        'geometry_changed': False,
    }
    if 'manual_fallback' in replacement:
        # Keep an operator fallback's original binding and all its other edits.
        # Validate its receipt before extending it with this automatic recovery.
        from kinbot.manual_fallback import logical_task, digest
        logical_task(task)
        replacement['manual_fallback']['patched_sha256'] = digest(replacement)
    spec['tasks'][spec['tasks'].index(task)] = replacement
    assign_partitions(spec)
    validate_spec(spec)

    # Journal before archival and fence concurrent allocated workers with the
    # manifest lock. A process interruption must not split workflow/state.
    reset = {ident} | {n for n in children if state['tasks'].get(n, {}).get('status')
                      in ('blocked', 'skipped')}
    reason = 'Gaussian VPT2 framework-group recovery with ' + keyword
    live_control.commit(run_dir, spec, state, reset, {n: reason for n in reset})
    _, spec, state = _load(run_dir)
    archive = Path(state['pending_retries'][ident]['previous_attempt'])
    _stage_task(run_dir, spec, state, replacement)
    state['tasks'][ident]['migration'] = \
        'Gaussian VPT2 framework-group recovery with ' + keyword
    _clear_blocked_tasks(state)
    _atomic_json(run_dir / 'state.json', state)
    return archive


def recover_vpt2_failures(run_dir, spec, state):
    """Apply known framework fallbacks for every Gaussian VPT2 dispatcher task.

    Selection is by parser/backend and native diagnostic, never task name,
    molecule, ANL variant or the controller that happens to own the graph.
    """
    from kinbot.anl.timeout_recovery import _job_id
    for task in list(spec['tasks']):
        ident = task['id']
        entry = state['tasks'].get(ident, {})
        if not can_recover_gaussian_vpt2_symmetry(run_dir, task, entry):
            continue
        job = _job_id(Path(run_dir)/'tasks'/ident, ident, entry)
        if job and _job_active(job):
            continue
        # Commit the latest allocation reconciliation before starting the
        # journaled migration, whose public API reloads its own manifest.
        _atomic_json(Path(run_dir)/'state.json', state)
        try:
            recover_gaussian_vpt2_symmetry(run_dir, ident)
        except VPT2RecoveryUnavailable:
            continue  # Unsupported user keywords remain available for manual recovery.
        except RuntimeError as error:
            if 'legacy workers' not in str(error):
                raise
            # A legacy allocation cannot take the manifest read lock. Keep it
            # running and reconsider this recovery at the next dispatcher poll.
            continue
        _, spec, state = _load(run_dir)
    return spec, state


def reroute_cfour_higher_order_to_mrcc(
        run_dir, ident, *, command='dmrcc', walltime=None,
        partition=None, max_cores=None):
    """Replace one failed CFOUR CCSDT(Q) probe with direct MRCC.

    CFOUR 2.1 selects its restricted closed-shell ``xncc`` solver for
    CCSDT(Q), even when a generated ZMAT requests ``CC_PROG=VCC``.  Such a
    result cannot supply the required RHF determinant plus unrestricted CC
    component.  This narrow migration archives the rejected attempt, rewrites
    only that task in the prepared workflow, and stages an RHF-UCCSDT(Q)
    direct-MRCC replacement on the identical geometry.
    """
    run_dir, spec, state = _load(run_dir)
    by_id = {task['id']: task for task in spec['tasks']}
    task = by_id.get(ident)
    entry = state['tasks'].get(ident)
    if task is None or entry is None:
        raise ValueError(f'{ident}: task is unavailable.')
    parser = task.get('result_parser', {})
    if (task.get('kind') != 'external'
            or task.get('backend', '').lower() != 'cfour'
            or parser.get('kind') != 'cfour_energy'
            or parser.get('method') != 'CCSDT(Q)'
            or parser.get('reference') != 'RHF'):
        raise ValueError(
            f'{ident}: rerouting requires a closed-shell CFOUR CCSDT(Q) task.')
    if entry.get('status') not in ('failed', 'unavailable'):
        raise ValueError(
            f'{ident}: only a failed or unavailable CFOUR task can be migrated.')
    if entry.get('job_id') and _job_active(entry['job_id']):
        raise RuntimeError(f'{ident}: Slurm job is still active.')
    _verify_stage_files(run_dir, task, entry)
    resources = deepcopy(task['resources'])
    selected_walltime = walltime or resources['walltime']
    selected_partition = (partition if partition is not None
                          else resources.get('partition'))
    selected_max_cores = (max_cores if max_cores is not None
                          else resources.get('max_cores',
                                             resources.get('cores')))
    if selected_max_cores == 'auto':
        selected_max_cores = 8
    _positive_integer(selected_max_cores, 'MRCC reroute max_cores')
    replacement = mrcc_task(
        ident, 'CCSDT(Q)', parser['basis'], multiplicity=1, reference='RHF',
        geometry_from=task.get('geometry_from', 'initial'),
        walltime=selected_walltime, max_cores=selected_max_cores,
        partition=selected_partition, command=command)
    replacement['resources'] = resources
    replacement['resources']['walltime'] = selected_walltime
    replacement['resources']['max_cores'] = selected_max_cores
    replacement['resources']['min_memory_mb_per_core'] = 4096
    if selected_partition is None:
        replacement['resources'].pop('partition', None)
    else:
        replacement['resources']['partition'] = selected_partition
    if 'depends_on' in task:
        replacement['depends_on'] = deepcopy(task['depends_on'])
    spec['tasks'][spec['tasks'].index(task)] = replacement
    validate_spec(spec)

    old = run_dir / 'tasks' / ident
    attempt = entry.get('attempt', 1)
    archive = run_dir / 'attempts' / ident / str(attempt)
    if archive.exists():
        raise RuntimeError(f'{ident}: retry archive already exists.')
    archive.parent.mkdir(parents=True, exist_ok=True)
    old_workflow = (run_dir / 'workflow.json').read_text()
    old.replace(archive)
    try:
        _atomic_json(run_dir / 'workflow.json', spec)
        state['spec_sha256'] = _file_hash(run_dir / 'workflow.json')
        _stage_task(run_dir, spec, state, replacement)
    except Exception:
        failed_stage = run_dir / 'tasks' / ident
        if failed_stage.exists():
            shutil.rmtree(failed_stage)
        archive.replace(old)
        (run_dir / 'workflow.json').write_text(old_workflow)
        state['spec_sha256'] = _file_hash(run_dir / 'workflow.json')
        state['tasks'][ident] = entry
        _atomic_json(run_dir / 'state.json', state)
        raise
    state['tasks'][ident]['attempt'] = attempt + 1
    state['tasks'][ident]['previous_attempt'] = str(archive)
    state['tasks'][ident]['migration'] = \
        'rejected CFOUR xncc CCSDT(Q) to direct MRCC RHF-UCCSDT(Q)'
    _clear_blocked_tasks(state)
    _atomic_json(run_dir / 'state.json', state)
    return archive


def migrate_cfour_vcc_keyword(run_dir, ident):
    """Backward-compatible alias for the corrected CFOUR-to-MRCC migration."""
    return reroute_cfour_higher_order_to_mrcc(run_dir, ident)


def resume_mrcc_failed(run_dir, ident, *, walltime='partition_max',
                       partition=None, max_cores=None):
    """Restage a time-limited direct-MRCC task without discarding amplitudes.

    MRCC documents ``rest=1`` for restarting canonical CC calculations from
    saved amplitudes. This path keeps the task's working memory and scratch
    files, archives its prior input and text outputs, and regenerates the
    immutable dispatcher records with a longer resource request.
    """
    run_dir, spec, state = _load(run_dir)
    by_id = {task['id']: task for task in spec['tasks']}
    task = by_id.get(ident)
    entry = state['tasks'].get(ident)
    if task is None or entry is None:
        raise ValueError(f'{ident}: task is unavailable.')
    if (task.get('kind') != 'external'
            or task.get('backend', '').lower() != 'mrcc'):
        raise ValueError(f'{ident}: resume-mrcc requires a direct MRCC task.')
    if entry.get('status') not in ('failed', 'unavailable'):
        raise ValueError(
            f'{ident}: only a failed or unavailable MRCC task can be resumed.')
    if entry.get('job_id') and _job_active(entry['job_id']):
        raise RuntimeError(f'{ident}: Slurm job is still active.')
    if any(other != ident and other in state['tasks']
           and state['tasks'][other].get('status') not in ('blocked', 'skipped')
           and ident in _task_dependencies(other_task)
           for other, other_task in by_id.items()):
        raise RuntimeError(f'{ident}: downstream task already exists; review run.')
    directory = run_dir / 'tasks' / ident
    _verify_stage_files(run_dir, task, entry)
    amplitude = directory / 'fort.16'
    if not amplitude.is_file() or amplitude.stat().st_size == 0:
        raise RuntimeError(f'{ident}: nonempty MRCC fort.16 is required to resume.')
    if (not isinstance(walltime, str)
            or (walltime != 'partition_max' and not _TIME.fullmatch(walltime))):
        raise ValueError('MRCC resume walltime must be partition_max, HH:MM:SS or D-HH:MM:SS.')
    if partition is not None:
        _basename(partition, 'MRCC resume partition')
    if max_cores is not None:
        _positive_integer(max_cores, 'MRCC resume max_cores')

    replacement = deepcopy(task)
    resources = replacement['resources']
    resources['walltime'] = walltime
    resources['cores'] = 'auto'
    if max_cores is not None:
        resources['max_cores'] = max_cores
    if partition is not None:
        resources['partition'] = partition
    template = replacement['input_template']
    if re.search(r'^\s*rest\s*=\s*1\s*$', template,
                 re.IGNORECASE | re.MULTILINE) is None:
        template, count = re.subn(
            r'(^\s*ccprog\s*=\s*mrcc\s*$)', r'\1\nrest=1', template,
            count=1, flags=re.IGNORECASE | re.MULTILINE)
        if count != 1:
            raise RuntimeError(f'{ident}: MRCC input has no unique ccprog=mrcc.')
        replacement['input_template'] = template
    replacement['restart'] = {
        'kind': 'mrcc_amplitudes', 'keyword': 'rest=1',
        'checkpoint': 'fort.16', 'sha256': _file_hash(amplitude),
    }
    if (directory/'fort.17').is_file():
        replacement['restart']['fort17_sha256'] = _file_hash(directory/'fort.17')
    spec['tasks'][spec['tasks'].index(task)] = replacement
    assign_partitions(spec)
    validate_spec(spec)

    attempt = entry.get('attempt', 1)
    archive = run_dir / 'attempts' / ident / str(attempt)
    if archive.exists():
        raise RuntimeError(f'{ident}: retry archive already exists.')
    archive.mkdir(parents=True)
    archive_names = {
        'task.json', 'job.slurm', task['input_name'], 'execution.json',
        'slurm.stdout', 'slurm.stderr', task.get('stdout', 'stdout.txt'),
        task.get('stderr', 'stderr.txt'), 'EXIT',
    }
    archive_names.update(task.get('required_outputs', []))
    for name in sorted(archive_names):
        source = directory / name
        if source.is_file():
            shutil.copy2(source, archive / name)
    _atomic_json(archive / 'resume_manifest.json', {
        'schema': 1, 'task_id': ident, 'attempt': attempt,
        'checkpoint': 'fort.16',
        'checkpoint_bytes': amplitude.stat().st_size,
        'working_memory_mb': task['resources']['memory_mb'],
    })

    atoms = read(directory / 'geometry.xyz')
    rendered = _render_input(
        replacement['input_template'], atoms, spec['molecule'],
        replacement['resources'])
    (directory / replacement['input_name']).write_text(rendered)
    initial = _atoms(spec['molecule'])
    _atomic_json(directory / 'task.json', {
        'schema': 1, 'task': replacement,
        'molecule': {'symbols': initial.get_chemical_symbols(),
                     'charge': spec['molecule'].get('charge', 0),
                     'multiplicity': spec['molecule'].get('multiplicity', 1)},
        'geometry_sha256': entry['geometry_sha256'],
    })
    (directory / 'job.slurm').write_text(
        _slurm_script(replacement, directory, state['python']))
    for name in {'execution.json', 'slurm.stdout', 'slurm.stderr',
                 replacement.get('stdout', 'stdout.txt'),
                 replacement.get('stderr', 'stderr.txt'), 'EXIT'}:
        path = directory / name
        if path.is_file():
            path.unlink()

    _atomic_json(run_dir / 'workflow.json', spec)
    state['spec_sha256'] = _file_hash(run_dir / 'workflow.json')
    state['tasks'][ident] = {
        'status': 'staged',
        'geometry_from': replacement.get('geometry_from', 'initial'),
        'geometry_sha256': entry['geometry_sha256'],
        'task_sha256': _file_hash(directory / 'task.json'),
        'job_sha256': _file_hash(directory / 'job.slurm'),
        'input_sha256': _file_hash(directory / replacement['input_name']),
        'attempt': attempt + 1,
        'previous_attempt': str(archive),
        'restart_checkpoint': 'fort.16',
        'restart_checkpoint_bytes': amplitude.stat().st_size,
    }
    _clear_blocked_tasks(state)
    _atomic_json(run_dir / 'state.json', state)
    return archive


def advance(run_dir, submit=False, submit_only=None):
    """Complete finished tasks, stage ready tasks, and optionally submit jobs."""
    from kinbot.anl import account_recovery, scratch_recovery, timeout_recovery, lock_recovery, optimization_recovery
    from kinbot import timeout_retry
    from kinbot.anl import live_control, interface_recovery
    live_control.resume(run_dir)
    account_recovery.resume(run_dir)
    scratch_recovery.resume(run_dir)
    timeout_recovery.resume(run_dir)
    lock_recovery.resume(run_dir)
    optimization_recovery.resume(run_dir)
    run_dir, spec, state = _load(run_dir)
    from kinbot.reaction_decisions import graph_rejected
    if graph_rejected(run_dir):
        return state
    from kinbot.anl.refinement import pending as pending_refinement
    refining = pending_refinement(run_dir)
    by_id = {task['id']: task for task in spec['tasks']}
    if submit_only is not None:
        submit_only = set(submit_only)
        unknown = submit_only - by_id.keys()
        if unknown:
            raise ValueError(f"Unknown submission task(s): {', '.join(sorted(unknown))}")

    # Reconcile two deterministic non-QC cases before dependency handling.
    # First, a newer parser may accept a native calculation that an older
    # KinBot rejected, or a site launcher may have reported nonzero after the
    # native program recorded normal termination.  Second, an atom has no
    # geometry coordinate to optimize, so its geometry nodes are exact
    # identity operations and do not belong in the scheduler at all.
    for ident, entry in state['tasks'].items():
        task = by_id[ident]
        outcome = run_dir / 'tasks' / ident / 'execution.json'
        previous = (json.loads(outcome.read_text())
                    if outcome.is_file() else None)
        if (entry.get('status') == 'failed' and previous is not None):
            try:
                previous = _recover_one_electron_f12_nan(
                    run_dir, spec, state, task, entry, previous)
            except ValueError:
                pass
        if (entry.get('status') in ('failed', 'unavailable') and previous is not None
                and _native_success_recovery_candidate(task, previous)
                and _native_success_present(
                    run_dir / 'tasks' / ident, task)):
            try:
                recovered, actual_hash = _reparse_execution(
                    run_dir, task, entry, previous)
            except ValueError:
                # The current parser still rejects the native result.  Keep
                # the original failure available for diagnosis or migration.
                pass
            else:
                entry.update(status='complete', execution='executed',
                             automatically_reparsed=True)
                entry.pop('error', None)
                entry.pop('job_id', None)
                entry.pop('unavailable_reason', None)
                entry.pop('failure_phase', None)
                if actual_hash is not None:
                    entry['final_geometry_sha256'] = actual_hash
                previous = recovered
        if (_monatomic_geometry_task(spec, task)
                and entry.get('status') in ('staged', 'failed')):
            if entry['status'] == 'staged' and previous is not None:
                # Let the ordinary reconciliation below accept an already
                # written execution record after a state-file interruption.
                continue
            _complete_monatomic_geometry(
                run_dir, spec, task, entry,
                previous=previous if entry['status'] == 'failed' else None)

    for ident, entry in state['tasks'].items():
        task = by_id[ident]
        if entry['status'] == 'submitting':
            raise RuntimeError(f'{ident}: submission outcome uncertain; reconcile Slurm job manually.')
        if entry['status'] in ('staged', 'submitted', 'complete'):
            _verify_stage_files(run_dir, task, entry)
        if entry['status'] == 'complete':
            outcome = run_dir / 'tasks' / ident / 'execution.json'
            if not outcome.is_file():
                raise RuntimeError(f'{ident}: completed execution record is missing.')
            actual_hash = _verify_execution(run_dir, task, entry,
                                            json.loads(outcome.read_text()))
            if actual_hash != entry.get('final_geometry_sha256'):
                raise RuntimeError(f'{ident}: completed geometry state changed.')
        if entry['status'] not in ('staged', 'submitted'):
            continue
        # The retry guard may replace execution.json while this allocation
        # is alive. Query first, then read the final record, never a failed
        # intermediate attempt captured before a slow scheduler response.
        if entry['status'] == 'submitted' and _job_active(entry['job_id']):
            timeout_retry.observe(run_dir / 'tasks' / ident, entry['job_id'])
            entry.pop('missing_since', None)
            continue
        outcome = run_dir / 'tasks' / ident / 'execution.json'
        if outcome.exists():
            execution = json.loads(outcome.read_text())
            entry.pop('missing_since', None)
            if (_native_success_recovery_candidate(task, execution)
                    and _native_success_present(
                        run_dir / 'tasks' / ident, task)):
                try:
                    execution, _ = _reparse_execution(
                        run_dir, task, entry, execution)
                    entry['automatically_reparsed'] = True
                except ValueError:
                    # Preserve a real parser rejection as a failed task.
                    pass
            actual_hash = _verify_execution(run_dir, task, entry, execution)
            entry['status'] = ('complete' if execution['status'] == 'executed'
                               else 'failed')
            entry['execution'] = execution['status']
            if actual_hash is not None:
                entry['final_geometry_sha256'] = actual_hash
        elif entry['status'] == 'submitted':
            # Queue exit can precede visibility of the atomic result rename
            # on the controller's shared-filesystem client. Retain ownership
            # and the node-cap slot throughout this bounded publication wait.
            missing_since = entry.setdefault('missing_since', time.time())
            if time.time() - missing_since <= _RESULT_PUBLICATION_GRACE:
                continue
            entry['status'] = 'failed'
            entry['error'] = 'Slurm job disappeared without execution.json.'
            entry['failure_phase'] = 'scheduler'

    # Old pinned workers may fail after replacement controllers start. Repair
    # only this recognized startup bug, once, before deriving component skips.
    if not spec.get('intent', {}).get('vrc_fragment_retry') and not refining:
        state = lock_recovery.recover(run_dir, spec, state)
        state = interface_recovery.recover(run_dir, spec, state)
        from kinbot.anl import molpro_mpi_recovery, molpro_solver_recovery, molpro_gradient_recovery
        spec, state = molpro_mpi_recovery.recover(run_dir, spec, state)
        spec, state = molpro_solver_recovery.recover(run_dir, spec, state)
        from kinbot.anl import cc_iteration_recovery
        spec, state = cc_iteration_recovery.recover(run_dir, spec, state)
        spec, state = recover_vpt2_failures(run_dir, spec, state)
        from kinbot.anl import vpt2_refinement
        spec, state = vpt2_refinement.recover(run_dir, spec, state)
        spec, state = molpro_gradient_recovery.recover(run_dir, spec, state)
        state = optimization_recovery.recover(run_dir, spec, state)
        state = account_recovery.recover(run_dir, spec, state)
        state = scratch_recovery.recover(run_dir, spec, state)
        spec, state = timeout_recovery.recover(run_dir, spec, state)
    by_id = {task['id']: task for task in spec['tasks']}

    # An optional component that cannot execute is unavailable at this
    # requested level. Keep the failed native record for provenance, stop
    # only explicitly related work that has not been submitted, and allow
    # independent components to finish so the CBH ladder can use a lower
    # complete composite. Parser/provenance failures remain hard failures.
    skip_queue = []
    for ident, entry in state['tasks'].items():
        task = by_id[ident]
        if (entry.get('status') != 'failed'
                or task.get('failure_policy') != 'skip_component'):
            continue
        from kinbot.live_control import holding
        if holding(run_dir/'tasks'/ident):
            continue  # Keep the retry window open before optional fallback.
        outcome = run_dir / 'tasks' / ident / 'execution.json'
        execution = (json.loads(outcome.read_text())
                     if outcome.is_file() else {})
        phase = entry.get('failure_phase') or execution.get('failure_phase')
        if phase not in ('execution', 'scheduler'):
            continue
        entry.update(
            status='unavailable',
            unavailable_reason=(execution.get('error')
                                or entry.get('error')
                                or 'native calculation unavailable'),
            component_families=list(task.get('component_families', ())),
        )
        entry.pop('job_id', None)
        skip_queue.extend((target, ident)
                          for target in task.get('failure_skips', ())
                          if target in by_id)

    skipped_from = {}
    while skip_queue:
        target, source = skip_queue.pop(0)
        target_task = by_id[target]
        entry = state['tasks'].get(target)
        if entry is not None and entry.get('status') in (
                'submitted', 'complete', 'unavailable', 'failed'):
            continue
        sources = skipped_from.setdefault(target, set())
        if source in sources:
            continue
        sources.add(source)
        replacement = dict(entry or {})
        replacement.update({
            'status': 'skipped',
            'geometry_from': target_task.get('geometry_from', 'initial'),
            'skipped_by': sorted(sources),
            'skip_reason': 'requested component is unavailable',
            'component_families': list(
                target_task.get('component_families', ())),
        })
        if entry is not None and entry.get('status') == 'staged':
            replacement['resume_status'] = 'staged'
        state['tasks'][target] = replacement
        skip_queue.extend((child, target)
                          for child in target_task.get('failure_skips', ())
                          if child in by_id)
    for task in spec['tasks']:
        ident = task['id']
        if ident in refining:
            continue  # Do not admit new descendants of a retiring geometry.
        dependencies = _task_dependencies(task, state)
        failed_dependencies = [
            dependency for dependency in dependencies
            if state['tasks'].get(dependency, {}).get('status') in
            ('failed', 'blocked', 'unavailable', 'skipped')]
        entry = state['tasks'].get(ident)
        from kinbot.live_control import holding
        sources = set(_task_dependencies(task)) | set(task.get('fallback_on', ()))
        held_sources = [n for n in sources if state['tasks'].get(n,{}).get('status')
                        in ('failed','unavailable') and not state['tasks'][n].get('vpt2_refinement_exhausted')
                        and holding(run_dir/'tasks'/n)]
        if held_sources and (entry is None or entry.get('status')=='blocked'):
            state['tasks'][ident] = dict(status='blocked',blocked_by=held_sources,
                                        reason='Upstream failure is held for operator review.')
            continue
        if (entry is not None and entry.get('status') == 'skipped'
                and entry.get('skipped_by')
                and all(state['tasks'].get(source, {}).get('status') == 'complete'
                        for source in entry['skipped_by'])):
            # A verified native reparse can recover an optional source that
            # was previously unavailable. Release only its derived skips;
            # unrelated unavailability and deliberate fallback skips remain.
            if entry.get('resume_status') == 'staged':
                entry['status'] = 'staged'
                for key in ('resume_status', 'skipped_by', 'skip_reason',
                            'component_families'):
                    entry.pop(key, None)
            else:
                state['tasks'].pop(ident)
                entry = None
        if entry is not None and entry.get('status') == 'blocked':
            if failed_dependencies:
                entry['blocked_by'] = failed_dependencies
                continue
            # A failed dependency may have been explicitly retried. Return the
            # child to the ordinary dependency gate without fabricating a
            # staged calculation.
            state['tasks'].pop(ident)
            entry = None
        if entry is not None:
            continue
        fallback_on = task.get('fallback_on', ())
        if fallback_on:
            fallback_states = {
                dependency: state['tasks'].get(dependency, {}).get('status')
                for dependency in fallback_on}
            if any(value == 'complete' for value in fallback_states.values()):
                state['tasks'][ident] = {
                    'status': 'skipped',
                    'geometry_from': task.get('geometry_from', 'initial'),
                    'skip_reason': 'higher tier completed; fallback not needed',
                    'fallback_on': list(fallback_on),
                }
                continue
            if not all(value in ('unavailable', 'skipped')
                       for value in fallback_states.values()):
                # A pending trigger must finish before fallback is decided.
                # A hard failed/blocked trigger remains visible to the caller.
                continue
        if failed_dependencies:
            soft = all(state['tasks'][dependency].get('status') in
                       ('unavailable', 'skipped')
                       for dependency in failed_dependencies)
            state['tasks'][ident] = {
                'status': 'skipped' if soft else 'blocked',
                'geometry_from': task.get('geometry_from', 'initial'),
                ('skipped_by' if soft else 'blocked_by'): failed_dependencies,
                **({'skip_reason': 'dependency is unavailable'} if soft else {}),
            }
        elif all(state['tasks'].get(dep, {}).get('status') == 'complete'
                 for dep in dependencies):
            _stage_task(run_dir, spec, state, task)
    _atomic_json(run_dir / 'state.json', state)
    if submit:
        from kinbot import run_policy
        active = sum(entry['status'] == 'submitted' for entry in state['tasks'].values())
        for task in spec['tasks']:
            ident = task['id']
            if ident in refining:
                continue
            entry = state['tasks'].get(ident)
            if active >= run_policy.limit(spec['limits']['max_nodes'], run_dir):
                break
            if submit_only is not None and ident not in submit_only:
                continue
            if not entry or entry['status'] != 'staged':
                continue
            # Migrate an older, verified staged script when it next submits;
            # completed/running allocations and native input hashes stay put.
            from kinbot.slurm_retry import wrap_script
            directory = run_dir / 'tasks' / ident
            script_path = directory / 'job.slurm'
            original_script = script_path.read_text()
            guarded = wrap_script(original_script, identity=ident,
                                  directories=[str(directory)], python=state['python'])
            guarded = scratch_recovery.apply_exclusions(guarded, state)
            if guarded != original_script:
                script_path.write_text(guarded)
                entry['job_sha256'] = _file_hash(script_path)
            entry['status'] = 'submitting'
            _atomic_json(run_dir / 'state.json', state)
            directory = run_dir / 'tasks' / ident
            submission_options = {}
            if spec.get('intent', {}).get('vrc_fragment_retry'):
                from kinbot.run_walltime import seconds, duration
                manual_time = task.get('manual_fallback', {}).get('options', {}).get('resources', {}).get('walltime')
                submission_options['walltime_ceiling'] = (manual_time or duration(min(
                    seconds(task['resources']['walltime']), 7200)))
            if state.get('timeout_retries', {}).get(ident):
                # SBATCH_* inherited from the controller must not shorten a
                # verified replacement budget or move it to another partition.
                submission_options['env'] = dict(
                    os.environ, SBATCH_TIMELIMIT=task['resources']['walltime'],
                    SBATCH_PARTITION=task['resources']['partition'])
            try:
                response = run_policy.submit(
                    ['sbatch', '--parsable', 'job.slurm'], cwd=directory,
                    capture_output=True, text=True, check=False,
                    timeout=_SCHEDULER_TIMEOUT, **submission_options)
            except SchedulerUnavailableError:
                # The shared pool could not count occupancy; no sbatch was
                # attempted. Preserve a runnable task for the next poll.
                entry['status'] = 'staged'
                _atomic_json(run_dir / 'state.json', state)
                raise
            except subprocess.TimeoutExpired as error:
                raise RuntimeError(
                    f'{ident}: sbatch timed out; submission outcome is '
                    'uncertain and state remains submitting for manual '
                    'reconciliation.') from error
            except OSError as error:
                if not getattr(error, 'kinbot_submission_not_started', False):
                    raise RuntimeError(
                        f'{ident}: submission outcome uncertain after I/O error; '
                        'state remains submitting for manual reconciliation.') from error
                entry['status'] = 'staged'
                _atomic_json(run_dir / 'state.json', state)
                if isinstance(error, FileNotFoundError):
                    raise
                raise SchedulerUnavailableError(
                    f'{ident}: sbatch could not start: {error}') from error
            if response is None:
                entry['status'] = 'staged'
                _atomic_json(run_dir / 'state.json', state)
                break
            if response.returncode:
                entry['status'] = 'staged'
                _atomic_json(run_dir / 'state.json', state)
                detail = response.stderr.strip() or response.stdout.strip()
                if _scheduler_unavailable(detail):
                    raise SchedulerUnavailableError(
                        f'{ident}: sbatch is unavailable: {detail}')
                raise RuntimeError(f'{ident}: sbatch failed: {response.stderr.strip()}')
            job_id = response.stdout.strip().split(';', 1)[0]
            if not job_id.isdecimal():
                raise RuntimeError(f'{ident}: unrecognized sbatch job id {response.stdout!r}; '
                                   'reconcile submission manually.')
            entry.update(status='submitted', job_id=job_id)
            _atomic_json(run_dir / 'state.json', state)
            active += 1
    return state


def _summary(spec, state):
    return {task['id']: state['tasks'].get(task['id'], {'status': 'waiting'})['status']
            for task in spec['tasks']}


def refresh_status(run_dir):
    """Reconcile finished jobs before reporting, without submitting QC work."""
    run_dir = Path(run_dir).resolve()
    if not (run_dir / 'workflow.json').is_file():
        return {'workflow': 'not_prepared'}
    with (run_dir / 'drive.lock').open('w') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            # A running driver owns reconciliation; report its last atomic
            # snapshot rather than racing it or waiting for its polling loop.
            _, spec, state = _load(run_dir)
        else:
            state = advance(run_dir, submit=False)
            _, spec, _ = _load(run_dir)
    return _summary(spec, state)


def drive_run(run_dir, *, interval=20., submit_only=None, hold_failures=True):
    """Drive a dispatcher graph without creating an orphanable controller."""
    run_dir = Path(run_dir).resolve()
    with (run_dir / 'drive.lock').open('w') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(
                'Another driver is already managing this run.') from exc
        scheduler_error = None
        while True:
            # A live PES review must yield graph ownership even when every
            # task is already submitted (so _submit is never reached). The
            # epoch exception releases drive.lock without touching Slurm jobs
            # or their ledger; the resumed branch reattaches to those IDs.
            from kinbot.pes_phase import checkpoint
            checkpoint()
            try:
                from kinbot.anl.live_control import process as process_controls
                process_controls(run_dir)
                state = advance(
                    run_dir, submit=True, submit_only=submit_only)
            except SchedulerUnavailableError as error:
                # A status-query failure says nothing about the submitted
                # job.  Preserve the atomic state and wait for Slurm rather
                # than failing the task or the scientific branch that owns
                # this in-process driver.
                message = str(error)
                if message != scheduler_error:
                    print(
                        'Scheduler status is unavailable; retaining submitted '
                        f'tasks and retrying: {message}',
                        file=sys.stderr, flush=True)
                    scheduler_error = message
                time.sleep(float(interval))
                continue
            if scheduler_error is not None:
                print(
                    'Scheduler status query recovered; resuming dispatcher.',
                    file=sys.stderr, flush=True)
                scheduler_error = None
            _, spec, _ = _load(run_dir)
            summary = _summary(spec, state)
            terminal = {
                'complete', 'failed', 'blocked', 'unavailable', 'skipped'}
            if all(value in terminal for value in summary.values()):
                from kinbot.live_control import holding
                if hold_failures and any(value in ('failed','unavailable')
                        and not state['tasks'].get(name, {}).get('vpt2_refinement_exhausted') and holding(
                        run_dir/'tasks'/name) for name,value in summary.items()):
                    time.sleep(float(interval))
                    continue
                successful = all(value in (
                    'complete', 'unavailable', 'skipped')
                    for value in summary.values())
                return (0 if successful else 1), summary
            time.sleep(float(interval))


def main(argv=None):
    parser = argparse.ArgumentParser(description='Stage and drive exclusive QC jobs')
    sub = parser.add_subparsers(dest='action', required=True)
    stage = sub.add_parser('prepare')
    stage.add_argument('spec')
    stage.add_argument('run_dir')
    drive = sub.add_parser('drive')
    drive.add_argument('run_dir')
    drive.add_argument('--once', action='store_true')
    drive.add_argument('--only', action='append', metavar='TASK_ID')
    drive.add_argument('--interval', type=int, default=20)
    status = sub.add_parser('status')
    status.add_argument('run_dir')
    status.add_argument('--cached', action='store_true',
                        help='show saved state without querying Slurm')
    check = sub.add_parser('preflight')
    check.add_argument('run_dir')
    retry = sub.add_parser('retry')
    retry.add_argument('run_dir')
    retry.add_argument('task_id')
    recover_vpt2 = sub.add_parser('recover-gaussian-vpt2-symmetry')
    recover_vpt2.add_argument('run_dir')
    recover_vpt2.add_argument('task_id')
    migrate_cfour = sub.add_parser('migrate-cfour-vcc')
    migrate_cfour.add_argument('run_dir')
    migrate_cfour.add_argument('task_id')
    reroute_cfour = sub.add_parser('reroute-cfour-mrcc')
    reroute_cfour.add_argument('run_dir')
    reroute_cfour.add_argument('task_id')
    reroute_cfour.add_argument('--mrcc-command', default='dmrcc')
    reroute_cfour.add_argument('--walltime')
    reroute_cfour.add_argument('--partition')
    reroute_cfour.add_argument('--max-cores', type=int)
    resume_mrcc = sub.add_parser('resume-mrcc')
    resume_mrcc.add_argument('run_dir')
    resume_mrcc.add_argument('task_id')
    resume_mrcc.add_argument('--walltime', default='partition_max')
    resume_mrcc.add_argument('--partition')
    resume_mrcc.add_argument('--max-cores', type=int)
    reparse = sub.add_parser('reparse')
    reparse.add_argument('run_dir')
    reparse.add_argument('task_id')
    run = sub.add_parser('run-task')
    run.add_argument('task_file')
    args = parser.parse_args(argv)
    if args.action == 'prepare':
        print(prepare(args.spec, args.run_dir))
    elif args.action == 'run-task':
        run_task(args.task_file)
    elif args.action == 'status':
        if args.cached and not (Path(args.run_dir) / 'workflow.json').is_file():
            summary = {'workflow': 'not_prepared'}
        elif args.cached:
            _, spec, state = _load(args.run_dir)
            summary = _summary(spec, state)
        else:
            summary = refresh_status(args.run_dir)
        print(json.dumps(summary, indent=2))
    elif args.action == 'preflight':
        print(json.dumps(preflight(args.run_dir), indent=2))
    elif args.action == 'retry':
        run_dir = Path(args.run_dir).resolve()
        with (run_dir / 'drive.lock').open('w') as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise RuntimeError('Another driver is managing this run.') from exc
            print(retry_failed(run_dir, args.task_id))
    elif args.action == 'recover-gaussian-vpt2-symmetry':
        run_dir = Path(args.run_dir).resolve()
        with (run_dir / 'drive.lock').open('w') as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise RuntimeError('Another driver is managing this run.') from exc
            print(recover_gaussian_vpt2_symmetry(run_dir, args.task_id))
    elif args.action == 'migrate-cfour-vcc':
        run_dir = Path(args.run_dir).resolve()
        with (run_dir / 'drive.lock').open('w') as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise RuntimeError('Another driver is managing this run.') from exc
            print(migrate_cfour_vcc_keyword(run_dir, args.task_id))
    elif args.action == 'reroute-cfour-mrcc':
        run_dir = Path(args.run_dir).resolve()
        with (run_dir / 'drive.lock').open('w') as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise RuntimeError('Another driver is managing this run.') from exc
            print(reroute_cfour_higher_order_to_mrcc(
                run_dir, args.task_id, command=args.mrcc_command,
                walltime=args.walltime, partition=args.partition,
                max_cores=args.max_cores))
    elif args.action == 'resume-mrcc':
        run_dir = Path(args.run_dir).resolve()
        with (run_dir / 'drive.lock').open('w') as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise RuntimeError('Another driver is managing this run.') from exc
            print(resume_mrcc_failed(
                run_dir, args.task_id, walltime=args.walltime,
                partition=args.partition, max_cores=args.max_cores))
    elif args.action == 'reparse':
        run_dir = Path(args.run_dir).resolve()
        with (run_dir / 'drive.lock').open('w') as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise RuntimeError('Another driver is managing this run.') from exc
            print(json.dumps(reparse_failed(run_dir, args.task_id), indent=2))
    else:
        if args.only and not args.once:
            parser.error('--only requires --once')
        if args.interval < 1 or args.interval > 60:
            parser.error('--interval must be 1 through 60 seconds')
        if args.once:
            run_dir = Path(args.run_dir).resolve()
            with (run_dir / 'drive.lock').open('w') as lock:
                try:
                    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError as exc:
                    raise RuntimeError(
                        'Another driver is already managing this run.') from exc
                state = advance(run_dir, submit=True, submit_only=args.only)
                _, spec, _ = _load(run_dir)
                summary = _summary(spec, state)
                print(json.dumps(summary, sort_keys=True), flush=True)
                return int(any(value in ('failed', 'blocked')
                               for value in summary.values()))
        returncode, summary = drive_run(
            args.run_dir, interval=args.interval, submit_only=args.only)
        print(json.dumps(summary, sort_keys=True), flush=True)
        return returncode


if __name__ == '__main__':
    sys.exit(main())
