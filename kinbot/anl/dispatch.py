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
from kinbot.theory import TheoryProfile


_SAFE_NAME = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*\Z')
_TOKEN = re.compile(r'\{\{([A-Z_]+)\}\}')
_TIME = re.compile(r'(?:\d+-)?\d{1,3}:[0-5]\d:[0-5]\d\Z')
_INPUT_TOKENS = {'CARTESIAN', 'XYZ', 'MRCC_XYZ', 'NATOMS', 'CHARGE',
                 'MULT', 'NELECTRONS', 'SPIN', 'CORES', 'MEMORY_MB',
                 'WORK_MEMORY_MB'}
_RESERVED_FILES = {'geometry.xyz', 'task.json', 'job.slurm', 'execution.json',
                   'slurm.stdout', 'slurm.stderr', 'drive.lock'}


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


def _file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


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
    dependencies = task.get('depends_on', [])
    if not isinstance(dependencies, list) or any(dep not in ids for dep in dependencies):
        raise ValueError(f'{ident}: depends_on must contain task ids.')
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
                set(options) - {'fmax', 'steps', 'sella_kwargs'}):
            raise ValueError(f'{ident}: invalid optimizer options.')
        if ('sella_kwargs' in options and not isinstance(options['sella_kwargs'], dict)):
            raise ValueError(f'{ident}: sella_kwargs must be an object.')
        if ('fmax' in options and
                (isinstance(options['fmax'], bool) or
                 not isinstance(options['fmax'], (int, float)) or
                 not np.isfinite(options['fmax']) or options['fmax'] <= 0)):
            raise ValueError(f'{ident}: fmax must be a positive finite number.')
        if 'steps' in options:
            _positive_integer(options['steps'], f'{ident} steps')
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
    _atoms(spec.get('molecule'))
    limits = spec.get('limits')
    if not isinstance(limits, dict):
        raise ValueError('limits must be an object.')
    _positive_integer(limits.get('max_nodes'), 'max_nodes')
    for key in ('max_cores_per_node', 'max_memory_mb_per_node'):
        if key in limits:
            _positive_integer(limits[key], key)
    tasks = spec.get('tasks')
    if not isinstance(tasks, list) or not tasks:
        raise ValueError('tasks must be a nonempty list.')
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
            pinned_methyl_uhf = (
                parser.get('kind') == 'mrcc_energy'
                and reference == 'UHF'
                and task.get('source_reference_validation') ==
                'methyl-qz-2017-UHF'
                and spec.get('intent', {}).get(
                    'source_reference_validation') ==
                'methyl-qz-2017-UHF'
                and spec.get('intent', {}).get(
                    'literature_benchmark', {}).get('name') ==
                'methyl-qz-2017'
                and multiplicity == 2)
            if (not pinned_methyl_uhf
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
            if (source_task is None or source_task.get('kind') != 'ase_optimize'
                    or not isinstance(profile, dict)
                    or not isinstance(keywords, dict)
                    or profile.get('calculator', '').lower() not in ('gaussian', 'gauss')
                    or profile.get('method', '').casefold() != parser['method'].casefold()
                    or profile.get('basis', '').casefold() != parser['basis'].casefold()
                    or keywords.get('EmpiricalDispersion', '').casefold()
                    != parser.get('dispersion', '').casefold()):
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


def _atomic_json(path, data):
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
    temp.replace(path)


def _load(run_dir):
    run_dir = Path(run_dir).resolve()
    spec = json.loads((run_dir / 'workflow.json').read_text())
    validate_spec(spec, allow_obsolete=True)
    state = json.loads((run_dir / 'state.json').read_text())
    checksum = hashlib.sha256((run_dir / 'workflow.json').read_bytes()).hexdigest()
    if state['spec_sha256'] != checksum:
        raise RuntimeError('Workflow changed after preparation; start a new run.')
    return run_dir, spec, state


def _task_dependencies(task):
    dependencies = set(task.get('depends_on', []))
    if task.get('geometry_from', 'initial') != 'initial':
        dependencies.add(task['geometry_from'])
    return dependencies


def _render_input(template, atoms, molecule, resources):
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

    return _TOKEN.sub(replace, template)


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
    lines.append(f'{shlex.quote(python)} -m kinbot.anl.dispatch run-task task.json')
    return '\n'.join(lines) + '\n'


def _stage_task(run_dir, spec, state, task):
    ident = task['id']
    directory = run_dir / 'tasks' / ident
    if directory.exists():
        raise RuntimeError(f'{ident}: existing task directory needs manual review.')
    source = task.get('geometry_from', 'initial')
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
                                 spec['molecule'], task['resources'])
        if _backend(task) == 'cfour':
            rendered = normalize_cfour_zmat(rendered)
        (directory / task['input_name']).write_text(rendered)
    _atomic_json(directory / 'task.json', {
        'schema': 1, 'task': task,
        'molecule': {'symbols': initial.get_chemical_symbols(),
                     'charge': spec['molecule'].get('charge', 0),
                     'multiplicity': spec['molecule'].get('multiplicity', 1)},
        'geometry_sha256': _geometry_hash(read(directory / 'geometry.xyz')),
    })
    (directory / 'job.slurm').write_text(
        _slurm_script(task, directory, state['python']))
    state['tasks'][ident] = {
        'status': 'staged', 'geometry_from': source,
        'geometry_sha256': _geometry_hash(atoms),
        'task_sha256': _file_hash(directory / 'task.json'),
        'job_sha256': _file_hash(directory / 'job.slurm'),
    }
    if task['kind'] == 'external':
        state['tasks'][ident]['input_sha256'] = _file_hash(
            directory / task['input_name'])


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
        if not _task_dependencies(task):
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
    stdin = (directory / task['input_name']).open('rb') if task.get('stdin') else None
    try:
        with (directory / task.get('stdout', 'stdout.txt')).open('wb') as stdout, \
                (directory / task.get('stderr', 'stderr.txt')).open('wb') as stderr:
            result = subprocess.run(command, cwd=directory, stdin=stdin,
                                    stdout=stdout, stderr=stderr, check=False,
                                    env=child_env)
    finally:
        if stdin:
            stdin.close()
        cleanup_qc_runtime(runtime)
    if result.returncode:
        error_file = directory / task.get('stderr', 'stderr.txt')
        detail = error_file.read_text(errors='replace').strip()[-1200:]
        suffix = f' Last stderr: {detail}' if detail else ''
        raise RuntimeError(f'Program exited with status {result.returncode}.{suffix}')
    _validate_external_outputs(directory, task)
    details = {'command': command, 'returncode': result.returncode, **runtime}
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


def _validate_external_outputs(directory, task):
    missing = [name for name in task['required_outputs']
               if not (directory / name).is_file()
               or not (directory / name).stat().st_size]
    if missing:
        raise RuntimeError(f'Missing or empty required output: {missing}')
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


def _parser_failure(task, execution):
    """Whether a failed record is eligible for a no-execution reparse."""
    return (task.get('kind') == 'external'
            and bool(task.get('result_parser'))
            and bool(task.get('success_marker'))
            and execution.get('status') == 'failed'
            and 'parse_result(' in execution.get('traceback', ''))


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
    """Reparse one hash-checked native success without changing state.json."""
    ident = task['id']
    directory = run_dir / 'tasks' / ident
    outcome = directory / 'execution.json'
    if not _parser_failure(task, previous):
        raise ValueError(f'{ident}: failure did not occur during result parsing.')
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
    details = {
        'command': _external_command(task),
        'returncode': 0,
        'returncode_source': 'original_parser_path_and_native_success_marker',
        'parsed_result': parsed,
        'reparsed_without_execution': True,
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
    """Recover a successful external calculation rejected only by its parser."""
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
        from kinbot.anl.results import parse_result
        parsed = parse_result(
            (directory / requested['file']).read_text(errors='replace'),
            requested)
        if parsed != previous.get('details', {}).get('parsed_result'):
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
    from sella import Sella
    from kinbot.ase_modules.calculators.factory import build_calculator

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
    if len(atoms) > 2 and profile.optimizer == 'sella':
        optimizer = Sella(atoms, order=0,
                          **common, **options.get('sella_kwargs', {}))
    else:
        optimizer = BFGS(atoms, **common)
    converged = optimizer.run(fmax=options.get('fmax', 0.03),
                              steps=options.get('steps', 100))
    if not converged:
        raise RuntimeError('ASE geometry optimization did not converge.')
    energy_ev = float(atoms.get_potential_energy())
    write(directory / 'final.xyz', atoms)
    return {'energy_ev': energy_ev, 'energy_hartree': energy_ev / Hartree,
            'optimizer': profile.optimizer if len(atoms) > 2 else 'bfgs_or_atom',
            'generated_files': getattr(atoms.calc, 'generated_files', [])}


def run_task(task_file):
    task_file = Path(task_file).resolve()
    directory = task_file.parent
    record = json.loads(task_file.read_text())
    task = record['task']
    outcome = directory / 'execution.json'
    if outcome.exists():
        raise RuntimeError('Task already has execution.json; review before rerunning.')
    result = {'schema': 1, 'task_id': task['id'],
              'geometry_sha256': record['geometry_sha256']}
    try:
        run_dir, spec, state = _load(directory.parents[1])
        task_spec = next(item for item in spec['tasks'] if item['id'] == task['id'])
        if task != task_spec or task['id'] != directory.name:
            raise RuntimeError('Task record does not match prepared workflow.')
        _verify_stage_files(run_dir, task, state['tasks'][task['id']])
        atoms = read(directory / 'geometry.xyz')
        if _geometry_hash(atoms) != record['geometry_sha256']:
            raise RuntimeError('Staged geometry changed after preparation.')
        details = (_run_external(directory, record) if task['kind'] == 'external'
                   else _run_ase_optimize(directory, record))
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
        else:
            names |= {'final.xyz', 'optimization.log'}
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
    result = subprocess.run(['squeue', '-h', '-j', str(job_id), '-o', '%i'],
                            capture_output=True, text=True, check=False)
    if result.returncode:
        # A completed job can disappear from slurmctld before its execution
        # record is reconciled. Some Slurm versions return an error rather
        # than an empty queue for this specific case.
        if re.search(r'\binvalid job id specified\b', result.stderr, re.I):
            return False
        raise RuntimeError(f'squeue failed: {result.stderr.strip()}')
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


def preflight(run_dir, _site_refresh_attempts=None):
    """Check login-node tools and the same site setup sourced by batch jobs."""
    run_dir, spec, state = _load(run_dir)
    refresh_attempts = (set() if _site_refresh_attempts is None
                        else _site_refresh_attempts)
    for name in ('sbatch', 'squeue', 'bash'):
        if not shutil.which(name):
            raise RuntimeError(f'Preflight: {name} is unavailable on PATH.')
    programs = set()
    env_files = set()
    for task in spec['tasks']:
        entry = state['tasks'].get(task['id'])
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
        result = subprocess.run(['bash', '-c', '\n'.join(lines)],
                                cwd=run_dir, capture_output=True, text=True,
                                check=False)
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
                        run_dir, _site_refresh_attempts=refresh_attempts)
                    refreshed['site_setup_refreshed'] = True
                    return refreshed
            raise RuntimeError(f"{task['id']}: preflight failed after site setup: "
                               + detail)
    checked = 0
    for task in spec['tasks']:
        entry = state['tasks'].get(task['id'])
        if not entry or entry['status'] != 'staged':
            continue
        result = subprocess.run(['sbatch', '--test-only', 'job.slurm'],
                                cwd=run_dir / 'tasks' / task['id'],
                                capture_output=True, text=True, check=False)
        if result.returncode:
            raise RuntimeError(f"{task['id']}: Slurm rejected test-only batch script: "
                               + (result.stderr.strip() or result.stdout.strip()))
        checked += 1
    return {'programs': sorted(programs), 'file_variables': sorted(env_files),
            'exclusive_jobs_checked': len(state['tasks']),
            'slurm_scripts_tested': checked}


def _clear_blocked_tasks(state):
    """Discard derived dependency blocks so the graph can be reevaluated."""
    for task_id in [
            task_id for task_id, entry in state['tasks'].items()
            if entry.get('status') == 'blocked']:
        state['tasks'].pop(task_id)


def retry_failed(run_dir, ident):
    """Archive one failed attempt and restage it without rerunning siblings."""
    run_dir, spec, state = _load(run_dir)
    by_id = {task['id']: task for task in spec['tasks']}
    if ident not in by_id or state['tasks'].get(ident, {}).get('status') != 'failed':
        raise ValueError(f'{ident}: only a failed, staged task can be retried.')
    entry = state['tasks'][ident]
    if entry.get('job_id') and _job_active(entry['job_id']):
        raise RuntimeError(f'{ident}: Slurm job is still active.')
    if any(other != ident and other in state['tasks']
           and state['tasks'][other].get('status') != 'blocked'
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


def recover_gaussian_vpt2_symmetry(run_dir, ident):
    """Retry the Gaussian framework-group VPT2 failure in a fixed C1 group.

    Gaussian 16 can classify the input and Eckart-oriented structures with
    different framework groups. Link 717 then aborts with an explicit
    framework-group inconsistency. This migration is intentionally narrow: it
    accepts only that native error with named, different old and new groups,
    replaces ``NoSymm`` (or an implicit symmetry setting) with the documented
    ``Symmetry=(PG=C1)`` cap, keeps the same accepted L2 geometry and method,
    and archives the entire failed attempt before staging the retry.
    """
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
    if entry.get('status') != 'failed':
        raise ValueError(f'{ident}: only a failed Gaussian VPT2 task can recover.')
    if entry.get('job_id') and _job_active(entry['job_id']):
        raise RuntimeError(f'{ident}: Slurm job is still active.')
    if any(other != ident and other in state['tasks']
           and state['tasks'][other].get('status') != 'blocked'
           and ident in _task_dependencies(other_task)
           for other, other_task in by_id.items()):
        raise RuntimeError(f'{ident}: downstream task already exists; review run.')
    _verify_stage_files(run_dir, task, entry)
    directory = run_dir / 'tasks' / ident
    native = directory / parser['file']
    if not native.is_file():
        raise RuntimeError(f'{ident}: native Gaussian output is missing.')
    output = native.read_text(errors='replace')
    required = (
        'ERROR: Inconsistency found in framework group definition:',
        'Error termination via Lnk1e',
    )
    framework_change = re.search(
        r'New:\s*([A-Za-z0-9]+)\s*-\s*Old:\s*([A-Za-z0-9]+)\b',
        output, re.IGNORECASE)
    if (not all(marker in output for marker in required)
            or framework_change is None
            or framework_change.group(1).casefold()
            == framework_change.group(2).casefold()
            or 'Normal termination of Gaussian' in output):
        raise RuntimeError(
            f'{ident}: native output is not the framework-group VPT2 failure.')
    replacement = deepcopy(task)
    original_template = replacement['input_template']
    if re.search(r'(?i)Symm(?:etry)?\s*=\s*\(\s*PG\s*=\s*C1\s*\)',
                 original_template):
        raise RuntimeError(
            f'{ident}: input already fixes the Gaussian framework to C1.')
    occurrences = len(re.findall(
        r'(?i)(?<!\S)NoSymm(?:etry)?(?!\S)', original_template))
    if occurrences > 1:
        raise RuntimeError(f'{ident}: input has multiple NoSymm keywords.')
    if occurrences == 1:
        template = re.sub(
            r'(?i)(?<!\S)NoSymm(?:etry)?(?!\S)', 'Symmetry=(PG=C1)',
            original_template, count=1)
    else:
        template, count = re.subn(
            r'(?i)(\bFreq\s*=\s*Anharmonic\b)',
            r'\1 Symmetry=(PG=C1)', original_template, count=1)
        if count != 1:
            raise RuntimeError(
                f'{ident}: input has no unique Freq=Anharmonic keyword.')
    replacement['input_template'] = template
    replacement['recovery'] = {
        'kind': 'gaussian_vpt2_framework_group',
        'change': ('fixed Gaussian framework at C1 with the documented '
                   'Symmetry=(PG=C1) option'),
        'observed_framework_change': framework_change.group(0),
        'geometry_changed': False,
    }
    spec['tasks'][spec['tasks'].index(task)] = replacement
    assign_partitions(spec)
    validate_spec(spec)

    old = directory
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
        'Gaussian VPT2 framework-group recovery with Symmetry=(PG=C1)'
    _clear_blocked_tasks(state)
    _atomic_json(run_dir / 'state.json', state)
    return archive


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
    if entry.get('status') != 'failed':
        raise ValueError(f'{ident}: only a failed CFOUR task can be migrated.')
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


def resume_mrcc_failed(run_dir, ident, *, walltime='7-00:00:00',
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
    if entry.get('status') != 'failed':
        raise ValueError(f'{ident}: only a failed MRCC task can be resumed.')
    if entry.get('job_id') and _job_active(entry['job_id']):
        raise RuntimeError(f'{ident}: Slurm job is still active.')
    if any(other != ident and other in state['tasks']
           and state['tasks'][other].get('status') != 'blocked'
           and ident in _task_dependencies(other_task)
           for other, other_task in by_id.items()):
        raise RuntimeError(f'{ident}: downstream task already exists; review run.')
    directory = run_dir / 'tasks' / ident
    _verify_stage_files(run_dir, task, entry)
    amplitude = directory / 'fort.16'
    if not amplitude.is_file() or amplitude.stat().st_size == 0:
        raise RuntimeError(f'{ident}: nonempty MRCC fort.16 is required to resume.')
    if (not isinstance(walltime, str) or not _TIME.fullmatch(walltime)):
        raise ValueError('MRCC resume walltime must be HH:MM:SS or D-HH:MM:SS.')
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
    if partition is None:
        resources.pop('partition', None)
    else:
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
        'checkpoint': 'fort.16',
    }
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
    run_dir, spec, state = _load(run_dir)
    by_id = {task['id']: task for task in spec['tasks']}
    if submit_only is not None:
        submit_only = set(submit_only)
        unknown = submit_only - by_id.keys()
        if unknown:
            raise ValueError(f"Unknown submission task(s): {', '.join(sorted(unknown))}")

    # Reconcile two deterministic non-QC cases before dependency handling.
    # First, a newer parser may accept a native calculation that an older
    # installed KinBot rejected.  Second, an atom has no geometry coordinate
    # to optimize, so its geometry nodes are exact identity operations and do
    # not belong in the scheduler at all.
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
        if (entry.get('status') == 'failed' and previous is not None
                and _parser_failure(task, previous)):
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
        outcome = run_dir / 'tasks' / ident / 'execution.json'
        if outcome.exists():
            execution = json.loads(outcome.read_text())
            if entry['status'] == 'submitted' and _job_active(entry['job_id']):
                continue
            if _parser_failure(task, execution):
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
        elif entry['status'] == 'submitted' and not _job_active(entry['job_id']):
            entry['status'] = 'failed'
            entry['error'] = 'Slurm job disappeared without execution.json.'
    for task in spec['tasks']:
        ident = task['id']
        dependencies = _task_dependencies(task)
        failed_dependencies = [
            dependency for dependency in dependencies
            if state['tasks'].get(dependency, {}).get('status') in
            ('failed', 'blocked')]
        entry = state['tasks'].get(ident)
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
        if failed_dependencies:
            state['tasks'][ident] = {
                'status': 'blocked',
                'geometry_from': task.get('geometry_from', 'initial'),
                'blocked_by': failed_dependencies,
            }
        elif all(state['tasks'].get(dep, {}).get('status') == 'complete'
                 for dep in dependencies):
            _stage_task(run_dir, spec, state, task)
    _atomic_json(run_dir / 'state.json', state)
    if submit:
        active = sum(entry['status'] == 'submitted' for entry in state['tasks'].values())
        for task in spec['tasks']:
            ident = task['id']
            entry = state['tasks'].get(ident)
            if active >= spec['limits']['max_nodes']:
                break
            if submit_only is not None and ident not in submit_only:
                continue
            if not entry or entry['status'] != 'staged':
                continue
            entry['status'] = 'submitting'
            _atomic_json(run_dir / 'state.json', state)
            directory = run_dir / 'tasks' / ident
            try:
                response = subprocess.run(['sbatch', '--parsable', 'job.slurm'],
                                          cwd=directory, capture_output=True,
                                          text=True, check=False)
            except OSError:
                entry['status'] = 'staged'
                _atomic_json(run_dir / 'state.json', state)
                raise
            if response.returncode:
                entry['status'] = 'staged'
                _atomic_json(run_dir / 'state.json', state)
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
    resume_mrcc.add_argument('--walltime', default='7-00:00:00')
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
        run_dir = Path(args.run_dir).resolve()
        with (run_dir / 'drive.lock').open('w') as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise RuntimeError('Another driver is already managing this run.') from exc
            while True:
                state = advance(run_dir, submit=True, submit_only=args.only)
                _, spec, _ = _load(run_dir)
                summary = _summary(spec, state)
                print(json.dumps(summary, sort_keys=True), flush=True)
                terminal = {'complete', 'failed', 'blocked'}
                if all(value in terminal for value in summary.values()):
                    return 0 if all(value == 'complete'
                                    for value in summary.values()) else 1
                if args.once:
                    return int(any(value in ('failed', 'blocked')
                                   for value in summary.values()))
                if args.interval < 1 or args.interval > 60:
                    parser.error('--interval must be 1 through 60 seconds')
                time.sleep(args.interval)


if __name__ == '__main__':
    sys.exit(main())
