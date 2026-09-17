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
    """Budget Molpro node memory after documented per-process overhead."""
    total_mw = int(resources['memory_mb'] * 0.85 / 8) - 200 * resources['cores']
    if total_mw < 32:
        raise ValueError('Molpro memory allocation leaves no usable stack/GA '
                         'after per-process overhead; request more node memory '
                         'or fewer cores.')
    return total_mw


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


def _validate_task(task, ids, limits):
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
    if cores > limits['max_cores_per_node'] or memory > limits['max_memory_mb_per_node']:
        raise ValueError(f'{ident}: task exceeds per-node resource limits.')
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
    else:
        _basename(task.get('backend'), f'{ident} backend')
        _basename(task.get('input_name'), f'{ident} input_name')
        if not isinstance(task.get('input_template'), str) or not task['input_template'].strip():
            raise ValueError(f'{ident}: input_template is required.')
        unknown_tokens = set(_TOKEN.findall(task['input_template'])) - _INPUT_TOKENS
        if unknown_tokens:
            raise ValueError(f'{ident}: unknown input tokens {sorted(unknown_tokens)}.')
        command = task.get('command')
        if not isinstance(command, list) or not command or any(
                not isinstance(arg, str) or not arg for arg in command):
            raise ValueError(f'{ident}: command must be an argument list.')
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
                or (task['backend'].lower() == 'molpro' and stdout_name in outputs)):
            raise ValueError(f'{ident}: input, stdout, and stderr filenames must not collide.')
        backend = task['backend'].lower()
        if backend == 'molpro':
            _molpro_total_mw(resources)
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
                    or filename in _RESERVED_FILES:
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


def validate_spec(spec):
    """Validate a declarative graph before creating or submitting any job."""
    if not isinstance(spec, dict) or spec.get('schema') != 1:
        raise ValueError('Workflow must be a schema 1 object.')
    _basename(spec.get('name'), 'Workflow name')
    _atoms(spec.get('molecule'))
    limits = spec.get('limits')
    if not isinstance(limits, dict):
        raise ValueError('limits must be an object.')
    for key in ('max_nodes', 'max_cores_per_node', 'max_memory_mb_per_node'):
        _positive_integer(limits.get(key), key)
    tasks = spec.get('tasks')
    if not isinstance(tasks, list) or not tasks:
        raise ValueError('tasks must be a nonempty list.')
    ids = [_basename(task.get('id') if isinstance(task, dict) else None, 'Task id')
           for task in tasks]
    if len(ids) != len(set(ids)):
        raise ValueError('Task ids must be unique.')
    by_id = dict(zip(ids, tasks))
    for task in tasks:
        _validate_task(task, by_id, limits)
        source = task.get('geometry_from', 'initial')
        if source != 'initial' and not by_id[source].get('geometry_output'):
            raise ValueError(f"{task['id']}: geometry source has no geometry output.")
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
    validate_spec(spec)
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
        'MRCC_XYZ': f'{len(atoms)}\n{coordinates}',
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
    if values.get('calculator', '').lower() in ('gaussian', 'gauss'):
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


def _slurm_script(task, directory, python):
    resources = task['resources']
    lines = [
        '#!/usr/bin/env bash',
        f"#SBATCH --job-name=kb-{task['id']}",
        '#SBATCH --nodes=1',
        f"#SBATCH --cpus-per-task={resources['cores']}",
        f"#SBATCH --mem={resources['memory_mb']}M",
        f"#SBATCH --time={resources['walltime']}",
        '#SBATCH --exclusive',
        '#SBATCH --output=slurm.stdout',
        '#SBATCH --error=slurm.stderr',
    ]
    if resources.get('partition'):
        lines.append(f"#SBATCH --partition={resources['partition']}")
    lines += ['', 'set -euo pipefail', f'cd {shlex.quote(str(directory))}',
              f"export OMP_NUM_THREADS={resources['cores']}",
              'source ../../site_setup.sh']
    lines += task.get('setup', [])
    lines.append(f'{shlex.quote(python)} -m kinbot.anl.dispatch run-task task.json')
    return '\n'.join(lines) + '\n'


def _stage_task(run_dir, spec, state, task):
    ident = task['id']
    directory = run_dir / 'tasks' / ident
    if directory.exists():
        raise RuntimeError(f'{ident}: existing task directory needs manual review.')
    source = task.get('geometry_from', 'initial')
    atoms = (_atoms(spec['molecule']) if source == 'initial' else
             read(run_dir / 'tasks' / source /
                  next(item['geometry_output'] for item in spec['tasks']
                       if item['id'] == source)))
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
        (directory / task['input_name']).write_text(rendered)
    _atomic_json(directory / 'task.json', {
        'schema': 1, 'task': task,
        'molecule': {'symbols': initial.get_chemical_symbols(),
                     'charge': spec['molecule'].get('charge', 0),
                     'multiplicity': spec['molecule'].get('multiplicity', 1)},
        'geometry_sha256': _geometry_hash(read(directory / 'geometry.xyz')),
    })
    (directory / 'job.slurm').write_text(_slurm_script(task, directory, sys.executable))
    state['tasks'][ident] = {'status': 'staged', 'geometry_from': source,
                             'geometry_sha256': _geometry_hash(atoms)}


def prepare(spec_path, run_dir):
    spec = validate_spec(json.loads(Path(spec_path).read_text()))
    run_dir = Path(run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=False)
    (run_dir / 'tasks').mkdir()
    (run_dir / 'site_setup.sh').write_text(
        '#!/usr/bin/env bash\n'
        '# Load licensed QC modules and set CFOUR_GENBAS on this site.\n')
    _atomic_json(run_dir / 'workflow.json', spec)
    state = {'schema': 1,
             'spec_sha256': hashlib.sha256((run_dir / 'workflow.json').read_bytes()).hexdigest(),
             'tasks': {}}
    for task in spec['tasks']:
        if not _task_dependencies(task):
            _stage_task(run_dir, spec, state, task)
    _atomic_json(run_dir / 'state.json', state)
    return run_dir


def _run_external(directory, record):
    task = record['task']
    for filename, env_name in task.get('files_from_env', {}).items():
        source = os.environ.get(env_name)
        if not source or not Path(source).is_file():
            raise RuntimeError(f'{filename} requires existing file in ${env_name}.')
        shutil.copyfile(source, directory / filename)
    # -M is a total-node limit. Molpro documents at least 200 MW of extra
    # program memory per compute process, beyond the stack/GA allocation.
    molpro_total_mw = (_molpro_total_mw(task['resources'])
                       if task['backend'].lower() == 'molpro' else 0)
    command = [arg.replace('{cores}', str(task['resources']['cores']))
               .replace('{input}', task['input_name'])
               .replace('{molpro_total_mw}', str(molpro_total_mw))
               for arg in task['command']]
    stdin = (directory / task['input_name']).open('rb') if task.get('stdin') else None
    try:
        with (directory / task.get('stdout', 'stdout.txt')).open('wb') as stdout, \
                (directory / task.get('stderr', 'stderr.txt')).open('wb') as stderr:
            result = subprocess.run(command, cwd=directory, stdin=stdin,
                                    stdout=stdout, stderr=stderr, check=False,
                                    env={**os.environ,
                                         'OMP_NUM_THREADS': str(task['resources']['cores'])})
    finally:
        if stdin:
            stdin.close()
    if result.returncode:
        raise RuntimeError(f'Program exited with status {result.returncode}.')
    missing = [name for name in task['required_outputs']
               if not (directory / name).is_file()
               or not (directory / name).stat().st_size]
    if missing:
        raise RuntimeError(f'Missing or empty required output: {missing}')
    if task.get('success_marker'):
        marker = task['success_marker']
        output = directory / marker['file']
        if not output.is_file() or marker['contains'] not in output.read_text(errors='replace'):
            raise RuntimeError(f"Success marker absent from {marker['file']}.")
    return {'command': command, 'returncode': result.returncode}


def _run_ase_optimize(directory, record):
    from ase.optimize import BFGS
    from sella import Sella
    from kinbot.ase_modules.calculators.factory import build_calculator

    task = record['task']
    atoms = read(directory / 'geometry.xyz')
    # UMA/OMol reads these from Atoms.info; XYZ does not preserve them.
    atoms.info['charge'] = record['molecule']['charge']
    atoms.info['spin'] = record['molecule']['multiplicity']
    profile = _runtime_profile(task)
    atoms.calc = build_calculator(profile, directory, {
        'name': task['id'], 'charge': record['molecule']['charge'],
        'mult': record['molecule']['multiplicity']})
    options = task.get('optimizer', {})
    if len(atoms) > 1:
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
            'optimizer': profile.optimizer if len(atoms) > 2 else 'bfgs_or_atom'}


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
        artifacts = {name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
                     for name in sorted(names) if (directory / name).is_file()}
        result.update(status='executed', details=details, artifacts=artifacts)
    except Exception as exc:
        result.update(status='failed', error=str(exc),
                      traceback=traceback.format_exc())
    _atomic_json(outcome, result)
    if result['status'] != 'executed':
        raise RuntimeError(f"{task['id']}: {result['error']}")
    return result


def _job_active(job_id):
    result = subprocess.run(['squeue', '-h', '-j', str(job_id), '-o', '%i'],
                            capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError(f'squeue failed: {result.stderr.strip()}')
    return str(job_id) in result.stdout.split()


def advance(run_dir, submit=False):
    """Complete finished tasks, stage ready tasks, and optionally submit jobs."""
    run_dir, spec, state = _load(run_dir)
    by_id = {task['id']: task for task in spec['tasks']}
    for ident, entry in state['tasks'].items():
        if entry['status'] == 'submitting':
            raise RuntimeError(f'{ident}: submission outcome uncertain; reconcile Slurm job manually.')
        if entry['status'] not in ('staged', 'submitted'):
            continue
        outcome = run_dir / 'tasks' / ident / 'execution.json'
        if outcome.exists():
            execution = json.loads(outcome.read_text())
            if execution['geometry_sha256'] != json.loads(
                    (run_dir / 'tasks' / ident / 'task.json').read_text())['geometry_sha256']:
                raise RuntimeError(f'{ident}: execution geometry hash mismatch.')
            if entry['status'] == 'submitted' and _job_active(entry['job_id']):
                continue
            geometry_output = by_id[ident].get('geometry_output')
            if execution['status'] == 'executed' and geometry_output:
                final = run_dir / 'tasks' / ident / geometry_output
                if not final.is_file():
                    raise RuntimeError(f'{ident}: accepted geometry file is missing.')
                actual_hash = _geometry_hash(read(final))
                if actual_hash != execution.get('details', {}).get('final_geometry_sha256'):
                    raise RuntimeError(f'{ident}: accepted geometry hash mismatch.')
            entry['status'] = ('complete' if execution['status'] == 'executed'
                               else 'failed')
            entry['execution'] = execution['status']
            if geometry_output and execution['status'] == 'executed':
                entry['final_geometry_sha256'] = actual_hash
        elif entry['status'] == 'submitted' and not _job_active(entry['job_id']):
            entry['status'] = 'failed'
            entry['error'] = 'Slurm job disappeared without execution.json.'
    if any(entry['status'] == 'failed' for entry in state['tasks'].values()):
        _atomic_json(run_dir / 'state.json', state)
        return state
    for task in spec['tasks']:
        ident = task['id']
        if ident in state['tasks']:
            continue
        if all(state['tasks'].get(dep, {}).get('status') == 'complete'
               for dep in _task_dependencies(task)):
            _stage_task(run_dir, spec, state, task)
    _atomic_json(run_dir / 'state.json', state)
    if submit:
        active = sum(entry['status'] == 'submitted' for entry in state['tasks'].values())
        for task in spec['tasks']:
            ident = task['id']
            entry = state['tasks'].get(ident)
            if active >= spec['limits']['max_nodes']:
                break
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


def main(argv=None):
    parser = argparse.ArgumentParser(description='Stage and drive exclusive QC jobs')
    sub = parser.add_subparsers(dest='action', required=True)
    stage = sub.add_parser('prepare')
    stage.add_argument('spec')
    stage.add_argument('run_dir')
    drive = sub.add_parser('drive')
    drive.add_argument('run_dir')
    drive.add_argument('--once', action='store_true')
    drive.add_argument('--interval', type=int, default=20)
    status = sub.add_parser('status')
    status.add_argument('run_dir')
    run = sub.add_parser('run-task')
    run.add_argument('task_file')
    args = parser.parse_args(argv)
    if args.action == 'prepare':
        print(prepare(args.spec, args.run_dir))
    elif args.action == 'run-task':
        run_task(args.task_file)
    elif args.action == 'status':
        _, spec, state = _load(args.run_dir)
        print(json.dumps(_summary(spec, state), indent=2))
    else:
        run_dir = Path(args.run_dir).resolve()
        with (run_dir / 'drive.lock').open('w') as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise RuntimeError('Another driver is already managing this run.') from exc
            while True:
                state = advance(run_dir, submit=True)
                _, spec, _ = _load(run_dir)
                summary = _summary(spec, state)
                print(json.dumps(summary, sort_keys=True), flush=True)
                if all(value == 'complete' for value in summary.values()):
                    return 0
                if 'failed' in summary.values():
                    return 1
                if args.once:
                    return 0
                if args.interval < 1 or args.interval > 60:
                    parser.error('--interval must be 1 through 60 seconds')
                time.sleep(args.interval)


if __name__ == '__main__':
    sys.exit(main())
