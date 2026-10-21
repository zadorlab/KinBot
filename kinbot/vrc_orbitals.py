"""Molpro fragment optimization and native orbital cubes for VRC pivot points.

These auxiliary orbitals do not supply the VRC interaction energies.  Their
own dispatcher graph retains the native input, output, geometry and cube.
"""

from copy import copy
import hashlib
import json
from pathlib import Path
import re
import shutil

import numpy as np
from ase import Atoms
from ase.io import read
from ase.units import Bohr, Hartree

from kinbot.anl.tasks import MOLPRO_HEADER
from kinbot.molpro import validate_molpro_convergence
from kinbot.orbital_cube import read_orbital_cube
from kinbot.theory import fortran_threshold


GAUSSIAN_EXPORT_FILES = ('formchk.stdout', 'formchk.stderr', 'cubegen.stdout',
                         'cubegen.stderr', 'gaussian_orbital_export.json')


def copy_fragment_geometries(products):
    """Give a VRC path independent coordinates, retaining read-only metadata.

    Products can own conformer searches, QC routers and native molecule
    handles. Those live services must not be deep-copied. The MR scan changes
    only ``geom`` after atom ordering has been established; connectivity and
    the accepted thermochemistry remain shared and are never modified here.
    """
    copies = []
    for product in products:
        fragment = copy(product)
        fragment.geom = np.array(product.geom, dtype=float, copy=True)
        copies.append(fragment)
    return copies


def validate_fragment_controls(parameters):
    if parameters.get('vrc_tst_fragment_geometry', 'legacy') not in ('legacy', 'asymptotic_scan'):
        raise ValueError('vrc_tst_fragment_geometry must be legacy or asymptotic_scan.')
    if parameters.get('vrc_tst_pivot_level', 'L2') not in ('L2', 'custom'):
        raise ValueError('vrc_tst_pivot_level must be L2 or custom.')


def asymptotic_fragments(parameters):
    validate_fragment_controls(parameters)
    return parameters.get('vrc_tst_fragment_geometry', 'legacy') == 'asymptotic_scan'


def fixed_pivot_profile(parameters):
    """Resolve orbital-only single points independently of MR geometry/energy."""
    from copy import deepcopy
    from dataclasses import asdict
    from kinbot.theory import resolve_profiles, _legacy_profile
    if parameters.get('vrc_tst_pivot_level', 'L2') == 'L2':
        profiles = resolve_profiles(parameters)
        profile = asdict(profiles['l2']) if profiles else _legacy_profile(parameters, 'l2')
    else:
        profile = {'calculator': pivot_backend(parameters),
                   'method': parameters['vrc_tst_pivot_method'],
                   'basis': parameters['vrc_tst_pivot_basis'],
                   'calculator_kwargs': deepcopy(parameters.get('vrc_tst_calc_kwargs') or {})}
    if profile['calculator'] not in ('gaussian', 'molpro'):
        raise ValueError('Fixed-geometry orbital cubes currently support Gaussian or Molpro; '
                         'set vrc_tst_pivot_level=custom to select either explicitly.')
    return profile


def fixed_pivot_spec(product, parameters):
    """One scheduled orbital calculation; never optimize the MR fragment again."""
    from io import StringIO
    import shlex
    from kinbot.ase_modules.io.gaussian import write_gaussian_in, _problem_methods
    from kinbot.vrc_tst_scan import vrc_dispatch_resources
    profile = fixed_pivot_profile(parameters)
    atoms = Atoms(product.atom, positions=product.geom)
    if len(atoms) < 2:
        raise ValueError('Atomic fragments use atom-centered pivots and need no orbital job.')
    if profile['calculator'] == 'molpro':
        controls = dict(parameters, vrc_tst_pivot_method=profile['method'],
                        vrc_tst_pivot_basis=profile['basis'])
        return molpro_pivot_spec(product, controls, optimize=False)
    kwargs = dict(profile.get('calculator_kwargs', {}))
    if profile['method'].lower() in _problem_methods:
        raise ValueError('Composite Gaussian methods cannot supply fixed-geometry L2 pivot orbitals; '
                         'select an explicit orbital method with vrc_tst_pivot_level=custom.')
    if re.search(r'\b(?:opt|freq|irc|geom|force)\b', str(kwargs.get('extra', '')), re.I):
        raise ValueError('L2 extra route contains a geometry/derivative instruction; '
                         'configure a fixed-geometry custom pivot method.')
    # Job-role/old-checkpoint options must not turn this into an optimization
    # or replace the accepted input geometry. Physical SCF/functional/grid
    # settings remain those of L2.
    for key in list(kwargs):
        if key.lower() in {'opt', 'freq', 'force', 'irc', 'geom', 'guess', 'chk',
                           'oldchk', 'label', 'directory', 'nprocshared', 'mem',
                           'command', 'scf_recovery'}:
            kwargs.pop(key)
    text = StringIO()
    write_gaussian_in(text, atoms, properties=['energy'], method=profile['method'],
        basis=profile['basis'], charge=int(product.charge), mult=int(product.mult),
        chk='pivot.chk', nprocshared='{{CORES}}', mem='{{WORK_MEMORY_MB}}MB', **kwargs)
    task = {'id': 'pivot', 'kind': 'external', 'backend': 'gaussian',
        'geometry_from': 'initial', 'resources': vrc_dispatch_resources(parameters),
        'input_name': 'pivot.com', 'input_template': text.getvalue(),
        'command': shlex.split(profile.get('command') or 'g16'), 'stdin': 'pivot.com',
        'stdout': 'pivot.log', 'stderr': 'launcher.stderr',
        'required_executables': ['formchk', 'cubegen'],
        'required_outputs': ['pivot.log', 'pivot.chk', 'pivot.fchk', 'homo.cube'],
        'success_marker': {'file': 'pivot.log', 'contains': 'Normal termination'},
        'gaussian_orbital_export': {'checkpoint': 'pivot.chk', 'formatted': 'pivot.fchk',
                                    'cube': 'homo.cube', 'log': 'pivot.log',
                                    'orbital': 'AMO=HOMO' if product.mult > 1 else 'MO=HOMO'}}
    return {'schema': 1, 'name': 'kinbot-vrc-fixed-fragment-orbital',
        'molecule': {'symbols': atoms.get_chemical_symbols(), 'positions': atoms.positions.tolist(),
                     'charge': int(product.charge), 'multiplicity': int(product.mult)},
        'limits': {'max_nodes': 1}, 'tasks': [task],
        'intent': {'cap_walltime_to_explicit_partition': True, 'fixed_geometry': True,
                   'orbital_only': True, 'pivot_profile': profile}}


def export_gaussian_cube(directory, settings, environment):
    """Run conversion on the compute node, inside the native local workspace."""
    import subprocess
    from kinbot.anl.dispatch import CalculationUnavailableError
    from kinbot.execution import atomic_json
    directory = Path(directory)
    log = directory/settings['log']
    text = log.read_text(errors='replace')
    if ('Normal termination' not in text or re.search(
            r'Convergence criterion not met|Convergence failure|Error termination', text, re.I)):
        raise ValueError('Gaussian orbital calculation did not converge normally.')
    # cubegen's nprocs argument requests one processor. Do not inherit the
    # preceding Gaussian SCF's OpenMP allocation or mutate its environment.
    environment = dict(environment, OMP_NUM_THREADS='1')
    commands = [['formchk', settings['checkpoint'], settings['formatted']],
                ['cubegen', '1', settings['orbital'], settings['formatted'], settings['cube'], '-4']]
    evidence = {'schema': 1, 'directory': str(directory), 'settings': dict(settings),
                'environment': {key: environment.get(key) for key in (
                    'GAUSS_MEMDEF', 'GAUSS_PDEF', 'GAUSS_CDEF', 'OMP_NUM_THREADS',
                    'GAUSS_SCRDIR', 'TMPDIR')}, 'commands': []}
    evidence_path = directory / 'gaussian_orbital_export.json'
    def excerpt(path):
        with path.open('rb') as stream:
            head = stream.read(2048)
            stream.seek(0, 2)
            size = stream.tell()
            if size > len(head):
                stream.seek(max(len(head), size - 1024))
                head += b'\n...\n' + stream.read(1024)
        return head.decode(errors='replace').strip()
    for command in commands:
        from kinbot.anl.site import _executable
        executable = (environment.get(f'KINBOT_{command[0].upper()}_EXECUTABLE')
                      or shutil.which(command[0], path=environment.get('PATH'))
                      or _executable(command[0], 'gaussian'))
        if executable is None:
            raise FileNotFoundError(f'Missing Gaussian cube utility {command[0]}.')
        stdout, stderr = (directory / (command[0] + suffix)
                          for suffix in ('.stdout', '.stderr'))
        record = {'command': [str(executable), *command[1:]],
                  'stdout': stdout.name, 'stderr': stderr.name}
        evidence['commands'].append(record)
        atomic_json(evidence_path, evidence)
        try:
            with stdout.open('wb') as out, stderr.open('wb') as err:
                result = subprocess.run(record['command'], cwd=directory, env=environment,
                                        stdout=out, stderr=err, check=False)
            record['returncode'] = result.returncode
        except OSError as error:
            record['error'] = str(error)
            raise CalculationUnavailableError(
                f'Could not start {command[0]}: {error}') from error
        finally:
            atomic_json(evidence_path, evidence)
        if result.returncode:
            raise CalculationUnavailableError(
                f'{command[0]} failed (exit {result.returncode}); '
                f'stdout: {excerpt(stdout)}; stderr: {excerpt(stderr)}')
    validate_cube_geometry(read_orbital_cube(directory/settings['cube']), read(directory/'geometry.xyz'))


def validate_cube_geometry(cube, atoms):
    if cube['atoms'].get_chemical_symbols() != atoms.get_chemical_symbols():
        raise ValueError('Orbital cube changed fragment atom order.')
    source = cube['atoms'].positions - cube['atoms'].positions.mean(axis=0)
    target = atoms.positions - atoms.positions.mean(axis=0)
    u, _, vh = np.linalg.svd(source.T @ target)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vh)
    if not np.allclose(source @ (u @ correction @ vh), target, atol=2.e-5, rtol=0):
        raise ValueError('Orbital cube differs from the fixed MR fragment geometry.')


def prepare_fixed_pivot(product, parameters, destination, geometry_evidence):
    from kinbot.anl.dispatch import _load, preflight_until_scheduler_ready
    from kinbot.vrc_tst_scan import _drive_vrc_dispatch
    from kinbot.execution import atomic_json, file_sha256
    spec = fixed_pivot_spec(product, parameters)
    destination = Path(destination)
    from kinbot.vrc_graph import select_graph
    graph = select_graph(Path('vrctst/molpro_pivots'), destination.name+'_fixed', spec)
    preflight_until_scheduler_ready(graph, interval=20.)
    code, summary = _drive_vrc_dispatch(graph, ['pivot'])
    if code or _load(graph)[2]['tasks']['pivot']['status'] != 'complete':
        raise RuntimeError(f'Fixed-geometry VRC orbital task failed at {graph}: {summary}')
    cube = graph/'tasks/pivot/homo.cube'
    validate_cube_geometry(read_orbital_cube(cube), Atoms(product.atom, positions=product.geom))
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix('.cube.tmp')
    shutil.copyfile(cube, temporary)
    temporary.replace(destination.with_suffix('.cube'))
    evidence = {'schema': 1, 'backend': spec['tasks'][0]['backend'],
        'run_dir': str(graph.resolve()), 'settings': spec['intent'],
        'fixed_geometry': True, 'geometry_source': geometry_evidence,
        'cube_sha256': file_sha256(cube)}
    atomic_json(destination.with_suffix('.orbitals.json'), evidence)
    return str(destination.with_suffix('.cube'))


def pivot_backend(parameters):
    backend = parameters.get('vrc_tst_pivot_backend', 'gaussian')
    if backend not in ('gaussian', 'molpro'):
        raise ValueError('vrc_tst_pivot_backend must be gaussian or molpro.')
    return backend


def molpro_pivot_settings(parameters):
    """Require native functional names; do not equate two B3LYP conventions."""
    method = parameters.get('vrc_tst_pivot_method', '')
    match = re.fullmatch(r'(UKS|RKS)/([A-Za-z][A-Za-z0-9_-]*)', method, re.I)
    if match is None:
        raise ValueError('Molpro vrc_tst_pivot_method must be UKS/<functional> '
                         'or RKS/<functional>, e.g. UKS/B3LYP3. Molpro B3LYP3 '
                         'uses the VWN3 convention used by Gaussian B3LYP.')
    basis = parameters.get('vrc_tst_pivot_basis', '')
    if not isinstance(basis, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_+*().,-]*', basis):
        raise ValueError('vrc_tst_pivot_basis must be a native Molpro basis name.')
    result = {'command': match[1].lower(), 'functional': match[2].lower(),
              'basis': basis,
              'scf': parameters.get('vrc_tst_scf_convergence', 1.e-8)}
    fortran_threshold(result['scf'], 'VRC pivot SCF convergence')
    for name, default in (('fmax', 5.e-4), ('grid_step', 0.2)):
        value = parameters.get(f'vrc_tst_pivot_{name}', default)
        if (isinstance(value, bool) or not isinstance(value, (int, float))
                or not np.isfinite(value) or value <= 0):
            raise ValueError(f'vrc_tst_pivot_{name} must be finite and positive.')
        result[name] = value
    for name, default, upper in (('maxit', 100, 100), ('scf_maxit', 200, None)):
        value = parameters.get(f'vrc_tst_pivot_{name}', default)
        if type(value) is not int or value < 1 or (upper and value > upper):
            raise ValueError(f'Invalid vrc_tst_pivot_{name}: {value}.')
        result[name] = value
    return result


def molpro_pivot_spec(product, parameters, *, optimize=True):
    from kinbot.vrc_tst_scan import vrc_dispatch_resources, _native_state_integer

    settings = molpro_pivot_settings(parameters)
    atoms = Atoms(list(product.atom), positions=np.asarray(product.geom))
    charge = _native_state_integer(product.charge, 'Fragment charge')
    mult = _native_state_integer(product.mult, 'Fragment multiplicity')
    electrons = int(sum(atoms.numbers)) - charge
    if mult < 1 or electrons < mult - 1 or (electrons + mult - 1) % 2:
        raise ValueError('Fragment electron count and multiplicity are incompatible.')
    if settings['command'] == 'rks' and mult != 1:
        raise ValueError('Open-shell Molpro VRC pivots require UKS.')
    homo = (electrons + mult - 1) // 2
    if homo < 1:
        raise ValueError('Fragment has no occupied orbital for a HOMO pivot.')
    threshold = fortran_threshold(settings['scf'], 'VRC pivot SCF convergence')
    scf = (f"{{{settings['command']},{settings['functional']},"
           f"energy={threshold},maxit={settings['scf_maxit']};"
           'wf,{{NELECTRONS}},1,{{SPIN}};orbital,2200.2}\n')
    body = (MOLPRO_HEADER + f"basis={settings['basis']}\n"
            f'gthresh,energy={threshold}\n' + scf)
    outputs = ['pivot.out', 'final.xyz']
    # Atomic fragments need no geometric pivot direction and no optimization.
    needs_cube = (len(atoms) > 1
                  and parameters.get('pp_orient', 'homo') == 'homo')
    if len(atoms) > 1 and optimize:
        # OPTG tests individual Cartesian components; this conservative
        # conversion implies a vector force norm below the requested fmax.
        gradient = settings['fmax'] * Bohr / Hartree / np.sqrt(3.)
        body += (f"{{optg,maxit={settings['maxit']},gradient={gradient:.12e},"
                 f'energy={threshold}}}\n' + scf
                 + 'forces\nput,xyzgrad,final.xyzgrad\n')
        outputs.append('final.xyzgrad')
    if needs_cube:
        selection = ',type=alpha' if settings['command'] == 'uks' else ''
        body += (f'{{cube,homo.cube;orbital,{homo}.1,record=2200.2{selection};'
                 f"step,{settings['grid_step']},{settings['grid_step']},"
                 f"{settings['grid_step']}}}\n")
        outputs.append('homo.cube')
    body += 'put,xyz,final.xyz\n'
    task = {
        'id': 'pivot', 'kind': 'external', 'backend': 'molpro',
        'geometry_from': 'initial', 'geometry_output': 'final.xyz',
        'resources': vrc_dispatch_resources(parameters),
        'input_name': 'pivot.inp', 'input_template': body,
        'command': ['molpro', '-n', '{cores}', '-m', '{molpro_stack_mw}', '{input}'],
        'stdout': 'launcher.stdout', 'stderr': 'launcher.stderr',
        'required_outputs': outputs,
        'success_marker': {'file': 'pivot.out', 'contains': 'Molpro calculation terminated'},
    }
    if needs_cube:
        task['output_aliases'] = {'homo.cube': [
            f'homo_orbital_{homo}.1.cube', f'homo.cube_orbital_{homo}.1.cube']}
    return {
        'schema': 1, 'name': 'kinbot-vrc-molpro-pivot',
        'molecule': {'symbols': atoms.get_chemical_symbols(),
                     'positions': atoms.positions.tolist(), 'charge': charge,
                     'multiplicity': mult},
        'limits': {'max_nodes': 1}, 'tasks': [task],
        'intent': {'cap_walltime_to_explicit_partition': True,
                   'pivot_settings': settings, 'homo_index': homo,
                   'fixed_geometry': not optimize, 'orbital_only': not optimize,
                   'orbital_type': 'alpha' if settings['command'] == 'uks' else 'canonical'},
    }


def validate_molpro_pivot(directory, spec):
    """Require native completion, a stationary geometry and a matching cube."""
    from kinbot.ase_modules.calculators.molpro import parse_xyzgrad
    directory = Path(directory)
    output = (directory/'pivot.out').read_text(errors='replace')
    if 'Molpro calculation terminated' not in output:
        raise ValueError('Molpro pivot calculation did not terminate normally.')
    validate_molpro_convergence(output)
    if re.search(r'(?i)no convergence|not converged|GLOBAL ERROR|ERROR EXIT', output):
        raise ValueError('Molpro pivot output contains a native solver failure.')
    atoms = read(directory/'final.xyz')
    if (atoms.get_chemical_symbols() != spec['molecule']['symbols']
            or not np.isfinite(atoms.positions).all()):
        raise ValueError('Molpro pivot geometry changed atom identity or is nonfinite.')
    if len(atoms) > 1 and not spec['intent'].get('fixed_geometry'):
        forces = parse_xyzgrad(directory/'final.xyzgrad', atoms)
        fmax = np.linalg.norm(forces, axis=1).max()
        if fmax > spec['intent']['pivot_settings']['fmax']:
            raise ValueError(f'Molpro pivot optimization is not converged: fmax={fmax}.')
    if 'homo.cube' in spec['tasks'][0]['required_outputs']:
        cube = read_orbital_cube(directory/'homo.cube')
        if (cube['atoms'].get_chemical_symbols() != atoms.get_chemical_symbols()
                or not np.allclose(cube['atoms'].positions, atoms.positions,
                                   atol=2.e-5, rtol=0)):
            raise ValueError('Molpro pivot cube does not match the final geometry.')
    return atoms


def prepare_molpro_fragment(product, parameters, destination):
    """Resume one immutable pivot task without requiring Gaussian utilities."""
    from kinbot.anl.dispatch import _load, preflight_until_scheduler_ready
    from kinbot.vrc_tst_scan import _drive_vrc_dispatch

    spec = molpro_pivot_spec(product, parameters)
    destination = Path(destination)
    from kinbot.vrc_graph import select_graph
    run_dir = select_graph(destination.parent/'molpro_pivots', destination.name, spec)
    preflight_until_scheduler_ready(run_dir, interval=20.)
    returncode, summary = _drive_vrc_dispatch(run_dir, ['pivot'])
    if returncode:
        raise RuntimeError(f'Molpro VRC pivot failed at {run_dir}: {summary}')
    _, stored, state = _load(run_dir)
    if state['tasks']['pivot']['status'] != 'complete':
        raise RuntimeError(f'Molpro VRC pivot is incomplete at {run_dir}.')
    task_dir = run_dir/'tasks'/'pivot'
    atoms = validate_molpro_pivot(task_dir, stored)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if (task_dir/'homo.cube').exists():
        temporary = destination.with_suffix('.cube.tmp')
        shutil.copyfile(task_dir/'homo.cube', temporary)
        temporary.replace(destination.with_suffix('.cube'))
    destination.with_suffix('.orbitals.json').write_text(json.dumps({
        'schema': 1, 'backend': 'molpro', 'run_dir': str(run_dir.resolve()),
        'settings': stored['intent'],
        'cube_sha256': (hashlib.sha256((task_dir/'homo.cube').read_bytes()).hexdigest()
                        if (task_dir/'homo.cube').exists() else None),
    }, indent=2) + '\n')
    return atoms
