"""Molpro fragment optimization and native HOMO cubes for VRC pivot points.

These auxiliary orbitals do not supply the VRC interaction energies.  Their
own dispatcher graph retains the native input, output, geometry and cube.
"""

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


def molpro_pivot_spec(product, parameters):
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
    if len(atoms) > 1:
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
    if len(atoms) > 1:
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
    from kinbot.anl.dispatch import _load, prepare, preflight_until_scheduler_ready
    from kinbot.vrc_tst_scan import _drive_vrc_dispatch

    spec = molpro_pivot_spec(product, parameters)
    digest = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:16]
    destination = Path(destination)
    run_dir = destination.parent/'molpro_pivots'/f'{destination.name}_{digest}'
    if not run_dir.exists():
        prepare(spec, run_dir)
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
