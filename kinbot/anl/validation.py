"""Build and audit the portable non-MRCC ANL interface validation graph.

This module deliberately calls the result an *interface validation*.  It runs
the available Gaussian, Molpro, and CFOUR calculation types and verifies their
native outputs, but it cannot label the result ANL1 or ANL1-F12 while the
MRCC-only CCSDTQ(P)/DZ component is disabled.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from ase.db import connect

from kinbot.anl.dispatch import _load, prepare
from kinbot.anl.recipes import recipe
from kinbot.anl.tasks import higher_order_task, molpro_task
from kinbot.anl.workflow import (_verified_task_result,
                                 cbs_task_component, task_component)


MOLPRO_HEADER = """***,KinBot ANL interface validation
symmetry,nosym
orient,noorient
geomtyp=xyz
geometry={
{{XYZ}}
}
set,charge={{CHARGE}}
set,spin={{SPIN}}
"""


def _resources(walltime, *, max_cores=8, min_stack_mw=None, partition=None):
    result = {'cores': 'auto', 'memory_mb': 'node', 'walltime': walltime,
              'max_cores': max_cores}
    if min_stack_mw is not None:
        result['min_stack_mw'] = min_stack_mw
    if partition is not None:
        result['partition'] = partition
    return result


def _molpro_task(ident, body, *, geometry_from='l3_geometry',
                  walltime='04:00:00', max_cores=8, parser=None,
                  partition=None):
    task = {
        'id': ident, 'kind': 'external', 'backend': 'molpro',
        'geometry_from': geometry_from,
        'resources': _resources(walltime, max_cores=max_cores,
                                min_stack_mw=1024, partition=partition),
        'input_name': f'{ident}.inp',
        'input_template': MOLPRO_HEADER + body,
        'command': ['molpro', '-n', '{cores}', '-m',
                    '{molpro_stack_mw}', '{input}'],
        'stdout': 'launcher.stdout', 'stderr': 'launcher.stderr',
        'required_outputs': [f'{ident}.out'],
        'success_marker': {'file': f'{ident}.out',
                           'contains': 'Molpro calculation terminated'},
    }
    if parser:
        task['result_parser'] = {'file': f'{ident}.out', **parser}
    return task


def interface_validation_spec(molecule, *, max_nodes=3, partition=None):
    """Return a molecule-independent Gaussian/Molpro/CFOUR task graph."""
    if not isinstance(molecule, dict):
        raise TypeError('molecule must be an object.')
    l2_profile = {
        'calculator': 'gaussian', 'method': 'B2PLYP',
        'basis': 'cc-pVTZ', 'command': 'g16',
        'calculator_kwargs': {
            'EmpiricalDispersion': 'GD3BJ', 'Symm': 'None',
            'scf': 'xqc', 'integral': 'UltraFine'},
        'optimizer': 'sella',
    }
    reference = 'RHF' if molecule.get('multiplicity', 1) == 1 else 'ROHF'
    conventional = 'uccsd(t),uhf_uccsd=1'
    f12 = 'uccsd(t)-f12b'
    tasks = [
        {
            'id': 'l2_geometry', 'kind': 'ase_optimize',
            'geometry_from': 'initial', 'geometry_output': 'final.xyz',
            'resources': _resources('06:00:00', max_cores=8,
                                    partition=partition),
            'profile': l2_profile,
            'optimizer': {'fmax': 0.0005, 'steps': 160,
                          'sella_kwargs': {'internal': True}},
        },
        {
            'id': 'l3_geometry', 'kind': 'ase_optimize',
            'geometry_from': 'l2_geometry', 'geometry_output': 'final.xyz',
            'resources': _resources('24:00:00', max_cores=12,
                                    min_stack_mw=1024, partition=partition),
            'profile': {
                'calculator': 'molpro', 'method': 'CCSD(T)',
                'basis': 'cc-pVTZ', 'command': 'molpro',
                'optimizer': 'sella'},
            'optimizer': {'fmax': 0.03, 'steps': 100,
                          'sella_kwargs': {'internal': True}},
        },
        _molpro_task(
            'harmonic',
            f'basis=cc-pVTZ\nrhf\n{conventional}\nfrequencies,numerical\n',
            walltime='24:00:00', max_cores=8,
            parser={'kind': 'molpro_harmonic', 'basis': 'cc-pVTZ',
                    'reference': reference},
            partition=partition),
        _molpro_task(
            'f12_tz',
            f'basis=cc-pVTZ-F12\nrhf\n{f12},scale_trip=1\n'
            'kb_f12b=energy\n',
            walltime='12:00:00', max_cores=12,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)-F12b',
                    'basis': 'cc-pVTZ-F12', 'reference': reference},
            partition=partition),
        _molpro_task(
            'f12_qz',
            f'basis=cc-pVQZ-F12\nrhf\n{f12},scale_trip=1\n'
            'kb_f12b=energy\n',
            walltime='24:00:00', max_cores=12,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)-F12b',
                    'basis': 'cc-pVQZ-F12', 'reference': reference},
            partition=partition),
        _molpro_task(
            'ccsdt_dz',
            f'basis=cc-pVDZ\nrhf\n{conventional}\nkb_dz_energy=energy\n',
            walltime='08:00:00', max_cores=8,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)',
                    'basis': 'cc-pVDZ', 'reference': reference},
            partition=partition),
        {
            'id': 'cfour_dboc', 'kind': 'external', 'backend': 'cfour',
            'geometry_from': 'l3_geometry',
            'resources': _resources('08:00:00', max_cores=4,
                                    partition=partition),
            'input_name': 'ZMAT',
            'input_template': (
                'KinBot DBOC interface validation\n{{CARTESIAN}}\n\n'
                '*CFOUR(CALC=SCF\nBASIS=cc-pVTZ\nDBOC=ON\n'
                'COORDINATES=CARTESIAN\nUNITS=ANGSTROM\nCHARGE={{CHARGE}}\n'
                'MULTIPLICITY={{MULT}}\nMEM_UNIT=MB\n'
                'MEMORY_SIZE={{WORK_MEMORY_MB}})\n'),
            'command': ['xcfour'], 'stdout': 'cfour.out',
            'stderr': 'cfour.err', 'required_outputs': ['cfour.out'],
            'files_from_env': {'GENBAS': 'CFOUR_GENBAS'},
            'success_marker': {
                'file': 'cfour.out',
                'contains': 'The total diagonal Born-Oppenheimer correction (DBOC) is:'},
            'result_parser': {'kind': 'cfour_dboc', 'file': 'cfour.out',
                              'level': 'HF', 'basis': 'cc-pVTZ'},
        },
        {
            'id': 'gaussian_vpt2', 'kind': 'external',
            'backend': 'gaussian', 'geometry_from': 'l2_geometry',
            'depends_on': ['l3_geometry'],
            'resources': _resources('24:00:00', max_cores=8,
                                    partition=partition),
            'input_name': 'vpt2.com',
            'input_template': (
                '%nprocshared={{CORES}}\n%mem={{WORK_MEMORY_MB}}MB\n'
                '#p B2PLYP/cc-pVTZ Freq=Anharmonic NoSymm SCF=XQC '
                'EmpiricalDispersion=GD3BJ Integral=UltraFine\n\n'
                'KinBot frequency-only VPT2 interface validation\n\n'
                '{{CHARGE}} {{MULT}}\n{{CARTESIAN}}\n\n'),
            'command': ['g16'], 'stdin': 'vpt2.com',
            'stdout': 'vpt2.log', 'stderr': 'vpt2.err',
            'required_outputs': ['vpt2.log'],
            'success_marker': {'file': 'vpt2.log',
                               'contains': 'Normal termination of Gaussian'},
            'result_parser': {'kind': 'gaussian_vpt2', 'file': 'vpt2.log',
                              'method': 'B2PLYP', 'basis': 'cc-pVTZ',
                              'dispersion': 'GD3BJ'},
        },
    ]
    return {
        'schema': 1, 'name': 'anl1-f12-non-mrcc-interface-validation',
        'molecule': molecule,
        'limits': {'max_nodes': max_nodes},
        'tasks': tasks,
        'intent': {
            'requested_ladder_head': 'ANL1-F12',
            'mrcc_enabled': False,
            'claim': 'interface-validation-only',
            'unavailable_components': [
                'ANL1-F12 pinned electronic equation',
                'CCSDTQ(P)/cc-pVDZ (MRCC)',
                'core-valence CBS provider',
                'scalar-relativistic provider',
                'state-specific spin-orbit provider',
                'complete recipe assembly and CBH/ATcT solve'],
        },
    }


def higher_order_validation_spec(molecule, *, max_nodes=3, partition=None):
    """Build the five-job ANL1 higher-order interface probe.

    This graph checks native input generation, execution, parsing, and the
    cross-program correction arithmetic on one supplied geometry.  It is a
    focused backend test and does not claim a complete composite energy.
    """
    if not isinstance(molecule, dict):
        raise TypeError('molecule must be an object.')
    multiplicity = molecule.get('multiplicity', 1)
    restricted_reference = 'RHF' if multiplicity == 1 else 'ROHF'
    conventional = 'uccsd(t),uhf_uccsd=1'
    common = dict(geometry_from='initial', partition=partition)
    tasks = [
        molpro_task(
            'ccsdt_tz',
            f'basis=cc-pVTZ\nrhf\n{conventional}\nkb_ccsdt=energy\n',
            walltime='08:00:00', max_cores=8,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)',
                    'basis': 'cc-pVTZ',
                    'reference': restricted_reference}, **common),
        molpro_task(
            'ccsdt_dz',
            f'basis=cc-pVDZ\nrhf\n{conventional}\nkb_ccsdt=energy\n',
            walltime='04:00:00', max_cores=8,
            parser={'kind': 'molpro_energy', 'method': 'CCSD(T)',
                    'basis': 'cc-pVDZ',
                    'reference': restricted_reference}, **common),
        higher_order_task(
            'ccsdtq_tz', 'CCSDT(Q)', 'cc-pVTZ',
            multiplicity=multiplicity, walltime='24:00:00',
            max_cores=8, **common),
        higher_order_task(
            'ccsdtq_dz', 'CCSDT(Q)', 'cc-pVDZ',
            multiplicity=multiplicity, walltime='12:00:00',
            max_cores=8, **common),
        higher_order_task(
            'ccsdtqp_dz', 'CCSDTQ(P)', 'cc-pVDZ',
            multiplicity=multiplicity, walltime='24:00:00',
            max_cores=4, **common),
    ]
    return {
        'schema': 1, 'name': 'anl1-higher-order-interface-validation',
        'molecule': molecule, 'limits': {'max_nodes': max_nodes},
        'tasks': tasks,
        'intent': {
            'claim': 'higher-order-interface-validation-only',
            'equation': ('CCSDT(Q)/TZ - CCSD(T)/TZ + CCSDTQ(P)/DZ '
                         '- CCSDT(Q)/DZ'),
            'backend_policy': {
                'closed_shell_ccsdt_q': 'direct-mrcc-rhf-unrestricted-cc',
                'open_shell_ccsdt_q': 'direct-mrcc-semicanonical-rohf',
                'closed_shell_ccsdtq_p': 'direct-mrcc-rhf-unrestricted-cc',
                'open_shell_ccsdtq_p': 'direct-mrcc-semicanonical-rohf',
            },
        },
    }


def audit_higher_order_run(run_dir):
    """Reparse and combine a completed higher-order interface probe."""
    _, spec, state = _load(run_dir)
    if spec.get('name') != 'anl1-higher-order-interface-validation':
        raise ValueError('Run is not an ANL1 higher-order validation graph.')
    statuses = {task['id']: state['tasks'].get(task['id'], {}).get(
        'status', 'waiting') for task in spec['tasks']}
    if any(value != 'complete' for value in statuses.values()):
        raise RuntimeError(f'Higher-order run is incomplete: {statuses}')
    state_id = 'higher-order-validation-state'
    components = {
        ident: task_component(run_dir, ident, key=ident, state_id=state_id)
        for ident in ('ccsdt_tz', 'ccsdt_dz', 'ccsdtq_tz',
                      'ccsdtq_dz', 'ccsdtqp_dz')
    }
    correction = math.fsum((
        components['ccsdtq_tz'].value_hartree,
        -components['ccsdt_tz'].value_hartree,
        components['ccsdtqp_dz'].value_hartree,
        -components['ccsdtq_dz'].value_hartree,
    ))
    return {
        'status': 'higher_order_interface_complete',
        'task_statuses': statuses,
        'correction_hartree': correction,
        'components': {
            key: {
                'energy_hartree': value.value_hartree,
                'method': value.method, 'basis': value.basis,
                'backend': value.backend,
                'reference': value.settings.get('reference'),
                'correlation': value.settings.get('correlation'),
                'driver': value.settings.get('driver'),
                'program': value.settings.get('program'),
                'source_sha256': value.source_sha256,
            } for key, value in components.items()
        },
        'claim': 'interface-validation-only',
    }


def molecule_from_database(database, job, *, charge, multiplicity):
    rows = list(connect(str(database)).select(name=job))
    if not rows:
        raise ValueError(f'No KinBot database row named {job!r}.')
    row = rows[-1]
    if row.data.get('status') != 'normal':
        raise ValueError(f'{job}: latest KinBot record is not normal.')
    missing = [key for key in ('energy', 'frequencies', 'zpe')
               if row.data.get(key) is None]
    if missing:
        raise ValueError(f'{job}: incomplete accepted L2 record: {missing}.')
    return {'symbols': list(row.symbols), 'positions': row.positions.tolist(),
            'charge': charge, 'multiplicity': multiplicity}


def audit_interface_run(run_dir):
    """Reparse every declared native result and report the honest boundary."""
    _, spec, state = _load(run_dir)
    statuses = {task['id']: state['tasks'].get(task['id'], {}).get(
        'status', 'waiting') for task in spec['tasks']}
    if any(value != 'complete' for value in statuses.values()):
        raise RuntimeError(f'Interface run is incomplete: {statuses}')
    parsed = {}
    for task in spec['tasks']:
        if task.get('result_parser'):
            *_, result = _verified_task_result(run_dir, task['id'])
            parsed[task['id']] = result
    equation = recipe('ANL0-F12', vpt2_method='B2PLYP-D3BJ')
    requirement = next(item for item in equation.requirements
                       if item.key == 'reference_cbs')
    reference = cbs_task_component(
        run_dir, 'f12_tz', 'f12_qz', requirement=requirement,
        state_id='validation-state', lower_basis='cc-pVTZ-F12',
        upper_basis='cc-pVQZ-F12')
    return {
        'status': 'interface_complete_recipe_incomplete',
        'requested_ladder_head': 'ANL1-F12',
        'mrcc_enabled': False,
        'task_statuses': statuses,
        'parsed_kinds': {key: value['kind'] for key, value in parsed.items()},
        'verified_f12_cbs_hartree': reference.value_hartree,
        'unavailable_components': spec['intent']['unavailable_components'],
    }


def audit_kinbot_run(run_dir, reaction, *, parent=None, hir_points=0,
                     require_rotdpy=False,
                     require_rotdpy_execution=False):
    """Gate downstream ANL work on an accepted KinBot reaction result."""
    run_dir = Path(run_dir).resolve()
    monitor = run_dir / 'kinbot_monitor.out'
    if not monitor.is_file():
        raise RuntimeError('KinBot did not write kinbot_monitor.out.')
    matches = []
    for line in monitor.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[2] == reaction:
            matches.append(fields)
    if len(matches) != 1 or matches[0][0] != '-1':
        raise RuntimeError(f'{reaction}: expected one accepted channel in '
                           f'kinbot_monitor.out, found {matches!r}.')
    products = matches[0][3:]
    if len(products) != 2:
        raise RuntimeError(f'{reaction}: expected two product entries, found '
                           f'{products!r}.')

    normal_hir = []
    if hir_points:
        database = run_dir / 'kinbot.db'
        if not database.is_file():
            raise RuntimeError('KinBot database is missing.')
        prefix = f'hir/{parent}_hir_' if parent else 'hir/'
        for row in connect(str(database)).select():
            if (getattr(row, 'name', '').startswith(prefix)
                    and row.data.get('status') == 'normal'):
                normal_hir.append(row.name)
        if len(set(normal_hir)) < hir_points:
            raise RuntimeError(f'Expected at least {hir_points} accepted '
                               f'hindered-rotor points for {parent}, found '
                               f'{len(set(normal_hir))}.')
        log = run_dir / 'kinbot.log'
        if (log.is_file()
                and 'will be treated as harmonic oscillators' in log.read_text()):
            raise RuntimeError('KinBot demoted a requested hindered rotor to '
                               'a harmonic oscillator.')

    rotdpy = None
    correction = None
    if require_rotdpy_execution:
        require_rotdpy = True
    if require_rotdpy:
        correction = run_dir / 'vrctst' / f'corr_{reaction}.json'
        rotdpy = run_dir / 'rotdPy' / f'{reaction}.py'
        if not correction.is_file() or not rotdpy.is_file():
            raise RuntimeError(f'{reaction}: VRC correction or rotdPy input '
                               'is missing.')
        payload = json.loads(correction.read_text())
        required = {'dist', 'e_samp', 'e_high', 'scan_ref', 'ra',
                    'e_inf_samp', 'e_inf_high', 'frags_atom', 'frags_geom',
                    'frags_mult'}
        if required - payload.keys():
            raise RuntimeError(f'{reaction}: incomplete VRC correction record.')
        if (len(payload['dist']) != len(payload['e_samp'])
                or len(payload['dist']) != len(payload['e_high'])):
            raise RuntimeError(f'{reaction}: inconsistent VRC correction arrays.')
        if not rotdpy.read_text().strip():
            raise RuntimeError(f'{reaction}: rotdPy input is empty.')
        try:
            compile(rotdpy.read_text(), str(rotdpy), 'exec')
        except SyntaxError as error:
            raise RuntimeError(f'{reaction}: rotdPy input is invalid: '
                               f'{error}') from error
        if require_rotdpy_execution:
            execution_file = rotdpy.with_name(f'{reaction}.execution.json')
            if not execution_file.is_file():
                raise RuntimeError(f'{reaction}: rotdPy was not executed.')
            execution = json.loads(execution_file.read_text())
            input_hash = hashlib.sha256(rotdpy.read_bytes()).hexdigest()
            if (execution.get('status') != 'complete'
                    or execution.get('returncode') != 0
                    or execution.get('input_sha256') != input_hash):
                raise RuntimeError(f'{reaction}: rotdPy execution is incomplete '
                                   'or does not match its input.')
            from kinbot.rotdpy import read_result
            rotdpy_result = read_result(rotdpy)
        else:
            execution_file = None
            rotdpy_result = None
    else:
        execution_file = None
        rotdpy_result = None

    return {
        'status': 'kinbot_reaction_complete',
        'reaction': reaction,
        'products': products,
        'normal_hir_points': len(set(normal_hir)),
        'vrc_correction': str(correction) if correction else None,
        'rotdpy_input': str(rotdpy) if rotdpy else None,
        'rotdpy_execution': (str(execution_file)
                             if execution_file else None),
        'rotdpy_surfaces': (rotdpy_result['surface_count']
                            if rotdpy_result else 0),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Prepare or audit ANL native-interface validations')
    commands = parser.add_subparsers(dest='action', required=True)
    build = commands.add_parser('from-db')
    build.add_argument('database', type=Path)
    build.add_argument('job')
    build.add_argument('spec', type=Path)
    build.add_argument('--charge', type=int, default=0)
    build.add_argument('--multiplicity', type=int, default=1)
    build.add_argument('--max-nodes', type=int, default=3)
    build.add_argument('--partition')
    stage = commands.add_parser('prepare-from-db')
    stage.add_argument('database', type=Path)
    stage.add_argument('job')
    stage.add_argument('run_dir', type=Path)
    stage.add_argument('--charge', type=int, default=0)
    stage.add_argument('--multiplicity', type=int, default=1)
    stage.add_argument('--max-nodes', type=int, default=3)
    stage.add_argument('--partition')
    audit = commands.add_parser('audit')
    audit.add_argument('run_dir', type=Path)
    higher = commands.add_parser('higher-order-from-db')
    higher.add_argument('database', type=Path)
    higher.add_argument('job')
    higher.add_argument('spec', type=Path)
    higher.add_argument('--charge', type=int, default=0)
    higher.add_argument('--multiplicity', type=int, default=1)
    higher.add_argument('--max-nodes', type=int, default=3)
    higher.add_argument('--partition')
    prepare_higher = commands.add_parser('prepare-higher-order-from-db')
    prepare_higher.add_argument('database', type=Path)
    prepare_higher.add_argument('job')
    prepare_higher.add_argument('run_dir', type=Path)
    prepare_higher.add_argument('--charge', type=int, default=0)
    prepare_higher.add_argument('--multiplicity', type=int, default=1)
    prepare_higher.add_argument('--max-nodes', type=int, default=3)
    prepare_higher.add_argument('--partition')
    audit_higher = commands.add_parser('audit-higher-order')
    audit_higher.add_argument('run_dir', type=Path)
    gate = commands.add_parser('gate-kinbot')
    gate.add_argument('run_dir', type=Path)
    gate.add_argument('reaction')
    gate.add_argument('--parent')
    gate.add_argument('--hir-points', type=int, default=0)
    gate.add_argument('--require-rotdpy', action='store_true')
    gate.add_argument('--require-rotdpy-execution', action='store_true')
    args = parser.parse_args(argv)
    if args.action == 'audit':
        print(json.dumps(audit_interface_run(args.run_dir), indent=2,
                         sort_keys=True))
        return 0
    if args.action == 'audit-higher-order':
        print(json.dumps(audit_higher_order_run(args.run_dir), indent=2,
                         sort_keys=True))
        return 0
    if args.action == 'gate-kinbot':
        print(json.dumps(audit_kinbot_run(
            args.run_dir, args.reaction, parent=args.parent,
            hir_points=args.hir_points,
            require_rotdpy=args.require_rotdpy,
            require_rotdpy_execution=args.require_rotdpy_execution),
            indent=2, sort_keys=True))
        return 0
    molecule = molecule_from_database(
        args.database, args.job, charge=args.charge,
        multiplicity=args.multiplicity)
    is_higher = args.action in ('higher-order-from-db',
                                'prepare-higher-order-from-db')
    spec = (higher_order_validation_spec(
        molecule, max_nodes=args.max_nodes, partition=args.partition)
        if is_higher else interface_validation_spec(
            molecule, max_nodes=args.max_nodes, partition=args.partition))
    if args.action in ('from-db', 'higher-order-from-db'):
        args.spec.write_text(json.dumps(spec, indent=2) + '\n')
        print(args.spec.resolve())
    else:
        print(prepare(spec, args.run_dir))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
