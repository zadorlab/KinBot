"""Pinned small-species targets for chemistry-level ANL validation.

These values are regression evidence, not replacement calculation results.
Only a task that completed normally and passed KinBot's native-output and
provenance checks is compared.  A published number is never inserted into an
assembled composite result.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from kinbot.anl.dispatch import _load
from kinbot.anl.workflow import _verified_task_result


SOURCE = {
    'doi': '10.1021/acs.jpca.7b05945',
    'supplement_doi': '10.1021/acs.jpca.7b05945.s001',
    'geometry_supplement_doi': '10.1021/acs.jpca.7b05945.s002',
    'title': ('Ab Initio Computations and Active Thermochemical Tables Hand '
              'in Hand: Heats of Formation of Core Combustion Species'),
}

# Energies are hartree and formation enthalpies are kcal mol-1.  The targets
# below are transcribed from the source workbook's ``stable,TZ`` and
# ``Stable,QZ`` sheets.  The geometry convention is part of the benchmark
# identity because an absolute-energy comparison across different optimized
# geometries is not meaningful at the tolerances used here.
BENCHMARKS = {
    'ethane-tz-2017': {
        'species': 'C2H6',
        'geometry': 'CCSD(T)/cc-pVTZ optimized geometry',
        'worksheet': 'stable,TZ!A33:CU33',
        'tasks': {
            'harmonic': {'path': ('zpe', 'hartree'),
                         'expected': 0.07472969, 'atol': 2e-6},
            'f12_tz': {'path': ('energy_hartree',),
                       'expected': -79.705837768001, 'atol': 2e-6},
            'f12_qz': {'path': ('energy_hartree',),
                       'expected': -79.709469984796, 'atol': 2e-6},
            'ccsdt_dz': {'path': ('energy_hartree',),
                         'expected': -79.582320541811, 'atol': 2e-6},
            'ccsdtq_dz': {'path': ('energy_hartree',),
                          'expected': -79.583342872557, 'atol': 2e-6},
            'cfour_dboc': {'path': ('selected', 'hartree'),
                           'expected': 0.0047235917, 'atol': 2e-7},
        },
        'derived': {
            'delta_q_dz': {
                'high': 'ccsdtq_dz', 'low': 'ccsdt_dz',
                'expected': -0.001022330746, 'atol': 3e-7},
        },
        'formation_enthalpy_0k_kcal_mol': {
            'ANL0': -16.458410745169683,
            'ANL0-F12': -16.453662494159296,
        },
    },
    'ethane-qz-2017': {
        'species': 'C2H6',
        'geometry': 'CCSD(T)/cc-pVQZ optimized geometry',
        'worksheet': 'Stable,QZ!A31:CF31',
        'tasks': {
            'ccsdt_tz': {'path': ('energy_hartree',),
                         'expected': -79.674438136905, 'atol': 2e-6},
            'ccsdtq_tz': {'path': ('energy_hartree',),
                          'expected': -79.675402085034, 'atol': 2e-6},
            'ccsdtq_dz': {'path': ('energy_hartree',),
                          'expected': -79.583194551738, 'atol': 2e-6},
            'ccsdtqp_dz': {'path': ('energy_hartree',),
                           'expected': -79.583204093642, 'atol': 2e-6},
        },
        'derived': {
            'delta_q_tz': {
                'high': 'ccsdtq_tz', 'low': 'ccsdt_tz',
                'expected': -0.000963948129, 'atol': 3e-7},
            'delta_p_dz': {
                'high': 'ccsdtqp_dz', 'low': 'ccsdtq_dz',
                'expected': -0.000009541904, 'atol': 3e-7},
        },
        'formation_enthalpy_0k_kcal_mol': {
            'ANL1': -16.486399886603884,
        },
    },
    'methane-qz-2017': {
        'species': 'CH4',
        'geometry': 'CCSD(T)/cc-pVQZ optimized geometry',
        'worksheet': 'Stable,QZ!A22:CF22',
        'molecule': {
            'symbols': ['C', 'H', 'H', 'H', 'H'],
            'positions': [
                [0.0, 0.0, -0.0000000012],
                [0.8882878215, 0.0, 0.6281143484],
                [-0.8882878215, 0.0, 0.6281143484],
                [0.0, 0.8882878285, -0.6281143410],
                [0.0, -0.8882878285, -0.6281143410],
            ],
            'charge': 0, 'multiplicity': 1,
        },
        'tasks': {
            'ccsdt_dz': {'path': ('energy_hartree',),
                         'expected': -40.38698142, 'atol': 2e-6,
                         'request': {
                             'backend': 'molpro', 'method': 'CCSD(T)',
                             'basis': 'cc-pVDZ', 'reference': 'RHF'}},
            'ccsdtq_dz': {'path': ('energy_hartree',),
                          'expected': -40.3875184, 'atol': 2e-6,
                          'request': {
                              'backend': 'cfour', 'method': 'CCSDT(Q)',
                              'basis': 'cc-pVDZ', 'reference': 'RHF',
                              'correlation': 'unrestricted',
                              'driver': 'VCC'}},
            'ccsdtqp_dz': {'path': ('energy_hartree',),
                           'expected': -40.387527542311, 'atol': 2e-6,
                           'request': {
                               'backend': 'mrcc', 'method': 'CCSDTQ(P)',
                               'basis': 'cc-pVDZ', 'reference': 'RHF',
                               'correlation': 'unrestricted'}},
        },
        'derived': {
            'delta_q_dz': {
                'high': 'ccsdtq_dz', 'low': 'ccsdt_dz',
                'expected': -0.00053698, 'atol': 3e-7},
            'delta_p_dz': {
                'high': 'ccsdtqp_dz', 'low': 'ccsdtq_dz',
                'expected': -0.000009142311, 'atol': 3e-7},
        },
        'formation_enthalpy_0k_kcal_mol': {
            'ANL0': -15.905831739961759,
            'ANL0-F12': -15.905831739961759,
            'ANL1': -15.905831739961759,
        },
    },
    'methyl-qz-2017': {
        'species': 'CH3',
        'geometry': 'CCSD(T)/cc-pVQZ optimized geometry',
        'worksheet': 'Stable,QZ!A21:CF21',
        'molecule': {
            'symbols': ['C', 'H', 'H', 'H'],
            'positions': [
                [0.0, 0.0, 0.0],
                [1.0777376714, 0.0, 0.0],
                [-0.5388688357, -0.9333482020, 0.0],
                [-0.5388688357, 0.9333482020, 0.0],
            ],
            'charge': 0, 'multiplicity': 2,
        },
        'tasks': {
            'ccsdt_dz': {'path': ('energy_hartree',),
                         'expected': -39.71574341, 'atol': 2e-6,
                         'request': {
                             'backend': 'molpro', 'method': 'CCSD(T)',
                             'basis': 'cc-pVDZ', 'reference': 'ROHF'}},
            'ccsdtq_dz': {'path': ('energy_hartree',),
                          'expected': -39.71626066, 'atol': 2e-6,
                          'request': {
                              'backend': 'mrcc', 'method': 'CCSDT(Q)',
                              'basis': 'cc-pVDZ', 'reference': 'ROHF',
                              'correlation': 'unrestricted'}},
            'ccsdtqp_dz': {'path': ('energy_hartree',),
                           'expected': -39.716271925574, 'atol': 2e-6,
                           'request': {
                               'backend': 'mrcc', 'method': 'CCSDTQ(P)',
                               'basis': 'cc-pVDZ', 'reference': 'ROHF',
                               'correlation': 'unrestricted'}},
        },
        'derived': {
            'delta_q_dz': {
                'high': 'ccsdtq_dz', 'low': 'ccsdt_dz',
                'expected': -0.00051725, 'atol': 3e-7},
            'delta_p_dz': {
                'high': 'ccsdtqp_dz', 'low': 'ccsdtq_dz',
                'expected': -0.000011265574, 'atol': 3e-7},
        },
        'formation_enthalpy_0k_kcal_mol': {
            'ANL0': 35.8358988626,
            'ANL0-F12': 35.86931404,
            'ANL1': 35.79702527592331,
        },
    },
}


def _nested(record, path):
    value = record
    for key in path:
        if not isinstance(value, dict) or key not in value:
            raise ValueError(f'Parsed result lacks {".".join(path)}.')
        value = value[key]
    if isinstance(value, bool) or not isinstance(value, (int, float)) \
            or not math.isfinite(value):
        raise ValueError(f'Parsed result {".".join(path)} is not finite.')
    return float(value)


def compare_values(name, observed):
    """Compare verified task values with one source and geometry convention."""
    if name not in BENCHMARKS:
        raise ValueError(f'Unknown literature benchmark {name!r}.')
    benchmark = BENCHMARKS[name]
    checks = {}
    for task_id, target in benchmark['tasks'].items():
        if task_id not in observed:
            continue
        value = observed[task_id]
        error = value - target['expected']
        checks[task_id] = {
            'observed_hartree': value,
            'expected_hartree': target['expected'],
            'error_hartree': error,
            'absolute_tolerance_hartree': target['atol'],
            'passed': abs(error) <= target['atol'],
        }
    for key, target in benchmark.get('derived', {}).items():
        if target['high'] not in observed or target['low'] not in observed:
            continue
        value = math.fsum((observed[target['high']],
                           -observed[target['low']]))
        error = value - target['expected']
        checks[key] = {
            'observed_hartree': value,
            'expected_hartree': target['expected'],
            'error_hartree': error,
            'absolute_tolerance_hartree': target['atol'],
            'passed': abs(error) <= target['atol'],
            'derived_from': [target['high'], target['low']],
        }
    if not checks:
        raise ValueError('Run and benchmark have no completed components in common.')
    return {
        'schema': 1,
        'status': ('passed' if all(item['passed'] for item in checks.values())
                   else 'failed'),
        'benchmark': name,
        'species': benchmark['species'],
        'geometry_convention': benchmark['geometry'],
        'source': SOURCE,
        'worksheet': benchmark['worksheet'],
        'checks': checks,
        'published_formation_enthalpy_0k_kcal_mol':
            benchmark['formation_enthalpy_0k_kcal_mol'],
        'note': ('Published values are comparison targets only and were not '
                 'inserted into any computed component or recipe.'),
    }


def compare_run(name, run_dir):
    """Verify completed native tasks, then compare their parsed values."""
    run_dir, spec, state = _load(run_dir)
    benchmark = BENCHMARKS.get(name)
    if benchmark is None:
        raise ValueError(f'Unknown literature benchmark {name!r}.')
    tasks = {task['id']: task for task in spec['tasks']}
    observed = {}
    for task_id, target in benchmark['tasks'].items():
        if (task_id not in tasks
                or state['tasks'].get(task_id, {}).get('status') != 'complete'):
            continue
        task = tasks[task_id]
        parser = task.get('result_parser', {})
        actual_request = {'backend': task.get('backend'), **parser}
        mismatch = {
            key: {'expected': value, 'observed': actual_request.get(key)}
            for key, value in target.get('request', {}).items()
            if actual_request.get(key) != value}
        if mismatch:
            raise ValueError(
                f'{task_id}: calculation request does not match the '
                f'literature benchmark: {mismatch}')
        *_, parsed = _verified_task_result(run_dir, task_id)
        observed[task_id] = _nested(parsed, target['path'])
    return compare_values(name, observed)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Compare hash-verified ANL tasks with pinned literature')
    commands = parser.add_subparsers(dest='action', required=True)
    show = commands.add_parser('show')
    show.add_argument('benchmark', choices=sorted(BENCHMARKS))
    compare = commands.add_parser('compare-run')
    compare.add_argument('benchmark', choices=sorted(BENCHMARKS))
    compare.add_argument('run_dir', type=Path)
    prepare_parser = commands.add_parser('prepare-higher-order')
    prepare_parser.add_argument(
        'benchmark', choices=tuple(
            name for name, value in BENCHMARKS.items() if 'molecule' in value))
    prepare_parser.add_argument('run_dir', type=Path)
    prepare_parser.add_argument('--max-nodes', type=int, default=2)
    prepare_parser.add_argument('--partition')
    prepare_parser.add_argument('--mrcc-command', default='dmrcc')
    prepare_parser.add_argument(
        '--task', dest='task_ids', action='append',
        choices=('ccsdt_tz', 'ccsdt_dz', 'ccsdtq_tz', 'ccsdtq_dz',
                 'ccsdtqp_dz'))
    args = parser.parse_args(argv)
    if args.action == 'show':
        print(json.dumps({
            'source': SOURCE,
            'benchmark': args.benchmark,
            **BENCHMARKS[args.benchmark],
        }, indent=2, sort_keys=True))
        return 0
    if args.action == 'prepare-higher-order':
        from kinbot.anl.dispatch import prepare
        from kinbot.anl.validation import higher_order_validation_spec
        benchmark = BENCHMARKS[args.benchmark]
        selected = args.task_ids or (
            'ccsdt_dz', 'ccsdtq_dz', 'ccsdtqp_dz')
        spec = higher_order_validation_spec(
            benchmark['molecule'], max_nodes=args.max_nodes,
            partition=args.partition, mrcc_command=args.mrcc_command,
            task_ids=selected)
        for task in spec['tasks']:
            if task['id'] == 'ccsdtqp_dz':
                # The small source-matched probes should finish on a day
                # partition. A genuine timeout can then exercise resume-mrcc.
                task['resources']['walltime'] = '24:00:00'
                task['resources']['max_cores'] = 4
        spec['intent']['literature_benchmark'] = {
            'name': args.benchmark,
            'source': SOURCE,
            'worksheet': benchmark['worksheet'],
            'geometry_convention': benchmark['geometry'],
        }
        print(prepare(spec, args.run_dir))
        return 0
    result = compare_run(args.benchmark, args.run_dir)
    print(json.dumps(result, indent=2, sort_keys=True))
    return int(result['status'] != 'passed')


if __name__ == '__main__':
    raise SystemExit(main())
