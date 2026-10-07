"""Pinned small-species targets for chemistry-level ANL validation.

These values are regression evidence, not replacement calculation results.
Only a task that completed normally and passed KinBot's native-output and
provenance checks is compared.  A published number is never inserted into an
assembled composite result.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
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
        'smiles': 'CC',
        'geometry': 'CCSD(T)/cc-pVTZ optimized geometry',
        'worksheet': 'stable,TZ!A33:CU33',
        'molecule': {
            'symbols': ['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'],
            'positions': [
                [0.0, 0.0, -0.7644921462],
                [0.0, 0.0, 0.7644921462],
                [1.0179776872, 0.0, -1.1592669740],
                [-0.5089888436, -0.8815945376, -1.1592669740],
                [-0.5089888436, 0.8815945376, -1.1592669740],
                [-1.0179776872, 0.0, 1.1592669740],
                [0.5089888436, 0.8815945376, 1.1592669740],
                [0.5089888436, -0.8815945376, 1.1592669740],
            ],
            'charge': 0, 'multiplicity': 1,
        },
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
        'smiles': 'CC',
        'geometry': 'CCSD(T)/cc-pVQZ optimized geometry',
        'worksheet': 'Stable,QZ!A31:CF31',
        'molecule': {
            'symbols': ['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'],
            'positions': [
                [0.0, 0.0, -0.7632756988],
                [0.0, 0.0, 0.7632756988],
                [1.0169667994, 0.0, -1.1578258559],
                [-0.5084833997, -0.8807190831, -1.1578258559],
                [-0.5084833997, 0.8807190831, -1.1578258559],
                [-1.0169667994, 0.0, 1.1578258559],
                [0.5084833997, 0.8807190831, 1.1578258559],
                [0.5084833997, -0.8807190831, 1.1578258559],
            ],
            'charge': 0, 'multiplicity': 1,
        },
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
        'smiles': 'C',
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
                              'backend': 'mrcc', 'method': 'CCSDT(Q)',
                              'basis': 'cc-pVDZ', 'reference': 'RHF',
                              'correlation': 'unrestricted',
                              'driver': 'direct'}},
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
        'smiles': '[CH3]',
        'geometry': 'CCSD(T)/cc-pVQZ optimized geometry',
        'worksheet': 'Stable,QZ!A21:CF21',
        'reference_note': (
            'The source RHF-UCCSD(T) component uses the restricted HF '
            'determinant. Its MRCC CCSDT(Q) and CCSDTQ(P) components use a '
            'UHF determinant. Current validation uses a semicanonical ROHF '
            'determinant with unrestricted coupled cluster and compares the '
            'results directly as an explicitly recorded reference variant.'),
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
                         'profiled_atol': 5e-5,
                         'comparison_profile': (
                             'published-MOLPRO-2012.1-to-current-MOLPRO-'
                             'RUCCSD(T)'),
                         'request': {
                             'backend': 'molpro', 'method': 'CCSD(T)',
                             'basis': 'cc-pVDZ', 'reference': 'ROHF'}},
            'ccsdtq_dz': {'path': ('energy_hartree',),
                          'expected': -39.71626066, 'atol': 2e-6,
                          'profiled_atol': 5e-5,
                          'request': {
                              'backend': 'mrcc', 'method': 'CCSDT(Q)',
                              'basis': 'cc-pVDZ', 'reference': 'UHF',
                              'correlation': 'unrestricted'},
                          'comparison_variants': {
                              'reference': ['ROHF']}},
            'ccsdtqp_dz': {'path': ('energy_hartree',),
                           'expected': -39.716271925574, 'atol': 2e-6,
                           'profiled_atol': 5e-5,
                           'request': {
                               'backend': 'mrcc', 'method': 'CCSDTQ(P)',
                               'basis': 'cc-pVDZ', 'reference': 'UHF',
                               'correlation': 'unrestricted'},
                           'comparison_variants': {
                               'reference': ['ROHF']}},
        },
        'derived': {
            'delta_q_dz': {
                'high': 'ccsdtq_dz', 'low': 'ccsdt_dz',
                'expected': -0.00051725, 'atol': 3e-7,
                'profiled_atol': 5e-5},
            'delta_p_dz': {
                'high': 'ccsdtqp_dz', 'low': 'ccsdtq_dz',
                'expected': -0.000011265574, 'atol': 3e-7,
                'profiled_atol': 5e-5},
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


def compare_values(name, observed, *, include_absolute=True,
                   profiled_tasks=()):
    """Compare verified task values with one source and geometry convention."""
    if name not in BENCHMARKS:
        raise ValueError(f'Unknown literature benchmark {name!r}.')
    benchmark = BENCHMARKS[name]
    checks = {}
    profiled_tasks = set(profiled_tasks)
    for task_id, target in benchmark['tasks'].items():
        if not include_absolute:
            continue
        if task_id not in observed:
            continue
        value = observed[task_id]
        error = value - target['expected']
        profiled = (task_id in profiled_tasks
                    or bool(target.get('comparison_profile')))
        tolerance = (target.get('profiled_atol', target['atol'])
                     if profiled else target['atol'])
        checks[task_id] = {
            'observed_hartree': value,
            'expected_hartree': target['expected'],
            'error_hartree': error,
            'absolute_tolerance_hartree': tolerance,
            'passed': abs(error) <= tolerance,
        }
        if profiled and tolerance != target['atol']:
            checks[task_id].update({
                'source_exact_tolerance_hartree': target['atol'],
                'tolerance_profile': target.get(
                    'comparison_profile',
                    'published-UHF-to-current-ROHF'),
            })
    for key, target in benchmark.get('derived', {}).items():
        if target['high'] not in observed or target['low'] not in observed:
            continue
        value = math.fsum((observed[target['high']],
                           -observed[target['low']]))
        error = value - target['expected']
        profiled = bool(
            {target['high'], target['low']} & profiled_tasks)
        tolerance = (target.get('profiled_atol', target['atol'])
                     if profiled else target['atol'])
        checks[key] = {
            'observed_hartree': value,
            'expected_hartree': target['expected'],
            'error_hartree': error,
            'absolute_tolerance_hartree': tolerance,
            'passed': abs(error) <= tolerance,
            'derived_from': [target['high'], target['low']],
        }
        if profiled and tolerance != target['atol']:
            checks[key].update({
                'source_exact_tolerance_hartree': target['atol'],
                'tolerance_profile': 'published-UHF-to-current-ROHF',
            })
    if not checks:
        raise ValueError('Run and benchmark have no completed components in common.')
    result = {
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
    if 'reference_note' in benchmark:
        result['reference_note'] = benchmark['reference_note']
    return result


def _distance_signature(molecule):
    """Return an orientation and same-element-permutation invariant geometry."""
    symbols = molecule.get('symbols')
    positions = molecule.get('positions')
    if (not isinstance(symbols, list) or not isinstance(positions, list)
            or len(symbols) != len(positions)):
        raise ValueError('Benchmark molecule has invalid symbols or positions.')
    groups = defaultdict(list)
    for left in range(len(symbols)):
        if len(positions[left]) != 3:
            raise ValueError('Benchmark molecule coordinates must be Cartesian.')
        for right in range(left + 1, len(symbols)):
            key = tuple(sorted((symbols[left], symbols[right])))
            distance = math.sqrt(math.fsum(
                (float(positions[left][axis]) -
                 float(positions[right][axis])) ** 2
                for axis in range(3)))
            groups[key].append(distance)
    return {key: sorted(values) for key, values in groups.items()}


def _geometry_deviation(molecule, reference):
    """Maximum pair-distance difference in angstrom, or infinity on mismatch."""
    if (molecule.get('charge', 0), molecule.get('multiplicity', 1)) != (
            reference.get('charge', 0), reference.get('multiplicity', 1)):
        return math.inf
    actual = _distance_signature(molecule)
    expected = _distance_signature(reference)
    if actual.keys() != expected.keys() or any(
            len(actual[key]) != len(expected[key]) for key in actual):
        return math.inf
    return max((abs(value - expected[key][index])
                for key, values in actual.items()
                for index, value in enumerate(values)), default=0.0)


def compare_run(name, run_dir):
    """Verify completed native tasks, then compare their parsed values."""
    run_dir, spec, state = _load(run_dir)
    benchmark = BENCHMARKS.get(name)
    if benchmark is None:
        raise ValueError(f'Unknown literature benchmark {name!r}.')
    tasks = {task['id']: task for task in spec['tasks']}
    geometry_deviation = _geometry_deviation(
        spec['molecule'], benchmark['molecule'])
    exact_geometry = geometry_deviation <= 5e-5
    observed = {}
    request_variants = {}
    comparison_profiles = {}
    for task_id, target in benchmark['tasks'].items():
        if (task_id not in tasks
                or state['tasks'].get(task_id, {}).get('status') != 'complete'):
            continue
        task = tasks[task_id]
        parser = task.get('result_parser', {})
        actual_request = {'backend': task.get('backend'), **parser}
        allowed_variants = target.get('comparison_variants', {})
        mismatch = {
            key: {'expected': value, 'observed': actual_request.get(key)}
            for key, value in target.get('request', {}).items()
            if (actual_request.get(key) != value
                and actual_request.get(key) not in
                allowed_variants.get(key, ()))}
        if mismatch:
            raise ValueError(
                f'{task_id}: calculation request does not match the '
                f'literature benchmark: {mismatch}')
        variants = {
            key: {'published': value, 'observed': actual_request.get(key)}
            for key, value in target.get('request', {}).items()
            if (actual_request.get(key) != value
                and actual_request.get(key) in
                allowed_variants.get(key, ()))}
        if variants:
            request_variants[task_id] = variants
        if target.get('comparison_profile'):
            comparison_profiles[task_id] = target['comparison_profile']
        *_, parsed = _verified_task_result(run_dir, task_id)
        observed[task_id] = _nested(parsed, target['path'])
    try:
        result = compare_values(
            name, observed, include_absolute=exact_geometry,
            profiled_tasks=set(request_variants) | set(comparison_profiles))
    except ValueError as exc:
        if exact_geometry or 'no completed components' not in str(exc):
            raise
        result = {
            'schema': 1, 'status': 'incomplete', 'benchmark': name,
            'species': benchmark['species'],
            'geometry_convention': benchmark['geometry'], 'source': SOURCE,
            'worksheet': benchmark['worksheet'], 'checks': {},
            'published_formation_enthalpy_0k_kcal_mol':
                benchmark['formation_enthalpy_0k_kcal_mol'],
            'note': ('No geometry-independent literature difference can be '
                     'formed from the completed tasks.'),
        }
        if 'reference_note' in benchmark:
            result['reference_note'] = benchmark['reference_note']
    result['geometry_match'] = exact_geometry
    result['maximum_pair_distance_error_angstrom'] = geometry_deviation
    result['calculation_request_variants'] = request_variants
    result['calculation_profiles'] = comparison_profiles
    result['profiled_variant'] = bool(request_variants or comparison_profiles)
    if not exact_geometry:
        result['absolute_checks_skipped'] = sorted(
            task_id for task_id in observed if task_id in benchmark['tasks'])
        result['note'] += (
            ' Absolute energies were skipped because the run geometry is not '
            'the pinned source geometry; completed same-geometry differences '
            'remain comparable.')
    return result


def compare_formation(name, source, *, method='ANL0-F12',
                      absolute_tolerance_kcal_mol=0.1):
    """Compare one completed CBH/ANL Hf(0 K) with the pinned workbook."""
    if name not in BENCHMARKS:
        raise ValueError(f'Unknown literature benchmark {name!r}.')
    if (not isinstance(absolute_tolerance_kcal_mol, (int, float))
            or isinstance(absolute_tolerance_kcal_mol, bool)
            or not math.isfinite(absolute_tolerance_kcal_mol)
            or absolute_tolerance_kcal_mol < 0.):
        raise ValueError('Formation-enthalpy tolerance must be nonnegative.')
    path = Path(source).resolve()
    try:
        record = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f'{path}: invalid CBH/ANL result.') from error
    formation = record.get('formation') if isinstance(record, dict) else None
    if (not isinstance(record, dict) or record.get('schema') != 1
            or record.get('status') != 'complete'
            or not isinstance(formation, dict)):
        raise ValueError(f'{path}: CBH/ANL result is incomplete.')
    benchmark = BENCHMARKS[name]
    from kinbot.anl.atct import _canonical_smiles
    if _canonical_smiles(formation.get('target_smiles')) != _canonical_smiles(
            benchmark['smiles']):
        raise ValueError('CBH/ANL target does not match the benchmark species.')
    if formation.get('method') != method:
        raise ValueError(
            f'CBH/ANL result uses {formation.get("method")!r}, expected {method!r}.')
    benchmark_method = (method.split(':', 2)[1]
                        if method.startswith('profiled:') else method)
    try:
        expected = float(
            benchmark['formation_enthalpy_0k_kcal_mol'][benchmark_method])
        observed = float(formation['formation_0k_kj_mol']) / 4.184
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError('Formation enthalpy or benchmark method is unavailable.') from error
    if not math.isfinite(observed):
        raise ValueError('CBH/ANL formation enthalpy is nonfinite.')
    error = observed - expected
    return {
        'schema': 1,
        'status': ('passed' if abs(error) <= absolute_tolerance_kcal_mol
                   else 'failed'),
        'benchmark': name, 'species': benchmark['species'],
        'method': method, 'benchmark_method_family': benchmark_method,
        'profiled_variant': method != benchmark_method,
        'observed_formation_enthalpy_0k_kcal_mol': observed,
        'expected_formation_enthalpy_0k_kcal_mol': expected,
        'error_kcal_mol': error,
        'absolute_tolerance_kcal_mol': float(
            absolute_tolerance_kcal_mol),
        'source': SOURCE, 'worksheet': benchmark['worksheet'],
        'computed_result': str(path),
        'note': ('The published value is an independent comparison target '
                 'and was not inserted into the CBH/ANL result.'),
    }


def methyl_uhf_source_validation_spec(*, max_nodes=1, partition=None,
                                      mrcc_command='dmrcc'):
    """Build the one-off published-UHF methyl CCSDT(Q) reproduction.

    Normal ANL generation deliberately remains RHF/ROHF-only.  This pinned
    graph exists solely to reproduce the historical UHF value in the 2017
    source workbook once, at its published geometry and method.
    """
    from kinbot.anl.validation import higher_order_validation_spec

    benchmark = BENCHMARKS['methyl-qz-2017']
    spec = higher_order_validation_spec(
        benchmark['molecule'], max_nodes=max_nodes, partition=partition,
        mrcc_command=mrcc_command, task_ids=('ccsdtq_dz',))
    task = spec['tasks'][0]
    task['input_template'] = task['input_template'].replace(
        'scftype=ROHF\n'
        'rohftype=semicanonical\n'
        'rohfcore=semicanonical\n',
        'scftype=UHF\n')
    task['result_parser']['reference'] = 'UHF'
    task['source_reference_validation'] = 'methyl-qz-2017-UHF'
    spec['intent'].update({
        'claim': 'pinned-literature-source-reference-validation-only',
        'source_reference_validation': 'methyl-qz-2017-UHF',
        'literature_benchmark': {
            'name': 'methyl-qz-2017',
            'source': SOURCE,
            'worksheet': benchmark['worksheet'],
            'geometry_convention': benchmark['geometry'],
        },
    })
    return spec


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Compare hash-verified ANL tasks with pinned literature')
    commands = parser.add_subparsers(dest='action', required=True)
    show = commands.add_parser('show')
    show.add_argument('benchmark', choices=sorted(BENCHMARKS))
    compare = commands.add_parser('compare-run')
    compare.add_argument('benchmark', choices=sorted(BENCHMARKS))
    compare.add_argument('run_dir', type=Path)
    compare_formation_parser = commands.add_parser('compare-formation')
    compare_formation_parser.add_argument(
        'benchmark', choices=sorted(BENCHMARKS))
    compare_formation_parser.add_argument('result', type=Path)
    compare_formation_parser.add_argument('--method', default='ANL0-F12')
    compare_formation_parser.add_argument(
        '--absolute-tolerance-kcal-mol', type=float, default=0.1)
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
    prepare_uhf = commands.add_parser('prepare-methyl-uhf')
    prepare_uhf.add_argument('run_dir', type=Path)
    prepare_uhf.add_argument('--partition')
    prepare_uhf.add_argument('--mrcc-command', default='dmrcc')
    args = parser.parse_args(argv)
    if args.action == 'show':
        print(json.dumps({
            'source': SOURCE,
            'benchmark': args.benchmark,
            **BENCHMARKS[args.benchmark],
        }, indent=2, sort_keys=True))
        return 0
    if args.action == 'compare-formation':
        result = compare_formation(
            args.benchmark, args.result, method=args.method,
            absolute_tolerance_kcal_mol=
                args.absolute_tolerance_kcal_mol)
        print(json.dumps(result, indent=2, sort_keys=True))
        return int(result['status'] != 'passed')
    if args.action == 'prepare-methyl-uhf':
        from kinbot.anl.dispatch import prepare
        spec = methyl_uhf_source_validation_spec(
            partition=args.partition, mrcc_command=args.mrcc_command)
        print(prepare(spec, args.run_dir))
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
