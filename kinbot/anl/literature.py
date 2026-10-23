"""Compare verified ANL results against explicitly supplied reference data.

Reference values are never used as computed results or implicit defaults.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path

from kinbot.anl.dispatch import _load
from kinbot.anl.workflow import _verified_task_result


def load_reference(reference):
    """Read an explicitly supplied reference catalog; no built-in targets."""
    if isinstance(reference, (str, Path)):
        reference = json.loads(Path(reference).read_text())
    if (not isinstance(reference, dict) or reference.get('schema') != 1
            or not isinstance(reference.get('benchmarks'), dict)
            or not reference['benchmarks']):
        raise ValueError('Reference must be a schema-1 catalog with benchmarks.')
    return reference


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


def compare_values(name, observed, *, reference, include_absolute=True,
                   profiled_tasks=()):
    """Compare verified task values with one source and geometry convention."""
    reference = load_reference(reference)
    benchmarks = reference['benchmarks']
    if name not in benchmarks:
        raise ValueError(f'Unknown literature benchmark {name!r}.')
    benchmark = benchmarks[name]
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
                    'declared_request_variant'),
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
                'tolerance_profile': 'declared_request_variant',
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
        'source': reference.get('source', {}),
        'worksheet': benchmark.get('worksheet'),
        'checks': checks,
        'published_formation_enthalpy_0k_kcal_mol':
            benchmark.get('formation_enthalpy_0k_kcal_mol', {}),
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


def compare_run(name, run_dir, *, reference):
    """Verify completed native tasks, then compare their parsed values."""
    run_dir, spec, state = _load(run_dir)
    reference = load_reference(reference)
    benchmarks = reference['benchmarks']
    benchmark = benchmarks.get(name)
    if benchmark is None:
        raise ValueError(f'Unknown literature benchmark {name!r}.')
    tasks = {task['id']: task for task in spec['tasks']}
    geometry_deviation = _geometry_deviation(
        spec['molecule'], benchmark['molecule'])
    exact_geometry = geometry_deviation <= 5e-5
    observed = {}
    request_variants = {}
    native_request_evidence = {}
    comparison_profiles = {}
    for task_id, target in benchmark['tasks'].items():
        if (task_id not in tasks
                or state['tasks'].get(task_id, {}).get('status') != 'complete'):
            continue
        task = tasks[task_id]
        parser = task.get('result_parser', {})
        *_, parsed = _verified_task_result(run_dir, task_id)
        actual_request = {'backend': task.get('backend'), **parser}
        # Some execution properties are established only by the verified
        # native result.  In particular, MRCC parser declarations identify
        # the program and method but the parser records ``driver=direct``
        # after checking the native output.  Use native evidence only for a
        # benchmark field absent from the staged declaration; an explicit
        # declaration remains authoritative and conflicts still fail.
        for key in target.get('request', {}):
            if actual_request.get(key) is None and parsed.get(key) is not None:
                actual_request[key] = parsed[key]
                native_request_evidence.setdefault(task_id, {})[key] = \
                    parsed[key]
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
        observed[task_id] = _nested(parsed, target['path'])
    try:
        result = compare_values(
            name, observed, reference=reference, include_absolute=exact_geometry,
            profiled_tasks=set(request_variants) | set(comparison_profiles))
    except ValueError as exc:
        if exact_geometry or 'no completed components' not in str(exc):
            raise
        result = {
            'schema': 1, 'status': 'incomplete', 'benchmark': name,
            'species': benchmark['species'],
            'geometry_convention': benchmark['geometry'], 'source': reference.get('source', {}),
            'worksheet': benchmark.get('worksheet'), 'checks': {},
            'published_formation_enthalpy_0k_kcal_mol':
                benchmark.get('formation_enthalpy_0k_kcal_mol', {}),
            'note': ('No geometry-independent literature difference can be '
                     'formed from the completed tasks.'),
        }
        if 'reference_note' in benchmark:
            result['reference_note'] = benchmark['reference_note']
    result['geometry_match'] = exact_geometry
    result['maximum_pair_distance_error_angstrom'] = geometry_deviation
    result['calculation_request_variants'] = request_variants
    result['calculation_request_native_evidence'] = native_request_evidence
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


def compare_formation(name, source, *, reference, method='ANL0-F12',
                      absolute_tolerance_kcal_mol=0.1):
    """Compare one completed CBH/ANL Hf(0 K) with the pinned workbook."""
    reference = load_reference(reference)
    benchmarks = reference['benchmarks']
    if name not in benchmarks:
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
    benchmark = benchmarks[name]
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
        'source': reference.get('source', {}), 'worksheet': benchmark.get('worksheet'),
        'computed_result': str(path),
        'note': ('The published value is an independent comparison target '
                 'and was not inserted into the CBH/ANL result.'),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Compare hash-verified ANL tasks with pinned literature')
    parser.add_argument('--reference', type=Path, required=True)
    commands = parser.add_subparsers(dest='action', required=True)
    show = commands.add_parser('show')
    show.add_argument('benchmark')
    compare = commands.add_parser('compare-run')
    compare.add_argument('benchmark')
    compare.add_argument('run_dir', type=Path)
    compare_formation_parser = commands.add_parser('compare-formation')
    compare_formation_parser.add_argument(
        'benchmark')
    compare_formation_parser.add_argument('result', type=Path)
    compare_formation_parser.add_argument('--method', default='ANL0-F12')
    compare_formation_parser.add_argument(
        '--absolute-tolerance-kcal-mol', type=float, default=0.1)
    prepare_parser = commands.add_parser('prepare-higher-order')
    prepare_parser.add_argument(
        'benchmark')
    prepare_parser.add_argument('run_dir', type=Path)
    prepare_parser.add_argument('--max-nodes', type=int, default=2)
    prepare_parser.add_argument('--partition')
    prepare_parser.add_argument('--mrcc-command', default='dmrcc')
    prepare_parser.add_argument(
        '--task', dest='task_ids', action='append',
        choices=('ccsdt_tz', 'ccsdt_dz', 'ccsdtq_tz', 'ccsdtq_dz',
                 'ccsdtqp_dz'))
    args = parser.parse_args(argv)
    reference = load_reference(args.reference)
    benchmarks = reference['benchmarks']
    if args.action == 'show':
        print(json.dumps({
            'source': reference.get('source', {}),
            'benchmark': args.benchmark,
            **benchmarks[args.benchmark],
        }, indent=2, sort_keys=True))
        return 0
    if args.action == 'compare-formation':
        result = compare_formation(
            args.benchmark, args.result, reference=reference, method=args.method,
            absolute_tolerance_kcal_mol=
                args.absolute_tolerance_kcal_mol)
        print(json.dumps(result, indent=2, sort_keys=True))
        return int(result['status'] != 'passed')
    if args.action == 'prepare-higher-order':
        from kinbot.anl.dispatch import prepare
        from kinbot.anl.validation import higher_order_validation_spec
        benchmark = benchmarks[args.benchmark]
        selected = args.task_ids or (
            'ccsdt_dz', 'ccsdtq_dz', 'ccsdtqp_dz')
        spec = higher_order_validation_spec(
            benchmark['molecule'], max_nodes=args.max_nodes,
            partition=args.partition, mrcc_command=args.mrcc_command,
            task_ids=selected)
        spec['intent']['literature_benchmark'] = {
            'name': args.benchmark,
            'source': reference.get('source', {}),
            'worksheet': benchmark.get('worksheet'),
            'geometry_convention': benchmark['geometry'],
        }
        print(prepare(spec, args.run_dir))
        return 0
    result = compare_run(args.benchmark, args.run_dir, reference=reference)
    print(json.dumps(result, indent=2, sort_keys=True))
    return int(result['status'] != 'passed')


if __name__ == '__main__':
    raise SystemExit(main())
