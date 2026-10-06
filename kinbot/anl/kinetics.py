"""Hash-verified MESS rate-table parsing and literature diagnostics."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import re


METHYL_RECOMBINATION_SOURCE = {
    'doi': '10.1021/jp951218o',
    'title': 'Recombination of Methyl Radicals. 2. Global Fits of the Rate Coefficient',
    'fit': '8.78e-11 * exp(-T/723) cm^3 s^-1',
}


@dataclass(frozen=True)
class HighPressureRate:
    reactant: str
    product: str
    temperature_k: float
    rate: float


def _text(source):
    if hasattr(source, 'read'):
        return source.read()
    if isinstance(source, Path):
        return source.read_text()
    if isinstance(source, str):
        return source
    return str(source)


def _number(word):
    try:
        value = float(word.replace('D', 'E').replace('d', 'e'))
    except ValueError as error:
        raise ValueError(f'Invalid MESS rate value {word!r}.') from error
    if not math.isfinite(value) or value < 0.:
        raise ValueError('MESS rate values must be finite and nonnegative.')
    return value


def parse_high_pressure_rates(source) -> tuple[HighPressureRate, ...]:
    """Parse MESS's temperature-species high-pressure rate tables.

    The table format is defined by the official MESS ``mess_driver`` output.
    Missing entries marked ``***`` are omitted. Duplicate keys and malformed
    rows fail closed rather than returning a partial table.
    """
    text = _text(source)
    marker = 'High Pressure Rate Coefficients (Temperature-Species Rate Tables):'
    if marker not in text:
        raise ValueError('MESS output has no high-pressure rate table.')
    section = text.split(marker, 1)[1]
    section = section.split('Capture/Escape Rate Coefficients:', 1)[0]
    lines = section.splitlines()
    reactant_pattern = re.compile(r'^\s*Reactant\s*=\s*(\S+)\s*$')
    points = {}
    index = 0
    while index < len(lines):
        match = reactant_pattern.match(lines[index])
        if match is None:
            index += 1
            continue
        reactant = match.group(1)
        index += 1
        while index < len(lines) and not lines[index].strip():
            index += 1
        if index >= len(lines):
            raise ValueError(f'{reactant}: missing MESS high-pressure table header.')
        header = lines[index].split()
        if not header or header[0] != 'T(K)' or len(header) < 2:
            raise ValueError(f'{reactant}: invalid MESS high-pressure table header.')
        products = header[1:]
        if len(products) != len(set(products)):
            raise ValueError(f'{reactant}: duplicate MESS product columns.')
        index += 1
        rows = 0
        while index < len(lines):
            line = lines[index]
            if not line.strip() or reactant_pattern.match(line):
                break
            words = line.split()
            if len(words) != len(products) + 1:
                raise ValueError(f'{reactant}: malformed MESS rate row {line!r}.')
            temperature = _number(words[0])
            if temperature <= 0.:
                raise ValueError('MESS rate temperatures must be positive.')
            for product, word in zip(products, words[1:]):
                if word == '***':
                    continue
                point = HighPressureRate(
                    reactant, product, temperature, _number(word))
                key = reactant, product, temperature
                if key in points:
                    raise ValueError(f'Duplicate MESS high-pressure rate {key}.')
                points[key] = point
            rows += 1
            index += 1
        if not rows:
            raise ValueError(f'{reactant}: MESS high-pressure table has no rows.')
    if not points:
        raise ValueError('MESS high-pressure rate tables contain no rates.')
    return tuple(points[key] for key in sorted(points))


def rate_series(points, reactant, product):
    series = sorted(
        (point for point in points
         if point.reactant == reactant and point.product == product),
        key=lambda point: point.temperature_k)
    if not series:
        raise ValueError(f'No MESS high-pressure rates for {reactant}->{product}.')
    temperatures = [point.temperature_k for point in series]
    if len(temperatures) != len(set(temperatures)):
        raise ValueError(f'Duplicate MESS temperatures for {reactant}->{product}.')
    return tuple(series)


def methyl_recombination_experimental_fit(temperature_k):
    """1996 global experimental high-pressure fit, in cm^3 s^-1."""
    temperature = float(temperature_k)
    if not math.isfinite(temperature) or temperature <= 0.:
        raise ValueError('Temperature must be finite and positive.')
    return 8.78e-11 * math.exp(-temperature / 723.)


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def audit_methyl_recombination(run_dir, *, job='mess_0000', max_factor=3.):
    """Audit a completed CH3 + CH3 MESS result against a published fit.

    Network metadata determines the separated bimolecular reactant and bound
    product names. The comparison is diagnostic for the complete ROTD/MESS
    setup; the empirical fit is never used to construct an input or rate.
    """
    if not math.isfinite(max_factor) or max_factor < 1.:
        raise ValueError('max_factor must be finite and at least one.')
    root = Path(run_dir).resolve()
    directory = root if root.name == 'me' else root / 'me'
    execution_path = directory / 'mess_execution.json'
    networks_path = directory / 'mess_networks.json'
    output_path = directory / f'{job}.out'
    for path in (execution_path, networks_path, output_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    execution = json.loads(execution_path.read_text())
    if execution.get('status') != 'complete':
        raise ValueError('MESS execution record is not complete.')
    jobs = execution.get('jobs', {})
    if job not in jobs:
        raise ValueError(f'{job}: missing from MESS execution record.')
    expected_hash = jobs[job].get('output_sha256')
    observed_hash = _sha256(output_path)
    if expected_hash != observed_hash:
        raise ValueError(f'{job}: MESS output hash disagrees with execution record.')
    networks = json.loads(networks_path.read_text())
    matches = [entry for entry in networks if entry.get('stem') == job]
    if len(matches) != 1:
        raise ValueError(f'{job}: expected one MESS network record.')
    network = matches[0]
    if len(network.get('wells', ())) != 1 or len(network.get('products', ())) != 1:
        raise ValueError(
            f'{job}: methyl recombination audit requires one well and one '
            'bimolecular product.')
    reactant = network['products'][0]
    product = network['wells'][0]
    series = rate_series(
        parse_high_pressure_rates(output_path), reactant, product)
    comparison = []
    for point in series:
        reference = methyl_recombination_experimental_fit(point.temperature_k)
        factor = (max(point.rate, reference) / min(point.rate, reference)
                  if point.rate > 0. else math.inf)
        comparison.append({
            **asdict(point), 'reference_rate': reference,
            'factor_difference': factor,
            'within_factor': factor <= max_factor,
        })
    monotone = all(
        right.rate <= left.rate
        for left, right in zip(series, series[1:]))
    return {
        'schema': 1,
        'status': ('passed' if monotone and all(
            item['within_factor'] for item in comparison) else 'failed'),
        'job': job,
        'output': str(output_path),
        'output_sha256': observed_hash,
        'reaction': f'{reactant}->{product}',
        'units': 'cm^3 s^-1',
        'maximum_accepted_factor': float(max_factor),
        'monotonically_nonincreasing_with_temperature': monotone,
        'source': METHYL_RECOMBINATION_SOURCE,
        'rates': comparison,
        'note': ('The literature fit is an independent diagnostic target and '
                 'was not used to generate the ROTD_py surface or MESS input.'),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Parse and audit hash-verified MESS rate tables')
    commands = parser.add_subparsers(dest='action', required=True)
    parse = commands.add_parser('parse-high-pressure')
    parse.add_argument('output', type=Path)
    compare = commands.add_parser('compare-methyl-recombination')
    compare.add_argument('run_dir', type=Path)
    compare.add_argument('--job', default='mess_0000')
    compare.add_argument('--max-factor', type=float, default=3.)
    args = parser.parse_args(argv)
    if args.action == 'parse-high-pressure':
        print(json.dumps([asdict(point) for point in
                          parse_high_pressure_rates(args.output)], indent=2))
        return 0
    result = audit_methyl_recombination(
        args.run_dir, job=args.job, max_factor=args.max_factor)
    print(json.dumps(result, indent=2, sort_keys=True))
    return int(result['status'] != 'passed')


if __name__ == '__main__':
    raise SystemExit(main())
