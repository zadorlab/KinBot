"""Hash-verified MESS rate-table parsing and execution audits."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import re

from kinbot.execution import file_sha256 as _sha256


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


def audit_rate_tables(run_dir, *, job='mess_0000'):
    """Verify a completed MESS network without assuming its chemistry or fit.

    Passing verifies execution/provenance and parseable native rates. It is
    deliberately not a claim of agreement with experiment or rate accuracy.
    """
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
    execution_job = job
    if execution_job not in jobs and job.startswith('mess_'):
        legacy_job = job.removeprefix('mess_')
        if legacy_job in jobs:
            execution_job = legacy_job
    if execution_job not in jobs:
        raise ValueError(f'{job}: missing from MESS execution record.')
    expected_hash = jobs[execution_job].get('output_sha256')
    observed_hash = _sha256(output_path)
    if expected_hash != observed_hash:
        raise ValueError(f'{job}: MESS output hash disagrees with execution record.')
    networks = json.loads(networks_path.read_text())
    matches = [entry for entry in networks if entry.get('stem') == job]
    if len(matches) != 1:
        raise ValueError(f'{job}: expected one MESS network record.')
    network = matches[0]
    points = parse_high_pressure_rates(output_path)
    return {'schema': 1, 'status': 'passed', 'job': job,
            'execution_job_key': execution_job, 'output': str(output_path),
            'output_sha256': observed_hash, 'network': network,
            'rates': [asdict(point) for point in points],
            'accuracy_status': 'not_assessed',
            'note': 'Native execution and rate-table audit; no empirical fit is applied.'}


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Parse and audit hash-verified MESS rate tables')
    commands = parser.add_subparsers(dest='action', required=True)
    parse = commands.add_parser('parse-high-pressure')
    parse.add_argument('output', type=Path)
    audit = commands.add_parser('audit-rates')
    audit.add_argument('run_dir', type=Path)
    audit.add_argument('--job', default='mess_0000')
    args = parser.parse_args(argv)
    if args.action == 'parse-high-pressure':
        print(json.dumps([asdict(point) for point in
                          parse_high_pressure_rates(args.output)], indent=2))
        return 0
    if args.action == 'audit-rates':
        print(json.dumps(audit_rate_tables(args.run_dir, job=args.job), indent=2))
        return 0


if __name__ == '__main__':
    raise SystemExit(main())
