"""Opt-in rate comparison using an explicit reference file and MESS endpoints."""
import argparse
from dataclasses import asdict
import json
import math
from pathlib import Path

from kinbot.anl.kinetics import audit_rate_tables, HighPressureRate, rate_series


def fit(temperature_k, coefficients):
    """Modified Arrhenius expression with an optional exponential T decay."""
    t = float(temperature_k)
    if not math.isfinite(t) or t <= 0:
        raise ValueError('Temperature must be finite and positive.')
    c = coefficients
    value = float(c['A']) * (t / float(c.get('T0', 1))) ** float(c.get('n', 0))
    value *= math.exp(-float(c.get('activation_temperature', 0)) / t)
    if 'temperature_decay' in c:
        value *= math.exp(-t / float(c['temperature_decay']))
    if not math.isfinite(value) or value <= 0:
        raise ValueError('Reference rate must be finite and positive.')
    return value


def compare(run_dir, *, reference, reactant, product, job='mess_0000'):
    if isinstance(reference, (str, Path)):
        reference = json.loads(Path(reference).read_text())
    if reference.get('schema') != 1:
        raise ValueError('Expected a schema-1 rate reference.')
    maximum_factor = float(reference['maximum_factor'])
    if not math.isfinite(maximum_factor) or maximum_factor < 1:
        raise ValueError('maximum_factor must be finite and at least one.')
    native = audit_rate_tables(run_dir, job=job)
    series = rate_series((HighPressureRate(**r) for r in native['rates']), reactant, product)
    primary = reference['primary']
    points, compared = [], []
    for point in series:
        comparisons = {}
        for name, target in reference['fits'].items():
            interval = target['source'].get('temperature_range_k')
            if interval and not interval[0] <= point.temperature_k <= interval[1]:
                continue
            expected = fit(point.temperature_k, target['coefficients'])
            factor = max(point.rate, expected) / min(point.rate, expected) if point.rate > 0 else None
            comparisons[name] = dict(reference_rate=expected, factor_difference=factor,
                                     within_factor=factor is not None and factor <= maximum_factor)
        points.append(dict(asdict(point), comparisons=comparisons))
        if primary in comparisons:
            compared.append(points[-1])
    if len(compared) < 2:
        raise ValueError('Fewer than two rates in the primary reference range.')
    trend = reference.get('trend')
    if trend not in (None, 'nonincreasing', 'nondecreasing'):
        raise ValueError('Unknown reference trend.')
    trend_ok = all((b['rate'] <= a['rate'] if trend == 'nonincreasing' else b['rate'] >= a['rate'])
                   for a, b in zip(compared, compared[1:])) if trend else True
    return dict(schema=1, status='passed' if trend_ok and all(
        r['comparisons'][primary]['within_factor'] for r in compared) else 'failed',
        job=job, execution_job_key=native['execution_job_key'],
        output=native['output'], output_sha256=native['output_sha256'],
        reaction=f'{reactant}->{product}', reference=reference, trend_matches=trend_ok,
        rates=points, note='Opt-in validation only; no rates or production inputs are changed.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('run_dir', type=Path)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--reactant', required=True)
    p.add_argument('--product', required=True)
    p.add_argument('--job', default='mess_0000')
    a = p.parse_args()
    print(json.dumps(compare(a.run_dir, reference=a.reference, reactant=a.reactant,
                             product=a.product, job=a.job), indent=2))


if __name__ == '__main__':
    main()
