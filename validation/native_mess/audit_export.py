"""Validate the experimental native MESS generator export, without changing it.

This audits the supplied energy-grained model. It does not validate the QC,
the species reduction, or the chemical accuracy of the model.
"""
import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path


def audit(path):
    path = Path(path)
    states, diagonal = {}, {}
    edges, collision = defaultdict(float), defaultdict(float)
    sinks = 0
    header = None
    ended = False
    digest = hashlib.sha256()

    def finite_positive(token, allow_zero=False):
        value = float(token)
        if not math.isfinite(value) or value < 0 or (not allow_zero and value == 0):
            raise ValueError(f'Invalid native number: {token}')
        return value

    with path.open('rb') as handle:
        for lineno, raw in enumerate(handle, 1):
            digest.update(raw)
            words = raw.decode().split()
            if not words:
                continue
            if ended:
                raise ValueError(f'Record after END at line {lineno}')
            kind = words[0]
            if kind == 'BEGIN':
                if header is not None or len(words) != 5 or words[1] != '1':
                    raise ValueError('Expected one schema-1 condition')
                header = [finite_positive(x) for x in words[2:]]
                continue
            if header is None:
                raise ValueError('Missing BEGIN')
            if kind == 'STATE':
                if len(words) != 7:
                    raise ValueError('Malformed STATE')
                i = int(words[1])
                if i in states or i < 0 or int(words[3]) < 0 or not math.isfinite(float(words[4])):
                    raise ValueError('Invalid/duplicate state identity or energy')
                states[i] = (words[2], finite_positive(words[5]), finite_positive(words[6]))
            elif kind == 'COLLISION_DIAGONAL':
                i = int(words[1])
                if len(words) != 3 or i in diagonal:
                    raise ValueError('Invalid/duplicate collision diagonal')
                diagonal[i] = finite_positive(words[2], True)
            elif kind == 'EDGE':
                if len(words) not in (5, 6) or words[4] not in ('collision', 'reaction'):
                    raise ValueError('Malformed/unsupported edge')
                i, j = map(int, words[1:3])
                if i == j:
                    raise ValueError('Self edge is not an off-diagonal transition')
                rate = finite_positive(words[3], True)
                edges[i, j] += rate
                if words[4] == 'collision':
                    collision[i] += rate
            elif kind == 'SINK':
                if len(words) not in (4, 5):
                    raise ValueError('Malformed sink')
                if int(words[1]) not in states:
                    raise ValueError('Sink references unknown state')
                finite_positive(words[3], True)
                sinks += 1
            elif kind == 'END':
                if len(words) != 1:
                    raise ValueError('Malformed END')
                ended = True
            else:
                raise ValueError(f'Unknown record {kind}')
    if not ended or not states or set(states) != set(range(len(states))) or set(diagonal) != set(states):
        raise ValueError('Incomplete export or state/diagonal coverage')
    collision_error = 0.
    for i, value in diagonal.items():
        outgoing = collision[i]
        if not math.isfinite(outgoing):
            raise ValueError('Nonfinite summed collision rate')
        scale = max(outgoing, value)
        collision_error = max(collision_error, abs(outgoing-value)/scale if scale else 0.)
    balance_error = 0.
    for (i, j), forward in edges.items():
        if i not in states or j not in states or not math.isfinite(forward):
            raise ValueError('Invalid edge identity or summed rate')
        reverse = edges.get((j, i), 0.)
        if not forward and not reverse:
            continue
        if not forward or not reverse:
            balance_error = 1.
            continue
        # Compare fluxes in log space, avoiding overflow of rate*b**2.
        delta = math.log(forward) + 2*math.log(states[i][1]) - math.log(reverse) - 2*math.log(states[j][1])
        balance_error = max(balance_error, -math.expm1(-abs(delta)))
    return dict(sha256=digest.hexdigest(), states=len(states), edges=len(edges),
                sink_records=sinks, temperature_K=header[0], pressure_atm=header[1],
                grain_kcal_mol=header[2], all_rates_finite_nonnegative=True,
                max_collision_diagonal_relative_error=collision_error,
                max_detailed_balance_relative_error=balance_error,
                scope='Native unreduced generator only; no output modification or rate certification.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('export', type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.export), indent=2))
