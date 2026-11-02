"""Read native MESS rate matrices without confusing execution with rate validity.

Signs are retained, including tiny negatives. A nonnegative table passes only
this necessary sign/format check, not convergence, detailed balance, or a
physical validation of the reduced kinetic model. No retry is requested here.
"""
import argparse
import json
import math
from pathlib import Path
import re

from kinbot.execution import file_sha256

_CONDITION = re.compile(
    r'^\s*Temperature\s*=\s*(\S+)\s+K'
    r'(?:\s+Pressure\s*=\s*(\S+)\s+(torr|bar|atm))?\s*$')
_SUMMARY = 'High Pressure Rate Coefficients (Temperature-Species Rate Tables):'


def _number(word):
    return float(word.replace('D', 'E').replace('d', 'e'))


def _scope():
    return dict(status='not_present', tables=0, offdiagonal_entries=0,
                negative_count=0, nonfinite_count=0, violations=[], issues=[])


def audit_text(text):
    """Audit species-species matrices, including escape columns and all units.

    Native MESS repeats tables in other layouts; inspect each condition once.
    Diagonals are loss summaries, not independent transitions. Missing ``***``
    entries remain missing, not zero. A pressure reduction may leave entire
    original species rows missing: that does not mean their chemistry vanished.
    """
    result = dict(schema=1, high_pressure=_scope(), pressure_dependent=_scope(),
                  accuracy_status='not_assessed', automatic_retry=False,
                  negative_values_clipped=False)
    lines = text.split(_SUMMARY, 1)[0].splitlines()
    start = next((i+1 for i, line in enumerate(lines)
                  if line.strip() == 'Species-Species Rate Tables:'), None)
    if start is None:
        for scope in (result['high_pressure'], result['pressure_dependent']):
            scope.update(status='not_assessed', issues=['No native species-species rate section.'])
        return result
    starts = [(i, _CONDITION.fullmatch(line)) for i, line in enumerate(lines[start:], start)
              if _CONDITION.fullmatch(line)]
    species = None
    seen = set()
    for block_number, (begin, match) in enumerate(starts):
        end = starts[block_number+1][0] if block_number+1 < len(starts) else len(lines)
        scope = result['pressure_dependent' if match[2] else 'high_pressure']
        scope['tables'] += 1
        context = dict(temperature_token=match[1], pressure_token=match[2],
                       pressure_unit=match[3], line=begin+1)
        try:
            temperature = _number(match[1]); pressure = _number(match[2]) if match[2] else None
            if not math.isfinite(temperature) or temperature <= 0 or (pressure is not None and
                    (not math.isfinite(pressure) or pressure <= 0)):
                raise ValueError('Invalid temperature or pressure.')
            condition = (temperature, pressure, match[3])
            if condition in seen:
                raise ValueError('Duplicate condition table.')
            seen.add(condition)
            header = next((i for i in range(begin+1, end)
                           if lines[i].split()[:1] == ['From\\To']), None)
            if header is None:
                raise ValueError('Missing rate-matrix header.')
            targets = lines[header].split()[1:]
            if not targets or len(set(targets)) != len(targets):
                raise ValueError('Missing or duplicate product columns.')
            if pressure is None:
                if species is not None and species != targets:
                    raise ValueError('Species columns changed between temperatures.')
                species = targets
            if species is None or targets[:len(species)] != species:
                raise ValueError('Pressure matrix lacks its matching high-pressure species columns.')
            for offset, source in enumerate(species, header+1):
                if offset >= end:
                    raise ValueError('Truncated rate matrix.')
                words = lines[offset].split()
                if len(words) != len(targets)+1 or words[0] != source:
                    raise ValueError(f'Missing, reordered, or malformed row for {source}.')
                values = {}
                for target, word in zip(targets, words[1:]):
                    if word == '***':
                        continue
                    value = _number(word)
                    if target != source or not math.isfinite(value):
                        values[target] = (value, word)
                scale = max([value for value, _ in values.values()
                             if math.isfinite(value) and value > 0] or [0.])
                for target, (value, word) in values.items():
                    if target != source:
                        scope['offdiagonal_entries'] += 1
                    if not math.isfinite(value) or value < 0:
                        finite = math.isfinite(value)
                        scope['negative_count' if finite else 'nonfinite_count'] += 1
                        ratio = abs(value)/scale if finite and scale else None
                        if ratio is not None and not math.isfinite(ratio):
                            ratio = None
                        scope['violations'].append(dict(
                            temperature_k=temperature, pressure=pressure,
                            pressure_unit=match[3], source=source, target=target,
                            value=value if finite else None, native_value=word,
                            diagonal=target == source,
                            largest_positive_same_row=scale,
                            negative_to_positive_row_ratio=ratio,
                            line=offset+1))
        except (ValueError, OverflowError) as error:
            scope['issues'].append(dict(**context, reason=str(error)))
    for scope in (result['high_pressure'], result['pressure_dependent']):
        if scope['issues'] or scope['violations']:
            scope['status'] = 'invalid'
        elif scope['offdiagonal_entries']:
            scope['status'] = 'passed_sign_check'
        elif scope['tables']:
            scope['status'] = 'no_rates'
    result['note'] = ('Sign/format audit only. No proof of numerical convergence or '
                      'physical accuracy; no negative rate is replaced by zero. '
                      'High-pressure and finite-pressure scopes are independent.')
    return result


def audit_output(output):
    path = Path(output)
    report = audit_text(path.read_text(errors='replace'))
    return dict(report, output=str(path), output_sha256=file_sha256(path))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args(argv)
    print(json.dumps(audit_output(args.output), indent=2, allow_nan=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
