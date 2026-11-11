"""Read native MESS rate matrices without confusing execution with rate validity.

Signs are retained, including tiny negatives. A nonnegative table passes only
this necessary sign/format check, not convergence, detailed balance, or a
physical validation of the reduced kinetic model. No retry is requested here.
"""
import argparse
import json
import math
from decimal import Decimal
from pathlib import Path
import re

from kinbot.execution import file_sha256

_CONDITION = re.compile(
    r'^\s*Temperature\s*=\s*(\S+)\s+K'
    r'(?:\s+Pressure\s*=\s*(\S+)\s+(torr|bar|atm))?\s*$')
_SUMMARY = 'High Pressure Rate Coefficients (Temperature-Species Rate Tables):'


def _number(word):
    return float(word.replace('D', 'E').replace('d', 'e'))


def _matrix_row(line, header):
    """Read native right-aligned fields that can touch at full precision.

    MESS uses setw without delimiters. A negative ten-digit mantissa with
    a three-digit exponent fills its whole field and touches the previous
    value or stars. Column ends in the header give the native boundaries.
    """
    words = line.split()
    columns = list(re.finditer(r'\S+', header))
    if len(words) == len(columns):
        try:
            for token in words[1:]:
                if token != '***':
                    _number(token)
            return words
        except ValueError:
            pass
    source = re.match(r'\s*(\S+)', line)
    if source is None or len(columns) < 2:
        return words
    fields = [source[1]]
    start = source.end()
    for column in columns[1:]:
        fields.append(line[start:column.end()].strip())
        start = column.end()
    if line[start:].strip():
        return words
    try:
        for token in fields[1:]:
            if token != '***':
                _number(token)
    except ValueError:
        return words
    return fields


def _scope():
    return dict(status='not_present', tables=0, offdiagonal_entries=0,
                negative_count=0, nonfinite_count=0, violations=[], issues=[])


def _printed_interval(token):
    """Rounding interval of a finite native %g-style printed value."""
    value = Decimal(token.replace('D', 'E').replace('d', 'e'))
    # Native %g prints exactly zero as '0'; it does not round small numbers
    # to an integer zero. Do not assign zero an arbitrary +/- 0.5 tolerance.
    radius = Decimal(0) if not value else Decimal(5).scaleb(value.as_tuple().exponent - 1)
    return value, radius


def _check_loss_balance(tokens, source, context, report):
    """Check native loss summaries, with printed rounding kept explicit.

    This is separate from the sign check. Missing entries are never zeroed;
    an incomplete row can only establish an excessive positive lower bound.
    The alternative driver has zero diagonal placeholders and is excluded.
    """
    if source not in tokens or tokens[source] == '***':
        report['not_assessed_rows'] += 1
        return
    finite = {k: v for k, v in tokens.items() if v != '***'}
    if not all(math.isfinite(_number(v)) for v in finite.values()):
        report['not_assessed_rows'] += 1
        return
    loss, loss_error = _printed_interval(finite[source])
    values = [_printed_interval(v) for k, v in finite.items() if k != source]
    outgoing = sum((v for v, _ in values), Decimal(0))
    error = loss_error + sum((e for _, e in values), Decimal(0))
    missing = [k for k, v in tokens.items() if k != source and v == '***']
    complete = not missing
    violation = abs(outgoing-loss) > error if complete else (
        all(v >= 0 for v, _ in values) and outgoing-loss > error)
    if complete:
        report['complete_rows'] += 1
    else:
        report['not_assessed_rows'] += 1
    if violation:
        report['violations'].append(dict(
            **context, source=source, native_loss=finite[source],
            printed_outgoing_sum=str(outgoing),
            printed_rounding_allowance=str(error),
            missing_targets=missing,
            reason=('Loss differs from the sum of outgoing rates beyond printed rounding.'
                    if complete else 'Known positive rates already exceed the loss summary.')))


def audit_text(text):
    """Audit species-species matrices, including escape columns and all units.

    Native MESS repeats tables in other layouts; inspect each condition once.
    Native diagonals are nonnegative loss summaries, not generator diagonals
    with a minus sign. Missing ``***``
    entries remain missing, not zero. A pressure reduction may leave entire
    original species rows missing: that does not mean their chemistry vanished.
    """
    if 'From\\to' in text and 'From\\To' not in text:
        return _audit_v2_text(text)
    result = dict(schema=1, high_pressure=_scope(), pressure_dependent=_scope(),
                  accuracy_status='not_assessed', automatic_retry=False,
                  negative_values_clipped=False)
    loss_balance = dict(complete_rows=0, not_assessed_rows=0, violations=[])
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
                words = _matrix_row(lines[offset], lines[header])
                if len(words) != len(targets)+1 or words[0] != source:
                    raise ValueError(f'Missing, reordered, or malformed row for {source}.')
                if pressure is not None:
                    _check_loss_balance(dict(zip(targets, words[1:])), source,
                                        dict(temperature_k=temperature, pressure=pressure,
                                             pressure_unit=match[3], line=offset+1), loss_balance)
                values = {}
                for target, word in zip(targets, words[1:]):
                    if word == '***':
                        continue
                    value = _number(word)
                    if target != source or not math.isfinite(value) or value < 0:
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
    loss_balance['status'] = ('inconsistent' if loss_balance['violations'] else
                              'consistent_with_printed_rounding' if loss_balance['complete_rows'] else
                              'not_assessed')
    result['pressure_dependent']['loss_balance'] = loss_balance
    return result


def _audit_v2_text(text):
    """MESS's alternative driver prints finite-P tables in column chunks.

    It has no preceding HPL table. Read every chunk rather than calling the
    absent HPL section an error, or silently missing columns after the first.
    Its printed zeros are preserved; they do not establish a convergence bound.
    """
    result = dict(schema=1, format='native_v2', high_pressure=_scope(),
                  pressure_dependent=_scope(), accuracy_status='not_assessed',
                  automatic_retry=False, negative_values_clipped=False)
    scope = result['pressure_dependent']
    lines = text.split('Temperature-Species Rate Tables:', 1)[0].splitlines()
    starts = [(i, _CONDITION.fullmatch(line)) for i, line in enumerate(lines)
              if _CONDITION.fullmatch(line)]
    seen, expected = set(), None
    for block, (begin, match) in enumerate(starts):
        end = starts[block+1][0] if block+1 < len(starts) else len(lines)
        scope['tables'] += 1
        try:
            if match[2] is None:
                raise ValueError('Alternative-driver table lacks a pressure.')
            t, p = _number(match[1]), _number(match[2])
            if not all(math.isfinite(v) and v > 0 for v in (t, p)):
                raise ValueError('Invalid temperature or pressure.')
            condition = (t, p, match[3])
            if condition in seen:
                raise ValueError('Duplicate condition table.')
            seen.add(condition)
            rows, columns, order = {}, [], None
            i = begin + 1
            while i < end:
                words = lines[i].split()
                if not words or words[0] != 'From\\to':
                    i += 1
                    continue
                targets = words[1:]
                if not targets or len(set(targets)) != len(targets) or set(targets) & set(columns):
                    raise ValueError('Missing or repeated column chunk.')
                columns.extend(targets)
                header = lines[i]
                i += 1
                sources = []
                while i < end and lines[i].split():
                    words = _matrix_row(lines[i], header)
                    if len(words) != len(targets)+1 or words[0] in sources:
                        raise ValueError('Malformed or duplicate row in column chunk.')
                    source = words[0]; sources.append(source)
                    row = rows.setdefault(source, {})
                    for target, token in zip(targets, words[1:]):
                        if token != '***':
                            row[target] = (_number(token), token, i+1)
                    i += 1
                if not sources or (order is not None and order != sources):
                    raise ValueError('Incomplete or reordered column chunk.')
                order = sources
            if not columns or [c for c in columns if c in order] != order:
                raise ValueError('Incomplete or reordered species columns.')
            # Native cool.cc omits source rows for Dummy bimolecular products,
            # while retaining their sink columns. Without the input inventory
            # we cannot distinguish these from consistently missing source rows.
            omitted = [c for c in columns if c not in order]
            if omitted:
                scope.setdefault('source_coverage_unverified', []).append(dict(
                    temperature_k=t, pressure=p, pressure_unit=match[3],
                    omitted_sources=omitted))
            if expected is not None and expected != columns:
                raise ValueError('Species changed between conditions.')
            expected = columns
            for source, row in rows.items():
                scale = max([v for target, (v, _, _) in row.items()
                             if target != source and math.isfinite(v) and v > 0] or [0.])
                for target, (v, token, line) in row.items():
                    if source != target:
                        scope['offdiagonal_entries'] += 1
                    if not math.isfinite(v) or v < 0:
                        finite = math.isfinite(v)
                        scope['negative_count' if finite else 'nonfinite_count'] += 1
                        ratio = abs(v)/scale if finite and scale else None
                        scope['violations'].append(dict(
                            temperature_k=t, pressure=p, pressure_unit=match[3],
                            source=source, target=target, value=v if finite else None,
                            native_value=token, diagonal=source == target,
                            largest_positive_same_row=scale, line=line,
                            negative_to_positive_row_ratio=ratio if ratio is None or math.isfinite(ratio) else None))
        except (TypeError, ValueError, OverflowError) as error:
            scope['issues'].append(dict(line=begin+1, reason=str(error)))
    if scope['issues'] or scope['violations']:
        scope['status'] = 'invalid'
    elif scope['offdiagonal_entries']:
        scope['status'] = 'passed_sign_check'
    elif scope['tables']:
        scope['status'] = 'no_rates'
    result['note'] = ('Native alternative-driver finite-pressure tables; sign/format '
                      'audit only. Zeros and missing entries are not convergence '
                      'bounds. Detailed balance, reduction and grid convergence '
                      'require separate validation. No values are clipped.')
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
