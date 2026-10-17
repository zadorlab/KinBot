"""Native multireference evidence, independent of molecule and energy acceptance.

No orbital/state identity is inferred from energy ordering. Coefficients below
the native printing threshold are absent, not zero; AO coefficients are not
populations. Retain the native output and WFU for quantitative overlap analysis.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile


DEFAULTS = {'orbital_coefficient_threshold': 0.05,
            'ci_coefficient_threshold': 0.05, 'virtual_orbitals': 2}
NUMBER = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?'
PROGRAM = re.compile(r'^\s*PROGRAM\s*\*\s*(\S+)', re.I)
MR_PROGRAM = re.compile(r'^(?:MULTI|MCSCF|MRCI(?:-F12)?|CI(?:-F12)?|RS2\S*)$', re.I)
MULTI = re.compile(r'\{multi\b[^}]*\}', re.I | re.S)


def diagnostic_options(options=None):
    if options is None:
        return dict(DEFAULTS)
    if not isinstance(options, dict) or set(options) - set(DEFAULTS):
        raise ValueError('Unknown orbital_diagnostics option(s).')
    result = {**DEFAULTS, **options}
    for key in ('orbital_coefficient_threshold', 'ci_coefficient_threshold'):
        value = result[key]
        if type(value) not in (int, float) or not math.isfinite(value) or not 0 < value < 1:
            raise ValueError(f'orbital_diagnostics.{key} must be between zero and one.')
    value = result['virtual_orbitals']
    if type(value) is not int or value < 0:
        raise ValueError('orbital_diagnostics.virtual_orbitals must be a nonnegative integer.')
    return result


def native_print_directives(options=None):
    options = diagnostic_options(options)
    return ('! KinBot native orbital diagnostics (printing only)\n'
            'gprint,orbitals,orben,civector,ref\n'
            f'gthresh,thrprint={options["orbital_coefficient_threshold"]:.12g},'
            f'printci={options["ci_coefficient_threshold"]:.12g}\n')


def ensure_native_printing(text):
    """Cover custom Molpro dispatcher decks, preserving explicit print settings.

    Generated MR decks already carry the marker. For custom decks, insert the
    defaults before MULTI (or stand-alone CI), without changing any orbital
    records, energies or numerical solver controls.
    """
    if 'KinBot native orbital diagnostics' in text:
        return text
    pattern = r'(?im)^(\s*)(\{?(?:multi|mrci(?:-f12)?|rs2(?:c|-f12)?)\b)'
    if not re.search(pattern, text):
        return text
    # Respect users' explicit global thresholds. GPRINT only enables evidence.
    directives = native_print_directives()
    for key in ('thrprint', 'printci'):
        if re.search(rf'(?im)^\s*gthresh\s*,[^\n]*\b{key}\s*=', text):
            directives = re.sub(rf'{key}=[^,\n]+,?', '', directives)
    directives = directives.replace(',\n', '\n').replace('gthresh\n', '')
    text = re.sub(pattern, lambda m: directives + m[0], text, count=1)
    def add_orbprint(m):
        block = m[0]
        if re.search(r'\borbprint\b', block, re.I):
            return block
        return block[:-1] + '\norbprint,2;}'
    return MULTI.sub(add_orbprint, text)


def _float(value):
    return float(value.replace('D', 'E').replace('d', 'e'))


def requested_models(text):
    """Read declarations from the input actually executed, not a manifest."""
    models = []
    for match in MULTI.finditer(text):
        block = match[0]
        clean = re.sub(r'!.*', '', block)
        def card(name):
            m = re.search(rf'(?:^|[;\n])\s*{name}\s*,([^;\n}}]+)', clean, re.I)
            return m[1].strip() if m else None
        def integers(name):
            value = card(name)
            if value is None:
                return None
            fields = value.split(',')
            return ([int(v.strip()) for v in fields]
                    if all(re.fullmatch(r'\s*\d+\s*', v) for v in fields) else None)
        occ, closed = integers('occ'), integers('closed')
        if occ is not None and closed is not None:
            length = max(len(occ), len(closed))
            occ += [0] * (length - len(occ))
            closed += [0] * (length - len(closed))
        model = {'input_line': text[:match.start()].count('\n') + 1,
                 'occ': occ, 'closed': closed, 'states': [],
                 'start': card('start'), 'natural_record': card('natorb'),
                 'canonical_record': card('canonical'), 'diab': card('diab'),
                 'rotations': re.findall(r'(?im)^\s*rotate\s*,([^;\n}]+)', clean)}
        source = re.search(r'!\s*Source active indices:\s*([^\n]+)', block, re.I)
        model['source_active_indices'] = source[1].strip() if source else None
        if occ is not None and closed is not None and len(occ) == len(closed):
            model['active_orbitals'] = sum(occ) - sum(closed)
            model['destination_active_indices'] = [f'{i}.{sym}'
                for sym, (c, o) in enumerate(zip(closed, occ), 1)
                for i in range(c + 1, o + 1)]
        for m in re.finditer(r'(?:^|[;\n])\s*wf\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)(.*?)(?=(?:[;\n]\s*wf\s*,)|\Z)', clean, re.I | re.S):
            states = re.search(r'(?:[;\n])\s*state\s*,\s*(\d+)', m[4], re.I)
            weights = re.search(r'(?:[;\n])\s*weight\s*,([^;\n}]+)', m[4], re.I)
            model['states'].append({'electrons': int(m[1]), 'symmetry': int(m[2]),
                'multiplicity': int(m[3]) + 1, 'roots': int(states[1]) if states else 1,
                'weights': weights[1].strip() if weights else None,
                'active_electrons': int(m[1]) - 2 * sum(closed) if closed is not None else None})
        models.append(model)
    return models


def _partition(lines, label):
    for line in lines:
        m = re.search(rf'Number of {label}\s+orbitals:\s*(\d+)\s*(?:\(([^)]+)\))?', line, re.I)
        if m:
            return {'total': int(m[1]), 'by_symmetry':
                    [int(v) for v in re.findall(r'\d+', m[2])] if m[2] else None}
    return None


def _orbital_tables(lines, offset, closed, active):
    tables = []
    current, orbital = None, None
    row = re.compile(rf'^\s*(\d+\.\d+)\s+({NUMBER})\s+({NUMBER})(.*)$')
    ao = re.compile(rf'(\d+)\s+(\d+)\s+([A-Za-z][A-Za-z0-9+\-]*)\s+({NUMBER})')
    for i, line in enumerate(lines):
        if re.match(r'^\s*(?:NATURAL|CANONICAL|DIABATIC) ORBITALS(?:\s*\([^)]*\))?\s*$', line, re.I):
            current = {'heading': line.strip(), 'line': offset + i,
                       'orbitals': [], 'format': 'native_labeled_coefficients'}
            tables.append(current)
            orbital = None
        elif current is not None:
            m = row.match(line)
            if m:
                index, sym = map(int, m[1].split('.'))
                c = (closed or {}).get('by_symmetry')
                a = (active or {}).get('by_symmetry')
                role = 'unknown'
                if c is not None and a is not None and sym <= min(len(c), len(a)):
                    role = ('closed' if index <= c[sym-1] else
                            'active' if index <= c[sym-1] + a[sym-1] else 'virtual')
                orbital = {'index': m[1], 'occupation': _float(m[2]),
                           'energy_Eh': _float(m[3]), 'role': role,
                           'line': offset + i, 'ao_coefficients': []}
                current['orbitals'].append(orbital)
                coefficients = m[4]
            elif orbital and re.match(r'^\s+\d+\s+\d+\s+[A-Za-z]', line):
                coefficients = line
            elif (not line.strip() or re.match(r'^\s*(?:=+|Orbital\s+Occupation)', line, re.I)):
                continue
            else:
                current = orbital = None
                continue
            for coeff in ao.finditer(coefficients):
                orbital['ao_coefficients'].append({'center': int(coeff[1]),
                    'function': int(coeff[2]), 'component': coeff[3],
                    'coefficient': _float(coeff[4]), 'line': offset + i})
    return tables


def _ci_tables(lines, offset):
    """Keep native coefficient vectors verbatim with their state headings."""
    result = []
    starts = re.compile(r'^\s*(?:CI Coefficients of|Reference coefficients greater|Coefficients of (?:fixed|rotated) reference functions)', re.I)
    for i, line in enumerate(lines):
        if not starts.search(line):
            continue
        end = i + 1
        began = False
        while end < len(lines):
            value = lines[end].strip()
            if began and not value:
                break
            if value and not re.fullmatch(r'=+', value):
                began = True
            end += 1
        result.append({'line_start': offset + i, 'line_end': offset + end - 1,
                       'text': '\n'.join(lines[i:end])})
    return result


def parse_native_output(text, input_text=''):
    lines = text.splitlines()
    starts = [(i, m[1]) for i, line in enumerate(lines) if (m := PROGRAM.match(line))]
    models = requested_models(input_text)
    report = {'schema': 1, 'requested_casscf': models, 'native_sections': [],
              'issues': [], 'identity_verdict': 'requires_orbital_and_state_overlap_review',
              'limitations': ['AO coefficients are thresholded amplitudes, not atomic populations.',
                  'Matching CAS dimensions/occupations does not prove orbital or state continuity.',
                  'Degenerate roots may rotate; equal energies alone do not establish collapse.',
                  'No calculation is accepted or rejected by this diagnostic report.']}
    report['input_evidence'] = [{'line': i+1, 'text': line.strip()}
        for i, line in enumerate(input_text.splitlines()) if re.search(
            r'(?i)\b(?:basis\s*[=,]|file,|sha-256|thrprint|printci|orbprint|KBVRC|KBGAPS|MYENERGY)', line)]
    report['requested_correlation_calls'] = [m[0] for m in re.finditer(
        r'\{(?:rs2(?:c|-f12)?|mrci(?:-f12)?)\b[^}]*\}', input_text, re.I | re.S)]
    report['native_basis_evidence'] = [{'line': i+1, 'text': line.strip()}
        for i, line in enumerate(lines) if re.search(
            r'(?i)\b(?:basis .*loaded|NUMBER OF (?:CONTRACTIONS|PRIMITIVE AOS|SYMMETRY AOS)|Library entry)', line)]
    report['native_geometries'] = []
    for i, line in enumerate(lines):
        if re.match(r'^\s*ATOMIC COORDINATES\s*$', line):
            end = i + 1
            began = False
            while end < len(lines):
                if began and not lines[end].strip():
                    break
                if re.match(r'^\s*\d+\s+\S+\s+', lines[end]):
                    began = True
                end += 1
            report['native_geometries'].append({'line_start': i+1, 'line_end': end,
                'text': '\n'.join(lines[i:end]), 'units': 'bohr (Molpro ATOMIC COORDINATES)'})
    state_number = 0
    for j, (start, program) in enumerate(starts):
        if not MR_PROGRAM.fullmatch(program):
            continue
        end = starts[j+1][0] if j+1 < len(starts) else len(lines)
        block = lines[start:end]
        section = {'program': program, 'line_start': start + 1, 'line_end': end,
                   'evidence': [], 'state_energies': [], 'states': [],
                   'ci_tables': _ci_tables(block, start + 1)}
        evidence = re.compile(r'(?:Number of (?:active|closed|frozen|external|optimized|reference|states|CSFs)|'
            r'Weight factors|Reference symmetry|orbital.*(?:record|dump|overlap)|'
            r'<\s*\d+\.\d+\s*\||ROTATION OF ORBITALS|CONVERGENCE REACHED|'
            r'IPEA shift|Level shift|Cluster corrected energies|Coefficient of reference function|'
            r'RECON[RC]\d+\s*=|KBVRCSEL\s*=|KB_GEOM_ENERGY\s*=|'
            r'\?\s*(?:Error|.*converged|.*problem)|NO CONVERGENCE|ERROR EXIT|GLOBAL ERROR)', re.I)
        for i, line in enumerate(block):
            if evidence.search(line):
                section['evidence'].append({'line': start + i + 1, 'text': line.strip()})
            m = re.search(rf'^\s*!(\S+) STATE\s+(\d+\.\d+)\s*(.*?)\s+Energy\s+({NUMBER})', line, re.I)
            if m:
                section['state_energies'].append({'method': m[1], 'state': m[2],
                    'spin_label': m[3], 'energy_Eh': _float(m[4]), 'line': start + i + 1})
            m = re.search(r'Number of active electrons:\s*(\d+)\s+Spin symmetry=(\S+)\s+Space symmetry=(\d+)', line, re.I)
            if m:
                section['states'].append({'active_electrons': int(m[1]), 'spin': m[2],
                                         'symmetry': int(m[3]), 'line': start + i + 1})
            m = re.search(r'Number of states:\s*(\d+)', line, re.I)
            if m and section['states']:
                section['states'][-1]['roots'] = int(m[1])
        if program.upper() in ('MULTI', 'MCSCF'):
            section['closed'] = _partition(block, 'closed-shell')
            section['active'] = _partition(block, 'active')
            section['orbital_tables'] = _orbital_tables(block, start + 1, section['closed'], section['active'])
            model = models[0] if len(models) == 1 else (models[state_number] if state_number < len(models) else None)
            state_number += 1
            section['requested_model'] = 0 if len(models) == 1 else (state_number - 1 if model else None)
            issues = []
            if not section['active'] or not section['closed'] or not section['states']:
                issues.append('Native CAS dimensions/state declarations missing.')
            natural = [t for t in section['orbital_tables'] if 'NATURAL' in t['heading'].upper()]
            if not natural or not any(t['orbitals'] for t in natural):
                issues.append('Native natural-orbital occupation table missing or unsupported format.')
            for table in natural:
                active_rows = [r for r in table['orbitals'] if r['role'] == 'active']
                if section['active'] and len(active_rows) != section['active']['total']:
                    issues.append('Natural-orbital table lacks some native active indices.')
                if any(not r['ao_coefficients'] for r in active_rows):
                    issues.append('Some active AO coefficients are absent at the print threshold or format is unsupported.')
                if active_rows and section['states']:
                    table['active_occupation_sum'] = sum(r['occupation'] for r in active_rows)
                    if abs(table['active_occupation_sum'] - section['states'][0]['active_electrons']) > max(1e-4, len(active_rows)*1e-5):
                        issues.append('Printed active occupations do not sum to the native active-electron count.')
            if model and section['active'] and model.get('active_orbitals') is not None:
                if model['active_orbitals'] != section['active']['total']:
                    issues.append('Requested and native active-orbital counts disagree.')
                def same_partition(left, right):
                    if right is None:
                        return False
                    size = max(len(left), len(right))
                    return left + [0]*(size-len(left)) == right + [0]*(size-len(right))
                if section['closed'] and not same_partition(model['closed'], section['closed']['by_symmetry']):
                    issues.append('Requested and native closed-orbital partitions disagree.')
                if not same_partition([o-c for o, c in zip(model['occ'], model['closed'])], section['active']['by_symmetry']):
                    issues.append('Requested and native active-orbital symmetry partitions disagree.')
                expected_e = {s['active_electrons'] for s in model['states']}
                if expected_e and any(s['active_electrons'] not in expected_e for s in section['states']):
                    issues.append('Requested and native active-electron counts disagree.')
                if model['states'] and [s['roots'] for s in model['states']] != [s.get('roots') for s in section['states']]:
                    issues.append('Requested and native state counts disagree.')
                spins = {'singlet': 1, 'doublet': 2, 'triplet': 3, 'quartet': 4,
                         'quintet': 5, 'sextet': 6, 'septet': 7, 'octet': 8, 'nonet': 9}
                observed = [spins.get(s['spin'].casefold()) for s in section['states']]
                if model['states'] and None not in observed and observed != [s['multiplicity'] for s in model['states']]:
                    issues.append('Requested and native spin manifolds disagree.')
                if not model['states']:
                    issues.append('Input WF/root declarations are inherited or symbolic; inspect native states.')
            if not section['state_energies']:
                issues.append('Native CASSCF state energies missing.')
            elif len({(s['state'], s['spin_label']) for s in section['state_energies']}) != sum(
                    s.get('roots', 0) for s in section['states']):
                issues.append('Native CASSCF energy list does not cover every declared state.')
            if not section['ci_tables']:
                issues.append('Native CASSCF CI-vector printout missing.')
            section['issues'] = issues
            report['issues'].extend(f'Line {start+1}: {issue}' for issue in issues)
        report['native_sections'].append(section)
    if not report['native_sections']:
        report['issues'].append('No executed multireference program section found.')
    elif state_number == 0:
        report['issues'].append('No native CASSCF orbital table in this output; inspect the source checkpoint calculation.')
    report['normal_termination_seen'] = bool(re.search(r'Molpro calculation terminated', text, re.I))
    report['native_error_lines'] = [{'line': i+1, 'text': line.strip()} for i, line in enumerate(lines)
        if re.search(r'^\s*(?:\?\s*Error|GLOBAL ERROR|ERROR EXIT|NO CONVERGENCE)', line, re.I)]
    if report['native_error_lines'] or not report['normal_termination_seen']:
        report['issues'].append('Native output has an error or no normal-termination marker (may be incomplete).')
    report['diagnostic_status'] = 'review_required' if report['issues'] else 'evidence_recorded'
    return report


def write_native_report(output, input_path=None):
    """Write <stem>.orbital_audit.json for native MR output, including failures.

    Return None for non-MR jobs. Each invocation hashes the exact bytes parsed;
    a live file is a snapshot and must be reread when it changes.
    """
    output = Path(output)
    input_path = Path(input_path) if input_path else output.with_suffix('.inp')
    source = input_path.read_bytes() if input_path.is_file() else b''
    native = output.read_bytes() if output.is_file() else b''
    text, declaration = native.decode(errors='replace'), source.decode(errors='replace')
    if not (MULTI.search(declaration) or re.search(r'\{(?:mrci|rs2)\b', declaration, re.I)
            or any(PROGRAM.match(line) and MR_PROGRAM.fullmatch(PROGRAM.match(line)[1]) for line in text.splitlines())):
        return None
    try:
        report = parse_native_output(text, declaration)
        json.dumps(report, allow_nan=False)
    except (ValueError, TypeError, OverflowError, IndexError) as exc:
        report = {'schema': 1, 'diagnostic_status': 'extraction_error',
                  'native_sections': [], 'issues': [f'{type(exc).__name__}: {exc}'],
                  'identity_verdict': 'unavailable'}
    report['sources'] = {'output': {'path': str(output.resolve()), 'bytes': len(native),
        'sha256': hashlib.sha256(native).hexdigest() if native else None},
        'input': {'path': str(input_path.resolve()), 'sha256': hashlib.sha256(source).hexdigest() if source else None}}
    path = output.with_suffix('.orbital_audit.json')
    # Unique temporary names also permit a read-only monitor and worker to
    # finish a report concurrently without sharing a partially written JSON.
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, prefix=path.name+'.', suffix='.tmp', delete=False) as stream:
        temporary = Path(stream.name)
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return path


def main():
    parser = argparse.ArgumentParser(description='Extract native CAS/MRCI orbital/state evidence; never modifies calculations.')
    parser.add_argument('paths', nargs='+', type=Path, help='Native .out file(s) or run directories (recursive).')
    args = parser.parse_args()
    for entry in args.paths:
        outputs = sorted(entry.rglob('*.out')) if entry.is_dir() else [entry]
        for output in outputs:
            native_input = output.with_suffix('.native.inp')
            report = write_native_report(output, native_input if native_input.exists() else None)
            if report:
                data = json.loads(report.read_text())
                print(f'{data["diagnostic_status"]:18s} {report} ({len(data["native_sections"])} native sections)')


if __name__ == '__main__':
    main()
