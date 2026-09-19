"""Method-aware native QC result extraction for ANL task outputs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re

from ase.units import Hartree, invcm, kJ, mol


_NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DdEe][+-]?\d+)?'
_DBOC_HEADER = re.compile(
    r'^\s*Summary of diagonal Born-Oppenheimer correction at\s+'
    r'(.+?)\s+level\s*$', re.IGNORECASE | re.MULTILINE)
_DBOC_VALUE = re.compile(
    rf'^\s*The total diagonal Born-Oppenheimer correction \(DBOC\) is:\s*'
    rf'({_NUMBER})\s*(a\.u\.|cm-1|kJ/mole)\s*$',
    re.IGNORECASE | re.MULTILINE)
_FINAL_ENERGY = re.compile(
    rf'^\s*The final electronic energy is\s+({_NUMBER})\s+a\.u\.\s*$',
    re.IGNORECASE | re.MULTILINE)
_SCF_WITH_DBOC = re.compile(
    rf'^\s*Total SCF energy including DBOC\s+({_NUMBER})\s*$',
    re.IGNORECASE | re.MULTILINE)
_MOLPRO_CCSD_T = r'^\s*ccsd\(t\)(?:\s*,[^\n]*)?\s*$'
_MOLPRO_F12B = r'^\s*ccsd\(t\)-f12b?\b[^\n]*\bscale_trip\s*=\s*1\b[^\n]*$'


def _number(value):
    number = float(value.replace('D', 'E').replace('d', 'e'))
    if not math.isfinite(number):
        raise ValueError('QC output contains a nonfinite number.')
    return number


def _level(name):
    folded = name.strip().casefold()
    return {'hartree-fock': 'HF', 'hf': 'HF', 'scf': 'HF',
            'mp1': 'MP1'}.get(folded, name.strip().upper())


def parse_cfour_dboc(output, *, level='HF', basis=None):
    """Select one named DBOC level, preserving separately reported levels.

    CFOUR 2.1 prints identical value labels in HF and MP1 summary sections;
    selecting the last label would silently return the wrong ANL component.
    """
    if 'ERROR ERROR ERROR' in output:
        raise ValueError('CFOUR reported an error flag.')
    if 'This computation required' not in output:
        raise ValueError('CFOUR output has no final completion line.')
    if basis is not None:
        if not isinstance(basis, str) or not basis:
            raise ValueError('CFOUR DBOC basis is invalid.')
        echoed = re.findall(r'^\s*BASIS\s*=\s*([^\s,)]+)\s*$', output,
                            re.IGNORECASE | re.MULTILINE)
        if len(echoed) != 1 or echoed[0].casefold() != basis.casefold():
            raise ValueError(f'CFOUR output does not echo basis {basis}.')
    headers = list(_DBOC_HEADER.finditer(output))
    if not headers:
        raise ValueError('CFOUR output has no named DBOC summary.')
    reported = {}
    for index, header in enumerate(headers):
        name = _level(header.group(1))
        if name in reported:
            raise ValueError(f'CFOUR repeated the {name} DBOC summary.')
        end = headers[index + 1].start() if index + 1 < len(headers) else len(output)
        values = {}
        for match in _DBOC_VALUE.finditer(output, header.end(), end):
            unit = {'a.u.': 'hartree', 'cm-1': 'cm_inverse',
                    'kj/mole': 'kj_mol'}[match.group(2).lower()]
            if unit in values:
                raise ValueError(f'CFOUR repeated the {name} DBOC {unit} value.')
            values[unit] = _number(match.group(1))
        if 'hartree' not in values or 'cm_inverse' not in values:
            raise ValueError(f'CFOUR {name} DBOC summary lacks a.u. or cm-1.')
        if values['hartree'] <= 0:
            raise ValueError(f'CFOUR {name} DBOC is not positive.')
        expected_cm = values['hartree'] * Hartree / invcm
        if abs(expected_cm - values['cm_inverse']) > 0.001:
            raise ValueError(f'CFOUR {name} DBOC a.u./cm-1 values disagree.')
        if 'kj_mol' in values:
            expected_kj = values['hartree'] * Hartree * mol / kJ
            if abs(expected_kj - values['kj_mol']) > 0.005:
                raise ValueError(f'CFOUR {name} DBOC a.u./kJ/mole values disagree.')
        reported[name] = values
    requested = _level(level)
    if requested not in reported:
        raise ValueError(f'CFOUR has no {requested} DBOC summary.')
    finals = _FINAL_ENERGY.findall(output)
    if len(finals) != 1:
        raise ValueError('CFOUR needs exactly one final electronic energy.')
    final = _number(finals[0])
    if requested == 'HF':
        scf_with_dboc = _SCF_WITH_DBOC.findall(output)
        if len(scf_with_dboc) != 1 or abs(_number(scf_with_dboc[0]) - final) > 1e-8:
            raise ValueError('CFOUR final energy disagrees with SCF plus DBOC.')
    return {
        'kind': 'cfour_dboc', 'selected_level': requested,
        'selected': reported[requested], 'reported_levels': reported,
        'final_energy_including_dboc_hartree': final,
    }


def _molpro_output(output, basis):
    if not re.search(r'^\s*Molpro calculation terminated\s*$', output,
                     re.IGNORECASE | re.MULTILINE):
        raise ValueError('Molpro output has no normal completion line.')
    if not re.search(rf'^\s*basis\s*=\s*{re.escape(basis)}\s*$', output,
                     re.IGNORECASE | re.MULTILINE):
        raise ValueError(f'Molpro output does not echo basis {basis}.')


def _one_number(output, pattern, label):
    matches = re.findall(pattern, output, re.IGNORECASE | re.MULTILINE)
    if len(matches) != 1:
        raise ValueError(f'Expected exactly one {label}; found {len(matches)}.')
    return _number(matches[0])


def parse_molpro_energy(output, *, method, basis):
    """Read the exact named total energy, rather than a rounded variable or F12a."""
    _molpro_output(output, basis)
    labels = {
        'CCSD(T)': r'CCSD\(T\)',
        'CCSD(T)-F12b': r'CCSD\(T\)-F12b',
    }
    if method not in labels:
        raise ValueError(f'Unsupported Molpro energy method {method!r}.')
    if method == 'CCSD(T)-F12b':
        if not re.search(_MOLPRO_F12B,
                         output, re.IGNORECASE | re.MULTILINE):
            raise ValueError('Molpro output does not echo scaled-triples F12 input.')
    else:
        if not re.search(_MOLPRO_CCSD_T, output,
                         re.IGNORECASE | re.MULTILINE):
            raise ValueError('Molpro output does not echo conventional CCSD(T).')
    energy = _one_number(
        output, rf'^\s*!{labels[method]} total energy\s+({_NUMBER})\s*$',
        f'Molpro {method} total energy')
    if method == 'CCSD(T)-F12b':
        summary = _one_number(
            output, rf'^\s*CCSD\(T\)-F12/{re.escape(basis)} energy\s*=\s*'
                    rf'({_NUMBER})\s*$', 'Molpro F12 summary')
        if abs(summary - energy) > 1e-9:
            raise ValueError('Molpro F12b total energy disagrees with the summary.')
    return {'kind': 'molpro_energy', 'method': method, 'basis': basis,
            'energy_hartree': energy}


def parse_molpro_harmonic(output, *, basis):
    """Read vibrational modes and ZPE, excluding rotations/translations."""
    _molpro_output(output, basis)
    if not re.search(r'^\s*frequencies\s*,\s*numerical\s*$', output,
                     re.IGNORECASE | re.MULTILINE):
        raise ValueError('Molpro output does not echo numerical frequencies.')
    frequency_sections = re.findall(
        r'^\s*PROGRAM \* FREQUENCIES \(Calculation of harmonic '
        r'vibrational spectra for (.+?)\)\s*$', output,
        re.IGNORECASE | re.MULTILINE)
    if len(frequency_sections) != 1 or frequency_sections[0].upper() != 'CCSD(T)':
        raise ValueError('Expected one CCSD(T) Molpro frequency section.')
    section = output.split('PROGRAM * FREQUENCIES', 1)[1]
    low_heading = re.search(r'^\s*Low Vibration\s+Wavenumber\s*$', section,
                            re.IGNORECASE | re.MULTILINE)
    low_modes = []
    if low_heading is not None:
        for line in section[low_heading.end():].splitlines():
            match = re.match(r'^\s*(\d+)\s+(\S+)\s*$', line)
            if match is None:
                if low_modes:
                    break
                continue
            try:
                low_modes.append(_number(match.group(2)))
            except ValueError as exc:
                raise ValueError('Molpro reports a nonreal low mode.') from exc
        if any(value < -1.0 for value in low_modes):
            raise ValueError('Molpro reports an imaginary low mode.')
    heading = re.search(r'^\s*Vibration\s+Wavenumber\s*$', section,
                        re.IGNORECASE | re.MULTILINE)
    if heading is None:
        raise ValueError('Molpro vibrational wavenumber table is absent.')
    modes = []
    for line in section[heading.end():].splitlines():
        match = re.match(r'^\s*(\d+)\s+(\S+)\s*$', line)
        if match is None:
            if modes:
                break
            continue
        if int(match.group(1)) != len(modes) + 1:
            raise ValueError('Molpro vibrational mode numbers are not consecutive.')
        try:
            value = _number(match.group(2))
        except ValueError as exc:
            raise ValueError('Molpro reports a nonreal vibrational mode.') from exc
        if value <= 0:
            raise ValueError('Molpro reports a nonpositive vibrational mode.')
        modes.append(value)
    if not modes:
        raise ValueError('Molpro has no positive vibrational modes.')
    zpe = re.findall(
        rf'^\s*Zero point energy:\s*({_NUMBER})\s*\[H\]\s*'
        rf'({_NUMBER})\s*\[1/CM\]\s*({_NUMBER})\s*\[KJ/MOL\]\s*$',
        output, re.IGNORECASE | re.MULTILINE)
    if len(zpe) != 1:
        raise ValueError('Expected exactly one Molpro harmonic ZPE.')
    hartree, cm_inverse, kj_mol = map(_number, zpe[0])
    if abs(hartree * Hartree / invcm - cm_inverse) > 0.02:
        raise ValueError('Molpro harmonic ZPE units disagree.')
    if abs(0.5 * sum(modes) - cm_inverse) > 0.02:
        raise ValueError('Molpro harmonic modes disagree with ZPE.')
    if abs(hartree * Hartree * mol / kJ - kj_mol) > 0.02:
        raise ValueError('Molpro harmonic ZPE kJ/mol disagrees.')
    return {'kind': 'molpro_harmonic', 'method': 'CCSD(T)', 'basis': basis,
            'wavenumbers_cm_inverse': modes,
            'low_modes_cm_inverse': low_modes,
            'review_required': any(abs(value) > 1.0 for value in low_modes),
            'zpe': {'hartree': hartree, 'cm_inverse': cm_inverse,
                    'kj_mol': kj_mol}}


def _gaussian_dispersion(route):
    match = re.search(r'\bEmpiricalDispersion\s*(?:=|\()\s*([A-Za-z0-9]+)',
                      route, re.IGNORECASE)
    return match.group(1).upper() if match else ''


def _gaussian_fundamental_bands(output):
    """Read Gaussian's mode-resolved harmonic/anharmonic fundamental table.

    Gaussian has used both a numeric-only table and a table with a status
    column.  Resonance flags can also precede the mode number.  Locate the
    ``n(1)`` mode token, skip any intervening text, and take the first two
    numeric fields as E(harm) and E(anharm), respectively.
    """
    marker = output.rfind('Fundamental Bands')
    if marker < 0:
        raise ValueError('Gaussian fundamental-band table is absent.')
    modes = []
    labels = set()
    started = False
    for line in output[marker:].splitlines()[1:]:
        if started and re.match(
                r'^\s*(?:Overtones|Combination Bands|Anharmonic Infrared|Thermochemistry)\b',
                line, re.IGNORECASE):
            break
        tokens = line.split()
        mode_index = None
        mode = None
        label = None
        for index, token in enumerate(tokens):
            match = re.fullmatch(r'(\d+)\(1(?:,[+-]?\d+)?\)', token)
            if match:
                mode_index, mode, label = index, int(match.group(1)), token
                break
        if mode_index is None:
            continue
        numbers = []
        for token in tokens[mode_index + 1:]:
            if re.fullmatch(_NUMBER, token):
                numbers.append(_number(token))
                if len(numbers) == 2:
                    break
        if len(numbers) != 2:
            raise ValueError(f'Gaussian fundamental mode {mode} lacks harmonic '
                             'or anharmonic energy.')
        if label in labels:
            raise ValueError(f'Gaussian repeats fundamental mode {label}.')
        labels.add(label)
        modes.append((mode, label, *numbers))
        started = True
    if not modes:
        raise ValueError('Gaussian fundamental-band table has no modes.')
    base_modes = sorted({row[0] for row in modes})
    if base_modes != list(range(1, max(base_modes) + 1)):
        raise ValueError('Gaussian fundamental mode numbers are not consecutive.')
    harmonic = [row[2] for row in modes]
    anharmonic = [row[3] for row in modes]
    if any(value <= 0. for value in harmonic):
        raise ValueError('Gaussian reports a nonpositive harmonic fundamental.')
    return [row[1] for row in modes], harmonic, anharmonic


def parse_gaussian_vpt2(output, *, method, basis, dispersion=''):
    """Read VPT2 ZPE and mode fundamentals, surfacing native warnings."""
    lines = output.strip().splitlines()
    if not lines or not lines[-1].lstrip().startswith('Normal termination of Gaussian'):
        raise ValueError('Gaussian output has no final normal termination.')
    route = output[:10000]
    if (not re.search(rf'\b{re.escape(method)}/{re.escape(basis)}(?![\w-])', route,
                      re.IGNORECASE)
            or not re.search(r'\bFreq\s*=\s*Anharmonic\b', route,
                             re.IGNORECASE)):
        raise ValueError('Gaussian output does not echo the requested Freq=Anharmonic route.')
    if _gaussian_dispersion(route) != dispersion.upper():
        raise ValueError('Gaussian output dispersion disagrees with the requested level.')
    marker = output.rfind('Anharmonic Zero Point Energy')
    if marker < 0:
        raise ValueError('Gaussian anharmonic ZPE section is absent.')
    section = output[marker:]
    components = {}
    for name, label in (('harmonic', 'Harmonic'),
                        ('anharmonic_potential', r'Anharmonic Pot\.'),
                        ('watson_coriolis', r'Watson\+Coriolis'),
                        ('total_anharmonic', 'Total Anharm')):
        components[name] = _one_number(
            section, rf'^\s*{label}\s*:\s*cm-1\s*=\s*({_NUMBER})\s*;',
            f'Gaussian {name} ZPE')
    expected = (components['harmonic'] + components['anharmonic_potential']
                + components['watson_coriolis'])
    if abs(expected - components['total_anharmonic']) > 0.02:
        raise ValueError('Gaussian anharmonic ZPE components disagree.')
    warnings = [line.strip() for line in output.splitlines()
                if re.match(r'^\s*WARNING:', line, re.IGNORECASE)]
    mode_labels, harmonic_modes, anharmonic_modes = \
        _gaussian_fundamental_bands(output)
    invalid_modes = [index for index, value in enumerate(anharmonic_modes, 1)
                     if value <= 0.]
    if invalid_modes:
        warnings.append('Nonpositive anharmonic fundamental mode(s): ' +
                        ', '.join(map(str, invalid_modes)))
    correction = components['total_anharmonic'] - components['harmonic']
    return {'kind': 'gaussian_vpt2', 'method': method, 'basis': basis,
            'dispersion': dispersion.upper(),
            'optimized_in_job': bool(re.search(r'\bOpt\s*(?:=|\()', route,
                                               re.IGNORECASE)),
            'zpe_cm_inverse': components,
            'anharmonic_correction_cm_inverse': correction,
            'anharmonic_correction_hartree': correction * invcm / Hartree,
            'harmonic_fundamentals_cm_inverse': harmonic_modes,
            'anharmonic_fundamentals_cm_inverse': anharmonic_modes,
            'fundamental_mode_labels': mode_labels,
            'warnings': warnings, 'review_required': bool(warnings)}


def validate_result_parser(request, *, backend, template, outputs):
    """Reject a parser whose declared method is inconsistent with its input."""
    if not isinstance(request, dict) or not isinstance(request.get('file'), str) \
            or request['file'] not in outputs:
        raise ValueError('invalid result_parser.')
    kind = request.get('kind')
    if kind == 'cfour_dboc':
        valid = (set(request) in ({'kind', 'file', 'level'},
                                  {'kind', 'file', 'level', 'basis'})
                 and backend == 'cfour' and request.get('level') in ('HF', 'MP1')
                 and re.search(r'\bDBOC\s*=\s*ON\b', template,
                               re.IGNORECASE) is not None)
        if valid and 'basis' in request:
            valid = (isinstance(request['basis'], str) and bool(request['basis'])
                     and re.search(rf'\bBASIS\s*=\s*{re.escape(request["basis"])}(?=\s*[,\n)])',
                                   template, re.IGNORECASE) is not None)
    elif kind == 'molpro_energy':
        valid = (set(request) == {'kind', 'file', 'method', 'basis'}
                 and backend == 'molpro'
                 and request.get('method') in ('CCSD(T)', 'CCSD(T)-F12b')
                 and isinstance(request.get('basis'), str) and bool(request['basis'])
                 and re.search(rf'^\s*basis\s*=\s*{re.escape(request["basis"])}\s*$',
                               template, re.IGNORECASE | re.MULTILINE) is not None)
        if valid:
            method_line = (_MOLPRO_CCSD_T if request['method'] == 'CCSD(T)'
                           else _MOLPRO_F12B)
            valid = re.search(method_line, template,
                              re.IGNORECASE | re.MULTILINE) is not None
    elif kind == 'molpro_harmonic':
        valid = (set(request) == {'kind', 'file', 'basis'}
                 and backend == 'molpro' and isinstance(request.get('basis'), str)
                 and bool(request['basis'])
                 and re.search(rf'^\s*basis\s*=\s*{re.escape(request["basis"])}\s*$',
                               template, re.IGNORECASE | re.MULTILINE) is not None
                 and re.search(_MOLPRO_CCSD_T, template,
                               re.IGNORECASE | re.MULTILINE) is not None
                 and bool(re.search(r'^\s*frequencies\s*,\s*numerical\s*$',
                                    template, re.IGNORECASE | re.MULTILINE)))
    elif kind == 'gaussian_vpt2':
        valid = (set(request) in ({'kind', 'file', 'method', 'basis'},
                                  {'kind', 'file', 'method', 'basis', 'dispersion'})
                 and backend == 'gaussian'
                 and isinstance(request.get('method'), str)
                 and isinstance(request.get('basis'), str)
                 and bool(request['method']) and bool(request['basis'])
                 and isinstance(request.get('dispersion', ''), str)
                 and re.search(rf'\b{re.escape(request["method"])}/'
                               rf'{re.escape(request["basis"])}(?![\w-])', template,
                               re.IGNORECASE) is not None
                 and re.search(r'\bFreq\s*=\s*Anharmonic\b', template,
                               re.IGNORECASE) is not None
                 and _gaussian_dispersion(template)
                 == request.get('dispersion', '').upper())
    else:
        valid = False
    if not valid:
        raise ValueError('invalid result_parser.')


def parse_result(output, request):
    kind = request['kind']
    if kind == 'cfour_dboc':
        return parse_cfour_dboc(output, level=request['level'],
                                basis=request.get('basis'))
    if kind == 'molpro_energy':
        return parse_molpro_energy(output, method=request['method'],
                                   basis=request['basis'])
    if kind == 'molpro_harmonic':
        return parse_molpro_harmonic(output, basis=request['basis'])
    if kind == 'gaussian_vpt2':
        return parse_gaussian_vpt2(output, method=request['method'],
                                   basis=request['basis'],
                                   dispersion=request.get('dispersion', ''))
    raise ValueError(f'Unsupported result parser {kind!r}.')


def main(argv=None):
    parser = argparse.ArgumentParser(description='Inspect native ANL QC output')
    subparsers = parser.add_subparsers(dest='kind', required=True)
    cfour = subparsers.add_parser('cfour-dboc')
    cfour.add_argument('output', type=Path)
    cfour.add_argument('--level', choices=('HF', 'MP1'), default='HF')
    energy = subparsers.add_parser('molpro-energy')
    energy.add_argument('output', type=Path)
    energy.add_argument('--method', choices=('CCSD(T)', 'CCSD(T)-F12b'), required=True)
    energy.add_argument('--basis', required=True)
    harmonic = subparsers.add_parser('molpro-harmonic')
    harmonic.add_argument('output', type=Path)
    harmonic.add_argument('--basis', required=True)
    vpt2 = subparsers.add_parser('gaussian-vpt2')
    vpt2.add_argument('output', type=Path)
    vpt2.add_argument('--method', required=True)
    vpt2.add_argument('--basis', required=True)
    vpt2.add_argument('--dispersion', default='')
    args = parser.parse_args(argv)
    if args.kind == 'cfour-dboc':
        print(json.dumps(parse_cfour_dboc(args.output.read_text(errors='replace'),
                                         level=args.level), indent=2))
    elif args.kind == 'molpro-energy':
        print(json.dumps(parse_molpro_energy(args.output.read_text(errors='replace'),
                                             method=args.method, basis=args.basis), indent=2))
    elif args.kind == 'molpro-harmonic':
        print(json.dumps(parse_molpro_harmonic(args.output.read_text(errors='replace'),
                                               basis=args.basis), indent=2))
    elif args.kind == 'gaussian-vpt2':
        print(json.dumps(parse_gaussian_vpt2(args.output.read_text(errors='replace'),
                                             method=args.method, basis=args.basis,
                                             dispersion=args.dispersion), indent=2))


if __name__ == '__main__':
    main()
