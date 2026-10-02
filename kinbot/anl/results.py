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
_MOLPRO_CCSD_T = (r'^\s*\{?\s*uccsd\(t\)\s*,[^\n]*'
                  r'\buhf_uccsd\s*=\s*1\b[^\n]*$')
_MOLPRO_F12B = (r'^\s*\{?\s*uccsd\(t\)-f12b\b[^\n]*'
                  r'\bscale_trip\s*=\s*1\b[^\n]*$')
_MOLPRO_LEGACY_CCSD_T = r'^\s*ccsd\(t\)(?:\s*,[^\n]*)?\s*$'
_MOLPRO_LEGACY_F12B = (r'^\s*ccsd\(t\)-f12b?\b[^\n]*'
                        r'\bscale_trip\s*=\s*1\b[^\n]*$')
_MOLPRO_RHF = r'^\s*\{?\s*rhf(?:\s*,[^\n]*)?\s*$'
_MOLPRO_ALL_ELECTRON = r'(?:^|;)\s*core\s*}'
_MRCC_METHODS = ('CCSDT(Q)', 'CCSDTQ(P)')


def legacy_molpro_parser(request, text):
    """Identify an old parser declaration from its paired input/output text.

    Some intermediate workflow archives omitted ``reference`` while already
    using the unrestricted command.  The command echo, rather than the absent
    field alone, therefore distinguishes those records from the older
    restricted implementation.
    """
    if request.get('reference') is not None:
        return False
    kind = request.get('kind')
    if kind == 'molpro_energy':
        method = request.get('method')
        modern = _MOLPRO_F12B if method == 'CCSD(T)-F12b' else _MOLPRO_CCSD_T
        old = (_MOLPRO_LEGACY_F12B if method == 'CCSD(T)-F12b'
               else _MOLPRO_LEGACY_CCSD_T)
    elif kind == 'molpro_harmonic':
        modern, old = _MOLPRO_CCSD_T, _MOLPRO_LEGACY_CCSD_T
    else:
        return False
    return (re.search(old, text, re.IGNORECASE | re.MULTILINE) is not None
            and re.search(modern, text, re.IGNORECASE | re.MULTILINE) is None)


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


def parse_cfour_energy(output, *, method, basis, reference, correlation, core,
                       program, driver):
    """Read a native CFOUR unrestricted higher-order single-point energy."""
    if (method != 'CCSDT(Q)' or reference != 'RHF'
            or correlation != 'unrestricted' or core != 'frozen'
            or program != 'cfour' or driver != 'VCC'):
        raise ValueError('CFOUR higher-order method settings are unsupported.')
    if ('ERROR ERROR ERROR' in output
            or 'Job has terminated with error flag' in output):
        raise ValueError('CFOUR reported an error flag.')
    if 'This computation required' not in output:
        raise ValueError('CFOUR output has no final completion line.')
    expected = {
        'CALC_LEVEL': method, 'BASIS': basis, 'REFERENCE': reference,
        'CC_PROGRAM': driver, 'FROZEN_CORE': 'ON',
    }
    for keyword, value in expected.items():
        matches = re.findall(
            rf'^\s*{keyword}\s*=\s*{re.escape(value)}\s*$', output,
            re.IGNORECASE | re.MULTILINE)
        if len(matches) != 1:
            raise ValueError(f'CFOUR output does not uniquely echo '
                             f'{keyword}={value}.')
    finals = _FINAL_ENERGY.findall(output)
    if len(finals) != 1:
        raise ValueError('CFOUR needs exactly one final electronic energy.')
    return {
        'kind': 'cfour_energy', 'method': method, 'basis': basis,
        'reference': reference, 'correlation': correlation, 'core': core,
        'program': program, 'driver': driver,
        'program_variant': 'RHF-UCCSDT(Q)',
        'energy_hartree': _number(finals[0]),
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


def _parse_legacy_molpro_energy(output, *, method, basis):
    """Reproduce parser records written by pre-unrestricted workflows.

    This path exists only so completed dispatcher archives remain readable.
    Its deliberately smaller result has no reference/correlation provenance,
    so current ANL recipes cannot mistake it for an unrestricted component.
    """
    command = (_MOLPRO_LEGACY_F12B if method == 'CCSD(T)-F12b'
               else _MOLPRO_LEGACY_CCSD_T)
    if not re.search(command, output, re.IGNORECASE | re.MULTILINE):
        raise ValueError('Legacy Molpro output does not echo its declared method.')
    label = (r'CCSD\(T\)-F12b' if method == 'CCSD(T)-F12b'
             else r'CCSD\(T\)')
    energy = _one_number(
        output, rf'^\s*!{label} total energy\s+({_NUMBER})\s*$',
        f'legacy Molpro {method} total energy')
    if method == 'CCSD(T)-F12b':
        summary = _one_number(
            output, rf'^\s*CCSD\(T\)-F12/{re.escape(basis)} energy\s*=\s*'
                    rf'({_NUMBER})\s*$', 'legacy Molpro F12 summary')
        if abs(summary - energy) > 1e-9:
            raise ValueError('Legacy Molpro F12b energy disagrees with its summary.')
    return {'kind': 'molpro_energy', 'method': method, 'basis': basis,
            'energy_hartree': energy}


def parse_molpro_energy(output, *, method, basis, reference=None, core=None,
                        relativistic=None, legacy=False):
    """Read the exact named total energy, rather than a rounded variable or F12a."""
    _molpro_output(output, basis)
    if method not in ('CCSD(T)', 'CCSD(T)-F12b'):
        raise ValueError(f'Unsupported Molpro energy method {method!r}.')
    if legacy:
        if reference is not None or core is not None or relativistic is not None:
            raise ValueError('Legacy Molpro parser cannot declare modern settings.')
        return _parse_legacy_molpro_energy(output, method=method, basis=basis)
    if method == 'CCSD(T)-F12b':
        if not re.search(_MOLPRO_F12B,
                         output, re.IGNORECASE | re.MULTILINE):
            raise ValueError('Molpro output does not echo scaled-triples F12 input.')
    else:
        if not re.search(_MOLPRO_CCSD_T, output,
                         re.IGNORECASE | re.MULTILINE):
            raise ValueError('Molpro output does not echo conventional '
                             'RHF-UCCSD(T) with UHF_UCCSD=1.')
    if reference is not None:
        if reference not in ('RHF', 'ROHF'):
            raise ValueError(f'Unsupported Molpro reference {reference!r}.')
        if (not re.search(_MOLPRO_RHF, output, re.IGNORECASE | re.MULTILINE)
                or not re.search(r'^\s*PROGRAMS?\s+\*.*\bRHF-SCF\b.*$',
                                 output, re.IGNORECASE | re.MULTILINE)):
            raise ValueError('Molpro output does not contain the requested '
                             'restricted HF reference calculation.')
    if core not in (None, 'frozen', 'all-electron'):
        raise ValueError(f'Unsupported Molpro core treatment {core!r}.')
    # Molpro's documented all-electron form is ``{ccsd(t);core}``.  Echoed
    # input can retain that one-line spelling or wrap the local directive onto
    # its own line, so accept either representation while still requiring the
    # directive to close the coupled-cluster command block.
    all_electron = re.search(_MOLPRO_ALL_ELECTRON, output,
                             re.IGNORECASE | re.MULTILINE) is not None
    if core == 'all-electron' and not all_electron:
        raise ValueError('Molpro output does not echo the all-electron core directive.')
    if core == 'frozen' and all_electron:
        raise ValueError('Molpro frozen-core output contains an all-electron directive.')
    if relativistic not in (None, 'none', 'DKH2'):
        raise ValueError(f'Unsupported Molpro relativistic setting {relativistic!r}.')
    dkh2 = re.search(r'^\s*set\s*,\s*dkho\s*=\s*2\s*$', output,
                     re.IGNORECASE | re.MULTILINE) is not None
    if relativistic == 'DKH2' and not dkh2:
        raise ValueError('Molpro output does not echo SET,DKHO=2.')
    if relativistic == 'none' and dkh2:
        raise ValueError('Molpro nonrelativistic output contains SET,DKHO=2.')
    method_pattern = (_MOLPRO_F12B if method == 'CCSD(T)-F12b'
                      else _MOLPRO_CCSD_T)
    if not re.search(method_pattern, output, re.IGNORECASE | re.MULTILINE):
        raise ValueError('Molpro output does not echo the requested '
                         'RHF-UCCSD(T) method.')
    output_label = (r'RHF-UCCSD\(T\)-F12'
                    if method == 'CCSD(T)-F12b'
                    else r'RHF-UCCSD\(T\)')
    energy = _one_number(
        output, rf'^\s*!{output_label}(?: total)? energy\s+({_NUMBER})\s*$',
        f'Molpro {method} total energy')
    if method == 'CCSD(T)-F12b':
        evidence = {
            'unrestricted coupled-cluster startup':
                r'^\s*Starting UCCSD calculation\s*$',
            'F12b unrestricted correlation result':
                rf'^\s*UCCSD-F12b correlation energy\s+{_NUMBER}\s*$',
            'unrestricted F12 program summary':
                r'^\s*PROGRAMS\s+\*.*\bUCCSD\(T\)(?=\s|$).*'
                r'\bRHF-SCF\b.*$',
        }
        for label, pattern in evidence.items():
            if not re.search(pattern, output, re.IGNORECASE | re.MULTILINE):
                raise ValueError(f'Molpro output lacks {label}.')
        restricted = (r'^\s*Starting RCCSD calculation\s*$|'
                      r'^\s*!RCCSD\(T\)-F12(?:[ab])?\s+energy\b|'
                      r'^\s*!RHF-RCCSD\(T\)-F12(?:[ab])?\s+energy\b')
        if re.search(restricted, output, re.IGNORECASE | re.MULTILINE):
            raise ValueError('Molpro output contains restricted F12 coupled '
                             'cluster evidence.')
        stored = re.findall(
            rf'^\s*(?:SETTING\s+)?KB_F12B\s*=\s*({_NUMBER})(?:\s+AU)?\s*$',
            output, re.IGNORECASE | re.MULTILINE)
        if stored and (len(stored) != 1
                       or abs(_number(stored[0]) - energy) > 1e-7):
            raise ValueError('Molpro stored F12b energy disagrees with the '
                             'native total energy.')
    result = {'kind': 'molpro_energy', 'method': method, 'basis': basis,
              'energy_hartree': energy}
    if reference is not None:
        result['reference'] = reference
    result['program_variant'] = ('RHF-UCCSD(T)-F12b'
                                 if method == 'CCSD(T)-F12b'
                                 else 'RHF-UCCSD(T)')
    result['correlation'] = 'unrestricted'
    if core is not None:
        result['core'] = core
    if relativistic is not None:
        result['relativistic'] = relativistic
    return result


def parse_mrcc_energy(output, *, method, basis, reference, correlation, core,
                      program):
    """Read one final total energy from a direct MRCC calculation."""
    if method not in _MRCC_METHODS:
        raise ValueError(f'Unsupported direct MRCC method {method!r}.')
    if (reference not in ('RHF', 'ROHF', 'UHF')
            or correlation != 'unrestricted' or core != 'frozen'
            or program != 'mrcc'):
        raise ValueError('MRCC reference or core treatment is unsupported.')
    normal = re.findall(r'^\s*Normal termination of mrcc\.\s*$', output,
                        re.IGNORECASE | re.MULTILINE)
    if len(normal) != 1:
        raise ValueError('MRCC output needs exactly one normal termination.')
    if re.search(r'Error at the termination of mrcc|\bFatal error\b', output,
                 re.IGNORECASE):
        raise ValueError('MRCC reported a fatal or termination error.')
    if not re.search(rf'^\s*basis\s*=\s*{re.escape(basis)}\s*$', output,
                     re.IGNORECASE | re.MULTILINE):
        raise ValueError(f'MRCC output does not echo basis {basis}.')
    if not re.search(rf'^\s*calc\s*=\s*{re.escape(method)}\s*$', output,
                     re.IGNORECASE | re.MULTILINE):
        raise ValueError(f'MRCC output does not echo calc={method}.')
    if not re.search(r'^\s*ccprog\s*=\s*mrcc\s*$', output,
                     re.IGNORECASE | re.MULTILINE):
        raise ValueError('MRCC output does not echo ccprog=mrcc.')
    if not re.search(rf'^\s*scftype\s*=\s*{reference}\s*$', output,
                     re.IGNORECASE | re.MULTILINE):
        raise ValueError(f'MRCC output does not echo scftype={reference}.')
    if (reference == 'ROHF'
            and not re.search(r'^\s*rohftype\s*=\s*semicanonical\s*$', output,
                              re.IGNORECASE | re.MULTILINE)):
        raise ValueError('MRCC output does not echo semicanonical ROHF orbitals.')
    if (reference == 'ROHF'
            and not re.search(r'^\s*rohfcore\s*=\s*semicanonical\s*$', output,
                              re.IGNORECASE | re.MULTILINE)):
        raise ValueError('MRCC output does not echo semicanonical ROHF core '
                         'orbitals.')
    if not re.search(rf'^\s*core\s*=\s*{core}\s*$', output,
                     re.IGNORECASE | re.MULTILINE):
        raise ValueError(f'MRCC output does not echo core={core}.')
    energy = _one_number(
        output,
        rf'^\s*Total\s+{re.escape(method)}\s+energy\s*\[au\]\s*:\s*'
        rf'({_NUMBER})\s*$', f'MRCC {method} total energy')
    return {'kind': 'mrcc_energy', 'method': method, 'basis': basis,
            'reference': reference, 'correlation': correlation, 'core': core,
            'program': program,
            'energy_hartree': energy, 'driver': 'direct',
            'program_variant': f'{reference}-U{method}'}


def parse_molpro_harmonic(output, *, basis, reference=None, legacy=False):
    """Read vibrational modes and ZPE, excluding rotations/translations."""
    _molpro_output(output, basis)
    command = _MOLPRO_LEGACY_CCSD_T if legacy else _MOLPRO_CCSD_T
    if not re.search(command, output,
                     re.IGNORECASE | re.MULTILINE):
        label = ('legacy CCSD(T)' if legacy else 'forced RHF-UCCSD(T)')
        raise ValueError(f'Molpro harmonic output does not echo the {label} method.')
    if legacy and reference is not None:
        raise ValueError('Legacy Molpro harmonic parser cannot declare a reference.')
    if reference is not None:
        if reference not in ('RHF', 'ROHF'):
            raise ValueError(f'Unsupported Molpro reference {reference!r}.')
        if (not re.search(_MOLPRO_RHF, output, re.IGNORECASE | re.MULTILINE)
                or not re.search(
                    r'^\s*PROGRAMS?\s+\*.*\bRHF-SCF\b.*$', output,
                    re.IGNORECASE | re.MULTILINE)):
            raise ValueError('Molpro harmonic output does not contain the '
                             'requested restricted HF reference calculation.')
    if not re.search(r'^\s*frequencies\s*,\s*numerical\s*$', output,
                     re.IGNORECASE | re.MULTILINE):
        raise ValueError('Molpro output does not echo numerical frequencies.')
    frequency_sections = re.findall(
        r'^\s*PROGRAM \* FREQUENCIES \(Calculation of harmonic '
        r'vibrational spectra for (.+?)\)\s*$', output,
        re.IGNORECASE | re.MULTILINE)
    expected_method = 'CCSD(T)' if legacy else 'UCCSD(T)'
    if len(frequency_sections) != 1 \
            or frequency_sections[0].upper() != expected_method:
        raise ValueError(f'Expected one {expected_method} Molpro frequency section.')
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
    result = {'kind': 'molpro_harmonic', 'method': 'CCSD(T)', 'basis': basis,
            'wavenumbers_cm_inverse': modes,
            'low_modes_cm_inverse': low_modes,
            'review_required': any(abs(value) > 1.0 for value in low_modes),
            'zpe': {'hartree': hartree, 'cm_inverse': cm_inverse,
                    'kj_mol': kj_mol}}
    if reference is not None:
        result['reference'] = reference
    if not legacy:
        result['program_variant'] = 'RHF-UCCSD(T)'
        result['correlation'] = 'unrestricted'
    return result


def _gaussian_dispersions(output):
    """Return every dispersion model echoed by Gaussian.

    Gaussian can omit this long keyword from the short route near the start
    while retaining it in the archive record written near the end. Search the
    complete, normally terminated output and require the requested model to be
    present rather than trusting only the first fixed-size text window.
    """
    return {match.upper() for match in re.findall(
        r'\bEmpiricalDispersion\s*(?:=|\()\s*([A-Za-z0-9]+)',
        output, re.IGNORECASE)}


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
    reported_dispersion = _gaussian_dispersions(output)
    requested_dispersion = dispersion.upper()
    if ((requested_dispersion and requested_dispersion not in reported_dispersion)
            or (not requested_dispersion and reported_dispersion)):
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
    elif kind == 'cfour_energy':
        method = request.get('method')
        basis = request.get('basis')
        valid = (
            set(request) == {'kind', 'file', 'method', 'basis', 'reference',
                             'correlation', 'core', 'program', 'driver'}
            and backend == 'cfour' and method == 'CCSDT(Q)'
            and isinstance(basis, str) and bool(basis)
            and request.get('reference') == 'RHF'
            and request.get('correlation') == 'unrestricted'
            and request.get('core') == 'frozen'
            and request.get('program') == 'cfour'
            and request.get('driver') == 'VCC'
            and re.search(r'\bCALC\s*=\s*CCSDT\(Q\)(?=\s*[,\n)])',
                          template, re.IGNORECASE)
            and re.search(rf'\bBASIS\s*=\s*{re.escape(basis)}(?=\s*[,\n)])',
                          template, re.IGNORECASE)
            and re.search(r'\bREFERENCE\s*=\s*RHF(?=\s*[,\n)])',
                          template, re.IGNORECASE)
            and re.search(r'\bCC_PROGRAM\s*=\s*VCC(?=\s*[,\n)])',
                          template, re.IGNORECASE)
            and re.search(r'\bFROZEN_CORE\s*=\s*ON(?=\s*[,\n)])',
                          template, re.IGNORECASE))
    elif kind == 'molpro_energy':
        required = {'kind', 'file', 'method', 'basis'}
        optional = {'reference', 'core', 'relativistic'}
        valid = (required <= set(request) <= required | optional
                 and backend == 'molpro'
                 and request.get('method') in ('CCSD(T)', 'CCSD(T)-F12b')
                 and isinstance(request.get('basis'), str) and bool(request['basis'])
                 and re.search(rf'^\s*basis\s*=\s*{re.escape(request["basis"])}\s*$',
                              template, re.IGNORECASE | re.MULTILINE) is not None)
        legacy = legacy_molpro_parser(request, template)
        if valid and legacy:
            valid = set(request) == required
        if valid and 'reference' in request:
            method_command = (_MOLPRO_F12B
                              if request['method'] == 'CCSD(T)-F12b'
                              else _MOLPRO_CCSD_T)
            valid = (request['reference'] in ('RHF', 'ROHF')
                     and re.search(_MOLPRO_RHF, template,
                                   re.IGNORECASE | re.MULTILINE) is not None
                     and re.search(method_command, template,
                                   re.IGNORECASE | re.MULTILINE) is not None)
        if valid:
            method_line = (
                (_MOLPRO_LEGACY_CCSD_T if request['method'] == 'CCSD(T)'
                 else _MOLPRO_LEGACY_F12B) if legacy else
                (_MOLPRO_CCSD_T if request['method'] == 'CCSD(T)'
                 else _MOLPRO_F12B))
            valid = re.search(method_line, template,
                              re.IGNORECASE | re.MULTILINE) is not None
        if valid and 'core' in request:
            core_line = re.search(_MOLPRO_ALL_ELECTRON, template,
                                  re.IGNORECASE | re.MULTILINE)
            valid = (request['core'] in ('frozen', 'all-electron')
                     and ((request['core'] == 'all-electron') == bool(core_line)))
        if valid and 'relativistic' in request:
            dkh2 = re.search(r'^\s*set\s*,\s*dkho\s*=\s*2\s*$', template,
                             re.IGNORECASE | re.MULTILINE)
            valid = (request['relativistic'] in ('none', 'DKH2')
                     and ((request['relativistic'] == 'DKH2') == bool(dkh2)))
    elif kind == 'molpro_harmonic':
        legacy = legacy_molpro_parser(request, template)
        valid = (set(request) in ({'kind', 'file', 'basis'},
                                  {'kind', 'file', 'basis', 'reference'})
                 and backend == 'molpro' and isinstance(request.get('basis'), str)
                 and bool(request['basis'])
                 and re.search(rf'^\s*basis\s*=\s*{re.escape(request["basis"])}\s*$',
                               template, re.IGNORECASE | re.MULTILINE) is not None
                 and re.search((_MOLPRO_LEGACY_CCSD_T if legacy
                                else _MOLPRO_CCSD_T), template,
                               re.IGNORECASE | re.MULTILINE) is not None
                 and bool(re.search(r'^\s*frequencies\s*,\s*numerical\s*$',
                                    template, re.IGNORECASE | re.MULTILINE)))
        if valid and 'reference' in request:
            command = _MOLPRO_CCSD_T
            valid = (request['reference'] in ('RHF', 'ROHF')
                     and re.search(_MOLPRO_RHF, template,
                                   re.IGNORECASE | re.MULTILINE)
                     and re.search(command, template,
                                   re.IGNORECASE | re.MULTILINE))
    elif kind == 'mrcc_energy':
        method = request.get('method')
        basis = request.get('basis')
        reference = request.get('reference')
        correlation = request.get('correlation')
        core = request.get('core')
        program = request.get('program')
        valid = (set(request) == {'kind', 'file', 'method', 'basis',
                                  'reference', 'correlation', 'core', 'program'}
                 and backend == 'mrcc' and method in _MRCC_METHODS
                 and reference in ('RHF', 'ROHF', 'UHF')
                 and correlation == 'unrestricted' and core == 'frozen'
                 and program == 'mrcc'
                 and isinstance(basis, str) and bool(basis)
                 and re.search(rf'^\s*calc\s*=\s*{re.escape(method)}\s*$',
                               template, re.IGNORECASE | re.MULTILINE)
                 and re.search(r'^\s*ccprog\s*=\s*mrcc\s*$', template,
                               re.IGNORECASE | re.MULTILINE)
                 and re.search(rf'^\s*basis\s*=\s*{re.escape(basis)}\s*$',
                               template, re.IGNORECASE | re.MULTILINE)
                 and re.search(rf'^\s*scftype\s*=\s*{reference}\s*$',
                               template, re.IGNORECASE | re.MULTILINE)
                 and re.search(r'^\s*core\s*=\s*frozen\s*$', template,
                               re.IGNORECASE | re.MULTILINE)
                 and re.search(r'^\s*geom\s*=\s*xyz\s*$', template,
                               re.IGNORECASE | re.MULTILINE)
                 and re.search(r'^\s*unit\s*=\s*angs\s*$', template,
                               re.IGNORECASE | re.MULTILINE)
                 and (reference != 'ROHF' or re.search(
                     r'^\s*rohftype\s*=\s*semicanonical\s*$', template,
                     re.IGNORECASE | re.MULTILINE))
                 and (reference != 'ROHF' or re.search(
                     r'^\s*rohfcore\s*=\s*semicanonical\s*$', template,
                     re.IGNORECASE | re.MULTILINE)))
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
                 and _gaussian_dispersions(template)
                 == ({request['dispersion'].upper()}
                     if request.get('dispersion') else set()))
    else:
        valid = False
    if not valid:
        raise ValueError('invalid result_parser.')


def parse_result(output, request):
    kind = request['kind']
    if kind == 'cfour_dboc':
        return parse_cfour_dboc(output, level=request['level'],
                                basis=request.get('basis'))
    if kind == 'cfour_energy':
        return parse_cfour_energy(
            output, method=request['method'], basis=request['basis'],
            reference=request['reference'], correlation=request['correlation'],
            core=request['core'], program=request['program'],
            driver=request['driver'])
    if kind == 'molpro_energy':
        return parse_molpro_energy(output, method=request['method'],
                                   basis=request['basis'],
                                   reference=request.get('reference'),
                                   core=request.get('core'),
                                   relativistic=request.get('relativistic'),
                                   legacy=legacy_molpro_parser(request, output))
    if kind == 'molpro_harmonic':
        return parse_molpro_harmonic(output, basis=request['basis'],
                                     reference=request.get('reference'),
                                     legacy=legacy_molpro_parser(request, output))
    if kind == 'mrcc_energy':
        return parse_mrcc_energy(output, method=request['method'],
                                 basis=request['basis'],
                                 reference=request['reference'],
                                 correlation=request['correlation'],
                                 core=request['core'],
                                 program=request['program'])
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
