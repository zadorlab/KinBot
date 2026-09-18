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


def _number(value):
    number = float(value.replace('D', 'E').replace('d', 'e'))
    if not math.isfinite(number):
        raise ValueError('CFOUR DBOC output contains a nonfinite number.')
    return number


def _level(name):
    folded = name.strip().casefold()
    return {'hartree-fock': 'HF', 'hf': 'HF', 'scf': 'HF',
            'mp1': 'MP1'}.get(folded, name.strip().upper())


def parse_cfour_dboc(output, *, level='HF'):
    """Select one named DBOC level, preserving separately reported levels.

    CFOUR 2.1 prints identical value labels in HF and MP1 summary sections;
    selecting the last label would silently return the wrong ANL component.
    """
    if 'ERROR ERROR ERROR' in output:
        raise ValueError('CFOUR reported an error flag.')
    if 'This computation required' not in output:
        raise ValueError('CFOUR output has no final completion line.')
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


def main(argv=None):
    parser = argparse.ArgumentParser(description='Inspect native ANL QC output')
    subparsers = parser.add_subparsers(dest='kind', required=True)
    cfour = subparsers.add_parser('cfour-dboc')
    cfour.add_argument('output', type=Path)
    cfour.add_argument('--level', choices=('HF', 'MP1'), default='HF')
    args = parser.parse_args(argv)
    if args.kind == 'cfour-dboc':
        print(json.dumps(parse_cfour_dboc(args.output.read_text(errors='replace'),
                                         level=args.level), indent=2))


if __name__ == '__main__':
    main()
