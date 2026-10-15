"""Reproduce the NH2O SI structural/counting audit without vendoring its data.

Run: python tests/verify_nh2o_si.py /path/to/downloaded/SI report.json
SI: Klippenstein et al., 10.1021/acs.jpca.3c05181.s002 (CC BY-NC 4.0).
This verifies input structure/algebra, not native MESS rates or Ne file aliases.
"""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import re
import tempfile

import numpy as np

from kinbot.vrc_multistate import assemble_states


NETWORK_SHA256 = 'f76e0cdb130618f55850ade3e61f8268ea38e95c1ce877fb209f3d2b9e648403'


def _levels(text):
    match = re.search(r'ElectronicLevels\[1/cm\]\s+(\d+)\s*\n', text)
    if match is None:
        raise ValueError('Missing SI electronic level block.')
    values = []
    for line in text[match.end():].splitlines():
        words = re.split('[!#]', line)[0].split()
        if not words:
            continue
        values.append({'energy_cm_inverse': float(words[0]), 'degeneracy': int(words[1])})
        if len(values) == int(match[1]):
            return values
    raise ValueError('Incomplete SI electronic level block.')


def verify(directory):
    directory = Path(directory)
    network = directory / 'jp3c05181_si_002.txt'
    digest = hashlib.sha256(network.read_bytes()).hexdigest()
    if digest != NETWORK_SHA256:
        raise ValueError('SI network differs from the publisher-checksummed reference.')
    text = network.read_text()
    headers = list(re.finditer(r'^\s*(Well|Bimolecular|Barrier)\s+([^\n]+)', text, re.M))
    blocks = {}
    for index, header in enumerate(headers):
        words = re.split('[!#]', header[2])[0].split()
        blocks[words[0]] = (header[1], words[1:], text[header.end():
            headers[index + 1].start() if index + 1 < len(headers) else len(text)])
    expected = {'B1v': ['W1v', 'P1'], 'B1': ['W1v', 'W1'],
                'B2v': ['W2v', 'P3'], 'B2': ['W2', 'W2v']}
    report = {}
    for name, endpoints in expected.items():
        kind, actual, block = blocks[name]
        if kind != 'Barrier' or actual != endpoints:
            raise ValueError('SI entrance network topology changed.')
        rrho = re.split(r'^\s*RRHO\b[^\n]*\n', block, flags=re.M)[1:]
        components = []
        for part in rrho:
            file = re.search(r'\bFile\s+(\S+)', part)
            zero = re.search(r'ZeroEnergy\[kcal/mol\]\s+([-+\d.]+)', part)
            symmetry = re.search(r'SymmetryFactor\s+([-+\d.]+)', part)
            components.append({'file': file[1], 'zero_energy_kcal_mol': float(zero[1]),
                'symmetry_factor': float(symmetry[1]), 'electronic_levels': _levels(part)})
        report[name] = {'endpoints': actual, 'components': components}
    if [len(report[name]['components']) for name in expected] != [3, 2, 2, 1]:
        raise ValueError('SI state component counts changed.')
    degeneracy = defaultdict(int)
    for component in report['B1v']['components']:
        for level in component['electronic_levels']:
            degeneracy[level['energy_cm_inverse']] += level['degeneracy']
    if dict(degeneracy) != {0.: 10, 158.5: 6, 226.5: 2}:
        raise ValueError('SI long-range level normalization changed.')
    # Test physical-level summation using an actual SI Ne table, reused for
    # all three components to isolate electronic populations from PES shape.
    # No claim is made about which original table filename this download maps to.
    table = directory / 'jp3c05181_si_003.txt'
    native = np.loadtxt(table)
    temperatures = [100., 300., 1000.]
    errors = []
    with tempfile.TemporaryDirectory() as folder:
        target = Path(folder) / 'sum.out'
        assemble_states([{'source': table, 'multiplicity': 2, 'energy_offset_cm_inverse': 0.,
                          'electronic_levels': component['electronic_levels']}
                         for component in report['B1v']['components']], target, baseline_cm_inverse=0.)
        summed = np.loadtxt(target)
        for temperature in temperatures:
            beta = 1.438776877 / temperature
            def integrate(data):
                integrand = data[:, 1] * np.exp(-beta * data[:, 0])
                return np.sum((integrand[1:] + integrand[:-1]) * np.diff(data[:, 0]) / 2.)
            ratio = integrate(summed) / integrate(native)
            expected_ratio = 2 * (5 + 3 * np.exp(-beta * 158.5) + np.exp(-beta * 226.5))
            errors.append(abs(ratio / expected_ratio - 1))
    if max(errors) > .003:
        raise ValueError('Electronic flux convolution disagrees with the SI level partition.')
    table_files = sorted(directory.glob('jp3c05181_si_*.txt'))
    for path in table_files:
        if path != network:
            values = np.loadtxt(path)
            if (values.ndim != 2 or values.shape[1] != 2 or not np.isfinite(values).all()
                    or np.any(values[:, 1] < 0.) or np.any(np.diff(values[:, 0]) <= 0.)):
                raise ValueError(f'Invalid SI Ne table: {path.name}')
    return {'source_doi': '10.1021/acs.jpca.3c05181.s002',
        'verification': 'SI topology, explicit level counting and actual-table convolution; no native rate claim',
        'network_sha256': digest, 'entrance_barriers': report,
        'outer_NH2_O_total_level_degeneracies': dict(degeneracy),
        'canonical_level_counting_relative_errors': errors, 'temperatures_kelvin': temperatures,
        'files_sha256': {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in table_files},
        'deferred_for_sulfur_examples': 'automatic homolytic vdW discovery/stabilization',
        'status': 'passed'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory')
    parser.add_argument('report')
    args = parser.parse_args()
    result = verify(args.directory)
    Path(args.report).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
