"""Reproduce the matched CASPT2 shift review from two collected native exports.

Inputs are extracted collect archives, with their export_inventory.json files.
All inventoried bytes are verified by review_export before comparison. This
script diagnoses a matched geometry; it does not certify a PES or a rate error.
"""
import argparse
import json
from pathlib import Path
import re

import numpy as np

from kinbot.execution import atomic_json, file_sha256
from review_export import KCAL, NUMBER, review, same_geometry_overlap


def require(condition, message):
    if not condition:
        raise ValueError(message)


def correlation_calls(text):
    calls = []
    for block in re.split(r'PROGRAM\s*\*\s*RS2-F12\b', text, flags=re.I)[1:]:
        # The last full iteration includes VAR(S/P); the shorter post-shift
        # correction table is deliberately excluded.
        rows = re.findall(
            r'(?m)^\s*(\d+)\s+1\s+(\d+)\s+('+NUMBER+r')\s+'
            r'('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')\s+'
            r'('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')\s*$', block)
        require(bool(rows), 'Missing full native iteration table')
        refs = re.findall(r'(?m)^\s*Reference energy\s+('+NUMBER+')', block)
        require(len(refs) == 1, 'Ambiguous reference energy')
        optimized = re.search(r'Number of optimized states:\s*(\d+)\s+Roots:\s*([^\n]+)', block)
        reference = re.search(r'Number of reference states:\s*(\d+)\s+Roots:\s*([^\n]+)', block)
        weights = re.search(r'Weight factors for SA-density in H0:\s*([^\n]+)', block)
        require(optimized and int(optimized[1]) == 1 and reference and int(reference[1]) == 6,
                'Unexpected optimized/reference manifold')
        require([int(v) for v in reference[2].split()] == list(range(1, 7)), 'Unexpected reference roots')
        require(weights and np.allclose([float(v) for v in weights[1].split()], [1/6]*6, atol=1e-6),
                'Unexpected reference-density weights')
        last = rows[-1]
        calls.append({'root': int(optimized[2]), 'iterations': int(last[0]),
                      'squared_norm': float(last[2]), 'reference_energy_hartree': float(refs[0])})
    require(len(calls) == 12, 'Expected twelve independently converged correlated calls')
    return calls


def compare(previous_root, current_root):
    old_root, new_root = Path(previous_root).resolve(), Path(current_root).resolve()
    old, new = review(old_root), review(new_root)
    report = {'schema': 1, 'status': 'matched_point_review_complete',
              'production_end_to_end_validated': False,
              'verified_files': {'previous': old['verified_files'], 'current': new['verified_files']},
              'inventory_sha256': {'previous': file_sha256(old_root/'export_inventory.json'),
                                   'current': file_sha256(new_root/'export_inventory.json')},
              'cases': [], 'interpretation':
              'The shift-0.1 matched points converge with unchanged CASSCF/reference states and active subspaces. '
              'Use 0.1 consistently in production; this does not establish shift independence or a rate uncertainty.'}
    gaps = {}
    for basis in ('adz', 'atz'):
        old_title = f'ch3ss_ss_{basis}_caspt2_shared_reference'
        new_title = f'ch3ss_ss_{basis}_caspt2_shift01'
        a = next(c for c in old['systems'] if c['title'] == old_title)
        b = next(c for c in new['systems'] if c['title'] == new_title)
        name = basis+'__inner'
        aj = next(j for j in a['jobs'] if j['name'] == name)
        bj = next(j for j in b['jobs'] if j['name'] == name)
        require(aj['status'] == bj['status'] == 'complete', 'Incomplete matched point')
        require(not aj['native_alerts'] and not bj['native_alerts'], 'Native diagnostic requires review')
        ad, bd = old_root/old_title/name, new_root/new_title/name
        native_a = (ad/(name+'.native.inp')).read_text()
        native_b = (bd/(name+'.native.inp')).read_text()
        def normalize(text):
            text = re.sub(r'(?m)^\*\*\*,.*$', '***,MATCHED', text)
            text = re.sub(r'(?m)^diab,2140\.2\s*\n', '', text)
            return text.replace('shift=0.3', 'shift=0.1')
        require(normalize(native_a) == normalize(native_b), 'Unexpected native-input change')
        # DIAB was removed by the recovery adapter. Record it, and validate its
        # actual effect rather than asserting the input differed only in SHIFT.
        ca, cb = correlation_calls((ad/(name+'.out')).read_text()), correlation_calls((bd/(name+'.out')).read_text())
        require([c['reference_energy_hartree'] for c in ca] == [c['reference_energy_hartree'] for c in cb],
                'Native reference energies changed')
        ea, eb = np.array(aj['correlated_hartree']), np.array(bj['correlated_hartree'])
        require(ea.shape == eb.shape == (12,), 'Wrong correlated root count')
        overlap = same_geometry_overlap(ad/(name+'.natural.orbitals.xml'), bd/(name+'.natural.orbitals.xml'), 17, 5)
        casscf_delta = float(np.max(abs(np.array(aj['casscf_hartree'])-bj['casscf_hartree'])))
        require(casscf_delta < 1e-7 and min(overlap['active_singular_values']) > .999999,
                'CASSCF or active subspace changed')
        gap_a, gap_b = (ea[:6]-ea[6])*KCAL, (eb[:6]-eb[6])*KCAL
        gaps[basis] = (gap_a, gap_b)
        report['cases'].append({'basis': basis, 'active_indices_one_based': [18,19,20,21,22],
            'native_outputs': [aj['native_output'], bj['native_output']],
            'native_output_sha256': [aj['output_sha256'], bj['output_sha256']],
            'native_input_sha256': [file_sha256(ad/(name+'.native.inp')), file_sha256(bd/(name+'.native.inp'))],
            'additional_input_difference': 'DIAB,2140.2 removed; title changed. Initial checkpoint hash identical.',
            'casscf_max_change_hartree': casscf_delta, 'active_overlap': overlap,
            'old_shift03_calls': ca, 'new_shift01_calls': cb,
            'old_shift03_energies_hartree': ea.tolist(), 'new_shift01_energies_hartree': eb.tolist(),
            'absolute_energy_changes_kcal': ((eb-ea)*KCAL).tolist(),
            'doublet_minus_quartet1_gap_kcal': {'shift03': gap_a.tolist(), 'shift01': gap_b.tolist(),
                                               'change': (gap_b-gap_a).tolist()},
            'equal_label_spin_gap_changes_kcal': (((eb[:6]-eb[6:])-(ea[:6]-ea[6:]))*KCAL).tolist(),
            'state_order_by_energy': {'shift03': [(np.argsort(v)+1).tolist() for v in (ea[:6],ea[6:])],
                                      'shift01': [(np.argsort(v)+1).tolist() for v in (eb[:6],eb[6:])]}})
    power = 3.7
    cbs_a = (3**power*gaps['atz'][0]-2**power*gaps['adz'][0])/(3**power-2**power)
    cbs_b = (3**power*gaps['atz'][1]-2**power*gaps['adz'][1])/(3**power-2**power)
    report['cbs_gap_estimate_kcal'] = {'power': power, 'shift03': cbs_a.tolist(),
        'shift01': cbs_b.tolist(), 'change': (cbs_b-cbs_a).tolist(),
        'qualification': 'Two-basis estimate, not a calibrated CASPT2-F12 extrapolation uncertainty.'}
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('previous_root', type=Path)
    parser.add_argument('current_root', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.previous_root, args.current_root)
    atomic_json(args.output, result)
    print(json.dumps({'status': result['status'], 'verified_files': result['verified_files'],
                      'lowest_CBS_gap_change_kcal': result['cbs_gap_estimate_kcal']['change'][0]}, indent=2))
