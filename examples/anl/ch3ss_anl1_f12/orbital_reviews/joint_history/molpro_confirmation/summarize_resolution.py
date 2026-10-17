"""Reproduce the recovery/gate review from the retained native exports."""
import argparse
import json
from pathlib import Path
import re

import numpy as np

from kinbot.execution import atomic_json, file_sha256
from review_export import review, KCAL, NUMBER


def mrci_calls(path):
    calls = []
    for block in re.split(r'PROGRAM\s*\*\s*MRCI-F12\b', path.read_text(), flags=re.I)[1:]:
        energy = re.search(r'!CI-F12 STATE\s+\d+\.\d+\s+Energy\s+('+NUMBER+')', block)
        dipole = re.search(r'!CI-F12 STATE\s+\d+\.\d+\s+Dipole moment\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+')', block)
        iterations = [line.split() for line in block.splitlines()
                      if re.match(r'^\s*\d+\s+\d+\s+\d+\s+\d+\.\d+',line)]
        calls.append({'variational_hartree': float(energy[1]),
                      'dipole_au': [float(x) for x in dipole.groups()],
                      'initial_internal_energy_hartree': float(iterations[0][5]),
                      'iterations': int(iterations[-1][0])})
    return calls


def summarize(recovery, gates, output):
    output.mkdir(parents=True, exist_ok=True)
    rr, gr = review(recovery), review(gates)
    atomic_json(output/'native_recovery_review.json', rr)
    atomic_json(output/'native_cabs3_review.json', gr)
    result = {'schema': 1, 'production_ready': False,
        'recovery_verified_files': rr['verified_files'], 'gates_verified_files': gr['verified_files'],
        'recovery_complete_points': sum(j['status']=='complete' for s in rr['systems'] for j in s['jobs']),
        'gate_complete_points': sum(j['status']=='complete' for s in gr['systems'] for j in s['jobs']),
        'caspt2': [], 'cabs3': [],
        'remaining': ['Simultaneous MRCI roots at the interacting geometries, with native state/energy checks.',
                      'Consistent CABS3 far/center/tilted interactions and matched CABS1 tilted control.',
                      'Install accepted basis-specific models in production JSON, then user-run end-to-end checks.']}
    returns={c['chain']: c for c in rr['return_checks']}
    for chain in rr['systems']:
        if 'caspt2' not in chain['title']:
            continue
        rows=[]
        far=next(j for j in chain['jobs'] if j['name'].endswith('__far'))
        ef=np.asarray(far['correlated_hartree'])
        n=len(ef)//2
        for job in chain['jobs']:
            e=np.asarray(job['correlated_hartree'])
            t=(recovery/job['native_output']).read_text()
            rows.append({'name': job['name'], 'output_sha256': job['output_sha256'],
                'wall_seconds': float(re.findall(r'REAL TIME\s*\*\s*([\d.]+) SEC', t)[-1]),
                'lowest_doublet_interaction_kcal': float((min(e[:n])-min(ef[:n]))*KCAL),
                'lowest_doublet_minus_lowest_quartet_kcal': float((min(e[:n])-min(e[n:]))*KCAL),
                'reference_counts': sorted(set(map(int,re.findall(r'Number of reference states:\s*(\d+)',t)))),
                'native_reference_weights': sorted(set(re.findall(r'Weight factors for SA-density in H0:\s*([^\n]+)',t))),
                'frozen_core_counts': sorted(set(map(int,re.findall(r'Number of core orbitals:\s*(\d+)',t))))})
        group_count = 4 if n==6 else 2
        result['caspt2'].append({'chain':chain['title'], 'active_labels': chain['active_labels'],
            'far_sulfur_component_spans_kcal': [float(np.ptp(a)*KCAL) for a in np.split(ef,group_count)],
            'far_fragment_gaps_kcal': ([float((np.mean(a[3:])-np.mean(a[:3]))*KCAL) for a in np.split(ef,2)] if n==6 else []),
            'return_check': returns[chain['title']], 'points': rows})
    for chain in gr['systems']:
        job=chain['jobs'][0]; e=np.asarray(job['correlated_hartree']); n=len(e)//2
        result['cabs3'].append({'chain':chain['title'], 'output_sha256': job['output_sha256'],
            'roots':len(e), 'native_alerts':job['native_alerts'],
            'same_component_D_minus_Q_kcal': ((e[:n]-e[n:])*KCAL).tolist(),
            'sulfur_component_spans_kcal': [float(np.ptp(a)*KCAL) for a in np.split(e,4 if n==6 else 2)],
            'energy_selector_checks':job['mrci_energy_checks']})
    native=recovery/'ch3s_cs_atz_mrci_maxdav20/atz__tilted/atz__tilted.out'
    calls=mrci_calls(native)
    pairs=[]
    for i,j in [(1,2),(3,5)]:
        a,b=calls[i],calls[j]
        pairs.append({'calls':[i+1,j+1], 'native_variational_difference_hartree': b['variational_hartree']-a['variational_hartree'],
            'dipole_max_component_difference_au':float(np.max(abs(np.array(a['dipole_au'])-b['dipole_au']))),
            'initial_internal_energy_separation_hartree':b['initial_internal_energy_hartree']-a['initial_internal_energy_hartree']})
    result['mrci_root_collapse']={'native_output':str(native), 'output_sha256':file_sha256(native),
        'calls':calls,'duplicate_pairs':pairs,
        'interpretation':'Separate requested roots reach effectively duplicate variational energies and dipoles after distinct initial energies. Numerical convergence does not establish distinct eigenstates; simultaneous same-spin roots are required for the next native check.'}
    atomic_json(output/'native_resolution_findings.json',result)
    print(json.dumps({k:result[k] for k in ['recovery_verified_files','gates_verified_files','recovery_complete_points','gate_complete_points','production_ready']},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('recovery',type=Path)
    p.add_argument('gates',type=Path)
    p.add_argument('--output',type=Path,default=Path(__file__).parent)
    a=p.parse_args()
    summarize(a.recovery,a.gates,a.output)
