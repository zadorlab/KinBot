"""Read-only review of a collected confirmation export; never accepts production.

Run against an extracted export. Absolute cluster paths are mapped using the
inventory, never followed on this host. NumPy is used only for same-geometry
MO-space comparisons, not to infer overlaps between displaced geometries.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np

from kinbot.execution import file_sha256, atomic_json
from kinbot.molpro_reconnaissance import convergence_diagnostics

KCAL = 627.5094740631
NUMBER = r'[-+]?\d+(?:\.\d*)?(?:[DEde][-+]?\d+)?'


def floats(text):
    return [float(x.replace('D', 'E').replace('d', 'e')) for x in text.split()]


def orbital_data(path):
    tree = ET.parse(path).getroot()
    sets = [e for e in tree.iter() if e.tag.endswith('}orbitals')]
    if len(sets) != 1 or sets[0].get('angular') != 'cartesian':
        raise ValueError('Expected one exported Cartesian MO set: ' + str(path))
    orbitals = list(sets[0])
    coefficients = np.asarray([floats(o.text) for o in orbitals]).T
    atoms = [dict(e.attrib) for e in tree.iter() if e.tag.endswith('}atom')]
    basis = next(e for e in tree.iter() if e.tag.endswith('}basisSet') and e.get('id') == 'ORBITAL')
    basis_identity = re.sub(r"molecule\[@id='[^']+'\]", "molecule[@id='MOLECULE']", ET.tostring(basis).decode())
    return {'atoms': atoms, 'basis': basis_identity, 'coefficients': coefficients,
            'orbitals': [dict(o.attrib) for o in orbitals], 'set': dict(sets[0].attrib)}


def same_geometry_overlap(left, right, inactive, count):
    """Expand right MOs in the complete orthonormal left MO basis.

    With C_right = C_left T and C_left^T S C_left = I, T is the MO overlap.
    This avoids inventing an AO metric for Cartesian XML. It assumes native
    normalization, checks reconstruction/unitarity, and forbids changed bases
    or coordinates. It does NOT independently validate Molpro's AO integrals.
    """
    a, b = orbital_data(left), orbital_data(right)
    if a['atoms'] != b['atoms'] or a['basis'] != b['basis']:
        raise ValueError('Overlap review requires identical coordinates and AO basis.')
    ca, cb = a['coefficients'], b['coefficients']
    transform, _, rank, _ = np.linalg.lstsq(ca, cb, rcond=1.e-12)
    residual = float(np.max(abs(ca @ transform - cb)))
    orth = float(np.max(abs(transform.T @ transform - np.eye(cb.shape[1]))))
    if rank != ca.shape[1] or residual > 1.e-7 or orth > 1.e-5:
        raise ValueError(f'Native MO comparison failed: rank={rank}, residual={residual}, orth={orth}')
    block = transform[inactive:inactive+count, inactive:inactive+count]
    return {'definition': 'Same-geometry MO expansion; assumes native MO normalization',
            'reconstruction_max_error': residual, 'unitarity_max_error': orth,
            'active_singular_values': np.linalg.svd(block, compute_uv=False).tolist(),
            'active_absolute_overlap': abs(block).tolist()}


def orbital_table(text, title, stop):
    match = re.search(title + r'\s*\(STATE AVERAGED\)', text)
    if not match:
        return []
    body = text[match.end():].split(stop, 1)[0]
    result = []
    for line in body.splitlines():
        head = re.match(r'\s*(\d+\.\d+)\s+('+NUMBER+r')\s+('+NUMBER+r')\s+(.*)', line)
        if head:
            result.append({'orbital': head[1], 'occupation': float(head[2]),
                           'energy_hartree': float(head[3]), 'printed_coefficients': head[4]})
        elif result and re.match(r'\s+\d+\s+\d+\s+[spdf]', line):
            result[-1]['printed_coefficients'] += ' ' + line.strip()
    return result


def review(root):
    root = Path(root).resolve()
    inventory = json.loads((root/'export_inventory.json').read_text())
    mapping = {}
    for row in inventory['files']:
        name = Path(row['archive_path'])
        if name.is_absolute() or '..' in name.parts:
            raise ValueError('Unsafe inventory member')
        path = root/name
        if path.is_symlink() or path.stat().st_size != row['bytes'] or file_sha256(path) != row['sha256']:
            raise ValueError('Export inventory mismatch: ' + str(path))
        mapping[row['source_path']] = path
    registry = json.loads((root/'confirmation.json').read_text())
    report = {'schema': 1, 'source_directory': inventory['source_directory'],
              'captured_utc': inventory['captured_utc'], 'verified_files': len(mapping),
              'production_approved': False, 'systems': [], 'return_checks': []}
    for system in registry['systems']:
        config = json.loads(mapping[system['configuration']].read_text())
        manifest = json.loads(mapping[system['manifest']].read_text())
        definition = config['probe_definition']
        ne = manifest['electron_count']
        active, count = definition['active_electrons'], definition['active_orbitals']
        inactive = (ne-active)//2
        chain = {'title': system['title'], 'active_space': [active, count],
                 'inactive': inactive, 'occupied': inactive+count,
                 'active_labels': [f'{i}.1' for i in range(inactive+1, inactive+count+1)],
                 'state_model': definition['state_model'], 'jobs': []}
        completed = {}
        for job in manifest['jobs']:
            directory = root/system['title']/job['name']
            output = directory/(job['name']+'.out')
            row = {'name': job['name'], 'status': job['status'], 'job_id': job.get('job_id'),
                   'predecessor': job.get('orbital_predecessor'),
                   'input_sha256': job['input_sha256'],
                   'wfu_sha256': job.get('wfu_sha256'),
                   'wfu': job.get('wfu'),
                   'orbital_xml': job.get('orbital_xml'),
                   'native_output': str(output.relative_to(root)),
                   'output_sha256': file_sha256(output) if output.exists() else None}
            chain['jobs'].append(row)
            if not output.exists():
                row['waiting_for'] = job.get('orbital_predecessor')
                continue
            text = output.read_text(errors='replace')
            native = (directory/(job['name']+'.native.inp')).read_text()
            execution = json.loads((directory/'execution.json').read_text())
            row.update(execution_status=execution['status'], returncode=execution.get('returncode'),
                       convergence=convergence_diagnostics(text, job),
                       normal_termination=bool(re.search(r'Molpro calculation terminated', text)),
                       completion_marker=bool(re.search(r'(?m)^\s*RECONDONE\s*=\s*1', text)),
                       version=re.findall(r'Version\s+[^\n]+', text)[:1],
                       native_alerts=[{'line': i, 'text': line.strip()} for i, line in enumerate(text.splitlines(), 1)
                           if re.search(r'^\s*\?\s|WARNING|NO CONVERGENCE|NOT CONVERGED|GLOBAL ERROR|negative norm', line, re.I)])
            multi = re.split(r'PROGRAM\s*\*\s*MULTI\b', text, flags=re.I)[-1]
            multi = re.split(r'PROGRAM\s*\*', multi, maxsplit=1, flags=re.I)[0]
            row['natural_active_table'] = [o for o in orbital_table(multi, 'NATURAL ORBITALS', 'Natural orbital dump')
                                           if o['orbital'] in chain['active_labels']]
            row['canonical_active_table'] = [o for o in orbital_table(multi, 'PSEUDO CANONICAL ORBITALS', 'Canonical orbital dump')
                                             if o['orbital'] in chain['active_labels']]
            row['native_counts'] = {key: re.findall(pattern, multi) for key, pattern in {
                'active_orbitals': r'Number of active\s+orbitals:\s+(\d+)',
                'active_electrons': r'Number of active electrons:\s+(\d+)',
                'closed': r'Number of closed-shell orbitals:\s+(\d+)'}.items()}
            row['native_states'] = job['root_states']
            row['diabatic_displaced_overlap_diagonal'] = re.findall(
                r'Diagonal elements of S\(R\+DR\|R\)[^:]+:\s*([^\n]+)', multi)
            # RESTART prints inherited variables before the actual calculation.
            # Those are not observations at this geometry.
            result_text = re.split(r'PROGRAM\s*\*\s*MULTI\b', text, flags=re.I)[-1]
            row['casscf_hartree'] = [float(x) for x in re.findall(r'(?m)^\s*RECONR\d+\s*=\s*('+NUMBER+')', result_text)]
            row['correlated_hartree'] = [float(x) for x in re.findall(r'(?m)^\s*RECONC\d+\s*=\s*('+NUMBER+')', result_text)]
            row['correlated_labels'] = job.get('correlated_energies_hartree', [])
            row['mrci_energy_checks'] = []
            for block in re.split(r'PROGRAM\s*\*\s*MRCI-F12\b', text, flags=re.I)[1:]:
                q = re.findall(r'Cluster corrected energies\s+('+NUMBER+r')\s*\(Davidson, relaxed reference\)', block)
                singles = re.findall(r'- singles correction\s+for state\s+\d+\s+('+NUMBER+')', block)
                selected = re.findall(r'(?m)^\s*SETTING KBVRCSEL\s*=\s*('+NUMBER+')', block)
                if q and singles and selected:
                    residual = float(selected[0])-float(q[0])-float(singles[0])
                    row['mrci_energy_checks'].append({'relaxed_davidson_before_cabs': float(q[0]),
                        'cabs_correction': float(singles[0]), 'selected_energy': float(selected[0]),
                        'rounding_residual_hartree': residual})
            if job['status'] == 'complete':
                if (not row['normal_termination'] or not row['completion_marker'] or
                    row['convergence']['errors'] or execution['status'] != 'complete' or execution['returncode'] != 0):
                    raise ValueError('Completed job failed independent review: '+str(output))
                if row['output_sha256'] != job['output_sha256'] or file_sha256(mapping[job['wfu']]) != job['wfu_sha256']:
                    raise ValueError('Completed artifact hash mismatch')
                if len(row['casscf_hartree']) != len(job['root_states']) or len(row['correlated_hartree']) != len(job['correlated_energies_hartree']):
                    raise ValueError('Root count mismatch')
                if any(abs(a-b['energy_hartree']) > 1.e-9 for a,b in zip(row['correlated_hartree'],job['correlated_energies_hartree'])):
                    raise ValueError('Independent energy extraction differs from manifest')
                if 'mrci-f12' in native:
                    if len(row['mrci_energy_checks']) != len(row['correlated_hartree']) or any(abs(v['rounding_residual_hartree'])>2.e-8 for v in row['mrci_energy_checks']):
                        raise ValueError('Relaxed-Davidson plus CABS selector differs')
                completed[job['anchor']['name'] if isinstance(job['anchor'],dict) else job['anchor']] = (row, directory)
        if 'center' in completed and 'return' in completed:
            (a, ad), (b, bd) = completed['center'], completed['return']
            comparison = same_geometry_overlap(ad/(a['name']+'.natural.orbitals.xml'),
                         bd/(b['name']+'.natural.orbitals.xml'), inactive, count)
            comparison.update(chain=system['title'],
                max_casscf_difference_kcal=max(abs(x-y)*KCAL for x,y in zip(a['casscf_hartree'], b['casscf_hartree'])),
                max_correlated_difference_kcal=max(abs(x-y)*KCAL for x,y in zip(a['correlated_hartree'], b['correlated_hartree'])))
            report['return_checks'].append(comparison)
        report['systems'].append(chain)
    report['status_counts'] = dict(Counter(j['status'] for c in report['systems'] for j in c['jobs']))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('extracted_export', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    atomic_json(args.output, review(args.extracted_export))
    print(args.output)
