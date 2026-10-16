"""Audit completed simultaneous-MRCI chains from a partial native export.

Unlike collect exports, a partial tar has no producer inventory. Check every
archived file against the extracted bytes, then check native artifacts against
their job/execution hashes. Never follow the absolute paths from Blodgett.
This read-only audit does not approve an entire production workflow.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import tarfile

import numpy as np

from kinbot.execution import atomic_json, file_sha256
from kinbot.molpro_reconnaissance import convergence_diagnostics
from review_export import KCAL, NUMBER, floats, mrci_energy_checks, same_geometry_overlap


def require(condition, message):
    if not condition:
        raise ValueError(message)


def audit(root, archive):
    root = root.resolve()
    members = {}
    with tarfile.open(archive, 'r:*') as tar:
        for member in tar:
            name = Path(member.name)
            require(not name.is_absolute() and '..' not in name.parts, 'Unsafe archive path')
            require(member.isdir() or member.isfile(), 'Unsupported archive member')
            path = root/name
            require(path.resolve().is_relative_to(root), 'Extracted member escapes root')
            if not member.isfile():
                continue
            require(str(name) not in members, 'Duplicate archive member')
            digest = hashlib.sha256(tar.extractfile(member).read()).hexdigest()
            require(path.stat().st_size == member.size and file_sha256(path) == digest,
                    'Archive/extraction mismatch: ' + str(name))
            members[str(name)] = {'bytes': member.size, 'sha256': digest}

    def bound(path, expected=None):
        name = str(path.relative_to(root))
        require(name in members, 'Artifact missing from archive: ' + name)
        if expected is not None:
            require(members[name]['sha256'] == expected, 'Bound hash mismatch: ' + name)
        return members[name]['sha256']

    report = {'schema': 1, 'production_approved': False,
              'provenance': {'archive': str(archive), 'archive_sha256': file_sha256(archive),
                             'extracted_root': str(root), 'producer_inventory_present': False,
                             'verified_archive_files': len(members)},
              'chains': [], 'same_geometry_controls': []}
    comparable = []
    for rel in sorted(members):
        if Path(rel).name != 'reconnaissance.json':
            continue
        manifest_path = root/rel
        manifest = json.loads(manifest_path.read_text())
        chain = {'name': manifest_path.parent.name, 'manifest_sha256': bound(manifest_path), 'jobs': []}
        for job in manifest['jobs']:
            require(job['status'] == 'complete', 'Partial chain is not complete: ' + rel)
            directory = manifest_path.parent/job['name']
            require(directory.resolve().is_relative_to(root), 'Unsafe job name')
            output = directory/(job['name']+'.out')
            native_path = directory/(job['name']+'.native.inp')
            execution_path = directory/'execution.json'
            bound(execution_path)
            execution = json.loads(execution_path.read_text())
            require(execution['status'] == 'complete' and execution['returncode'] == 0,
                    'Execution did not complete')
            bound(directory/(job['name']+'.inp'), job['input_sha256'])
            require(execution['input_sha256'] == job['input_sha256'], 'Execution input mismatch')
            bound(native_path, execution['native_input_sha256'])
            bound(output, job['output_sha256'])
            bound(directory/Path(job['wfu']).name, job['wfu_sha256'])
            for orbital in job['orbital_xml'].values():
                bound(directory/Path(orbital['path']).name, orbital['sha256'])
            text, native = output.read_text(), native_path.read_text()
            result = re.split(r'PROGRAM\s*\*\s*MULTI\b', text, flags=re.I)[-1]
            alerts = [line.strip() for line in text.splitlines()
                      if re.search(r'^\s*\?\s|WARNING|NO CONVERGENCE|NOT CONVERGED|GLOBAL ERROR|negative norm', line, re.I)]
            require(not alerts and not convergence_diagnostics(text, job)['errors'], 'Native diagnostic found')
            require('Molpro calculation terminated' in text and
                    re.search(r'(?m)^\s*RECONDONE\s*=\s*1', result), 'Missing completion markers')
            energies = [float(x) for x in re.findall(r'(?m)^\s*RECONC\d+\s*=\s*('+NUMBER+')', result)]
            casscf = [float(x) for x in re.findall(r'(?m)^\s*RECONR\d+\s*=\s*('+NUMBER+')', result)]
            declared = job['correlated_energies_hartree']
            require(len(energies) == len(declared) == len(casscf) == len(job['root_states']), 'Root count mismatch')
            require(all(abs(e-d['energy_hartree']) < 1.e-9 for e, d in zip(energies, declared)), 'Energy mismatch')
            require(all(d['root_convention'] == 'correlated_variational_energy_order' for d in declared), 'Wrong root convention')
            checks = mrci_energy_checks(text, native)
            require(len(checks) == len(energies) and
                    all(abs(c['rounding_residual_hartree']) < 2.e-8 for c in checks), 'Davidson selector mismatch')
            commands = re.findall(r'\{(mrci-f12\b[^{}]*)\}', native, re.I)
            blocks = re.split(r'PROGRAM\s*\*\s*MRCI-F12\b', text, flags=re.I)[1:]
            require(len(blocks) == len(commands) == len(job['state_model']['state_average']), 'Native call count mismatch')
            calls = []
            offset = 0
            for command, block, state in zip(commands, blocks, job['state_model']['state_average']):
                count = state['roots']
                require(re.search(r'(?:,|;)noswap(?:,|;)', command, re.I), 'Missing energy-order policy')
                for label in ('optimized', 'reference'):
                    match = re.search(r'Number of '+label+r' states:\s*(\d+)\s+Roots:\s*([^\n]+)', block)
                    require(match and int(match[1]) == count and
                            [int(x) for x in match[2].split()] == list(range(1, count+1)), 'Native manifold mismatch')
                variational = [float(x) for x in re.findall(r'!CI-F12 STATE\s+\d+\.\d+\s+Energy\s+('+NUMBER+')', block)]
                require(len(variational) == count and all(a <= b+1.e-9 for a,b in zip(variational,variational[1:])), 'Unexpected eigenvalue ordering')
                section = block.split('Coefficients of fixed reference functions:', 1)[1].split('Coefficients of rotated reference functions:', 1)[0]
                coeff = [floats(m[1]) for m in re.finditer(r'(?m)^\s*\d+\s+([-+\d. ]+)$', section)]
                require(np.shape(coeff) == (count, count), 'Reference coefficient dimensions differ')
                # A truncated reference projection is not a full CI overlap.
                singular = np.linalg.svd(coeff, compute_uv=False)
                overlaps = [{'pair': [int(m[0]), int(m[1])], 'value': float(m[2])}
                            for m in re.findall(r'!MRCI overlap\s+<(\d+)\.\d+\|\|(\d+)\.\d+>\s*('+NUMBER+')', block)]
                calls.append({'multiplicity': state['multiplicity'], 'roots': count,
                              'variational_hartree': variational,
                              'relaxed_davidson_hartree': energies[offset:offset+count],
                              'fixed_reference_coefficients': coeff,
                              'reference_projection_singular_values': singular.tolist(),
                              'dominant_reference_indices': (np.argmax(abs(np.array(coeff)), axis=1)+1).tolist(),
                              'printed_state_overlaps': overlaps})
                offset += count
            active = job['probe_active_space']
            inactive = (manifest['electron_count']-active['active_electrons'])//2
            multi = re.split(r'PROGRAM\s*\*', result, maxsplit=1, flags=re.I)[0]
            for pattern, expected in [(r'Number of active\s+orbitals:\s*(\d+)', active['active_orbitals']),
                                      (r'Number of active electrons:\s*(\d+)', active['active_electrons']),
                                      (r'Number of closed-shell orbitals:\s*(\d+)', inactive)]:
                counts = [int(x) for x in re.findall(pattern, multi)]
                require(counts and all(x == expected for x in counts), 'Native active-space count mismatch')
            row = {'name': job['name'], 'anchor': job['anchor'], 'job_id': job['job_id'], 'native_output': str(output.relative_to(root)),
                   'output_sha256': bound(output), 'native_input_sha256': bound(native_path),
                   'wfu_sha256': job['wfu_sha256'], 'orbital_xml': job['orbital_xml'],
                   'elapsed_seconds': execution['elapsed_seconds'], 'native_alerts': alerts,
                   'active_space': [active['active_electrons'], active['active_orbitals']],
                   'active_labels': [f'{i}.1' for i in range(inactive+1,inactive+active['active_orbitals']+1)],
                   'casscf_hartree': casscf, 'correlated_hartree': energies,
                   'energy_selector_checks': checks, 'calls': calls}
            chain['jobs'].append(row)
            comparable.append((chain['name'], job, row, directory, inactive))
        far = next((r for r in chain['jobs'] if r['anchor'] == 'far'), None)
        if far:
            for row in chain['jobs']:
                for call, asymptote in zip(row['calls'], far['calls']):
                    call['lowest_interaction_kcal'] = (min(call['relaxed_davidson_hartree'])-
                                                      min(asymptote['relaxed_davidson_hartree']))*KCAL
        report['chains'].append(chain)
    for i, (an, a, ar, ad, inactive) in enumerate(comparable):
        for bn, b, br, bd, bi in comparable[i+1:]:
            if (a['symbols'], a['geometry_angstrom'], a['basis'], ar['active_space']) != (b['symbols'], b['geometry_angstrom'], b['basis'], br['active_space']):
                continue
            control = same_geometry_overlap(ad/(a['name']+'.natural.orbitals.xml'),
                                            bd/(b['name']+'.natural.orbitals.xml'), inactive, ar['active_space'][1])
            control.update(chains=[an,bn], points=[a['name'],b['name']],
                           max_casscf_difference_hartree=max(abs(x-y) for x,y in zip(ar['casscf_hartree'],br['casscf_hartree'])),
                           selected_energy_changes_kcal=[(y-x)*KCAL for x,y in zip(ar['correlated_hartree'],br['correlated_hartree'])])
            minima_a = [min(c['relaxed_davidson_hartree']) for c in ar['calls']]
            minima_b = [min(c['relaxed_davidson_hartree']) for c in br['calls']]
            control['lowest_spin_gap_changes_kcal'] = [
                {'multiplicities': [ar['calls'][j]['multiplicity'], ar['calls'][k]['multiplicity']],
                 'change': ((minima_b[j]-minima_b[k])-(minima_a[j]-minima_a[k]))*KCAL}
                for j in range(len(minima_a)) for k in range(j+1,len(minima_a))]
            report['same_geometry_controls'].append(control)
    report['completed_points'] = len(comparable)
    report['selected_root_energies_checked'] = sum(len(r['correlated_hartree']) for _,_,r,_,_ in comparable)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('extracted_root', type=Path)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.extracted_root, args.archive)
    atomic_json(args.output, result)
    print(json.dumps({k: result[k] for k in ('completed_points', 'selected_root_energies_checked', 'production_approved')}, indent=2))
