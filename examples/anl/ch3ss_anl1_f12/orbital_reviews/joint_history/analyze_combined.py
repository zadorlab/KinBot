"""Review the sulfur reference panel across the five October 15 campaigns.

This case analysis verifies exported evidence and compares converged attempts.
It does not select or install a production active space. WFU bytes remain on
Blodgett; printed CI projections are diagnostics, not energy-error bounds.
"""
import argparse
from collections import Counter
import json
from pathlib import Path, PurePosixPath
import re
import tarfile
import tempfile

import numpy as np

from analyze import CANDIDATES, KCAL, ci_tables, plot, sha256
from kinbot.molpro_reconnaissance import convergence_diagnostics


def analyze(root, archive):
    report = json.loads((root / 'combined_review.json').read_text())
    verified, native, manifests, execution_counts = {}, {}, {}, Counter()

    def local(remote_path, campaign):
        return root / PurePosixPath(campaign).name / PurePosixPath(remote_path).relative_to(campaign)

    def verify(path, expected):
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f'Artifact hash mismatch: {path}')
        verified[str(path.relative_to(root))] = expected

    for source in report['campaigns']:
        campaign_path = source['campaign']
        campaign = json.loads(local(campaign_path + '/campaign.json', campaign_path).read_text())
        for system in campaign['systems']:
            config_path = local(system['configuration'], campaign_path)
            verify(config_path, system['sha256'])
            config = json.loads(config_path.read_text())
            if config.get('kinbot_input'):
                verify(local(config['kinbot_input'], campaign_path), system['context_sha256'])
            for key in ('manifest', 'pilot_manifest'):
                manifest = json.loads(local(system[key], campaign_path).read_text())
                for job in manifest['jobs']:
                    directory = local(job['directory'], campaign_path)
                    verify(directory / (job['name'] + '.inp'), job['input_sha256'])
                    manifests[job['directory']] = job
                    output = directory / (job['name'] + '.out')
                    if not output.is_file():
                        if job['status'] == 'complete':
                            raise ValueError(f'Accepted output missing: {directory}')
                        execution_counts['unrun'] += 1
                        continue
                    execution = json.loads((directory / 'execution.json').read_text())
                    if (execution['input_sha256'] != job['input_sha256']
                            or execution['status'] != 'complete' or execution['returncode'] != 0):
                        raise ValueError(f'Unexpected execution status/input: {directory}')
                    verify(local(execution['native_input'], campaign_path), execution['native_input_sha256'])
                    text = output.read_text()
                    diagnostics = convergence_diagnostics(text, job)
                    native[job['directory']] = (text, diagnostics)
                    execution_counts['native_outputs'] += 1
                    if job['status'] != 'complete':
                        execution_counts['rejected'] += 1
                        continue
                    if diagnostics['errors']:
                        raise ValueError(f'Accepted result failed native convergence: {directory}')
                    execution_counts['accepted_' + job['probe_method']] += 1
                    verify(output, job['output_sha256'])
                    for evidence in job['orbital_xml'].values():
                        verify(local(evidence['path'], campaign_path), evidence['sha256'])
                    space = job['probe_active_space']
                    if (int(re.search(r'Number of active\s+orbitals:\s+(\d+)', text)[1]) != space['active_orbitals']
                            or int(re.search(r'Number of active electrons:\s+(\d+)', text)[1]) != space['active_electrons']
                            or abs(sum(job['natural_occupations']) - manifest['electron_count']) > .0002):
                        raise ValueError(f'Native active-space/occupation mismatch: {directory}')

    if (len(report['jobs']) != len(manifests)
            or {j['directory'] for j in report['jobs']} != set(manifests)):
        raise ValueError('Combined report must include every pilot and panel attempt exactly once.')
    accepted = []
    for job in report['jobs']:
        original = manifests[job['directory']]
        for key in ('status', 'input_sha256', 'output_sha256', 'wfu_sha256', 'state_model', 'probe_active_space'):
            if job.get(key) != original.get(key):
                raise ValueError(f'Combined report differs from manifest for {key}: {job["directory"]}')
        if job['status'] != 'complete' or job['probe_method'] != 'casscf':
            continue
        text, diagnostics = native[job['directory']]
        if job.get('correlation_checks'):
            raise ValueError('This reference-only analysis must not mix multiple correlated MULTI blocks.')
        full_precision = re.findall(r'!MCSCF STATE\s+(\d+)\.\d+\s+(\w+)\s+Energy\s+([-0-9.]+)', text)
        spin_names = {2: 'Doublet', 4: 'Quartet'}
        labels = [(str(s['root']), spin_names[s['multiplicity']]) for s in job['root_states']]
        if [(r, spin) for r, spin, energy in full_precision] != labels:
            raise ValueError(f'Native root labels differ from the requested average: {job["directory"]}')
        energies = [float(energy) for _, _, energy in full_precision]
        if max(abs(a - b) for a, b in zip(energies, job['root_energies_hartree'])) > 5.1e-9:
            raise ValueError('Native precise energies disagree with the rounded manifest.')
        natural = text.rsplit('NATURAL ORBITALS (STATE AVERAGED)', 1)[1]
        tables = ci_tables(natural, job['probe_active_space']['active_orbitals'])
        expected_roots = Counter(spin for _, spin in labels)
        if (set(tables) != set(expected_roots)
                or any(len(coefficients) != expected_roots[spin]
                       for spin, entries in tables.items() for _, coefficients in entries)):
            raise ValueError(f'Printed CI tables do not cover every requested root: {job["directory"]}')
        projections = {}
        for electrons, orbitals, freeze, discard in CANDIDATES[job['review']['model']]:
            if (job['probe_active_space']['active_electrons'] - 2 * freeze != electrons
                    or job['probe_active_space']['active_orbitals'] - freeze - discard != orbitals):
                raise ValueError('Candidate truncation has inconsistent electron/orbital counts.')
            projections[f'{electrons},{orbitals}'] = {
                spin: np.sum([np.square(coefficients) for occupation, coefficients in entries
                              if occupation[:freeze] != '2' * freeze
                              or discard and occupation[-discard:] != '0' * discard], axis=0).tolist()
                if any(occupation[:freeze] != '2' * freeze or discard and occupation[-discard:] != '0' * discard
                       for occupation, _ in entries) else [0.] * len(entries[0][1])
                for spin, entries in tables.items()}
        row = {key: job[key] for key in ('anchor', 'basis', 'campaign', 'phase', 'directory',
            'distance_angstrom', 'orientation', 'geometry_angstrom', 'review', 'natural_occupations',
            'root_states', 'probe_active_space', 'state_model', 'elapsed_seconds', 'wfu', 'wfu_sha256',
            'solver_controls', 'initialization')}
        row.update(root_energies_hartree=energies, convergence=diagnostics,
                   fixed_orbital_truncations=projections,
                   natural_orbital_and_ci_text=natural)
        accepted.append(row)

    # First validated attempt in archive order; never choose by lowest energy.
    representatives = {}
    for row in accepted:
        key = (row['review']['model'], row['review']['variant']['name'], row['anchor'], row['basis'])
        representatives.setdefault(key, row)
    summaries = []
    for model in CANDIDATES:
        rows = [r for r in accepted if r['review']['model'] == model]
        expected = {j['anchor'] for j in report['jobs'] if j['review']['model'] == model}
        coverage = {variant: sorted({r['anchor'] for r in rows if r['review']['variant']['name'] == variant})
                    for variant in ('casscf_fresh', 'casscf_propagated')}
        pairs = []
        for fresh in rows:
            if fresh['review']['variant']['name'] != 'casscf_fresh':
                continue
            for seeded in rows:
                if (seeded['review']['variant']['name'] != 'casscf_propagated'
                        or fresh['anchor'] != seeded['anchor'] or fresh['basis'] != seeded['basis']):
                    continue
                if (fresh['state_model'] != seeded['state_model']
                        or fresh['probe_active_space'] != seeded['probe_active_space']
                        or fresh['root_states'] != seeded['root_states']
                        or not np.allclose(fresh['geometry_angstrom'], seeded['geometry_angstrom'], rtol=0., atol=1e-12)):
                    raise ValueError('Initialization comparison changes more than the starting orbitals.')
                delta = np.asarray(seeded['root_energies_hartree']) - fresh['root_energies_hartree']
                pairs.append({'anchor': fresh['anchor'], 'fresh': fresh['directory'], 'seeded': seeded['directory'],
                              'differences_hartree': delta.tolist(),
                              'maximum_absolute_kcal_mol': float(np.max(np.abs(delta))) * KCAL,
                              'maximum_occupation_difference': float(np.max(np.abs(
                                  np.asarray(fresh['natural_occupations']) - seeded['natural_occupations'])))})
        occupations = np.asarray([r['natural_occupations'] for r in rows])
        truncations = {}
        for candidate in rows[0]['fixed_orbital_truncations']:
            weight, directory, spin, state = max(
                (weight, r['directory'], spin, i + 1) for r in rows
                for spin, values in r['fixed_orbital_truncations'][candidate].items()
                for i, weight in enumerate(values))
            truncations[candidate] = {'maximum_excluded_printed_weight': weight,
                                     'directory': directory, 'spin': spin, 'root': state}
        summaries.append({'model': model, 'expected_anchors': sorted(expected), 'coverage': coverage,
            'missing': {variant: sorted(expected - set(found)) for variant, found in coverage.items()},
            'comparisons': pairs,
            'maximum_initialization_difference_kcal_mol': max(p['maximum_absolute_kcal_mol'] for p in pairs),
            'maximum_occupation_difference': max(p['maximum_occupation_difference'] for p in pairs),
            'natural_occupation_ranges': [{'index_one_based': i + 1, 'minimum': float(column.min()),
                                          'maximum': float(column.max())} for i, column in enumerate(occupations.T)],
            'fixed_orbital_truncation_summary': truncations})
    return {'archive_sha256': sha256(archive), 'verified_hash_count': len(verified), 'verified_files': verified,
            'counts': dict(execution_counts), 'accepted_casscf_attempts': len(accepted),
            'unique_casscf_model_variant_anchor_basis': len(representatives),
            'models': summaries, 'jobs': list(representatives.values()),
            'all_accepted_casscf_attempts': accepted,
            'limitations': ['WFU bytes remain on Blodgett; only their recorded hashes are available locally.',
                'Agreement of converged starts at these anchors is not a proof for every sampling geometry.',
                'Nearly degenerate states/orbitals can rotate; equal NO spectra do not prove identical orbitals.',
                'Printed CI truncations omit small coefficients and do not bound correlated or kinetic errors.',
                'No candidate reoptimization or correlated basis/cost comparison has been run.',
                'The user must select the final active spaces and inner/outer state models.']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='kinbot-combined-analysis-') as temporary:
        root = Path(temporary)
        with tarfile.open(args.archive) as bundle:
            bundle.extractall(root, filter='data')
        for path in root.glob('*.tar.gz'):
            with tarfile.open(path) as bundle:
                bundle.extractall(root / path.name.removesuffix('.tar.gz'), filter='data')
        result = analyze(root, args.archive)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    plot(result, args.output / 'diagnostics.png', title='Converged TZ-F12 CASSCF reference panel\n'
         'Eight geometries per channel; repaired points included. No correlated energies.')
    print(json.dumps({key: result[key] for key in ('archive_sha256', 'verified_hash_count', 'counts',
          'accepted_casscf_attempts', 'unique_casscf_model_variant_anchor_basis')}, indent=2))
    for model in result['models']:
        print(model['model'], 'missing:', model['missing'],
              'maximum initialization difference:', model['maximum_initialization_difference_kcal_mol'], 'kcal/mol')
