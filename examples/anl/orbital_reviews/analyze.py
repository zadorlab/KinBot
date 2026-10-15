"""Reproduce the October 15 sulfur orbital review without running chemistry.

Usage: python analyze.py ARCHIVE OUTPUT_DIRECTORY
The proposed truncations below are case-specific review hypotheses, not an
automatic active-space selector. Printed CI weights are incomplete and do
not bound correlated energy errors. All orbital numbers are local NO ranks.
"""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import tarfile
import tempfile

import numpy as np

from kinbot.molpro_reconnaissance import convergence_diagnostics


KCAL = 627.509474
CANDIDATES = {
    'system_1__cs_valence_reference': [(5, 4, 4, 0), (9, 6, 2, 0)],
    'system_2__cs_valence_reference': [(3, 3, 8, 1), (9, 7, 5, 0), (11, 8, 4, 0), (13, 9, 3, 0)],
    'system_2__ss_valence_reference': [(7, 5, 6, 1), (9, 7, 5, 0), (13, 9, 3, 0)],
}
LABELS = ['CH3S: C–S', 'CH3SS: C–S', 'CH3SS: S–S']


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def ci_tables(text, orbitals):
    parts = re.split(r'CI Coefficients of symmetry \d+ \((\w+)\)', text)
    tables = {}
    for spin, block in zip(parts[1::2], parts[2::2]):
        rows = []
        for line in block.split('Energy:', 1)[0].splitlines():
            if re.match(r'^\s*[20ab]{' + str(orbitals) + r'}\s+', line):
                fields = line.split()
                rows.append((fields[0], [float(x) for x in fields[1:]]))
        if rows:
            tables[spin] = rows
    return tables


def analyze(root, archive):
    reports = sorted(root.glob('review_report*.json'))
    report = json.loads(reports[-1].read_text())
    remote = PurePosixPath(report['campaign'])
    campaign = json.loads((root / 'campaign.json').read_text())
    checks = {}

    def local(path):
        return root / PurePosixPath(path).relative_to(remote)

    def verify(path, expected):
        path = Path(path)
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f'Artifact hash mismatch: {path}')
        checks[str(path.relative_to(root))] = expected

    native, diagnostics = {}, {}
    for system in campaign['systems']:
        source = local(system['configuration'])
        verify(source, system['sha256'])
        config = json.loads(source.read_text())
        verify(local(config['kinbot_input']), system['context_sha256'])
        for key in ('pilot_manifest', 'manifest'):
            manifest = json.loads(local(system[key]).read_text())
            for job in manifest['jobs']:
                assert job['status'] == 'complete', job['directory']
                directory = local(job['directory'])
                verify(directory / (job['name'] + '.inp'), job['input_sha256'])
                output = directory / (job['name'] + '.out')
                verify(output, job['output_sha256'])
                for evidence in job['orbital_xml'].values():
                    verify(local(evidence['path']), evidence['sha256'])
                execution = json.loads((directory / 'execution.json').read_text())
                assert execution['input_sha256'] == job['input_sha256']
                assert execution['returncode'] == 0 and execution['status'] == 'complete'
                verify(local(execution['native_input']), execution['native_input_sha256'])
                text = output.read_text()
                space = job['probe_active_space']
                assert int(re.search(r'Number of active\s+orbitals:\s+(\d+)', text)[1]) == space['active_orbitals']
                assert int(re.search(r'Number of active electrons:\s+(\d+)', text)[1]) == space['active_electrons']
                assert abs(sum(job['natural_occupations']) - manifest['electron_count']) < .0002
                diagnostics[job['directory']] = convergence_diagnostics(text, job)
                native[job['directory']] = text

    rows = []
    for job in report['jobs']:
        verify(local(job['directory']) / (job['name'] + '.out'), job['output_sha256'])
        text = native[job['directory']]
        natural = text.rsplit('NATURAL ORBITALS (STATE AVERAGED)', 1)[1]
        tables = ci_tables(natural, job['probe_active_space']['active_orbitals'])
        truncations = {}
        for electrons, orbitals, freeze, discard in CANDIDATES[job['review']['model']]:
            assert job['probe_active_space']['active_electrons'] - 2 * freeze == electrons
            assert job['probe_active_space']['active_orbitals'] - freeze - discard == orbitals
            excluded = {}
            for spin, entries in tables.items():
                weights = np.zeros(len(entries[0][1]))
                for configuration, coefficients in entries:
                    if configuration[:freeze] != '2' * freeze or (
                            discard and configuration[-discard:] != '0' * discard):
                        weights += np.square(coefficients)
                excluded[spin] = weights.tolist()
            truncations[f'{electrons},{orbitals}'] = {
                'freeze_first_active_natural_orbitals': freeze,
                'discard_last_active_natural_orbitals': discard,
                'excluded_printed_squared_coefficients': excluded}
        rows.append({key: job[key] for key in (
            'anchor', 'basis', 'directory', 'distance_angstrom', 'geometry_angstrom',
            'orientation', 'probe_method', 'review', 'natural_occupations',
            'root_energies_hartree', 'root_states', 'probe_active_space',
            'elapsed_seconds', 'wfu', 'wfu_sha256')})
        rows[-1].update(natural_orbital_table=natural.split('Natural orbital dump', 1)[0].strip(),
                        fixed_orbital_truncations=truncations,
                        convergence=diagnostics[job['directory']])

    grouped = defaultdict(list)
    for row in rows:
        grouped[row['review']['model']].append(row)
    summaries = []
    for model, jobs in grouped.items():
        fresh = [j for j in jobs if j['review']['variant']['name'] == 'casscf_fresh']
        seeded = {j['anchor']: j for j in jobs if j['review']['variant']['name'] == 'casscf_propagated'}
        pairs = [{'anchor': j['anchor'], 'both_converged': not (
            j['convergence']['errors'] or seeded[j['anchor']]['convergence']['errors']),
            'state_differences_hartree': (
            np.array(seeded[j['anchor']]['root_energies_hartree']) - j['root_energies_hartree']).tolist()}
            for j in fresh]
        far = max((j for j in fresh if j['orientation'] is None), key=lambda j: j['distance_angstrom'])
        occupations = np.asarray([j['natural_occupations'] for j in fresh])
        truncation_summary = {}
        for space in fresh[0]['fixed_orbital_truncations']:
            values = [(weight, j['anchor'], spin, i + 1) for j in fresh
                for spin, weights in j['fixed_orbital_truncations'][space]['excluded_printed_squared_coefficients'].items()
                for i, weight in enumerate(weights)]
            worst = max(values)
            truncation_summary[space] = {'maximum_excluded_printed_weight': worst[0],
                'anchor': worst[1], 'spin': worst[2], 'root': worst[3]}
        timing = {}
        for variant in sorted({j['review']['variant']['name'] for j in jobs}):
            times = [j['elapsed_seconds'] for j in jobs if j['review']['variant']['name'] == variant]
            timing[variant] = {'minimum_seconds': min(times), 'median_seconds': float(np.median(times)),
                               'maximum_seconds': max(times), 'total_seconds': sum(times)}
        summaries.append({'model': model, 'comparisons': pairs,
            'maximum_raw_initialization_difference_kcal_mol': max(abs(v) for p in pairs
                for v in p['state_differences_hartree']) * KCAL,
            'converged_pairs': sum(p['both_converged'] for p in pairs),
            'maximum_converged_initialization_difference_kcal_mol': max((abs(v) for p in pairs
                if p['both_converged'] for v in p['state_differences_hartree']), default=0.) * KCAL,
            'natural_occupation_ranges': [{'index_one_based': i + 1, 'minimum': float(column.min()),
                'maximum': float(column.max())} for i, column in enumerate(occupations.T)],
            'far_anchor': far['anchor'], 'far_state_splittings_kcal_mol': (
                (np.asarray(far['root_energies_hartree']) - far['root_energies_hartree'][0]) * KCAL).tolist(),
            'fixed_orbital_truncation_summary': truncation_summary, 'timing': timing})
    return {'archive_sha256': sha256(archive), 'verified_files': checks,
        'verified_hash_count': len(checks), 'native_jobs_checked_including_pilots': len(native),
        'rejected_native_jobs': {path: d for path, d in diagnostics.items() if d['errors']},
        'native_warnings': {path: d['warnings'] for path, d in diagnostics.items() if d['warnings']},
        'optimized_jobs_checked_including_pilots': sum('CONVERGENCE REACHED!' in t for t in native.values()),
        'full_panel_jobs': len(rows), 'models': summaries, 'jobs': rows,
        'limitations': ['WFU bytes remain on Blodgett; their recorded hashes cannot be reverified locally.',
            'Initialization agreement is established only at these geometries, basis and state averages.',
            'Natural orbital rank is not an orbital identity across geometries.',
            'Printed CI coefficients omit small terms; zero excluded printed weight is not zero truncation error.',
            'Orbital occupations and fixed-orbital projections do not bound correlated energy or rate errors.',
            'No correlated checks were requested in this campaign.']}


def plot(report, destination):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for col, (model, label) in enumerate(zip(CANDIDATES, LABELS)):
        rows = [j for j in report['jobs'] if j['review']['model'] == model
                and j['review']['variant']['name'] == 'casscf_fresh']
        active = rows[0]['probe_active_space']['active_orbitals']
        occupations = np.asarray([j['natural_occupations'][-active:] for j in rows])
        start = len(rows[0]['natural_occupations']) - active + 1
        ranks = np.arange(start, start + active)
        ax = axes[0, col]
        ax.vlines(ranks, occupations.min(axis=0), occupations.max(axis=0), color='#2a6f97', lw=5)
        far = max((j for j in rows if j['orientation'] is None), key=lambda j: j['distance_angstrom'])
        ax.scatter(ranks, far['natural_occupations'][-active:], color='#d95f02', s=25, zorder=3, label='30 Å')
        ax.set(title=label, xlabel='Natural orbital rank (local ordering)', ylabel='Occupation',
               xticks=ranks, ylim=(-.05, 2.08))
        ax.grid(axis='y', alpha=.2)
        ax.legend(fontsize=9)
        ax = axes[1, col]
        radial = sorted((j for j in rows if j['orientation'] is None and j['distance_angstrom'] < 10),
                        key=lambda j: j['distance_angstrom'])
        for i, state in enumerate(far['root_states']):
            spin = state['multiplicity']
            ax.plot([j['distance_angstrom'] for j in radial],
                [(j['root_energies_hartree'][i] - far['root_energies_hartree'][i]) * KCAL for j in radial],
                marker='o', ls='-' if spin == 2 else '--',
                color=plt.get_cmap('tab10')((state['root'] - 1) % 10),
                label=('D' if spin == 2 else 'Q') + str(state['root']))
            bad = [j for j in radial if j['convergence']['errors']]
            ax.scatter([j['distance_angstrom'] for j in bad],
                [(j['root_energies_hartree'][i] - far['root_energies_hartree'][i]) * KCAL for j in bad],
                marker='x', s=65, color='black', zorder=4)
        ax.axhline(0, color='gray', lw=.8)
        ax.set(xlabel='Reactive atom distance (Å)', ylabel='Energy − own 30 Å root (kcal/mol)')
        ax.legend(fontsize=8, ncol=2 if len(far['root_states']) <= 6 else 4)
        ax.grid(alpha=.2)
    fig.suptitle('TZ-F12 CASSCF diagnostic references: occupations and sparse radial energies\n'
                 'Black crosses: iteration-limit failures. Ranges include provisional points. No correlated energies.',
                 fontsize=12)
    fig.savefig(destination, dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='kinbot-review-analysis-') as temporary:
        root = Path(temporary)
        with tarfile.open(args.archive) as bundle:
            bundle.extractall(root, filter='data')
        result = analyze(root, args.archive)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    plot(result, args.output / 'diagnostics.png')
    print(json.dumps({key: result[key] for key in (
        'verified_hash_count', 'native_jobs_checked_including_pilots',
        'optimized_jobs_checked_including_pilots', 'full_panel_jobs')}, indent=2))
