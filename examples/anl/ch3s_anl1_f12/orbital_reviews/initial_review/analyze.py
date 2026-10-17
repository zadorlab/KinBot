"""Reproduce the CH3S review from its supplied archive; no QC or submissions.

Usage: python analyze.py ARCHIVE OUTPUT_DIRECTORY
Requires numpy and matplotlib for the diagnostic figure. Read archive members
directly rather than extracting them. CI weights below are lower bounds from
the printed coefficients, not complete CI norms or production error estimates.
"""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import tarfile
import xml.etree.ElementTree as ET

import numpy as np


def analyze(archive_path):
    with tarfile.open(archive_path) as archive:
        # Read each compressed member once; repeated gzip seeks are quadratic.
        members = {member.name: archive.extractfile(member).read()
                   for member in archive if member.isfile()}
        manifest_bytes = members['01_recon/reconnaissance.json']
        manifest = json.loads(manifest_bytes)
        remote = PurePosixPath(manifest['jobs'][0]['directory']).parent.parent

        def read(path):
            return members[str(PurePosixPath(path).relative_to(remote))]

        verified = 0
        profiles = []
        for job in manifest['jobs']:
            directory = PurePosixPath(job['directory'])
            output = read(directory / (job['name'] + '.out'))
            checks = [(directory / (job['name'] + '.inp'), job['input_sha256']),
                      (directory / (job['name'] + '.out'), job['output_sha256'])]
            checks.extend((entry['path'], entry['sha256'])
                          for entry in job.get('orbital_xml', {}).values())
            execution = json.loads(read(directory / 'execution.json'))
            if execution.get('native_input'):
                checks.append((execution['native_input'], execution['native_input_sha256']))
            for path, expected in checks:
                if hashlib.sha256(read(path)).hexdigest() != expected:
                    raise ValueError(f'Archive hash mismatch: {path}')
                verified += 1
            text = output.decode()
            natural = text.rsplit('NATURAL ORBITALS (STATE AVERAGED)', 1)[1]
            table = natural.split('Natural orbital dump', 1)[0].strip()
            blocks = re.split(r'CI Coefficients of symmetry \d+ \((Doublet|Quartet)\)', natural)
            omitted = {}
            for spin, block in zip(blocks[1::2], blocks[2::2]):
                rows = re.findall(r'^\s*([20ab]{8})\s+([-+\d.]+)\s+([-+\d.]+)\s+([-+\d.]+)\s*$',
                                  block.split('Energy:', 1)[0], re.M)
                for cas, freeze in [('5,4', 4), ('9,6', 2)]:
                    omitted[f'{cas}:{spin}'] = [sum(float(row[root + 1]) ** 2 for row in rows
                        if row[0][:freeze] != '2' * freeze) for root in range(3)]
            xml = ET.fromstring(read(job['orbital_xml']['natural']['path']))
            atoms = xml.findall('.//{http://www.xml-cml.org/schema}atom')
            geometry = [[float(atom.attrib[axis]) for axis in ('x3', 'y3', 'z3')] for atom in atoms]
            profiles.append({
                'name': job['name'], 'distance_angstrom': job['distance_angstrom'],
                'orientation': job.get('orientation'), 'predecessor': job.get('orbital_predecessor'),
                'states': job['state_energies_hartree'], 'occupations': job['natural_occupations'],
                'geometry_angstrom': geometry, 'natural_orbital_table': table,
                'printed_ci_weight_excluded_by_freezing': omitted,
                'source_wfu': job['wfu'], 'source_wfu_sha256': job['wfu_sha256'],
            })
    equivalences = []
    for index, first in enumerate(profiles):
        a = np.array(first['geometry_angstrom'])
        da = np.linalg.norm(a[:, None] - a[None, :], axis=-1)
        for second in profiles[index + 1:]:
            b = np.array(second['geometry_angstrom'])
            db = np.linalg.norm(b[:, None] - b[None, :], axis=-1)
            error = float(np.max(np.abs(da - db)))
            if error < 1.e-7:
                equivalences.append({'first': first['name'], 'second': second['name'],
                    'maximum_pair_distance_difference_angstrom': error,
                    'root_energy_differences_hartree': [right['energy_hartree'] - left['energy_hartree']
                        for left, right in zip(first['states'], second['states'])]})
    occupations = np.array([p['occupations'] for p in profiles])
    return {
        'archive_sha256': hashlib.sha256(Path(archive_path).read_bytes()).hexdigest(),
        'manifest_sha256': hashlib.sha256(manifest_bytes).hexdigest(),
        'verified_file_hashes': verified,
        'status': manifest['status'], 'coverage': manifest['coverage'],
        'occupation_rank_ranges': [{'index_one_based': i + 1, 'minimum': float(column.min()),
                                   'maximum': float(column.max())}
                                  for i, column in enumerate(occupations.T)],
        'isometric_geometry_comparisons': equivalences, 'profiles': profiles,
        'limitations': ['Occupation rank is not a tracked orbital identity.',
                        'Printed CI weights omit unprinted coefficients and do not bound energy errors.',
                        'WFU bytes are not in this archive; saved WFU hashes cannot be reverified here.',
                        'No CASSCF optimization or correlated validation is performed by this analysis.'],
    }


def plot(report, destination):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    radial = sorted((p for p in report['profiles'] if not p['orientation']),
                    key=lambda p: p['distance_angstrom'])
    angular = [p for p in report['profiles'] if p['orientation']]
    fresh = next(p for p in angular if p['distance_angstrom'] == 30 and p['orientation']['axis'] == 'bond')
    offset = fresh['states'][0]['energy_hartree']
    hartree_to_kcal = 627.509474
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)
    radii = [p['distance_angstrom'] for p in radial]
    for i, label in [(8, 'NO rank 9'), (9, 'NO rank 10')]:
        line, = axes[0].plot(radii, [p['occupations'][i] for p in radial], label=label + ', radial')
        axes[0].scatter([p['distance_angstrom'] for p in angular],
                        [p['occupations'][i] for p in angular], color=line.get_color(), marker='x')
    axes[0].set(xlim=(1.5, 7), ylim=(1.8, 2.01), xlabel='C–S distance (Å)', ylabel='Occupation',
                title='Extra occupied orbitals matter off-axis')
    axes[0].legend(fontsize=8)
    axes[0].text(.36, .22, '× = independent orientation probe\nRanks do not track orbital identity',
                 transform=axes[0].transAxes, fontsize=8)
    for i in range(3):
        axes[1].plot(radii, [(p['states'][i]['energy_hartree'] - offset) * hartree_to_kcal for p in radial],
                     marker='.', label=f'Radial doublet {i + 1}')
        axes[1].scatter([30], [(fresh['states'][i]['energy_hartree'] - offset) * hartree_to_kcal],
                        marker='x', color='black', zorder=5)
    axes[1].set(xlim=(13, 85), ylim=(-1, 40), xlabel='C–S distance (Å)',
                ylabel='Energy − fresh 30 Å reference (kcal/mol)', title='Propagated CASCI retains a different subspace')
    axes[1].legend(fontsize=8)
    axes[1].text(.04, .08, '× = fresh 30 Å roots\nQuartets coincide with doublets at long range',
                 transform=axes[1].transAxes, fontsize=8)
    worst = next(p for p in angular if p['name'].endswith('perpendicular_2_a60_r1p776838825'))
    for x, (cas, color) in enumerate([('5,4', 'tab:red'), ('9,6', 'tab:blue')]):
        axes[2].bar(np.arange(3) + x * .35, worst['printed_ci_weight_excluded_by_freezing'][cas + ':Quartet'],
                    width=.35, label=f'CAS({cas})', color=color)
    axes[2].set(xticks=np.arange(3) + .175, xticklabels=['Quartet 1', 'Quartet 2', 'Quartet 3'],
                ylim=(0, 1.02), ylabel='Excluded printed CI weight (lower bound)',
                title='Tilted 1.777 Å probe: fixed-orbital truncation')
    axes[2].legend(fontsize=8)
    fig.suptitle('CH3S / avtz-f12: diagnostic CASCI evidence, not production potentials', fontsize=13)
    fig.savefig(destination, dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    report = analyze(args.archive)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    plot(report, args.output / 'diagnostics.png')
    print(f'Analyzed {len(report["profiles"])} probes; verified {report["verified_file_hashes"]} file hashes.')
    print(args.output.resolve())
