"""Read the completed detached panels; write evidence, never a model approval.

python audit.py EXTRACTED_DETACHED EXTRACTED_ROOT_SPECTRUM --output audit.json
All MO comparisons use identical geometries and full XML coefficient matrices.
"""
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np

from kinbot.molpro_diagnostics import parse_native_output, _orbital_tables


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(data)
    return h.hexdigest()


def checked(path, expected):
    actual = digest(path)
    if actual != expected:
        raise ValueError(f'Changed evidence: {path}')
    return actual


def orbital_excerpt(orbital):
    """Keep a readable amplitude excerpt; full native tables stay hash-bound."""
    result = {k: v for k, v in orbital.items() if k != 'ao_coefficients'}
    coefficients = sorted(orbital['ao_coefficients'],
                          key=lambda c: abs(c['coefficient']), reverse=True)
    result['largest_printed_AO_amplitudes'] = [
        f"{c['center']}:{c['function']}{c['component']} {c['coefficient']:+.5f}"
        for c in coefficients[:12]]
    result['other_printed_amplitudes'] = max(0, len(coefficients) - 12)
    return result


def main(detached, spectrum):
    module_path = (Path(__file__).resolve().parent.parent /
                   'r7_compact_union_20261025/compare_orbitals.py')
    spec = importlib.util.spec_from_file_location('same_geometry_mos', module_path)
    mo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mo)
    pool = detached / 'runs/ch3s_ch3ss_orbital_reviews_20261027'
    jobs = []
    for manifest_path in sorted(pool.glob('*/*/*/reconnaissance.json')):
        campaign = manifest_path.parent
        for entry in json.loads(manifest_path.read_text())['jobs']:
            name = entry['name']
            directory = campaign / name
            output = directory / (name + '.out')
            hashes = {'output': checked(output, entry['output_sha256']),
                      'wfu': checked(directory / (name + '.wfu'), entry['wfu_sha256'])}
            for kind, item in entry['orbital_xml'].items():
                hashes[kind + '_xml'] = checked(directory / Path(item['path']).name,
                                                item['sha256'])
            native_input = directory / (name + '.native.inp')
            hashes['native_input'] = digest(native_input)
            text = output.read_text()
            native = parse_native_output(text, native_input.read_text())
            if not native['normal_termination_seen'] or native['native_error_lines']:
                raise ValueError(f'Native completion failed: {output}')
            sections = [s for s in native['native_sections'] if s['program'] == 'MULTI']
            if entry.get('probe_active_space'):
                if len(sections) != 1 or native['issues']:
                    raise ValueError(f'Native CAS validation: {output}: {native["issues"]}')
                section = sections[0]
                orbitals = section['orbital_tables'][-1]['orbitals']
                ci = section['ci_tables']
                # Compare energies parsed from the native text to the manifest,
                # not merely the scalar RECONENERGY (which can select another spin).
                for e in entry['state_energies_hartree']:
                    matches = [v for v in section['state_energies']
                               if int(v['state'].split('.')[0]) == e['root']
                               and (not v['spin_label'].strip() or
                                    v['spin_label'].strip().lower() ==
                                    {2: 'doublet', 4: 'quartet'}[e['multiplicity']])]
                    if len(matches) != 1 or abs(matches[0]['energy_Eh'] - e['energy_hartree']) > 2.e-8:
                        raise ValueError(f'Native/manifest state mismatch: {output}: {e}')
            else:
                lines = text.replace('ELECTRON ORBITALS', 'CANONICAL ORBITALS').splitlines()
                orbitals = _orbital_tables(lines, 1, None, None)[-1]['orbitals']
                ci = []
            relative = campaign.relative_to(pool)
            job = dict(panel=relative.parts[0], screen=relative.parts[1],
                       candidate=relative.parts[2], job=name,
                       model=entry.get('probe_active_space'), hashes=hashes,
                       states=entry.get('state_energies_hartree', []), ci=ci,
                       orbitals=[orbital_excerpt(o) for o in orbitals
                                 if 14 <= int(o['index'].split('.')[0]) <= 24],
                       source_manifest=str(manifest_path.relative_to(detached)))
            jobs.append(job)

    comparisons = []

    def compare(left, right, selections, label):
        value = mo.compare(Path('/'), left, right, selections)
        a, _, _ = mo.read(left)
        b, _, _ = mo.read(right)
        transform = b @ np.linalg.pinv(a)
        # Projection of each left orbital into the selected right subspace.
        # This is an MO overlap, not a sum of squared printed AO coefficients.
        value['left_projection_into_right'] = {
            key: {str(j): float(np.sum(transform[np.array(i) - 1, j - 1] ** 2))
                  for j in left_indices}
            for key, (i, left_indices) in selections.items()}
        def portable(path):
            for prefix, root in [('detached', detached), ('spectrum', spectrum)]:
                if path.is_relative_to(root):
                    return prefix + '/' + str(path.relative_to(root))
            raise ValueError(path)
        value['sha256'] = {portable(left): digest(left), portable(right): digest(right)}
        value.update(left=portable(left), right=portable(right), label=label)
        comparisons.append(value)

    support = pool / 'state_support'
    old = spectrum / ('production/orbital_recon/root_spectrum_20261026/'
                      'csh2_sh_exit/csh2_sh_exit_fixed_spectrum')
    old_manifest = {j['name']: j for j in json.loads((old / 'reconnaissance.json').read_text())['jobs']}
    for basis in ('adz', 'atz'):
        for anchor in ('exit_l2', 'exit_r_30p0'):
            name = basis + '__' + anchor
            source = old / name / (name + '.natural.orbitals.xml')
            checked(source, old_manifest[name]['orbital_xml']['natural']['sha256'])
            for candidate, first in [('csh2_sh_six_state_frontier', 19),
                                     ('csh2_sh_six_state_union', 17)]:
                target = support / 'csh2_sh_exit' / candidate / name / source.name
                compare(source, target, {
                    'SH_pair_in_active': (list(range(first, 23)), [19, 20]),
                    'frontier_in_active': (list(range(first, 23)), [19, 20, 21, 22]),
                }, candidate + '/' + name)

    pairs = [
        ('new_cs', 'new_cs_ch_frontier', 'new_cs_sulfur_donor', 20, 22),
        ('new_ss', 'new_ss_sh_pi', 'new_ss_carbon_sigma_union', 19, 22),
        ('new_hcs_h2s_exit', 'new_hcs_h2s_pi', 'new_hcs_h2s_donor', 20, 22),
        ('new_hcssh_h_exit', 'new_hcssh_h_somo', 'new_hcssh_h_s_donor', 21, 21),
        ('new_hcssh_h_exit', 'new_hcssh_h_somo', 'new_hcssh_h_c_donor', 21, 21),
    ]
    for screen, small, large, first, last in pairs:
        for directory in sorted((support / screen / small).glob('*__*')):
            if not directory.is_dir():
                continue
            name = directory.name
            compare(directory / (name + '.natural.orbitals.xml'),
                    support / screen / large / name / (name + '.natural.orbitals.xml'), {
                        'small_active_in_large_active':
                            (list(range(first - 1, last + 1)), list(range(first, last + 1))),
                        'closed_17_to_20_in_large_active':
                            (list(range(first - 1, last + 1)), list(range(17, 21))),
                        'closed_17_to_20_in_large_closed':
                            (list(range(1, first - 1)), list(range(17, 21))),
                    }, large + '/' + name)

    return dict(
        schema=1,
        scope='Native completion and orbital/state evidence; does not approve or install a model.',
        method='Full same-geometry MO transforms, native CI and energies. No cross-geometry overlap; '
               'printed AO amplitudes are not populations; occupation thresholds do not select a CAS. '
               'Orbital excerpts list at most 12 largest printed amplitudes; full tables remain in '
               'the hash-bound native outputs and full coefficients in the XML used for overlaps.',
        archive_sha256='e4c5bbe1dc5346b95db594f41d9d59319173bbc3fc8a0790892adaba2bcfe2a8',
        counts=dict(jobs=len(jobs), casscf=sum(bool(j['model']) for j in jobs),
                    scf=sum(not j['model'] for j in jobs),
                    by_panel=dict(Counter(j['panel'] for j in jobs))),
        jobs=jobs, comparisons=comparisons)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('detached', type=Path)
    parser.add_argument('spectrum', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = main(args.detached.resolve(), args.spectrum.resolve())
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['counts'], sort_keys=True))
