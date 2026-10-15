"""Prepare finite example tests through the general orbital-review machinery.

No submission, automatic CAS selection, cross-basis checkpoint reuse, or
modification of old evidence. Chemical choices live in the JSON recipes.
Each output campaign shares one job cap and is run with --pilot-only.
"""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import tempfile

from kinbot import orbital_campaign as campaign, orbital_review as review
from kinbot import molpro_reconnaissance as recon
from kinbot.execution import atomic_json
from prepare_validation import plan as initial_plan


def limits(item, recipe):
    item['casscf_options'] = {'max_iterations': recipe['casscf_max_iterations']}
    item['correlation_options'] = {
        'max_reference_iterations': recipe['reference_max_iterations'],
        'max_iterations': recipe['correlated_max_iterations']}
    return item


def checks(electrons, orbitals, *, diagnostic=False, mrci=False):
    method = f'caspt2-f12({electrons},{orbitals})'
    values = [{'method': method, 'shift': .3, 'ipea': 0.},
              {'method': method, 'shift': .3, 'ipea': 0., 'reference_policy': 'state_average'}]
    if diagnostic:
        values = [
            {'method': method, 'shift': .3, 'ipea': .25}, values[0],
            {'method': method, 'shift': .1, 'ipea': 0.},
            {'method': method, 'shift': .3, 'ipea': .25, 'reference_policy': 'state_average'},
            {'method': method, 'shift': .3, 'ipea': 0., 'reference_policy': 'state_average'}]
    if mrci:
        values.append({'method': f'mrci-f12+q({electrons},{orbitals})',
                       'reference_policy': 'state_average'})
    return values


def reference_items(original, cases, recipe):
    """Independent-method, state-treatment, and distance/orientation tests."""
    result = []
    for case in cases:
        for space in case['spaces']:
            group = f'{case["label"]}_cas{space["electrons"]}_{space["orbitals"]}'
            for point in ('center', 'far', 'tilted'):
                # Large MRCI references initially need only center and far.
                if point == 'tilted' and not space['radial_neighbors']:
                    continue
                item = deepcopy(next(i for i in original if i['name'] == group + '_' + point))
                item.update(name=item['name'] + '_methods', model=group + '_methods',
                            correlation_checks=checks(space['electrons'], space['orbitals'],
                                                      diagnostic=point == 'far', mrci=True))
                result.append(limits(item, recipe))
        small = case['spaces'][0]
        group = f'{case["label"]}_cas{small["electrons"]}_{small["orbitals"]}'
        for point, distance in (('far', recipe['long_distance_angstrom']),
                                ('tilted', recipe['far_distance_angstrom'])):
            item = deepcopy(next(i for i in original if i['name'] == group + '_' + point))
            item.update(name=group + '_' + point + '_separation_check',
                        model=group + '_separation_check', target_distance_angstrom=distance,
                        correlation_checks=checks(small['electrons'], small['orbitals'], diagnostic=True))
            result.append(limits(item, recipe))
    return result


def fresh_configuration(source_directory, name, *, basis=None, distance=None, cache=None):
    cache = cache if cache is not None else {}
    directory = Path(source_directory).expanduser().resolve()
    root = directory.parents[1]
    if root not in cache:
        cache[root] = review._review_rows(root)  # Validates native convergence and hashes.
    saved, rows = cache[root]
    prior = next((row for row in rows if row['directory'] == str(directory)), None)
    if prior is None or prior['status'] != 'complete' or prior['probe_method'] != 'casscf':
        raise ValueError(f'Need a completed CASSCF source: {directory}')
    record = next(s for s in saved['systems'] if s['title'] == prior['review']['title'])
    loaded = recon._load(record['configuration'])
    anchor = deepcopy(next(a for a in loaded['anchors'] if a['name'] == prior['anchor']))
    if distance is not None:
        anchor['distance_angstrom'] = distance
    anchor['name'] = name
    data = {k: deepcopy(v) for k, v in loaded.items()
            if not k.startswith('_') and k not in ('coordinates', 'orientation_probes', 'scan_policy', 'initial_records')}
    if data.get('probe_definition', {}).get('orbital_partition'):
        raise ValueError('New-basis references require an AVAS source, not a truncated checkpoint.')
    data.update(title=name, anchors=[anchor], bases=[basis] if basis else
                [b for b in data['bases'] if b['basis'] == prior['basis']],
                probe_method='casscf', orbital_initialization='independent', correlation_checks=[],
                scf_options={'solver': 'so_sci', 'max_iterations': 150},
                casscf_options={'solver': 'so_sci', 'energy_convergence': 1e-8,
                                'max_iterations': 150, 'quasi_newton': True})
    metadata = {**prior['review'], 'title': name, 'model': name,
                'variant': {'name': 'casscf_fresh', 'probe_method': 'casscf',
                            'orbital_initialization': 'independent'}}
    return data, metadata


def prepare_references(output, old, original, cases, recipe):
    """One campaign: unfinished old points + new independent tests + DZ seeds."""
    configurations, metadata, preserved = [], [], []
    cache = {}
    with tempfile.TemporaryDirectory(prefix='kinbot-followup-') as temporary:
        tmp = Path(temporary)

        def gather(root):
            saved = json.loads((root / 'campaign.json').read_text())
            configurations.extend(Path(s['configuration']) for s in saved['systems'])
            metadata.extend(json.loads((root / 'review_index.json').read_text())['systems'])

        try:
            review.prepare_incomplete(old, tmp / 'remaining',
                casscf_max_iterations=recipe['casscf_max_iterations'],
                reference_max_iterations=recipe['reference_max_iterations'])
        except ValueError as error:
            if str(error) != 'All validation pilots are complete; no new campaign needed.':
                raise
        else:
            gather(tmp / 'remaining')
            preserved = json.loads((tmp / 'remaining' / 'review_plan.json').read_text())['completed_preserved']
        atomic_json(tmp / 'partitions.json', {'schema': 1, 'models': reference_items(original, cases, recipe)})
        review.prepare_partitions(tmp / 'partitions.json', tmp / 'methods')
        gather(tmp / 'methods')

        def add(data, meta):
            path = tmp / (data['title'] + '.json')
            atomic_json(path, data)
            recon._load(path)
            configurations.append(path)
            metadata.append(meta)

        for case in cases:
            space = case['spaces'][0]
            group = f'{case["label"]}_cas{space["electrons"]}_{space["orbitals"]}'
            for point in ('center', 'far', 'tilted', 'left', 'right'):
                item = next(i for i in original if i['name'] == group + '_' + point)
                name = case['label'] + '_basis_' + point
                data, meta = fresh_configuration(item['source_directory'], name,
                    basis=recipe['basis'], distance=item.get('target_distance_angstrom'), cache=cache)
                meta.update(purpose='basis_reference', case_label=case['label'], candidate_point=point)
                add(data, meta)
        for control in recipe.get('fragment_controls', []):
            name = control['name'] + '_term_check'
            data, meta = fresh_configuration(original[0]['source_directory'], name, cache=cache)
            data.pop('kinbot_input', None)
            data.update({k: control[k] for k in ('symbols', 'geometry_angstrom', 'charge', 'multiplicity')})
            data.update(anchors=[{'name': name}],
                probe_definition={
                    **{k: control[k] for k in ('active_electrons', 'active_orbitals', 'chemical_note')},
                    'state_model': {'state_average': control['state_average'],
                        'orbital_selection': {'mode': 'avas', 'centers': control['centers'],
                                             'threshold': .3, 'minao_basis': 'minao'},
                        'correlation_options': {'max_reference_iterations': recipe['reference_max_iterations'],
                                                'max_iterations': recipe['correlated_max_iterations']}}},
                correlation_checks=checks(control['active_electrons'], control['active_orbitals'],
                                          diagnostic=True, mrci=True))
            meta.update(channel=name, bond=None, region='diagnostic', chemical_note=control['chemical_note'],
                        rotd_target={'status': 'isolated_fragment_control'}, purpose='fragment_control')
            add(data, meta)
        result = campaign.prepare(output, configurations)
    atomic_json(result.parent / 'review_index.json', {'schema': 1, 'systems': metadata})
    atomic_json(result.parent / 'review_plan.json', {'schema': 1, 'stage': 'references',
        'old_campaign': str(Path(old).resolve()), 'completed_preserved': preserved, 'recipe': recipe})
    return result


def prepare_potentials(output, source, original, cases, recipe):
    """Matched DZ potentials, TZ shift tests, and explicit inner-average tests."""
    _, rows = review._review_rows(source)
    items = []
    for case in cases:
        for space in case['spaces']:
            group = f'{case["label"]}_cas{space["electrons"]}_{space["orbitals"]}'
            for item in [i for i in original if i['model'] == group]:
                point = item['name'][len(group) + 1:]
                matches = [j for j in rows if j['status'] == 'complete' and j['phase'] == 'pilot'
                           and j['review'].get('purpose') == 'basis_reference'
                           and j['review']['case_label'] == case['label']
                           and j['review']['candidate_point'] == point]
                if len(matches) != 1:
                    raise ValueError(f'Need one completed basis reference: {case["label"]}/{point}.')
                dz = deepcopy(item)
                dz.update(name=item['name'] + '_dz', model=group + '_dz',
                          source_directory=matches[0]['directory'])
                dz.pop('target_distance_angstrom', None)
                dz['correlation_checks'] = [deepcopy(item['correlation_checks'][0]),
                                            *checks(space['electrons'], space['orbitals'])]
                if space == case['spaces'][0]:
                    high_spin = max(s['multiplicity'] for s in matches[0]['state_model']['state_average'])
                    dz['correlation_checks'].append({'method': 'uccsd(t)-f12b', 'multiplicity': high_spin})
                items.append(limits(dz, recipe))
                # Center/far/tilted IPEA=0 already accompany the MRCI checks.
                # Repeat only when a new CC anchor is needed at that geometry.
                if (space == case['spaces'][0] or point in ('left', 'right')
                        or (point == 'tilted' and not space['radial_neighbors'])):
                    tz = deepcopy(item)
                    tz.update(name=item['name'] + '_zero_ipea', model=group + '_zero_ipea',
                              correlation_checks=checks(space['electrons'], space['orbitals']))
                    if space == case['spaces'][0]:
                        tz['correlation_checks'].append({'method': 'uccsd(t)-f12b', 'multiplicity': high_spin})
                    items.append(limits(tz, recipe))
            if case['label'] in recipe['inner_average_cases'] and space['radial_neighbors']:
                for index, average in enumerate(recipe['inner_averages'], 1):
                    for point in ('center', 'far', 'tilted'):
                        item = deepcopy(next(i for i in original if i['name'] == group + '_' + point))
                        item.update(name=item['name'] + f'_average{index}', model=group + f'_average{index}',
                            region='inner', state_average=average,
                            correlation_checks=checks(space['electrons'], space['orbitals'], mrci=True))
                        items.append(limits(item, recipe))
    with tempfile.TemporaryDirectory(prefix='kinbot-followup-') as temporary:
        specification = Path(temporary) / 'partitions.json'
        atomic_json(specification, {'schema': 1, 'models': items})
        return review.prepare_partitions(specification, output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['references', 'potentials'])
    parser.add_argument('source', type=Path, help='Old validation campaign, or completed reference follow-up.')
    parser.add_argument('output', type=Path)
    parser.add_argument('--recipe', type=Path, default=Path(__file__).with_name('followup_recipe.json'))
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a fresh output campaign; prior evidence is preserved.')
    recipe = json.loads(args.recipe.read_text())
    path = args.recipe.resolve().parent / recipe['validation_recipe']
    definition = json.loads(path.read_text())
    original = initial_plan(definition, path.parent)['models']
    prepare = prepare_references if args.stage == 'references' else prepare_potentials
    result = prepare(args.output, args.source, original, definition['cases'], recipe)
    saved = json.loads(result.read_text())
    print(f'Prepared {len(saved["systems"])} independent pilots: {result}')
    print('Run with --pilot-only --probe-method casscf --max-jobs 20 --max-refinements 0')


if __name__ == '__main__':
    main()
