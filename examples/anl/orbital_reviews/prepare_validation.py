"""Prepare a finite, editable example comparison; submit nothing.

All molecule/space/anchor choices are data in validation_candidates.json.
Methods, basis and radial stencil points come from the linked KinBot inputs.
The resulting single-geometry systems run once with --pilot-only.
"""
import argparse
import json
from pathlib import Path

from kinbot.molpro import vrc_multireference_energy_contract
from kinbot.parameters import Parameters
from kinbot.vrc_correction import resolve_vrc_level


def plan(recipe, base):
    summary = json.loads((base / recipe['summary']).read_text())
    models = []
    for case in recipe['cases']:
        parameters = Parameters(str((base / case['kinbot_input']).resolve())).par
        grid = sorted(set(float(d) for d in parameters['vrc_tst_scan_points']))

        def source(anchor):
            matches = [j for j in summary['jobs'] if j['review']['model'] == case['model']
                       and j['anchor'] == anchor
                       and j['review']['variant']['name'] == 'casscf_fresh']
            if len(matches) != 1 or matches[0]['convergence']['errors']:
                raise ValueError(f'Need one validated fresh seed: {case["model"]} {anchor}')
            return matches[0]

        radial = source(case['radial_seed'])
        # Interior requested scan point nearest the already reviewed seed.
        # These are atom-pair correction coordinates, not inferred COM pivots.
        if len(grid) < 3:
            raise ValueError('A radial stencil needs at least three requested scan distances.')
        center = min(range(1, len(grid) - 1), key=lambda i: abs(grid[i] - radial['distance_angstrom']))
        for space in case['spaces']:
            electrons, orbitals = space['electrons'], space['orbitals']
            name = f'{case["label"]}_cas{electrons}_{orbitals}'
            points = [('center', radial, grid[center]),
                      ('far', source(case['far']), None),
                      ('tilted', source(case['orientation_seed']), grid[center])]
            if space['radial_neighbors']:
                points += [('left', radial, grid[center-1]), ('right', radial, grid[center+1])]
            level = resolve_vrc_level({'method': parameters['vrc_tst_sample_method'],
                                       'basis': parameters['vrc_tst_sample_basis']},
                                      {'active_electrons': electrons, 'active_orbitals': orbitals})
            if not vrc_multireference_energy_contract(level['method']):
                raise ValueError('This example comparison needs a supported correlated MR sampling method.')
            for position, seed, distance in points:
                choice = space.get('partition_seeds', {}).get(position, {})
                if choice.get('anchor'):
                    seed = source(choice['anchor'])
                if level['basis'] != seed['basis']:
                    raise ValueError('Prepare basis-specific reviewed seeds before changing the sampling basis.')
                item = {'name': name + '_' + position, 'model': name,
                    'source_directory': seed['directory'], 'region': 'diagnostic',
                    'chemical_note': case['note'] + ' Reoptimize the reviewed natural-orbital partition; '
                        'test interaction accuracy, state character and continuity before selecting it.',
                    'active_electrons': electrons, 'active_orbitals': orbitals,
                    'close_lowest': space['close_lowest'], 'discard_highest': space.get('discard_highest', 0),
                    'correlation_checks': [{'method': level['method'],
                        'shift': parameters['vrc_tst_scan_shift'], 'ipea': parameters['vrc_tst_scan_ipea']}]}
                if 'orbital_swaps' in choice:
                    item['orbital_swaps'] = choice['orbital_swaps']
                if space.get('bind_reviewed_checkpoint'):
                    item['source_wfu_sha256'] = seed['wfu_sha256']
                if distance is not None:
                    item['target_distance_angstrom'] = distance
                models.append(item)
    return {'schema': 1, 'models': models}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--recipe', type=Path, default=Path(__file__).with_name('validation_candidates.json'))
    args = parser.parse_args()
    definition = plan(json.loads(args.recipe.read_text()), args.recipe.resolve().parent)
    with args.output.open('x') as stream:
        json.dump(definition, stream, indent=2)
        stream.write('\n')
    print(f'Prepared {len(definition["models"])} single-geometry comparisons: {args.output}')
