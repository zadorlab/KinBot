"""Prepare the requested five-point CAS(7,5) example addition; submit nothing.

All orbital/geometry choices are in cas75_candidates.json. Source checkpoints
are hash-bound to the reviewed native orbitals, not selected by occupations.
"""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import tempfile

from kinbot.execution import atomic_json
from kinbot import orbital_review
from prepare_validation import plan
from prepare_followup import checks, limits


def definition(recipe_path):
    recipe_path = Path(recipe_path).resolve()
    recipe = json.loads(recipe_path.read_text())
    controls = json.loads(recipe_path.with_name('followup_recipe.json').read_text())
    result = plan(recipe, recipe_path.parent)
    for item in result['models']:
        far = item['name'].endswith('_far')
        item['correlation_checks'] = checks(item['active_electrons'], item['active_orbitals'],
                                           diagnostic=far, mrci=True)
        # Match the original finite-stencil IPEA=0.25 potential too.
        if not far:
            original = deepcopy(item['correlation_checks'][0])
            original['ipea'] = .25
            item['correlation_checks'].insert(0, original)
        limits(item, controls)
    return result


def prepare(output, recipe_path):
    requested = definition(recipe_path)
    output = Path(output).resolve()
    if output.exists():
        saved = output / 'review_plan.json'
        if not (output / 'campaign.json').is_file() or not saved.is_file() or json.loads(saved.read_text()) != requested:
            raise ValueError('Existing addition is incomplete or differs; preserve it and choose a new path.')
        return output / 'campaign.json'
    with tempfile.TemporaryDirectory(prefix='kinbot-cas-addition-') as temporary:
        source = Path(temporary) / 'partitions.json'
        atomic_json(source, requested)
        return orbital_review.prepare_partitions(source, output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--recipe', type=Path, default=Path(__file__).with_name('cas75_candidates.json'))
    args = parser.parse_args()
    print(prepare(args.output, args.recipe))
