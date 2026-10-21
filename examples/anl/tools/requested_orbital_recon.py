"""Prepare/run a finite orbital screen for an existing VRC review request.

Scientific choices live in the supplied recipe. This tool never accepts a
production model, edits a PES input, or starts/restarts a production driver.
"""
import argparse
from copy import deepcopy
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

from kinbot import molpro_reconnaissance as recon, run_policy
from kinbot.execution import atomic_json, file_sha256


def _read(path):
    return json.loads(Path(path).read_text())


def prepare(request, case_input, recipe, output):
    request, case_input, recipe = [Path(p).expanduser().resolve(strict=True)
                                    for p in (request, case_input, recipe)]
    output = Path(output).expanduser().resolve()
    evidence, controls = _read(request), _read(recipe)
    if evidence.get('schema') != 1 or evidence.get('status') != 'needs_review':
        raise ValueError('Expected a current, unaccepted VRC orbital review request.')
    if controls.get('schema') != 1 or not controls.get('request_match'):
        raise ValueError('Recipe must identify its intended review request.')
    for key, value in controls['request_match'].items():
        if evidence.get(key) != value:
            raise ValueError(f'Request differs from the reviewed recipe: {key}.')
    if case_input.parent not in request.parents:
        raise ValueError('The case input must belong to the requested calculation tree.')
    policy_file = run_policy.find(request.parent)
    policy = run_policy.read(request.parent) or {}
    slots = policy.get('orbital_review_reserved_slots', 0)
    if slots < 1:
        raise ValueError('No orbital-review capacity is configured/currently reserved.')
    requested_cap = controls.get('max_jobs', slots)
    if type(requested_cap) is not int or requested_cap < 1:
        raise ValueError('Recipe max_jobs must be a positive integer.')
    if policy_file.parent == output or policy_file.parent in output.parents:
        raise ValueError('Keep reconnaissance outside the production run directory.')
    context = _read(case_input)
    config = deepcopy(controls['reconnaissance'])
    config.update(symbols=evidence['symbols'], geometry_angstrom=evidence['positions'],
                  charge=evidence['charge'], partition=context['queue_name'])
    # The recipe may explicitly choose a high-spin ROHF initialization; the
    # parent and its production multiplicity remain in the immutable request.
    config.setdefault('multiplicity', evidence['multiplicity'])
    config['schema'] = 1
    config['kinbot_input'] = str(case_input)
    with tempfile.TemporaryDirectory(prefix='kinbot-request-preflight-') as temporary:
        source = Path(temporary)/'input.json'
        atomic_json(source, config)
        loaded = recon._load(source)
        # Validate every generated native deck before publishing the directory.
        for basis in loaded['bases']:
            for anchor in loaded['anchors']:
                recon.validate_input_blocks(recon._deck(
                    loaded, basis['basis'], anchor, basis['name']+'__'+anchor['name']), 'screen')
    output.mkdir(parents=True, exist_ok=False)  # Repeating preparation cannot duplicate jobs.
    atomic_json(output/'request.json', evidence)
    atomic_json(output/'electronic_context.json', context)
    atomic_json(output/'recipe.json', controls)
    config['kinbot_input'] = str(output/'electronic_context.json')
    atomic_json(output/'input.json', config)
    manifest = recon.prepare(output/'input.json', output/'native')
    state = _read(manifest)
    state.update(slurm_excluded_nodes=policy['excluded_nodes'], startup_retries=2)
    atomic_json(manifest, state)
    record = dict(schema=1, status='prepared', request=str(request), manifest=str(manifest),
                  max_jobs=min(requested_cap, slots), production_policy=str(policy_file),
                  files={name: file_sha256(output/name) for name in
                         ('request.json', 'recipe.json', 'electronic_context.json', 'input.json')},
                  systems=[dict(title=config['title'], configuration=str(output/'input.json'),
                                manifest=str(manifest))])
    atomic_json(output/'review_launch.json', record)
    return dict(directory=str(output), tasks=len(state['jobs']), max_jobs=record['max_jobs'],
                active_space=loaded['_probe_space'], states=state['jobs'][0]['root_states'],
                partition=config['partition'], excluded_nodes=policy['excluded_nodes'])


def run(output):
    root = Path(output).expanduser().resolve(strict=True)
    record = _read(root/'review_launch.json')
    for name, digest in record['files'].items():
        if file_sha256(root/name) != digest:
            raise ValueError(f'Prepared evidence changed: {name}.')
    current, saved = _read(record['request']), _read(root/'request.json')
    if current != saved:
        raise ValueError('Live review request changed; inspect before running this screen.')
    policy = run_policy.read(Path(record['request']).parent) or {}
    slots = policy.get('orbital_review_reserved_slots', 0)
    if slots < 1:
        raise ValueError('No review capacity remains reserved; production may have resumed.')
    # Check setup and scheduler, without starting native work on the login node.
    setup = Path(record['manifest']).parent/'site_setup.sh'
    env = dict(os.environ, KINBOT_BACKEND='molpro')
    subprocess.run(['bash', '-c', 'set -eo pipefail; source '+shlex.quote(str(setup))+
                    '; command -v molpro'], env=env, check=True, timeout=60)
    config = _read(root/'input.json')
    from kinbot.run_walltime import partition_walltime
    print('Partition walltime:', partition_walltime(config['partition']), flush=True)
    # Hold the same manifest lock used by every submission path while updating
    # exclusions, then let run_many take it for the complete controller lifetime.
    manifest = Path(record['manifest'])
    with recon._manifest_lock(manifest.parent):
        state = _read(manifest)
        state['slurm_excluded_nodes'] = sorted(set(state.get('slurm_excluded_nodes', []))
                                               | set(policy['excluded_nodes']))
        atomic_json(manifest, state)
    previous = Path.cwd()
    try:
        os.chdir(root)
        result = recon.run_many([manifest], max_jobs=min(record['max_jobs'], slots),
                               allocation_retries=2)
    finally:
        os.chdir(previous)
    return dict(status=result[0]['status'], tasks=len(result[0]['jobs']),
                production_model_accepted=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='action', required=True)
    prepare_cmd = commands.add_parser('prepare')
    for name in ('request', 'case_input', 'recipe', 'output'):
        prepare_cmd.add_argument(name)
    for name in ('run', 'monitor'):
        commands.add_parser(name).add_argument('output')
    args = parser.parse_args()
    if args.action == 'prepare':
        result = prepare(args.request, args.case_input, args.recipe, args.output)
    elif args.action == 'run':
        result = run(args.output)
    else:
        from kinbot.orbital_monitor import snapshot, render
        root = Path(args.output).expanduser().resolve()
        data = snapshot(root, registry=root/'review_launch.json', include_all=True)
        data['status'] = _read(_read(root/'review_launch.json')['manifest'])['status']
        print(render(data, details=True))
        return
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
