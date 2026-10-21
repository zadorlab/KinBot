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

import numpy as np

from kinbot import molpro_reconnaissance as recon, run_policy
from kinbot.execution import atomic_json, file_sha256


def _read(path):
    return json.loads(Path(path).read_text())


def _saved(root):
    record = _read(root/'review_launch.json')
    for name, digest in record['files'].items():
        if file_sha256(root/name) != digest:
            raise ValueError(f'Prepared evidence changed: {name}.')
    if _read(record['request']) != _read(root/'request.json'):
        raise ValueError('Live review request changed; inspect before continuing this review.')
    return record


def candidates(source, recipe, output):
    """Prepare explicit CASSCF partitions from completed CASCI/CASSCF seeds.

    This does not select orbitals by occupation or accept a production model.
    Every partition, seed geometry and native check comes from the recipe.
    """
    source = Path(source).expanduser().resolve(strict=True)
    output = Path(output).expanduser().resolve()
    previous = _saved(source)
    controls = _read(recipe)
    request = _read(source/'request.json')
    if controls.get('schema') != 1 or not controls.get('request_match'):
        raise ValueError('Candidate recipe must identify its intended review request.')
    for key, value in controls['request_match'].items():
        if request.get(key) != value:
            raise ValueError(f'Request differs from candidate recipe: {key}.')
    policy_path = run_policy.find(Path(previous['request']).parent)
    policy = run_policy.read(Path(previous['request']).parent) or {}
    slots = policy.get('orbital_review_reserved_slots', 0)
    cap = controls.get('max_jobs', previous['max_jobs'])
    if type(cap) is not int or cap < 1 or slots < 1:
        raise ValueError('Candidate review needs positive reserved capacity.')
    if output.exists() or output == source or source in output.parents:
        raise FileExistsError('Preserve the source review; choose a fresh sibling directory.')
    if policy_path.parent == output or policy_path.parent in output.parents:
        raise ValueError('Keep reconnaissance outside the production run directory.')
    # A single completed source manifest is explicit. Multi-model successors
    # require choosing a source, rather than silently picking one candidate.
    if len(previous['systems']) != 1:
        raise ValueError('Choose one source model for a new partition comparison.')
    manifest_path = Path(previous['systems'][0]['manifest'])
    original = _read(manifest_path)
    if not original['jobs'] or any(j['status'] != 'complete' for j in original['jobs']):
        raise ValueError('Finish the source screen before preparing candidates.')
    loaded = recon._load(previous['systems'][0]['configuration'])
    anchors = {a['name']: a for a in loaded['anchors']}
    jobs = {j['name']: j for j in original['jobs']}
    definitions = controls.get('candidates')
    if not isinstance(definitions, list) or not definitions:
        raise ValueError('Candidate recipe needs a nonempty candidates list.')
    configurations, names = [], set()
    with tempfile.TemporaryDirectory(prefix='kinbot-candidate-preflight-') as tmp:
        for item in definitions:
            name = item['name']
            if not recon._SAFE.fullmatch(name) or name in names:
                raise ValueError('Candidate names must be safe and distinct.')
            names.add(name)
            data = {k: deepcopy(v) for k, v in loaded.items()
                    if not k.startswith('_') and k not in
                    ('coordinates', 'anchors', 'scan_policy', 'orientation_probes', 'initial_records')}
            definition = deepcopy(item['probe_definition'])
            partition = definition.get('orbital_partition', {})
            for key in ('active_electrons', 'active_orbitals'):
                if partition.get('source_'+key) != loaded['_probe_space'][key]:
                    raise ValueError('Partition counts differ from the source active space.')
            data.update(title=name, probe_method='casscf', orbital_initialization='independent',
                        probe_definition=definition, anchors=[], initial_records={},
                        casscf_options=deepcopy(item['casscf_options']),
                        correlation_checks=deepcopy(item['correlation_checks']))
            for choice in item['anchors']:
                if set(choice) - {'source', 'name', 'distance_angstrom', 'rotation'}:
                    raise ValueError('Candidate anchors accept source, name, distance and rotation.')
                anchor = deepcopy(anchors[choice['source']])
                anchor['name'] = choice.get('name', choice['source'])
                if 'distance_angstrom' in choice:
                    if 'bond' not in anchor:
                        raise ValueError('Distance changes require a dissociation coordinate.')
                    anchor['distance_angstrom'] = choice['distance_angstrom']
                if 'rotation' in choice:
                    rotate = choice['rotation']
                    if (not isinstance(rotate, dict) or set(rotate) !=
                            {'atoms', 'pivot_atom', 'axis_bond', 'angle_degrees'}):
                        raise ValueError('Rotation needs atoms, pivot_atom, axis_bond and angle_degrees.')
                    atoms, axis_bond, pivot = rotate['atoms'], rotate['axis_bond'], rotate['pivot_atom']
                    angle = rotate['angle_degrees']
                    if (not isinstance(atoms, list) or not atoms
                            or not isinstance(axis_bond, list) or len(axis_bond) != 2
                            or any(type(i) is not int or not 1 <= i <= len(data['symbols'])
                                   for i in [*atoms, *axis_bond, pivot])
                            or len(set(atoms)) != len(atoms) or pivot not in atoms
                            or axis_bond[0] == axis_bond[1] or type(angle) not in (int, float)
                            or not np.isfinite(angle)):
                        raise ValueError('Invalid one-based fragment rotation.')
                    xyz = recon._anchor_geometry(data, anchor)
                    axis = xyz[axis_bond[1]-1]-xyz[axis_bond[0]-1]
                    length = np.linalg.norm(axis)
                    if length < 1e-8:
                        raise ValueError('Rotation axis has zero length.')
                    anchor['_rotation'] = dict(atoms=[i-1 for i in atoms], pivot=pivot-1,
                                                axis=(axis/length).tolist(), angle_degrees=angle)
                    anchor['_orientation'] = dict(channel=anchor.get('_channel'), pivot_atom=pivot,
                        axis='bond_'+'_'.join(map(str, axis_bond)), angle_degrees=angle)
                    if np.allclose(xyz, recon._anchor_geometry(data, anchor), atol=1e-10, rtol=0.):
                        raise ValueError('Rotation would duplicate the unchanged geometry.')
                data['anchors'].append(anchor)
                for basis in data['bases']:
                    old = jobs[basis['name']+'__'+choice['source']]
                    if (old['probe_active_space']['active_electrons'] != loaded['_probe_space']['active_electrons']
                            or old['probe_active_space']['active_orbitals'] != loaded['_probe_space']['active_orbitals']
                            or old['basis'] != basis['basis']
                            or file_sha256(old['wfu']) != old['wfu_sha256']
                            or file_sha256(Path(old['directory'])/(old['name']+'.out')) != old['output_sha256']):
                        raise ValueError('Source checkpoint/output or model differs from validated evidence.')
                    data['initial_records'][basis['name']+'__'+anchor['name']] = dict(
                        path=old['wfu'], sha256=old['wfu_sha256'], basis=old['basis'])
            path = Path(tmp)/(name+'.json')
            atomic_json(path, data)
            check = recon._load(path)
            for basis in check['bases']:
                for anchor in check['anchors']:
                    stem = basis['name']+'__'+anchor['name']
                    recon.validate_input_blocks(recon._deck(check, basis['basis'], anchor, stem,
                        previous=check['initial_records'][stem]), stem)
            configurations.append(data)
    output.mkdir(parents=True, exist_ok=False)
    for name in ('request.json', 'electronic_context.json'):
        atomic_json(output/name, _read(source/name))
    atomic_json(output/'recipe.json', controls)
    atomic_json(output/'source_manifest.json', original)
    files = ['request.json', 'electronic_context.json', 'recipe.json', 'source_manifest.json']
    systems = []
    for data in configurations:
        config = output/(data['title']+'.json')
        data['kinbot_input'] = str(output/'electronic_context.json')
        atomic_json(config, data)
        manifest = recon.prepare(config, output/data['title'])
        state = _read(manifest)
        state.update(slurm_excluded_nodes=policy['excluded_nodes'], startup_retries=2)
        atomic_json(manifest, state)
        files.append(config.name)
        systems.append(dict(title=data['title'], configuration=str(config), manifest=str(manifest)))
    record = dict(schema=1, status='prepared', request=previous['request'],
                  max_jobs=min(cap, slots), production_policy=str(policy_path),
                  files={name: file_sha256(output/name) for name in files}, systems=systems)
    atomic_json(output/'review_launch.json', record)
    return dict(directory=str(output), models=len(systems), max_jobs=record['max_jobs'],
                tasks=sum(len(_read(s['manifest'])['jobs']) for s in systems),
                production_model_accepted=False)


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
    record = _saved(root)
    policy = run_policy.read(Path(record['request']).parent) or {}
    slots = policy.get('orbital_review_reserved_slots', 0)
    if slots < 1:
        raise ValueError('No review capacity remains reserved; production may have resumed.')
    # Check setup and scheduler, without starting native work on the login node.
    manifests = [Path(s['manifest']) for s in record['systems']]
    setup = manifests[0].parent/'site_setup.sh'
    env = dict(os.environ, KINBOT_BACKEND='molpro')
    subprocess.run(['bash', '-c', 'set -eo pipefail; source '+shlex.quote(str(setup))+
                    '; command -v molpro'], env=env, check=True, timeout=60)
    config = _read(record['systems'][0]['configuration'])
    from kinbot.run_walltime import partition_walltime
    print('Partition walltime:', partition_walltime(config['partition']), flush=True)
    # Hold the same manifest lock used by every submission path while updating
    # exclusions, then let run_many take it for the complete controller lifetime.
    for manifest in manifests:
        with recon._manifest_lock(manifest.parent):
            state = _read(manifest)
            state['slurm_excluded_nodes'] = sorted(set(state.get('slurm_excluded_nodes', []))
                                                   | set(policy['excluded_nodes']))
            atomic_json(manifest, state)
    previous = Path.cwd()
    try:
        os.chdir(root)
        result = recon.run_many(manifests, max_jobs=min(record['max_jobs'], slots),
                               allocation_retries=2)
    finally:
        os.chdir(previous)
    statuses = sorted({r['status'] for r in result})
    return dict(status=statuses[0] if len(statuses) == 1 else statuses,
                tasks=sum(len(r['jobs']) for r in result),
                production_model_accepted=False)


def evidence(output):
    """Saved evidence only; no audit mutation, scheduler call or acceptance."""
    root = Path(output).expanduser().resolve(strict=True)
    record = _read(root/'review_launch.json')
    lines = ['Recorded native results; no production model has been accepted.']
    for system in record['systems']:
        state = _read(system['manifest'])
        lines.append(system['title']+': '+state['status'])
        config = _read(system['configuration'])
        for j in state['jobs']:
            lines.append(j['name']+' '+j.get('status', 'prepared'))
            for field, label in (('state_energies_hartree', config.get('probe_method', 'casci').upper()),
                                 ('correlated_energies_hartree', 'CORRELATED')):
                rows = j.get(field, [])
                if rows:
                    lines.append(label+' Eh '+', '.join(
                        f'{r.get("method", "")}:m{r["multiplicity"]}r{r["root"]}={r["energy_hartree"]:.10f}'
                        for r in rows))
            n = (state['electron_count']-j['probe_active_space']['active_electrons'])//2
            stop = n+j['probe_active_space']['active_orbitals']
            occ = j.get('natural_occupations', [])
            if occ:
                lines.append('ACTIVE '+', '.join(f'{i+1}:{occ[i]:.6f}'
                             for i in range(n, min(stop, len(occ)))))
            for field in ('validation_errors', 'validation_warnings', 'terminal_failure'):
                if j.get(field):
                    lines.append(field+': '+json.dumps(j[field]))
        for key, value in state.get('coverage', {}).items():
            lines.append('COVERAGE '+key+' '+json.dumps({k:value.get(k) for k in
                ('status', 'asymptote_plateau', 'suggested_midpoints_angstrom',
                 'suggested_far_distance_angstrom')}))
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='action', required=True)
    prepare_cmd = commands.add_parser('prepare')
    for name in ('request', 'case_input', 'recipe', 'output'):
        prepare_cmd.add_argument(name)
    candidate_cmd = commands.add_parser('candidates')
    for name in ('source', 'recipe', 'output'):
        candidate_cmd.add_argument(name)
    for name in ('run', 'monitor', 'evidence'):
        commands.add_parser(name).add_argument('output')
    args = parser.parse_args()
    if args.action == 'prepare':
        result = prepare(args.request, args.case_input, args.recipe, args.output)
    elif args.action == 'candidates':
        result = candidates(args.source, args.recipe, args.output)
    elif args.action == 'run':
        result = run(args.output)
    elif args.action == 'evidence':
        print(evidence(args.output))
        return
    else:
        from kinbot.orbital_monitor import snapshot, render
        root = Path(args.output).expanduser().resolve()
        data = snapshot(root, registry=root/'review_launch.json', include_all=True)
        states = sorted({_read(s['manifest'])['status']
                         for s in _read(root/'review_launch.json')['systems']})
        data['status'] = ', '.join(states)
        print(render(data, details=True))
        return
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
