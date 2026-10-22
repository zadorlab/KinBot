"""Prepare/run a finite orbital screen for an existing VRC review request.

Scientific choices live in the supplied recipe. Only the explicit ``review
--apply`` command installs a human-reviewed selection in an existing PES
input. It uses live refresh and never starts/restarts a production driver.
"""
import argparse
from copy import deepcopy
import fcntl
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


def _saved(root, *, allow_accepted=False):
    record = _read(root/'review_launch.json')
    for name, digest in record['files'].items():
        if file_sha256(root/name) != digest:
            raise ValueError(f'Prepared evidence changed: {name}.')
    if record.get('reviews'):
        for child in record['reviews']:
            saved = _saved(Path(child), allow_accepted=allow_accepted)
            if saved['production_policy'] != record['production_policy']:
                raise ValueError('Batch reviews no longer share one resource policy.')
        return record
    if record.get('provisional'):
        original = _read(root/'provisional_source.json')
        current = _read(record['request'])
        if any(current.get(k) != original.get(k) for k in ('structure', 'charge', 'mult')):
            raise ValueError('Provisional parent structure/spin changed; inspect the saved screen.')
        return record
    live, original = _read(record['request']), _read(root/'request.json')
    if allow_accepted and live.get('status') == 'model_accepted':
        # Acceptance changes only bookkeeping. A repeated review must still
        # match the saved chemistry/geometry and the installed channel below.
        live.pop('surface_model', None)
        live['status'] = original.get('status')
    if live != original:
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


def prepare(request, case_input, recipe, output, *, provisional=False):
    request, case_input, recipe = [Path(p).expanduser().resolve(strict=True)
                                    for p in (request, case_input, recipe)]
    output = Path(output).expanduser().resolve()
    evidence, controls = _read(request), _read(recipe)
    preliminary = None
    if provisional:
        if request.name != request.parent.name+'.json':
            raise ValueError('Provisional screening requires the existing well input.')
        preliminary = deepcopy(evidence)
        structure = evidence.get('structure', [])
        if not structure or len(structure) % 4:
            raise ValueError('Provisional well input needs explicit atom coordinates.')
        match = controls.get('request_match', {})
        evidence = dict(schema=1, status='needs_review',
            parent=request.parent.name, reaction=match.get('reaction'), bond=match.get('bond'),
            symbols=structure[::4], positions=[structure[i+1:i+4] for i in range(0,len(structure),4)],
            charge=evidence.get('charge'), multiplicity=evidence.get('mult'),
            products=match.get('products'), smiles=evidence.get('smiles', ''),
            provisional=True, note='Orbital screen only; reaction/product acceptance remains pending.')
        if (evidence['reaction'] != evidence['parent']+'_hom_sci_'+'_'.join(map(str,evidence['bond'] or []))
                or not evidence['products']):
            raise ValueError('Provisional recipe must explicitly identify the cleavage and fragment states.')
    if evidence.get('schema') != 1 or evidence.get('status') != 'needs_review':
        raise ValueError('Expected a current, unaccepted VRC orbital review request.')
    if controls.get('schema') != 1 or not controls.get('request_match'):
        raise ValueError('Recipe must identify its intended review request.')
    for key, value in controls['request_match'].items():
        if evidence.get(key) != value:
            raise ValueError(f'Request differs from the reviewed recipe: {key}.')
    if provisional:
        monitor = request.parent/'kinbot_monitor.out'
        rows = [s.split() for s in monitor.read_text().splitlines()] if monitor.exists() else []
        if not any(len(s)>=3 and s[2]==evidence['reaction'] and s[0]!='-999' for s in rows):
            raise ValueError('Provisional cleavage is not an active/accepted search in this well.')
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
    if preliminary is not None:
        atomic_json(output/'provisional_source.json', preliminary)
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
    if provisional:
        record['provisional'] = True
        record['files']['provisional_source.json'] = file_sha256(output/'provisional_source.json')
    atomic_json(output/'review_launch.json', record)
    return dict(directory=str(output), tasks=len(state['jobs']), max_jobs=record['max_jobs'],
                active_space=loaded['_probe_space'], states=state['jobs'][0]['root_states'],
                partition=config['partition'], excluded_nodes=policy['excluded_nodes'])


def batch(production, plan, output):
    """Stage independent orbital-only screens with one shared submission cap.

    Paths and scientific settings are explicit in the plan/recipes. Existing
    review selections are skipped; no rejected or imaginary future channel is
    manufactured. A provisional source can never directly approve production.
    """
    production = Path(production).expanduser().resolve(strict=True)
    plan = Path(plan).expanduser().resolve(strict=True)
    output = Path(output).expanduser().resolve()
    spec = _read(plan)
    entries = spec.get('screens', [])
    if (spec.get('schema') != 1 or not entries or type(spec.get('max_jobs')) is not int
            or spec['max_jobs'] < 1):
        raise ValueError('Batch plan requires schema 1, screens and positive max_jobs.')
    if output.exists():
        raise FileExistsError('Batch already exists; use run to resume without duplicate preparation.')
    if production == output or production in output.parents:
        raise ValueError('Keep reconnaissance outside the production run directory.')
    pending, skipped, names, policies = [], [], set(), set()
    for item in entries:
        name = item['name']
        if not recon._SAFE.fullmatch(name) or name in names:
            raise ValueError('Batch screen names must be safe and distinct.')
        names.add(name)
        case = (production/item['case']).resolve(strict=True)
        if case.parent != production:
            raise ValueError('Batch case must be directly inside the production directory.')
        recipe = (plan.parent/item['recipe']).resolve(strict=True)
        control = _read(recipe)
        config = control['reconnaissance']
        if config.get('correlation_checks') or config.get('probe_method') != 'casci':
            raise ValueError('This batch command stages orbital-only CASCI screens.')
        match = control['request_match']
        parent, reaction = match['parent'], match['reaction']
        if not all(recon._SAFE.fullmatch(x) for x in (parent,reaction)):
            raise ValueError('Unsafe parent/reaction in batch recipe.')
        context = _read(case/'input.json')
        model = context.get('vrc_tst_surface_models', {}).get(parent, {})
        if any(sorted(c['bond'])==sorted(match['bond']) for c in model.get('channels', [])):
            skipped.append(dict(name=name, reason='channel already supplied in case input'))
            continue
        request = case/parent/'vrctst/orbital_review'/reaction/'request.json'
        provisional = not request.exists()
        if provisional:
            if not item.get('allow_provisional', False):
                raise ValueError(f'{name}: no live orbital review request yet.')
            request = case/parent/(parent+'.json')
        policy_file = run_policy.find(case)
        if policy_file is None:
            raise ValueError('Every batch case needs a production resource policy.')
        policies.add(str(policy_file))
        pending.append((name,request,case/'input.json',recipe,provisional))
    if not pending:
        return dict(status='nothing_to_prepare', skipped=skipped, tasks=0)
    if len(policies)!=1:
        raise ValueError('A batch requires one shared production resource policy.')
    # Preflight the entire batch first, so a bad later recipe cannot publish a
    # partially usable campaign. This performs no scheduler/native calls.
    with tempfile.TemporaryDirectory(prefix='kinbot-batch-preflight-') as tmp:
        for name,request,case_input,recipe,provisional in pending:
            prepare(request,case_input,recipe,Path(tmp)/name,provisional=provisional)
    output.mkdir(parents=True, exist_ok=False)
    systems, reviews, caps = [], [], []
    for name,request,case_input,recipe,provisional in pending:
        child = output/name
        prepare(request,case_input,recipe,child,provisional=provisional)
        saved = _saved(child)
        systems.extend(saved['systems']); reviews.append(str(child)); caps.append(saved['max_jobs'])
    atomic_json(output/'plan.json',spec)
    record = dict(schema=1,status='prepared', request=str(pending[0][1]),
        production_policy=next(iter(policies)),max_jobs=min(spec['max_jobs'],*caps),
        systems=systems,reviews=reviews,skipped=skipped,
        files={'plan.json':file_sha256(output/'plan.json')})
    atomic_json(output/'review_launch.json',record)
    return dict(directory=str(output),screens=len(reviews),skipped=skipped,
        tasks=sum(len(_read(s['manifest'])['jobs']) for s in systems),max_jobs=record['max_jobs'],
        production_model_accepted=False)


def bind(source, request):
    """Attach a completed provisional screen to a subsequently accepted channel.

    This only binds evidence. The separate reviewed selection and live refresh
    remain required. Geometry changes are recorded; atom order/connectivity,
    fragment spin/charge, parent and bond must agree.
    """
    from kinbot.electronic_states import _canonical_smiles
    from kinbot.stationary_pt import StationaryPoint
    source = Path(source).expanduser().resolve(strict=True)
    request = Path(request).expanduser().resolve(strict=True)
    with recon._manifest_lock(source/'native'):
        record = _saved(source)
        if not record.get('provisional'):
            raise ValueError('Only a provisional screen needs binding.')
        if not all(j['status']=='complete' for j in _read(record['manifest'])['jobs']):
            raise ValueError('Finish the orbital screen before binding its evidence.')
        old, current = _read(source/'request.json'), _read(request)
        if (current.get('status')!='needs_review' or
                any(current.get(k)!=old.get(k) for k in ('parent','reaction','symbols','bond','charge','multiplicity'))
                or request != Path(record['request']).parent/'vrctst/orbital_review'/old['reaction']/'request.json'):
            raise ValueError('Live review request does not match the provisional channel.')
        def products(m):
            return sorted((_canonical_smiles(p['smiles']),p['charge'],p['multiplicity']) for p in m['products'])
        if products(old)!=products(current):
            raise ValueError('Accepted fragments differ from the provisional fragment states.')
        molecules=[]
        for m in (old,current):
            p=StationaryPoint('screen',m['charge'],m['multiplicity'],
                             atom=np.array(m['symbols']),geom=np.array(m['positions']))
            p.characterize(); molecules.append(p)
        if not np.array_equal(molecules[0].bond01,molecules[1].bond01):
            raise ValueError('Parent atom connectivity/order changed; review cannot reuse this screen.')
        atomic_json(source/'provisional_request.json',old)
        atomic_json(source/'request.json',current)
        record.update(provisional=False,request=str(request),
                      binding=dict(geometry_changed=old['positions']!=current['positions']))
        for name in ('provisional_request.json','request.json'):
            record['files'][name]=file_sha256(source/name)
        atomic_json(source/'review_launch.json',record)
    return dict(status='bound_for_review',request=str(request),production_model_accepted=False)


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
        options = {'keep_going': True} if record.get('reviews') else {}
        result = recon.run_many(manifests, max_jobs=min(record['max_jobs'], slots),
                               allocation_retries=2, **options)
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


def review(source, case_input, recipe, *, apply=False):
    """Publish an explicitly reviewed orbital selection to the existing owner.

    Selection is supplied by a human-reviewed recipe, never inferred from
    occupations or correlation energies. No reconnaissance jobs are launched.
    """
    case_input = Path(case_input).expanduser().resolve(strict=True)
    if not apply:
        return _review(source, case_input, recipe, apply=False)
    # Serialize helper invocations before reading/modifying the shared input;
    # otherwise two simultaneous channel additions could overwrite each other.
    with (case_input.parent/'.kinbot_orbital_review.lock').open('a+') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError('Another orbital review is editing this input; retry after it finishes.') from error
        return _review(source, case_input, recipe, apply=True)


def _review(source, case_input, recipe, *, apply):
    from kinbot import input_refresh
    from kinbot.molpro import validate_vrc_state_model
    from kinbot.vrc_multistate import validate_surface_model, extend_surface_model
    source = Path(source).expanduser().resolve(strict=True)
    case_input = Path(case_input).expanduser().resolve(strict=True)
    saved = _saved(source, allow_accepted=True)
    if saved.get('provisional'):
        raise ValueError('Bind the completed screen to the actual VRC review request before selection.')
    request, chosen = _read(source/'request.json'), _read(recipe)
    if (chosen.get('schema') != 1 or chosen.get('reviewed') is not True
            or not chosen.get('review_note') or not chosen.get('request_match')):
        raise ValueError('An explicit reviewed selection and chemistry note are required.')
    if any(request.get(k) != v for k,v in chosen['request_match'].items()):
        raise ValueError('Selected orbital recipe does not match this review request.')
    if case_input.parent not in Path(saved['request']).parents:
        raise ValueError('Case input does not belong to this review request.')
    if len(saved['systems']) != 1:
        raise ValueError('Select one source screen with basis-specific orbital records.')
    manifest = _read(saved['systems'][0]['manifest'])
    jobs = {j['name']: j for j in manifest['jobs']}
    active = deepcopy(chosen['active_space'])
    basis_models = {}
    for basis, selection in chosen['basis_seeds'].items():
        job = jobs[selection['job']]
        if (job['status'] != 'complete' or job['basis'] != basis
                or file_sha256(job['wfu']) != job['wfu_sha256']
                or file_sha256(Path(job['directory'])/(job['name']+'.out')) != job['output_sha256']):
            raise ValueError('Selected native seed is not unchanged, complete and basis matched.')
        model = deepcopy(chosen['state_model'])
        model['orbital_selection'] = dict(mode='molpro_record', path=job['wfu'],
            sha256=job['wfu_sha256'], record=selection['record'], orbital_type='natural',
            align_active_orbitals=False, active_indices=selection['active_indices'])
        basis_models[basis] = validate_vrc_state_model(model,
            active_electrons=active['active_electrons'], active_orbitals=active['active_orbitals'],
            nelectron=manifest['electron_count'], atom_count=len(request['symbols']))
    previous = _read(case_input)
    parameters = deepcopy(chosen.get('parameters', {}))
    parameters.update(vrc_tst_active_space=active, vrc_tst_state_model_file='', vrc_tst_state_model=dict(
        schema=2, molecule=dict(symbols=request['symbols'], charge=request['charge']),
        active_space=active, state_average_policy='reviewed', basis_models=basis_models))
    stage = dict(name='single', coordinate=chosen['coordinate'],
        distances_angstrom=deepcopy(previous[chosen['distances_parameter']]),
        states=deepcopy(chosen['states']), parameters=parameters)
    # Match JSON's string keys on repeat application (selected_roots is
    # normalized to integer keys by the state-model validator).
    channel = json.loads(json.dumps(dict(bond=request['bond'], combination='single', bottlenecks=[stage])))
    parent = request['parent']
    models = deepcopy(previous.get('vrc_tst_surface_models', {}))
    prior = models.get(parent)
    model = deepcopy(prior) if prior else dict(schema=1, reviewed=True,
        review_note=chosen['review_note'], channels=[])
    matches = [c for c in model['channels'] if sorted(c['bond']) == sorted(channel['bond'])]
    if matches and matches != [channel]:
        raise ValueError('This channel already has different settings; use a scientific migration.')
    if not matches and _read(saved['request']).get('status') == 'model_accepted':
        raise ValueError('Accepted review has no matching channel in this input; inspect its provenance.')
    if not matches:
        model['channels'].append(channel)
    models[parent] = extend_surface_model(prior, validate_surface_model(model))
    proposed = dict(previous, vrc_tst_surface_models=models)
    input_refresh.validate_models(previous, proposed, case_input.parent)
    if input_refresh.owner_directory(case_input) is None:
        raise ValueError('No live PES owner for this input; no driver was launched.')
    summary = dict(status='reviewed', parent=parent, active_space=active,
        basis_seeds={b:m['orbital_selection'] for b,m in basis_models.items()},
        drivers_started=0, reconnaissance_jobs_submitted=0)
    if not apply:
        return summary
    if previous != proposed:
        import hashlib
        digest = hashlib.sha256(json.dumps(proposed,sort_keys=True).encode()).hexdigest()
        archive = case_input.parent/'orbital_reviews/live_refresh'/digest
        archive.mkdir(parents=True, exist_ok=True)
        # Preserve the exact pre-edit input once. A retry may find the proposal
        # already installed; submit() reconciles its existing receipt then.
        for name, payload in [('input_before.json',previous), ('input_after.json',proposed),
                              ('selection.json',chosen)]:
            path = archive/name
            if path.exists() and _read(path) != payload:
                raise ValueError('Existing review archive differs; nothing was overwritten.')
            atomic_json(path,payload)
        if _read(case_input) != previous:
            raise ValueError('Case input changed during review preparation; retry from its new contents.')
        atomic_json(case_input,proposed)
    # Uses the existing journaled live-refresh path, which validates the owner,
    # accepted input, worker state and additive models. It never starts a driver.
    # A blocked submission leaves a recoverable input proposal and its backup.
    refresh = input_refresh.submit(case_input,apply=True)
    return dict(summary, refresh=refresh, status=refresh['status'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='action', required=True)
    prepare_cmd = commands.add_parser('prepare')
    for name in ('request', 'case_input', 'recipe', 'output'):
        prepare_cmd.add_argument(name)
    candidate_cmd = commands.add_parser('candidates')
    for name in ('source', 'recipe', 'output'):
        candidate_cmd.add_argument(name)
    batch_cmd = commands.add_parser('batch',help='Prepare one shared-cap orbital-only campaign.')
    for name in ('production','plan','output'):
        batch_cmd.add_argument(name)
    bind_cmd = commands.add_parser('bind',help='Bind preliminary evidence to an actual review request.')
    for name in ('source','request'):
        bind_cmd.add_argument(name)
    review_cmd = commands.add_parser('review', help='Apply an explicitly reviewed orbital selection without QC submissions.')
    for name in ('source', 'case_input', 'recipe'):
        review_cmd.add_argument(name)
    review_cmd.add_argument('--apply', action='store_true')
    for name in ('run', 'monitor', 'evidence'):
        commands.add_parser(name).add_argument('output')
    args = parser.parse_args()
    if args.action == 'prepare':
        result = prepare(args.request, args.case_input, args.recipe, args.output)
    elif args.action == 'batch':
        result = batch(args.production,args.plan,args.output)
    elif args.action == 'bind':
        result = bind(args.source,args.request)
    elif args.action == 'candidates':
        result = candidates(args.source, args.recipe, args.output)
    elif args.action == 'review':
        result = review(args.source, args.case_input, args.recipe, apply=args.apply)
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
