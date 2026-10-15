"""Chemically annotated, finite orbital comparisons before user selection.

This is a screening and validation aid, not an automatic CAS or kinetic-model
selector. Native jobs use orbital_campaign's locks, pilots and shared job cap.
"""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import tempfile

import numpy as np

from kinbot import molpro_reconnaissance as recon
from kinbot.execution import atomic_json, file_sha256


def _channel(anchor):
    return anchor.get('_channel') or anchor.get('_orientation', {}).get('channel')


def _fresh(path):
    path = Path(path).expanduser().resolve()
    if path.exists():
        raise FileExistsError(f'Preserving existing evidence: choose a fresh output: {path}')
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def plan(source, output):
    """Expose chemistry, available geometries and editable region/model choices."""
    data = recon._load(source)
    if not data.get('probe_active_space'):
        raise ValueError('Orbital review requires probe_active_space=true.')
    channels = sorted({_channel(anchor) for anchor in data['anchors']} - {None})
    if not channels:
        raise ValueError('Use dissociation coordinates to provide fragment/channel context.')
    models = []
    for channel in channels:
        anchors = [a for a in data['anchors'] if _channel(a) == channel]
        radial = sorted((a for a in anchors if a.get('_channel')),
                        key=lambda a: a['distance_angstrom'])
        # A finite panel, including off-axis geometries and independent far starts.
        # These are review locations, never an inferred kinetic bottleneck.
        indices = {0, len(radial) // 2, len(radial) - 1}
        bond = radial[0]['bond']
        geom = np.asarray(data['geometry_angstrom'])
        bonded = np.linalg.norm(geom[bond[0] - 1] - geom[bond[1] - 1])
        indices.add(min(range(len(radial)), key=lambda i: abs(radial[i]['distance_angstrom'] - bonded)))
        panel = [a['name'] for i, a in enumerate(radial) if i in indices]
        oriented = [a for a in anchors if a.get('_orientation')]
        # One nearest probe per rotation axis plus a far independent orientation.
        # All generated anchors remain available for explicit expansion.
        axes = sorted({a['_orientation']['axis'] for a in oriented})
        panel += [min((a for a in oriented if a['_orientation']['axis'] == axis),
                      key=lambda a: abs(a['distance_angstrom'] - bonded))['name'] for axis in axes]
        if oriented:
            panel.append(max(oriented, key=lambda a: a['distance_angstrom'])['name'])
        panel = list(dict.fromkeys(panel))
        models.append({
            'name': channel + '_valence_reference', 'channel': channel,
            'region': 'diagnostic',
            'chemical_note': 'Valence reference from configured atomic targets; inspect radical, '
                'forming/breaking bond, resonance and degenerate fragment manifolds. '
                'This is not a selected production CAS or state average.',
            'active_electrons': data['_probe_space']['active_electrons'],
            'active_orbitals': data['_probe_space']['active_orbitals'],
            'state_model': recon._probe_model(data, radial[0]),
            'anchors': panel, 'correlation_checks': [],
            'fragment_context': radial[0].get('_context'), 'bond': bond,
        })
    source = Path(source).expanduser().resolve()
    context = data.get('kinbot_input')
    context_path = (source.parent / context).resolve() if context else None
    payload = {'schema': 1, 'source': str(source), 'source_sha256': file_sha256(source),
        'context_sha256': file_sha256(context_path) if context_path else None,
        'purpose': 'User-aided chemical orbital/state selection; no automatic production approval.',
        'models': models,
        'variants': [
            {'name': 'casci_fresh', 'probe_method': 'casci', 'orbital_initialization': 'independent'},
            {'name': 'casscf_fresh', 'probe_method': 'casscf', 'orbital_initialization': 'independent'},
            {'name': 'casscf_propagated', 'probe_method': 'casscf', 'orbital_initialization': 'propagated'}],
        'available_anchors': data['anchors'],
        'inner_outer_guidance': {
            'inner': 'Reactive-atom/pivot dividing surfaces; identify attractive contributing states.',
            'outer': 'Fragment center-of-mass dividing surfaces; retain required fragment degeneracies.',
            'region_boundaries': 'User-defined from potentials and flux minimization; not inferred from occupations.',
            'multiple_spaces': 'Duplicate a model with region=inner or outer, choose anchors and state averages independently.',
            'vdw': 'A potential minimum is neither necessary nor sufficient for two variational bottlenecks.'}}
    if context_path:
        context_data = json.loads(context_path.read_text())
        surface_name = context_data.get('vrc_tst_surface_model_file')
        if surface_name:
            surface_path = (context_path.parent / surface_name).resolve()
            payload['surface_model_context'] = {'path': str(surface_path),
                'available': surface_path.is_file(),
                'model': json.loads(surface_path.read_text()) if surface_path.is_file() else None}
    target = _fresh(output)
    atomic_json(target, payload)
    return target


def prepare(plan_file, output, *, additional_plans=(), reuse_from=None):
    """Prepare comparisons without changing or rerunning an existing campaign."""
    from kinbot import orbital_campaign
    if additional_plans:
        return _prepare_many([plan_file, *additional_plans], output, reuse_from=reuse_from)
    plan_file = Path(plan_file).expanduser().resolve()
    definition = json.loads(plan_file.read_text())
    if definition.get('schema') != 1:
        raise ValueError('Orbital review plan requires schema 1.')
    source = (plan_file.parent / definition['source']).resolve()
    if file_sha256(source) != definition['source_sha256']:
        raise ValueError('Recon source changed since planning; regenerate the plan.')
    data = recon._load(source)
    if data.get('kinbot_input'):
        context = (source.parent / data['kinbot_input']).resolve()
        if file_sha256(context) != definition.get('context_sha256'):
            raise ValueError('Electronic context changed since planning; regenerate the plan.')
        data['kinbot_input'] = str(context)
    anchors = {a['name']: a for a in data['anchors']}
    models, variants = definition.get('models'), definition.get('variants')
    if not isinstance(models, list) or not models or not isinstance(variants, list) or not variants:
        raise ValueError('Review plan needs models and variants.')
    names, configurations, index = set(), [], []
    reusable = []
    if reuse_from:
        old = Path(reuse_from).expanduser().resolve()
        previous = json.loads((old / 'campaign.json').read_text())
        for system in previous['systems']:
            for key in ('manifest', 'pilot_manifest'):
                manifest = json.loads(Path(system[key]).read_text())
                reusable.extend(j for j in manifest['jobs'] if j.get('status') == 'complete'
                                and j.get('natural_orbital_probe'))
    with tempfile.TemporaryDirectory(prefix='kinbot-orbital-review-') as temporary:
        for model in models:
            if (set(model) - {'name', 'channel', 'region', 'chemical_note', 'active_electrons',
                    'active_orbitals', 'state_model', 'anchors', 'correlation_checks', 'fragment_context', 'bond'}
                    or not recon._SAFE.fullmatch(str(model.get('name', '')))
                    or model.get('region') not in ('diagnostic', 'single', 'inner', 'outer')
                    or not isinstance(model.get('anchors'), list) or not model['anchors']
                    or len(set(model['anchors'])) != len(model['anchors'])
                    or any(name not in anchors or _channel(anchors[name]) != model.get('channel')
                           for name in model['anchors'])):
                raise ValueError('Each model needs a safe name, region and unique anchors from its channel.')
            selected = [anchors[name] for name in model['anchors']]
            # Radial propagation follows distance, independently of plan display order.
            selected.sort(key=lambda a: (bool(a.get('_orientation')), a['distance_angstrom'], a['name']))
            for variant in variants:
                if (set(variant) != {'name', 'probe_method', 'orbital_initialization'}
                        or not recon._SAFE.fullmatch(str(variant.get('name', '')))):
                    raise ValueError('Invalid comparison variant.')
                title = model['name'] + '__' + variant['name']
                if title in names:
                    raise ValueError('Duplicate model/variant name.')
                names.add(title)
                config = {key: deepcopy(value) for key, value in data.items()
                          if not key.startswith('_') and key not in ('coordinates', 'scan_policy', 'orientation_probes')}
                config.update(title=title, anchors=selected, avas_reference='rhf',
                              probe_method=variant['probe_method'],
                              orbital_initialization=variant['orbital_initialization'])
                config['probe_definition'] = {key: model[key] for key in (
                    'active_electrons', 'active_orbitals', 'state_model', 'chemical_note')}
                config['multiplicity'] = model['state_model']['state_average'][0]['multiplicity']
                config['correlation_checks'] = (model.get('correlation_checks', [])
                                                if variant['probe_method'] == 'casscf' else [])
                config.pop('initial_records', None)
                # Reoptimize old CASCI WFU records; never relabel their energies
                # as CASSCF results. Only the matching reference space can reuse them.
                seeds = {}
                if reusable and variant['orbital_initialization'] == 'propagated' and variant['probe_method'] == 'casscf':
                    for basis in config['bases']:
                        for anchor in selected:
                            matches = [job for job in reusable
                                if job['basis'] == basis['basis'] and job['anchor'] == anchor['name']
                                and job['symbols'] == config['symbols']
                                and all((job.get('probe_active_space') or {}).get(key) == model[key]
                                        for key in ('active_electrons', 'active_orbitals'))
                                and (job.get('probe_active_space') or {}).get('orbital_selection') ==
                                    model['state_model'].get('orbital_selection')]
                            if matches:
                                prior = matches[0]
                                # Older manifests lack geometry: check it from the
                                # immutable prepared deck instead of assuming names match.
                                original = Path(prior['directory']) / (prior['name'] + '.inp')
                                expected_rows = [f'{symbol} {x:.12f} {y:.12f} {z:.12f}' for symbol, (x, y, z)
                                                 in zip(config['symbols'], recon._anchor_geometry(data, anchor))]
                                if (not original.is_file() or file_sha256(original) != prior['input_sha256']
                                        or not all(row in original.read_text().splitlines() for row in expected_rows)
                                        or f'wf,{data["_electrons"]},1,' not in original.read_text().lower()
                                        or not Path(prior['wfu']).is_file()
                                        or file_sha256(prior['wfu']) != prior['wfu_sha256']):
                                    raise ValueError('Reusable checkpoint input, geometry or WFU changed.')
                                seeds[basis['name'] + '__' + anchor['name']] = {
                                    'path': prior['wfu'], 'sha256': prior['wfu_sha256']}
                    config['initial_records'] = seeds
                target = Path(temporary) / (title + '.json')
                atomic_json(target, config)
                recon._load(target)  # Validate all definitions before creating the campaign.
                configurations.append(target)
                index.append({'title': title, 'model': model['name'], 'channel': model['channel'],
                              'region': model['region'], 'chemical_note': model['chemical_note'],
                              'variant': variant, 'bond': selected[0]['bond'],
                              'reused_checkpoint_count': len(seeds)})
        result = orbital_campaign.prepare(output, configurations)
    root = result.parent
    atomic_json(root / 'review_plan.json', definition)
    atomic_json(root / 'review_index.json', {'schema': 1, 'systems': index})
    return result


def _prepare_many(plans, output, *, reuse_from=None):
    """Combine prepared definitions into one cap/lock without submitting jobs."""
    from kinbot import orbital_campaign
    configurations, indices, definitions = [], [], []
    with tempfile.TemporaryDirectory(prefix='kinbot-review-plans-') as temporary:
        for number, source in enumerate(plans, 1):
            prepared = prepare(source, Path(temporary) / str(number), reuse_from=reuse_from)
            campaign = json.loads(prepared.read_text())
            index = json.loads((prepared.parent / 'review_index.json').read_text())['systems']
            definitions.append(json.loads((prepared.parent / 'review_plan.json').read_text()))
            for record, metadata in zip(campaign['systems'], index):
                data = json.loads(Path(record['configuration']).read_text())
                data['title'] = f'system_{number}__' + data['title']
                metadata['title'] = data['title']
                # Model identity includes the source system, so unrelated
                # geometries cannot be compared merely because channels share a name.
                metadata['model'] = f'system_{number}__' + metadata['model']
                config = Path(temporary) / (data['title'] + '.json')
                atomic_json(config, data)
                configurations.append(config)
                indices.append(metadata)
        result = orbital_campaign.prepare(output, configurations)
    atomic_json(result.parent / 'review_plan.json', {'schema': 1, 'plans': definitions})
    atomic_json(result.parent / 'review_index.json', {'schema': 1, 'systems': indices})
    return result


def _signature(job):
    return [(e['multiplicity'], e['root']) for e in job.get('state_energies_hartree', [])]


def report(root, output):
    """Read atomic snapshots; retain native tables and label diagnostic limits."""
    root = Path(root).expanduser().resolve()
    campaign = json.loads((root / 'campaign.json').read_text())
    index = {item['title']: item for item in json.loads((root / 'review_index.json').read_text())['systems']}
    rows, comparisons = [], []
    for system in campaign['systems']:
        metadata = index[system['title']]
        manifest_path = Path(system['manifest'])
        manifest = json.loads(manifest_path.read_text())
        for job in manifest['jobs']:
            row = {**deepcopy(job), 'review': metadata, 'manifest': str(manifest_path)}
            if job.get('status') == 'complete':
                directory = Path(job['directory'])
                evidence = [(directory / (job['name'] + '.inp'), job['input_sha256']),
                            (directory / (job['name'] + '.out'), job['output_sha256']),
                            (Path(job['wfu']), job['wfu_sha256'])]
                if any(not path.is_file() or file_sha256(path) != digest for path, digest in evidence):
                    raise ValueError(f'{job["name"]}: completed evidence changed; audit before review.')
                native = directory / (job['name'] + '.out')
                text = native.read_text(errors='replace')
                # Preserve complete native character/CI tables, avoiding assumptions
                # about AO labels, symmetry, orbital rank identities or elements.
                row['native_output'] = str(native)
                row['native_orbital_and_ci_text'] = text[text.find('NATURAL ORBITALS'):] if 'NATURAL ORBITALS' in text else text
                execution = directory / 'execution.json'
                row['elapsed_seconds'] = (json.loads(execution.read_text()).get('elapsed_seconds')
                                          if execution.is_file() else None)
            rows.append(row)
    completed = [row for row in rows if row.get('status') == 'complete']
    for i, left in enumerate(completed):
        for right in completed[i + 1:]:
            if (left['basis'] != right['basis'] or left['review']['model'] != right['review']['model']
                    or left.get('probe_method') != right.get('probe_method')
                    or left['symbols'] != right['symbols'] or _signature(left) != _signature(right)):
                continue
            x, y = np.asarray(left['geometry_angstrom']), np.asarray(right['geometry_angstrom'])
            delta = np.max(np.abs(np.linalg.norm(x[:, None] - x, axis=2) -
                                  np.linalg.norm(y[:, None] - y, axis=2)))
            if delta < 1.e-7:
                differences = [b['energy_hartree'] - a['energy_hartree'] for a, b in
                               zip(left['state_energies_hartree'], right['state_energies_hartree'])]
                comparisons.append({'left': left['directory'], 'right': right['directory'],
                    'method': left['probe_method'], 'maximum_distance_matrix_difference_angstrom': float(delta),
                    'state_energy_differences_hartree': differences,
                    'limitation': 'Root indices are labels, not a demonstrated cross-geometry state assignment.'})
    groups = {}
    for row in completed:
        if row.get('channel'):
            groups.setdefault((row['review']['title'], row['basis'], row['channel']), []).append(row)
    profiles = []
    for key, group in groups.items():
        group.sort(key=lambda row: row['distance_angstrom'])
        far = group[-1]
        profiles.append({'system': key[0], 'basis': key[1], 'channel': key[2],
            'reference_job': far['name'], 'reference_distance_angstrom': far['distance_angstrom'],
            'reference_status': 'finite separation; asymptotic adequacy requires review',
            'points': [{'job': row['name'], 'distance_angstrom': row['distance_angstrom'],
                'states': [{**state, 'relative_to_far_hartree': state['energy_hartree'] - reference['energy_hartree']}
                    for state, reference in zip(row.get('state_energies_hartree', []),
                                                 far.get('state_energies_hartree', []))],
                'correlated': [{**state, 'relative_to_far_hartree': state['energy_hartree'] - reference['energy_hartree']}
                    for state, reference in zip(row.get('correlated_energies_hartree', []),
                                                 far.get('correlated_energies_hartree', []))]}
                for row in group]})
    result = {'schema': 1, 'campaign': str(root), 'complete': len(completed), 'total': len(rows),
        'selection_status': 'user_review_required', 'jobs': rows,
        'radial_profiles': profiles,
        'equivalent_geometry_comparisons': comparisons,
        'limitations': [
            'A finite panel does not prove coverage of the VRC bottleneck.',
            'CASCI occupations and CI weights are screening diagnostics, not correlated error bounds.',
            'Compare interaction energies relative to consistent separated limits, not total energies.',
            'Root character, intruder sensitivity and rate sensitivity require review.',
            'Inner/outer state contributions, spin-orbit treatment and equilibration are user choices.']}
    output = _fresh(output)
    atomic_json(output, result)
    return output


def install_selection(input_file, selection_file, output, *, surface_file=None, bond=None, stage=None):
    """Copy a reviewed selection into an input or one explicit bottleneck.

    New files only. This registers a user's choice, never selects orbitals,
    invents electronic weights, or claims that a kinetics run has passed.
    """
    from kinbot.molpro import validate_vrc_state_model
    from kinbot.vrc_multistate import validate_surface_model
    input_file, selection_file = (Path(p).expanduser().resolve() for p in (input_file, selection_file))
    data = json.loads(input_file.read_text())
    selection = json.loads(selection_file.read_text())
    if selection.get('schema') != 2 or selection.get('state_average_policy') != 'reviewed':
        raise ValueError('Use a basis-specific selection with state_average_policy=reviewed.')
    space = selection['active_space']
    if space.get('mode') != 'selected' or not str(space.get('source', '')).strip():
        raise ValueError('A selected active space and review note are required.')
    models = selection.get('basis_models', {})
    if not models:
        raise ValueError('Selection has no basis-specific orbital models.')
    for model in models.values():
        validated = validate_vrc_state_model(model, active_electrons=space['active_electrons'],
                                             active_orbitals=space['active_orbitals'])
        if validated.get('orbital_selection', {}).get('mode') != 'molpro_record':
            raise ValueError('Final selection requires a hash-bound reviewed WFU.')
    # Keep pre-existing input-relative external files valid in the new location.
    for key in ('vrc_tst_state_model_file', 'vrc_tst_surface_model_file'):
        if data.get(key):
            data[key] = str((input_file.parent / data[key]).resolve())
    updates = {'vrc_tst_active_space': space, 'vrc_tst_state_model_file': str(selection_file)}
    output = _fresh(output)
    if surface_file is None:
        if bond is not None or stage is not None:
            raise ValueError('A stage selection requires surface_file, bond and stage.')
        if data.get('vrc_tst_surface_model_file'):
            raise ValueError('Input has a surface model; select an explicit bond/stage.')
        data.update(updates)
        data.pop('vrc_tst_state_model', None)
        data.pop('vrc_tst_state_average_roots', None)
        data.pop('vrc_tst_selected_roots', None)
    else:
        if bond is None or stage is None:
            raise ValueError('A stage selection requires bond and stage.')
        surface_file = Path(surface_file).expanduser().resolve()
        surface = validate_surface_model(json.loads(surface_file.read_text()), reviewed=False)
        matches = []
        for channel in surface['channels']:
            for bottleneck in channel['bottlenecks']:
                par = bottleneck['parameters']
                if par.get('vrc_tst_state_model_file'):
                    par['vrc_tst_state_model_file'] = str((surface_file.parent / par['vrc_tst_state_model_file']).resolve())
                if sorted(channel['bond']) == sorted(bond) and bottleneck['name'] == stage:
                    matches.append(bottleneck)
        if len(matches) != 1:
            raise ValueError('Selection must match exactly one declared bond/stage.')
        bottleneck = matches[0]
        for model in models.values():
            counts = {entry['multiplicity']: entry['roots'] for entry in model['state_average']}
            if any(state['root'] > counts.get(state['multiplicity'], 0) for state in bottleneck['states']):
                raise ValueError('Selected average cannot supply a requested bottleneck surface.')
        bottleneck['parameters'].update(updates)
        for key in ('vrc_tst_state_model', 'vrc_tst_state_average_roots', 'vrc_tst_selected_roots'):
            bottleneck['parameters'].pop(key, None)
        # Empty leaf overrides prevent unrelated global root selections/counts
        # leaking into this independently reviewed region.
        bottleneck['parameters']['vrc_tst_state_average_roots'] = {}
        bottleneck['parameters']['vrc_tst_selected_roots'] = {}
        surface_output = _fresh(output.with_suffix('.surfaces.json'))
        atomic_json(surface_output, surface)
        data['vrc_tst_surface_model_file'] = str(surface_output)
    atomic_json(output, data)
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(prog='python -m kinbot.orbital_review')
    commands = parser.add_subparsers(dest='command', required=True)
    for name in ('plan', 'prepare', 'report'):
        command = commands.add_parser(name)
        command.add_argument('source')
        command.add_argument('output')
        if name == 'prepare':
            command.add_argument('--additional-plan', action='append', default=[])
            command.add_argument('--reuse-from')
    install = commands.add_parser('install-selection')
    install.add_argument('input')
    install.add_argument('selection')
    install.add_argument('output')
    install.add_argument('--surface-file')
    install.add_argument('--bond', type=int, nargs=2)
    install.add_argument('--stage', choices=['single', 'inner', 'outer'])
    args = parser.parse_args(argv)
    if args.command == 'install-selection':
        print(install_selection(args.input, args.selection, args.output,
              surface_file=args.surface_file, bond=args.bond, stage=args.stage))
    elif args.command == 'prepare':
        print(prepare(args.source, args.output, additional_plans=args.additional_plan, reuse_from=args.reuse_from))
    else:
        print({'plan': plan, 'prepare': prepare, 'report': report}[args.command](args.source, args.output))


if __name__ == '__main__':
    main()
