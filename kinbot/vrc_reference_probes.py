"""Prepare finite native probes from immutable accepted VRC scan artifacts."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from ase.data import atomic_numbers
from ase.io import read as read_atoms

from kinbot.execution import file_sha256
from kinbot.molpro import continued_vrc_state_model, vrc_multireference_energy_contract
from kinbot.vrc_correction import correction_components, spin_splitting_settings
from kinbot.vrc_multistate import load_surface_model, state_parameters


def read(path):
    return json.loads(Path(path).read_text())


def context_fragment_multiplicities(parameters, saved):
    """Resolve explicit fragment states without waiting for ANL finalization.

    Use the same geometry perception/state-model functions as production, but
    refuse undeclared spins or charged fragmentation rather than guessing.
    """
    from kinbot.electronic_states import apply_electronic_state_model
    from kinbot.stationary_pt import StationaryPoint
    from kinbot.vrc_scan_recovery import digest
    context = saved['context']
    if digest(context) != saved.get('context_sha256'):
        raise ValueError('Fragment context hash is missing or changed.')
    if context['charge'] != 0 or context.get('fragment_geometry') is None:
        raise ValueError('Need explicit fragment records for charged/unseparated channels.')
    parent = StationaryPoint('reference_fragments', context['charge'], context['multiplicity'],
        atom=np.array(context['symbols']), geom=np.array(context['fragment_geometry']))
    parent.characterize()
    products, maps = parent.start_multi_molecular()
    actual = {tuple(sorted(map(int, m))): p for p, m in zip(products, maps)}
    expected = [tuple(sorted(m)) for m in context['maps']]
    if len(products) != 2 or len(actual) != 2 or set(actual) != set(expected):
        raise ValueError('Separated fragment connectivity differs from the saved atom maps.')
    result = []
    for mapping in expected:
        product = actual[mapping]
        if apply_electronic_state_model(product, parameters) is None:
            raise ValueError('No explicit electronic-state model for fragment '+product.smiles)
        result.append(int(product.mult))
    return result


def leaf_parameters(well, leaf, saved=None):
    """Resolve actual accepted products and the current declared scalar leaf."""
    source = well/(well.name+'.json')
    par = read(source)
    par['input_file'] = str(source)
    override = par.get('vrc_tst_surface_models', {}).get(well.name)
    if override:
        par.update(vrc_tst_surface_model=override, vrc_tst_surface_model_file='')
    match = re.fullmatch(r'(.+_hom_sci_(\d+)_(\d+))__(single|inner|outer)_m(\d+)_r(\d+)', leaf)
    if not match:
        raise ValueError('No unambiguous scalar leaf identity: '+leaf)
    reaction, a, b, bottleneck, mult, root = match.groups()
    requested = well/'vrctst/orbital_review'/reaction/'request.json'
    if requested.is_file():
        fragments = [p['multiplicity'] for p in read(requested)['products']]
    elif (well/'anl_automatic/network_request.json').is_file():
        network = read(well/'anl_automatic/network_request.json')['network']
        rx = next(r for r in network['reactions'] if r['name'] == reaction)
        species = {s['routing_name']: s for s in network['stable_species']}
        fragments = [species[p]['molecule']['multiplicity'] for p in rx['products']]
    else:
        fragments = context_fragment_multiplicities(par, saved or {})
    surface = load_surface_model(par)
    if surface is None:
        raise ValueError('No declared production surface model.')
    channel = next(c for c in surface[0]['channels'] if sorted(c['bond']) == sorted([int(a), int(b)]))
    stage = next(s for s in channel['bottlenecks'] if s['name'] == bottleneck)
    state = next(s for s in stage['states'] if (s['multiplicity'], s['root']) == (int(mult), int(root)))
    return state_parameters(par, stage, state, fragments), fragments, int(mult), int(root)


def prepare_configs(production, template):
    root = Path(production).resolve(strict=True)
    configs, skipped = [], []
    for request in sorted(root.glob('*/[0-9]*/vrctst/scan_review/*/scan_review.json')):
        saved = read(request)
        well = request.parents[3]
        leaf = saved['reaction']
        points = {**saved.get('history', {}), **saved.get('points', {})}
        accepted = [p for p in points.values() if p.get('status') == 'accepted' and p.get('graph')]
        far = saved['context']['separation_angstrom']
        inner = [p for p in accepted if p['distance_angstrom'] < far]
        chosen = ([min(inner, key=lambda p: (abs(p['distance_angstrom']-3.5), p['distance_angstrom']))]
                  if inner else [])
        chosen += [p for p in accepted if p['distance_angstrom'] == far]
        # Older/in-progress scans have not yet recorded the final grid's far
        # point, but already bind their optimized fragment reference in context.
        if len(chosen) == 1:
            try:
                selector = saved['context']['state_model']['orbital_selection']
                wfu = Path(selector['path']).resolve(strict=True)
                graph = wfu.parents[2]
                if (not graph.is_relative_to(root) or '_vts_fragments_asymptote_' not in graph.name
                        or file_sha256(wfu) != selector['sha256']):
                    raise ValueError('Invalid bound fragment checkpoint.')
                directory = graph/'tasks/relaxed_scan'
                atoms = read_atoms(directory/'final.xyz')
                if not np.allclose(atoms.positions, saved['context']['fragment_geometry'], atol=1.e-7, rtol=0.):
                    raise ValueError('Bound fragment geometry differs from final.xyz.')
                chosen.append(dict(distance_angstrom=far, graph=str(graph),
                    source='accepted_fragment_context', hashes={str(p):file_sha256(p) for p in
                    (graph/'workflow.json', directory/'execution.json', directory/'final.xyz')}))
            except (KeyError, OSError, ValueError, TypeError, IndexError):
                pass
        if len(chosen) != 2:
            skipped.append(dict(channel=leaf, reason='Need accepted interior and separated-reference points.'))
            continue
        try:
            par, fragments, target, target_root = leaf_parameters(well, leaf, saved)
            if par.get('vrc_tst_correction_model') != 'spin_splitting':
                raise ValueError('This finite CC/MR panel requires a spin-splitting production contract.')
            spin = spin_splitting_settings(par, parent_multiplicity=target, fragment_multiplicities=fragments)
            reference = spin['reference_multiplicity']
        except (ValueError, KeyError, OSError, StopIteration) as exc:
            skipped.append(dict(channel=leaf, reason='Cannot resolve live contract: '+str(exc)))
            continue
        for point in chosen:
            try:
                graph = Path(point['graph']).resolve(strict=True)
                directory = graph/'tasks/relaxed_scan'
                if not graph.is_relative_to(root) or (directory/'scan_rejection.json').exists():
                    raise ValueError('Graph is outside run or marked rejected.')
                hashes = point.get('hashes', {})
                required = [graph/'workflow.json', directory/'execution.json', directory/'final.xyz']
                if any(hashes.get(str(p)) != file_sha256(p) for p in required):
                    raise ValueError('Accepted-point artifact hashes missing or changed.')
                execution = read(directory/'execution.json')
                if execution['status'] != 'executed':
                    raise ValueError('Production calculation is not complete.')
                task = read(directory/'task.json')['task']
                profile = task['profile']; options = profile['calculator_kwargs']
                molecule = read(graph/'workflow.json')['molecule']
                atoms = read_atoms(directory/'final.xyz')
                if atoms.get_chemical_symbols() != molecule['symbols']:
                    raise ValueError('Final atom order differs from the native molecule.')
                generated = execution['details']['generated_files']
                inputs = [directory/Path(n).name for n in generated if re.search(r'_step_\d+\.inp$', n)]
                inp = inputs[-1]
                last = read_atoms(inp.with_suffix('.xyz'))
                if (last.get_chemical_symbols() != atoms.get_chemical_symbols()
                        or not np.allclose(last.positions, atoms.positions, atol=1.e-7, rtol=0.)):
                    raise ValueError('Final force and accepted geometry differ.')
                for p in (directory/'task.json', inp, inp.with_suffix('.out'), inp.with_suffix('.xyz')):
                    if execution.get('artifacts', {}).get(p.name) != file_sha256(p):
                        raise ValueError('Native execution artifact missing or changed: '+p.name)
                contract = vrc_multireference_energy_contract(profile['method'])
                ne, no = contract['active_electrons'], contract['active_orbitals']
                electrons = sum(atomic_numbers[s] for s in molecule['symbols'])-molecule['charge']
                model = execution['details'].get('orbital_checkpoint_model')
                if model is None:
                    wfu = [directory/Path(n).name for n in generated if n.endswith('.wfu')][-1]
                    if execution['artifacts'].get(wfu.name) != file_sha256(wfu):
                        raise ValueError('Final checkpoint hash missing or changed.')
                    model = continued_vrc_state_model(options['state_model'], wfu, nelectron=electrons,
                                                     active_electrons=ne, active_orbitals=no)
                selector = model['orbital_selection']
                checkpoint = Path(selector['path']).resolve(strict=True)
                closed = (electrons-ne)//2
                if (not checkpoint.is_relative_to(directory) or file_sha256(checkpoint) != selector['sha256']
                        or selector['record'] != '2140.2' or selector['orbital_type'] != 'natural'
                        or selector['active_indices'] != list(range(closed+1, closed+no+1))):
                    raise ValueError('Expected a hash-bound, already packed final natural checkpoint.')
                roots = {e['multiplicity']: e['roots'] for e in model['state_average']}
                components = correction_components(par, symbols=molecule['symbols'], charge=molecule['charge'],
                                                   settings=spin, state_roots=roots)
                basis = profile['basis']
                cc = [c for c in components if c['role'].startswith('reference:') and c['basis'] == basis]
                mr = [c for c in components if c['role'].startswith('gap:')
                      and ':reference:' in c['role'] and c['basis'] == basis]
                if len(cc) != 1 or len(mr) != 1 or cc[0]['method'].lower() != 'uccsd(t)-f12b':
                    raise ValueError('Expected one production CC/MR reference at the saved scan basis.')
                other = json.loads(json.dumps(mr[0]['state_model']))
                if (mr[0]['method'].lower() != profile['method'].lower()
                        or any(model.get(k) != other.get(k) for k in
                               ('state_average', 'selected_roots', 'energy_selection', 'correlation_options',
                                'casscf_options', 'canonical_orbitals'))):
                    raise ValueError('Current correction contract differs from saved scan model; inspect before replay.')
                shift, ipea = mr[0].get('shift'), mr[0].get('ipea')
                if (shift != options.get('caspt2_shift') or ipea != options.get('caspt2_ipea')):
                    raise ValueError('Correction shifts differ from saved scan profile.')
                data = {k: deepcopy(v) for k, v in template.items() if not k.startswith('_') and k not in
                        ('coordinates', 'anchors', 'scan_policy', 'orientation_probes', 'initial_records',
                         'probe_definition', 'probe_space', 'casscf_options', 'correlation_checks')}
                title = 'production_'+hashlib.sha256(str(graph).encode()).hexdigest()[:16]
                clean = deepcopy(model); clean.pop('orbital_selection')
                data.update(title=title, symbols=molecule['symbols'], geometry_angstrom=atoms.positions.tolist(),
                    charge=molecule['charge'], multiplicity=target, probe_method='casscf', probe_active_space=True,
                    orbital_initialization='independent', bases=[dict(name='native', basis=basis)],
                    anchors=[dict(name='point', geometry_angstrom=atoms.positions.tolist())],
                    initial_records={'native__point':dict(path=str(checkpoint),sha256=selector['sha256'],basis=basis)},
                    scf_convergence=options['scf_convergence'], scf_options={},
                    correlation_convergence=options['correlation_convergence'],
                    casscf_options=clean.get('casscf_options', {}), reference_identity_check=True,
                    correlation_checks=[dict(method=profile['method'],shift=shift,ipea=ipea,
                        reference_policy=clean.get('correlation_options', {}).get('reference_policy','selected_root')),
                        dict(method=cc[0]['method'], multiplicity=reference)],
                    probe_definition=dict(active_electrons=ne,active_orbitals=no,state_model=clean,
                        chemical_note='Native fixed-geometry production reference check; preserve the saved CAS and SA. '
                                      'Use the final packed checkpoint without repeating initial orbital swaps.',
                        orbital_partition=dict(source_active_electrons=ne,source_active_orbitals=no,
                                               close_lowest=0,discard_highest=0)))
                from kinbot.vrc_reference_review import contract, digest
                data['reference_review_contract_sha256'] = digest(contract(par))
                data['reference_review_distance_angstrom'] = point['distance_angstrom']
                metadata = dict(channel=leaf, graph=str(graph), distance_angstrom=point['distance_angstrom'],
                    target_multiplicity=target, target_root=target_root, reference_multiplicity=reference,
                    reference_root=model['selected_roots'].get(str(reference), model['selected_roots'].get(reference)),
                    reference_selection=model.get('energy_selection', {}).get(str(reference),
                                          model.get('energy_selection', {}).get(reference, 'selected_root')),
                    production_review_status=saved.get('status'), source_hashes=hashes,
                    profile=profile, checkpoint=selector)
                configs.append((data, metadata))
            except (ValueError, KeyError, OSError, IndexError, TypeError) as exc:
                skipped.append(dict(channel=leaf, distance_angstrom=point['distance_angstrom'], reason=str(exc)))
    return configs, skipped
