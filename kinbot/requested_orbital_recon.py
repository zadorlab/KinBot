"""General finite orbital screens for existing or provisional VRC requests.

Scientific choices live in the supplied recipe. Only the explicit ``review
--apply`` command installs a human-reviewed selection in an existing PES
input. It uses live refresh and never starts/restarts a production driver.
"""
import argparse
from contextlib import ExitStack, contextmanager
from copy import deepcopy
import fcntl
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

import numpy as np

from kinbot import molpro_reconnaissance as recon, run_policy
from kinbot.execution import atomic_json, file_sha256


def _read(path):
    return json.loads(Path(path).read_text())


@contextmanager
def _new_campaign(output):
    """Roll back only this invocation's unsubmitted staging on a late failure.

    mkdir is the ownership claim: a competing preparer cannot remove an
    existing campaign. Preparation performs no scheduler/native calls.
    """
    output.mkdir(parents=True, exist_ok=False)
    try:
        yield
    except BaseException:
        import shutil
        shutil.rmtree(output)
        raise


def _require_available_recipe(payload):
    # An absent status retains compatibility with explicitly defined recipes.
    # Withdrawn/proposed recipes can preserve evidence without authorizing QC.
    if payload.get('status', 'ready') != 'ready':
        raise ValueError('Orbital review recipe/plan is not ready: '
                         + str(payload.get('status')) + '. '
                         + str(payload.get('status_reason', 'Inspect its scientific selection.')))


def _saved(root, *, allow_accepted=False):
    record = _read(root/'review_launch.json')
    for name, digest in record['files'].items():
        if file_sha256(root/name) != digest:
            raise ValueError(f'Prepared evidence changed: {name}.')
    if record.get('reviews'):
        for child in record['reviews']:
            saved = _saved(Path(child), allow_accepted=allow_accepted)
            if (saved['production_policy'] != record['production_policy']
                    or saved.get('shared_pool', False) != record.get('shared_pool', False)):
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
        for key in ('surface_model', 'installed_channel_sha256', 'acknowledgement'):
            live.pop(key, None)
            original.pop(key, None)
        live['status'] = original.get('status')
    if live != original:
        raise ValueError('Live review request changed; inspect before continuing this review.')
    return record


def _review_policy(request, output, resource_root=None, *, preflight=False):
    """Separate immutable scientific evidence from a new scheduling owner."""
    source_policy = run_policy.find(Path(request).parent)
    if resource_root is None:
        return source_policy, run_policy.read(Path(request).parent) or {}
    root = Path(resource_root).expanduser().resolve(strict=True)
    policy = run_policy.find(root)
    if policy != root/run_policy.NAME:
        raise ValueError('The explicit review resource root needs its own scheduling policy.')
    if source_policy is None or (root.is_relative_to(source_policy.parent)
                                or source_policy.parent.is_relative_to(root)):
        raise ValueError('A detached review pool must be disjoint from its evidence policy tree.')
    if not preflight and (not output.is_relative_to(root) or output == root
                         or run_policy.find(output.parent) != policy):
        raise ValueError('Detached reviews must stay beneath their explicit resource root.')
    return policy, run_policy.read(root)


def candidates(source, recipe, output, *, preflight=False, resume=False, resource_root=None):
    """Prepare compact CASSCF or fixed-orbital CASCI from validated seeds.

    This does not select orbitals by occupation or accept a production model.
    Every partition, seed geometry and native check comes from the recipe.
    """
    source = Path(source).expanduser().resolve(strict=True)
    output = Path(output).expanduser().resolve()
    previous = _saved(source)
    controls = _read(recipe)
    _require_available_recipe(controls)
    request = _read(source/'request.json')
    if controls.get('schema') != 1 or not controls.get('request_match'):
        raise ValueError('Candidate recipe must identify its intended review request.')
    for key, value in controls['request_match'].items():
        if request.get(key) != value:
            raise ValueError(f'Request differs from candidate recipe: {key}.')
    policy_path, policy = _review_policy(previous['request'], output, resource_root)
    shared = resource_root is not None or previous.get('shared_pool', False)
    slots = (policy.get('configured_max_jobs', policy.get('max_jobs', 0)) if shared
             else policy.get('orbital_review_reserved_slots', 0))
    cap = controls.get('max_jobs', previous['max_jobs'])
    if type(cap) is not int or cap < 1 or slots < 1:
        raise ValueError('Candidate review needs positive reserved capacity.')
    if output.exists() and resume:
        saved = _saved(output)
        if (saved.get('source_review') != str(source) or _read(output/'recipe.json') != controls
                or saved.get('shared_pool', False) != shared
                or saved['production_policy'] != str(policy_path)):
            raise ValueError('Existing candidate plan or ownership differs.')
        for system in saved['systems']:
            loaded = recon._load(system['configuration'])
            for seed in loaded['initial_records'].values():
                if file_sha256(seed['path']) != seed['sha256']:
                    raise ValueError('Candidate source checkpoint changed.')
        return dict(status='resume_ready', directory=str(output), models=len(saved['systems']),
                    tasks=sum(len(_read(s['manifest'])['jobs']) for s in saved['systems']))
    if output.exists() or output == source or source in output.parents:
        raise FileExistsError('Preserve the source review; choose a fresh sibling directory.')
    if shared:
        if run_policy.find(output.parent) != policy_path:
            raise ValueError('Shared candidates must remain in the same production policy tree.')
    elif policy_path.parent == output or policy_path.parent in output.parents:
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
    scf_source = not bool(loaded.get('_probe_space'))
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
            if scf_source:
                if definition.get('orbital_source') != 'scf' or partition:
                    raise ValueError('SCF candidates need explicit orbital_source=scf and selected indices.')
            else:
                for key in ('active_electrons', 'active_orbitals'):
                    if partition.get('source_'+key) != loaded['_probe_space'][key]:
                        raise ValueError('Partition counts differ from the source active space.')
            initialization = item.get('orbital_initialization', 'independent')
            method = item.get('probe_method', 'casscf')
            if method not in ('casscf', 'casci'):
                raise ValueError('Candidate probe_method must be casscf or casci.')
            if method == 'casci' and (scf_source or initialization != 'independent'):
                raise ValueError('Seeded CASCI diagnostics require independent natural checkpoints.')
            if initialization not in ('independent', 'propagated') or (
                    scf_source and initialization != 'independent'):
                raise ValueError('Propagation requires a completed natural-orbital source.')
            data.update(title=name, probe_method=method, orbital_initialization=initialization,
                        probe_active_space=True, orbital_character_only=False,
                        probe_definition=definition, anchors=[], initial_records={},
                        casscf_options=deepcopy(item['casscf_options']),
                        reference_identity_check=item.get('reference_identity_check', False),
                        correlation_checks=deepcopy(item['correlation_checks']))
            for key in ('scf_convergence', 'correlation_convergence', 'scf_options'):
                if key in item:
                    data[key] = deepcopy(item[key])
            chains = {}
            for choice in item['anchors']:
                if set(choice) - {'source', 'name', 'distance_angstrom', 'rotation',
                                  'seed', 'active_indices', 'seed_sha256',
                                  'propagate', 'orbital_chain'}:
                    raise ValueError('Unsupported candidate anchor field.')
                if not scf_source and 'active_indices' in choice:
                    raise ValueError('Explicit SCF indices cannot alter a natural-orbital partition.')
                propagate = choice.get('propagate', False)
                if type(propagate) is not bool or (propagate and (
                        initialization != 'propagated' or
                        set(choice) & {'seed', 'active_indices', 'seed_sha256'})):
                    raise ValueError('A propagated anchor must consume its predecessor, not another seed.')
                anchor = deepcopy(anchors[choice['source']])
                anchor['name'] = choice.get('name', choice['source'])
                if 'orbital_chain' in choice:
                    anchor['orbital_chain'] = choice['orbital_chain']
                if initialization == 'propagated':
                    chain = anchor.get('orbital_chain')
                    if not isinstance(chain, str) or not recon._SAFE.fullmatch(chain):
                        raise ValueError('Propagated candidates require explicit safe orbital_chain names.')
                    # A chain follows one mapped fragment separation. Geometry
                    # changes remain explicit; never silently join two exits.
                    identity = {k: anchor.get(k) for k in
                                ('bond', 'moving_atoms', '_channel', '_context', '_rotation')}
                    if propagate and chains.get(chain) != identity:
                        raise ValueError('A propagated anchor needs a preceding seed on the same coordinate.')
                    if not propagate and chain in chains:
                        raise ValueError('An orbital chain has one seed; use a distinct chain for another start.')
                    if not propagate and not {'seed', 'seed_sha256'} <= set(choice):
                        raise ValueError('A continuation seed needs explicit job and basis-specific hashes.')
                    if 'rotation' in choice:
                        raise ValueError('Use independent candidates for new fragment orientations.')
                    chains[chain] = identity
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
                if propagate:
                    continue
                for basis in data['bases']:
                    old = jobs[basis['name']+'__'+choice.get('seed', choice['source'])]
                    if method == 'casci' and not np.allclose(
                            old['geometry_angstrom'], recon._anchor_geometry(data, anchor),
                            atol=1.e-12, rtol=0.):
                        raise ValueError('Seeded CASCI diagnostics require the source geometry unchanged.')
                    wrong_space = (not scf_source and any(
                        old['probe_active_space'][key] != loaded['_probe_space'][key]
                        for key in ('active_electrons', 'active_orbitals')))
                    if (wrong_space or old.get('validation_errors')
                            or not scf_source and not old.get('natural_orbital_probe')
                            or old['basis'] != basis['basis']
                            or file_sha256(old['wfu']) != old['wfu_sha256']
                            or file_sha256(Path(old['directory'])/(old['name']+'.out')) != old['output_sha256']):
                        raise ValueError('Source checkpoint/output or model differs from validated evidence.')
                    if 'seed_sha256' in choice and choice['seed_sha256'].get(basis['name']) != old['wfu_sha256']:
                        raise ValueError('Source checkpoint differs from the explicitly reviewed seed hash.')
                    data['initial_records'][basis['name']+'__'+anchor['name']] = dict(
                        path=old['wfu'], sha256=old['wfu_sha256'], basis=old['basis'])
                    if scf_source:
                        inp = Path(old['directory'])/(old['name']+'.inp')
                        if (old.get('validation_errors') or old.get('probe_active_space')
                                or file_sha256(inp) != old['input_sha256']):
                            raise ValueError('Source is not an intact validated SCF probe.')
                        if choice.get('seed_sha256', {}).get(basis['name']) != old['wfu_sha256']:
                            raise ValueError('SCF source checkpoint differs from the explicitly reviewed seed hash.')
                        data['initial_records'][basis['name']+'__'+anchor['name']].update(
                            active_indices=choice.get('active_indices', {}).get(basis['name']),
                            source_output=str(Path(old['directory'])/(old['name']+'.out')),
                            source_output_sha256=old['output_sha256'], source_input=str(inp),
                            source_input_sha256=old['input_sha256'])
            path = Path(tmp)/(name+'.json')
            atomic_json(path, data)
            check = recon._load(path)
            for basis in check['bases']:
                predecessors = {}
                for anchor in check['anchors']:
                    stem = basis['name']+'__'+anchor['name']
                    chain = anchor.get('orbital_chain')
                    orbital_start = check['initial_records'].get(stem) or predecessors.get(chain)
                    recon.validate_input_blocks(recon._deck(check, basis['basis'], anchor, stem,
                        previous=orbital_start), stem)
                    if initialization == 'propagated':
                        predecessors[chain] = stem
            configurations.append(data)
    if preflight:
        return dict(status='preflight_passed', models=len(configurations),
                    tasks=sum(len(d['anchors'])*len(d['bases']) for d in configurations))
    with _new_campaign(output):
        for name in ('request.json', 'electronic_context.json'):
            atomic_json(output/name, _read(source/name))
        atomic_json(output/'recipe.json', controls)
        atomic_json(output/'source_manifest.json', original)
        files = ['request.json', 'electronic_context.json', 'recipe.json', 'source_manifest.json']
        if (source/'geometry_source.json').exists():
            atomic_json(output/'geometry_source.json', _read(source/'geometry_source.json'))
            files.append('geometry_source.json')
        systems = []
        for data in configurations:
            config = output/(data['title']+'.json')
            data['kinbot_input'] = str(output/'electronic_context.json')
            atomic_json(config, data)
            manifest = recon.prepare(config, output/data['title'])
            state = _read(manifest)
            state.update(slurm_excluded_nodes=policy['excluded_nodes'], startup_retries=2)
            if shared:
                state.update(shared_policy_file=str(policy_path), slurm_excluded_nodes=[])
            atomic_json(manifest, state)
            files.append(config.name)
            systems.append(dict(title=data['title'], configuration=str(config), manifest=str(manifest)))
        record = dict(schema=1, status='prepared', shared_pool=shared, request=previous['request'],
                      max_jobs=min(cap, slots), production_policy=str(policy_path),
                      files={name: file_sha256(output/name) for name in files}, systems=systems,
                      candidate_comparison=True, source_review=str(source))
        if resource_root is not None:
            record['resource_root'] = str(policy_path.parent)
        atomic_json(output/'review_launch.json', record)
        return dict(directory=str(output), models=len(systems), max_jobs=record['max_jobs'],
                    tasks=sum(len(_read(s['manifest'])['jobs']) for s in systems),
                    production_model_accepted=False)


def product_complex_geometry(request, evidence, selection):
    """Read an accepted path's mapped product geometry for exit-channel probes."""
    import sqlite3
    import networkx as nx
    from ase.db import connect
    from rdkit import Chem
    from kinbot.stationary_pt import StationaryPoint
    worker=Path(request).parents[3]
    network=_read(worker/'accepted_l2_network.json')['network']['reactions']
    rows=[r for r in network if r['name']==evidence['reaction']]
    if len(rows)!=1 or rows[0].get('product_complex')!=selection['product_complex']:
        raise ValueError('Product-side screen is not bound to the accepted reaction complex.')
    name=selection['record_name']
    if name not in [evidence['reaction']+s for s in ('_IRC_F_prod_high','_IRC_R_prod_high')]:
        raise ValueError('Product-side screen requires an explicit L2 product record.')
    path=worker/'kinbot.db'
    db=connect(str(path),use_lock_file=False)
    db.connection=sqlite3.connect(path.as_uri()+'?mode=ro',uri=True);db.change_count=0
    try:
        found=list(db.select(name=name))
        if not found:raise ValueError('Missing saved L2 product geometry: '+name)
        row=found[-1];atoms=row.toatoms()
    finally:db.connection.close()
    if row.data.get('status')!='normal':
        raise ValueError('Saved L2 product geometry has no normal completion: '+name)
    if atoms.get_chemical_symbols()!=evidence['symbols'] or not np.isfinite(atoms.positions).all():
        raise ValueError('Product-side atom identity/coordinates differ.')
    fragments=selection['fragments'];maps=[f['atoms'] for f in fragments]
    if len(maps)!=2 or sorted(i for m in maps for i in m)!=list(range(1,len(atoms)+1)):
        raise ValueError('Explicit one-based fragment maps must partition all atoms.')
    species=StationaryPoint('product',evidence['charge'],evidence['multiplicity'],
        atom=np.array(evidence['symbols']),geom=atoms.positions)
    species.characterize()
    parts=list(nx.connected_components(nx.from_numpy_array(species.bond01)))
    if {frozenset(p) for p in parts}!={frozenset(i-1 for i in m) for m in maps}:
        raise ValueError('Saved product geometry does not match the declared separated fragments.')
    from kinbot.electronic_states import _canonical_smiles
    product_key=lambda p:(_canonical_smiles(p['smiles']),p['charge'],p['multiplicity'])
    if sorted(map(product_key,fragments))!=sorted(map(product_key,evidence['products'])):
        raise ValueError('Mapped fragment states differ from the live review request.')
    for fragment,mapping in zip(fragments,maps):
        mol=Chem.AddHs(Chem.MolFromSmiles(fragment['smiles']))
        expected=nx.Graph();expected.add_nodes_from((a.GetIdx(),dict(element=a.GetSymbol())) for a in mol.GetAtoms())
        expected.add_edges_from((b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in mol.GetBonds())
        actual=nx.from_numpy_array(species.bond01).subgraph([i-1 for i in mapping]).copy()
        nx.set_node_attributes(actual,{i:atoms[i].symbol for i in actual},'element')
        if not nx.is_isomorphic(actual,expected,node_match=lambda a,b:a['element']==b['element']):
            raise ValueError('Mapped fragment connectivity differs from its declared identity.')
    pair=evidence['bond'];moving=next((m for m in maps if pair[1] in m and pair[0] not in m),None)
    if moving is None:raise ValueError('VRC coordinate does not separate the mapped fragments.')
    anchors=[dict(name='exit_r_'+str(d).replace('.','p'),bond=pair,moving_atoms=moving,
                  distance_angstrom=d) for d in selection['distances_angstrom']]
    include_reference = selection.get('include_reference_geometry', False)
    if not isinstance(include_reference, bool):
        raise ValueError('include_reference_geometry must be boolean.')
    if include_reference:
        distance = float(np.linalg.norm(atoms.positions[pair[0]-1]-atoms.positions[pair[1]-1]))
        anchors.insert(0, dict(name='exit_l2', bond=pair, moving_atoms=moving,
                              distance_angstrom=distance))
    proof=dict(database=str(path),row_id=row.id,record_name=name,network=rows[0],
               symbols=evidence['symbols'],positions=atoms.positions.tolist(),fragments=fragments)
    return atoms.positions.tolist(),anchors,proof


def prepare(request, case_input, recipe, output, *, provisional=False, shared_pool=False,
            _preflight=False, resource_root=None):
    request, case_input, recipe = [Path(p).expanduser().resolve(strict=True)
                                    for p in (request, case_input, recipe)]
    output = Path(output).expanduser().resolve()
    evidence, controls = _read(request), _read(recipe)
    _require_available_recipe(controls)
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
    from kinbot.reaction_decisions import rejected
    worker = request.parent if provisional else request.parents[3]
    if rejected(worker, evidence['reaction']):
        raise ValueError('An operator-rejected reaction cannot start a new orbital screen.')
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
    if resource_root is not None and not shared_pool:
        raise ValueError('An explicit resource root requires shared-pool scheduling.')
    policy_file, policy = _review_policy(request, output, resource_root, preflight=_preflight)
    slots = (policy.get('configured_max_jobs', policy.get('max_jobs', 0)) if shared_pool
             else policy.get('orbital_review_reserved_slots', 0))
    if slots < 1:
        raise ValueError('No orbital-review capacity is configured/currently reserved.')
    requested_cap = controls.get('max_jobs', slots)
    if type(requested_cap) is not int or requested_cap < 1:
        raise ValueError('Recipe max_jobs must be a positive integer.')
    inside = policy_file.parent in output.parents
    if shared_pool:
        if not _preflight and (not inside or run_policy.find(output.parent) != policy_file):
            raise ValueError('Shared-pool reviews must be inside their production policy tree.')
    elif policy_file.parent == output or inside:
        raise ValueError('Keep reconnaissance outside the production run directory.')
    context = _read(case_input)
    config = deepcopy(controls['reconnaissance'])
    recon.require_explicit_probe_space(config)
    config.update(symbols=evidence['symbols'], geometry_angstrom=evidence['positions'],
                  charge=evidence['charge'], partition=context['queue_name'])
    geometry_proof=None
    if controls.get('geometry_source'):
        if provisional:raise ValueError('Product-side screens require an accepted live request.')
        positions,anchors,geometry_proof=product_complex_geometry(request,evidence,controls['geometry_source'])
        config.pop('coordinates',None)
        config.update(geometry_angstrom=positions,anchors=anchors)
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
    with _new_campaign(output):
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
        if shared_pool:
            state['shared_policy_file'] = str(policy_file)
            state['slurm_excluded_nodes'] = []  # Pool exclusions are read at each admission.
        atomic_json(manifest, state)
        record = dict(schema=1, status='prepared', shared_pool=shared_pool, request=str(request), manifest=str(manifest),
                      max_jobs=min(requested_cap, slots), production_policy=str(policy_file),
                      files={name: file_sha256(output/name) for name in
                             ('request.json', 'recipe.json', 'electronic_context.json', 'input.json')},
                      systems=[dict(title=config['title'], configuration=str(output/'input.json'),
                                    manifest=str(manifest))])
        if resource_root is not None:
            record['resource_root'] = str(policy_file.parent)
        if provisional:
            record['provisional'] = True
            record['files']['provisional_source.json'] = file_sha256(output/'provisional_source.json')
        atomic_json(output/'review_launch.json', record)
        if geometry_proof:
            atomic_json(output/'geometry_source.json',geometry_proof)
            record['files']['geometry_source.json']=file_sha256(output/'geometry_source.json')
            atomic_json(output/'review_launch.json',record)
        return dict(directory=str(output), tasks=len(state['jobs']), max_jobs=record['max_jobs'],
                    active_space=loaded.get('_probe_space'), states=state['jobs'][0]['root_states'],
                    partition=config['partition'], excluded_nodes=policy['excluded_nodes'])


def batch(production, plan, output, *, shared_pool=False, preflight=False, resume=False,
          resource_root=None):
    """Stage explicit SCF, CASCI or seeded compact probes under one review cap.

    Paths and scientific settings are explicit in the plan/recipes. Existing
    review selections are skipped; no rejected or imaginary future channel is
    manufactured. A provisional source can never directly approve production.
    """
    production = Path(production).expanduser().resolve(strict=True)
    plan = Path(plan).expanduser().resolve(strict=True)
    output = Path(output).expanduser().resolve()
    if resource_root is not None and not shared_pool:
        raise ValueError('An explicit resource root requires shared-pool scheduling.')
    policy_file, _ = _review_policy(production/'input.json', output, resource_root)
    spec = _read(plan)
    _require_available_recipe(spec)
    entries = spec.get('screens', [])
    if (spec.get('schema') != 1 or not entries or type(spec.get('max_jobs')) is not int
            or spec['max_jobs'] < 1):
        raise ValueError('Batch plan requires schema 1, screens and positive max_jobs.')
    if output.exists():
        if resume:
            saved = _saved(output)
            if (_read(output/'plan.json') != spec or saved.get('shared_pool', False) != shared_pool
                    or saved['production_policy'] != str(policy_file)):
                raise ValueError('Existing batch plan or resource ownership differs.')
            for child in saved.get('reviews', []):
                child = Path(child)
                entry = next(e for e in entries if e['name'] == child.name)
                if _read(child/'recipe.json') != _read(plan.parent/entry['recipe']):
                    raise ValueError('Existing batch recipe changed: '+child.name)
                if 'source' in entry:
                    candidates(production/entry['source'], plan.parent/entry['recipe'], child,
                               resume=True, resource_root=resource_root)
            return dict(status='resume_ready', directory=str(output),
                        tasks=sum(len(_read(s['manifest'])['jobs']) for s in saved['systems']))
        raise FileExistsError('Batch already exists; use run to resume without duplicate preparation.')
    if shared_pool and resource_root is None:
        if production not in output.parents or run_policy.find(output.parent) != run_policy.find(production):
            raise ValueError('Shared-pool reviews must be inside their production policy tree.')
    elif production == output or production in output.parents:
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
        _require_available_recipe(control)
        candidate = 'source' in item
        source = (production/item['source']).resolve(strict=True) if candidate else None
        if candidate and (production not in source.parents
                          or not _saved(source).get('shared_pool') or not shared_pool):
            raise ValueError('Batch candidates require a completed source in the shared production tree.')
        config = control.get('reconnaissance', {})
        if not candidate and (config.get('correlation_checks')
                or (config.get('probe_active_space') and config.get('probe_method', 'casci') != 'casci')):
            raise ValueError('This batch command stages orbital-only SCF or CASCI screens.')
        match = control['request_match']
        parent, reaction = match['parent'], match['reaction']
        if not all(recon._SAFE.fullmatch(x) for x in (parent,reaction)):
            raise ValueError('Unsafe parent/reaction in batch recipe.')
        context = _read(case/'input.json')
        model = context.get('vrc_tst_surface_models', {}).get(parent, {})
        if any(sorted(c['bond'])==sorted(match['bond'])
               and c.get('reaction', reaction) == reaction for c in model.get('channels', [])):
            skipped.append(dict(name=name, reason='channel already supplied in case input'))
            continue
        request = case/parent/'vrctst/orbital_review'/reaction/'request.json'
        provisional = not request.exists()
        if provisional:
            if not item.get('allow_provisional', False):
                raise ValueError(f'{name}: no live orbital review request yet.')
            request = case/parent/(parent+'.json')
        policy_file, _ = _review_policy(case/'input.json', output, resource_root)
        if policy_file is None:
            raise ValueError('Every batch case needs a production resource policy.')
        policies.add(str(policy_file))
        if candidate and (provisional or Path(_saved(source)['request']) != request):
            raise ValueError('Candidate batch source belongs to a different review request.')
        pending.append((name,request,case/'input.json',recipe,provisional,source))
    if not pending:
        return dict(status='nothing_to_prepare', skipped=skipped, tasks=0)
    if len(policies)!=1:
        raise ValueError('A batch requires one shared production resource policy.')
    # Preflight the entire batch first, so a bad later recipe cannot publish a
    # partially usable campaign. This performs no scheduler/native calls.
    with tempfile.TemporaryDirectory(prefix='kinbot-batch-preflight-') as tmp:
        counts = []
        for name,request,case_input,recipe,provisional,source in pending:
            result = (candidates(source,recipe,output/name,preflight=True,
                                 resource_root=resource_root) if source else
                      prepare(request,case_input,recipe,Path(tmp)/name,provisional=provisional,
                              shared_pool=shared_pool,_preflight=True,resource_root=resource_root))
            counts.append(result['tasks'])
    if preflight:
        return dict(status='preflight_passed', screens=len(pending), tasks=sum(counts), skipped=skipped)
    with _new_campaign(output):
        systems, reviews, caps = [], [], []
        for name,request,case_input,recipe,provisional,source in pending:
            child = output/name
            if source:
                candidates(source,recipe,child,resource_root=resource_root)
            else:
                prepare(request,case_input,recipe,child,provisional=provisional,
                        shared_pool=shared_pool,resource_root=resource_root)
            saved = _saved(child)
            systems.extend(saved['systems']); reviews.append(str(child)); caps.append(saved['max_jobs'])
        atomic_json(output/'plan.json',spec)
        record = dict(schema=1,status='prepared', shared_pool=shared_pool, request=str(pending[0][1]),
            production_policy=next(iter(policies)),max_jobs=min(spec['max_jobs'],*caps),
            systems=systems,reviews=reviews,skipped=skipped,
            files={'plan.json':file_sha256(output/'plan.json')})
        if resource_root is not None:
            record['resource_root'] = str(policy_file.parent)
        atomic_json(output/'review_launch.json',record)
        return dict(directory=str(output),screens=len(reviews),skipped=skipped,
            tasks=sum(len(_read(s['manifest'])['jobs']) for s in systems),max_jobs=record['max_jobs'],
            production_model_accepted=False)


def bind(source, request, *, apply=True):
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
        if apply:
            atomic_json(source/'provisional_request.json',old)
            atomic_json(source/'request.json',current)
            record.update(provisional=False,request=str(request),
                          binding=dict(geometry_changed=old['positions']!=current['positions']))
            for name in ('provisional_request.json','request.json'):
                record['files'][name]=file_sha256(source/name)
            atomic_json(source/'review_launch.json',record)
    return dict(status='bound_for_review',request=str(request),production_model_accepted=False)


def _idle_panel(record):
    """Call with the campaign and all existing manifest locks held."""
    live = set(subprocess.run(['squeue', '-h', '-o', '%A'], capture_output=True,
        text=True, check=True, timeout=30).stdout.split())
    for system in record['systems']:
        manifest = _read(system['manifest'])
        for job in manifest['jobs'] + manifest.get('excluded_jobs', []):
            if str(job.get('job_id')) in live or job.get('submission_status') == 'submitting':
                raise RuntimeError('Wait for panel allocations/submissions to finish before repair or extension.')
    return live


def retry_reference_exports(output):
    """Archive only known density-selector failures; leave scientific failures alone."""
    from kinbot.molpro_reference_identity import repair_legacy_selectors
    root = Path(output).expanduser().resolve(strict=True)
    with recon._manifest_lock(root), ExitStack() as locks:
        record = _saved(root)
        paths = [Path(s['manifest']) for s in record['systems']]
        for parent in sorted({p.parent for p in paths}):
            locks.enter_context(recon._manifest_lock(parent))
        live = _idle_panel(record)
        selected = {}
        for path in paths:
            selected[path] = set()
            for job in _read(path)['jobs']:
                d = Path(job['directory']); inp = d/(job['name']+'.inp')
                if file_sha256(inp) != job['input_sha256']:
                    raise ValueError('Prepared native input changed: '+str(inp))
                _, repairs = repair_legacy_selectors(inp.read_text())
                receipt = d/'execution.json'
                if not repairs or not receipt.exists() or not job.get('job_id'):
                    continue
                execution = _read(receipt)
                out = d/(job['name']+'.out')
                text = out.read_text(errors='replace') if out.exists() else ''
                if (execution.get('status') == 'failed'
                        and re.search(r'REQUESTED DENSITY NOT FOUND IN RECORD\s+300\d+\.2', text)):
                    selected[path].add(str(job['job_id']))
        results = [recon._retry_failed(p, job_ids=ids, live_job_ids=live)
                   for p, ids in selected.items() if ids]
        if results:
            # Gate this new export layout on a repaired case, not a previously
            # passing case whose WF order never exercised the failure.
            singleton = next((p for p, ids in selected.items()
                              if ids and len(_read(p)['jobs']) == 1), None)
            if singleton is not None:
                record['reference_export_preflight'] = str(singleton)
            atomic_json(root/'review_launch.json', record)
    return dict(status='retry_ready', manifests=results, production_model_accepted=False)


def production_checks(production, output, *, extend=False):
    """Append fixed-geometry native checks to one unstarted shared-cap panel."""
    from kinbot.vrc_reference_probes import prepare_configs
    root = Path(output).expanduser().resolve(strict=True)
    production = Path(production).expanduser().resolve(strict=True)
    if root.is_relative_to(production):
        raise ValueError('Keep reference checks outside production.')
    with recon._manifest_lock(root), ExitStack() as locks:
        record = _saved(root)
        if Path(record['production_policy']).parent != production:
            raise ValueError('The panel and production must share one resource policy.')
        if record.get('production_checks') and not extend:
            return dict(status='already_prepared', **record['production_checks'])
        if extend:
            for parent in sorted({Path(s['manifest']).parent for s in record['systems']}):
                locks.enter_context(recon._manifest_lock(parent))
            _idle_panel(record)
        else:
            for system in record['systems']:
                state = _read(system['manifest'])
                if any(j.get('job_id') or j.get('status', 'prepared') != 'prepared' for j in state['jobs']):
                    raise ValueError('Append production checks before starting this panel.')
        template = _read(record['systems'][0]['configuration'])
        configurations, skipped = prepare_configs(production, template)
        if not configurations and not extend:
            raise ValueError('No eligible production points: '+json.dumps(skipped))
        # Preserve the original interior/far selection even if production has
        # advanced. An extension fills missing roles, never adds a second near
        # point or changes a previously prepared model/checkpoint.
        def key(p):
            return p['channel'], '_vts_fragments_asymptote_' in Path(p['graph']).name
        existing = {key(s['production_probe']) for s in record['systems'] if s.get('production_probe')}
        configurations = [(d, p) for d, p in configurations if key(p) not in existing]
        policy = run_policy.read(production) or {}
        # Preflight every native deck before publishing the new systems.
        with tempfile.TemporaryDirectory(prefix='kinbot-production-probe-') as tmp:
            for data, provenance in configurations:
                path = Path(tmp)/(data['title']+'.json')
                atomic_json(path, data)
                loaded = recon._load(path)
                basis, anchor = loaded['bases'][0], loaded['anchors'][0]
                recon.validate_input_blocks(recon._deck(loaded, basis['basis'], anchor, 'native__point',
                    previous=loaded['initial_records']['native__point']), data['title'])
        for data, provenance in configurations:
            path = root/(data['title']+'.json')
            if path.exists() or (root/data['title']).exists():
                raise FileExistsError('Partial production staging exists; use a fresh panel directory.')
            atomic_json(path, data)
            manifest = recon.prepare(path, root/data['title'])
            state = _read(manifest)
            state.update(slurm_excluded_nodes=policy['excluded_nodes'], startup_retries=2)
            atomic_json(manifest, state)
            record['files'][path.name] = file_sha256(path)
            record['systems'].append(dict(title=data['title'], configuration=str(path),
                                         manifest=str(manifest), production_probe=provenance))
        record['max_jobs'] = min(3, record['max_jobs'])
        record['production_checks'] = dict(added_tasks=sum(bool(s.get('production_probe')) for s in record['systems']), skipped=skipped,
            total_tasks=sum(len(_read(s['manifest'])['jobs']) for s in record['systems']),
            max_jobs=record['max_jobs'], production_model_accepted=False)
        atomic_json(root/'review_launch.json', record)
    return dict(status='prepared', **record['production_checks'])


def _recheck_seed_jobs(root, identifiers):
    """Explicit one-time rechecks of unpublished, independent legacy probes.

    All targets are checked before archival; the ordinary runner subsequently
    submits them through the same admission policy. Old IDs in attempt history
    make a repeated handoff a no-op, including after a partial application.
    """
    from kinbot.recon_scf_seed import projection_failure
    identifiers = {str(i) for i in identifiers}
    if not identifiers or any(not i.isdigit() for i in identifiers):
        raise ValueError('Seed recheck requires numeric Slurm job IDs.')
    record = _saved(root)  # A changed/accepted live request must not be replayed.
    paths = [Path(s['manifest']) for s in record['systems']]
    with ExitStack() as locks:
        for parent in sorted({p.parent for p in paths}):
            locks.enter_context(recon._manifest_lock(parent))
        queued = subprocess.run(['squeue', '-h', '-o', '%A'], capture_output=True,
                                text=True, check=True, timeout=30)
        live = set(queued.stdout.split())
        snapshots = {p: recon._audit(p, live_job_ids=live) for p in paths}
        found, selected = set(), {}
        for p, m in snapshots.items():
            for j in m['jobs']:
                for h in j.get('attempt_history', []):
                    if (str(h.get('job_id')) in identifiers
                            and h.get('recovery', {}).get('kind') == 'reviewed_seed_recheck'):
                        found.add(str(h['job_id']))
                identifier = str(j.get('job_id'))
                if identifier not in identifiers:
                    continue
                if identifier in found or identifier in live:
                    raise ValueError('Seed recheck target is live or has ambiguous attempt identity.')
                if (not projection_failure(j, recheck_completed=True)
                        or any(k.get('orbital_predecessor') == j['name'] for k in m['jobs'])
                        or any(k.get('submission_status') == 'submitting' for k in m['jobs'])):
                    raise ValueError('Seed recheck target is not an intact independent legacy completion: '+identifier)
                scheduler = recon.slurm_job_record(identifier)
                if scheduler.get('JobId') == identifier and (
                        scheduler.get('WorkDir') != j['directory']
                        or scheduler.get('JobState') not in recon.SLURM_TERMINAL_STATES):
                    raise RuntimeError('Slurm still owns the recheck target, or its identity differs.')
                found.add(identifier)
                selected.setdefault(p, set()).add(identifier)
        if found != identifiers:
            raise ValueError('Unknown or inapplicable seed recheck job IDs: '+','.join(sorted(identifiers-found)))
        # All participating paths stay locked through validation and archival.
        for p, ids in selected.items():
            result = recon._retry_failed(p, job_ids=ids, live_job_ids=live, recheck_scf_seed=True)
            print('Archived completed seed probes for explicit recheck: '
                  + ', '.join(result['retry_ready']), flush=True)


def run(output, *, report_failures=False, queue_all=False, allow_nodes=(), wait_owner=0.,
        recheck_seed_jobs=()):
    root = Path(output).expanduser().resolve(strict=True)
    record = _read(root/'review_launch.json')
    if report_failures and not (record.get('reviews') or record.get('candidate_comparison')):
        raise ValueError('Failure reporting requires independent screens or a candidate comparison.')
    with recon._manifest_lock(root, wait=wait_owner):
        record = _read(root/'review_launch.json')
        if recheck_seed_jobs:
            _recheck_seed_jobs(root, recheck_seed_jobs)
        try:
            result = (_run(root, queue_all=queue_all, allow_nodes=allow_nodes)
                      if queue_all or allow_nodes else _run(root))
        except recon.ProbeFailure as error:
            if not report_failures:
                raise
            # Batch run_many(keep_going=True) raises this only after auditing
            # all allocations to termination. No still-live jobs can overlap
            # the next independent campaign using these reserved slots.
            result = dict(status='failed', error=str(error), production_model_accepted=False)
        if record.get('shared_pool'):
            # Retire the last verified claims even if production owners loaded
            # older code and no further review submission will occur.
            run_policy.reconcile_allocations(root)
        return result


def _reference_preflight(record, manifests):
    """Validate one complete native export before allowing panel fan-out.

    Reuse one existing singleton production probe, never stage extra QC. The
    other manifests still count against capacity and are not submitted here.
    Re-audit on resume rather than trusting a stale 'passed' receipt.
    """
    from kinbot.molpro_reference_identity import report, state_records
    preferred = record.get('reference_export_preflight')
    system = next((s for s in record['systems'] if str(s['manifest']) == preferred), None)
    if system is None:
        system = next((s for s in record['systems'] if s.get('production_probe')), None)
    if system is None:
        raise ValueError('Reference panel needs a singleton production probe for native preflight.')
    path = Path(system['manifest'])
    manifest = _read(path)
    if len(manifest['jobs']) != 1:
        raise ValueError('Reference preflight must contain exactly one existing probe.')
    print('Native reference-export preflight: '+system['title']+'; remaining probes wait.', flush=True)
    result = recon.run_many([path], max_jobs=1, allocation_retries=2,
                            capacity_manifests=[p for p in manifests if p != path])
    manifest = result[0]; job = manifest['jobs'][0]
    if job.get('status') != 'complete':
        raise RuntimeError('Native reference preflight has not completed validation.')
    mult = next(c['multiplicity'] for c in job['correlation_checks']
                if c['method'].lower() == 'uccsd(t)-f12b')
    receipt = _read(Path(job['directory'])/'execution.json')
    keys = {'s', 'cas', 'hf', 'dhf'} | {f'dr{r}' for r in state_records(job['state_model'], mult)}
    expected = {job['name']+'.identity_'+k+'.out' for k in keys}
    if set(receipt.get('reference_identity_matrices', {})) != expected:
        raise RuntimeError('Native reference preflight lacks the complete set of root-density exports.')
    result = report(job['directory'], job['name'], job['probe_active_space'],
                    manifest['electron_count'], mult)
    if result.get('verdict') == 'unavailable':
        raise RuntimeError('Native reference preflight has no valid matrices.')
    print('Native export/energy preflight passed; state identity still requires scientific review.', flush=True)


def _run(output, *, queue_all=False, allow_nodes=()):
    root = Path(output).expanduser().resolve(strict=True)
    record = _saved(root)
    policy = run_policy.read(record.get('resource_root', Path(record['request']).parent)) or {}
    shared = record.get('shared_pool', False)
    if shared and (queue_all or allow_nodes):
        raise ValueError('Shared-pool reviews use the live pool limit and node exclusions.')
    slots = record['max_jobs'] if shared else policy.get('orbital_review_reserved_slots', 0)
    if shared and run_policy.find(root) != Path(record['production_policy']):
        raise ValueError('Shared review has moved outside its production policy.')
    if slots < 1 and not queue_all:
        raise ValueError('No review capacity remains reserved; production may have resumed.')
    for system in record['systems']:
        config = _read(system['configuration'])
        if (config.get('probe_definition') or {}).get('orbital_source') == 'scf':
            # Refuse changed source evidence before allocation; the executor
            # repeats this check to close the submission-to-execution interval.
            recon._load(system['configuration'])
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
            if shared and (state.get('shared_policy_file') != record['production_policy']
                           or run_policy.find(manifest.parent) != Path(record['production_policy'])):
                raise ValueError('Shared review manifest policy changed.')
            if queue_all:
                allowed = set(run_policy._nodes(allow_nodes))
                previous = set(state.get('slurm_inherited_exclusion_allow', []))
                state['slurm_excluded_nodes'] = sorted(
                    (set(state.get('slurm_excluded_nodes', [])) - (allowed - previous))
                    | (set(policy['excluded_nodes']) - allowed))
                state['slurm_inherited_exclusion_allow'] = sorted(allowed)
            else:
                if allow_nodes:
                    raise ValueError('Node overrides require explicit --queue-all review scheduling.')
                if not shared:
                    state['slurm_excluded_nodes'] = sorted(set(state.get('slurm_excluded_nodes', []))
                                                           | set(policy['excluded_nodes']))
            atomic_json(manifest, state)
    previous = Path.cwd()
    try:
        os.chdir(root)
        if record.get('production_checks'):
            _reference_preflight(record, manifests)
        options = {'keep_going': True} if (record.get('reviews') or record.get('production_checks')
                                         or record.get('candidate_comparison')) else {}
        cap = sum(len(_read(p)['jobs']) for p in manifests) if queue_all else min(record['max_jobs'], slots)
        if queue_all:
            record.setdefault('queue_all_previous_max_jobs',record['max_jobs'])
            record.update(max_jobs=cap,queue_all_scheduling=dict(allowed_nodes=sorted(allowed),
                          production_model_accepted=False))
            atomic_json(root/'review_launch.json',record)
            print(f'Queueing all {cap} prepared probes; orbital predecessor dependencies remain enforced.', flush=True)
        result = recon.run_many(manifests, max_jobs=cap, allocation_retries=2,
                               **({'startup_checked': True} if queue_all else {}), **options)
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
            space = j.get('probe_active_space')
            if not space and 'energy_hartree' in j:
                lines.append(f'SCF Eh {j["energy_hartree"]:.10f}')
            occ = j.get('natural_occupations', [])
            if occ and space:
                n = (state['electron_count']-space['active_electrons'])//2
                stop = n+space['active_orbitals']
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


def character(output):
    """Hash-bound printed orbital character; no population or state assignment."""
    from kinbot.molpro_diagnostics import _orbital_tables, parse_native_output
    root = Path(output).expanduser().resolve(strict=True)
    record = _read(root/'review_launch.json')
    lines = [evidence(root),
        'AO entries are amplitudes, not populations. SCF occupations are not correlated natural occupations.',
        'SCF character alone does not establish fragment charge/spin, roots, weights or an accepted CAS.']
    for system in record['systems']:
        manifest = _read(system['manifest'])
        for job in manifest['jobs']:
            if job.get('status') != 'complete' or job.get('validation_errors'):
                continue
            output = Path(job['directory'])/(job['name']+'.out')
            if file_sha256(output) != job.get('output_sha256'):
                raise ValueError(f'Native output changed: {output}')
            lines.append('\nNATIVE '+str(output))
            for kind in ('wfu', 'molden'):
                if not job.get(kind):
                    continue
                path = Path(job[kind])
                if file_sha256(path) != job.get(kind+'_sha256'):
                    raise ValueError(f'Native {kind} changed: {path}')
                lines.append(kind.upper()+' '+str(path)+' SHA256 '+job[kind+'_sha256'])
            for kind, export in job.get('orbital_xml', {}).items():
                if file_sha256(export['path']) != export['sha256']:
                    raise ValueError('Native orbital XML changed: '+export['path'])
                lines.append('XML '+kind+' '+export['path']+' SHA256 '+export['sha256'])
            text = output.read_text(errors='replace')
            if job.get('probe_active_space'):
                audit = parse_native_output(text)
                tables = [t for s in audit['native_sections'] for t in s.get('orbital_tables', [])
                          if 'NATURAL' in t['heading'].upper()]
                label = 'NATURAL'
            else:
                normalized = re.sub(r'(?im)^\s*ELECTRON ORBITALS\s*$', 'CANONICAL ORBITALS', text)
                tables = _orbital_tables(normalized.splitlines(), 1, None, None)
                label = 'SCF'
            if not tables:
                lines.append('ORBITAL TABLE MISSING; inspect saved Molden/WFU before selection.')
                continue
            for orbital in tables[-1]['orbitals']:
                ao = sorted(orbital['ao_coefficients'], key=lambda v: abs(v['coefficient']), reverse=True)[:6]
                lines.append(f'{label} {orbital["index"]} occ={orbital.get("occupation")} '
                    f'energy_Eh={orbital.get("energy_Eh")} AO='+str([
                        (v['center'], v['function'], v['component'], v['coefficient']) for v in ao]))
    return '\n'.join(lines)


def reference_report(output):
    """Read the finite reference panel; report evidence, never approve science."""
    from kinbot.molpro_diagnostics import _orbital_tables, parse_native_output
    root = Path(output).expanduser().resolve(strict=True)
    record = _saved(root)
    policy = _read(root/'recipe.json')['reference_check']
    target, reference = policy['target_multiplicity'], policy['reference_multiplicity']
    rows, lines = {}, [evidence(root),
        'Reference audit: convergence is not state-identity or CC-quality acceptance.']

    def orbitals(label, table, low, high):
        for orbital in table['orbitals']:
            if not low <= int(orbital['index'].split('.')[0]) <= high:
                continue
            ao = sorted(orbital['ao_coefficients'],
                        key=lambda v: abs(v['coefficient']), reverse=True)[:5]
            lines.append(f'{label} {orbital["index"]} occ={orbital["occupation"]} AO='
                         + str([(v['center'], v['component'], v['coefficient']) for v in ao]))

    for system in record['systems']:
        manifest = _read(system['manifest'])
        production = system.get('production_probe')
        reference = production['reference_multiplicity'] if production else policy['reference_multiplicity']
        target = production['target_multiplicity'] if production else policy['target_multiplicity']
        for job in manifest['jobs']:
            if job['status'] != 'complete' or job.get('validation_errors'):
                continue
            path = Path(job['directory'])/(job['name']+'.out')
            if file_sha256(path) != job['output_sha256']:
                raise ValueError(f'Native output changed: {path}')
            text = path.read_text(errors='replace')
            lines.append('\nNATIVE '+system['title']+'/'+job['name'])
            if production:
                lines.append('PRODUCTION CONTRACT '+json.dumps(production))
            audit = parse_native_output(text)
            tables = [t for section in audit['native_sections']
                      for t in section.get('orbital_tables', [])
                      if 'NATURAL' in t['heading'].upper()]
            closed = (manifest['electron_count']-job['probe_active_space']['active_electrons'])//2
            if tables:
                orbitals('CAS', tables[-1], closed+1,
                         closed+job['probe_active_space']['active_orbitals'])
            else:
                lines.append('CAS ORBITAL TABLE MISSING: inspect native output.')
            for section in audit['native_sections']:
                if section['program'].upper() not in ('MULTI', 'MCSCF'):
                    continue
                for table in section.get('ci_tables', []):
                    lines.append('CI '+table['text'][:1600])
            # The last RHF section is the independent reference immediately
            # before CC, not the initial doublet used to enter MULTI.
            blocks = re.split(r'(?im)^\s*PROGRAM\s*\*\s*', text)[1:]
            scf = [b for b in blocks if re.match(r'RHF(?:-SCF)?\b', b, re.I)]
            if scf:
                normalized = re.sub(r'(?im)^\s*ELECTRON ORBITALS\s*$',
                                    'CANONICAL ORBITALS', scf[-1])
                tables = _orbital_tables(normalized.splitlines(), 1, None, None)
                if tables:
                    nclosed = (manifest['electron_count']-reference+1)//2
                    orbitals('REFERENCE ROHF', tables[-1], nclosed+1, nclosed+reference-1)
                else:
                    lines.append('REFERENCE ORBITAL TABLE MISSING: inspect reference Molden/WFU.')
            diagnostics = [s.strip() for s in text.splitlines() if re.search(
                r'T1 diagnostic|D1 diagnostic|Norm of t[12]|Expectation value of S|'
                r'\(T\).*contribution|triples.*(?:energy|contribution)', s, re.I)]
            lines.extend('CC '+s for s in diagnostics[-12:])
            from kinbot.molpro_reference_identity import report as identity_report
            try:
                identity = identity_report(job['directory'], job['name'],
                                           job['probe_active_space'], manifest['electron_count'], reference)
            except (ValueError, OSError, KeyError) as exc:
                identity = dict(verdict='unavailable', reason=str(exc))
            reference_energies = [e for e in job.get('correlated_energies_hartree', [])
                if e['multiplicity'] == reference and e['method'].lower() != 'uccsd(t)-f12b']
            if reference_energies:
                root = (min(reference_energies, key=lambda e:e['energy_hartree'])['root']
                        if production and production['reference_selection'] == 'lowest'
                        else production['reference_root'] if production else policy['root'])
                identity['selected_mr_reference_root'] = root
                identity['correlated_reference_roots'] = reference_energies
                if 'closest_density_root' in identity:
                    identity['closest_density_matches_selected_root'] = identity['closest_density_root'] == root
                    identity['flags'] = []
                    if identity['closest_density_root'] != root:
                        identity['flags'].append('Selected MR reference is not the closest ROHF density; inspect root identity/degeneracy.')
                    if min(identity['rohf_somo_active_projection_eigenvalues'], default=0.) < .9:
                        identity['flags'].append('ROHF open-shell subspace projects below 90% into the MR active span.')
                    identity['threshold_note'] = '90% is a diagnostic screening threshold, not a production acceptance criterion.'
            lines.append('NATIVE IDENTITY '+json.dumps(identity))
            if production:
                continue  # DZ spot checks cannot stand in for a new CBS curve.
            energies = job.get('correlated_energies_hartree', [])
            values = []
            for mult, cc in ((reference, True), (target, False), (reference, False)):
                found = [v['energy_hartree'] for v in energies
                         if v['multiplicity'] == mult and v['root'] == policy['root']
                         and (v['method'].lower() == 'uccsd(t)-f12b') == cc]
                if len(found) != 1 or not np.isfinite(found[0]):
                    raise ValueError(f'{job["name"]}: reference panel energies are missing or ambiguous.')
                values.append(found[0])
            key = (job['distance_angstrom'], job['basis'])
            if key in rows:
                raise ValueError('Reference report requires one geometry per distance and basis.')
            rows[key] = np.array(values)
    # A deliberately independent DZ/TZ calculation of the production formula.
    # Report the finite 30 A reference explicitly; do not claim a proven limit.
    low, high = policy['bases']
    power = policy['extrapolation_power']
    complete = {}
    for distance in sorted({r for r, _ in rows}):
        if (distance, low) in rows and (distance, high) in rows:
            complete[distance] = (3.**power*rows[distance, high]
                                  - 2.**power*rows[distance, low])/(3.**power-2.**power)
    far = policy['far_distance_angstrom']
    lines.append('\nCBS columns: R/A, CCref/Eh, MRtarget/Eh, MRref/Eh, '
                 'MRgap/Eh, anchored interaction/kcal mol-1, direct MR interaction/kcal mol-1')
    if far not in complete:
        lines.append('Waiting for both far-reference results; no interaction curve yet.')
    else:
        anchor = complete[far][0]+complete[far][1]-complete[far][2]
        for distance, (cc, low_spin, high_spin) in complete.items():
            lines.append(f'{distance:g} {cc:.10f} {low_spin:.10f} {high_spin:.10f} '
                         f'{low_spin-high_spin:.10f} '
                         f'{(cc+low_spin-high_spin-anchor)*627.509474:.6f} '
                         f'{(low_spin-complete[far][1])*627.509474:.6f}')
    lines.append('No production input changed. Inspect CC/MR reference character, '
                 'active span and diagnostics before refresh; three distances do not certify continuity.')
    return '\n'.join(lines)


def production_evidence(production):
    """Inspect actual accepted scan artifacts without replaying native QC.

    Select the latest accepted interior point and the separated reference of
    every saved leaf. Never substitute a proposed recipe or a rejected point.
    """
    from ase.io import read as read_atoms
    from kinbot.molpro_diagnostics import parse_native_output
    root = Path(production).expanduser().resolve(strict=True)
    lines = ['READ ONLY: native production state evidence, not scientific approval.']
    for request in sorted(root.glob('*/[0-9]*/vrctst/scan_review/*/scan_review.json')):
        saved = _read(request)
        lines.append('\nCHANNEL '+str(request.relative_to(root)))
        lines.append('REVIEW '+str(saved.get('status'))+' '+str(saved.get('error')))
        points = {**saved.get('history', {}), **saved.get('points', {})}
        accepted = [v for v in points.values() if v.get('status') == 'accepted' and v.get('graph')]
        if not accepted:
            lines.append('No accepted native points available.')
            continue
        far = saved['context']['separation_angstrom']
        inner = [v for v in accepted if v['distance_angstrom'] < far]
        chosen = ([max(inner, key=lambda p: p['distance_angstrom'])] if inner else [])
        chosen += [v for v in accepted if v['distance_angstrom'] == far]
        for point in chosen:
            graph = Path(point['graph'])
            if not graph.resolve().is_relative_to(root):
                raise ValueError(f'Production graph lies outside this run: {graph}')
            for name, expected in point.get('hashes', {}).items():
                if file_sha256(name) != expected:
                    raise ValueError(f'Accepted production evidence changed: {name}')
            directory = graph/'tasks/relaxed_scan'
            if (directory/'scan_rejection.json').exists():
                lines.append('REJECTED marker present; excluding '+str(graph))
                continue
            execution = _read(directory/'execution.json')
            if execution.get('status') != 'executed':
                lines.append('Incomplete execution; excluding '+str(graph))
                continue
            task = _read(directory/'task.json')['task']
            profile = task['profile']
            options = profile['calculator_kwargs']
            lines.append('POINT '+str(point['distance_angstrom'])+' A '+str(graph.relative_to(root)))
            lines.append('PROFILE '+json.dumps({k: profile.get(k) for k in ('method', 'basis')}))
            model = options['state_model']
            lines.append('MODEL '+json.dumps({k: model.get(k) for k in
                ('state_average', 'selected_roots', 'energy_selection', 'correlation_options')}))
            generated = execution.get('details', {}).get('generated_files', [])
            inputs = [directory/Path(n).name for n in generated
                      if re.search(r'_step_\d+\.inp$', str(n))]
            if not inputs:
                lines.append('No recorded final force evaluation; native identity unresolved.')
                continue
            inp = inputs[-1]
            actual = read_atoms(inp.with_suffix('.xyz'))
            final = read_atoms(directory/'final.xyz')
            if (actual.get_chemical_symbols() != final.get_chemical_symbols()
                    or not np.allclose(actual.positions, final.positions, atol=1.e-7, rtol=0.)):
                lines.append('Last force geometry differs from final.xyz; inspect manually.')
                continue
            out = inp.with_suffix('.out')
            audit = parse_native_output(out.read_text(errors='replace'), inp.read_text())
            # Numerical-gradient output repeats MULTI at displaced geometries.
            # The first MULTI is the undisplaced final-force evaluation.
            sections = [s for s in audit['native_sections'] if s['program'].upper() in ('MULTI', 'MCSCF')]
            lines.append('NATIVE '+str(out.relative_to(root))+' sha256='+file_sha256(out))
            if not sections:
                lines.append('No native MULTI section; identity unresolved.')
                continue
            section = sections[0]
            lines.append('NATIVE STATES '+json.dumps(section['states']))
            lines.append('NATIVE ROOTS '+json.dumps(section['state_energies']))
            tables = [t for t in section['orbital_tables'] if 'NATURAL' in t['heading'].upper()]
            if tables:
                for orbital in tables[-1]['orbitals']:
                    if orbital['role'] != 'active':
                        continue
                    ao = sorted(orbital['ao_coefficients'], key=lambda a: abs(a['coefficient']), reverse=True)[:4]
                    lines.append(f'ACTIVE {orbital["index"]} occ={orbital["occupation"]} AO='
                        +str([(v['center'], v['component'], v['coefficient']) for v in ao]))
            for table in section['ci_tables']:
                lines.append('CI '+table['text'][:1000])
        correction = request.parents[2]/('corr_'+saved['reaction']+'.json')
        if correction.exists():
            payload = _read(correction)
            lines.append('CORRECTION '+json.dumps({k: payload.get(k) for k in
                ('correction_model', 'levels', 'e_inf_high', 'e_inf_samp')}))
            lines.append('Native correction components: '+str(list(payload.get('native_components', {}))))
            lines.append('CC/MR reference identity still requires inspection of those native components.')
        else:
            lines.append('CC/MR correction components not yet published; their state agreement is NOT verified.')
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
    previous = _read(case_input)
    proposed, summary, chosen = _selection(source, case_input, recipe, previous)
    if not apply:
        return summary
    return _publish_selection(case_input, previous, proposed, chosen, summary)


def proposal(source, case_input, recipe, output):
    """Stage a hash-bound additive model for a single controller migration.

    No live input, review request, job or controller is changed here.
    """
    case_input = Path(case_input).expanduser().resolve(strict=True)
    source = Path(source).expanduser().resolve(strict=True)
    recipe = Path(recipe).expanduser().resolve(strict=True)
    previous = _read(case_input)
    proposed, summary, chosen = _selection(source, case_input, recipe, previous)
    saved = _saved(source, allow_accepted=True)
    system = _selected_system(saved, chosen)
    manifest_path = Path(system['manifest'])
    manifest = _read(manifest_path)
    files = [recipe, source/'request.json', Path(saved['request']), manifest_path,
             source/'review_launch.json', Path(system['configuration'])]
    required = set(chosen.get('required_complete_jobs', [])) | {
        v['job'] for v in chosen['basis_seeds'].values()}
    for job in manifest['jobs']:
        if job['name'] in required:
            files.extend([Path(job['wfu']), Path(job['directory'])/(job['name']+'.out')])
    record = dict(schema=1, kind='orbital_model_proposal', case_input=str(case_input),
                  input_sha256=file_sha256(case_input),
                  models=proposed['vrc_tst_surface_models'], selection=chosen, summary=summary,
                  evidence={str(p):file_sha256(p) for p in files})
    output = Path(output).expanduser().resolve()
    if output in files or output == case_input:
        raise ValueError('Proposal must not overwrite an input or evidence file.')
    atomic_json(output, record)
    return dict(summary, status='staged_for_controller_refresh', proposal=str(output),
                sha256=file_sha256(output))


def _selected_system(saved, chosen):
    """Resolve an explicitly named candidate, never the first candidate by default."""
    systems = saved['systems']
    title = chosen.get('source_model')
    matches = ([s for s in systems if s.get('title') == title]
               if title is not None else systems)
    if len(matches) != 1:
        raise ValueError('Select one source screen with basis-specific orbital records; '
                         'source_model must identify exactly one candidate title.')
    return matches[0]


def _selection(source, case_input, recipe, previous, *, bind_request=None):
    from kinbot import input_refresh
    from kinbot.molpro import validate_vrc_state_model
    from kinbot.vrc_multistate import validate_surface_model, extend_surface_model
    source = Path(source).expanduser().resolve(strict=True)
    case_input = Path(case_input).expanduser().resolve(strict=True)
    saved = _saved(source, allow_accepted=True)
    if saved.get('provisional') and bind_request is None:
        raise ValueError('Bind the completed screen to the actual VRC review request before selection.')
    request = _read(bind_request or source/'request.json')
    chosen = _read(recipe)
    _require_available_recipe(chosen)
    if chosen.get('diagnostic_only'):
        raise ValueError('A diagnostic-only orbital seed cannot be published to production.')
    if (chosen.get('schema') != 1 or chosen.get('reviewed') is not True
            or not chosen.get('review_note') or not chosen.get('request_match')):
        raise ValueError('An explicit reviewed selection and chemistry note are required.')
    from kinbot.electronic_states import _canonical_smiles
    def match_value(k, value):
        if k == 'products':
            return sorted((_canonical_smiles(p['smiles']), p['charge'], p['multiplicity']) for p in value)
        return value
    if any(match_value(k, request.get(k)) != match_value(k, v)
           for k,v in chosen['request_match'].items()):
        raise ValueError('Selected orbital recipe does not match this review request.')
    if case_input.parent not in Path(saved['request']).parents:
        raise ValueError('Case input does not belong to this review request.')
    system = _selected_system(saved, chosen)
    manifest = _read(system['manifest'])
    if any(j.get('probe_method') == 'casci' and
           (j.get('probe_active_space') or {}).get('orbital_partition')
           for j in manifest['jobs']):
        raise ValueError('Seeded fixed-orbital CASCI is diagnostic only; '
                         'optimize and review the production state model first.')
    jobs = {j['name']: j for j in manifest['jobs']}
    if any(not j.get('natural_orbital_probe') for j in manifest['jobs']):
        raise ValueError('SCF character is not a reviewed multistate natural-orbital selection.')
    for name in chosen.get('required_complete_jobs', []):
        job = jobs[name]
        if (job['status'] != 'complete' or job.get('validation_errors')
                or file_sha256(job['wfu']) != job['wfu_sha256']
                or file_sha256(Path(job['directory'])/(name+'.out')) != job['output_sha256']):
            raise ValueError(f'Required orbital evidence is incomplete or changed: {name}.')
    active = deepcopy(chosen['active_space'])
    basis_models = {}
    for basis, selection in chosen['basis_seeds'].items():
        job = jobs[selection['job']]
        if (job['status'] != 'complete' or job.get('validation_errors') or job['basis'] != basis
                or file_sha256(job['wfu']) != job['wfu_sha256']
                or file_sha256(Path(job['directory'])/(job['name']+'.out')) != job['output_sha256']):
            raise ValueError('Selected native seed is not unchanged, complete and basis matched.')
        for key in ('wfu_sha256', 'output_sha256'):
            if key in selection and selection[key] != job[key]:
                raise ValueError('Selected native seed differs from the explicitly reviewed hash: '
                                 + job['name'] + ' ' + key)
        model = deepcopy(chosen['state_model'])
        model['orbital_selection'] = dict(mode='molpro_record', path=job['wfu'],
            sha256=job['wfu_sha256'], record=selection['record'], orbital_type='natural',
            align_active_orbitals=False, active_indices=selection['active_indices'])
        basis_models[basis] = validate_vrc_state_model(model,
            active_electrons=active['active_electrons'], active_orbitals=active['active_orbitals'],
            nelectron=manifest['electron_count'], atom_count=len(request['symbols']))
    parameters = deepcopy(chosen.get('parameters', {}))
    parameters.update(vrc_tst_active_space=active, vrc_tst_state_model_file='', vrc_tst_state_model=dict(
        schema=2, molecule=dict(symbols=request['symbols'], charge=request['charge']),
        active_space=active, state_average_policy='reviewed', basis_models=basis_models))
    stage = dict(name='single', coordinate=chosen['coordinate'],
        distances_angstrom=deepcopy(previous[chosen['distances_parameter']]),
        states=deepcopy(chosen['states']), parameters=parameters)
    # Match JSON's string keys on repeat application (selected_roots is
    # normalized to integer keys by the state-model validator).
    channel = json.loads(json.dumps(dict(bond=request['bond'], reaction=request['reaction'],
                                         combination='single', bottlenecks=[stage])))
    parent = request['parent']
    models = deepcopy(previous.get('vrc_tst_surface_models', {}))
    prior = models.get(parent)
    model = deepcopy(prior) if prior else dict(schema=1, reviewed=True,
        review_note=chosen['review_note'], channels=[])
    matches = [c for c in model['channels'] if sorted(c['bond']) == sorted(channel['bond'])
               and c.get('reaction', request['reaction']) == request['reaction']]
    # A repeat of an already published legacy selection must remain byte-for-byte
    # additive. It cannot convert that bond model into a scoped model implicitly.
    if len(matches) == 1 and 'reaction' not in matches[0]:
        channel.pop('reaction')
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
    config = _read(system['configuration'])
    a = np.asarray(config['geometry_angstrom'], dtype=float)
    b = np.asarray(request['positions'], dtype=float)
    distances = lambda xyz: np.linalg.norm(xyz[:, None]-xyz[None, :], axis=2)
    from kinbot.vrc_correction import asymptote_distance
    effective = {**previous, **parameters}
    summary = dict(status='reviewed', parent=parent, active_space=active,
        basis_seeds={b:m['orbital_selection'] for b,m in basis_models.items()},
        source_screen=str(source), source_model=system.get('title'),
        request=str(bind_request or saved['request']),
        geometry_provenance=dict(coordinates_equal=bool(np.array_equal(a, b)),
            max_pair_distance_change_angstrom=float(np.max(abs(distances(a)-distances(b)))),
            note='Source is a fixed-geometry orbital screen; production reoptimizes the current geometry.'),
        scan_plan=dict(anchor_atoms_one_based=request['bond'],
            correction_distances_angstrom=effective.get('vrc_tst_scan_points'),
            asymptote_distance_angstrom=asymptote_distance(effective),
            fragment_geometry=effective.get('vrc_tst_fragment_geometry'),
            unspecified_parameters='Unspecified scan/fragment settings retain the owning driver defaults.',
            sampling_coordinate=stage['coordinate'],
            sampling_distances_angstrom=stage['distances_angstrom'],
            geometry_source=(_read(source/'geometry_source.json')
                if (source/'geometry_source.json').is_file() else None)),
        drivers_started=0, reconnaissance_jobs_submitted=0)
    return proposed, summary, chosen


def _publish_selection(case_input, previous, proposed, chosen, summary):
    from kinbot import input_refresh
    # Validate with the real publisher before making even an archived proposal.
    # Its apply call repeats validation under locks before writing the input.
    input_refresh.submit(case_input, proposal=proposed, expected=previous)
    if previous != proposed:
        import hashlib
        digest = hashlib.sha256(json.dumps(proposed,sort_keys=True).encode()).hexdigest()
        archive = case_input.parent/'orbital_reviews/live_refresh'/digest
        archive.mkdir(parents=True, exist_ok=True)
        # Preserve the exact pre-edit input once. A retry may find the proposal
        # already installed; submit() reconciles its existing receipt then.
        for name, payload in [('input_before.json',previous), ('input_after.json',proposed),
                              ('selection.json',chosen), ('review_summary.json',summary)]:
            path = archive/name
            if path.exists() and _read(path) != payload:
                raise ValueError('Existing review archive differs; nothing was overwritten.')
            atomic_json(path,payload)
    # Uses the existing journaled live-refresh path, which validates the owner,
    # accepted input, worker state and additive models. It never starts a driver.
    refresh = input_refresh.submit(case_input,apply=True,proposal=proposed,expected=previous)
    return dict(summary, refresh=refresh, status=refresh['status'])


def _review_source(campaign, relative):
    relative = Path(relative)
    source = (campaign/relative).resolve(strict=True)
    if relative.is_absolute() or not source.is_relative_to(campaign) or source == campaign:
        raise ValueError('Selection sources must be relative paths beneath campaign.')
    return source


def _supporting_review(source):
    """Hash-check a completed comparison used in a human selection rationale."""
    saved = _saved(source, allow_accepted=True)
    evidence = {}
    for system in saved['systems']:
        manifest = _read(system['manifest'])
        if not manifest['jobs'] or any(j['status'] != 'complete' for j in manifest['jobs']):
            return None
        for job in manifest['jobs']:
            files = [(Path(job['wfu']), job['wfu_sha256']),
                     (Path(job['directory'])/(job['name']+'.out'), job['output_sha256'])]
            if job.get('validation_errors'):
                raise ValueError('Supporting comparison failed native validation: '+job['name'])
            for path, expected in files:
                if file_sha256(path) != expected:
                    raise ValueError('Supporting comparison evidence changed: '+str(path))
                evidence[str(path)] = expected
    return dict(source=str(source), files=evidence)


def review_batch(production, campaign, plan, *, apply=False, require_all=False):
    """Install ready, explicitly selected channels once per case, without QC.

    Provisional screens can bind only after actual accepted review requests
    appear. Waiting cases do not block independent ready cases. Re-running
    reconciles receipts and never starts a driver or a reconnaissance job.
    ``require_all`` requires one case and publishes its entire selection list
    in a single input revision, or none of it. Sources may be relative paths
    beneath campaign, so completed screens from separate campaigns can join.
    """
    from contextlib import ExitStack, contextmanager
    from kinbot import input_refresh
    production = Path(production).expanduser().resolve(strict=True)
    campaign = Path(campaign).expanduser().resolve(strict=True)
    plan = Path(plan).expanduser().resolve(strict=True)
    controls = _read(plan)
    _require_available_recipe(controls)
    if controls.get('schema') != 1 or not controls.get('selections'):
        raise ValueError('Expected an explicit selection plan.')
    cases = {}
    names = set()
    for item in controls['selections']:
        name, case = item['name'], item['case']
        if not recon._SAFE.fullmatch(name) or not recon._SAFE.fullmatch(case) or name in names:
            raise ValueError('Selection names/cases must be safe and distinct.')
        names.add(name)
        cases.setdefault(case, []).append(item)
    if require_all and len(cases) != 1:
        raise ValueError('An all-or-none publication must target one case input.')
    sources = set()
    for item in controls['selections']:
        source = _review_source(campaign, item.get('source', item['name']))
        if source in sources:
            raise ValueError('Selection sources must be distinct.')
        sources.add(source)
        for relative in item.get('supporting_reviews', []):
            _review_source(campaign, relative)
    results = {}
    for case, items in cases.items():
        case_input = production/case/'input.json'
        with ExitStack() as stack:
            if apply:
                handle = stack.enter_context((case_input.parent/'.kinbot_orbital_review.lock').open('a+'))
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            previous = _read(case_input)
            proposed = deepcopy(previous)
            selections, summaries, bindings, waiting = [], [], [], []
            for item in items:
                source = _review_source(campaign, item.get('source', item['name']))
                saved = _saved(source, allow_accepted=True)
                request = _read(source/'request.json')
                actual = case_input.parent/request['parent']/'vrctst/orbital_review'/request['reaction']/'request.json'
                if saved.get('provisional') and not actual.is_file():
                    waiting.append(dict(channel=item['name'], reason='waiting_for_accepted_channel', request=str(actual)))
                    continue
                chosen_recipe = _read(plan.parent/item['recipe'])
                _require_available_recipe(chosen_recipe)
                manifest = _read(_selected_system(saved, chosen_recipe)['manifest'])
                if any(j['status'] != 'complete' for j in manifest['jobs']):
                    waiting.append(dict(channel=item['name'], reason='waiting_for_native_completion'))
                    continue
                supporting = [_supporting_review(_review_source(campaign, relative))
                              for relative in item.get('supporting_reviews', [])]
                if any(e is None for e in supporting):
                    waiting.append(dict(channel=item['name'], reason='waiting_for_supporting_comparison'))
                    continue
                binding = actual if saved.get('provisional') else None
                if binding:
                    bind(source, binding, apply=False)
                    bindings.append((source, binding))
                proposed, summary, chosen = _selection(source, case_input,
                    plan.parent/item['recipe'], proposed, bind_request=binding)
                if supporting:
                    summary['supporting_evidence'] = supporting
                    chosen = dict(chosen, supporting_evidence=supporting)
                selections.append(chosen); summaries.append(summary)
            if not summaries or (require_all and waiting):
                results[case] = dict(status='waiting', waiting=waiting)
                continue
            summary = dict(status='reviewed', channels=summaries, waiting=waiting,
                           drivers_started=0, reconnaissance_jobs_submitted=0)
            # A preceding revision must remain byte-for-byte intact until
            # consumed. Do not edit its input and then discover it was pending.
            try:
                input_refresh.submit(case_input, proposal=proposed, expected=previous)
            except ValueError as error:
                if 'input revision is still pending' not in str(error):
                    raise
                results[case] = dict(summary, status='waiting_for_previous_refresh', error=str(error))
                continue
            if apply:
                for source, actual in bindings:
                    bind(source, actual)
                summary = _publish_selection(case_input, previous, proposed,
                                              dict(selections=selections), summary)
            results[case] = summary
    return dict(cases=results, drivers_started=0, reconnaissance_jobs_submitted=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='action', required=True)
    prepare_cmd = commands.add_parser('prepare')
    for name in ('request', 'case_input', 'recipe', 'output'):
        prepare_cmd.add_argument(name)
    prepare_cmd.add_argument('--shared-pool', action='store_true')
    candidate_cmd = commands.add_parser('candidates')
    for name in ('source', 'recipe', 'output'):
        candidate_cmd.add_argument(name)
    candidate_cmd.add_argument('--preflight', action='store_true')
    candidate_cmd.add_argument('--resume', action='store_true')
    production_cmd = commands.add_parser('production-checks', help='Append native production CC/MR identity probes.')
    for name in ('production', 'output'):
        production_cmd.add_argument(name)
    production_cmd.add_argument('--extend', action='store_true', help='Add missing roles only after every panel allocation has stopped.')
    commands.add_parser('retry-reference-exports', help='Archive known legacy density-selector failures for retry.').add_argument('output')
    batch_cmd = commands.add_parser('batch',help='Prepare one shared-cap orbital-only campaign.')
    for name in ('production','plan','output'):
        batch_cmd.add_argument(name)
    batch_cmd.add_argument('--shared-pool', action='store_true',
        help='Stage within production and admit all review jobs through its existing resource policy.')
    batch_cmd.add_argument('--preflight', action='store_true', help='Validate in temporary storage; create no campaign.')
    batch_cmd.add_argument('--resume', action='store_true', help='Verify an identical existing campaign without restaging.')
    for command in (prepare_cmd, candidate_cmd, batch_cmd):
        command.add_argument('--resource-root', type=Path,
            help='Use a separate existing scheduling pool; source requests/checkpoints remain read-only.')
    bind_cmd = commands.add_parser('bind',help='Bind preliminary evidence to an actual review request.')
    for name in ('source','request'):
        bind_cmd.add_argument(name)
    review_cmd = commands.add_parser('review', help='Apply an explicitly reviewed orbital selection without QC submissions.')
    for name in ('source', 'case_input', 'recipe'):
        review_cmd.add_argument(name)
    review_cmd.add_argument('--apply', action='store_true')
    proposal_cmd = commands.add_parser('proposal', help='Stage a reviewed model for one controller refresh; no live publication.')
    for name in ('source', 'case_input', 'recipe', 'output'):
        proposal_cmd.add_argument(name)
    reviews_cmd = commands.add_parser('review-batch', help='Publish ready selections once per case, without new QC.')
    for name in ('production', 'campaign', 'plan'):
        reviews_cmd.add_argument(name)
    reviews_cmd.add_argument('--apply', action='store_true')
    reviews_cmd.add_argument('--require-all', action='store_true',
        help='Require every selection for one case; publish one input revision or none.')
    for name in ('run', 'monitor', 'evidence', 'character', 'reference-report', 'production-evidence'):
        command = commands.add_parser(name)
        command.add_argument('output')
        if name == 'run':
            command.add_argument('--recheck-seed-job', action='append', default=[], metavar='JOB_ID',
                help='Archive/recheck an explicitly reviewed completed legacy forced-OVL SCF-seeded probe once.')
            command.add_argument('--wait-owner', type=float, default=0., metavar='SECONDS',
                help='Wait for the current review owner to finish before resuming; never stop jobs.')
            command.add_argument('--report-failures', action='store_true',
                help='Report terminal independent-batch probe failures without aborting a later campaign.')
            command.add_argument('--queue-all', action='store_true',
                help='Queue the finite prepared inventory, overriding reserved-slot and pilot throttles.')
            command.add_argument('--allow-node', action='append', default=[])
    args = parser.parse_args()
    if args.action == 'prepare':
        result = prepare(args.request, args.case_input, args.recipe, args.output,
                         shared_pool=args.shared_pool, resource_root=args.resource_root)
    elif args.action == 'batch':
        result = batch(args.production,args.plan,args.output,shared_pool=args.shared_pool,
                       preflight=args.preflight, resume=args.resume, resource_root=args.resource_root)
    elif args.action == 'bind':
        result = bind(args.source,args.request)
    elif args.action == 'candidates':
        result = candidates(args.source, args.recipe, args.output,
                            preflight=args.preflight, resume=args.resume, resource_root=args.resource_root)
    elif args.action == 'production-checks':
        result = production_checks(args.production, args.output, extend=args.extend)
    elif args.action == 'retry-reference-exports':
        result = retry_reference_exports(args.output)
    elif args.action == 'review':
        result = review(args.source, args.case_input, args.recipe, apply=args.apply)
    elif args.action == 'proposal':
        result = proposal(args.source, args.case_input, args.recipe, args.output)
    elif args.action == 'review-batch':
        result = review_batch(args.production, args.campaign, args.plan, apply=args.apply,
                              require_all=args.require_all)
    elif args.action == 'run':
        result = run(args.output, report_failures=args.report_failures,
                     queue_all=args.queue_all, allow_nodes=args.allow_node, wait_owner=args.wait_owner,
                     recheck_seed_jobs=args.recheck_seed_job)
    elif args.action == 'evidence':
        print(evidence(args.output))
        return
    elif args.action == 'character':
        print(character(args.output))
        return
    elif args.action == 'reference-report':
        print(reference_report(args.output))
        return
    elif args.action == 'production-evidence':
        print(production_evidence(args.output))
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
    if args.action == 'review-batch' and args.require_all and any(
            row['status'].startswith('waiting') for row in result['cases'].values()):
        raise SystemExit(2)


if __name__ == '__main__':
    main()
