"""Reviewed, owner-applied changes to automatic ANL calculation protocols.

The request is bound to the existing scientific graph and attempt generations.
Completed artifacts are never edited. Independent allocations can finish while
the owner archives and replaces only the affected dependency chains.
"""
from copy import deepcopy
import json
from pathlib import Path
import uuid
import time

from kinbot.execution import atomic_json, file_sha256

QUEUE = '.kinbot_controls/anl_protocol_requests'
TIERS = ('ANL1-F12', 'ANL1', 'ANL0-F12', 'ANL0', 'L3', 'L3-SP@L2', 'L2')
FAMILIES = {'core_valence': ('cv_',), 'scalar_relativistic': ('rel_',),
            'higher_order': ('ccsdt_', 'ccsdtq_', 'ccsdtqp_'), 'dboc': ('cfour_dboc',),
            'higher_order_tz': ('ccsdt_tz', 'ccsdtq_tz'),
            'higher_order_dz': ('ccsdtq_dz', 'ccsdtqp_dz')}


def validate(options):
    if not isinstance(options, dict) or not options or set(options) - {
            'tier', 'geometry', 'harmonic', 'omit_tasks', 'omit_components', 'anharmonic_qz'}:
        raise ValueError('Supply tier, geometry, harmonic, omit_tasks, omit_components or anharmonic_qz.')
    if options.get('tier', 'L2') not in TIERS:
        raise ValueError('Unknown ANL tier.')
    if options.get('geometry', 'default') not in ('default', 'tz', 'qz', 'l2'):
        raise ValueError('Geometry must be default, tz, qz or l2.')
    if options.get('harmonic', 'default') not in ('default', 'tz', 'qz', 'cbs', 'l2'):
        raise ValueError('Harmonic treatment must be default, tz, qz, cbs or l2.')
    if 'anharmonic_qz' in options and type(options['anharmonic_qz']) is not bool:
        raise ValueError('anharmonic_qz must be boolean.')
    for key in ('omit_tasks', 'omit_components'):
        if key in options and (not isinstance(options[key], list) or
                               any(not isinstance(x, str) for x in options[key])):
            raise ValueError(f'{key} must be a list of names.')
    if set(options.get('omit_components', [])) - FAMILIES.keys():
        raise ValueError('Unknown optional component family.')
    return deepcopy(options)


def sources(spec):
    p = spec.get('intent', {}).get('live_protocol', {})
    high = spec.get('intent', {}).get('highest_tier', '').startswith('ANL1')
    return (p.get('geometry_task', 'l3_geometry'),
            p.get('harmonic_tasks', ['harmonic_tz', 'harmonic_qz'] if high else ['harmonic']))


def omitted(spec):
    return set(spec.get('intent', {}).get('live_protocol', {}).get('omitted_tasks', []))


def equation(recipe, spec, *, geometry_fallback=False, harmonic_fallback=False):
    """Label the actual mixed protocol and require its native method/bases."""
    from dataclasses import replace
    geom, harmonic = sources(spec)
    tasks = {t['id']:t for t in spec['tasks']}
    requirements, labels = [], []
    for r in recipe.requirements:
        if r.key == 'reference_cbs' and not geometry_fallback:
            p = tasks[geom]['profile']
            if (r.settings.get('geometry_method'), r.settings.get('geometry_basis')) != (p['method'], p['basis']):
                settings = {**r.settings, 'geometry_method':p['method'], 'geometry_basis':p['basis']}
                basis = r.basis.split('//')[0]+'//'+p['basis'] if '//' in r.basis else r.basis
                r = replace(r, settings=settings, basis=basis)
                labels.append('geometry='+p['method']+'/'+p['basis'])
        if r.key == 'harmonic_zpe' and not harmonic_fallback:
            if len(harmonic) == 1:
                p = tasks[harmonic[0]]['result_parser']
                method = p.get('method', 'CCSD(T)')
                if (r.method, r.basis) != (method, p['basis']):
                    labels.append('harmonic='+method+'/'+p['basis'])
                settings = {k:v for k,v in r.settings.items() if k not in
                            ('extrapolation_power', 'upper_cardinal', 'geometry_mode', 'native_basis_pair')}
                r = replace(r, method=method, basis=p['basis'], backends=('molpro', 'known_zero'), settings=settings)
            r = replace(r, geometry_role='harmonic')
        requirements.append(r)
    return replace(recipe, name=recipe.name + (' ['+';'.join(labels)+']' if labels else ''),
                   requirements=tuple(requirements))


def _rename(value, mapping):
    if isinstance(value, dict):
        return {k: _rename(v, mapping) for k, v in value.items()}
    if isinstance(value, list):
        return [_rename(v, mapping) for v in value]
    if isinstance(value, str):
        return mapping.get(value, value)
    return value


def _meaning(task):
    from kinbot.manual_fallback import logical_task
    task = deepcopy(logical_task(task))
    for k in ('resources', 'runtime_setup', 'setup'):
        task.pop(k, None)
    return task


def replacement(spec, options):
    """Construct a full replacement; unused tasks remain as archived evidence."""
    from kinbot.anl.validation import automatic_composite_spec
    from kinbot.anl import dispatch, anharmonic
    options = validate(options)
    if spec.get('name') != 'anl-automatic-post-l2':
        raise ValueError('Protocol selection requires an automatic ANL graph.')
    intent = spec['intent']; old = {t['id']: t for t in spec['tasks']}
    chosen = {**intent.get('live_protocol', {}).get('options', {}), **options}
    tier = chosen.get('tier', intent['highest_tier'])
    if tier != intent['highest_tier']:
        conv = intent['convergence']
        vpt = old.get('gaussian_vpt2', {}).get('result_parser', {})
        resource = next((t['resources'] for t in spec['tasks']), {})
        recovery = next((t['optimizer'].get('minimum_recovery', {}) for t in old.values()
                         if t.get('kind') == 'ase_optimize'
                         and not t.get('anharmonic_basis_extension')), {})
        generated = automatic_composite_spec(spec['molecule'], highest_tier=tier,
            max_nodes=spec['limits']['max_nodes'], partition=resource.get('partition'),
            stationary_point=intent['stationary_point'], vpt2_method=vpt.get('method', 'B2PLYP'),
            vpt2_basis=vpt.get('basis', 'cc-pVTZ'), vpt2_dispersion=vpt.get('dispersion', 'GD3BJ'),
            scf_convergence=conv['scf'], cc_convergence=conv['coupled_cluster_or_ci'],
            geometry_fmax=conv['geometry_fmax_ev_angstrom'],
            minimum_recovery=recovery,
            f12_augmented=intent.get('f12_augmented', False), augmented=intent.get('anl_augmented', False),
            component_version=intent.get('component_version', 1),
            cv_f12_auxiliary=intent.get('cv_f12_auxiliary'),
            require_l3_geometry=intent.get('require_l3_geometry', False),
            l3_max_retries=intent.get('l3_max_retries', 0),
            anharmonic_qz=chosen.get('anharmonic_qz', True))
        # Site launch commands are not a protocol choice. Preserve same-task
        # launch configuration; otherwise require a unique backend precedent.
        def backend(t):
            return (t.get('profile', {}).get('calculator') if t.get('kind') == 'ase_optimize'
                    else t.get('result_parser', {}).get('kind', '').split('_')[0])
        for t in generated['tasks']:
            previous = old.get(t['id'])
            candidates = ([previous] if previous and previous['kind'] == t['kind']
                          and backend(previous) == backend(t) else
                          [v for v in old.values() if v['kind'] == t['kind']
                           and backend(v) == backend(t)])
            launch = []
            for v in candidates:
                payload = {k: deepcopy(v[k]) for k in ('setup', 'runtime_setup') if k in v}
                payload['command'] = (v['profile']['command'] if v['kind'] == 'ase_optimize'
                                      else v.get('command'))
                if payload not in launch:
                    launch.append(payload)
            if len(launch) > 1:
                raise ValueError(f'No unique site launch configuration for new task {t["id"]}.')
            if launch:
                payload = launch[0]
                if t['kind'] == 'ase_optimize':
                    t['profile']['command'] = payload['command']
                elif payload['command'] is not None:
                    t['command'] = payload['command']
                t.update({k:v for k,v in payload.items() if k != 'command'})
        # The TZ optimization in ANL1 is the same operation as the ANL0
        # primary geometry. Keep its original native task ID and provenance.
        mapping = {}
        for t in generated['tasks']:
            if t.get('kind') != 'ase_optimize':
                continue
            for n, previous in old.items():
                if (previous.get('kind') == 'ase_optimize' and
                        previous.get('geometry_from') == t.get('geometry_from') and
                        previous.get('profile') == t.get('profile') and
                        previous.get('optimizer') == t.get('optimizer')):
                    mapping[t['id']] = n
                    break
        # A lower tier's primary TZ ID can be reused by the new TZ task;
        # give the new QZ primary a distinct ID instead of colliding with it.
        reserved = set(old) | {t['id'] for t in generated['tasks']} | set(mapping.values())
        for t in generated['tasks']:
            ident = t['id']
            if ident not in mapping and ident in mapping.values():
                candidate = ident+'_protocol'
                while candidate in reserved:
                    candidate += '_new'
                mapping[ident] = candidate; reserved.add(candidate)
        generated['tasks'] = [_rename(t, mapping) for t in generated['tasks']]
        primary = mapping.get('l3_geometry', 'l3_geometry')
        harmonic_mapping = {}
        for t in generated['tasks']:
            parser = t.get('result_parser', {})
            if parser.get('kind') != 'molpro_harmonic':
                continue
            for n, previous in old.items():
                pp = previous.get('result_parser', {})
                if ({k:v for k,v in pp.items() if k != 'file'} ==
                        {k:v for k,v in parser.items() if k != 'file'}
                        and previous.get('geometry_from') == t.get('geometry_from')):
                    harmonic_mapping[t['id']] = n
                    t.clear(); t.update(deepcopy(previous))
                    break
        for t in generated['tasks']:
            previous = old.get(t['id'])
            if previous and _meaning(t) == _meaning(previous):
                t.clear(); t.update(deepcopy(previous))
        result = deepcopy(generated)
        result['limits'] = deepcopy(spec['limits'])
        selected = {t['id'] for t in result['tasks']}
        result['tasks'].extend(deepcopy(t) for n, t in old.items() if n not in selected)
    else:
        result = deepcopy(spec)
        policy = intent.get('live_protocol', {})
        primary = policy.get('default_geometry_task', 'l3_geometry')
        selected = set(policy.get('eligible_tasks', old))
        harmonic_mapping = {}
    if chosen.get('anharmonic_qz', True):
        before_extension = {t['id'] for t in result['tasks']}
        result = anharmonic.extend(result)
        selected.update(t['id'] for t in result['tasks'] if t['id'] not in before_extension)
    tasks = {t['id']: t for t in result['tasks']}
    high = tier.startswith('ANL1')
    harmonic = ([harmonic_mapping.get(n, n) for n in
                 (['harmonic_tz', 'harmonic_qz'] if high else ['harmonic'])]
                if tier != intent['highest_tier'] else
                intent.get('live_protocol', {}).get('default_harmonic_tasks',
                    ['harmonic_tz', 'harmonic_qz'] if high else ['harmonic']))
    if tier in ('L2', 'L3-SP@L2'):
        primary, harmonic = 'initial', []
    defaults = dict(default_geometry_task=primary, default_harmonic_tasks=harmonic[:], eligible_tasks=sorted(selected))
    gchoice = chosen.get('geometry', 'default')
    if gchoice in ('tz', 'qz'):
        candidates = [n for n, t in tasks.items() if t.get('kind') == 'ase_optimize'
            and not t.get('anharmonic_basis_extension')
            and t.get('profile', {}).get('calculator') == 'molpro'
            and ('pV' + gchoice[0].upper() + 'Z').casefold() in t['profile']['basis'].casefold()]
        if not candidates:
            raise ValueError(f'This graph has no {gchoice.upper()} CC geometry; choose an ANL tier providing it.')
        primary = candidates[0]
    elif gchoice == 'l2':
        primary = 'initial'
    hchoice = chosen.get('harmonic', 'default')
    if hchoice == 'default' and gchoice in ('tz', 'l2'):
        hchoice = gchoice
    if hchoice in ('tz', 'qz'):
        harmonic = [n for n, t in tasks.items()
                    if t.get('result_parser', {}).get('kind') == 'molpro_harmonic'
                    and ('pV' + hchoice[0].upper() + 'Z').casefold() in t['result_parser']['basis'].casefold()]
        if len(harmonic) != 1:
            raise ValueError(f'This graph has no unique {hchoice.upper()} harmonic calculation.')
    elif hchoice == 'cbs':
        harmonic = []
        for cardinal in ('T', 'Q'):
            matches = [n for n,t in tasks.items()
                       if t.get('result_parser', {}).get('kind') == 'molpro_harmonic'
                       and ('pV'+cardinal+'Z').casefold() in t['result_parser']['basis'].casefold()]
            if len(matches) != 1:
                raise ValueError('CBS harmonic treatment needs a unique declared TZ/QZ pair.')
            harmonic.extend(matches)
    elif hchoice == 'l2' or (hchoice == 'default' and primary == 'initial'):
        harmonic = []
    # Rewire electronic corrections only. A basis-specific Hessian always
    # retains its own basis-specific optimized geometry.
    previous_primary = sources(spec)[0] if tier == intent['highest_tier'] else defaults['default_geometry_task']
    for n, t in tasks.items():
        if (n in selected and t.get('kind') != 'ase_optimize'
                and t.get('result_parser', {}).get('kind') not in ('molpro_harmonic', 'gaussian_vpt2')
                and t.get('geometry_from') in (previous_primary, 'l3_geometry')):
            t['geometry_from'] = primary
            t['depends_on'] = [primary if d in (previous_primary, 'l3_geometry') else d
                               for d in t.get('depends_on', []) if primary != 'initial' or d not in (previous_primary, 'l3_geometry')]
    for n, t in tasks.items():
        if ((t.get('kind') == 'ase_optimize' and not t.get('anharmonic_basis_extension'))
                or t.get('result_parser', {}).get('kind') == 'molpro_harmonic'):
            selected.discard(n)
    selected.update(harmonic)
    if primary != 'initial':
        selected.add(primary)
    # Geometry predecessors of the selected Hessians/optimization are needed.
    def dependencies(n):
        t = tasks[n]; values = set(t.get('depends_on', [])) | {t.get('geometry_from', 'initial')}
        for dep in values - {'initial'}:
            if dep in tasks and dep not in selected:
                selected.add(dep); dependencies(dep)
    for n in list(selected):
        dependencies(n)
    explicit = set(chosen.get('omit_tasks', []))
    if explicit - tasks.keys():
        raise ValueError('Unknown omitted tasks: ' + ', '.join(sorted(explicit - tasks.keys())))
    if any(tasks[n].get('kind') == 'ase_optimize' or n in harmonic for n in explicit):
        raise ValueError('Select geometry/harmonic sources with their dedicated options.')
    for family in chosen.get('omit_components', []):
        explicit.update(n for n in tasks if n.startswith(FAMILIES[family]))
    if not chosen.get('anharmonic_qz', True):
        explicit.update({anharmonic.GEOMETRY, anharmonic.HARMONIC, anharmonic.VPT2} & tasks.keys())
    # Removing a prerequisite also removes its consumers. A retained task
    # must never wait forever for a deliberately disabled native calculation.
    from kinbot.anl.salvage import descendants
    for n in list(explicit):
        explicit.update(descendants(result, n))
    selected -= explicit
    result['intent']['live_protocol'] = dict(schema=1, options=chosen, **defaults, geometry_task=primary,
        harmonic_tasks=harmonic, omitted_tasks=sorted(tasks.keys() - selected))
    result['intent']['require_l3_geometry'] = primary != 'initial'
    dispatch.assign_partitions(result)
    dispatch.validate_spec(result)
    return result


def plan(graph, options, *, allow_closed=False):
    from kinbot.anl import dispatch, refinement
    from kinbot.anl.salvage import descendants
    graph, spec, state = dispatch._load(graph)
    reopen = None
    if not allow_closed:
        refinement.require_open(graph)
    else:
        try:
            refinement.require_open(graph)
        except ValueError:
            from kinbot import live_control as control
            worker = next((p for p in graph.parents if (p/(p.name+'.json')).is_file()), None)
            if worker is None:
                raise ValueError('Closed standalone graph needs its downstream owner reopened before editing.')
            reopen = dict(worker=str(worker), expected=control.generation(worker, 'network'))
    revised = replacement(spec, options)
    before = {t['id']: t for t in spec['tasks']}
    after = {t['id']: t for t in revised['tasks']}
    reset = {n for n in before if before[n] != after.get(n)}
    for n in list(reset):
        reset.update(descendants(spec, n))
    # Newly excluded running/failed/staged tasks are archived. Completed
    # evidence remains intact, but the assembler excludes it explicitly.
    for n in omitted(revised) - omitted(spec):
        if state['tasks'].get(n, {}).get('status') != 'complete':
            reset.add(n)
    reset.update(n for n in omitted(spec) - omitted(revised)
                 if state['tasks'].get(n, {}).get('protocol_omitted'))
    if any(state['tasks'].get(n, {}).get('status') == 'submitting' for n in reset):
        raise RuntimeError('Reconcile ambiguous submissions before changing this protocol.')
    return dict(graph=str(graph), options=validate(options), spec_sha256=state['spec_sha256'], reopen=reopen,
        reset=sorted(reset), expected={n: refinement._generation(state['tasks'].get(n, {})) for n in reset},
        affected=[dict(task=n, **{k: state['tasks'].get(n, {}).get(k) for k in ('status', 'job_id')}) for n in sorted(reset)],
        replacement=revised)


def requests(graph):
    from kinbot import live_control as control
    root = control.find(graph)
    if root is None:
        return []
    result = []
    for p in sorted((root.parent/QUEUE).glob('*.json')):
        r = json.loads(p.read_text())
        if r['graph'] == str(Path(graph).resolve()) and r['status'] in ('queued', 'applying'):
            result.append((p, r))
    return result


def pending(graph):
    return {n for _, r in requests(graph) for n in r['plan']['reset']}


def publish(root, plans, *, reason, cancel_running=False):
    from kinbot import live_control as control
    root = Path(root).resolve(strict=True)
    if control.find(root) != root/control.NAME or not reason.strip():
        raise ValueError('Enable live controls and supply a reason.')
    result = []
    with control.locked(root):
        # Validate the whole selection before publishing any request.
        for p in plans:
            if not Path(p['graph']).is_relative_to(root) or requests(p['graph']):
                raise ValueError('Graph is outside this root or already has a pending protocol request.')
            if plan(p['graph'], p['options'], allow_closed=True) != p:
                raise ValueError('Protocol preview changed; repeat it before publication.')
        for p in plans:
            ident = uuid.uuid4().hex
            r = dict(schema=1, id=ident, graph=p['graph'], plan=p, reason=reason,
                     cancel_running=bool(cancel_running), status='queued', created_at=control.now())
            atomic_json(root/QUEUE/(ident+'.json'), r); result.append(r)
    # Publish the graph fence first, then ask the existing PES owner to
    # invalidate derived assembly and reopen it at its normal epoch boundary.
    # A failed wakeup leaves a resumable fenced request, never stale assembly.
    wake_owners(root)
    return result


def wake_owners(root):
    from kinbot import live_control as control
    root = Path(root)
    for path in sorted((root/QUEUE).glob('*.json')):
        r = json.loads(path.read_text())
        reopening = r['plan'].get('reopen')
        if r['status'] not in ('queued', 'applying') or not reopening or r.get('owner_request'):
            continue
        worker = Path(reopening['worker'])
        existing = [v for _,v in control.requests(worker, kinds={'network'})
                    if v['target'] == str(worker)]
        try:
            receipt = existing[0] if existing else control.request(root, worker, kind='network',
                options={}, reason=r['reason'], expected=reopening['expected'])
            r['owner_request'] = receipt['id']; r.pop('owner_wakeup_error', None)
        except (ValueError, RuntimeError) as error:
            r['owner_wakeup_error'] = str(error)
        atomic_json(path, r)


def process(graph):
    """Only called by the owner while holding the graph drive lock."""
    from kinbot.anl import sharing
    with sharing.locked(graph):
        _process(graph)


def _process(graph):
    from kinbot import live_control as control
    from kinbot.anl import dispatch, live_control, refinement
    graph = Path(graph).resolve()
    for path, r in requests(graph):
        try:
            _, spec, state = dispatch._load(graph)
            p = r['plan']; reset = set(p['reset'])
            if state['spec_sha256'] != p['spec_sha256'] or any(
                    refinement._generation(state['tasks'].get(n, {})) != e for n, e in p['expected'].items()):
                raise ValueError('Protocol attempt generation changed; preview a new request.')
            try:
                refinement.require_open(graph)
            except ValueError:
                continue  # Existing PES owner must invalidate downstream assembly first.
            current = plan(graph, p['options'])
            if current['replacement'] != p['replacement']:
                raise ValueError('Protocol specification changed after preview.')
            by_id = {t['id']: t for t in spec['tasks']}
            live_control.check_commit(graph, current['replacement'], state, reset)
            from kinbot.anl import sharing
            if sharing.consumers(graph):
                raise RuntimeError('Shared canonical graph has receivers; detach/revise its receivers first.')
            for n in reset:
                e = state['tasks'].get(n, {})
                if e.get('task_sha256'):
                    dispatch._verify_stage_files(graph, by_id[n], e)
                    execution = graph/'tasks'/n/'execution.json'
                    if execution.exists():
                        dispatch._verify_execution(graph, by_id[n], e, json.loads(execution.read_text()))
            ready = True
            for n in reset:
                job = state['tasks'].get(n, {}).get('job_id')
                if not job or not dispatch._job_active(job):
                    continue
                if not r['cancel_running']:
                    ready = False; continue
                child_path = path.parent.parent/'protocol_cancellations'/(r['id']+'_'+n+'.json')
                child = (json.loads(child_path.read_text()) if child_path.exists() else
                         dict(r, target=str(graph/'tasks'/n), pending_only=False))
                ready &= live_control.stop_attempt(graph/'tasks'/n, job, child_path, child, allow_completed=True)
            if not ready:
                control.acknowledge(path, r, 'applying', detail='Waiting for affected jobs; unrelated work retained.')
                continue
            state.setdefault('protocol_changes', []).append(dict(request=r['id'], options=p['options'], reason=r['reason']))
            sharing.detach(graph, reason=r['reason'])
            atomic_json(graph/'protocol_preflight.json', dict(status='pending'))
            control.acknowledge(path, r, 'applying')
            live_control.commit(graph, current['replacement'], state, reset,
                                {n:r['reason'] for n in reset}, path)
        except (ValueError, RuntimeError) as error:
            transaction = graph/'.live_control_transaction.json'
            if transaction.exists() and json.loads(transaction.read_text()).get('status') == 'applying':
                raise
            control.acknowledge(path, r, 'blocked' if isinstance(error, ValueError) else r['status'], error=str(error))


def migrate(graph):
    """Add the new default branch without resetting existing native tasks."""
    from kinbot.anl import dispatch, anharmonic, live_control
    graph, spec, state = dispatch._load(graph)
    from kinbot.anl.refinement import require_open
    try:
        require_open(graph)
    except ValueError:
        return
    if spec.get('intent', {}).get('live_protocol', {}).get('options', {}).get('anharmonic_qz') is False:
        return
    revised = anharmonic.extend(spec)
    if revised != spec and revised['intent'].get('live_protocol'):
        eligible = revised['intent']['live_protocol']['eligible_tasks']
        eligible.extend(t['id'] for t in revised['tasks']
                        if t['id'] not in {x['id'] for x in spec['tasks']})
    if revised != spec:
        try:
            atomic_json(graph/'protocol_preflight.json', dict(status='pending'))
            live_control.commit(graph, revised, state, set(), {})
        except RuntimeError as error:
            if 'legacy workers' not in str(error):
                raise
            atomic_json(graph/'anharmonic_extension_status.json', dict(status='waiting_for_legacy_workers', error=str(error)))


def preflight_ready(graph):
    """New live branches receive the same site/Slurm checks as initial work."""
    from kinbot.anl import dispatch
    graph = Path(graph); path = graph/'protocol_preflight.json'
    if not path.exists():
        return True
    receipt = json.loads(path.read_text()); digest = file_sha256(graph/'workflow.json')
    if receipt.get('status') == 'verified' and receipt.get('spec_sha256') == digest:
        return True
    if receipt.get('retry_after', 0) > time.time():
        return False
    try:
        dispatch.preflight(graph, skip_optional=True)
    except dispatch.SchedulerUnavailableError:
        raise
    except (ValueError, RuntimeError) as error:
        atomic_json(path, dict(status='blocked', error=str(error), retry_after=time.time()+60))
        return False
    atomic_json(path, dict(status='verified', spec_sha256=digest))
    return True
