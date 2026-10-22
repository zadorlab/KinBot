"""Targeted VRC model migration used by the controller-refresh transaction.

Old scan directories stay at their original paths. Only owned relaxed-scan
allocations are retired; ANL, pivot and correction allocations are untouched.
"""
from copy import deepcopy
import json
from pathlib import Path
import re

from kinbot.execution import atomic_json, file_sha256


def settings(basis=None, fragments=None, pivots=None, reference_bases=None, reference_power=None,
             *, direct_caspt2_cbs=False, reviewed_models=None):
    values = {}
    if direct_caspt2_cbs:
        if reference_bases is not None or reference_power is not None:
            raise ValueError('Direct CASPT2 CBS cannot also change a CC reference.')
        values['direct_caspt2_cbs'] = True
    if reviewed_models is not None:
        path = Path(reviewed_models).expanduser().resolve(strict=True)
        values['reviewed_models'] = dict(path=str(path), sha256=file_sha256(path))
    if basis is not None:
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_+*().,-]*', basis):
            raise ValueError('Scan basis must be a basis name.')
        values['vrc_tst_scan_basis'] = basis
    if fragments is not None:
        if fragments not in ('legacy', 'asymptotic_scan'):
            raise ValueError('Invalid VRC fragment geometry mode.')
        values['vrc_tst_fragment_geometry'] = fragments
    if pivots is not None:
        if pivots not in ('L2', 'custom'):
            raise ValueError('Invalid VRC pivot level.')
        values['vrc_tst_pivot_level'] = pivots
    if reference_bases is not None:
        if len(reference_bases) not in (1, 2):
            raise ValueError('Reference needs one basis or an adjacent pair.')
        reference = {'basis': (reference_bases[0] if len(reference_bases) == 1
                               else list(reference_bases))}
        if reference_power is not None:
            if len(reference_bases) != 2:
                raise ValueError('Extrapolation power requires two reference bases.')
            reference['extrapolation_power'] = reference_power
        values['vrc_tst_spin_splitting'] = {'reference': reference}
    elif reference_power is not None:
        raise ValueError('Specify reference bases with their extrapolation power.')
    return values


def revised_input(payload, updates):
    revised = deepcopy(payload)
    def apply(target, selected):
        for key, value in selected.items():
            if key == 'vrc_tst_spin_splitting':
                spin = deepcopy(target.get(key) or {})
                if 'triplet_reference' in spin and 'reference' not in spin:
                    spin['reference'] = spin.pop('triplet_reference')
                reference = spin.setdefault('reference', {})
                reference.update(deepcopy(value['reference']))
                if isinstance(reference.get('basis'), str):
                    reference.pop('extrapolation_power', None)
                spin['reference'] = reference
                target[key] = spin
            else:
                target[key] = deepcopy(value)
    apply(revised, {k:v for k,v in updates.items() if k not in ('direct_caspt2_cbs', 'reviewed_models')})
    if revised.get('vrc_tst_surface_model_file'):
        raise ValueError('VRC migration requires an inline surface model; external files are not edited.')
    surfaces = [revised.get('vrc_tst_surface_model'),
                *revised.get('vrc_tst_surface_models', {}).values()]
    for surface in filter(None, surfaces):
        for channel in surface['channels']:
            for stage in channel['bottlenecks']:
                apply(stage['parameters'], {key: value for key, value in updates.items()
                                           if key in stage['parameters']})
    if updates.get('direct_caspt2_cbs'):
        # Derive each method from its own CAS. Never propagate the initial
        # well's CAS(5,4) into a later CAS(3,4) or CAS(7,5) channel.
        def direct(target, inherited):
            from kinbot.vrc_correction import apply_selected_active_space_methods
            effective = {**inherited, **target}
            method = effective.get('vrc_tst_scan_method', '')
            if not re.fullmatch(r'caspt2-f12\(\d+,\d+\)', method, re.I):
                raise ValueError('Direct CASPT2 migration requires explicit CASPT2-F12 scan methods.')
            effective['vrc_tst_high_method'] = method
            space = effective.get('vrc_tst_active_space')
            if isinstance(space, dict) and space.get('mode') == 'selected':
                apply_selected_active_space_methods(effective, space)
            target.update(vrc_tst_correction_model='direct',
                          vrc_tst_high_method=effective['vrc_tst_high_method'],
                          vrc_tst_high_basis=['avdz-f12', 'avtz-f12'],
                          vrc_tst_high_extrapolation_power=3.7)
            target.pop('vrc_tst_spin_splitting', None)
            target.pop('vrc_tst_reference_review', None)
        direct(revised, {})
        for surface in filter(None, surfaces):
            for channel in surface['channels']:
                for stage in channel['bottlenecks']:
                    direct(stage['parameters'], revised)

    return revised


def validate_input(path, updates, *, payload=None):
    """Resolve actual reviewed orbital records before stopping controllers."""
    from kinbot.parameters import Parameters
    from kinbot.vrc_multistate import load_surface_model
    from kinbot.vrc_correction import (correction_model, resolved_state_model,
                                       sampling_settings, spin_splitting_settings,
                                       trusted_settings, validate_vrc_levels)
    from kinbot.vrc_orbitals import fixed_pivot_profile, asymptotic_fragments
    from kinbot.molpro import vrc_multireference_energy_contract
    from kinbot.anl.site import _executable
    original = Parameters(str(path)).par
    if payload is not None:
        original.update(deepcopy(payload))
    par = revised_input(original, updates)
    surfaces = [load_surface_model(par)]
    surfaces.extend(load_surface_model(dict(par, vrc_tst_surface_model=model,
                    vrc_tst_surface_model_file=''))
                    for model in par.get('vrc_tst_surface_models', {}).values())
    leaves = []
    for surface in filter(None, surfaces):
        for channel in surface[0]['channels']:
            for stage in channel['bottlenecks']:
                for state in stage['states']:
                    leaf = {**par, **deepcopy(stage['parameters'])}
                    leaf['mult'] = state['multiplicity']
                    leaf['vrc_tst_selected_roots'] = {str(state['multiplicity']): state['root']}
                    leaves.append(leaf)
    if not leaves:
        leaves.append(par)
    for leaf in leaves:
        validate_vrc_levels(leaf)
        if vrc_multireference_energy_contract(leaf['vrc_tst_scan_method']) is None:
            raise ValueError('This migration requires a multireference relaxed scan.')
        envelope = leaf.get('vrc_tst_state_model') or {}
        if updates.get('direct_caspt2_cbs'):
            if leaf.get('vrc_tst_state_model_file'):
                model_path = Path(leaf['vrc_tst_state_model_file'])
                if not model_path.is_absolute():
                    model_path = Path(path).parent/model_path
                envelope = json.loads(model_path.read_text())
            if envelope.get('state_average_policy') != 'reviewed':
                raise ValueError('Preserving scans requires an explicit reviewed state average in every channel.')
        molecule = envelope.get('molecule', {})
        resolved_state_model(leaf, molecule.get('symbols'), molecule.get('charge', par['charge']),
                             None, basis=leaf['vrc_tst_scan_basis'])
        # Resolve every correction basis now too: a valid scan model alone
        # does not establish that the requested CBS pair has reviewed orbitals.
        levels = [sampling_settings(leaf)]
        if correction_model(leaf) == 'spin_splitting':
            spin = spin_splitting_settings(leaf)
            levels.extend([spin['reference'], *spin['gap_levels']])
        else:
            levels.append(trusted_settings(leaf))
        bases = {basis for level in levels for basis in
                 (level['basis'] if isinstance(level['basis'], list) else [level['basis']])}
        for basis in sorted(bases):
            # Accepted product spins are established by the live reaction
            # workflow. Do not impose default singlet/triplet coupling here.
            resolved_state_model(leaf, molecule.get('symbols'),
                                 molecule.get('charge', par['charge']), None, basis=basis)
        if asymptotic_fragments(leaf) and leaf.get('pp_orient', 'homo') == 'homo':
            profile = fixed_pivot_profile(leaf)
            if profile['calculator'] == 'gaussian':
                for program in ('formchk', 'cubegen'):
                    if _executable(program, 'gaussian') is None:
                        raise FileNotFoundError(f'Fixed L2 pivot orbitals require {program}.')


def plan(root, updates):
    from kinbot import run_salvage
    allowed = {'vrc_tst_scan_basis', 'vrc_tst_fragment_geometry', 'vrc_tst_pivot_level',
               'vrc_tst_spin_splitting', 'direct_caspt2_cbs', 'reviewed_models'}
    if 'direct_caspt2_cbs' in updates and updates['direct_caspt2_cbs'] is not True:
        raise ValueError('Direct CASPT2 migration must be explicitly enabled.')
    if updates.get('direct_caspt2_cbs') and 'vrc_tst_spin_splitting' in updates:
        raise ValueError('Cannot request direct and spin-splitting migrations together.')
    if not updates or set(updates) - allowed:
        raise ValueError('Unsupported VRC migration controls.')
    settings(updates.get('vrc_tst_scan_basis'), updates.get('vrc_tst_fragment_geometry'),
             updates.get('vrc_tst_pivot_level'))
    if 'vrc_tst_spin_splitting' in updates:
        spin = updates['vrc_tst_spin_splitting']
        if (not isinstance(spin, dict) or set(spin) != {'reference'}
                or not isinstance(spin['reference'], dict)
                or set(spin['reference']) - {'basis', 'extrapolation_power'}):
            raise ValueError('Reference migration can change only basis and extrapolation power.')
    record, inputs = run_salvage._cases(root)
    from kinbot import live_control, input_refresh
    if live_control.requests(root, kinds={'network', 'reaction'}):
        raise ValueError('Finish pending live input/reaction requests before a scientific migration.')
    for path in inputs:
        marker = path.parent/input_refresh.JOURNAL
        if marker.exists() and json.loads(marker.read_text()).get('status') in ('queued','publishing','blocked'):
            raise ValueError('Finish or resolve the pending input refresh before this migration.')
    proposal = None
    if 'reviewed_models' in updates:
        binding = updates['reviewed_models']
        proof = Path(binding['path'])
        if file_sha256(proof) != binding['sha256']:
            raise ValueError('Reviewed model proposal changed after staging.')
        proposal = json.loads(proof.read_text())
        if (proposal.get('schema') != 1 or proposal.get('kind') != 'orbital_model_proposal'
                or Path(proposal['case_input']) not in inputs):
            raise ValueError('Orbital proposal does not belong to this run.')
        if file_sha256(proposal['case_input']) != proposal['input_sha256']:
            raise ValueError('Case input changed; regenerate the orbital proposal before refreshing.')
        for name, digest in proposal['evidence'].items():
            if file_sha256(name) != digest:
                raise ValueError('Orbital proposal evidence changed: '+name)
        from kinbot.input_refresh import validate_models
        previous = json.loads(Path(proposal['case_input']).read_text())
        validate_models(previous, dict(previous, vrc_tst_surface_models=proposal['models']),
                        Path(proposal['case_input']).parent)
    paths = []
    def with_proposal(path, payload):
        if proposal and Path(proposal['case_input']).parent in (path.parent, path.parent.parent):
            payload = dict(payload, vrc_tst_surface_models=deepcopy(proposal['models']))
        return payload
    for path in inputs:
        if proposal and path == Path(proposal['case_input']):
            validate_input(path, updates, payload=with_proposal(path, json.loads(path.read_text())))
        else:
            validate_input(path, updates)
        paths.append(path)
        if (path.parent/'resolved_parameters.json').is_file():
            paths.append(path.parent/'resolved_parameters.json')
        for worker in sorted(path.parent.glob('[0-9]*')):
            if worker.is_dir():
                worker_input = worker/(worker.name+'.json')
                if worker_input.is_file():
                    paths.append(worker_input)
                    overrides=worker/'.kinbot_search_options.json'
                    if overrides.is_file():
                        paths.append(overrides)
    edits = []
    replace_geometry = False
    def geometry_controls(payload):
        keys = ('vrc_tst_scan_basis', 'vrc_tst_fragment_geometry')
        result = [tuple(payload.get(k) for k in keys)]
        # Added channels have no old scan to invalidate. Existing geometry
        # controls are unchanged by the direct-energy migration.
        for channel in (payload.get('vrc_tst_surface_model') or {}).get('channels', []):
            for stage in channel['bottlenecks']:
                result.append(tuple(stage['parameters'].get(k, payload.get(k)) for k in keys))
        return result
    for path in paths:
        if path.is_symlink() or not path.resolve().is_relative_to(root):
            raise ValueError(f'Input must be a regular file inside the run: {path}')
        previous = json.loads(path.read_text())
        sparse = path.name == '.kinbot_search_options.json'
        effective = previous
        if sparse:
            effective = {**json.loads((path.parent/(path.parent.name+'.json')).read_text()), **previous}
        if proposal:
            candidate = with_proposal(path, effective)
            if candidate != effective:
                validate_models(effective, candidate, path.parent)
        if path.name == path.parent.name+'.json' or sparse:
            native_input = path.parent/(path.parent.name+'.json')
            validate_input(native_input, updates, payload=with_proposal(path, effective))
        revised = revised_input(with_proposal(path, effective), updates)
        replace_geometry |= geometry_controls(effective) != geometry_controls(revised)
        if sparse:
            # kb.py accepts only live-search controls here, not a full input.
            revised = {key:revised[key] for key in previous}
            if proposal and with_proposal(path, {}) != {}:
                revised['vrc_tst_surface_models'] = revised_input(
                    with_proposal(path, effective), updates)['vrc_tst_surface_models']
        if revised != previous:
            edits.append({'path': str(path.relative_to(root)), 'before': file_sha256(path),
                          'after': revised})
    return {'settings': updates, 'inputs': edits, 'replace_geometry': replace_geometry,
            'prepared_before': record,
            'prepared_sha256': file_sha256(root/'prepared.json')}


def scan_inventory(root, jobs):
    """Only exact dispatcher-owned relaxed_scan jobs can be cancelled."""
    from kinbot.anl import dispatch
    live = {str(Path(j['directory']).resolve()): j for j in jobs}
    scans, selected = [], []
    for path in sorted(root.glob('*/[0-9]*/vrctst/molpro_scan/*/state.json')):
        graph, spec, state = dispatch._load(path.parent)
        for task in spec['tasks']:
            if task['id'] != 'relaxed_scan':
                continue
            entry = state['tasks'].get(task['id'], {})
            item = {'graph': str(graph), 'workflow_sha256': file_sha256(graph/'workflow.json'),
                    'status_at_retirement': entry.get('status'),
                    'distance_angstrom': spec.get('intent', {}).get('vrc_scan_distance_angstrom'),
                    'job_id': entry.get('job_id')}
            scans.append(item)
            job = live.get(str((graph/'tasks/relaxed_scan').resolve()))
            if job:
                if entry.get('status') != 'submitted' or str(entry.get('job_id')) != job['job_id']:
                    raise RuntimeError(f'{graph}: relaxed-scan job ownership changed.')
                selected.append(job)
    return scans, selected


def cancel_scans(plan, journal, marker):
    """Journal cancellations and require terminal scheduler evidence.

    Queue absence alone never authorizes new work. Repeating after a timeout
    rechecks the same IDs and paths; it never broadens the cancellation set.
    """
    from kinbot.execution import slurm_job_record
    from kinbot.run_refresh import _pending_request
    outcomes = plan.setdefault('cancellations', {})
    terminal = {'CANCELLED', 'COMPLETED', 'FAILED', 'TIMEOUT', 'NODE_FAIL',
                'OUT_OF_MEMORY', 'PREEMPTED', 'BOOT_FAIL', 'DEADLINE', 'REVOKED'}
    for job in plan['selected_jobs']:
        ident = job['job_id']
        if outcomes.get(ident, {}).get('phase') == 'terminal':
            continue
        current = slurm_job_record(ident)
        if (not ident.isdecimal() or current.get('JobId') != ident
                or Path(current.get('WorkDir', '/')).resolve() != Path(job['directory']).resolve()):
            raise RuntimeError(f'Selected job {ident}: scheduler identity unavailable; repeat this refresh.')
        outcome = outcomes.setdefault(ident, {'before': current, 'phase': 'inspected'})
        if current.get('JobState') in terminal:
            outcome.update(phase='terminal', record=current)
            atomic_json(marker, journal)
            continue
        if current.get('JobState') not in ('PENDING', 'RUNNING', 'CONFIGURING', 'SUSPENDED', 'COMPLETING'):
            raise RuntimeError(f'Selected job {ident}: unrecognized state; no cancellation attempted.')
        if current['JobState'] != 'COMPLETING':
            outcome['phase'] = 'requested'
            atomic_json(marker, journal)
            outcome['request_diagnostic'] = _pending_request(['scancel', ident])
            atomic_json(marker, journal)
            current = slurm_job_record(ident)
        # Slurm normally needs a short cleanup interval. Do not make a
        # successful cancellation look like an immediate handoff failure.
        import time
        deadline = time.monotonic() + 30
        waiting = False
        while (current.get('JobId') == ident
               and Path(current.get('WorkDir', '/')).resolve() == Path(job['directory']).resolve()
               and current.get('JobState') in ('RUNNING', 'PENDING', 'COMPLETING', 'CONFIGURING')
               and time.monotonic() < deadline):
            if not waiting:
                print(f'Waiting for superseded job {ident} to release its allocation...', flush=True)
                waiting = True
            time.sleep(2)
            current = slurm_job_record(ident)
        if (current.get('JobId') != ident or current.get('JobState') not in terminal
                or Path(current.get('WorkDir', '/')).resolve() != Path(job['directory']).resolve()):
            raise RuntimeError(f'Selected job {ident}: cancellation/finish is not yet confirmed. '
                               'Other allocations remain running; repeat this same refresh command.')
        outcome.update(phase='terminal', record=current)
        atomic_json(marker, journal)


def apply_inputs(root, plan, archive, *, record_name='vrc_restart.json'):
    """Idempotent bytes-first migration; native graph files are never edited."""
    from kinbot import run_salvage
    import hashlib
    for entry in plan['inputs']:
        path = root/entry['path']
        after = hashlib.sha256((json.dumps(entry['after'], indent=2, sort_keys=True)+'\n').encode()).hexdigest()
        actual = file_sha256(path)
        if actual == after:
            continue
        if actual != entry['before']:
            raise RuntimeError(f'{path}: input changed during VRC migration; refusing overwrite.')
        run_salvage._save(path, root, archive)
        atomic_json(path, entry['after'])
    prepared = deepcopy(plan['prepared_before'])
    prepared['sha256'] = {name: file_sha256(root/name) for name in prepared['sha256']}
    marker = root/'prepared.json'
    if json.loads(marker.read_text()) != prepared:
        if file_sha256(marker) != plan['prepared_sha256']:
            raise RuntimeError('Prepared manifest changed during VRC migration.')
        run_salvage._save(marker, root, archive)
        atomic_json(marker, prepared)
    # Graphs/checkpoints remain at their original locations, including failed
    # last evaluations. They are comparison evidence, not new-basis results.
    atomic_json(root/record_name, {
        'schema': 1, 'settings': plan['settings'], 'archive': str(archive),
        'retired_scans': plan['scans'], 'cancellations': plan.get('cancellations', {}),
        'note': 'Old scan artifacts retained in place. Cancelled jobs lose in-memory optimizer state.'})


def audit_direct_inputs(root):
    """Verify current input contracts after handoff, without reading QC tails."""
    from kinbot import run_salvage
    root=Path(root).resolve()
    _,inputs=run_salvage._cases(root)
    paths=[]
    for path in inputs:
        paths.append(path)
        resolved=path.parent/'resolved_parameters.json'
        if resolved.is_file(): paths.append(resolved)
        paths.extend(p for p in path.parent.glob('[0-9]*/*.json')
                     if p.name==p.parent.name+'.json')
    leaves=0
    for path in paths:
        payload=json.loads(path.read_text())
        override=path.parent/'.kinbot_search_options.json'
        if path.name == path.parent.name+'.json' and override.exists():
            payload.update(json.loads(override.read_text()))
        settings=[payload]
        for surface in [payload.get('vrc_tst_surface_model'),
                        *payload.get('vrc_tst_surface_models',{}).values()]:
            if surface:
                settings.extend({**payload,**s['parameters']} for c in surface['channels'] for s in c['bottlenecks'])
        for par in settings:
            if (par.get('vrc_tst_correction_model')!='direct'
                    or par.get('vrc_tst_spin_splitting')
                    or par.get('vrc_tst_high_basis')!=['avdz-f12','avtz-f12']
                    or par.get('vrc_tst_high_extrapolation_power')!=3.7
                    or not re.fullmatch(r'caspt2-f12\(\d+,\d+\)',par.get('vrc_tst_high_method',''),re.I)):
                raise RuntimeError(f'Direct CBS input audit failed: {path}')
            leaves+=1
    return dict(input_files=len(paths),parameter_blocks=leaves,
                correction='direct CASPT2-F12 avdz-f12/avtz-f12 CBS, power 3.7')
