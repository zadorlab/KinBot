"""Targeted VRC model migration used by the controller-refresh transaction.

Old scan directories stay at their original paths. Only owned relaxed-scan
allocations are retired; ANL, pivot and correction allocations are untouched.
"""
from copy import deepcopy
import json
from pathlib import Path
import re

from kinbot.execution import atomic_json, file_sha256


def settings(basis=None, fragments=None, pivots=None):
    values = {}
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
    return values


def revised_input(payload, updates):
    revised = deepcopy(payload)
    revised.update(updates)
    if revised.get('vrc_tst_surface_model_file'):
        raise ValueError('VRC migration requires an inline surface model; external files are not edited.')
    for channel in (revised.get('vrc_tst_surface_model') or {}).get('channels', []):
        for stage in channel['bottlenecks']:
            # Explicit stage overrides otherwise defeat top-level controls.
            for key, value in updates.items():
                if key in stage['parameters']:
                    stage['parameters'][key] = value
    return revised


def validate_input(path, updates):
    """Resolve actual reviewed orbital records before stopping controllers."""
    from kinbot.parameters import Parameters
    from kinbot.vrc_multistate import load_surface_model
    from kinbot.vrc_correction import resolved_state_model, validate_vrc_levels
    from kinbot.vrc_orbitals import fixed_pivot_profile, asymptotic_fragments
    from kinbot.molpro import vrc_multireference_energy_contract
    from kinbot.anl.site import _executable
    par = revised_input(Parameters(str(path)).par, updates)
    surface = load_surface_model(par)
    leaves = []
    if surface:
        for channel in surface[0]['channels']:
            for stage in channel['bottlenecks']:
                for state in stage['states']:
                    leaf = {**par, **deepcopy(stage['parameters'])}
                    leaf['mult'] = state['multiplicity']
                    leaf['vrc_tst_selected_roots'] = {str(state['multiplicity']): state['root']}
                    leaves.append(leaf)
    else:
        leaves.append(par)
    for leaf in leaves:
        validate_vrc_levels(leaf)
        if vrc_multireference_energy_contract(leaf['vrc_tst_scan_method']) is None:
            raise ValueError('This migration requires a multireference relaxed scan.')
        envelope = leaf.get('vrc_tst_state_model') or {}
        molecule = envelope.get('molecule', {})
        resolved_state_model(leaf, molecule.get('symbols'), molecule.get('charge', par['charge']),
                             None, basis=leaf['vrc_tst_scan_basis'])
        if asymptotic_fragments(leaf) and leaf.get('pp_orient', 'homo') == 'homo':
            profile = fixed_pivot_profile(leaf)
            if profile['calculator'] == 'gaussian':
                for program in ('formchk', 'cubegen'):
                    if _executable(program, 'gaussian') is None:
                        raise FileNotFoundError(f'Fixed L2 pivot orbitals require {program}.')


def plan(root, updates):
    from kinbot import run_salvage
    allowed = {'vrc_tst_scan_basis', 'vrc_tst_fragment_geometry', 'vrc_tst_pivot_level'}
    if not updates or set(updates) - allowed:
        raise ValueError('Only scan basis, fragment geometry and pivot level may change in this handoff.')
    settings(updates.get('vrc_tst_scan_basis'), updates.get('vrc_tst_fragment_geometry'),
             updates.get('vrc_tst_pivot_level'))
    record, inputs = run_salvage._cases(root)
    paths = []
    for path in inputs:
        validate_input(path, updates)
        paths.append(path)
        if (path.parent/'resolved_parameters.json').is_file():
            paths.append(path.parent/'resolved_parameters.json')
        for worker in sorted(path.parent.glob('[0-9]*')):
            if worker.is_dir():
                worker_input = worker/(worker.name+'.json')
                if worker_input.is_file():
                    paths.append(worker_input)
    edits = []
    replace_geometry = False
    def geometry_controls(payload):
        keys = ('vrc_tst_scan_basis', 'vrc_tst_fragment_geometry')
        result = [tuple(payload.get(k) for k in keys)]
        for channel in (payload.get('vrc_tst_surface_model') or {}).get('channels', []):
            for stage in channel['bottlenecks']:
                result.append(tuple(stage['parameters'].get(k, payload.get(k)) for k in keys))
        return result
    for path in paths:
        if path.is_symlink() or not path.resolve().is_relative_to(root):
            raise ValueError(f'Input must be a regular file inside the run: {path}')
        previous = json.loads(path.read_text())
        revised = revised_input(previous, updates)
        replace_geometry |= geometry_controls(previous) != geometry_controls(revised)
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
            raise RuntimeError(f'VRC job {ident}: scheduler identity unavailable; repeat this refresh.')
        outcome = outcomes.setdefault(ident, {'before': current, 'phase': 'inspected'})
        if current.get('JobState') in terminal:
            outcome.update(phase='terminal', record=current)
            atomic_json(marker, journal)
            continue
        if current.get('JobState') not in ('PENDING', 'RUNNING', 'CONFIGURING', 'SUSPENDED', 'COMPLETING'):
            raise RuntimeError(f'VRC job {ident}: unrecognized state; no cancellation attempted.')
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
                print(f'Waiting for superseded VRC job {ident} to release its allocation...', flush=True)
                waiting = True
            time.sleep(2)
            current = slurm_job_record(ident)
        if (current.get('JobId') != ident or current.get('JobState') not in terminal
                or Path(current.get('WorkDir', '/')).resolve() != Path(job['directory']).resolve()):
            raise RuntimeError(f'VRC job {ident}: cancellation/finish is not yet confirmed. '
                               'Other allocations remain running; repeat this same refresh command.')
        outcome.update(phase='terminal', record=current)
        atomic_json(marker, journal)


def apply_inputs(root, plan, archive):
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
    atomic_json(root/'vrc_restart.json', {
        'schema': 1, 'settings': plan['settings'], 'archive': str(archive),
        'retired_scans': plan['scans'], 'cancellations': plan.get('cancellations', {}),
        'note': 'Old scan artifacts retained in place. Cancelled jobs lose in-memory optimizer state.'})
