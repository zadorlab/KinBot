"""Live, dependency-aware geometry and harmonic-frequency refinement."""
from copy import deepcopy
import json
import math
from pathlib import Path
import re

from kinbot.execution import atomic_json, file_sha256


def require_open(graph):
    marker = Path(graph)/'refinement_window.json'
    if marker.exists() and json.loads(marker.read_text()).get('status') == 'closed':
        raise ValueError('Graph results have been handed downstream; reopen their owning workflow before refinement.')


def window(graph, *, close=False):
    """Serialize final handoff with publication of a last-moment live request."""
    from contextlib import nullcontext
    from kinbot import live_control as control
    graph = Path(graph)
    root = control.find(graph)
    with control.locked(root.parent) if root else nullcontext():
        from kinbot.anl import protocol
        if close and protocol.requests(graph):
            return False
        if close and control.requests(graph/'tasks', kinds={'task', 'refinement'}):
            return False
        atomic_json(graph/'refinement_window.json', dict(status='closed' if close else 'open'))
    return True


def validate(options):
    from kinbot.theory import convergence_exponent
    if not isinstance(options, dict) or set(options) != {'fmax', 'scf', 'cc'}:
        raise ValueError('Refinement requires fmax, scf and cc.')
    if (type(options['fmax']) not in (int, float)
            or not math.isfinite(options['fmax']) or options['fmax'] <= 0):
        raise ValueError('Refinement fmax must be positive and finite (eV/angstrom).')
    convergence_exponent(options['scf'], 'Refinement SCF convergence')
    convergence_exponent(options['cc'], 'Refinement CC convergence')
    return deepcopy(options)


def _generation(entry):
    # Status may change when an already observed job finishes. A new
    # submission, geometry or task generation must never expand this request.
    return {k: entry.get(k) for k in
            ('attempt', 'job_id', 'task_sha256', 'geometry_sha256')}


def _frequency_task(task, options):
    from kinbot.manual_fallback import task_options
    template = task['input_template']
    # Explicit native controls are required; do not guess defaults or alter
    # arbitrary user procedures. Preserve every stricter existing threshold.
    global_pattern = (r'(?im)^gthresh,energy=1\.d-(\d+),coeff=1\.d-(\d+),'
                      r'civec=1\.d-(\d+)[ \t]*$')
    cc_pattern = r'(?i)thrden=1\.d-(\d+),thrvar=1\.d-(\d+)'
    if len(re.findall(global_pattern, template)) != 1 or len(re.findall(cc_pattern, template)) != 1:
        raise ValueError('Frequency refinement requires explicit generated Molpro convergence controls.')
    se, ce = (round(-math.log10(options[k])) for k in ('scf', 'cc'))
    result = re.sub(global_pattern, lambda m: 'gthresh,' + ','.join(
        f'{k}=1.d-{max(int(m[i]), e)}' for i, k, e in
        ((1, 'energy', se), (2, 'coeff', ce), (3, 'civec', ce))), template)
    result = re.sub(cc_pattern, lambda m:
        f'thrden=1.d-{max(int(m[1]), ce)},thrvar=1.d-{max(int(m[2]), ce, 10)}', result)
    return task_options(task, {'schema': 1, 'input_edits': [
        {'old': template, 'new': result}]})


def plan(graph, ident, options):
    """Read-only exact scope; no source or scheduler mutation."""
    from kinbot.anl import dispatch
    from kinbot.anl.salvage import descendants
    from kinbot.manual_fallback import task_options
    options = validate(options)
    graph, spec, state = dispatch._load(graph)
    require_open(graph)
    tasks = {t['id']: t for t in spec['tasks']}
    target = tasks.get(ident, {})
    if target.get('kind') == 'ase_optimize':
        geometry = ident
    elif target.get('geometry_from') in tasks:
        geometry = target.get('geometry_from')
    else:
        raise ValueError('Select an ASE geometry or a task consuming an optimized geometry.')
    if geometry not in tasks or tasks[geometry]['kind'] != 'ase_optimize':
        raise ValueError('Frequency task has no optimizable geometry dependency.')
    reset = {geometry} | descendants(spec, geometry)
    if not any(state['tasks'].get(n, {}).get('status') in ('failed', 'unavailable')
               for n in reset):
        raise ValueError('Refinement requires a failed/unavailable task in the affected dependency chain.')
    if any(state['tasks'].get(n, {}).get('status') == 'submitting' for n in reset):
        raise ValueError('Reconcile ambiguous submissions before refinement.')
    replacement = deepcopy(spec)
    changed = []
    for i, task in enumerate(spec['tasks']):
        if task['id'] not in reset:
            continue
        if task['kind'] == 'ase_optimize':
            if task['profile']['calculator'] != 'molpro':
                raise ValueError('Paired refinement currently supports Molpro ASE geometry/harmonic graphs.')
            old = task.get('optimizer', {}).get('fmax', .03)
            if options['fmax'] > old:
                raise ValueError('Refinement cannot loosen an existing geometry force criterion.')
            kwargs = task['profile'].get('calculator_kwargs', {})
            values = {'scf_convergence': min(options['scf'], kwargs.get('scf_convergence', 1.e-10)),
                      'correlation_convergence': min(options['cc'], kwargs.get('correlation_convergence', 1.e-8))}
            patched = task_options(task, {'schema': 1, 'calculator_kwargs': values,
                                          'optimizer': {'fmax': options['fmax']}})
        elif task.get('result_parser', {}).get('kind') == 'molpro_harmonic':
            patched = _frequency_task(task, options)
        else:
            continue  # Recompute dependent energies at their declared method/settings.
        replacement['tasks'][i] = patched
        changed.append(task['id'])
    dispatch.validate_spec(replacement)
    return dict(graph=str(graph), target=ident, geometry=geometry, options=options,
                spec_sha256=state['spec_sha256'], reset=sorted(reset), changed=changed,
                expected={n: _generation(state['tasks'].get(n, {})) for n in sorted(reset)},
                affected=[dict(task=n, status=state['tasks'].get(n, {}).get('status', 'waiting'),
                               job_id=state['tasks'].get(n, {}).get('job_id')) for n in sorted(reset)],
                replacement=replacement)


def pending(graph):
    """Fence only affected tasks while their owner stops/waits for old jobs."""
    from kinbot import live_control as control
    from kinbot.anl import protocol
    result = protocol.pending(graph)
    for _, request in control.requests(Path(graph)/'tasks', kinds={'refinement'}):
        planned = request.get('options', {}).get('refinement_plan')
        if planned and planned.get('graph') == str(Path(graph).resolve()):
            result.update(planned['reset'])
    return result


def process(graph, path, request, spec, state):
    """Called by the owner with the graph lock; journal before replacing tasks."""
    from kinbot import live_control as control
    from kinbot.anl import dispatch, live_control
    planned = request['options']['refinement_plan']
    if str(graph) != planned['graph'] or state['spec_sha256'] != planned['spec_sha256']:
        raise ValueError('Workflow changed after refinement review; preview a new request.')
    reset = set(planned['reset'])
    if any(_generation(state['tasks'].get(n, {})) != e for n, e in planned['expected'].items()):
        raise ValueError('Affected task generation changed after refinement review; no replacement applied.')
    # Rebuild from validated user settings, not a serialized replacement spec.
    current = plan(graph, planned['target'], planned['options'])
    if current['reset'] != planned['reset']:
        raise ValueError('Refinement dependency scope changed.')
    # Verify all existing source/task files before sending any cancellation.
    for n in reset:
        entry = state['tasks'].get(n, {})
        if entry.get('task_sha256'):
            task = next(t for t in spec['tasks'] if t['id'] == n)
            dispatch._verify_stage_files(graph, task, entry)
            outcome = graph/'tasks'/n/'execution.json'
            if outcome.exists():
                dispatch._verify_execution(graph, task, entry, json.loads(outcome.read_text()))
    for n, entry in state['tasks'].items():
        if n not in reset and entry.get('status') in ('submitted', 'submitting'):
            record = json.loads((graph/'tasks'/n/'task.json').read_text())
            if record.get('live_control_schema') != 1:
                raise RuntimeError('Refinement deferred until legacy sibling workers finish.')
    ready = True
    for n in sorted(reset):
        entry = state['tasks'].get(n, {})
        job = entry.get('job_id')
        if not job or not dispatch._job_active(job):
            continue
        if not request['cancel_running']:
            ready = False
            continue
        # Each job gets its own cancellation receipt; never use one job's
        # cancellation flag to skip another job in the affected DAG.
        child = request.setdefault('dependent_cancellations', {}).setdefault(
            n, dict(request, target=str(graph/'tasks'/n)))
        child.pop('dependent_cancellations', None)
        child_path = path.parent/(request['id']+'_'+n+'.cancel.json')
        # Keep scheduler subreceipts outside the queue read by requests().
        child_path = path.parent.parent/'refinement_cancellations'/child_path.name
        if not live_control.stop_attempt(graph/'tasks'/n, job, child_path, child, allow_completed=True):
            ready = False
        atomic_json(path, request)
    if not ready:
        control.acknowledge(path, request, 'applying', detail='Waiting for affected allocations; independent work retained.')
        return
    updated = deepcopy(state)
    for n in reset:
        updated.get('optimizer_continuations', {}).pop(n, None)
    # Reuse a verified converged starting point, but recompute forces/Hessian
    # with the requested accuracy. Descendant geometries start from the NEW
    # parent, never their old continuation seeds.
    geom = planned['geometry']
    ge = state['tasks'].get(geom, {})
    seed = graph/'tasks'/geom/'final.xyz'
    if ge.get('status') == 'complete' and seed.is_file():
        from ase.io import read
        if dispatch._geometry_hash(read(seed)) != ge.get('final_geometry_sha256'):
            raise RuntimeError('Refinement seed geometry changed.')
        saved = graph/'refinement_seeds'/(request['id']+'.xyz')
        saved.parent.mkdir(exist_ok=True)
        saved.write_bytes(seed.read_bytes())
        updated.setdefault('optimizer_continuations', {})[geom] = dict(
            source=str(saved), sha256=file_sha256(saved), chunks_used=0,
            reason='Explicit tighter geometry/frequency refinement', request=request['id'])
    updated.setdefault('refinements', []).append(dict(request=request['id'],
        options=planned['options'], affected=planned['reset'], reason=request['reason']))
    control.acknowledge(path, request, 'applying')
    live_control.commit(graph, current['replacement'], updated, reset,
                        {n: request['reason'] for n in reset}, path)
