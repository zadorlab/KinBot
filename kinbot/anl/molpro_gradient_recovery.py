"""One numerical-accuracy retry for an exhausted, stagnating optimization.

Retain the selected energy, electronic model, constraint targets and force
criterion. A small trust radius alone is not evidence of this failure class:
require an exhausted step budget and a persistent force plateau in Sella's
native progress log. Unknown failures stay available for review.
"""
from copy import deepcopy
import json
import math
from pathlib import Path
import uuid


def evidence(task, directory, outcome):
    profile = task.get('profile', {})
    kwargs = profile.get('calculator_kwargs', {})
    from kinbot.molpro import vrc_multireference_energy_contract
    contract = vrc_multireference_energy_contract(profile.get('method', ''))
    analytic = (contract and contract['method'] == 'CASPT2'
                and not contract['explicitly_correlated']
                and kwargs.get('gradient_mode', 'auto') != 'numerical')
    if (task.get('kind') != 'ase_optimize' or profile.get('calculator') != 'molpro'
            or profile.get('optimizer') != 'sella' or analytic
            or outcome.get('status') != 'failed'
            or outcome.get('failure_phase') != 'execution'
            or outcome.get('error') != 'ASE geometry optimization did not converge.'):
        return None
    try:
        budget = json.loads((directory/'optimization_budget.json').read_text())
        with (directory/'optimization.log').open('rb') as stream:
            stream.seek(max(0, (directory/'optimization.log').stat().st_size-32_000))
            lines = stream.read().decode().splitlines()
        rows = [s.split() for s in lines if s.startswith('Sella ')][-20:]
        if len(rows) != 20 or any(len(row) != 9 for row in rows):
            return None
        steps = [int(row[1]) for row in rows]
        forces, radii = ([float(row[i]) for row in rows] for i in (4, 6))
        target = task['optimizer']['fmax']
        if (budget.get('status') != 'exhausted' or budget.get('fmax') != target
                or type(budget.get('chunks_used')) is not int or budget['chunks_used'] < 1
                or steps != list(range(steps[0], steps[0]+20))
                or not all(math.isfinite(v) for v in [*forces, *radii])
                or min(forces) <= target or not 0 < min(radii) <= max(radii) <= .001
                or min(forces[-5:]) < .8 * min(forces[:5])):
            return None
    except (OSError, ValueError, KeyError, TypeError):
        return None
    return dict(budget=budget, first_step=steps[0], last_step=steps[-1],
                force_range=[min(forces), max(forces)], maximum_trust_radius=max(radii))


def accuracy_options(task):
    """Tighten only numerical tolerances; never replace a model definition."""
    before = task['profile'].get('calculator_kwargs', {})
    after = deepcopy(before)
    after['scf_convergence'] = min(before.get('scf_convergence', 1.e-10), 1.e-11)
    after['correlation_convergence'] = min(before.get('correlation_convergence', 1.e-8), 1.e-10)
    after['numerical_gradient'] = {**(before.get('numerical_gradient') or {}),
                                   'stencil': 'fourpoint'}
    model = after.get('state_model')
    if model is not None:
        for group, limits in (
                ('casscf_options', {'energy_convergence': 1.e-11, 'gradient_convergence': 1.e-7}),
                ('correlation_options', {'energy_convergence': 1.e-11, 'coefficient_convergence': 1.e-8})):
            values = model.setdefault(group, {})
            for key, limit in limits.items():
                values[key] = min(values.get(key, limit), limit)
    return None if before == after else {'schema': 1, 'calculator_kwargs': after}


def recover(graph, spec, state):
    from ase.io import read, write
    from kinbot.anl import dispatch, live_control
    from kinbot.anl.optimization_recovery import _last_point
    from kinbot.anl.salvage import descendants
    from kinbot.anl.timeout_recovery import _job_id
    from kinbot.execution import file_sha256
    from kinbot.manual_fallback import task_options
    graph = Path(graph)
    if spec.get('intent', {}).get('vrc_fragment_retry'):
        return spec, state
    replacement, updated = deepcopy(spec), deepcopy(state)
    reset, reasons = set(), {}
    for index, task in enumerate(spec['tasks']):
        ident = task['id']
        entry = state['tasks'].get(ident, {})
        if (entry.get('status') not in ('failed', 'unavailable')
                or state.get('molpro_gradient_retries', {}).get(ident)):
            continue
        directory = graph/'tasks'/ident
        path = directory/'execution.json'
        if not path.is_file():
            continue
        outcome = json.loads(path.read_text())
        proof = evidence(task, directory, outcome)
        if not proof:
            continue
        options = accuracy_options(task)
        if not options:
            continue  # These settings have already been exhausted.
        job = _job_id(directory, ident, entry)
        if not job or dispatch._job_active(job):
            continue
        dispatch._verify_stage_files(graph, task, entry)
        dispatch._verify_execution(graph, task, entry, outcome)
        children = descendants(spec, ident)
        if any(state['tasks'].get(n, {}).get('status') in
               ('complete', 'submitted', 'submitting', 'staged') for n in children):
            continue
        # Resume positions only: recalculate energy/forces with a fresh
        # optimizer. The original geometry still defines constraint targets.
        last = _last_point(directory, task, read(directory/'geometry.xyz'))
        replacement['tasks'][index] = task_options(task, options)
        dispatch.validate_spec(replacement)
        seed = graph/'attempts'/('gradient_'+uuid.uuid4().hex)/'continuation.xyz'
        seed.parent.mkdir(parents=True)
        write(seed, last)
        updated.setdefault('optimizer_continuations', {})[ident] = dict(
            source=str(seed), sha256=file_sha256(seed), chunks_used=0,
            source_trajectory_sha256=file_sha256(directory/'optimization.traj'))
        proof.update(previous_job_id=job, numerical_options=options,
                     optimizer_log_sha256=file_sha256(directory/'optimization.log'))
        updated.setdefault('molpro_gradient_retries', {})[ident] = [proof]
        reset.add(ident)
        reset.update(n for n in children if state['tasks'].get(n, {}).get('status') in ('blocked', 'skipped'))
        reasons[ident] = ('One tighter electronic/four-point gradient retry after exhausted '
                         'small-trust-radius stagnation; model, constraints and fmax unchanged.')
    if reset:
        live_control.commit(graph, replacement, updated, reset, reasons)
        _, spec, state = dispatch._load(graph)
    return spec, state
