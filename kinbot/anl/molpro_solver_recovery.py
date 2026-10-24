"""Bounded native-diagnostic solver retries, preserving the electronic model.

PSPACE changes the primary space used by the CI eigensolver, not the full
configuration/reference space, requested states, or CASPT2 Hamiltonian.
Unknown failures and an explicitly configured P-space are left for review.
"""
from copy import deepcopy
import json
import re


def solver_options(model, text, attempt):
    """Shared two-step hsdel ladder for optimization and sampling evaluations."""
    if (attempt not in (0, 1) or not model.get('state_average')
            or ('pspace' in model.get('correlation_options', {}) and attempt == 0)
            or not re.search(r'Catastrophic failure in diagonalization\s*\(hsdel\)', text, re.I)
            or not re.search(r'The expansion set has become singular', text, re.I)):
        return None
    result = deepcopy(model)
    controls = result.setdefault('correlation_options', {})
    roots = max(s['roots'] for s in result['state_average'])
    counts = [int(x) for x in re.findall(r'Number of p-space configurations:\s*(\d+)', text, re.I)]
    if attempt == 0:
        count = 2 * max([roots, *counts])
        if count > 64:
            return None
        controls['pspace'] = dict(energy_threshold=0., minimum_configurations=count)
    else:
        floor = 2 * roots + 2
        old = controls.get('max_reference_vectors')
        controls['max_reference_vectors'] = max(floor, old // 2) if old else floor
        if old is not None:
            controls['max_reference_vectors'] = min(old, controls['max_reference_vectors'])
        controls['energy_convergence'] = min(controls.get('energy_convergence', 1.e-7), 1.e-10)
        controls['coefficient_convergence'] = min(controls.get('coefficient_convergence', 1.e-5), 1.e-7)
    return result if result != model else None


def evidence(task, directory, outcome, *, own_pspace=False):
    from kinbot.molpro import vrc_multireference_energy_contract
    profile = task.get('profile', {})
    model = profile.get('calculator_kwargs', {}).get('state_model') or {}
    if (task.get('kind') != 'ase_optimize' or profile.get('calculator') != 'molpro'
            or vrc_multireference_energy_contract(profile.get('method', '')) is None
            or not model.get('state_average')
            or ('pspace' in model.get('correlation_options', {}) and not own_pspace)
            or outcome.get('status') != 'failed'
            or outcome.get('failure_phase') != 'execution'):
        return None
    match = re.fullmatch(r'Molpro force evaluation (\w+) exited with status ([1-9]\d*); '
                         r'diagnostics: (\w+\.failure\.json)\.', outcome.get('error', ''))
    if not match or match[3] != match[1] + '.failure.json':
        return None
    try:
        record = json.loads((directory/match[3]).read_text())
    except (OSError, ValueError):
        return None
    if record.get('returncode') != int(match[2]) or record.get('explicit_io_failure'):
        return None
    text = '\n'.join(v.get('error_context', '') + '\n' + v.get('tail', '')
                     for v in record.get('files', {}).values())
    # Read only this failed evaluation. Native logs may retain the P-space
    # count outside the diagnostic excerpt; never infer it from another step.
    for suffix in ('.out', '.log'):
        path = directory/(match[1]+suffix)
        if path.is_file():
            with path.open('rb') as stream:
                stream.seek(max(0, path.stat().st_size-8_000_000))
                text += '\n' + stream.read().decode(errors='replace')
    if not (re.search(r'Catastrophic failure in diagonalization\s*\(hsdel\)', text, re.I)
            and re.search(r'The expansion set has become singular', text, re.I)):
        return None
    counts = [int(x) for x in re.findall(r'Number of p-space configurations:\s*(\d+)', text, re.I)]
    roots = max(s['roots'] for s in model['state_average'])
    # Small bounded primary solve. Do not silently enlarge a large native
    # eigensystem, nor override a deliberate user P-space choice.
    count = 2 * max([roots, *counts])
    if count > 64:
        return None
    return dict(file=match[3], host=record.get('host'), returncode=record['returncode'],
                observed_pspace_counts=sorted(set(counts)),
                pspace={'energy_threshold': 0., 'minimum_configurations': count})


def recover(graph, spec, state):
    from kinbot.anl import dispatch, live_control
    from kinbot.anl.salvage import descendants
    from kinbot.anl.timeout_recovery import _job_id
    from kinbot.manual_fallback import task_options
    replacement, updated = deepcopy(spec), deepcopy(state)
    reset, reasons = set(), {}
    for index, task in enumerate(replacement['tasks']):
        ident = task['id']
        entry = state['tasks'].get(ident, {})
        history = state.get('molpro_pspace_retries', {}).get(ident, [])
        if (entry.get('status') not in ('failed', 'unavailable') or len(history) >= 2):
            continue
        directory = graph/'tasks'/ident
        path = directory/'execution.json'
        if not path.is_file():
            continue
        outcome = json.loads(path.read_text())
        proof = evidence(task, directory, outcome, own_pspace=bool(history))
        if not proof:
            continue
        job = _job_id(directory, ident, entry)
        if not job or dispatch._job_active(job):
            continue
        dispatch._verify_stage_files(graph, task, entry)
        dispatch._verify_execution(graph, task, entry, outcome)
        children = descendants(spec, ident)
        if any(state['tasks'].get(n, {}).get('status') in
               ('complete', 'submitted', 'submitting', 'staged') for n in children):
            continue
        old_model = task['profile']['calculator_kwargs']['state_model']
        # evidence already validated the native signature/count. Use the same
        # numerical policy as the sampling bridge, preserving all other fields.
        native = ('Catastrophic failure in diagonalization(hsdel)\nThe expansion set has become singular\n'
                  + '\n'.join(f'Number of p-space configurations: {n}' for n in proof['observed_pspace_counts']))
        model = solver_options(old_model, native, len(history))
        if model is None:
            continue
        controls = model['correlation_options']
        kwargs = {'state_model': model}
        if not history:
            reason = 'P-space eigensolver retry after native hsdel singularity'
        else:
            # MAXVI is the iterative workspace, not NSTATI or the reference
            # configuration space. A smaller workspace limits dependence in
            # accumulated vectors; tighter tolerances reduce numerical noise.
            kwargs['scf_convergence'] = min(task['profile']['calculator_kwargs'].get('scf_convergence', 1.e-10), 1.e-10)
            reason = 'Final hsdel retry with reduced internal expansion workspace and tighter convergence'
            proof['pspace'] = controls['pspace']
        proof['numerical_controls'] = deepcopy(controls)
        replacement['tasks'][index] = task_options(task, {'schema': 1,
            'calculator_kwargs': kwargs})
        reset.add(ident)
        reset.update(n for n in children if state['tasks'].get(n, {}).get('status')
                     in ('blocked', 'skipped'))
        updated.setdefault('molpro_pspace_retries', {})[ident] = [
            *history, dict(proof, previous_job_id=job)]
        reasons[ident] = reason + '; active space, roots, weights and Hamiltonian unchanged.'
    if reset:
        live_control.commit(graph, replacement, updated, reset, reasons)
        _, spec, state = dispatch._load(graph)
    return spec, state
