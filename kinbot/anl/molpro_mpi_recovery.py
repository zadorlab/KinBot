"""One serial diagnostic retry for a recorded Molpro communication abort.

This is a resource fallback, not a diagnosis of the MPI implementation or a
license to accept a partial energy/gradient. Scientific settings are retained.
"""
from copy import deepcopy
import json
import re


def evidence(task, directory, outcome):
    if (task.get('kind') != 'ase_optimize'
            or task.get('profile', {}).get('calculator') != 'molpro'
            or task.get('resources', {}).get('cores', 1) <= 1
            or outcome.get('status') != 'failed'
            or outcome.get('failure_phase') != 'execution'):
        return None
    match = re.fullmatch(r'Molpro force evaluation (\w+) exited with status (\d+); '
                         r'diagnostics: (\w+\.failure\.json)\.', outcome.get('error', ''))
    if not match or match[3] != match[1]+'.failure.json':
        return None
    path = directory/match[3]
    try:
        record = json.loads(path.read_text())
    except (OSError, ValueError):
        return None
    if record.get('returncode') != int(match[2]) or record.get('explicit_io_failure'):
        return None
    text = '\n'.join(v.get('error_context', '') + '\n' + v.get('tail', '')
                     for v in record.get('files', {}).values())
    # A generic MPI_Abort/SIGKILL can be a consequence of a scientific error,
    # timeout, or OOM. Require the actual communication error and exclude the
    # recognizable alternate causes rather than retrying every exit 255.
    if (not re.search(r'Received an Error in Communication', text, re.I)
            or not re.search(r'MPI_Abort', text, re.I)
            or re.search(r'out.of.memory|oom.kill|time limit|insufficient memory|'
                         r'not enough memory|no space left|no convergence|'
                         r'catastrophic failure in diagonalization|expansion set has become singular|'
                         r'failed to converge|\?\s*(?:ERROR|INSUFFICIENT)', text, re.I)):
        return None
    return {'file': path.name, 'host': record.get('host'), 'returncode': record['returncode']}


def recover(graph, spec, state):
    from kinbot.anl import dispatch, live_control
    from kinbot.anl.salvage import descendants
    from kinbot.anl.timeout_recovery import _job_id
    from kinbot.manual_fallback import task_options
    if spec.get('intent', {}).get('vrc_fragment_retry'):
        return spec, state
    replacement, updated = deepcopy(spec), deepcopy(state)
    reset, reasons = set(), {}
    for index, task in enumerate(replacement['tasks']):
        ident = task['id']
        entry = state['tasks'].get(ident, {})
        if (entry.get('status') not in ('failed', 'unavailable')
                or state.get('molpro_serial_retries', {}).get(ident)):
            continue
        directory = graph/'tasks'/ident
        path = directory/'execution.json'
        if not path.is_file():
            continue
        outcome = json.loads(path.read_text())
        proof = evidence(task, directory, outcome)
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
        patched = task_options(task, {'schema': 1, 'resources': {'cores': 1}})
        replacement['tasks'][index] = patched
        reset.add(ident)
        reset.update(n for n in children if state['tasks'].get(n, {}).get('status')
                     in ('blocked', 'skipped'))
        updated.setdefault('molpro_serial_retries', {})[ident] = dict(proof, previous_job_id=job)
        reasons[ident] = 'One serial retry of recorded Molpro communication abort; scientific model unchanged.'
    if reset:
        live_control.commit(graph, replacement, updated, reset, reasons)
        _, spec, state = dispatch._load(graph)
    return spec, state
