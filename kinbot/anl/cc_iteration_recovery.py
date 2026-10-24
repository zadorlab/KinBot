"""Bounded native CC iteration retries; no convergence test is relaxed.

Molpro CC/CC-F12 supports MAXIT; CFOUR VCC supports CC_MAXCYC. MRCC's
ccprog=mrcc higher-order engine has a compiled iteration ceiling, so applying
ccmaxit to that engine would be an ineffective retry and is not attempted.
"""
from copy import deepcopy
import json
import re

from kinbot.execution import file_sha256


def _tail(path):
    try:
        with path.open('rb') as stream:
            stream.seek(max(0, path.stat().st_size-2_000_000))
            return stream.read().decode(errors='replace')
    except OSError:
        return ''


def options(task, directory, outcome):
    if outcome.get('status') != 'failed' or outcome.get('failure_phase') != 'execution':
        return None
    ase = task['kind'] == 'ase_optimize'
    backend = task.get('profile', {}).get('calculator') if ase else task.get('backend')
    if backend not in ('molpro', 'cfour'):
        return None
    if ase:
        profile = task['profile']
        if profile['method'].upper() not in ('CCSD(T)', 'CCSD(T)-F12B'):
            return None
        match = re.fullmatch(r'Molpro force evaluation (\w+) exited with status ([1-9]\d*); '
                            r'diagnostics: (\w+\.failure\.json)\.', outcome.get('error', ''))
        if not match or match[3] != match[1]+'.failure.json':
            return None
        try:
            record = json.loads((directory/match[3]).read_text())
        except (OSError, ValueError):
            return None
        if record.get('explicit_io_failure') or record.get('returncode') != int(match[2]):
            return None
        files = [directory/(match[1]+suffix) for suffix in ('.out', '.log')]
        old = profile.get('calculator_kwargs', {}).get('cc_max_iterations', 0)
        template = ''
    else:
        parser = task.get('result_parser', {})
        if not str(parser.get('method', '')).upper().startswith(('CC', 'UCC', 'RCC')):
            return None
        native = directory/parser.get('file', task.get('stdout', 'stdout.txt'))
        files = [native]
        if backend == 'molpro':
            files.append(native.with_suffix('.log'))
        template = task['input_template']
        if backend == 'molpro':
            blocks = list(re.finditer(r'\{\s*[ru]?ccsd(?:\(t\))?(?:-f12[ab])?(?=[,;\s}])[^{}]*\}', template, re.I))
            if len(blocks) != 1:
                return None  # Ambiguous composite deck; do not edit another solve.
            old_block = blocks[0][0]
            limits = re.findall(r'\bmaxit\s*(?:=|,)\s*(\d+)', old_block, re.I)
        else:
            if not re.search(r'(?im)^CC_PROG=VCC\s*$', template):
                return None
            limits = re.findall(r'\bCC_MAXCYC\s*=\s*(\d+)', template, re.I)
        if len(limits) > 1:
            return None
        old = int(limits[0]) if limits else 0
    text = '\n'.join(_tail(p) for p in files)
    if re.search(r'out.of.memory|no space left|time limit|insufficient memory', text, re.I):
        return None
    if backend == 'molpro':
        if not (re.search(r'CONVERGENCE NOT REACHED AFTER MAX\. ITERATIONS', text, re.I)
                and re.search(r'The problem occurs in ccsd', text, re.I)):
            return None
        limit = next((v for v in (100, 200) if v > old), None)
    else:
        # Require a CC-specific error, not an SCF/geometry iteration message.
        if not re.search(r'(?:CC(?:SDT?Q?)? (?:iterations|equations|amplitudes)[^\n]*'
                         r'(?:not converg|failed to converg|maximum)|'
                         r'coupled.cluster[^\n]*(?:not converg|failed to converg))', text, re.I):
            return None
        limit = next((v for v in (600, 1000) if v > old), None)
    if limit is None:
        return None
    if ase:
        change = {'calculator_kwargs': {'cc_max_iterations': limit}}
    elif backend == 'molpro':
        new = (re.sub(r'\bmaxit\s*(?:=|,)\s*\d+', f'maxit={limit}', old_block, flags=re.I)
               if limits else old_block[:-1]+f';maxit,{limit}'+'}')
        change = {'input_edits': [{'old': old_block, 'new': new}]}
    else:
        old_text = re.search(r'\bCC_MAXCYC\s*=\s*\d+', template, re.I)
        if old_text is None:
            return None
        change = {'input_edits': [{'old': old_text[0], 'new': f'CC_MAXCYC={limit}'}]}
    proof = dict(backend=backend, previous_limit=old or 'native_default', limit=limit,
                 native={p.name:file_sha256(p) for p in files if p.is_file()})
    return dict(schema=1, **change), proof


def recover(graph, spec, state):
    from kinbot.anl import dispatch, live_control
    from kinbot.anl.salvage import descendants
    from kinbot.anl.timeout_recovery import _job_id
    from kinbot.manual_fallback import task_options
    if spec.get('intent', {}).get('vrc_fragment_retry'):
        return spec, state
    replacement, updated = deepcopy(spec), deepcopy(state)
    reset, reasons = set(), {}
    for index, task in enumerate(spec['tasks']):
        ident = task['id']; entry = state['tasks'].get(ident, {})
        history = state.get('cc_iteration_retries', {}).get(ident, [])
        directory = graph/'tasks'/ident
        if (entry.get('status') not in ('failed', 'unavailable') or len(history) >= 2
                or not (directory/'execution.json').is_file()):
            continue
        outcome = json.loads((directory/'execution.json').read_text())
        planned = options(task, directory, outcome)
        if planned is None:
            continue
        change, proof = planned
        if any(h['limit'] >= proof['limit'] for h in history):
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
        replacement['tasks'][index] = task_options(task, change)
        reset.add(ident)
        reset.update(n for n in children if state['tasks'].get(n, {}).get('status') in ('blocked', 'skipped'))
        updated.setdefault('cc_iteration_retries', {})[ident] = [*history, dict(proof, previous_job_id=job)]
        reasons[ident] = f"Native CC nonconvergence: increase iteration limit to {proof['limit']}; thresholds and theory unchanged."
    if reset:
        live_control.commit(graph, replacement, updated, reset, reasons)
        _, spec, state = dispatch._load(graph)
    return spec, state
