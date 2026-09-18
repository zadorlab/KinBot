"""Read method-aware components from verified dispatcher task records."""

from __future__ import annotations

import json

from kinbot.anl.dispatch import _load, _verify_execution, _verify_stage_files
from kinbot.anl.model import ComponentResult, IncompleteRecipeError
from kinbot.anl.results import parse_result


def task_component(run_dir, task_id, *, key, state_id):
    """Return one validated native QC component from a completed task.

    The source is reparsed after checking all staged and output artifact hashes.
    Derived CBS, core-valence, relativistic, and higher-order providers are
    still separate work; a dispatch success alone cannot create those terms.
    """
    run_dir, spec, state = _load(run_dir)
    tasks = {task['id']: task for task in spec['tasks']}
    if task_id not in tasks:
        raise KeyError(task_id)
    task = tasks[task_id]
    entry = state['tasks'].get(task_id)
    if entry is None or entry['status'] != 'complete':
        raise IncompleteRecipeError(f'{task_id}: task is not complete.')
    request = task.get('result_parser')
    if request is None:
        raise IncompleteRecipeError(f'{task_id}: no method-aware parser was declared.')
    _verify_stage_files(run_dir, task, entry)
    directory = run_dir / 'tasks' / task_id
    execution = json.loads((directory / 'execution.json').read_text())
    _verify_execution(run_dir, task, entry, execution)
    if execution['status'] != 'executed':
        raise IncompleteRecipeError(f'{task_id}: QC execution failed.')
    native = directory / request['file']
    parsed = parse_result(native.read_text(errors='replace'), request)
    if parsed != execution.get('details', {}).get('parsed_result'):
        raise ValueError(f'{task_id}: saved parsed result differs from native output.')
    kind = parsed['kind']
    settings = {}
    review_required = parsed.get('review_required', False)
    if kind == 'molpro_energy':
        quantity = 'electronic'
        value = parsed['energy_hartree']
        method = parsed['method']
        basis = parsed['basis']
        if method == 'CCSD(T)-F12b':
            settings['scale_trip'] = 1
    elif kind == 'molpro_harmonic':
        quantity = 'zpe'
        value = parsed['zpe']['hartree']
        method = parsed['method']
        basis = parsed['basis']
    elif kind == 'gaussian_vpt2':
        if parsed['optimized_in_job']:
            raise IncompleteRecipeError(
                f'{task_id}: VPT2 optimized in the frequency job; '
                'use a separate accepted L2 geometry.')
        quantity = 'correction'
        value = parsed['anharmonic_correction_hartree']
        method = (parsed['method'] + '-D3BJ'
                  if parsed['dispersion'] == 'GD3BJ' else parsed['method'])
        basis = parsed['basis']
        settings['dispersion'] = parsed['dispersion']
    elif kind == 'cfour_dboc':
        if 'basis' not in request:
            raise IncompleteRecipeError(
                f'{task_id}: CFOUR DBOC task lacks a basis-specific parser.')
        quantity = 'correction'
        value = parsed['selected']['hartree']
        method = parsed['selected_level']
        basis = request['basis']
    else:
        raise ValueError(f'{task_id}: unsupported component parser {kind!r}.')
    return ComponentResult(
        key=key, value_hartree=value, quantity=quantity, method=method,
        basis=basis, backend=task['backend'].lower(), state_id=state_id,
        charge=spec['molecule'].get('charge', 0),
        multiplicity=spec['molecule'].get('multiplicity', 1),
        geometry_sha256=execution['geometry_sha256'],
        source_sha256=execution['artifacts'][request['file']],
        source=str(native), settings=settings, review_required=review_required,
    )
