"""Read method-aware components from verified dispatcher task records."""

from __future__ import annotations

import json
import hashlib
import math
import re

from kinbot.anl.dispatch import _load, _verify_execution, _verify_stage_files
from kinbot.anl.extrapolation import two_point_cbs
from kinbot.anl.model import (ComponentRequirement, ComponentResult,
                              IncompleteRecipeError)
from kinbot.anl.results import parse_result


def task_component(run_dir, task_id, *, key, state_id):
    """Return one validated native QC component from a completed task.

    The source is reparsed after checking all staged and output artifact hashes.
    Derived core-valence, relativistic, and higher-order providers are
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


def cbs_task_component(run_dir, lower_task_id, upper_task_id, *,
                       requirement: ComponentRequirement, state_id: str,
                       lower_basis: str, upper_basis: str) -> ComponentResult:
    """Extrapolate two independently verified Molpro task results.

    Both native outputs are reparsed and checked against their execution
    records by ``task_component``. The result records a digest of the two
    source hashes and the exact extrapolation parameters.
    """
    if lower_task_id == upper_task_id:
        raise ValueError('A CBS pair needs two distinct tasks.')
    if not isinstance(requirement, ComponentRequirement):
        raise TypeError('CBS requirement must be a ComponentRequirement.')
    lower = task_component(run_dir, lower_task_id, key=lower_task_id,
                           state_id=state_id)
    upper = task_component(run_dir, upper_task_id, key=upper_task_id,
                           state_id=state_id)
    return _cbs_components(lower, upper, requirement=requirement,
                           lower_basis=lower_basis, upper_basis=upper_basis)


def _cbs_components(lower, upper, *, requirement, lower_basis, upper_basis):
    basis_pair = f'CBS({lower_basis},{upper_basis})'
    if (not lower_basis or not upper_basis or lower_basis == upper_basis
            or not (requirement.basis == basis_pair
                    or requirement.basis.startswith(basis_pair + '//'))):
        raise ValueError('CBS basis pair disagrees with the recipe.')
    if requirement.quantity not in ('electronic', 'zpe') \
            or 'composite' not in requirement.backends:
        raise ValueError('Recipe requirement is not a composite CBS term.')
    if (lower.basis, upper.basis) != (lower_basis, upper_basis):
        raise ValueError('CBS native basis order differs from the recipe.')
    if (lower.method != requirement.method or upper.method != requirement.method
            or lower.quantity != requirement.quantity
            or upper.quantity != requirement.quantity):
        raise ValueError('CBS method or quantity differs between inputs.')
    if lower.backend != 'molpro' or upper.backend != 'molpro':
        raise ValueError('CBS inputs must come from Molpro native parsers.')
    if ((lower.state_id, lower.charge, lower.multiplicity,
         lower.geometry_sha256) !=
            (upper.state_id, upper.charge, upper.multiplicity,
             upper.geometry_sha256)):
        raise ValueError('CBS inputs have different state or geometry.')
    if lower.geometry_sha256 is None:
        raise ValueError('CBS inputs lack a geometry hash.')
    if lower.review_required or upper.review_required:
        raise IncompleteRecipeError('CBS input requires native quality review.')
    if lower.settings != upper.settings:
        raise ValueError('CBS input calculation settings differ.')
    settings = dict(requirement.settings)
    if 'upper_cardinal' not in settings or 'extrapolation_power' not in settings:
        raise ValueError('CBS recipe lacks cardinal number or exponent.')
    for key, value in settings.items():
        if key not in ('extrapolation_power', 'upper_cardinal') \
                and lower.settings.get(key) != value:
            raise ValueError(f'CBS input {key} setting differs from the recipe.')
    for item in (lower, upper):
        if (not item.source or not isinstance(item.source_sha256, str)
                or re.fullmatch(r'[0-9a-f]{64}', item.source_sha256) is None
                or not math.isfinite(item.value_hartree)):
            raise ValueError('CBS input lacks valid native provenance.')
    cardinal = settings['upper_cardinal']
    power = settings['extrapolation_power']
    value = two_point_cbs(lower.value_hartree, upper.value_hartree,
                          upper_cardinal=cardinal, power=power)
    provenance = {
        'formula': 'E_n + alpha_n*(E_n-E_(n-1))',
        'lower_sha256': lower.source_sha256,
        'upper_sha256': upper.source_sha256,
        'lower_basis': lower_basis, 'upper_basis': upper_basis,
        'upper_cardinal': cardinal, 'power': power,
    }
    digest = hashlib.sha256(json.dumps(provenance, sort_keys=True,
                                       separators=(',', ':')).encode()).hexdigest()
    return ComponentResult(
        key=requirement.key, value_hartree=value,
        quantity=requirement.quantity, method=requirement.method,
        basis=requirement.basis, backend='composite',
        state_id=lower.state_id, charge=lower.charge,
        multiplicity=lower.multiplicity,
        geometry_sha256=lower.geometry_sha256,
        source_sha256=digest,
        source=f'CBS({lower.source},{upper.source})', settings=settings,
    )
