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


def _verified_task_result(run_dir, task_id):
    """Reparse one complete task after verifying every recorded artifact."""
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
    return run_dir, spec, task, execution, native, parsed


def task_component(run_dir, task_id, *, key, state_id):
    """Return one validated native QC component from a completed task.

    The source is reparsed after checking all staged and output artifact hashes.
    Derived core-valence, relativistic, and higher-order providers are
    still separate work; a dispatch success alone cannot create those terms.
    """
    run_dir, spec, task, execution, native, parsed = _verified_task_result(
        run_dir, task_id)
    request = task['result_parser']
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


def attach_task_vpt2_frequencies(species, run_dir, task_id, **match_options):
    """Attach verified mode-resolved VPT2 corrections to a KinBot species."""
    from kinbot.energy import attach_anharmonic_frequencies

    _, _, task, execution, native, parsed = _verified_task_result(run_dir, task_id)
    if parsed.get('kind') != 'gaussian_vpt2' or task['backend'].lower() != 'gaussian':
        raise ValueError(f'{task_id}: expected a Gaussian VPT2 task.')
    frequencies = attach_anharmonic_frequencies(
        species, parsed, **match_options)
    species.anl_thermochemistry_frequency_source.update({
        'task_id': task_id,
        'native_output': str(native),
        'source_sha256': execution['artifacts'][task['result_parser']['file']],
    })
    return frequencies


def cbs_task_component(run_dir, lower_task_id, upper_task_id, *,
                       requirement: ComponentRequirement, state_id: str,
                       lower_basis: str, upper_basis: str) -> ComponentResult:
    """Extrapolate two verified electronic or zero-point task results.

    Both native outputs are reparsed and checked against their execution
    records by ``task_component``. Electronic single points must share one
    geometry. Harmonic and VPT2 pairs must each trace to a completed,
    basis-matched optimization, with the higher-basis geometry retained as
    the composite result's geometry. The source digest records both inputs,
    geometry sources, and exact extrapolation parameters.
    """
    if lower_task_id == upper_task_id:
        raise ValueError('A CBS pair needs two distinct tasks.')
    if not isinstance(requirement, ComponentRequirement):
        raise TypeError('CBS requirement must be a ComponentRequirement.')
    lower = task_component(run_dir, lower_task_id, key=lower_task_id,
                           state_id=state_id)
    upper = task_component(run_dir, upper_task_id, key=upper_task_id,
                           state_id=state_id)
    run_dir, spec, state = _load(run_dir)
    tasks = {task['id']: task for task in spec['tasks']}
    if requirement.quantity == 'electronic':
        settings = requirement.settings
        if (not settings.get('geometry_method')
                or not settings.get('geometry_basis')):
            raise ValueError('Electronic CBS recipe lacks a geometry level.')
        if (tasks[lower_task_id].get('geometry_from') !=
                tasks[upper_task_id].get('geometry_from')):
            raise ValueError('Electronic CBS inputs have different geometry sources.')
        geometry_sources = tuple(
            _verified_optimized_geometry(
                run_dir, tasks, state, tasks[task_id], component,
                expected_method=settings['geometry_method'],
                expected_basis=settings['geometry_basis'])
            for task_id, component in ((lower_task_id, lower),
                                       (upper_task_id, upper)))
    else:
        geometry_sources = (
            _verified_optimized_geometry(run_dir, tasks, state,
                                         tasks[lower_task_id], lower),
            _verified_optimized_geometry(run_dir, tasks, state,
                                         tasks[upper_task_id], upper),
        )
    return _cbs_components(lower, upper, requirement=requirement,
                           lower_basis=lower_basis, upper_basis=upper_basis,
                           geometry_sources=geometry_sources)


def _verified_optimized_geometry(run_dir, tasks, state, task, component, *,
                                 expected_method=None, expected_basis=None):
    source_id = task.get('geometry_from', 'initial')
    source_task = tasks.get(source_id)
    if source_task is None or source_task.get('kind') != 'ase_optimize':
        raise ValueError(f"{task['id']}: CBS input needs an optimized geometry.")
    profile = source_task.get('profile', {})
    backend = ('gaussian' if component.quantity == 'correction' else 'molpro')
    method = expected_method or (
        component.method.removesuffix('-D3BJ')
        if component.quantity == 'correction' else component.method)
    basis = expected_basis or component.basis
    if (not isinstance(profile, dict)
            or profile.get('calculator', '').lower() not in (
                ('gaussian', 'gauss') if backend == 'gaussian' else ('molpro',))
            or profile.get('method', '').casefold() != method.casefold()
            or profile.get('basis', '').casefold() != basis.casefold()):
        raise ValueError(f"{task['id']}: optimized geometry level differs from CBS input.")
    if component.quantity == 'correction':
        keywords = profile.get('calculator_kwargs', {})
        if (not isinstance(keywords, dict) or
                keywords.get('EmpiricalDispersion', '').upper()
                != component.settings.get('dispersion', '')):
            raise ValueError(f"{task['id']}: optimized geometry dispersion differs.")
    source_entry = state['tasks'].get(source_id)
    if source_entry is None or source_entry.get('status') != 'complete':
        raise IncompleteRecipeError(f'{source_id}: optimized geometry is not complete.')
    _verify_stage_files(run_dir, source_task, source_entry)
    execution = json.loads((run_dir / 'tasks' / source_id /
                            'execution.json').read_text())
    geometry_hash = _verify_execution(run_dir, source_task, source_entry,
                                      execution)
    if (execution['status'] != 'executed' or geometry_hash is None
            or geometry_hash != source_entry.get('final_geometry_sha256')
            or geometry_hash != component.geometry_sha256):
        raise ValueError(f"{task['id']}: accepted optimized geometry differs.")
    return {'task_id': source_id, 'geometry_sha256': geometry_hash,
            'final_xyz_sha256': execution['artifacts']['final.xyz']}


def _cbs_components(lower, upper, *, requirement, lower_basis, upper_basis,
                    geometry_sources=None):
    allowed = (requirement.quantity == 'electronic'
               or requirement.quantity == 'zpe'
               or (requirement.quantity == 'correction'
                   and requirement.key == 'vpt2_correction'))
    if not allowed or 'composite' not in requirement.backends:
        raise ValueError('Recipe requirement is not a supported CBS term.')
    basis_pair = f'CBS({lower_basis},{upper_basis})'
    if (not lower_basis or not upper_basis or lower_basis == upper_basis
            or not (requirement.basis == basis_pair
                    or requirement.basis.startswith(basis_pair + '//'))):
        raise ValueError('CBS basis pair disagrees with the recipe.')
    if (lower.basis, upper.basis) != (lower_basis, upper_basis):
        raise ValueError('CBS native basis order differs from the recipe.')
    if (lower.method != requirement.method or upper.method != requirement.method
            or lower.quantity != requirement.quantity
            or upper.quantity != requirement.quantity):
        raise ValueError('CBS method or quantity differs between inputs.')
    backend = 'gaussian' if requirement.quantity == 'correction' else 'molpro'
    if lower.backend != backend or upper.backend != backend:
        raise ValueError(f'CBS inputs must come from {backend} native parsers.')
    if ((lower.state_id, lower.charge, lower.multiplicity) !=
            (upper.state_id, upper.charge, upper.multiplicity)):
        raise ValueError('CBS inputs have different electronic states.')
    if lower.geometry_sha256 is None or upper.geometry_sha256 is None:
        raise ValueError('CBS inputs lack geometry hashes.')
    if requirement.quantity == 'electronic':
        if lower.geometry_sha256 != upper.geometry_sha256:
            raise ValueError('Electronic CBS inputs have different geometries.')
        if (not isinstance(geometry_sources, tuple)
                or len(geometry_sources) != 2
                or geometry_sources[0]['task_id'] != geometry_sources[1]['task_id']
                or any(source['geometry_sha256'] != lower.geometry_sha256
                       for source in geometry_sources)):
            raise ValueError('Electronic CBS lacks one verified highest-level geometry.')
    else:
        if (requirement.settings.get('geometry_mode') != 'basis_optimized'
                or not isinstance(geometry_sources, tuple)
                or len(geometry_sources) != 2
                or geometry_sources[0]['geometry_sha256'] != lower.geometry_sha256
                or geometry_sources[1]['geometry_sha256'] != upper.geometry_sha256
                or geometry_sources[0]['task_id'] == geometry_sources[1]['task_id']):
            raise ValueError('CBS zero-point inputs lack basis-matched optimized geometries.')
    if lower.review_required or upper.review_required:
        raise IncompleteRecipeError('CBS input requires native quality review.')
    if lower.settings != upper.settings:
        raise ValueError('CBS input calculation settings differ.')
    settings = dict(requirement.settings)
    if 'upper_cardinal' not in settings or 'extrapolation_power' not in settings:
        raise ValueError('CBS recipe lacks cardinal number or exponent.')
    for key, value in settings.items():
        if key not in ('extrapolation_power', 'upper_cardinal', 'geometry_mode',
                       'geometry_method', 'geometry_basis') \
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
        'geometry_sources': geometry_sources,
    }
    digest = hashlib.sha256(json.dumps(provenance, sort_keys=True,
                                       separators=(',', ':')).encode()).hexdigest()
    return ComponentResult(
        key=requirement.key, value_hartree=value,
        quantity=requirement.quantity, method=requirement.method,
        basis=requirement.basis, backend='composite',
        state_id=lower.state_id, charge=lower.charge,
        multiplicity=lower.multiplicity,
        geometry_sha256=upper.geometry_sha256,
        source_sha256=digest,
        source=f'CBS({lower.source},{upper.source})', settings=settings,
    )
