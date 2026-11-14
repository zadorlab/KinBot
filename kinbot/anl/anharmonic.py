"""L3-owned, basis-optimized TZ/QZ double-hybrid anharmonic corrections.

The accepted L2 point remains immutable. Only anharmonic *differences* are
extrapolated; neither a second harmonic ZPE nor rotor corrections are added.
"""
from copy import deepcopy
from dataclasses import replace

from kinbot.anl.extrapolation import ORIGINAL_ANL_POWER

GEOMETRY = 'l3_vpt2_qz_geometry'
HARMONIC = 'l3_vpt2_qz_harmonic'
VPT2 = 'gaussian_vpt2_qz'


def extend(spec):
    """Add the independent branch to a new or owner-migrated automatic graph."""
    result = deepcopy(spec)
    intent = result.get('intent', {})
    if (result.get('name') != 'anl-automatic-post-l2'
            or intent.get('highest_tier') not in ('ANL1', 'ANL1-F12')
            or intent.get('stationary_point') != 'minimum'
            or len(result['molecule']['symbols']) == 1):
        return result
    by_id = {t['id']: t for t in result['tasks']}
    source = by_id.get('gaussian_vpt2')
    if source is None:
        return result
    if intent.get('anharmonic_basis_extension') or VPT2 in by_id:
        return result
    from kinbot.anl.vpt2_refinement import _tasks
    from kinbot.manual_fallback import logical_task
    original = deepcopy(logical_task(source))
    parser = original['result_parser']
    from kinbot.anl.basis_reuse import canonical_basis, basis_pair, select
    main_basis = canonical_basis(parser.get('basis'))
    if (parser.get('method') != 'B2PLYP' or main_basis is None
            or parser.get('dispersion') != 'GD3BJ'):
        return result
    bases = basis_pair(main_basis)
    missing_basis = next(b for b in bases if b != main_basis)
    upper_missing = missing_basis == bases[1]
    ident = VPT2 if upper_missing else 'gaussian_vpt2_tz'
    geometry_id = GEOMETRY if upper_missing else 'l3_vpt2_tz_geometry'
    harmonic_id = HARMONIC if upper_missing else 'l3_vpt2_tz_harmonic'
    original.pop('manual_fallback', None)
    original['id'] = ident
    original['geometry_from'] = 'initial'
    original.pop('geometry_record', None)
    original.pop('depends_on', None)
    original['input_template'] = original['input_template'].replace(parser['basis'], missing_basis)
    original['result_parser']['basis'] = missing_basis
    point = select(result['molecule'], missing_basis)
    if point is not None:
        original['geometry_record'] = point
        additions = [original]
    else:
        geometry, harmonic, vpt2 = _tasks(result, original)
        rename = {geometry['id']: geometry_id, harmonic['id']: harmonic_id}
        additions = [geometry, harmonic, vpt2]
        for task in additions:
            task['id'] = rename.get(task['id'], task['id'])
            task['geometry_from'] = rename.get(task['geometry_from'], task['geometry_from'])
            if 'depends_on' in task:
                task['depends_on'] = [rename.get(n, n) for n in task['depends_on']]
        vpt2.pop('vpt2_geometry_refinement', None)
    for task in additions:
        task['anharmonic_basis_extension'] = True
    result['tasks'].extend(additions)
    config = dict(schema=2, owner='L3', method='B2PLYP', dispersion='GD3BJ',
        bases=bases, lower_task='gaussian_vpt2' if upper_missing else ident,
        upper_task=ident if upper_missing else 'gaussian_vpt2',
        extension_tasks=[t['id'] for t in additions],
        required_tasks=['gaussian_vpt2', *[t['id'] for t in additions]],
        upper_cardinal=4, extrapolation_power=ORIGINAL_ANL_POWER,
        hindered_rotors='accepted L2 only')
    if point is None:
        config.update(geometry_task=geometry_id, harmonic_task=harmonic_id)
    else:
        config['reused_geometry_sha256'] = point['record_sha256']
    intent['anharmonic_basis_extension'] = config
    return result


def available(spec, state):
    from kinbot.anl.protocol import omitted
    config = spec.get('intent', {}).get('anharmonic_basis_extension')
    if not config:
        return False
    required = config.get('required_tasks') or [config[k] for k in
        ('lower_task', 'upper_task', 'geometry_task', 'harmonic_task')]
    return all(n not in omitted(spec) and state['tasks'].get(n, {}).get('status') == 'complete'
               for n in required)


def component(run, *, state_id, requirement=None):
    """Use the existing native-verified, basis-optimized CBS evaluator."""
    from kinbot.anl.dispatch import _load
    from kinbot.anl.recipes import recipe
    from kinbot.anl.workflow import cbs_task_component
    _, spec, state = _load(run)
    if not available(spec, state):
        raise ValueError('TZ/QZ anharmonic branch is not complete.')
    config = spec['intent']['anharmonic_basis_extension']
    if requirement is None:
        requirement = next(r for r in recipe('ANL1', vpt2_method='B2PLYP-D3BJ',
                                             vpt2_cbs=True, vpt2_bases=config['bases']).requirements
                           if r.key == 'vpt2_correction')
    requirement = replace(requirement, settings={**requirement.settings,
        'extrapolation_power': config['extrapolation_power']})
    return cbs_task_component(run, config['lower_task'], config['upper_task'],
        requirement=requirement, state_id=state_id,
        lower_basis=config['bases'][0], upper_basis=config['bases'][1])
