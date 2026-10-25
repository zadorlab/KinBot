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
    if VPT2 in by_id:
        return result
    from kinbot.anl.vpt2_refinement import _tasks
    from kinbot.manual_fallback import logical_task
    original = deepcopy(logical_task(source))
    parser = original['result_parser']
    if (parser.get('method') != 'B2PLYP' or parser.get('basis') != 'cc-pVTZ'
            or parser.get('dispersion') != 'GD3BJ'):
        return result  # A custom lower-level model is not this declared pair.
    # Use the declared L2 method, not an incidental recovery or resource edit.
    original.pop('manual_fallback', None)
    original['id'] = VPT2
    original['geometry_from'] = 'initial'
    original.pop('depends_on', None)
    original['result_parser']['basis'] = 'cc-pVQZ'
    original['input_template'] = original['input_template'].replace('cc-pVTZ', 'cc-pVQZ')
    geometry, harmonic, vpt2 = _tasks(result, original)
    rename = {geometry['id']: GEOMETRY, harmonic['id']: HARMONIC}
    for task in (geometry, harmonic, vpt2):
        task['id'] = rename.get(task['id'], task['id'])
        task['geometry_from'] = rename.get(task['geometry_from'], task['geometry_from'])
        if 'depends_on' in task:
            task['depends_on'] = [rename.get(n, n) for n in task['depends_on']]
        task['anharmonic_basis_extension'] = True
    # This is the requested basis extension, not an exhausted-recovery attempt.
    vpt2.pop('vpt2_geometry_refinement', None)
    result['tasks'].extend((geometry, harmonic, vpt2))
    intent['anharmonic_basis_extension'] = dict(schema=1, owner='L3',
        method='B2PLYP', dispersion=source['result_parser'].get('dispersion', ''),
        bases=['cc-pVTZ', 'cc-pVQZ'], lower_task='gaussian_vpt2', upper_task=VPT2,
        geometry_task=GEOMETRY, harmonic_task=HARMONIC,
        upper_cardinal=4, extrapolation_power=ORIGINAL_ANL_POWER,
        hindered_rotors='accepted L2 only')
    return result


def available(spec, state):
    from kinbot.anl.protocol import omitted
    config = spec.get('intent', {}).get('anharmonic_basis_extension')
    return bool(config and all(config[key] not in omitted(spec) and state['tasks'].get(config[key], {}).get('status') == 'complete'
                               for key in ('lower_task', 'upper_task', 'geometry_task', 'harmonic_task')))


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
                                             vpt2_cbs=True).requirements
                           if r.key == 'vpt2_correction')
    requirement = replace(requirement, settings={**requirement.settings,
        'extrapolation_power': config['extrapolation_power']})
    return cbs_task_component(run, config['lower_task'], config['upper_task'],
        requirement=requirement, state_id=state_id,
        lower_basis=config['bases'][0], upper_basis=config['bases'][1])
