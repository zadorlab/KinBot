"""Shared, provenance-checked ROTD to MESS staging for direct and PES runs."""
from pathlib import Path
import hashlib
import json
import logging
import re

from kinbot import constants
from kinbot.execution import atomic_json, file_sha256


logger = logging.getLogger('KinBot')
_MARKER = re.compile(r'! kinbot_rotdpy_begin ([0-9a-f]{64})\n(.*?)'
                     r'! kinbot_rotdpy_end \1\n', re.S)


def selected_reaction(parameters, name):
    """None retains legacy global required policy when no selectors exist."""
    from kinbot.species_routing import matches_name
    selections = [choices for key in ('vrc_tst_scan', 'vrc_tst_noscan')
                  for choices in parameters.get(key, {}).values()]
    return (any(matches_name(name, choices) for choices in selections)
            if selections else None)


def locate_result(parameters, name):
    mode = str(parameters.get('rotdpy_mess_mode', 'auto')).casefold()
    if mode == 'off' or selected_reaction(parameters, name) is False:
        return None
    directory = Path(parameters.get('rotdpy_directory', 'rotdPy'))
    inp = directory / f'{name}.py'
    multistate = directory / f'{name}.multistate.json'
    if parameters.get('vrc_tst_surface_model_file') and not multistate.is_file():
        raise FileNotFoundError(f'{name}: state-resolved ROTD results are required.')
    if not multistate.is_file() and not (inp.is_file() and inp.with_suffix('.rotdpy.json').is_file()):
        if mode in ('required', 'production'):
            raise FileNotFoundError(f'{name}: completed ROTD_py input and manifest are required under {directory}.')
        return None
    return inp, multistate if multistate.is_file() else None


def build_model(parameters, context, ground_energy, destination='me'):
    """Stage flux once and convolve the worker's own conserved modes/weights."""
    name = context['reaction']
    located = locate_result(parameters, name)
    if located is None:
        return None, None
    inp, multistate = located
    mode = str(parameters.get('rotdpy_mess_mode', 'auto')).casefold()
    destination = Path(destination)
    from kinbot.rotdpy import (number_of_states_file, number_of_states_ground,
                               stage_number_of_states_file)
    try:
        if multistate:
            from kinbot.vrc_multistate import stage_channel
            provenance = stage_channel(multistate, destination,
                required_profile='production' if mode == 'production' else None)
        else:
            source, provenance = number_of_states_file(inp,
                parameters.get('rotdpy_mess_energy_index', -1),
                required_profile='production' if mode == 'production' else None)
    except (OSError, RuntimeError, ValueError) as error:
        if multistate or mode in ('required', 'production'):
            raise RuntimeError(f'{name}: cannot use ROTD_py result: {error}') from error
        logger.warning('%s: invalid ROTD_py result; using phase-space fallback: %s', name, error)
        return None, None
    provenance = dict(provenance)
    zero = round(float(ground_energy) - context['reference_shift'], 2)
    if not multistate:
        safe_name = re.sub(r'[^A-Za-z0-9_.-]+', '_', name)
        target = destination / f'rotd_{safe_name}_Ne_{provenance["energy_index"]}.out'
        preprocessing = stage_number_of_states_file(source, target)
        provenance.update(mess_file=target.name,
            mess_file_sha256=preprocessing.pop('mess_file_sha256'),
            mess_preprocessing=preprocessing)
        threshold = number_of_states_ground(target)
        zero += threshold * constants.AUtoKCAL / constants.AUtoCM
        provenance['transitional_ground_cm_inverse'] = threshold
    provenance.update({key: context[key] for key in (
        'rotational_symmetry_product', 'indistinguishability_factor',
        'uncertainty_factor', 'mess_symmetry_factor', 'fragment_multiplicities',
        'reactive_surface_multiplicity')})
    if provenance.get('dynamical_correction') is not None:
        provenance['raw_flux_equivalent_mess_symmetry_factor'] = (
            context['mess_symmetry_factor'] / provenance['dynamical_correction'])
    common = {key: context[key] for key in ('stoich', 'nfreq', 'freq', 'hinderedrotor')}
    common['symmetry'] = context['mess_symmetry_factor']
    if multistate:
        from kinbot.vrc_multistate import render_mess
        zero += provenance['baseline_cm_inverse'] * constants.AUtoKCAL / constants.AUtoCM
        def rrho(stage):
            return format_rotd_rrho(context['template'], **common, file=stage['file'],
                nelec=1, electronic_levels='            0 1',
                zeroenergy=zero + stage['ground_cm_inverse'] * constants.AUtoKCAL / constants.AUtoCM,
                provenance_sha256=stage['source_sha256'])
        return render_mess(provenance, rrho), provenance
    return format_rotd_rrho(context['template'], **common, **context['electronic'],
        file=target.name, zeroenergy=zero,
        provenance_sha256=provenance['source_sha256']), provenance


def format_rotd_rrho(template, **fields):
    """MESS requires omission, not a zero count, for vibration-free fragments."""
    text = template.format(**fields)
    if fields['nfreq'] == 0:
        text = re.sub(r'(?m)^\s*Frequencies\[1/cm\]\s+0\s*$', '', text)
    return text


def defer_model(context, fallback, directory='me/rotdpy_models'):
    """Bind each worker/UQ convolution model without its unfinished energy axis."""
    payload = (json.dumps(context, sort_keys=True, indent=2) + '\n').encode()
    digest = hashlib.sha256(payload).hexdigest()
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    (root / f'{digest}.json').write_bytes(payload)
    return f'! kinbot_rotdpy_begin {digest}\n{fallback}! kinbot_rotdpy_end {digest}\n'


def resolve_models(text, worker, parameters, ground_energy, sources, reaction):
    """Resolve only after global ROTD execution and final network/UQ energies."""
    matches = list(_MARKER.finditer(text))
    required = (parameters.get('rotdpy_mess_mode') in ('required', 'production')
                or bool(parameters.get('vrc_tst_surface_model_file')))
    if (not matches and required and parameters.get('rotdpy_mess_mode') != 'off'
            and selected_reaction(parameters, reaction) is not False):
        raise ValueError(f'{reaction}: missing deferred ROTD context; regenerate the PES worker MESS artifacts.')
    def replace(match):
        if (parameters.get('rotdpy_mess_mode') == 'off'
                or selected_reaction(parameters, reaction) is False):
            return match[2]
        path = Path(worker) / 'me' / 'rotdpy_models' / f'{match[1]}.json'
        if file_sha256(path) != match[1]:
            raise ValueError('PES ROTD convolution context changed after worker export.')
        context = json.loads(path.read_text())
        if context.get('schema') != 1 or context.get('reaction') != reaction:
            raise ValueError('PES ROTD context does not match the exported reaction.')
        model, provenance = build_model(parameters, context, ground_energy)
        if model is None:
            return match[2]
        sources[reaction] = {**provenance, 'worker_context': str(path.resolve()),
                             'worker_context_sha256': match[1],
                             'network_ground_energy_kcal_mol': float(ground_energy)}
        return model
    return _MARKER.sub(replace, text)


def write_sources(destination, parameters, sources):
    # An empty record replaces stale provenance when a rewrite disables ROTD
    # or returns to the explicitly permitted phase-space fallback.
    atomic_json(Path(destination), {'schema': 1, 'selection_policy': {
        'mode': parameters.get('rotdpy_mess_mode', 'auto'),
        'energy_index': parameters.get('rotdpy_mess_energy_index', -1)},
        'reactions': sources})
