"""VRC correction-potential models and their calculation provenance."""

from copy import deepcopy
from numbers import Real
import re

import numpy as np

from kinbot.anl.extrapolation import (ORIGINAL_ANL_POWER,
                                      extrapolation_coefficient,
                                      two_point_cbs)


_SPIN_SPLITTING_DEFAULTS = {
    'target_multiplicity': 1,
    'reference_multiplicity': 3,
    'fragment_multiplicities': None,
    'reference': {
        'method': 'uccsd(t)-f12b',
        'basis': ['avtz-f12', 'avqz-f12'],
        'extrapolation_power': ORIGINAL_ANL_POWER,
    },
    'gap_levels': [
        {
            'label': 'caspt2_f12',
            'method': 'caspt2-f12b(2,2)',
            'basis': ['avtz-f12', 'avqz-f12'],
            'extrapolation_power': ORIGINAL_ANL_POWER,
            'shift': 0.2,
            'ipea': 0.25,
        },
        {
            'label': 'mrci_f12_q',
            'method': 'mrci-f12b+q(2,2)',
            'basis': ['avtz-f12', 'avqz-f12'],
            'extrapolation_power': ORIGINAL_ANL_POWER,
        },
    ],
    # Ram et al. used an arithmetic mean.  This range-dependent refinement
    # gives the size-consistent CASPT2 gap more weight at long range and the
    # more complete MRCI gap more weight at short range.  The 2.5--5.0 A
    # interval spans the radical-radical transition region discussed there.
    'weighting': {
        'model': 'distance_switch',
        'distance_knots_angstrom': [2.5, 3.125, 3.75, 4.375, 5.0],
        'weights': {
            'caspt2_f12': [0.2, 0.4, 0.5, 0.6, 0.8],
            'mrci_f12_q': [0.8, 0.6, 0.5, 0.4, 0.2],
        },
    },
}


def correction_model(parameters):
    """Return the normalized correction model name."""
    model = str(parameters.get(
        'vrc_tst_correction_model', 'direct')).strip().casefold()
    aliases = {
        'direct': 'direct',
        'spin_splitting': 'spin_splitting',
        'klippenstein_spin_splitting': 'spin_splitting',
    }
    try:
        return aliases[model]
    except KeyError as exc:
        raise ValueError(
            'vrc_tst_correction_model must be direct or spin_splitting.') \
            from exc


def coupled_multiplicities(first, second):
    """Return all total-spin multiplicities allowed for two fragments."""
    if (isinstance(first, bool) or not isinstance(first, int) or first < 1
            or isinstance(second, bool) or not isinstance(second, int)
            or second < 1):
        raise ValueError('Fragment multiplicities must be positive integers.')
    twice_s_min = abs((first - 1) - (second - 1))
    twice_s_max = (first - 1) + (second - 1)
    return tuple(value + 1 for value in range(
        twice_s_min, twice_s_max + 1, 2))


def spin_splitting_settings(parameters):
    """Return validated spin-splitting settings with concise user overrides."""
    settings = deepcopy(_SPIN_SPLITTING_DEFAULTS)
    override = parameters.get('vrc_tst_spin_splitting')
    if override is None:
        override = {}
    if not isinstance(override, dict):
        raise ValueError('vrc_tst_spin_splitting must be an object.')
    override = deepcopy(override)
    legacy = {
        'singlet_multiplicity': 'target_multiplicity',
        'triplet_multiplicity': 'reference_multiplicity',
        'triplet_reference': 'reference',
    }
    for old, new in legacy.items():
        if old in override:
            if new in override:
                raise ValueError(
                    f'vrc_tst_spin_splitting cannot define both {old} and '
                    f'{new}.')
            override[new] = override.pop(old)
    allowed = set(settings)
    unknown = set(override) - allowed
    if unknown:
        raise ValueError('Unknown vrc_tst_spin_splitting option(s): '
                         + ', '.join(sorted(unknown)))
    for key, value in override.items():
        if key in ('reference', 'weighting'):
            if not isinstance(value, dict):
                raise ValueError(f'vrc_tst_spin_splitting.{key} must be an '
                                 'object.')
            settings[key].update(value)
        else:
            settings[key] = value
    if isinstance(settings.get('reference'), dict) \
            and isinstance(settings['reference'].get('basis'), str):
        settings['reference'].pop('extrapolation_power', None)
    for level in settings.get('gap_levels', []):
        if isinstance(level, dict) and isinstance(level.get('basis'), str):
            level.pop('extrapolation_power', None)
    if (isinstance(settings.get('weighting'), dict)
            and str(settings['weighting'].get('model', '')).casefold()
            == 'arithmetic_mean'):
        settings['weighting'] = {'model': 'arithmetic_mean'}
    _validate_spin_splitting_settings(settings)
    # Keep archived ethane manifests and callers readable while new records
    # use state-neutral names.  These aliases are views of the same validated
    # values and are never accepted alongside their canonical input names.
    settings['singlet_multiplicity'] = settings['target_multiplicity']
    settings['triplet_multiplicity'] = settings['reference_multiplicity']
    settings['triplet_reference'] = settings['reference']
    return settings


def _validate_spin_splitting_settings(settings):
    target = settings['target_multiplicity']
    reference_multiplicity = settings['reference_multiplicity']
    if (isinstance(target, bool) or not isinstance(target, int) or target < 1
            or isinstance(reference_multiplicity, bool)
            or not isinstance(reference_multiplicity, int)
            or reference_multiplicity < 1
            or target == reference_multiplicity):
        raise ValueError('Spin-splitting target and reference multiplicities '
                         'must be distinct positive integers.')
    fragments = settings['fragment_multiplicities']
    if fragments is not None:
        if (not isinstance(fragments, list) or len(fragments) != 2
                or any(isinstance(value, bool) or not isinstance(value, int)
                       or value < 1 for value in fragments)):
            raise ValueError('fragment_multiplicities must contain two '
                             'positive integers.')
        allowed = coupled_multiplicities(*fragments)
        if target not in allowed or reference_multiplicity not in allowed:
            raise ValueError('Spin-splitting multiplicities are not allowed '
                             'by the two fragment spins.')
    reference = settings['reference']
    if (not isinstance(reference, dict)
            or set(reference) - {'method', 'basis', 'extrapolation_power'}
            or not isinstance(reference.get('method'), str)
            or not reference['method'].strip()):
        raise ValueError('spin-splitting reference requires a method '
                         'and one basis or an adjacent basis pair.')
    basis = reference.get('basis')
    if isinstance(basis, str):
        if not basis.strip():
            raise ValueError('Reference basis cannot be empty.')
    elif isinstance(basis, list) and len(basis) == 2 \
            and all(isinstance(item, str) and item.strip() for item in basis):
        lower, upper = map(_basis_cardinal, basis)
        if upper != lower + 1:
            raise ValueError('Reference CBS bases must be adjacent.')
        power = reference.get('extrapolation_power', ORIGINAL_ANL_POWER)
        extrapolation_coefficient(upper, power=power)
    else:
        raise ValueError('Reference basis must be a string or a '
                         'two-item list of adjacent cardinal bases.')
    levels = settings['gap_levels']
    if not isinstance(levels, list) or not levels:
        raise ValueError('spin-splitting gap_levels must contain at least one '
                         'level.')
    labels = []
    for level in levels:
        required = {'label', 'method', 'basis'}
        if (not isinstance(level, dict) or not required <= set(level)
                or set(level) - required
                - {'shift', 'ipea', 'extrapolation_power'}):
            raise ValueError('Each spin-splitting gap level requires label, '
                             'method, and basis; shift, ipea, and the CBS '
                             'power are optional.')
        if (any(not isinstance(level[key], str) or not level[key].strip()
                for key in ('label', 'method'))):
            raise ValueError('Spin-splitting gap-level strings cannot be '
                             'empty.')
        basis = level['basis']
        if isinstance(basis, str):
            if not basis.strip():
                raise ValueError('Gap-level basis cannot be empty.')
        elif isinstance(basis, list) and len(basis) == 2 \
                and all(isinstance(item, str) and item.strip()
                        for item in basis):
            lower, upper = map(_basis_cardinal, basis)
            if upper != lower + 1:
                raise ValueError('Gap-level CBS bases must be adjacent.')
            extrapolation_coefficient(
                upper, power=level.get('extrapolation_power',
                                       ORIGINAL_ANL_POWER))
        else:
            raise ValueError('Gap-level basis must be a string or a '
                             'two-item list of adjacent cardinal bases.')
        labels.append(level['label'])
        for option in ('shift', 'ipea'):
            if option in level and (isinstance(level[option], bool)
                                    or not isinstance(level[option], Real)
                                    or not np.isfinite(level[option])):
                raise ValueError(f'Spin-splitting {option} must be finite.')
    if len(set(labels)) != len(labels):
        raise ValueError('Spin-splitting gap-level labels must be unique.')

    weighting = settings['weighting']
    if not isinstance(weighting, dict):
        raise ValueError('spin-splitting weighting must be an object.')
    model = str(weighting.get('model', '')).casefold()
    if model == 'arithmetic_mean':
        return
    if model != 'distance_switch':
        raise ValueError('spin-splitting weighting model must be '
                         'distance_switch or arithmetic_mean.')
    knots = weighting.get('distance_knots_angstrom')
    weights = weighting.get('weights')
    try:
        knot_array = np.asarray(knots, dtype=float)
    except (TypeError, ValueError):
        knot_array = np.asarray([])
    if (not isinstance(knots, list) or len(knots) < 2
            or knot_array.shape != (len(knots),)
            or not np.isfinite(knot_array).all()
            or np.any(np.diff(knot_array) <= 0.)):
        raise ValueError('Spin-splitting distance knots must be a strictly '
                         'increasing finite list.')
    if not isinstance(weights, dict) or set(weights) != set(labels):
        raise ValueError('Spin-splitting weights must match gap-level labels.')
    arrays = []
    for label in labels:
        values = weights[label]
        try:
            values_array = np.asarray(values, dtype=float)
        except (TypeError, ValueError):
            values_array = np.asarray([])
        if (not isinstance(values, list) or len(values) != len(knots)
                or values_array.shape != (len(knots),)
                or not np.isfinite(values_array).all()
                or np.any(values_array < 0.)):
            raise ValueError('Each spin-splitting weight series must be '
                             'finite, nonnegative, and match the knots.')
        arrays.append(values_array)
    if not np.allclose(sum(arrays), 1., rtol=0., atol=1.e-12):
        raise ValueError('Spin-splitting weights must sum to one at every '
                         'distance knot.')


def correction_components(parameters):
    """Describe the independent energies required at every VRC scan point."""
    sampling = sampling_settings(parameters)
    sample_bases = (sampling['basis'] if isinstance(sampling['basis'], list)
                    else [sampling['basis']])
    sample_components = []
    for index, basis in enumerate(sample_bases):
        sample_components.append({
            'role': ('sampling' if len(sample_bases) == 1
                     else f'sampling:{index}'),
            'suffix': ('_fr' if len(sample_bases) == 1
                       else f'_fr_sampling_{index}'),
            'geometry': 'frozen',
            'method': sampling['method'],
            'basis': basis,
            'multiplicity': None,
        })
    if correction_model(parameters) == 'direct':
        return [*sample_components, {
                'role': 'trusted', 'suffix': '', 'geometry': 'relaxed',
                'method': parameters['vrc_tst_high_method'],
                'basis': parameters['vrc_tst_high_basis'],
                'multiplicity': None,
            }]

    settings = spin_splitting_settings(parameters)
    legacy_states = (settings['target_multiplicity'] == 1
                     and settings['reference_multiplicity'] == 3)
    components = []
    for component in sample_components:
        component = dict(component)
        component['multiplicity'] = settings['target_multiplicity']
        components.append(component)
    reference = settings['reference']
    bases = reference['basis'] if isinstance(reference['basis'], list) \
        else [reference['basis']]
    for index, basis in enumerate(bases):
        stem = 'triplet_reference' if legacy_states else 'reference'
        role = (stem if len(bases) == 1 else f'{stem}:{index}')
        suffix_stem = 'triplet_cc' if legacy_states else 'reference_cc'
        suffix = (f'_ss_{suffix_stem}' if len(bases) == 1
                  else f'_ss_{suffix_stem}_{index}')
        components.append({
            'role': role, 'suffix': suffix, 'geometry': 'relaxed',
            'method': reference['method'], 'basis': basis,
            'multiplicity': settings['reference_multiplicity'],
        })
    for level in settings['gap_levels']:
        common = {key: value for key, value in level.items()
                  if key in ('method', 'shift', 'ipea')}
        label = level['label']
        bases = level['basis'] if isinstance(level['basis'], list) \
            else [level['basis']]
        state_pairs = ((('singlet', settings['target_multiplicity']),
                        ('triplet', settings['reference_multiplicity']))
                       if legacy_states else
                       (('target', settings['target_multiplicity']),
                        ('reference', settings['reference_multiplicity'])))
        for state, multiplicity in state_pairs:
            for index, basis in enumerate(bases):
                role = (f'gap:{label}:{state}' if len(bases) == 1
                        else f'gap:{label}:{state}:{index}')
                suffix = (f'_ss_{label}_{state}' if len(bases) == 1
                          else f'_ss_{label}_{state}_{index}')
                components.append({
                    'role': role, 'suffix': suffix,
                    'geometry': 'relaxed', 'multiplicity': multiplicity,
                    'basis': basis, **common,
                })
    return components


def expected_levels(parameters):
    """Return the immutable level declaration stored with a correction."""
    sampling = sampling_settings(parameters)
    if correction_model(parameters) == 'direct':
        trusted = {
            'method': parameters['vrc_tst_high_method'],
            'basis': parameters['vrc_tst_high_basis'],
        }
    else:
        settings = spin_splitting_settings(parameters)
        trusted = {
            'model': 'spin_splitting',
            'target_multiplicity': settings['target_multiplicity'],
            'reference_multiplicity': settings['reference_multiplicity'],
            'fragment_multiplicities': settings['fragment_multiplicities'],
            'reference': settings['reference'],
            'spin_gaps': settings['gap_levels'],
            'weighting': settings['weighting'],
        }
    return {'sampling': sampling, 'trusted_correction': trusted}


def sampling_settings(parameters):
    """Return the frozen sampling level and optional ANL-style CBS pair."""
    method = parameters.get('vrc_tst_sample_method')
    basis = parameters.get('vrc_tst_sample_basis')
    if not isinstance(method, str) or not method.strip():
        raise ValueError('vrc_tst_sample_method must be a nonempty string.')
    if isinstance(basis, str):
        if not basis.strip():
            raise ValueError('vrc_tst_sample_basis cannot be empty.')
        return {'method': method, 'basis': basis}
    if (not isinstance(basis, list) or len(basis) != 2
            or any(not isinstance(item, str) or not item.strip()
                   for item in basis)):
        raise ValueError('vrc_tst_sample_basis must be one basis or an '
                         'adjacent two-basis CBS pair.')
    lower, upper = map(_basis_cardinal, basis)
    if upper != lower + 1:
        raise ValueError('VRC sampling CBS bases must be adjacent.')
    power = parameters.get('vrc_tst_sample_extrapolation_power',
                           ORIGINAL_ANL_POWER)
    extrapolation_coefficient(upper, power=power)
    return {'method': method, 'basis': list(basis),
            'extrapolation_power': power}


def assemble_sampling_reference(distances, energies, settings):
    """Return the sampling curve, applying its optional CBS extrapolation."""
    basis = settings['basis']
    if isinstance(basis, str):
        values = np.asarray(energies.get('sampling', []), dtype=float)
        if len(values) != len(distances) or not np.isfinite(values).all():
            raise ValueError('Incomplete or nonfinite VRC sampling series.')
        return values.tolist(), None
    lower = energies.get('sampling:0', [])
    upper = energies.get('sampling:1', [])
    if len(lower) != len(distances) or len(upper) != len(distances):
        raise ValueError('Incomplete VRC sampling CBS series.')
    lower_cardinal, upper_cardinal = map(_basis_cardinal, basis)
    power = settings['extrapolation_power']
    values = np.asarray([
        two_point_cbs(float(lo), float(hi),
                      upper_cardinal=upper_cardinal, power=power)
        for lo, hi in zip(lower, upper)], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError('VRC sampling CBS energies must be finite.')
    return values.tolist(), {
        'formula': 'E_n + alpha_n*(E_n-E_(n-1))',
        'basis': list(basis),
        'lower_cardinal': lower_cardinal,
        'upper_cardinal': upper_cardinal,
        'power': power,
        'alpha': extrapolation_coefficient(upper_cardinal, power=power),
        'lower_hartree': list(lower),
        'upper_hartree': list(upper),
        'extrapolated_hartree': values.tolist(),
    }


def spin_splitting_weights(distances, settings):
    """Return one normalized gap weight per level and distance."""
    labels = [level['label'] for level in settings['gap_levels']]
    model = settings['weighting']['model'].casefold()
    if model == 'arithmetic_mean':
        return {label: [1. / len(labels)] * len(distances)
                for label in labels}
    knots = settings['weighting']['distance_knots_angstrom']
    return {
        label: np.interp(distances, knots,
                         settings['weighting']['weights'][label]).tolist()
        for label in labels
    }


def assemble_spin_splitting_reference(distances, energies, settings):
    """Reconstruct a target-spin surface from a reference and MR gaps."""
    distances = np.asarray(distances, dtype=float)
    reference_settings = settings.get(
        'reference', settings.get('triplet_reference'))
    if reference_settings is None:
        raise ValueError('Spin-splitting reference settings are missing.')
    cbs_metadata = None
    if isinstance(reference_settings['basis'], list):
        lower = energies.get(
            'reference:0', energies.get('triplet_reference:0', []))
        upper = energies.get(
            'reference:1', energies.get('triplet_reference:1', []))
        if len(lower) != len(distances) or len(upper) != len(distances):
            raise ValueError('Incomplete spin-splitting reference CBS series.')
        lower_cardinal, upper_cardinal = map(
            _basis_cardinal, reference_settings['basis'])
        power = reference_settings.get('extrapolation_power',
                                       ORIGINAL_ANL_POWER)
        reference_energy = np.asarray([
            two_point_cbs(float(lo), float(hi),
                          upper_cardinal=upper_cardinal, power=power)
            for lo, hi in zip(lower, upper)], dtype=float)
        cbs_metadata = {
            'formula': 'E_n + alpha_n*(E_n-E_(n-1))',
            'basis': list(reference_settings['basis']),
            'lower_cardinal': lower_cardinal,
            'upper_cardinal': upper_cardinal,
            'power': power,
            'alpha': extrapolation_coefficient(upper_cardinal, power=power),
            'lower_hartree': list(lower),
            'upper_hartree': list(upper),
            'extrapolated_hartree': reference_energy.tolist(),
        }
    else:
        reference_energy = np.asarray(
            energies.get('reference', energies.get('triplet_reference', [])),
            dtype=float)
    if len(reference_energy) != len(distances):
        raise ValueError('Incomplete spin-splitting reference series.')
    weights = spin_splitting_weights(distances, settings)
    result = reference_energy.copy()
    reconstructed = {}
    gap_cbs = {}
    for level in settings['gap_levels']:
        label = level['label']
        target, target_cbs = _state_series(
            distances, energies, level, label, 'target')
        reference_state, reference_cbs = _state_series(
            distances, energies, level, label, 'reference')
        if len(target) != len(distances) or len(reference_state) != len(distances):
            raise ValueError(f'Incomplete spin-splitting gap series: {label}.')
        candidate = reference_energy + target - reference_state
        reconstructed[label] = candidate.tolist()
        result += np.asarray(weights[label]) * (target - reference_state)
        if target_cbs is not None:
            gap_cbs[label] = {
                'target': target_cbs,
                'reference': reference_cbs,
                'gap_hartree': (target - reference_state).tolist(),
            }
    if not np.isfinite(result).all():
        raise ValueError('Spin-splitting reference energies must be finite.')
    metadata = {
        'formula': ('E_target = E_reference_CC + sum_i w_i(R) '
                    '[E_target_i - E_reference_i]'),
        'target_multiplicity': settings.get(
            'target_multiplicity', settings.get('singlet_multiplicity')),
        'reference_multiplicity': settings.get(
            'reference_multiplicity', settings.get('triplet_multiplicity')),
        'weights': weights,
        'reconstructed_target_hartree': reconstructed,
    }
    if cbs_metadata is not None:
        metadata['reference_cbs'] = cbs_metadata
        metadata['triplet_cbs'] = cbs_metadata
    if gap_cbs:
        metadata['gap_cbs'] = gap_cbs
    return result.tolist(), metadata


def _state_series(distances, energies, level, label, state):
    """Read one MR state series, applying the ANL CBS formula when requested."""
    basis = level['basis']
    legacy_state = {'target': 'singlet', 'reference': 'triplet'}.get(state)

    def series(key):
        value = energies.get(key)
        if value is None and legacy_state is not None:
            value = energies.get(key.replace(f':{state}',
                                             f':{legacy_state}'))
        return [] if value is None else value

    if isinstance(basis, str):
        return np.asarray(series(f'gap:{label}:{state}'),
                          dtype=float), None
    lower = series(f'gap:{label}:{state}:0')
    upper = series(f'gap:{label}:{state}:1')
    if len(lower) != len(distances) or len(upper) != len(distances):
        raise ValueError(f'Incomplete spin-splitting CBS series: '
                         f'{label} {state}.')
    lower_cardinal, upper_cardinal = map(_basis_cardinal, basis)
    power = level.get('extrapolation_power', ORIGINAL_ANL_POWER)
    values = np.asarray([
        two_point_cbs(float(lo), float(hi),
                      upper_cardinal=upper_cardinal, power=power)
        for lo, hi in zip(lower, upper)], dtype=float)
    return values, {
        'formula': 'E_n + alpha_n*(E_n-E_(n-1))',
        'basis': list(basis),
        'lower_cardinal': lower_cardinal,
        'upper_cardinal': upper_cardinal,
        'power': power,
        'alpha': extrapolation_coefficient(upper_cardinal, power=power),
        'lower_hartree': list(lower),
        'upper_hartree': list(upper),
        'extrapolated_hartree': values.tolist(),
    }


def _basis_cardinal(basis):
    """Read D/T/Q/5/... cardinality from a correlation-consistent basis."""
    compact = re.sub(r'[^a-z0-9]', '', str(basis).casefold())
    match = re.search(r'(?:aug)?ccpv([dtq]|[3-9])z', compact)
    if match is None:
        match = re.search(r'a?v([dtq]|[3-9])z', compact)
    if match is None:
        raise ValueError(f'Cannot infer cardinal number from basis {basis!r}.')
    return {'d': 2, 't': 3, 'q': 4}.get(
        match.group(1), int(match.group(1)) if match.group(1).isdigit()
        else None)
