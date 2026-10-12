"""VRC correction-potential models and their calculation provenance."""

from copy import deepcopy
from numbers import Real
import re

import numpy as np

from kinbot.anl.extrapolation import (ORIGINAL_ANL_POWER,
                                      extrapolation_coefficient,
                                      two_point_cbs)


_SPIN_SPLITTING_DEFAULTS = {
    'singlet_multiplicity': 1,
    'triplet_multiplicity': 3,
    'triplet_reference': {
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


def spin_splitting_settings(parameters):
    """Return validated spin-splitting settings with concise user overrides."""
    settings = deepcopy(_SPIN_SPLITTING_DEFAULTS)
    override = parameters.get('vrc_tst_spin_splitting')
    if override is None:
        override = {}
    if not isinstance(override, dict):
        raise ValueError('vrc_tst_spin_splitting must be an object.')
    allowed = set(settings)
    unknown = set(override) - allowed
    if unknown:
        raise ValueError('Unknown vrc_tst_spin_splitting option(s): '
                         + ', '.join(sorted(unknown)))
    for key, value in override.items():
        if key in ('triplet_reference', 'weighting'):
            if not isinstance(value, dict):
                raise ValueError(f'vrc_tst_spin_splitting.{key} must be an '
                                 'object.')
            settings[key].update(value)
        else:
            settings[key] = value
    if isinstance(settings.get('triplet_reference'), dict) \
            and isinstance(settings['triplet_reference'].get('basis'), str):
        settings['triplet_reference'].pop('extrapolation_power', None)
    for level in settings.get('gap_levels', []):
        if isinstance(level, dict) and isinstance(level.get('basis'), str):
            level.pop('extrapolation_power', None)
    if (isinstance(settings.get('weighting'), dict)
            and str(settings['weighting'].get('model', '')).casefold()
            == 'arithmetic_mean'):
        settings['weighting'] = {'model': 'arithmetic_mean'}
    _validate_spin_splitting_settings(settings)
    return settings


def _validate_spin_splitting_settings(settings):
    singlet = settings['singlet_multiplicity']
    triplet = settings['triplet_multiplicity']
    if (isinstance(singlet, bool) or not isinstance(singlet, int)
            or isinstance(triplet, bool) or not isinstance(triplet, int)
            or singlet != 1 or triplet != 3):
        raise ValueError('The implemented spin-splitting model requires '
                         'singlet multiplicity 1 and triplet multiplicity 3.')
    reference = settings['triplet_reference']
    if (not isinstance(reference, dict)
            or set(reference) - {'method', 'basis', 'extrapolation_power'}
            or not isinstance(reference.get('method'), str)
            or not reference['method'].strip()):
        raise ValueError('spin-splitting triplet_reference requires a method '
                         'and one basis or an adjacent basis pair.')
    basis = reference.get('basis')
    if isinstance(basis, str):
        if not basis.strip():
            raise ValueError('Triplet-reference basis cannot be empty.')
    elif isinstance(basis, list) and len(basis) == 2 \
            and all(isinstance(item, str) and item.strip() for item in basis):
        lower, upper = map(_basis_cardinal, basis)
        if upper != lower + 1:
            raise ValueError('Triplet-reference CBS bases must be adjacent.')
        power = reference.get('extrapolation_power', ORIGINAL_ANL_POWER)
        extrapolation_coefficient(upper, power=power)
    else:
        raise ValueError('Triplet-reference basis must be a string or a '
                         'two-item list of adjacent cardinal bases.')
    levels = settings['gap_levels']
    if not isinstance(levels, list) or len(levels) != 2:
        raise ValueError('spin-splitting gap_levels must contain exactly two '
                         'levels.')
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
    if correction_model(parameters) == 'direct':
        return [
            {
                'role': 'sampling', 'suffix': '_fr', 'geometry': 'frozen',
                'method': parameters['vrc_tst_sample_method'],
                'basis': parameters['vrc_tst_sample_basis'],
                'multiplicity': None,
            },
            {
                'role': 'trusted', 'suffix': '', 'geometry': 'relaxed',
                'method': parameters['vrc_tst_high_method'],
                'basis': parameters['vrc_tst_high_basis'],
                'multiplicity': None,
            },
        ]

    settings = spin_splitting_settings(parameters)
    components = [{
        'role': 'sampling', 'suffix': '_fr', 'geometry': 'frozen',
        'method': parameters['vrc_tst_sample_method'],
        'basis': parameters['vrc_tst_sample_basis'],
        'multiplicity': settings['singlet_multiplicity'],
    }]
    triplet = settings['triplet_reference']
    bases = triplet['basis'] if isinstance(triplet['basis'], list) \
        else [triplet['basis']]
    for index, basis in enumerate(bases):
        role = ('triplet_reference' if len(bases) == 1
                else f'triplet_reference:{index}')
        suffix = ('_ss_triplet_cc' if len(bases) == 1
                  else f'_ss_triplet_cc_{index}')
        components.append({
            'role': role, 'suffix': suffix, 'geometry': 'relaxed',
            'method': triplet['method'], 'basis': basis,
            'multiplicity': settings['triplet_multiplicity'],
        })
    for level in settings['gap_levels']:
        common = {key: value for key, value in level.items()
                  if key in ('method', 'shift', 'ipea')}
        label = level['label']
        bases = level['basis'] if isinstance(level['basis'], list) \
            else [level['basis']]
        for state, multiplicity in (
                ('singlet', settings['singlet_multiplicity']),
                ('triplet', settings['triplet_multiplicity'])):
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
    sampling = {
        'method': parameters['vrc_tst_sample_method'],
        'basis': parameters['vrc_tst_sample_basis'],
    }
    if correction_model(parameters) == 'direct':
        trusted = {
            'method': parameters['vrc_tst_high_method'],
            'basis': parameters['vrc_tst_high_basis'],
        }
    else:
        settings = spin_splitting_settings(parameters)
        trusted = {
            'model': 'spin_splitting',
            'triplet_reference': settings['triplet_reference'],
            'singlet_triplet_gaps': settings['gap_levels'],
            'weighting': settings['weighting'],
        }
    return {'sampling': sampling, 'trusted_correction': trusted}


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
    """Reconstruct singlet energies from CC triplets and weighted MR gaps."""
    distances = np.asarray(distances, dtype=float)
    triplet_settings = settings['triplet_reference']
    cbs_metadata = None
    if isinstance(triplet_settings['basis'], list):
        lower = energies.get('triplet_reference:0', [])
        upper = energies.get('triplet_reference:1', [])
        if len(lower) != len(distances) or len(upper) != len(distances):
            raise ValueError('Incomplete spin-splitting triplet CBS series.')
        lower_cardinal, upper_cardinal = map(
            _basis_cardinal, triplet_settings['basis'])
        power = triplet_settings.get('extrapolation_power',
                                     ORIGINAL_ANL_POWER)
        triplet_reference = np.asarray([
            two_point_cbs(float(lo), float(hi),
                          upper_cardinal=upper_cardinal, power=power)
            for lo, hi in zip(lower, upper)], dtype=float)
        cbs_metadata = {
            'formula': 'E_n + alpha_n*(E_n-E_(n-1))',
            'basis': list(triplet_settings['basis']),
            'lower_cardinal': lower_cardinal,
            'upper_cardinal': upper_cardinal,
            'power': power,
            'alpha': extrapolation_coefficient(upper_cardinal, power=power),
            'lower_hartree': list(lower),
            'upper_hartree': list(upper),
            'extrapolated_hartree': triplet_reference.tolist(),
        }
    else:
        triplet_reference = np.asarray(
            energies.get('triplet_reference', []), dtype=float)
    if len(triplet_reference) != len(distances):
        raise ValueError('Incomplete spin-splitting triplet-reference series.')
    weights = spin_splitting_weights(distances, settings)
    reference = triplet_reference.copy()
    reconstructed = {}
    gap_cbs = {}
    for level in settings['gap_levels']:
        label = level['label']
        singlet, singlet_cbs = _state_series(
            distances, energies, level, label, 'singlet')
        triplet, triplet_cbs = _state_series(
            distances, energies, level, label, 'triplet')
        if len(singlet) != len(distances) or len(triplet) != len(distances):
            raise ValueError(f'Incomplete spin-splitting gap series: {label}.')
        candidate = triplet_reference + singlet - triplet
        reconstructed[label] = candidate.tolist()
        reference += np.asarray(weights[label]) * (singlet - triplet)
        if singlet_cbs is not None:
            gap_cbs[label] = {
                'singlet': singlet_cbs,
                'triplet': triplet_cbs,
                'gap_hartree': (singlet - triplet).tolist(),
            }
    if not np.isfinite(reference).all():
        raise ValueError('Spin-splitting reference energies must be finite.')
    metadata = {
        'formula': ('E_singlet = E_triplet_CC + sum_i w_i(R) '
                    '[E_singlet_i - E_triplet_i]'),
        'weights': weights,
        'reconstructed_singlet_hartree': reconstructed,
    }
    if cbs_metadata is not None:
        metadata['triplet_cbs'] = cbs_metadata
    if gap_cbs:
        metadata['gap_cbs'] = gap_cbs
    return reference.tolist(), metadata


def _state_series(distances, energies, level, label, state):
    """Read one MR state series, applying the ANL CBS formula when requested."""
    basis = level['basis']
    if isinstance(basis, str):
        return np.asarray(energies.get(f'gap:{label}:{state}', []),
                          dtype=float), None
    lower = energies.get(f'gap:{label}:{state}:0', [])
    upper = energies.get(f'gap:{label}:{state}:1', [])
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
