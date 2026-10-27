"""VRC correction-potential models and their calculation provenance."""

from copy import deepcopy
import json
from numbers import Real
from pathlib import Path
import re

import numpy as np
from ase.data import atomic_numbers

from kinbot.anl.extrapolation import (ORIGINAL_ANL_POWER,
                                      extrapolation_coefficient,
                                      two_point_cbs)
from kinbot.molpro import (automatic_vrc_state_model,
                           full_valence_active_space,
                           method_with_active_space,
                           validate_vrc_single_reference_root,
                           validate_vrc_state_model,
                           vrc_multireference_energy_contract)


# Use the requested ANL exponent for the two-point gap estimate as well.
# Retain its basis sensitivity; this is not a calibrated error certificate.
CASPT2_GAP_CBS_POWER = ORIGINAL_ANL_POWER
_SPIN_SPLITTING_DEFAULTS = {
    'target_multiplicity': 1,
    'reference_multiplicity': 3,
    'fragment_multiplicities': None,
    'reference': {'method': 'uccsd(t)-f12b',
                  'basis': ['avdz-f12', 'avtz-f12'],
                  'extrapolation_power': ORIGINAL_ANL_POWER},
    'gap_levels': [{
        'label': 'caspt2_f12', 'method': 'caspt2-f12(2,2)',
        'basis': ['avdz-f12', 'avtz-f12'],
        'extrapolation_power': CASPT2_GAP_CBS_POWER,
        'shift': 0.1, 'ipea': 0.25,
    }],
    'weighting': {'model': 'arithmetic_mean'},
}

_VRC_ACCURACY_PRESETS = {
    'legacy': {},
    # CAS(2,2) is a method template, replaced by the reviewed reaction CAS.
    # No chemical active space, roots or physical surfaces are inferred here.
    'anl1-f12': {
        'vrc_tst_scan_method': 'caspt2-f12(2,2)',
        'vrc_tst_scan_basis': 'avtz-f12',
        'vrc_tst_scan_gradient': 'numerical',
        'vrc_tst_scan_shift': 0.1,
        'vrc_tst_scan_ipea': 0.25,
        'vrc_tst_sample_method': 'caspt2-f12(2,2)',
        'vrc_tst_sample_basis': 'avdz-f12',
        'vrc_tst_high_method': 'caspt2-f12(2,2)',
        'vrc_tst_high_basis': ['avdz-f12', 'avtz-f12'],
        'vrc_tst_high_extrapolation_power': CASPT2_GAP_CBS_POWER,
        'vrc_tst_correction_model': 'spin_splitting',
        'vrc_tst_spin_splitting': {
            key: deepcopy(_SPIN_SPLITTING_DEFAULTS[key])
            for key in ('reference', 'gap_levels', 'weighting')},
    },
}


def apply_vrc_accuracy_preset(parameters, user_keys=()):
    """Fill VRC production defaults without replacing expert input.

    ``legacy`` leaves the historical settings untouched. ``anl1-f12`` is a
    compact request for augmented-basis, spin-split VRC levels. It requires an
    explicit active-space policy or orbital/surface model; it does not decide
    the chemical space. Explicit full_valence remains supported.
    Every explicitly supplied top-level option wins over the preset.
    """
    name = str(parameters.get(
        'vrc_tst_accuracy_preset', 'legacy')).strip().casefold()
    aliases = {'legacy': 'legacy', 'default': 'legacy',
               'anl1-f12': 'anl1-f12', 'anl1_f12': 'anl1-f12'}
    try:
        name = aliases[name]
    except KeyError as exc:
        raise ValueError('vrc_tst_accuracy_preset must be legacy or '
                         'anl1-f12.') from exc
    explicit = set(user_keys)
    choices = ('vrc_tst_active_space', 'vrc_tst_state_model_file',
               'vrc_tst_state_model', 'vrc_tst_surface_model_file', 'vrc_tst_surface_model')
    if name != 'legacy' and not any(key in explicit and parameters.get(key) for key in choices):
        raise ValueError('The anl1-f12 VRC preset specifies energy levels, not a chemical active space. '
                         'Supply vrc_tst_active_space explicitly or install an orbital/state or '
                         'surface model after orbital review.')
    parameters['vrc_tst_accuracy_preset'] = name
    for key, value in _VRC_ACCURACY_PRESETS[name].items():
        if key not in explicit:
            parameters[key] = deepcopy(value)
    if correction_model(parameters) == 'direct' and 'vrc_tst_spin_splitting' not in explicit:
        parameters['vrc_tst_spin_splitting'] = None
    space = active_space_settings(parameters)
    if space is not None and space.get('selection') == 'selected_orbital_subspace':
        apply_selected_active_space_methods(parameters, space)
    return parameters


def apply_selected_active_space_methods(parameters, active_space):
    """Keep human-readable method labels consistent with the selected CAS."""
    for key in ('vrc_tst_scan_method', 'vrc_tst_sample_method', 'vrc_tst_high_method'):
        if parameters.get(key):
            parameters[key] = method_with_active_space(parameters[key], active_space)
    if isinstance(parameters.get('vrc_tst_spin_splitting'), dict):
        parameters['vrc_tst_spin_splitting'] = deepcopy(parameters['vrc_tst_spin_splitting'])
        for level in parameters['vrc_tst_spin_splitting'].get('gap_levels', []):
            if isinstance(level, dict) and 'method' in level:
                level['method'] = method_with_active_space(level['method'], active_space)
    return parameters


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


def asymptote_distance(parameters):
    """Resolve the separated-fragment reference independently of job names.

    Thirty angstrom is the legacy default, not a convergence certificate.
    Reconnaissance can establish a different value for the production input.
    """
    value = parameters.get('vrc_tst_asymptote_distance', 30.)
    if (isinstance(value, bool) or not isinstance(value, Real)
            or not np.isfinite(value) or value <= 0.):
        raise ValueError('vrc_tst_asymptote_distance must be positive and finite.')
    points = parameters.get('vrc_tst_scan_points', [])
    if any(float(point) >= value for point in points):
        raise ValueError('vrc_tst_asymptote_distance must exceed every VRC scan point.')
    return float(value)


def spin_splitting_settings(parameters, *, parent_multiplicity=None,
                            fragment_multiplicities=None):
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
    explicit = set(override)
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
    if parent_multiplicity is not None \
            and 'target_multiplicity' not in explicit:
        settings['target_multiplicity'] = int(parent_multiplicity)
    if fragment_multiplicities is not None \
            and 'fragment_multiplicities' not in explicit:
        settings['fragment_multiplicities'] = [
            int(value) for value in fragment_multiplicities]
    declared_fragments = settings['fragment_multiplicities']
    if declared_fragments is not None and (
            not isinstance(declared_fragments, list) or len(declared_fragments) != 2
            or any(isinstance(value, bool) or not isinstance(value, int) or value < 1
                   for value in declared_fragments)):
        raise ValueError('fragment_multiplicities must contain two positive integers.')
    if (parent_multiplicity is not None
            and settings['target_multiplicity'] != parent_multiplicity):
        raise ValueError('Spin-splitting target multiplicity differs from the '
                         'accepted parent surface.')
    if (fragment_multiplicities is not None
            and sorted(settings['fragment_multiplicities'])
            != sorted(fragment_multiplicities)):
        raise ValueError('Spin-splitting fragment multiplicities differ from '
                         'the accepted products.')
    if (settings['fragment_multiplicities'] is not None
            and 'reference_multiplicity' not in explicit):
        allowed = coupled_multiplicities(
            *settings['fragment_multiplicities'])
        alternatives = [value for value in allowed
                        if value != settings['target_multiplicity']]
        if not alternatives:
            raise ValueError(
                'The accepted fragments provide only one coupled spin '
                'surface; specify reference_multiplicity explicitly for '
                'a separately reviewed energy-reference state.')
        # The single-reference anchor is the maximum-spin coupling. Falling
        # back to a lower-spin anchor for an already-high-spin target reverses
        # the physical justification for this correction. Expert references
        # remain explicit; ordinary high-spin channels should use direct QC.
        if settings['target_multiplicity'] == max(allowed):
            raise ValueError('The target is already the highest coupled spin; '
                             'use a direct correction or explicitly justify '
                             'a different reference multiplicity.')
        settings['reference_multiplicity'] = max(allowed)
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
    # Keep archived manifests and callers readable while new records
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
    if target % 2 != reference_multiplicity % 2:
        raise ValueError('Spin-splitting target and reference multiplicities '
                         'must have the same electron-spin parity.')
    fragments = settings['fragment_multiplicities']
    if fragments is not None:
        if (not isinstance(fragments, list) or len(fragments) != 2
                or any(isinstance(value, bool) or not isinstance(value, int)
                       or value < 1 for value in fragments)):
            raise ValueError('fragment_multiplicities must contain two '
                             'positive integers.')
        allowed = coupled_multiplicities(*fragments)
        if target not in allowed:
            raise ValueError('Spin-splitting target multiplicity is not allowed '
                             'by the two fragment spins.')
        # The reference is an auxiliary energy calculation, not another
        # reactive surface. An explicitly selected reference may correlate
        # with excited fragment terms. Its CC and MR states must match; the
        # actual products and target-spin statistical weights stay unchanged.
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


def active_space_settings(parameters, symbols=None, charge=0):
    """Return the requested active-space policy for one VRC molecule.

    ``full_valence`` uses atom-labelled AVAS targets. ``radical_pair`` is the
    deliberately narrow CAS(2,2) shortcut for two-doublet recombination; its
    fragment-spin contract is checked when the accepted reaction is known.
    """
    request = parameters.get('vrc_tst_active_space', 'explicit')
    additions = []
    if isinstance(request, dict):
        unknown = set(request) - {
            'mode', 'additional_orbitals', 'active_electrons',
            'active_orbitals', 'source',
        }
        if unknown:
            raise ValueError('Unknown vrc_tst_active_space option(s): '
                             + ', '.join(sorted(unknown)))
        mode = str(request.get('mode', '')).strip().casefold()
        additions = request.get('additional_orbitals', [])
        if not isinstance(additions, list):
            raise ValueError('additional_orbitals must be a list.')
        normalized = []
        for entry in additions:
            if (not isinstance(entry, dict)
                    or set(entry) != {'atom', 'orbitals'}
                    or isinstance(entry['atom'], bool)
                    or not isinstance(entry['atom'], int)
                    or entry['atom'] < 1
                    or not isinstance(entry['orbitals'], list)
                    or not entry['orbitals']):
                raise ValueError('Each additional orbital request requires a '
                                 'positive atom index and an orbital list.')
            normalized.append({
                'atom': entry['atom'],
                'orbitals': [str(value).strip().casefold()
                             for value in entry['orbitals']],
            })
        additions = normalized
    else:
        mode = str(request).strip().casefold()
    if mode == 'explicit':
        if isinstance(request, dict) and (
                'active_electrons' in request
                or 'active_orbitals' in request):
            raise ValueError(
                'Counted active spaces use mode=selected; string explicit '
                'retains the CAS written in each method name.')
        if additions:
            raise ValueError('additional_orbitals requires full_valence mode.')
        return None
    if mode in ('radical_pair', 'auto_radical_pair'):
        if additions:
            raise ValueError('additional_orbitals requires full_valence mode.')
        return {
            'active_electrons': 2,
            'active_orbitals': 2,
            'orbital_selection': None,
            'selection': ('verified_radical_pair' if mode == 'auto_radical_pair'
                          else 'two_doublet_radical_pair'),
            'validation': ('native_verification_required' if mode == 'auto_radical_pair'
                           else 'accepted_products_must_be_two_doublets'),
        }
    if mode in ('selected', 'targeted'):
        if not isinstance(request, dict):
            raise ValueError(
                'vrc_tst_active_space=selected requires an object with '
                'active_electrons and active_orbitals.')
        if additions:
            raise ValueError(
                'additional_orbitals is an AVAS full_valence option.')
        electrons = request.get('active_electrons')
        orbitals = request.get('active_orbitals')
        if (isinstance(electrons, bool) or not isinstance(electrons, int)
                or electrons < 1 or isinstance(orbitals, bool)
                or not isinstance(orbitals, int) or orbitals < 1
                or electrons > 2 * orbitals):
            raise ValueError(
                'A selected active space requires compatible positive '
                'active_electrons and active_orbitals.')
        source = str(request.get('source', 'user-selected')).strip()
        if not source:
            raise ValueError('A selected active-space source cannot be empty.')
        return {
            'active_electrons': electrons,
            'active_orbitals': orbitals,
            'orbital_selection': None,
            'selection': 'selected_orbital_subspace',
            'source': source,
            'validation': 'explicit_orbital_identity_required',
        }
    if mode != 'full_valence':
        raise ValueError(
            'vrc_tst_active_space must be explicit, radical_pair, auto_radical_pair, '
            'selected, full_valence, or an object with one of those modes.')
    if symbols is None:
        return {'selection': 'whole_molecule_main_group_valence',
                'additional_orbitals': additions}
    return full_valence_active_space(
        symbols, charge, additional_centers=additions)


def _augmented_basis(basis):
    token = str(basis).strip().casefold().replace('_', '-').replace(' ', '')
    return (token.startswith('aug-')
            or re.fullmatch(r'a(?:v|cv|wcv)[dtq2-9]z(?:-f12)?', token)
            is not None)


def _resolve_level(level, active_space):
    result = deepcopy(level)
    if active_space is not None:
        if {'active_electrons', 'active_orbitals'} <= set(active_space):
            result['method'] = method_with_active_space(
                result['method'], active_space)
        contract = vrc_multireference_energy_contract(result['method'])
        if contract is not None:
            bases = (result['basis'] if isinstance(result['basis'], list)
                     else [result['basis']])
            if not all(_augmented_basis(basis) for basis in bases):
                raise ValueError(
                    'Automatic active-space multireference VRC levels '
                    'require augmented orbital bases.')
    return result


def resolve_vrc_level(level, active_space):
    """Resolve a VRC level against an automatic molecular active space."""
    return _resolve_level(level, active_space)


def validate_vrc_levels(parameters):
    """Validate every configured VRC level before any scheduler submission."""
    active_space = active_space_settings(parameters)
    state_file = parameters.get('vrc_tst_state_model_file', '')
    if (not isinstance(state_file, str)
            or state_file.strip() and (active_space is None
                                      or 'active_electrons' not in active_space)):
        raise ValueError('vrc_tst_state_model_file requires a counted selected '
                         'active space.')
    if (active_space is not None
            and active_space.get('selection') == 'selected_orbital_subspace'):
        state_model = parameters.get('vrc_tst_state_model')
        state_file = parameters.get('vrc_tst_state_model_file', '')
        inline = (isinstance(state_model, dict)
                  and ((isinstance(state_model.get('orbital_selection'), dict)
                        and state_model['orbital_selection'].get('mode') in ('molpro_record', 'avas'))
                       or (state_model.get('schema') == 2 and bool(state_model.get('basis_models')))
                       or (state_model.get('schema') == 3 and bool(state_model.get('spin_models')))))
        external = isinstance(state_file, str) and bool(state_file.strip())
        if inline and external:
            raise ValueError('Configure either vrc_tst_state_model or '
                             'vrc_tst_state_model_file, not both.')
        if not inline and not external and not (parameters.get('vrc_tst_surface_model_file') or parameters.get('vrc_tst_surface_model')):
            raise ValueError(
                'A selected VRC active space requires '
                'vrc_tst_state_model.orbital_selection or '
                'vrc_tst_state_model_file with a reviewed Molpro record, '
                'or per-channel selections in vrc_tst_surface_model_file.')
    method = parameters.get('vrc_tst_scan_method')
    basis = parameters.get('vrc_tst_scan_basis')
    if not isinstance(method, str) or not method.strip():
        raise ValueError('vrc_tst_scan_method must be a nonempty string.')
    if not isinstance(basis, str) or not basis.strip():
        raise ValueError('vrc_tst_scan_basis must be a nonempty string.')
    _resolve_level({'method': method, 'basis': basis}, active_space)
    sampling_settings(parameters, active_space=active_space)
    model = correction_model(parameters)
    if (active_space is not None
            and active_space.get('selection') in ('two_doublet_radical_pair', 'verified_radical_pair')):
        if model != 'spin_splitting':
            raise ValueError(
                'The radical_pair active space requires the spin_splitting '
                'VRC correction model.')
        override = parameters.get('vrc_tst_spin_splitting') or {}
        fragments = override.get('fragment_multiplicities')
        if fragments is not None and fragments != [2, 2]:
            raise ValueError(
                'The radical_pair active space requires two doublet '
                'fragments.')
        spin = spin_splitting_settings(parameters)
        if {spin['target_multiplicity'], spin['reference_multiplicity']} \
                != {1, 3}:
            raise ValueError(
                'The radical_pair active space requires coupled singlet and '
                'triplet surfaces.')
    if model == 'direct':
        trusted_settings(parameters, active_space=active_space)
    else:
        settings = spin_splitting_settings(parameters)
        _resolve_level(settings['reference'], active_space)
        for level in settings['gap_levels']:
            _resolve_level(level, active_space)


def active_space_feasibility(active_space, multiplicities):
    """Report exact determinant dimensions for an automatic CAS.

    This is provenance, not a hidden cutoff. A site or user can decide that
    the resulting MRCI expansion is too expensive without KinBot silently
    changing the requested chemistry.
    """
    if active_space is None or 'active_electrons' not in active_space:
        return None
    from math import comb

    electrons = int(active_space['active_electrons'])
    orbitals = int(active_space['active_orbitals'])
    result = {}
    for multiplicity in dict.fromkeys(map(int, multiplicities)):
        twice_spin = multiplicity - 1
        if (electrons - twice_spin) % 2 or twice_spin > electrons:
            raise ValueError(
                f'CAS({electrons},{orbitals}) cannot represent multiplicity '
                f'{multiplicity}.')
        beta = (electrons - twice_spin) // 2
        alpha = electrons - beta
        if alpha > orbitals or beta > orbitals:
            raise ValueError(
                f'CAS({electrons},{orbitals}) cannot represent multiplicity '
                f'{multiplicity}.')
        result[str(multiplicity)] = {
            'alpha_electrons': alpha,
            'beta_electrons': beta,
            'determinants': comb(orbitals, alpha) * comb(orbitals, beta),
        }
    return result


def root_overrides(parameters, key):
    value = parameters.get(key, {})
    if not isinstance(value, dict):
        raise ValueError(f'{key} must map multiplicities to positive roots.')
    result = {}
    for mult, count in value.items():
        try:
            integer = int(mult)
        except (TypeError, ValueError) as error:
            raise ValueError(f'{key} needs positive integer multiplicities.') from error
        if (isinstance(mult, bool) or str(integer) != str(mult) or integer < 1
                or isinstance(count, bool) or not isinstance(count, int) or count < 1):
            raise ValueError(f'{key} needs positive integer multiplicities and roots.')
        result[integer] = count
    return result


def resolved_state_model(parameters, symbols, charge, settings, *, roots=None,
                         basis=None, multiplicity=None):
    """Resolve an explicit model or build the automatic whole-molecule model."""
    roots = {**({} if roots is None else {int(key): value for key, value in roots.items()}),
             **root_overrides(parameters, 'vrc_tst_state_average_roots')}
    selected_overrides = root_overrides(parameters, 'vrc_tst_selected_roots')
    explicit = deepcopy(parameters.get('vrc_tst_state_model'))
    external = False
    average_policy = 'channel_derived'
    state_file = parameters.get('vrc_tst_state_model_file', '')
    if explicit is not None and isinstance(state_file, str) and state_file.strip():
        raise ValueError('Configure either vrc_tst_state_model or '
                         'vrc_tst_state_model_file, not both.')
    payload = explicit if isinstance(explicit, dict) and explicit.get('schema') in (2, 3) else None
    path = Path(parameters.get('input_file') or Path.cwd() / 'input.json').expanduser().resolve()
    if explicit is None and isinstance(state_file, str) and state_file.strip():
        path = Path(state_file).expanduser()
        if not path.is_absolute():
            input_file = parameters.get('input_file')
            base = (Path(input_file).expanduser().resolve().parent
                    if input_file else Path.cwd())
            path = base / path
        path = path.resolve()
        if not path.is_file():
            raise FileNotFoundError(
                f'Reviewed VRC state-model file is unavailable: {path}. '
                'Run and audit the Molpro orbital reconnaissance first.')
        payload = json.loads(path.read_text())
    spin_specific = payload is not None and payload.get('schema') == 3
    if spin_specific:
        chosen_spin = int(multiplicity or (settings or {}).get('target_multiplicity')
                          or parameters.get('mult', 1))
        models = payload.get('spin_models', {})
        if str(chosen_spin) not in models:
            raise ValueError('Reviewed spin-specific model omits the requested multiplicity.')
        payload = deepcopy(models[str(chosen_spin)])
        if payload.get('schema') != 2 or payload.get('state_average_policy') != 'reviewed':
            raise ValueError('Spin-specific orbital models require reviewed basis-specific children.')
        roots = {key: value for key, value in roots.items() if key == chosen_spin}
        selected_overrides = {key: value for key, value in selected_overrides.items()
                              if key == chosen_spin}
    if payload is not None:
        if (not isinstance(payload, dict) or payload.get('schema') not in (1, 2)
                or not isinstance(payload.get('active_space'), dict)):
            raise ValueError('Reviewed VRC state-model file has an invalid '
                             'schema.')
        molecule = payload.get('molecule')
        if molecule is not None and symbols is not None:
            if (molecule.get('symbols') != list(map(str, symbols))
                    or molecule.get('charge') != charge):
                raise ValueError('Reviewed orbital record atom order or charge differs from the calculation.')
        average_policy = payload.get('state_average_policy', 'channel_derived')
        if average_policy not in ('reviewed', 'channel_derived'):
            raise ValueError('Unknown reviewed state_average_policy.')
        active_space = active_space_settings(parameters, symbols, charge)
        declared = payload['active_space']
        if (declared.get('active_electrons')
                != active_space.get('active_electrons')
                or declared.get('active_orbitals')
                != active_space.get('active_orbitals')):
            raise ValueError('Reviewed VRC state-model active space differs '
                             'from vrc_tst_active_space.')
        if payload['schema'] == 1:
            if not isinstance(payload.get('state_model'), dict):
                raise ValueError('Reviewed VRC state-model file has an '
                                 'invalid schema-1 state_model.')
            explicit = payload['state_model']
        else:
            models = payload.get('basis_models')
            if not isinstance(models, dict) or not models:
                raise ValueError('Reviewed VRC state-model file has no '
                                 'basis_models.')
            if not isinstance(basis, str) or not basis.strip():
                raise ValueError('A basis is required to select a reviewed '
                                 'Molpro orbital record.')
            key = basis.strip().casefold()
            matches = {str(name).strip().casefold(): model
                       for name, model in models.items()}
            if len(matches) != len(models):
                raise ValueError('Reviewed VRC basis names are ambiguous.')
            if key not in matches:
                raise ValueError(
                    f'Reviewed VRC state-model file has no orbital record '
                    f'for basis {basis!r}.')
            if not isinstance(matches[key], dict):
                raise ValueError('Reviewed VRC basis model must be an object.')
            explicit = matches[key]
        from kinbot.molpro import resolve_orbital_record_path
        explicit = resolve_orbital_record_path(explicit, path.parent)
        external = True
    if explicit is not None:
        active_space = active_space_settings(parameters, symbols, charge)
        # A reviewed external file certifies the orbital subspace.  The
        # accepted reaction still determines which coupled spin manifolds
        # and spatial components must be averaged for that channel.  This
        # lets one reviewed orbital model cover, for example, both C-S and
        # S-S cleavage without imposing the atomic-P three-root average on a
        # nondegenerate molecular fragment channel.
        if external and average_policy == 'channel_derived':
            multiplicities = ([settings['target_multiplicity'], settings['reference_multiplicity']]
                              if settings is not None else [int(parameters.get('mult', 1))])
            counts = {} if roots is None else {
                int(key): int(value) for key, value in roots.items()}
            explicit['state_average'] = [
                {
                    'multiplicity': multiplicity,
                    'roots': counts.get(multiplicity, 1),
                    'weights': ([1. / counts.get(multiplicity, 1)]
                                * counts.get(multiplicity, 1)),
                }
                for multiplicity in multiplicities
            ]
            selected = explicit.get('selected_roots', {})
            explicit['selected_roots'] = {
                multiplicity: selected.get(str(multiplicity),
                                          selected.get(multiplicity, 1))
                for multiplicity in multiplicities
            }
        elif external:
            counts = {entry['multiplicity']: entry['roots'] for entry in explicit['state_average']}
            overrides = root_overrides(parameters, 'vrc_tst_state_average_roots')
            if spin_specific:
                overrides = {key: value for key, value in overrides.items() if key == chosen_spin}
            if any(counts.get(mult) != count for mult, count in overrides.items()):
                raise ValueError('Root-count overrides conflict with the reviewed state average.')
            required = ([chosen_spin] if spin_specific else
                        [settings['target_multiplicity'], settings['reference_multiplicity']]
                        if settings is not None else [int(parameters.get('mult', 1))])
            if any(mult not in counts for mult in required):
                raise ValueError('Reviewed state average omits a required target/reference multiplicity.')
        if selected_overrides:
            chosen = {int(key): value for key, value in explicit.get('selected_roots', {}).items()}
            chosen.update(selected_overrides)
            explicit['selected_roots'] = chosen
        return validate_vrc_state_model(
            explicit,
            active_electrons=(active_space or {}).get('active_electrons'),
            active_orbitals=(active_space or {}).get('active_orbitals'),
            nelectron=(None if symbols is None else
                       sum(atomic_numbers[str(symbol)]
                           for symbol in symbols) - charge),
            atom_count=None if symbols is None else len(symbols))
    active_space = active_space_settings(parameters, symbols, charge)
    if active_space is not None and active_space.get('selection') == 'verified_radical_pair':
        raise ValueError('Automatic radical-pair verification has not supplied an accepted orbital model.')
    if active_space is None:
        return None
    if settings is None:
        multiplicities = [int(parameters.get('mult', 1) or 1)]
    else:
        multiplicities = [settings['target_multiplicity'],
                          settings['reference_multiplicity']]
    return automatic_vrc_state_model(
        symbols, charge, multiplicities, roots=roots,
        selected_roots=selected_overrides,
        active_space=active_space)


def correction_components(parameters, *, symbols=None, charge=0,
                          settings=None, state_roots=None):
    """Describe the independent energies required at every VRC scan point."""
    active_space = active_space_settings(parameters, symbols, charge)
    sampling = sampling_settings(parameters, active_space=active_space)
    if settings is None and correction_model(parameters) == 'spin_splitting':
        settings = spin_splitting_settings(parameters)
    def state_model_for(basis, multiplicity=None):
        return (deepcopy(parameters.get('vrc_tst_state_model'))
                if active_space is None else
                resolved_state_model(
                    parameters, symbols, charge, settings,
                    roots=state_roots, basis=basis, multiplicity=multiplicity))
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
            'state_model': state_model_for(basis),
        })
    if correction_model(parameters) == 'direct':
        trusted = trusted_settings(parameters, active_space=active_space)
        bases = trusted['basis'] if isinstance(trusted['basis'], list) else [trusted['basis']]
        return [*sample_components, *[{
                'role': 'trusted' if len(bases) == 1 else f'trusted:{index}',
                'suffix': '' if len(bases) == 1 else f'_trusted_{index}', 'geometry': 'relaxed',
                'method': trusted['method'],
                'basis': basis,
                'multiplicity': None,
                'state_model': state_model_for(basis),
            } for index, basis in enumerate(bases)]]

    settings = settings or spin_splitting_settings(parameters)
    legacy_states = (settings['target_multiplicity'] == 1
                     and settings['reference_multiplicity'] == 3)
    components = []
    for component in sample_components:
        component = dict(component)
        component['multiplicity'] = settings['target_multiplicity']
        components.append(component)
    reference = _resolve_level(settings['reference'], active_space)
    bases = reference['basis'] if isinstance(reference['basis'], list) \
        else [reference['basis']]
    for index, basis in enumerate(bases):
        reference_model = state_model_for(basis, settings['reference_multiplicity'])
        validate_vrc_single_reference_root(reference['method'], reference_model,
                                           settings['reference_multiplicity'])
        stem = 'triplet_reference' if legacy_states else 'reference'
        role = (stem if len(bases) == 1 else f'{stem}:{index}')
        suffix_stem = 'triplet_cc' if legacy_states else 'reference_cc'
        suffix = (f'_ss_{suffix_stem}' if len(bases) == 1
                  else f'_ss_{suffix_stem}_{index}')
        components.append({
            'role': role, 'suffix': suffix, 'geometry': 'relaxed',
            'method': reference['method'], 'basis': basis,
            'multiplicity': settings['reference_multiplicity'],
            'state_model': reference_model,
        })
    for original_level in settings['gap_levels']:
        level = _resolve_level(original_level, active_space)
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
                    'basis': basis, 'state_model': state_model_for(basis, multiplicity),
                    **common,
                })
    return components


def expected_levels(parameters, *, symbols=None, charge=0, settings=None,
                    state_roots=None):
    """Return the immutable level declaration stored with a correction."""
    active_space = active_space_settings(parameters, symbols, charge)
    sampling = sampling_settings(parameters, active_space=active_space)
    if correction_model(parameters) == 'direct':
        trusted = trusted_settings(parameters, active_space=active_space)
        multiplicities = [int(parameters.get('mult', 1) or 1)]
    else:
        settings = settings or spin_splitting_settings(parameters)
        multiplicities = [settings['target_multiplicity'],
                          settings['reference_multiplicity']]
        trusted = {
            'model': 'spin_splitting',
            'target_multiplicity': settings['target_multiplicity'],
            'reference_multiplicity': settings['reference_multiplicity'],
            'fragment_multiplicities': settings['fragment_multiplicities'],
            'reference': _resolve_level(settings['reference'], active_space),
            'spin_gaps': [_resolve_level(level, active_space)
                          for level in settings['gap_levels']],
            'weighting': settings['weighting'],
        }
    external = parameters.get('vrc_tst_state_model_file', '')
    inline_basis_model = (parameters.get('vrc_tst_state_model') or {}).get('schema') in (2, 3)
    state_models = None
    if active_space is not None and (inline_basis_model or (isinstance(external, str) and external.strip())):
        bases = (list(sampling['basis']) if isinstance(sampling['basis'], list)
                 else [sampling['basis']])
        if correction_model(parameters) == 'direct':
            bases.extend(trusted['basis'] if isinstance(trusted['basis'], list) else [trusted['basis']])
        else:
            reference_basis = trusted['reference']['basis']
            bases.extend(reference_basis if isinstance(reference_basis, list)
                         else [reference_basis])
            for level in trusted['spin_gaps']:
                level_basis = level['basis']
                bases.extend(level_basis if isinstance(level_basis, list)
                             else [level_basis])
        state_models = {
            basis: resolved_state_model(
                parameters, symbols, charge, settings,
                roots=state_roots, basis=basis)
            for basis in dict.fromkeys(bases)
        }
        state_model = None
    else:
        state_model = (deepcopy(parameters.get('vrc_tst_state_model'))
                       if active_space is None
                       else resolved_state_model(
                           parameters, symbols, charge, settings,
                           roots=state_roots))
    result = {
        'sampling': sampling, 'trusted_correction': trusted,
        'active_space': active_space, 'state_model': state_model,
        'active_space_feasibility': active_space_feasibility(
            active_space, multiplicities),
    }
    if state_models is not None:
        result['state_models_by_basis'] = state_models
        payload = parameters.get('vrc_tst_state_model') or {}
        if external:
            path = Path(external).expanduser()
            if not path.is_absolute():
                path = Path(parameters.get('input_file') or Path.cwd()/'input.json').resolve().parent/path
            payload = json.loads(path.read_text())
        if payload.get('schema') == 3:
            result['state_models_by_spin'] = {
                str(mult): {basis: resolved_state_model(parameters, symbols, charge, settings,
                    roots=state_roots, basis=basis, multiplicity=mult) for basis in state_models}
                for mult in multiplicities}
    return result


def trusted_settings(parameters, *, active_space=None):
    return sampling_settings({
        'vrc_tst_sample_method': parameters['vrc_tst_high_method'],
        'vrc_tst_sample_basis': parameters['vrc_tst_high_basis'],
        'vrc_tst_sample_extrapolation_power': parameters.get(
            'vrc_tst_high_extrapolation_power', ORIGINAL_ANL_POWER)}, active_space=active_space)


def assemble_trusted_reference(distances, energies, parameters, *, active_space=None):
    renamed = {key.replace('trusted', 'sampling', 1): value
               for key, value in energies.items() if key.startswith('trusted')}
    return assemble_sampling_reference(distances, renamed,
        trusted_settings(parameters, active_space=active_space))


def sampling_settings(parameters, *, active_space=None):
    """Return the frozen sampling level and optional ANL-style CBS pair."""
    method = parameters.get('vrc_tst_sample_method')
    basis = parameters.get('vrc_tst_sample_basis')
    if not isinstance(method, str) or not method.strip():
        raise ValueError('vrc_tst_sample_method must be a nonempty string.')
    if isinstance(basis, str):
        if not basis.strip():
            raise ValueError('vrc_tst_sample_basis cannot be empty.')
        return _resolve_level({'method': method, 'basis': basis}, active_space)
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
    return _resolve_level({'method': method, 'basis': list(basis),
                           'extrapolation_power': power}, active_space)


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
    if distances.ndim != 1 or not len(distances) or not np.isfinite(distances).all():
        raise ValueError('Spin-splitting distances must be a nonempty finite vector.')
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
                'upper_basis_gap_hartree': (np.asarray(target_cbs['upper_hartree'])
                    - np.asarray(reference_cbs['upper_hartree'])).tolist(),
                'cbs_minus_upper_gap_hartree': (target - reference_state
                    - np.asarray(target_cbs['upper_hartree'])
                    + np.asarray(reference_cbs['upper_hartree'])).tolist(),
                'interpretation': 'two-point basis estimate; not an uncertainty bound',
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
        'method_sensitivity': {
            'energy_unit': 'hartree',
            'asymptote_reference_distance_angstrom': float(distances[-1]),
            'absolute_spread': np.ptp(np.asarray(list(reconstructed.values())), axis=0).tolist(),
            'asymptote_referenced_spread': np.ptp(np.asarray([
                np.asarray(curve) - curve[-1] for curve in reconstructed.values()]), axis=0).tolist(),
            'interpretation': ('method disagreement diagnostic, not a statistical error bar'
                               if len(reconstructed) > 1 else
                               'one method only: zero spread supplies no uncertainty estimate'),
        },
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
