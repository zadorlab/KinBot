"""Opt-in, per-level theory profiles for KinBot calculations.

This module only resolves configuration. Job routing is deliberately separate so
the legacy global ``qc`` setting keeps its meaning for existing inputs.
"""

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any, Mapping


ANL_AUTOMATIC_METHODS = frozenset(('ANL0', 'ANL0-F12'))

DEFAULT_ANL_SIZE_POLICY = {
    'full_anl0_f12': {
        'max_atoms': 16, 'max_heavy_atoms': 4, 'max_heavy_bonds': 4},
    'l3_geometry': {
        'max_atoms': 30, 'max_heavy_atoms': 8, 'max_heavy_bonds': 8},
    'l3_single_point': {
        'max_atoms': 50, 'max_heavy_atoms': 12, 'max_heavy_bonds': 12},
}


def automatic_cbh_anl_requested(parameters):
    """Return whether an ordinary input requests the production solver.

    A generated reference calculation retains the parent's composite theory
    profile so its L1/L2 routing is identical, but it writes a standalone
    thermochemistry model and must never recursively start another CBH graph.
    Existing handoff consumers likewise already have a completed graph.
    """
    method = str(parameters.get('composite_method', '')).upper()
    return (method in ANL_AUTOMATIC_METHODS
            and not parameters.get('anl_handoff_manifest')
            and not parameters.get('thermochemistry_model_output'))


def validate_anl_size_policy(policy):
    """Validate the monotone connectivity limits used by automatic ANL."""
    if (not isinstance(policy, dict)
            or set(policy) != set(DEFAULT_ANL_SIZE_POLICY)):
        raise ValueError(
            'anl_size_policy must define full_anl0_f12, l3_geometry, and '
            'l3_single_point.')
    normalized = {}
    previous = None
    keys = ('max_atoms', 'max_heavy_atoms', 'max_heavy_bonds')
    for tier in DEFAULT_ANL_SIZE_POLICY:
        entry = policy[tier]
        if not isinstance(entry, dict) or set(entry) != set(keys):
            raise ValueError(f'anl_size_policy.{tier} has invalid limits.')
        if any(isinstance(value, bool) or not isinstance(value, int)
               or value < 1 for value in entry.values()):
            raise ValueError(
                f'anl_size_policy.{tier} limits must be positive integers.')
        normalized[tier] = dict(entry)
        if previous is not None and any(
                normalized[tier][key] < previous[key] for key in keys):
            raise ValueError(
                'ANL size tiers must be monotonically nondecreasing.')
        previous = normalized[tier]
    return normalized


@dataclass(frozen=True)
class TheoryProfile:
    name: str
    calculator: str
    method: str = ''
    basis: str = ''
    command: str = ''
    calculator_kwargs: Mapping[str, Any] = field(default_factory=dict)
    optimizer: str = 'sella'
    frequency_mode: str = 'auto'
    label: str = ''
    model_path: str = ''
    task_name: str = 'omol'
    device: str = 'cpu'

    @classmethod
    def from_dict(cls, name, values):
        if not isinstance(values, dict):
            raise ValueError(f'{name}_profile must be an object.')
        allowed = set(cls.__dataclass_fields__) - {'name'}
        unknown = set(values) - allowed
        if unknown:
            raise ValueError(f'Unknown {name}_profile options: {sorted(unknown)}')
        calculator = values.get('calculator')
        if not isinstance(calculator, str) or not calculator.strip():
            raise ValueError(f'{name}_profile requires a calculator.')
        kwargs = values.get('calculator_kwargs', {})
        if not isinstance(kwargs, dict):
            raise ValueError(f'{name}_profile calculator_kwargs must be an object.')
        if values.get('frequency_mode', 'auto') not in ('auto', 'ase_forces', 'native_hessian'):
            raise ValueError(f'Invalid {name}_profile frequency_mode.')
        if values.get('optimizer', 'sella') not in ('sella', 'bfgs'):
            raise ValueError(f'Invalid {name}_profile optimizer.')
        return cls(name=name, **deepcopy(values))


THEORY_PRESETS = {
    'uma-b2plyp-anl': {
        'l1': {
            'calculator': 'fairchem',
            # A model name is portable across installations. The FairChem
            # loader resolves an existing Hugging Face cache snapshot before
            # attempting network access, so users never paste cache paths.
            'model_path': 'uma-s-1p2',
            'task_name': 'omol',
            'frequency_mode': 'ase_forces',
        },
        'l2': {
            'calculator': 'gaussian',
            'method': 'B2PLYP',
            'basis': 'cc-pVTZ',
            'calculator_kwargs': {
                'EmpiricalDispersion': 'GD3BJ',
                'Symm': 'None',
                'scf': 'xqc',
                'integral': 'UltraFine',
            },
            'optimizer': 'sella',
            'frequency_mode': 'native_hessian',
            'label': 'B2PLYP-D3(BJ)/cc-pVTZ',
        },
    },
    'uma-b3lyp-anl-low': {
        'l1': {
            'calculator': 'fairchem',
            'model_path': 'uma-s-1p2',
            'task_name': 'omol',
            'frequency_mode': 'ase_forces',
        },
        'l2': {
            'calculator': 'gaussian',
            'method': 'B3LYP',
            'basis': 'cc-pVTZ',
            'calculator_kwargs': {
                'Symm': 'None',
                'scf': 'xqc',
                'integral': 'UltraFine',
            },
            'optimizer': 'sella',
            'frequency_mode': 'native_hessian',
            'label': 'B3LYP/cc-pVTZ',
        },
    },
}

_COMPOSITE_DEFAULT_PRESET = {
    'ANL0': 'uma-b3lyp-anl-low',
    'ANL0-F12': 'uma-b3lyp-anl-low',
    'ANL1': 'uma-b2plyp-anl',
    # ANL1-F12 is a ladder request used by the CBH workflow. Its exact
    # electronic expression remains a separately provenance-labelled recipe;
    # L1/L2 nevertheless use the higher B2PLYP surface.
    'ANL1-F12': 'uma-b2plyp-anl',
    'ANL1-QZF': 'uma-b2plyp-anl',
}

_ALIASES = {'uma': 'fairchem', 'fc': 'fairchem', 'gauss': 'gaussian'}


def profiled_requested(par):
    """Return whether an explicit input opted in to per-level routing."""
    return bool(par.get('profiled_theory') or par.get('theory_preset')
                or par.get('l1') or par.get('l2')
                or par.get('l1_profile') or par.get('l2_profile')
                or par.get('composite_method') or par.get('l3_overrides')
                or par.get('l3_resource_overrides'))


def _merge(base, override):
    merged = deepcopy(base)
    for key, value in override.items():
        if key == 'calculator_kwargs' and isinstance(value, dict):
            merged[key] = {**merged.get(key, {}), **deepcopy(value)}
        else:
            merged[key] = deepcopy(value)
    return merged


def _legacy_profile(par, level):
    calculator = _ALIASES.get(par['qc'].lower(), par['qc'].lower())
    method_key, basis_key = (('method', 'basis') if level == 'l1'
                             else ('high_level_method', 'high_level_basis'))
    return {
        'calculator': calculator,
        'method': par[method_key],
        'basis': par[basis_key],
        'command': par['qc_command'],
    }


def resolve_profiles(par):
    """Resolve L1/L2 with precedence: legacy < preset < explicit profiles.

    No external program or optional package is loaded during resolution.
    """
    if not profiled_requested(par):
        return {}

    for key in ('l3_overrides', 'l3_resource_overrides'):
        if not isinstance(par.get(key, {}), dict):
            raise ValueError(f'{key} must be an object.')
        if par.get(key) and not par.get('composite_method'):
            raise ValueError(f'{key} requires composite_method.')

    preset_name = par.get('theory_preset')
    if not preset_name and par.get('composite_method'):
        method = par['composite_method'].upper()
        if method not in _COMPOSITE_DEFAULT_PRESET:
            raise ValueError(f'{method}: select an explicit theory_preset for '
                             'this composite method.')
        preset_name = _COMPOSITE_DEFAULT_PRESET[method]
    if preset_name and preset_name not in THEORY_PRESETS:
        raise ValueError(f'Unknown theory_preset: {preset_name}')
    preset = THEORY_PRESETS.get(preset_name, {})

    profiles = {}
    for level in ('l1', 'l2'):
        legacy = _legacy_profile(par, level)
        values = deepcopy(legacy)
        preset_values = preset.get(level, {})
        if (preset_values.get('calculator')
                and preset_values['calculator'] != values['calculator']):
            values['command'] = ''
        values = _merge(values, preset_values)
        alias = par.get(level)
        if alias:
            if not isinstance(alias, str):
                raise ValueError(f'{level} must be a calculator alias.')
            calculator = _ALIASES.get(alias.lower(), alias.lower())
            if calculator != values['calculator']:
                # A backend switch must not retain a different backend's
                # method, basis, frequency policy, or private keywords.
                values = deepcopy(legacy)
                values['command'] = ''
            values['calculator'] = calculator
        override = par.get(f'{level}_profile') or {}
        if not isinstance(override, dict):
            raise ValueError(f'{level}_profile must be an object.')
        if 'calculator' in override and not isinstance(override['calculator'], str):
            raise ValueError(f'{level}_profile calculator must be a string.')
        if ('calculator' in override
                and _ALIASES.get(override['calculator'].lower(),
                                 override['calculator'].lower()) != values['calculator']):
            values = deepcopy(legacy)
            values['command'] = ''
        values = _merge(values, override)
        values['calculator'] = _ALIASES.get(values['calculator'].lower(),
                                            values['calculator'].lower())
        if values['calculator'] == 'fairchem':
            # Explicit profile > legacy/global option > portable preset.
            # Keeping this precedence preserves existing site-local inputs
            # while allowing a hands-off named model in the preset.
            values['model_path'] = (
                override.get('model_path')
                or par.get('fc_model_path', '')
                or values.get('model_path', ''))
            values['task_name'] = values.get('task_name') or par.get('fc_task_name', 'omol')
            values['device'] = values.get('device') or par.get('fc_device', 'cpu')
            if not values['model_path']:
                raise ValueError(f'{level}_profile requires a FairChem model_path '
                                 'or fc_model_path.')
        profiles[level] = TheoryProfile.from_dict(level, values)
    return profiles
