"""Opt-in, per-level theory profiles for KinBot calculations.

This module only resolves configuration. Job routing is deliberately separate so
the legacy global ``qc`` setting keeps its meaning for existing inputs.
"""

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any, Mapping


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
            values['model_path'] = values.get('model_path') or par.get('fc_model_path', '')
            values['task_name'] = values.get('task_name') or par.get('fc_task_name', 'omol')
            values['device'] = values.get('device') or par.get('fc_device', 'cpu')
            if not values['model_path']:
                raise ValueError(f'{level}_profile requires a FairChem model_path '
                                 'or fc_model_path.')
        profiles[level] = TheoryProfile.from_dict(level, values)
    return profiles
