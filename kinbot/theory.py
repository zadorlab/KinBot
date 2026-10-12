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
        if values.get('optimizer', 'sella') not in ('sella', 'bfgs', 'native'):
            raise ValueError(f'Invalid {name}_profile optimizer.')
        return cls(name=name, **deepcopy(values))


_ALIASES = {'uma': 'fairchem', 'fc': 'fairchem', 'gauss': 'gaussian'}


def profiled_requested(par):
    """Return whether an explicit input opted in to per-level routing."""
    return bool(par.get('l0') or par.get('l0_profile')
                or par.get('l0_qc') or par.get('l0_qc_command')
                or par.get('l0_basis') or par.get('l0_calc_kwargs')
                or par.get('profiled_theory')
                or par.get('l1') or par.get('l2')
                or par.get('l1_profile') or par.get('l2_profile'))


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
    values = {
        'calculator': calculator,
        'method': par[method_key],
        'basis': par[basis_key],
        'command': par['qc_command'],
        'calculator_kwargs': deepcopy(par.get('calc_kwargs', {})),
        'optimizer': 'sella' if par['use_sella'] or calculator in ('fairchem', 'orca', 'nn_pes') else 'native',
    }
    if level == 'l0':
        # Preliminary searches keep the user's backend and command.
        l0_calculator = par.get('l0_qc') or par['qc']
        values.update(calculator=_ALIASES.get(l0_calculator.lower(), l0_calculator.lower()),
                      method=par['l0_method'], basis=par.get('l0_basis', ''),
                      calculator_kwargs=deepcopy(par.get('l0_calc_kwargs') or
                                                 par.get('calc_kwargs', {})))
        values['command'] = par.get('l0_qc_command') or (
            par['qc_command'] if values['calculator'] == calculator else '')
        if values['calculator'] == 'fairchem' and values['method'].lower() != 'am1':
            # Retain the old flat L0 model-name override.
            values['model_path'] = values['method']
    return values


def resolve_profiles(par):
    """Resolve L0/L1/L2 by applying explicit profiles over legacy settings.

    No external program or optional package is loaded during resolution.
    """
    if par.get('theory_preset', '') != '':
        raise ValueError('theory_preset is not supported in this branch. Specify '
                         'l0_profile, l1_profile and l2_profile explicitly; '
                         'a restart must keep the previously recorded methods and models.')
    if not profiled_requested(par):
        return {}

    profiles = {}
    for level in ('l0', 'l1', 'l2'):
        legacy = _legacy_profile(par, level)
        values = deepcopy(legacy)
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
                values['optimizer'] = 'native' if calculator == 'qchem' else 'sella'
                values['calculator_kwargs'] = {}
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
            values['optimizer'] = ('native' if override['calculator'].lower() == 'qchem'
                                   else 'sella')
            values['calculator_kwargs'] = {}
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
        if values['calculator'] == 'fairchem':
            values['method'], values['basis'] = '', ''
        profile = TheoryProfile.from_dict(level, values)
        if profile.calculator == 'nwchem':
            raise ValueError(f'{level}_profile: NWChem calculation workers are not '
                             'supported by the profile workflow. Select Gaussian, '
                             'FairChem, or Q-Chem with optimizer=native instead.')
        if profile.calculator not in ('gaussian', 'qchem', 'orca', 'nn_pes', 'fairchem'):
            raise ValueError(f'{level}_profile has unsupported calculator {profile.calculator!r}.')
        if profile.calculator == 'qchem' and profile.optimizer != 'native':
            raise ValueError(f'{level}_profile Q-Chem requires optimizer=native; '
                             'its Sella/BFGS workers do not support force-only calculations.')
        if profile.frequency_mode == 'native_hessian' and profile.calculator not in ('gaussian', 'qchem'):
            raise ValueError(f'{level}_profile native_hessian requires Gaussian or Q-Chem.')
        if profile.frequency_mode == 'ase_forces' and profile.optimizer == 'native':
            raise ValueError(f'{level}_profile ase_forces requires a Sella or BFGS optimizer.')
        profiles[level] = profile
    return profiles


def calculation_label(par):
    """Describe the selected L1 or L2 calculation in MESS output and PES logs."""
    profiles = resolve_profiles(par)
    if profiles:
        profile = profiles['l2' if par['high_level'] else 'l1']
        if profile.label:
            return profile.label
        if profile.calculator == 'fairchem':
            return f'FairChem {profile.model_path} ({profile.task_name})'
        label = f'{profile.method}/{profile.basis}'
        dispersion = next((value for key, value in profile.calculator_kwargs.items()
                           if key.lower() == 'empiricaldispersion'), None)
        if dispersion:
            label += f' (EmpiricalDispersion={dispersion})'
        return label
    if par['qc'] == 'fc':
        return f"FairChem {par['fc_model_path']} ({par['fc_task_name']})"
    if par['qc'] == 'nn_pes':
        return f"nn_pes {par['nn_model']}"
    prefix = 'high_level_' if par['high_level'] else ''
    return f"{par[prefix + 'method']}/{par[prefix + 'basis']}"
