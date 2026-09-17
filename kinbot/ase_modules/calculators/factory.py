"""Calculator registry for the opt-in profiled job path.

Imports are lazy so legacy KinBot and users without FairChem can import this
module without installing every possible calculator backend.
"""

from dataclasses import dataclass
from importlib import import_module
from pathlib import Path
import re

from kinbot.theory import TheoryProfile


@dataclass(frozen=True)
class CalculatorCapabilities:
    energy: bool = True
    forces: bool = False
    native_hessian: bool = False
    dipole: bool = False
    rohf: bool = False
    uhf: bool = False
    f12: bool = False
    numerical_forces: bool = False


@dataclass(frozen=True)
class CalculatorSpec:
    module: str
    class_name: str
    capabilities: CalculatorCapabilities


CALCULATORS = {
    'gaussian': CalculatorSpec(
        'kinbot.ase_modules.calculators.gaussian', 'Gaussian',
        CalculatorCapabilities(forces=True, dipole=True)),
    'qchem': CalculatorSpec(
        'kinbot.ase_modules.calculators.qchem', 'QChem',
        CalculatorCapabilities(forces=True)),
    'orca': CalculatorSpec(
        'kinbot.ase_modules.calculators.orca', 'ORCA',
        CalculatorCapabilities(forces=True, dipole=True)),
    'fairchem': CalculatorSpec(
        'fairchem.core', 'FAIRChemCalculator',
        CalculatorCapabilities(forces=True)),
}

_ALIASES = {'gauss': 'gaussian', 'fc': 'fairchem'}


def calculator_spec(name):
    name = _ALIASES.get(name.lower(), name.lower())
    try:
        return CALCULATORS[name]
    except KeyError as exc:
        raise ValueError(f'No profiled ASE calculator registered for {name!r}.') from exc


def capabilities(name):
    return calculator_spec(name).capabilities


def calculator_class(name):
    spec = calculator_spec(name)
    return getattr(import_module(spec.module), spec.class_name)


def build_calculator(profile, directory, task):
    """Construct a calculator for a job without changing global backend state.

    ``task`` is a mapping with ``name``, ``charge`` and ``mult``. The latter two
    default to a neutral singlet for simple calculator-layer tests.
    """
    if isinstance(profile, dict):
        profile = TheoryProfile.from_dict('job', profile)
    if not isinstance(profile, TheoryProfile):
        raise TypeError('profile must be a TheoryProfile or profile object.')
    if not isinstance(task, dict):
        raise TypeError('task must be an object with a name.')
    name = task.get('name')
    if (not isinstance(name, str)
            or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name)):
        raise ValueError('task.name must be a safe nonempty basename.')
    backend = _ALIASES.get(profile.calculator.lower(), profile.calculator.lower())
    calculator_spec(backend)
    if backend == 'fairchem' and not profile.model_path:
        raise ValueError('FairChem profile requires model_path.')
    if backend == 'orca' and not profile.command:
        raise ValueError('ORCA profile requires a command.')
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    calc_class = calculator_class(backend)
    kwargs = dict(profile.calculator_kwargs)
    charge = task.get('charge', 0)
    mult = task.get('mult', 1)

    if backend == 'gaussian':
        kwargs.update(method=profile.method, basis=profile.basis,
                      charge=charge, mult=mult)
        calc = calc_class(label=str(directory / name), **kwargs)
        if profile.command:
            calc.command = f'{profile.command} < PREFIX.com > PREFIX.log'
        return calc
    if backend == 'qchem':
        kwargs.update(method=profile.method, basis=profile.basis,
                      charge=charge, multiplicity=mult)
        calc = calc_class(label=str(directory / name), **kwargs)
        if profile.command:
            calc.command = f'{profile.command} PREFIX.inp PREFIX.out'
        return calc
    if backend == 'orca':
        from kinbot.ase_modules.calculators.orca import OrcaProfile
        kwargs.setdefault('orcasimpleinput', f'{profile.method} {profile.basis}'.strip())
        kwargs.update(charge=charge, mult=mult)
        return calc_class(profile=OrcaProfile(command=profile.command),
                          directory=directory / name, **kwargs)
    if backend == 'fairchem':
        from kinbot.fairchem_utils import load_predictor
        return calc_class(load_predictor(profile.model_path, profile.device),
                          task_name=profile.task_name, **kwargs)
    raise AssertionError(f'Unhandled calculator {backend}')
