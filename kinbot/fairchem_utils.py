"""Load FairChem models with KinBot's compilation setting."""

import ast
import os
import re


def load_predictor(model, device, *, compile_model=False):
    """Load a local checkpoint or a named pretrained model."""
    from fairchem.core.units.mlip_unit.api.inference import inference_settings_default

    # FairChem 2.22 enabled compilation by default. Short KinBot jobs often
    # finish before faster force evaluations repay that startup cost.
    settings = inference_settings_default()
    settings.compile = compile_model
    if compile_model:
        # FairChem 2.14 can fail in compiled UMA force evaluations without
        # precomputed composition weights. A KinBot job keeps its composition,
        # charge and spin fixed, as required by this preparation.
        settings.merge_mole = True
    if os.path.isfile(model):
        from fairchem.core.units.mlip_unit import load_predict_unit
        return load_predict_unit(model, device=device, inference_settings=settings)
    from fairchem.core import pretrained_mlip
    return pretrained_mlip.get_predict_unit(model, device=device,
                                          inference_settings=settings)


def uma_model(model, calculator=None):
    """Recognize named UMA models and local checkpoints by their backbone."""
    predictor = getattr(calculator, 'predictor', None)
    module = getattr(getattr(predictor, 'model', None), 'module', None)
    backbone = getattr(module, 'backbone', None)
    if backbone is not None:
        return type(backbone).__module__.startswith('fairchem.core.models.uma.')
    # Only names passed to the pretrained-model registry identify the model
    # family. Users can give a local checkpoint any filename.
    return (not os.path.isfile(model) and os.path.basename(str(model)) == str(model)
            and str(model).lower().startswith('uma-'))


def h2_cation_reference(mol, label, model, fmax, steps):
    """Replace only an isolated UMA H2+ minimum, using neutral H2 as reference.

    The measured adiabatic ionization energy includes ZPE. Subtract the fixed
    ion ZPE when storing its electronic energy so KinBot does not add it twice.
    The neutral reference uses this worker's UMA calculator and molecular task.
    """
    if (mol.get_chemical_symbols() != ['H', 'H']
            or mol.info.get('charge') != 1 or mol.info.get('spin') != 2
            or getattr(mol.calc, 'task_name', None) != 'omol'
            or not uma_model(model, mol.calc)):
        return None

    import shutil
    import numpy as np
    from ase import Atoms
    from ase.optimize import BFGS
    from kinbot import constants
    from kinbot.frequencies import calc_vibrations

    # NIST: https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Mask=20
    # Harmonic stretch: https://cccbdb.nist.gov/expdiatomicsx.asp
    ionization_ev, stretch_cm, distance = 15.42593, 2321.7, 1.058
    neutral = Atoms('H2', positions=[[0., 0., 0.], [.74, 0., 0.]])
    neutral.info.update(charge=0, spin=1)
    neutral.calc = mol.calc
    reference_label = label + '_neutral_h2'
    opt = BFGS(neutral, logfile=reference_label + '.log')
    if not opt.run(fmax=fmax, steps=steps):
        raise RuntimeError('Neutral H2 reference optimization did not converge.')
    energy = neutral.get_potential_energy()
    try:
        frequencies, zpe, _ = calc_vibrations(neutral, reference_label)
    finally:
        shutil.rmtree(reference_label + '_vib', ignore_errors=True)
    if (frequencies is None or len(frequencies) != 1
            or not np.isfinite(frequencies[0]) or frequencies[0] <= 0.
            or zpe is None or not np.isfinite(zpe) or not np.isfinite(energy)):
        raise RuntimeError('Neutral H2 reference has no valid minimum frequencies/ZPE.')

    axis = mol.positions[1] - mol.positions[0]
    axis /= np.linalg.norm(axis)
    center = mol.positions.mean(axis=0)
    mol.positions = center + np.array([-.5, .5])[:, None] * distance * axis
    # One longitudinal stretch and five zero rigid motions. Units match
    # calc_vibrations: unweighted Hessian in Hartree/bohr^2, ZPE in Hartree.
    conversion = 2. * np.pi * constants.SPEEDofLIGHT * constants.BOHRtoCM
    force_constant = .5 * constants.exact_mass['H'] * constants.MEtoAMU * (conversion * stretch_cm)**2
    block = force_constant * np.outer(axis, axis)
    ion_zpe = .5 * stretch_cm / constants.AUtoCM
    mol.calc = None
    return {'status': 'normal',
            'charge': 1, 'multiplicity': 2, 'uma_model': True,
            'energy': energy + (zpe - ion_zpe) / constants.EVtoHARTREE + ionization_ev,
            'frequencies': np.array([stretch_cm]), 'zpe': ion_zpe,
            'hess': np.block([[block, -block], [-block, block]]),
            'h2_cation_reference': {'ionization_energy_ev': ionization_ev,
                                    'neutral_h2_energy_ev': energy,
                                    'neutral_h2_zpe_hartree': zpe}}


class SavedHydrogenError(ValueError):
    """The saved electronic state or energy cannot be used for this hydrogen fragment."""


class UncorrectedH2Result(SavedHydrogenError):
    """A saved FairChem calculation may contain an obsolete H2+ energy."""


def saved_hydrogen_state(db, row, job, directory='.'):
    """Recover charge and multiplicity without guessing a charge from a species number."""
    sources = [job]
    if row.data.get('copied_from_job'):
        sources.append(row.data['copied_from_job'])
    # Old selected conformers often omit both charge and spin. Their ordinary
    # well has the same electronic state, even if no conformer worker survives.
    name_spin = None
    for source in tuple(sources):
        match = re.fullmatch(r'((?:1000000000000000000|2002000000000000000)([1-9])(?:-s[0-9a-f]{16})?)'
                             r'_(?:low|well(?:_high|_mp2)?|product|(?:semi_emp_)?[0-9]+(?:_high)?)',
                             os.path.basename(source))
        if match:
            name_spin = int(match[2])
            sources.append(match[1] + '_well')
            break
    sources = list(dict.fromkeys(sources))
    result_rows = [row]
    copied = row.data.get('copied_from_job')
    if copied:
        selector = {'name': copied}
        if row.data.get('copied_from_row_id') is not None:
            selector['id'] = row.data['copied_from_row_id']
        original = next(db.select(**selector, sort='-id', limit=1), None)
        if original is not None:
            result_rows.append(original)
    records = []
    for result in result_rows:
        records.extend((result.data, result.data.get('identity', {}),
                        getattr(result, 'calculator_parameters', {}) or {}))
        if result.data.get('proton_reference'):
            records.append({'charge': 1, 'mult': 1})

    def recorded_value(fields):
        values = {int(record[key]) for record in records for key in fields if record.get(key) is not None}
        if len(values) > 1:
            raise SavedHydrogenError(f'{os.path.join(directory, job)}: conflicting saved {fields[0]} '
                                     f'values {sorted(values)} for hydrogen. Use a fresh directory.')
        return next(iter(values), None)

    charge_fields = ('charge', 'input_charge')
    spin_fields = ('multiplicity', 'mult', 'input_multiplicity')
    need_worker = recorded_value(charge_fields) is None or recorded_value(spin_fields) is None
    for source in sources:
        reference = next(db.select(name='stereochemistry/' + source, sort='-id', limit=1), None)
        if reference is not None:
            records.append(reference.data)
        if not need_worker:
            continue
        try:
            with open(os.path.join(directory, source + '.py')) as handle:
                nodes = ast.parse(handle.read()).body
            for node in nodes:
                if (isinstance(node, ast.Assign)
                        and any(isinstance(target, ast.Name) and target.id == 'kwargs'
                                for target in node.targets)):
                    kwargs = ast.literal_eval(node.value)
                    if isinstance(kwargs, dict):
                        records.append(kwargs)
                    break
        except (OSError, SyntaxError, ValueError, TypeError):
            continue
    charge = recorded_value(charge_fields)
    mult = recorded_value(spin_fields)
    if list(row.symbols) == ['H']:
        if charge is None:
            # ASE can store an explicit initial charge without a calculator.
            # An absent array supplies no information about H+ versus H-.
            for result in result_rows:
                initial = result.get('initial_charges')
                if initial is not None and float(initial[0]).is_integer():
                    records.append({'charge': initial[0]})
            charge = recorded_value(charge_fields)
        if name_spin is not None and mult is not None and name_spin != mult:
            raise SavedHydrogenError(f'{os.path.join(directory, job)}: saved hydrogen multiplicity {mult} '
                                     f'differs from {name_spin} in the job name. Use a fresh directory.')
    return charge, name_spin if mult is None else mult


def check_h2_cation_cache(qc, job, directory='.', *, row=None, state=None):
    """Check saved hydrogen energies affected by earlier FairChem workers."""
    if qc.par.get('fc_task_name') != 'omol':
        return
    # Check the saved record even when its model weights are unavailable.
    if row is None:
        row = next(qc.db.select(name=job, sort='-id', limit=1), None)
    if row is None:
        return
    atomic_zero = (list(row.symbols) == ['H'] and row.data.get('energy') == 0.
                   and not row.data.get('atomic_energy_evaluated')
                   and not row.data.get('proton_reference'))
    if not atomic_zero and (list(row.symbols) != ['H', 'H'] or row.data.get('h2_cation_reference')
                            or row.data.get('uma_model') is False):
        return
    charge, mult = state if state is not None else saved_hydrogen_state(qc.db, row, job, directory)
    if atomic_zero:
        if charge == -1:
            raise SavedHydrogenError(f'{os.path.join(directory, job)}: saved FairChem hydride (H-) '
                                     'has the old zero proton energy. Recalculate in a fresh directory.')
        if charge is None and mult == 1:
            raise SavedHydrogenError(f'{os.path.join(directory, job)}: singlet hydrogen has no saved '
                                     'charge to distinguish H+ from H-. Use a fresh directory.')
        return
    if ((charge is not None and charge != 1)
            or (mult is not None and mult != 2)):
        return
    raise UncorrectedH2Result(f'{os.path.join(directory, job)}: saved FairChem H2 result lacks the corrected UMA H2+ '
                     'reference (or its model/charge/spin cannot be verified). Start in '
                     'a fresh directory; the existing calculations were preserved.')
