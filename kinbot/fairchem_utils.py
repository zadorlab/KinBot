"""Load FairChem models with KinBot's compilation setting."""

import os


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
    if os.path.basename(str(model)).lower().startswith('uma-'):
        return True
    predictor = getattr(calculator, 'predictor', None)
    module = getattr(getattr(predictor, 'model', None), 'module', None)
    backbone = getattr(module, 'backbone', None)
    return type(backbone).__module__.startswith('fairchem.core.models.uma.')


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


def check_h2_cation_cache(qc, job, directory='.'):
    """Do not reuse the old, uncorrected UMA H2+ minimum on restart."""
    if qc.par.get('fc_task_name') != 'omol':
        return
    # Check the saved record even when its model weights are unavailable.
    # A missing checkpoint cannot establish that an old H2+ result is non-UMA.
    row = next(qc.db.select(name=job, sort='-id', limit=1), None)
    if (row is None or list(row.symbols) != ['H', 'H']
            or row.data.get('h2_cation_reference')
            or row.data.get('uma_model') is False):
        return
    state = dict(row.data)
    charge = state.get('charge', state.get('input_charge'))
    mult = state.get('multiplicity', state.get('input_multiplicity'))
    sources = [job]
    if row.data.get('copied_from_job'):
        sources.append(row.data['copied_from_job'])
    for source in sources:
        if charge is not None and mult is not None:
            break
        reference = next(qc.db.select(name='stereochemistry/' + source, sort='-id', limit=1), None)
        if reference is not None:
            charge = reference.data.get('input_charge') if charge is None else charge
            mult = reference.data.get('input_multiplicity') if mult is None else mult
        # Old conformer rows lack state metadata. Read the literal worker
        # input, including the source of a published _low result, without
        # executing it; its filename does not encode charge.
        import ast
        try:
            with open(os.path.join(directory, source + '.py')) as handle:
                nodes = ast.parse(handle.read()).body
            for node in nodes:
                if (isinstance(node, ast.Assign)
                        and any(isinstance(target, ast.Name) and target.id == 'kwargs'
                                for target in node.targets)):
                    kwargs = ast.literal_eval(node.value)
                    charge = kwargs.get('charge') if charge is None else charge
                    mult = kwargs.get('mult') if mult is None else mult
                    break
        except (OSError, SyntaxError, ValueError, TypeError):
            pass
    if ((charge is not None and charge != 1)
            or (mult is not None and mult != 2)):
        return
    raise ValueError(f'{os.path.join(directory, job)}: saved FairChem H2 result lacks the corrected UMA H2+ '
                     'reference (or its model/charge/spin cannot be verified). Start in '
                     'a fresh directory; the existing calculations were preserved.')
