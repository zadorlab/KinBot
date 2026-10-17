import os
import sys
import pickle
import shutil

import numpy as np
from ase import Atoms
from ase.vibrations import Vibrations
from ase.optimize import BFGS
from ase.io import read, write
from sella import Sella

from kinbot.constants import EVtoHARTREE
from kinbot.ase_modules.calculators.{code} import {Code}
from kinbot.stationary_pt import StationaryPoint
from kinbot.frequencies import calc_vibrations, get_frequencies, frequency_validation_error

scratch_dir = os.getcwd()

mol = Atoms(symbols={atom}, 
            positions={geom})

kwargs = {kwargs}
if '{Code}' == 'Gaussian':
    # Also support a live older controller that omitted a conformer checkpoint.
    kwargs.setdefault('chk', os.path.basename('{label}'))
if '{Code}' == 'ORCA':
    from kinbot.ase_modules.calculators.orca import OrcaProfile
    kwargs['profile'] = OrcaProfile(command=kwargs['profile'])

mol.calc = {Code}(**kwargs)
if '{Code}' == 'Gaussian':
    mol.get_potential_energy()
    kwargs['guess'] = 'Read'
    mol.calc = {Code}(**kwargs)

basename = os.path.basename('{label}')
frequency_mode = '{frequency_mode}'
pkl_file = '{label}.pkl'
if os.path.isfile(pkl_file):
    os.remove(pkl_file)


def publish(data):
    """Atomically hand one worker result to the KinBot driver."""
    payload = {{'sym': mol.symbols, 'pos': mol.positions,
               'calc': '{code}', 'name': '{label}', 'data': data}}
    temporary = f'{{pkl_file}}.tmp.{{os.getpid()}}'
    with open(temporary, 'wb') as stream:
        pickle.dump(payload, stream)
    os.replace(temporary, pkl_file)
    # Human-readable status for the monitor, including isolated native checks
    # with no driver to ingest the pickle. The pickle remains the QC handoff.
    from kinbot.execution import atomic_json
    atomic_json('{label}.result.json', {{
        'status': data['status'], 'error': data.get('error'),
        'minimum_recovery_attempts': data.get('minimum_recovery_attempts', 0)}})

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')


def report_failure(reason):
    message = 'Worker failure: ' + reason
    print(message, file=sys.stderr, flush=True)
    with open('{label}_sella.log', 'a') as stream:
        stream.write(message + '\n')
    return reason

# For monoatomic wells, just calculate the energy and exit. 
if len(mol) == 1:
    e = mol.get_potential_energy()
    data = {{'energy': e, 'frequencies': np.array([]), 'zpe': 0.0,
            'hess': np.zeros([3, 3]), 'status': 'normal'}}

    if os.path.isdir(f'{{basename}}'):
        shutil.rmtree(f'{{basename}}')

    with open('{label}_sella.log', 'a') as f:
        f.write('Sella optimization is not needed for atoms.\ndone\n')
else:
    data = {{'status': 'error'}}
    order = {order}
    sella_kwargs = {sella_kwargs}
    from kinbot.minimum_recovery import pop_worker_controls
    minimum_options = pop_worker_controls(sella_kwargs)
    imagfreq_threshold = minimum_options['imagfreq_threshold']
    sella_kwargs.setdefault('internal', True)
    if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
        sella_kwargs['internal'] = False

    recovery_limit = (minimum_options['minimum_recovery_attempts'] if order == 0
                      and frequency_mode == 'native_hessian'
                      and '{Code}' == 'Gaussian' else 0)
    recovery_count = 0
    for minimum_attempt in range(recovery_limit + 1):
        if len(mol.symbols) > 2:
            opt = Sella(mol,
                        order=order,
                        trajectory='{label}.traj',
                        logfile='{label}_sella.log',
                        **sella_kwargs)
        else:
            opt = BFGS(mol,
                       trajectory='{label}.traj',
                       logfile='{label}_sella.log')
        freqs = []
        mol.calc.label = '{label}'
        rejection = None
        converged = False
        failure_reason = ''
        try:
            converged = opt.run(fmax={fmax}, steps={steps})
            traj = read('{label}.traj', index=':')
            write('{label}.xyz', traj, format='xyz')
            if not converged:
                failure_reason = report_failure('Optimization did not converge within {steps} steps.')
        except Exception as error:
            failure_reason = report_failure(
                f'Optimization failed: {{type(error).__name__}}: {{error}}')
            converged = False
        if converged:
            try:
                if frequency_mode == 'native_hessian' and '{Code}' == 'Gaussian':
                    # Sella owns the optimization. Run one fixed-geometry native
                    # Gaussian frequency calculation afterward; no Gaussian
                    # optimization keyword is present.
                    from kinbot import reader_gauss
                    from kinbot.ase_modules.calculators.gaussian import \
                        fixed_frequency_parameters
                    from kinbot.utils import iowait

                    frequency_kwargs = fixed_frequency_parameters(kwargs)
                    frequency_kwargs['label'] = '{label}'
                    # Rendered checkpoint option: {{'chk': '{checkpoint}'}}
                    frequency_kwargs['chk'] = os.path.basename('{label}')
                    mol.calc = Gaussian(**frequency_kwargs)
                    e = mol.get_potential_energy()
                    iowait('{label}.log', 'gauss')
                    freqs = reader_gauss.read_freq('{label}.log', {atom})
                    zpe = reader_gauss.read_zpe('{label}.log')
                    hessian = None
                else:
                    freqs, zpe, hessian = calc_vibrations(mol, f'{{basename}}')
                rejection = frequency_validation_error(freqs, order, imagfreq_threshold)
                if rejection:
                    failure_reason = report_failure(rejection)
                    converged = False
                else:
                    e = mol.get_potential_energy()
                    if not np.isfinite(e) or not np.isfinite(zpe):
                        raise ValueError('Nonfinite energy or zero-point energy.')
                    data = {{'energy': e, 'frequencies': freqs, 'zpe': zpe,
                            'status': 'normal'}}
                    if hessian is not None:
                        data['hess'] = hessian
                    pass
            except Exception as error:
                failure_reason = report_failure(
                    f'Frequency evaluation failed: {{type(error).__name__}}: {{error}}')
                converged = False


        if converged or not rejection or minimum_attempt >= recovery_limit:
            break
        # A parser/optimizer failure is not grounds for a geometry change.
        # Only a finite, significantly imaginary minimum spectrum is eligible.
        if (freqs is None or not np.isfinite(freqs).all()
                or not any(value < -imagfreq_threshold for value in freqs)):
            break
        try:
            from kinbot.minimum_recovery import gaussian_minimum_restart
            report = gaussian_minimum_restart(
                mol, kwargs, '{label}', minimum_attempt + 1,
                minimum_options['minimum_recovery_displacement'] * 2 ** minimum_attempt)
            recovery_count += 1
            message = (f'Minimum recovery {{recovery_count}}/{{recovery_limit}}: '
                       f'selected {{report["selected"]}} displacement; '
                       f'energy change {{report["energy_change_eV"]:.8g}} eV. '
                       'Repeating optimization and Hessian at unchanged theory. '
                       f'Evidence: {{report["archive"]}}')
            print(message, flush=True)
            with open('{label}_sella.log', 'a') as stream:
                stream.write(message + '\n')
        except Exception as error:
            failure_reason = report_failure(
                f'Minimum recovery failed: {{type(error).__name__}}: {{error}}')
            break

    os.chdir(scratch_dir)

    if not converged:
        data = {{'status': 'error', 'error': failure_reason, 'frequencies': freqs}}
    
    data['minimum_recovery_attempts'] = recovery_count

    if os.path.isdir(f'{{basename}}'):
        shutil.rmtree(f'{{basename}}')

    if os.path.isdir(f'{{basename}}_vib'):
        shutil.rmtree(f'{{basename}}_vib')

    with open('{label}_sella.log', 'a') as f:
        f.write('done\n')

publish(data)
