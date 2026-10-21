from kinbot.optimization_budget import run_optimizer
import os
import shutil
import sys

import numpy as np
from ase import Atoms
from ase.db import connect
from ase.io import read, write
from ase.vibrations import Vibrations
from sella import Sella

from kinbot.constants import EVtoHARTREE
from kinbot.ase_modules.calculators.{code} import {Code}
from kinbot.stationary_pt import StationaryPoint
from kinbot.frequencies import calc_vibrations, get_frequencies, frequency_validation_error

db = connect('{working_dir}/kinbot.db')
mol = Atoms(symbols={atom},
            positions={geom})

kwargs = {kwargs}
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

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')


def report_failure(reason):
    message = 'Worker failure: ' + reason
    print(message, file=sys.stderr, flush=True)
    with open('{label}_sella.log', 'a') as stream:
        stream.write(message + '\n')
    return reason

sella_kwargs = {sella_kwargs}
from kinbot.minimum_recovery import pop_worker_controls
imagfreq_threshold = pop_worker_controls(sella_kwargs)['imagfreq_threshold']
sella_kwargs.setdefault('internal', True)
if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
    sella_kwargs['internal'] = False
opt = Sella(mol, order=1, 
            trajectory='{label}.traj',
            logfile='{label}_sella.log',
            **sella_kwargs)
freqs = []

fmax = {fmax}
steps = {steps}
mol.calc.label = '{label}'

converged = False
failure_reason = ''
try:
    converged = run_optimizer(opt, fmax=fmax, steps=steps, report='{label}.optimization_budget.json')
    traj = read('{label}.traj', index=':')
    write('{label}.xyz', traj, format='xyz')
    if not converged:
        failure_reason = report_failure(f'Optimization did not converge within {{opt.nsteps}} steps.')
except Exception as error:
    failure_reason = report_failure(
        f'Optimization failed: {{type(error).__name__}}: {{error}}')
if converged:
    try:
        if frequency_mode == 'native_hessian' and '{Code}' == 'Gaussian':
            # Sella owns the TS optimization. Evaluate the native Gaussian
            # Hessian once at the fixed, optimized saddle geometry.
            from kinbot import reader_gauss
            from kinbot.ase_modules.calculators.gaussian import \
                fixed_frequency_parameters
            from kinbot.utils import iowait

            frequency_kwargs = fixed_frequency_parameters(kwargs)
            frequency_kwargs['label'] = '{label}'
            frequency_kwargs['chk'] = os.path.basename('{label}')
            mol.calc = Gaussian(**frequency_kwargs)
            e = mol.get_potential_energy()
            iowait('{label}.log', 'gauss')
            freqs = reader_gauss.read_freq('{label}.log', {atom})
            zpe = reader_gauss.read_zpe('{label}.log')
            hessian = None
        else:
            freqs, zpe, hessian = calc_vibrations(mol, f'{{basename}}')
        rejection = frequency_validation_error(freqs, 1, imagfreq_threshold)
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
            db.write(mol, name='{label}', data=data)
    except Exception as error:
        failure_reason = report_failure(
            f'Frequency evaluation failed: {{type(error).__name__}}: {{error}}')
        converged = False

if 'vib' in os.getcwd():
    os.chdir("..")

if not converged:
    data = {{'status': 'error', 'error': failure_reason}}
    data['frequencies'] = freqs
    db.write(mol, name='{label}', data=data)

if os.path.isdir(f'{{basename}}'):
    shutil.rmtree(f'{{basename}}')

if os.path.isdir(f'{{basename}}_vib'):
    if os.path.exists(os.path.join(f'{{basename}}_vib', 'vib.xyz')):
        shutil.copy2(os.path.join(f'{{basename}}_vib', 'vib.xyz'), 
                     os.path.join(os.getcwd(), f'{{basename}}_vib.xyz'))
    shutil.rmtree(f'{{basename}}_vib')

with open('{label}_sella.log', 'a') as f:
    f.write('done\n')
