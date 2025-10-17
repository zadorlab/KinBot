import os
import sys
import shutil

import numpy as np
from ase import Atoms
from ase.db import connect
from ase.vibrations import Vibrations
from ase.optimize import BFGS
from ase.io import read, write
from sella import Sella

from kinbot.constants import EVtoHARTREE
from kinbot.ase_modules.calculators.{code} import {Code}
from kinbot.stationary_pt import StationaryPoint
from kinbot.frequencies import calc_vibrations, get_frequencies

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

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

# For monoatomic wells, just calculate the energy and exit. 
if len(mol) == 1:
    e = mol.get_potential_energy()
    db.write(mol, name='{label}',
             data={{'energy': e, 'frequencies': np.array([]), 'zpe': 0.0,
                    'hess': np.zeros([3, 3]), 'status': 'normal'}})
    with open('{label}_sella.log', 'a') as f:
        f.write('done\n')
    sys.exit(0)

order = {order}
sella_kwargs = {sella_kwargs}
if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
    sella_kwargs['internal'] = False

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
attempts = 1
converged = False
while attempts <= 2:
    fmax = {fmax}
    steps = {steps}
    mol.calc.label = '{label}'
    try:
        converged = opt.run(fmax=fmax, steps=steps)
        traj = read('{label}.traj', index=':')
        write('{label}.xyz', traj, format='xyz')
        freqs, zpe, hessian = calc_vibrations(mol, '{label}')
        if order == 0 and (np.count_nonzero(np.array(freqs) < 0) > 1
                       or np.count_nonzero(np.array(freqs) < -50) >= 1):
            print(f'Found one or more imaginary frequencies. {{freqs[1:6]}}')
            converged = False
        elif order == 1 and (np.count_nonzero(np.array(freqs) < 0) > 2  # More than two imag frequencies
                         or np.count_nonzero(np.array(freqs) < -50) >= 2  # More than one imag frequency larger than 50i
                         or np.count_nonzero(np.array(freqs) < 0) == 0):  # No imaginary frequencies
            print(f'Wrong number of imaginary frequencies: {{freqs[6:]}}')
            converged = False
        else:
            converged = True
            e = mol.get_potential_energy()
            db.write(mol, name='{label}', 
                 data={{'energy': e, 'frequencies': freqs, 'zpe': zpe, 
                        'hess': hessian, 'status': 'normal'}})
            break
        if not converged:
            raise RuntimeError
    except:
        attempts += 1
        if len(mol.symbols) > 2 and attempts == 2:
            sella_kwargs['internal'] = 1 - sella_kwargs['internal']
            opt = Sella(mol,
                order=order,
                trajectory='{label}.traj',
                logfile='{label}_sella.log',
                **sella_kwargs)
        else:
            break

if not converged:
    data = {{'status': 'error'}}
    if freqs:
        data['frequencies'] = freqs
    db.write(mol, name='{label}', data=data)

if os.path.isdir('{label}'):
    shutil.rmtree('{label}')

if os.path.isdir('{label}_vib'):
    shutil.rmtree('{label}_vib')

with open('{label}_sella.log', 'a') as f:
    f.write('done\n')
