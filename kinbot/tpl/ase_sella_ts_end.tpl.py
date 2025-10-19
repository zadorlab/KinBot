import os
import shutil

import numpy as np
from ase import Atoms
from ase.db import connect
from ase.io import read, write
from ase.vibrations import Vibrations
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

basename = os.path.basename('{label}')

if os.path.isfile(f'{{basename}}_sella.log'):
    os.remove(f'{{basename}}_sella.log')

sella_kwargs = {sella_kwargs}
if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
    sella_kwargs['internal'] = False
opt = Sella(mol, order=1, 
            trajectory=f'{{basename}}.traj',
            logfile=f'{{basename}}_sella.log',
            **sella_kwargs)
freqs = []

fmax = {fmax}
steps = {steps}
mol.calc.label = '{label}'

converged = opt.run(fmax=fmax, steps=steps)
traj = read(f'{{basename}}.traj', index=':')
write(f'{{basename}}.xyz', traj, format='xyz')
if converged:
    freqs, zpe, hessian = calc_vibrations(mol, f'{{basename}}')
    if (np.count_nonzero(np.array(freqs) < 0) > 2  # More than two imag frequencies
        or np.count_nonzero(np.array(freqs) < -50) >= 2  # More than one frequency smaller than 50i
        or np.count_nonzero(np.array(freqs) < 0) == 0):  # No imaginary frequencies

        data = {{'status': 'error'}}
        data['frequencies'] = freqs
        db.write(mol, name='{label}', data=data)

    else:
        e = mol.get_potential_energy()
        db.write(mol, name='{label}', 
             data={{'energy': e, 'frequencies': freqs, 'zpe': zpe, 
                    'hess': hessian, 'status': 'normal'}})            

if os.path.isdir(f'{{basename}}'):
    shutil.rmtree(f'{{basename}}')

if os.path.isdir(f'{{basename}}_vib'):
    shutil.copy2(os.path.join(f'{{basename}}_vib', f'vib.xyz'), os.getcwd())
    shutil.rmtree(f'{{basename}}_vib')

with open(f'{{basename}}_sella.log', 'a') as f:
    f.write('done\n')
