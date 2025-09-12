import os
import random

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.db import connect
from sella import Sella

from fairchem.core import pretrained_mlip, FAIRChemCalculator
#from kinbot.ase_modules.calculators.{code} import {Code}
from kinbot.frequencies import calc_vibrations

db = connect('{working_dir}/kinbot.db')

mol = Atoms(symbols={atom},
            positions={geom})

kwargs = {kwargs}
mol.info.update({{"charge": kwargs['charge'], "spin": kwargs['mult']}})

mol.calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit("uma-s-1", device="cpu"), task_name="omol")

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

sella_kwargs = {sella_kwargs}
if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
    sella_kwargs['internal'] = False
opt = Sella(mol, order=1, 
            trajectory='{label}.traj',
            logfile='{label}_sella.log',
            **sella_kwargs)
freqs = []
converged = False
attempts = 1
while attempts <= 2:
    fmax = 1e-4
    steps = 250
    mol.calc.label = '{label}'
    try:
        converged = opt.run(fmax=fmax, steps=steps)
        traj = read('{label}.traj', index=':')
        write('{label}.xyz', traj, format='xyz')
        freqs, zpe, hessian = calc_vibrations(mol, '{label}')
        if (np.count_nonzero(np.array(freqs) < 0) > 2  # More than two imag frequencies
            or np.count_nonzero(np.array(freqs) < -50) >= 2  # More than one frequency smaller than 50i
            or np.count_nonzero(np.array(freqs) < 0) == 0):  # No imaginary frequencies
            converged = False
            mol.calc.label = '{label}'
        else:
            converged = True
            e = mol.get_potential_energy()

            random.seed()
            db.write(mol, name='{label}', 
                 data={{'energy': e, 'frequencies': freqs, 'zpe': zpe, 
                     'hess': hessian, 'status': 'normal'}})            
            break
        if not converged:
            raise RuntimeError
    except:
        attempts += 1
        if attempts == 2:
            sella_kwargs['internal'] = 1 - sella_kwargs['internal']
            opt = Sella(mol, order=1,
                trajectory='{label}.traj',
                logfile='{label}_sella.log',
                **sella_kwargs)
        else:
            break
if not converged:
    data = {{'status': 'error'}}
    if freqs:
        data['frequencies'] = freqs
    random.seed()
    db.write(mol, name='{label}', data=data)

with open('{label}_sella.log', 'a') as f:
    f.write('done\n')
