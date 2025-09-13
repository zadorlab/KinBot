import os

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.db import connect
from sella import Sella

from fairchem.core import pretrained_mlip, FAIRChemCalculator
#from kinbot.ase_modules.calculators.{code} import {Code}
from kinbot.frequencies import calc_vibrations
from kinbot.utils import sella_freq_check

db = connect('{working_dir}/kinbot.db')
if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

# molecule
mol = Atoms(symbols={atom},
            positions={geom})
kwargs = {kwargs}
mol.info.update({{"charge": kwargs['charge'], "spin": kwargs['mult']}})
mol.calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit("uma-s-1", device="cpu"), task_name="omol")
mol.calc.label = '{label}'
freqs = []

# sella
sella_kwargs = {sella_kwargs}
fmax = 0.0001
steps = 250
if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
    sella_kwargs['internal'] = False
opt = Sella(mol, order=1, 
            trajectory='{label}.traj',
            logfile='{label}_sella.log',
            **sella_kwargs)

# run
try:
    converged = opt.run(fmax=fmax, steps=steps)
except RuntimeError:
    converged = False
traj = read('{label}.traj', index=':')
write('{label}.xyz', traj, format='xyz')
e = mol.get_potential_energy()
del mol.calc.results['forces']

if converged:
    freqs, zpe, hessian = calc_vibrations(mol, '{label}')
    if sella_freq_check(freqs, 1):
        data = {{'energy': e, 'frequencies': freqs, 'zpe': zpe,
                 'hess': hessian, 'status': 'normal'}}
    else: 
        data = {{'status': 'error'}}
else:
    sella_kwargs['internal'] = 1 - sella_kwargs['internal']
    mol.positions = {geom}
    opt = Sella(mol, order=1,
        trajectory='{label}.traj',
        logfile='{label}_sella.log',
        **sella_kwargs)

    try:
        converged = opt.run(fmax=fmax, steps=steps)
    except RuntimeError:
        converged = False
    traj = read('{label}.traj', index=':')
    write('{label}.xyz', traj, format='xyz')
    e = mol.get_potential_energy()
    del mol.calc.results['forces']
    if converged:
        freqs, zpe, hessian = calc_vibrations(mol, '{label}')
        if sella_freq_check(freqs, 1):
            data = {{'energy': e, 'frequencies': freqs, 'zpe': zpe,
                     'hess': hessian, 'status': 'normal'}}
        else: 
            data = {{'status': 'error'}}
    else: 
        data = {{'status': 'error'}}
 
db.write(mol, name='{label}', data=data)
with open('{label}_sella.log', 'a') as f:
    f.write('done\n')
