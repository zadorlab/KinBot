import os
import random
import sys

import numpy as np
from ase import Atoms
from ase.db import connect
from ase.io import read, write
from ase.optimize import BFGS
from sella import Sella

#from kinbot.ase_modules.calculators.{code} import {Code}
from fairchem.core import pretrained_mlip, FAIRChemCalculator
from kinbot.frequencies import calc_vibrations

db = connect('{working_dir}/kinbot.db')
mol = Atoms(symbols={atom}, 
            positions={geom})

kwargs = {kwargs}
mol.info.update({{"charge": kwargs['charge'], "spin": kwargs['mult']}})

mol.calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit("uma-s-1", device="cpu"), task_name="omol")

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

# For monoatomic wells, just calculate the energy and exit. 
if len(mol) == 1:
    e = mol.get_potential_energy()
    forces = mol.calc.results['forces']
    del mol.calc.results['forces']
    random.seed()
    db.write(mol, name='{label}',
             data={{'energy': e, 'frequencies': np.array([]), 'zpe': 0.0,
                 'hess': np.zeros([3, 3]), 'status': 'normal'}})
    with open('{label}_sella.log', 'a') as f:
        f.write('Sella optimization is not needed for atoms.\ndone\n')
    sys.exit(0)

order = {order}
sella_kwargs = {sella_kwargs}
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
fmax = 1e-4
steps = 250
converged = False
try:
    mol.calc.label = '{label}'
    try:
        converged = opt.run(fmax=fmax, steps=steps)
        traj = read('{label}.traj', index=':')
        write('{label}.xyz', traj, format='xyz')
    except ValueError:
        mol.set_positions({geom} + np.random.normal(scale=0.05, size=(len(mol), 3)))
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
        converged = opt.run(fmax=fmax, steps=steps)
        traj = read('{label}.traj', index=':')
        write('{label}.xyz', traj, format='xyz')

    freqs, zpe, hessian = calc_vibrations(mol, '{label}')

    if order == 0 and (np.count_nonzero(np.array(freqs) < 0) > 1
                        or np.count_nonzero(np.array(freqs) < -50) >= 1):
        converged = False
    elif order == 1 and (np.count_nonzero(np.array(freqs) < 0) > 2  # More than two imag frequencies
                         or np.count_nonzero(np.array(freqs) < -50) >= 2  # More than one imag frequency larger than 50i
                         or np.count_nonzero(np.array(freqs) < 0) == 0):  # No imaginary frequencies
        converged = False
    else:
        converged = True
        e = mol.get_potential_energy()
        forces = mol.calc.results['forces']
        del mol.calc.results['forces']
        random.seed()
        db.write(mol, name='{label}', 
                    data={{'energy': e, 'frequencies': freqs, 'zpe': zpe, 
                    'hess': hessian, 'forces': forces, 'status': 'normal'}})
    if not converged:
        raise RuntimeError("Did not converge")

except (RuntimeError, ValueError):
    data = {{'status': 'error'}}
    if freqs:
        data['frequencies'] = freqs
    db.write(mol, name='{label}', data=data)

with open('{label}_sella.log', 'a') as f:
    f.write('done\n')
