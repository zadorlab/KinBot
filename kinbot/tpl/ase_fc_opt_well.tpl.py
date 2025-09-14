import os
import sys
import pickle

import numpy as np
from ase import Atoms
#from ase.db import connect
from ase.io import read, write
from ase.optimize import BFGS
from sella import Sella

#from kinbot.ase_modules.calculators.{code} import {Code}
from fairchem.core import pretrained_mlip, FAIRChemCalculator
from kinbot.frequencies import calc_vibrations
from kinbot.utils import sella_freq_check

#db = connect('{working_dir}/kinbot.db')
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

# For monoatomic wells, just calculate the energy and exit. 
if len(mol) == 1:
    e = mol.get_potential_energy()
    del mol.calc.results['forces']
    data={{'energy': e, 'frequencies': np.array([]), 'zpe': 0.0,
                 'hess': np.zeros([3, 3]), 'status': 'normal'}}
    mol_pkl = {{'mol': mol, 'name': '{label}', 'data': data}}
    with open('{label}.pickle', 'wb') as f:
        pickle.dump(mol_pkl, f)
    #db.write(mol, name='{label}',
    #         data={{'energy': e, 'frequencies': np.array([]), 'zpe': 0.0,
    #             'hess': np.zeros([3, 3]), 'status': 'normal'}})
    with open('{label}_sella.log', 'a') as f:
        f.write('Sella optimization is not needed for atoms.\ndone\n')
    sys.exit(0)

# sella
sella_kwargs = {sella_kwargs}
fmax = 0.0001
steps = 250
if len(mol.symbols) > 2:
    if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
        sella_kwargs['internal'] = False
    opt = Sella(mol, 
                order={order}, 
                trajectory='{label}.traj', 
                logfile='{label}_sella.log',
                **sella_kwargs)
else:
    opt = BFGS(mol, 
               trajectory='{label}.traj',
               logfile='{label}_sella.log')

#run
try:
    converged = opt.run(fmax=fmax, steps=steps)
except:
    converged = False
e = mol.get_potential_energy()
del mol.calc.results['forces']
traj = read('{label}.traj', index=':')
write('{label}.xyz', traj, format='xyz')

if converged:
    freqs, zpe, hessian = calc_vibrations(mol, '{label}')
    if sella_freq_check(freqs, {order}):
        data={{'energy': e, 'frequencies': freqs, 'zpe': zpe,
            'hess': hessian, 'status': 'normal'}}
elif len(mol.symbols) > 2 and (not converged):
    mol.positions = {geom}
    sella_kwargs['internal'] = 1 - sella_kwargs['internal']
    opt = Sella(mol,
        order={order},
        trajectory='{label}.traj',
        logfile='{label}_sella.log',
        **sella_kwargs)
    try:
        converged = opt.run(fmax=fmax, steps=steps)
    except RuntimeError:
        converged = False
    e = mol.get_potential_energy()
    del mol.calc.results['forces']
    traj = read('{label}.traj', index=':')
    write('{label}.xyz', traj, format='xyz')
    if converged:
        freqs, zpe, hessian = calc_vibrations(mol, '{label}')
        if sella_freq_check(freqs, {order}):
            data={{'energy': e, 'frequencies': freqs, 'zpe': zpe,
                'hess': hessian, 'status': 'normal'}}
        else:
            data = {{'status': 'error'}}
    else:
        data = {{'status': 'error'}}
else:
    data = {{'status': 'error'}}

mol_pkl = {{'mol': mol, 'name': '{label}', 'data': data}}
with open('{label}.pickle', 'wb') as f:
    pickle.dump(mol_pkl, f)

#db.write(mol, name='{label}', data=data)

with open('{label}_sella.log', 'a') as f:
    f.write('done\n')
