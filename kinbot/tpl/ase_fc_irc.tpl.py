import os
import random
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.db import connect

from sella import Sella, IRC

#from kinbot.ase_modules.calculators.{code} import {Code}
from kinbot.frequencies import calc_vibrations
from fairchem.core import pretrained_mlip, FAIRChemCalculator

db = connect('{working_dir}/kinbot.db')
mol = Atoms(symbols={atom}, 
            positions={geom})

kwargs = {kwargs}
mol.info.update({{"charge": kwargs['charge'], "spin": kwargs['mult']}})

mol.calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit("uma-s-1", device="cpu"), task_name="omol")

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

irc = IRC(mol, 
         trajectory='{label}.traj', 
         dx=0.1, 
         eta=1e-4, 
         gamma=0, 
         keep_going=True,
         logfile='{label}_sella.log')

if '{label}'.endswith('F'):
    direction = 'forward'
elif '{label}'.endswith('R'):
    direction = 'reverse'
    
converged_irc = irc.run(fmax=0.01, steps=100, direction=direction)
e = mol.get_potential_energy()
forces = mol.calc.results['forces']
del mol.calc.results['forces']
random.seed()
db.write(mol, name='{label}', data={{'energy': e, 'forces': forces, 'status': 'normal'}})

with open('{label}_sella.log', 'a') as f:
    f.write('done\n')

# product optimization
prod_kwargs = {prod_kwargs}
mol.calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit("uma-s-1", device="cpu"), task_name="omol")

sella_kwargs = {sella_kwargs}
opt = Sella(mol, 
            order=0, 
            trajectory='{label}_prod.traj', 
            logfile='{label}_prod_sella.log',
            **sella_kwargs)
try:
    converged_opt = opt.run(fmax=0.0001, steps=300)
    traj = read('{label}.traj', index=':')
    write('{label}.xyz', traj, format='xyz')
    if converged_opt:
        e = mol.get_potential_energy()
        freqs, zpe, hessian = calc_vibrations(mol, '{label}')
        forces = mol.calc.results['forces']
        del mol.calc.results['forces']
        random.seed()
        db.write(mol, name='{label}_prod', 
                    data={{'energy': e, 'frequencies': freqs, 'zpe': zpe, 
                    'hess': hessian, 'forces': forces, 'status': 'normal'}})
    else:
        raise RuntimeError
except (RuntimeError, ValueError):
    forces = mol.calc.results['forces']
    del mol.calc.results['forces']
    random.seed()
    db.write(mol, name='{label}_prod', data={{'status': 'error'}})    

with open('{label}_prod_sella.log', 'a') as f:
    f.write('done\n')
