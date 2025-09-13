import os
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.db import connect

from sella import Sella, IRC

#from kinbot.ase_modules.calculators.{code} import {Code}
from kinbot.frequencies import calc_vibrations
from fairchem.core import pretrained_mlip, FAIRChemCalculator

db = connect('{working_dir}/kinbot.db')
if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

# molecule
mol = Atoms(symbols={atom}, 
            positions={geom})
kwargs = {kwargs}
mol.info.update({{"charge": kwargs['charge'], "spin": kwargs['mult']}})
mol.calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit("uma-s-1", device="cpu"), task_name="omol")

# irc
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

# run
try:
    converged_irc = irc.run(fmax=0.01, steps=100, direction=direction)
except RuntimeError:
    pass
e = mol.get_potential_energy()
del mol.calc.results['forces']

db.write(mol, name='{label}', data={{'energy': e, 'status': 'normal'}})
with open('{label}_sella.log', 'a') as f:
    f.write('done\n')

# product optimization
#prod_kwargs = {prod_kwargs}

# sella
sella_kwargs = {sella_kwargs}
opt = Sella(mol, 
            order=0, 
            trajectory='{label}_prod.traj', 
            logfile='{label}_prod_sella.log',
            **sella_kwargs)

try:
    converged_opt = opt.run(fmax=0.0001, steps=250)
except RuntimeError:
    pass
traj = read('{label}.traj', index=':')
write('{label}.xyz', traj, format='xyz')
e = mol.get_potential_energy()
del mol.calc.results['forces']
freqs, zpe, hessian = calc_vibrations(mol, '{label}')

# product optimizations will not always converge, it's okay
data={{'energy': e, 'frequencies': freqs, 'zpe': zpe,
            'hess': hessian, 'status': 'normal'}}
db.write(mol, name='{label}_prod', data=data)

with open('{label}_prod_sella.log', 'a') as f:
    f.write('done\n')
