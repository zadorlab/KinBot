import os
import numpy as np
import pickle

from ase import Atoms
from ase.db import connect
from ase.io import read, write
from sella import Sella, Constraints
#from kinbot.ase_modules.calculators.{code} import {Code}

db = connect('{working_dir}/kinbot.db')
if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

# molecule
mol = Atoms(symbols={atom}, 
            positions={geom})
kwargs = {kwargs}
mol.info.update({{"charge": kwargs['charge'], "spin": kwargs['mult']}})
with open('fc_model.pkl', 'rb') as f:
    mol.calc = pickle.load(f)

# constraints 
const = Constraints(mol)
base_0_fix = [idx - 1 for idx in {fix}]
const.fix_dihedral(base_0_fix)

# optimizer
sella_kwargs = {sella_kwargs}
if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
    sella_kwargs['internal'] = False
opt = Sella(mol, 
            order={order}, 
            constraints=const,
            trajectory='{label}.traj', 
            logfile='{label}_sella.log',
            **sella_kwargs)
fmax = 0.001
steps = 250

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
    data={{'energy': e, 'status': 'normal'}}
else:
    data={{'status': 'error'}}

mol_pkl = {{'sym': mol.symbols, 
            'pos': mol.positions, 
            'calc': 'fairchemcalculator', 
            'name': '{label}', 
            'data': data}}
with open('{label}.pkl', 'wb') as f:
    pickle.dump(mol_pkl, f)

with open('{label}_sella.log', 'a') as f:
    f.write('done\n')
