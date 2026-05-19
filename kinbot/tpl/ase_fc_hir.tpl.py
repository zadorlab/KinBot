import os
import numpy as np
import pickle

from ase import Atoms
from ase.io import read, write
from sella import Sella, Constraints

from fairchem.core.units.mlip_unit import load_predict_unit
from fairchem.core import FAIRChemCalculator

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

# molecule
mol = Atoms(symbols={atom},
            positions={geom})
kwargs = {kwargs}
mol.info.update({{"charge": kwargs['charge'], "spin": kwargs['mult']}})
mol.calc = FAIRChemCalculator(load_predict_unit('{fc_model_path}', device='{fc_device}'), task_name='{fc_task_name}')

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
fmax = {fmax}
steps = {steps}

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
