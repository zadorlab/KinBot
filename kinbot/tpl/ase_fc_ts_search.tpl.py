import os
import numpy as np
import pickle

from ase import Atoms
from ase.io import read, write
from sella import Sella, Constraints

from kinbot.utils import too_far
#from kinbot.ase_modules.calculators.{code} import {Code}

#db = connect('{working_dir}/kinbot.db')
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
base_0_fix = [[idx - 1 for idx in fix] for fix in {fix}]
for fix in base_0_fix:
    if len(fix) == 2:
        const.fix_bond(fix)
    elif len(fix) == 3:
        const.fix_angle(fix)
    elif len(fix) == 4:
        const.fix_dihedral(fix)
    else:
        raise ValueError(f'Unexpected length of fix: {{fix}}')

# sella
fmax = 0.1
steps = 100
sella_kwargs = {sella_kwargs}
if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
    sella_kwargs['internal'] = False

opt = Sella(mol, 
            order=0,
            constraints=const,
            trajectory='{label}.traj', 
            logfile='{label}_sella.log',
            **sella_kwargs)

# run
try:
    converged = opt.run(fmax=fmax, steps=steps)
except:
    converged = False
traj = read('{label}.traj', index=':')
write('{label}.xyz', traj, format='xyz')
e = mol.get_potential_energy()
del mol.calc.results['forces']

if not mol.positions.any() or too_far(mol.positions):  # If all coordinates are 0 or too far
    mol.positions = {geom}   # Reset to the original geometry
data = {{'energy': e, 'status': 'normal'}}

mol_pkl = {{'sym': mol.symbols,
            'pos': mol.positions,
            'calc': 'fairchemcalculator',
            'name': '{label}',
            'data': data}}
with open('{label}.pkl', 'wb') as f:
    pickle.dump(mol_pkl, f)
with open('{label}_sella.log', 'a') as f:
    f.write('am1\ndone\n')  # is am1 is there, it'll be deleted on restart
