import os
import sys
import pickle

import numpy as np
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

const = Constraints(mol)
# make it zero indexed
base_0_fixes = [[idx - 1 for idx in fix] for fix in {fix}]
for fix in base_0_fixes:
    if len(fix) == 2:
        const.fix_bond(fix)
    elif len(fix) == 3:
        const.fix_angle(fix)
    elif len(fix) == 4:
        const.fix_dihedral(fix)
    else:
        raise ValueError(f'Unexpected length of fix: {{fix}}.')

base_0_changes = [[idx - 1 for idx in change] for change in {change}]
for c in base_0_changes:
    const.fix_dihedral(c[:-1], target=c[-1])

# sella
sella_kwargs = {sella_kwargs}
fmax = {fmax}
steps = {steps}
opt = Sella(mol, 
            order=0, 
            constraints=const,
            trajectory='{label}.traj', 
            logfile='{label}_sella.log',
            **sella_kwargs)

try:
    converged = opt.run(fmax=fmax, steps=steps)
except RuntimeError:
    pass
    
traj = read('{label}.traj', index=':')
write('{label}.xyz', traj, format='xyz')
e = mol.get_potential_energy()
del mol.calc.results['forces']
data={{'energy': e, 'status': 'normal'}}

# write even if no converged, this is an intermediate
mol_pkl = {{'sym': mol.symbols,
            'pos': mol.positions,
            'calc': 'fairchemcalculator',
            'name': '{label}',
            'data': data}}
with open('{label}.pkl', 'wb') as f:
    pickle.dump(mol_pkl, f)
with open('{label}_sella.log', 'a') as f:
    f.write('done\n')
