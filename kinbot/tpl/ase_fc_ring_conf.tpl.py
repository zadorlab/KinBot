import os
import sys
import shutil

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.db import connect
from sella import Sella, Constraints

#from kinbot.ase_modules.calculators.{code} import {Code}
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
mol.calc.label = '{label}'

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
fmax = 0.0001
steps = 250
if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
    sella_kwargs['internal'] = False
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

# write even if no converged, this is an intermediate
db.write(mol, name='{label}', 
         data={{'energy': e, 'status': 'normal'}})
with open('{label}_sella.log', 'a') as f:
    f.write('done\n')
