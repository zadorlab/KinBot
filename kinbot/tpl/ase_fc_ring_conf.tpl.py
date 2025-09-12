import os
import random
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

mol = Atoms(symbols={atom}, 
            positions={geom})

kwargs = {kwargs}
mol.info.update({{"charge": kwargs['charge'], "spin": kwargs['mult']}})

mol.calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit("uma-s-1", device="cpu"), task_name="omol")

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

    if os.path.isfile('{label}_sella.log'):
        os.remove('{label}_sella.log')

    sella_kwargs = {sella_kwargs}
    if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
        sella_kwargs['internal'] = False
    opt = Sella(mol, 
                order=0, 
                constraints=const,
                trajectory='{label}.traj', 
                logfile='{label}_sella.log',
                **sella_kwargs)
    fmax = 1e-4
    steps = 250
    mol.calc.label = '{label}'
    try:
        converged = opt.run(fmax=fmax, steps=steps)
    except:
        sella_kwargs['internal'] = 1 - sella_kwargs['internal']
        opt = Sella(mol,
                order=0,
                constraints=const,
                trajectory='{label}.traj',
                logfile='{label}_sella.log',
                **sella_kwargs)
        converged = opt.run(fmax=fmax, steps=steps)
    traj = read('{label}.traj', index=':')
    write('{label}.xyz', traj, format='xyz')
    e = mol.get_potential_energy()

if converged:
    random.seed()
    db.write(mol, name='{label}', 
             data={{'energy': e, 'status': 'normal'}})
else:
    data = {{'status': 'error'}}
    db.write(mol, name='{label}', data=data)
with open('{label}_sella.log', 'a') as f:
    f.write('done\n')
