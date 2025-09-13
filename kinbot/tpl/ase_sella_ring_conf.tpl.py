import os
import sys
import shutil

import numpy as np
from ase import Atoms
from ase.db import connect
from sella import Sella, Constraints

from kinbot.ase_modules.calculators.{code} import {Code}
from kinbot.stationary_pt import StationaryPoint


db = connect('{working_dir}/kinbot.db')
mol = Atoms(symbols={atom}, 
            positions={geom})

kwargs = {kwargs}
if '{Code}' == 'ORCA':
    from kinbot.ase_modules.calculators.orca import OrcaProfile
    kwargs['profile'] = OrcaProfile(command=kwargs['profile'])

mol.calc = {Code}(**kwargs)
if '{Code}' == 'Gaussian':
    mol.get_potential_energy()
    mol.calc = {Code}(**kwargs)

const = Constraints(mol)
fix_these = [[idx - 1 for idx in fix] for fix in {fix}]
for fix in fix_these:
    if len(fix) == 2:
        const.fix_bond(fix)
    elif len(fix) == 4:
        const.fix_dihedral(fix)
    else:
        raise ValueError(f'Unexpected length of fix: {{fix}}.')

for c in {change}:
    const.fix_dihedral((c[0]-1, c[1]-1, c[2]-1, c[3]-1), target=c[4])

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
            **sella_kwargs,
            )
mol.calc.label = '{label}'

try:
    opt.run(fmax=1e-4, steps=250)
    e = mol.get_potential_energy()
    db.write(mol, name='{label}', 
             data={{'energy': e, 'status': 'normal'}})
except (RuntimeError, ValueError):
    try:
        sella_kwargs['internal'] = 1 - sella_kwargs['internal']
        opt.run(fmax=1e-4, steps=250)
        e = mol.get_potential_energy()
        db.write(mol, name='{label}',
             data={{'energy': e, 'status': 'normal'}})
    except:
        data = {{'status': 'error'}}
        db.write(mol, name='{label}', data=data)

if os.path.isdir('{label}'):
    shutil.rmtree('{label}')

with open('{label}_sella.log', 'a') as f:
    f.write('done\n')
