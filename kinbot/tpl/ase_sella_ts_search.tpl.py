import os
import numpy as np
import shutil

from ase import Atoms
from ase.db import connect
from sella import Sella, Constraints

from kinbot.ase_modules.calculators.{code} import {Code}

db = connect('{working_dir}/kinbot.db')
mol = Atoms(symbols={atom}, 
            positions={geom})

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

kwargs = {kwargs}
if '{Code}' == 'ORCA':
    from kinbot.ase_modules.calculators.orca import OrcaProfile
    kwargs['profile'] = OrcaProfile(command=kwargs['profile'])

mol.calc = {Code}(**kwargs)
if '{Code}' == 'Gaussian':
    mol.get_potential_energy()
    kwargs['guess'] = 'Read'
    mol.calc = {Code}(**kwargs)

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

# intermediate steps don't need to fully converge
try:
    converged = opt.run(fmax=0.1, steps=100)
    e = mol.get_potential_energy()
    #if not converged:
        #raise RuntimeError
except:
    sella_kwargs['internal'] = 1 - sella_kwargs['internal']
    opt = Sella(mol,
            order=0,
            constraints=const,
            trajectory='{label}.traj',
            logfile='{label}_sella.log',
            **sella_kwargs)
    converged = opt.run(fmax=0.1, steps=100)
    e = mol.get_potential_energy()

if not mol.positions.any():  # If all coordinates are 0
    mol.positions = {geom}   # Reset to the original geometry
db.write(mol, name='{label}', data={{'energy': e, 'status': 'normal'}})

if os.path.isdir('{label}'):
    shutil.rmtree('{label}')

with open('{label}_sella.log', 'a') as f:
    f.write('done\n')
