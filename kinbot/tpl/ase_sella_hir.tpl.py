import os
import shutil
import numpy as np

from ase import Atoms
from ase.db import connect
from ase.io import read, write
from sella import Sella, Constraints

from kinbot.ase_modules.calculators.{code} import {Code}

db = connect('{working_dir}/kinbot.db')
mol = Atoms(symbols={atom}, 
            positions={geom})

const = Constraints(mol)
base_0_fix = [idx - 1 for idx in {fix}]
const.fix_dihedral(base_0_fix)

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
            order={order}, 
            constraints=const,
            trajectory='{label}.traj', 
            logfile='{label}_sella.log',
            **sella_kwargs)

converged = opt.run(fmax={fmax}, steps={steps})
traj = read('{label}.traj', index=':')
write('{label}.xyz', traj, format='xyz')
if converged:
    e = mol.get_potential_energy()
    db.write(mol, name='{label}', data={{'energy': e, 'status': 'normal'}})

else:
    db.write(mol, name='{label}', data={{'status': 'error'}})

if os.path.isdir('{label}'):
    shutil.rmtree('{label}')

with open('{label}_sella.log', 'a') as f:
    f.write('done\n')
