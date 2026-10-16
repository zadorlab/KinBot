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
if '{Code}' == 'Gaussian' and kwargs.get('chk'):
    # Converge the SCF once and reuse it. Without a checkpoint file
    # Guess=Read stops Gaussian in l301 ("No data on chk file").
    mol.get_potential_energy()
    kwargs['guess'] = 'Read'
    mol.calc = {Code}(**kwargs)

# ORCA runs in a scratch directory whose contents are copied back into the
# job's directory. The other codes run from the working directory, where
# KinBot reads {label}_sella.log, so conf/ and hir/ jobs keep their path.
basename = os.path.basename('{label}') if '{Code}' == 'ORCA' else '{label}'

if os.path.isfile(f'{{basename}}_sella.log'):
    os.remove(f'{{basename}}_sella.log')

sella_kwargs = {sella_kwargs}
if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
    sella_kwargs['internal'] = False
opt = Sella(mol, 
            order={order}, 
            constraints=const,
            trajectory=f'{{basename}}.traj', 
            logfile=f'{{basename}}_sella.log',
            **sella_kwargs)

try:
    converged = opt.run(fmax={fmax}, steps={steps})
except:
    converged = False
traj = read(f'{{basename}}.traj', index=':')
write(f'{{basename}}.xyz', traj, format='xyz')
if converged:
    e = mol.get_potential_energy()
    db.write(mol, name='{label}', data={{'energy': e, 'status': 'normal'}})

else:
    db.write(mol, name='{label}', data={{'status': 'error'}})

if os.path.isdir(f'{{basename}}'):
    shutil.rmtree(f'{{basename}}')

with open(f'{{basename}}_sella.log', 'a') as f:
    f.write('done\n')
