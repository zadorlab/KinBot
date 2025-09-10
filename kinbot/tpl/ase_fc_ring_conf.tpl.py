import os
import random
import sys
import shutil

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.db import connect
from sella import Sella, Constraints

from kinbot.modify_geom import modify_coordinates
#from kinbot.ase_modules.calculators.{code} import {Code}
from fairchem.core import pretrained_mlip, FAIRChemCalculator

from kinbot.stationary_pt import StationaryPoint
from kinbot.frequencies import get_frequencies

db = connect('{working_dir}/kinbot.db')

mol = Atoms(symbols={atom}, 
            positions={geom})

kwargs = {kwargs}
mol.info.update({{"charge": kwargs['charge'], "spin": kwargs['mult']}})

mol.calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit("uma-s-1", device="cpu"), task_name="omol")

const = Constraints(mol)
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

st_pt = StationaryPoint.from_ase_atoms(mol)
st_pt.characterize()
base_0_changes = []
for c in {change}:
    c_new = [ci - 1 for ci in c[:-1]]
    c_new.append(c[-1])
    base_0_changes.append(c_new)
if len(base_0_changes) > 0:
    _, mol.positions = modify_coordinates(st_pt, '{label}', mol.positions, 
                                           base_0_changes, st_pt.bond)
for change in base_0_changes:
    if len(change[:-1]) == 2:
        const.fix_bond(change[:-1])
    elif len(change[:-1]) == 3:
        const.fix_angle(change[:-1])
    elif len(change[:-1]) == 4:
        const.fix_dihedral(change[:-1])
    else:
        raise ValueError(f'Unexpected length of changes: {{change}}.')

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

order = {order}
sella_kwargs = {sella_kwargs}
opt = Sella(mol, 
            order=order, 
            trajectory='{label}.traj', 
            logfile='{label}_sella.log',
            **sella_kwargs)
freqs = []
fmax = 1e-4
steps = 250
mol.calc.label = '{label}'
converged = opt.run(fmax=fmax, steps=steps)
traj = read('{label}.traj', index=':')
write('{label}.xyz', traj, format='xyz')
e = mol.get_potential_energy()
forces = mol.calc.results['forces']
del mol.calc.results['forces']
random.seed()
db.write(mol, name='{label}', 
         data={{'energy': e, 'forces': forces, 'status': 'normal'}})
with open('{label}_sella.log', 'a') as f:
    f.write('done\n')
