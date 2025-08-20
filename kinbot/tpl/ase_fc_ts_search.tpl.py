import os
import random
import numpy as np
from ase import Atoms
from ase.db import connect
from sella import Sella, Constraints

#from kinbot.ase_modules.calculators.{code} import {Code}
from fairchem.core import pretrained_mlip, FAIRChemCalculator

with open('fairchem.log', 'a') as f:
    f.write('{label} | Starting transition state search...\n')

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
        with open('fairchem.log', 'a') as f:
            f.write('{label} | Transition state search failed (1).\n')

        raise ValueError(f'Unexpected length of fix: {{fix}}')

kwargs = {kwargs}
mol.info.update({{"charge": kwargs['charge'], "spin": kwargs['mult']}})

mol.calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit("uma-s-1", device="cpu"), task_name="omol")
if '{Code}' == 'Gaussian':
    mol.get_potential_energy()
    kwargs['guess'] = 'Read'
    mol.calc = {Code}(**kwargs)

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

sella_kwargs = {sella_kwargs}
opt = Sella(mol, 
            order=0,
            constraints=const,
            trajectory='{label}.traj', 
            logfile='{label}_sella.log',
            **sella_kwargs)
try:
    cvgd = opt.run(fmax=0.1, steps=300)
    if cvgd:
        e = mol.get_potential_energy()
    else:  # TODO Eventually we might want to correct something in case it fails.
        with open('fairchem.log', 'a') as f:
            f.write('{label} | Transition state search did not converge.\n')
        raise RuntimeError
except (RuntimeError, ValueError):
    e = 0.0

if not mol.positions.any():  # If all coordinates are 0
    mol.positions = {geom}   # Reset to the original geometry

forces = mol.calc.results['forces']
del mol.calc.results['forces']
random.seed()
db.write(mol, name='{label}', data={{'energy': e, 'forces': forces, 'status': 'normal'}})

with open('fairchem.log', 'a') as f:
    f.write('{label} | Transition state search successful!\n')

with open('{label}.log', 'a') as f:
    f.write('done\n')
