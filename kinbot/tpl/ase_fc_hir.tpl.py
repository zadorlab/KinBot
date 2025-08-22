import os
import numpy as np
from ase import Atoms
from ase.db import connect
from ase.io import read, write
from sella import Sella, Constraints
import random
#from kinbot.ase_modules.calculators.{code} import {Code}
from fairchem.core import pretrained_mlip, FAIRChemCalculator

with open('fairchem.log', 'a') as f:
    f.write('{label} | Performing hindered rotors\n')


db = connect('{working_dir}/kinbot.db')

mol = Atoms(symbols={atom}, 
            positions={geom})

const = Constraints(mol)
base_0_fix = [idx - 1 for idx in {fix}]
const.fix_dihedral(base_0_fix)

kwargs = {kwargs}
mol.info.update({{"charge": kwargs['charge'], "spin": kwargs['mult']}})

mol.calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit("uma-s-1", device="cpu"), task_name="omol")

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

sella_kwargs = {sella_kwargs}
opt = Sella(mol, 
            order={order}, 
            constraints=const,
            trajectory='{label}.traj', 
            logfile='{label}_sella.log',
            **sella_kwargs)
try:
    converged = opt.run(fmax=0.001, steps=300)
    traj = read('{label}.traj', index=':')
    write('{label}.xyz', traj, format='xyz')
    if converged:
        e = mol.get_potential_energy()
        forces = mol.calc.results['forces']
        del mol.calc.results['forces']
        random.seed()
        db.write(mol, name='{label}', data={{'energy': e, 'forces': forces, 'status': 'normal'}})
    else:  # TODO Eventually we might want to correct something in case it fails.
        raise RuntimeError
except (RuntimeError, ValueError):
    del mol.calc.results['forces']
    random.seed()
    db.write(mol, name='{label}', data={{'status': 'error'}})

with open('{label}.log', 'a') as f:
    f.write('done\n')

with open('fairchem.log', 'a') as f:
    f.write('{label} | Hindered rotors success!\n')

