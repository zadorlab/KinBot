import os
import random
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.db import connect

from sella import Sella, IRC

#from kinbot.ase_modules.calculators.{code} import {Code}
from fairchem.core import pretrained_mlip, FAIRChemCalculator

db = connect('{working_dir}/kinbot.db')
mol = Atoms(symbols={atom}, 
            positions={geom})
with open('fairchem.log', 'a') as f:
    f.write('{label} | Starting IRC Calculation...\n')

kwargs = {kwargs}
mol.info.update({{"charge": kwargs['charge'], "spin": kwargs['mult']}})

mol.calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit("uma-s-1", device="cpu"), task_name="omol")

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')
irc = IRC(mol, trajectory='{label}.traj', dx=0.1, eta=1e-4, gamma=0, 
          logfile='{label}_sella.log')
if '{label}'.endswith('F'):
    direction = 'forward'
elif '{label}'.endswith('R'):
    direction = 'reverse'
else:
    raise ValueError('Unexpected IRC name: {label}.')
try:
    converged_irc = irc.run(fmax=0.01, steps=100, direction=direction)
    if converged_irc:
        with open('fairchem.log', 'a') as f:
            f.write('{label} | IRC converged!\n')

        e = mol.get_potential_energy()
        forces = mol.calc.results['forces']
        del mol.calc.results['forces']
        random.seed()
        db.write(mol, name='{label}', data={{'energy': e, 'forces': forces, 'status': 'normal'}})
        success = True
    elif mol.positions is not None and mol.positions.any():
        # although there is an error, continue from the final geometry
        with open('fairchem.log', 'a') as f:
            f.write('{label} | Continuing from final geometry\n')
        forces = mol.calc.results['forces']
        del mol.calc.results['forces']
        random.seed()
        db.write(mol, name='{label}', data={{'forces': forces, 'status': 'normal'}})
        success = True
    else:
        raise RuntimeError
except (RuntimeError, ValueError):
    if mol.positions is not None and mol.positions.any():
        # although there is an error, continue from the final geometry
        with open('fairchem.log', 'a') as f:
            f.write('{label} | Contuining from final geometry_2\n')
        forces = mol.calc.results['forces']
        del mol.calc.results['forces']
        random.seed()
        db.write(mol, name='{label}', data={{'forces': forces, 'status': 'normal'}})
        success = True
    else:
        success = False
        with open('fairchem.log', 'a') as f:
            f.write('{label} | IRC Failed\n')
        
        del mol.calc.results['forces']
        random.seed()
        db.write(mol, name='{label}', data={{'status': 'error'}})

with open('{label}.log', 'a') as f:
    f.write('done\n')

if success:
    prod_kwargs = {prod_kwargs}
    mol.calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit("uma-s-1", device="cpu"), task_name="omol")
    with open('fairchem.log', 'a') as f:
        f.write('{label} | Optimizing IRC\n')

    sella_kwargs = {sella_kwargs}
    opt = Sella(mol, 
                order=0, 
                trajectory='{label}_prod.traj', 
                logfile='{label}_prod_sella.log',
                **sella_kwargs)
    try:
        converged_opt = opt.run(fmax=0.0001, steps=300)
        traj = read('{label}.traj', index=':')
        write('{label}.xyz', traj, format='xyz')
        if converged_opt:
            with open('fairchem.log', 'a') as f:
                f.write('{label} | IRC Optimization Successfully Converged!\n')

            e = mol.get_potential_energy()
            forces = mol.calc.results['forces']
            del mol.calc.results['forces']
            random.seed()
            db.write(mol, name='{label}_prod', data={{'energy': e, 'forces': forces, 'status': 'normal'}})
        else:
            raise RuntimeError
    except (RuntimeError, ValueError):
        with open('fairchem.log', 'a') as f:
            f.write('{label} | IRC optimization failed\n')
        forces = mol.calc.results['forces']
        del mol.calc.results['forces']
        random.seed()
        db.write(mol, name='{label}_prod', data={{'status': 'error'}})    
    
    with open('fairchem.log', 'a') as f:
        f.write('{label} | IRC Complete\n')

    with open('{label}_prod.log', 'a') as f:
        f.write('done\n')
