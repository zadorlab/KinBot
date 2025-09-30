import os
import shutil
import numpy as np

from ase import Atoms
from ase.db import connect

from sella import Sella, IRC

from kinbot.ase_modules.calculators.{code} import {Code}

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
    kwargs['guess'] = 'Read'
    mol.calc = {Code}(**kwargs)

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
        e = mol.get_potential_energy()
        db.write(mol, name='{label}', data={{'energy': e, 'status': 'normal'}})
        success = True
    elif mol.positions is not None and mol.positions.any():
        # although there is an error, continue from the final geometry
        db.write(mol, name='{label}', data={{'status': 'normal'}})
        success = True
    else:
        raise RuntimeError
except (RuntimeError, ValueError):
    if mol.positions is not None and mol.positions.any():
        # although there is an error, continue from the final geometry
        db.write(mol, name='{label}', data={{'status': 'normal'}})
        success = True
    else:
        success = False
        db.write(mol, name='{label}', data={{'status': 'error'}})

with open('{label}_sella.log', 'a') as f:
    f.write('done\n')

if success:
    prod_kwargs = {prod_kwargs}
    if '{Code}' == 'ORCA':
        prod_kwargs['profile'] = kwargs['profile']

    mol.calc = {Code}(**prod_kwargs)
    if '{Code}' == 'Gaussian':
        mol.get_potential_energy()
        kwargs['guess'] = 'Read'
        mol.calc = {Code}(**prod_kwargs)
    sella_kwargs = {sella_kwargs}
    if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
        sella_kwargs['internal'] = False
    opt = Sella(mol, 
                order=0, 
                trajectory='{label}_prod.traj', 
                logfile='{label}_prod_sella.log',
                **sella_kwargs)
    try:
        converged_opt = opt.run(fmax={fmax}, steps={steps})
        if converged_opt:
            e = mol.get_potential_energy()
            db.write(mol, name='{label}_prod', 
                     data={{'energy': e, 'status': 'normal'}})
        else:
            raise RuntimeError
    except (RuntimeError, ValueError):
        try:
            sella_kwargs['internal'] = 1 - sella_kwargs['internal']
            opt = Sella(mol,
                order=0,
                trajectory='{label}_prod.traj',
                logfile='{label}_prod_sella.log',
                **sella_kwargs)
            converged_opt = opt.run(fmax={fmax}, steps={steps})
            e = mol.get_potential_energy()
            db.write(mol, name='{label}_prod',
                     data={{'energy': e, 'status': 'normal'}})
        except:
            db.write(mol, name='{label}_prod', data={{'status': 'error'}})    
    

    with open('{label}_prod_sella.log', 'a') as f:
        f.write('done\n')
if os.path.isdir('{label}'):
    shutil.rmtree('{label}')
