import os
import shutil
import numpy as np

from ase import Atoms
from ase.db import connect
from ase.io import read, write

from sella import Sella, IRC

from kinbot.ase_modules.calculators.{code} import {Code}
from kinbot.frequencies import profile_vibrations, calc_vibrations

db = connect('{working_dir}/kinbot.db')
mol = Atoms(symbols={atom},
            positions={geom})

kwargs = {kwargs}
if '{Code}' == 'ORCA':
    from kinbot.ase_modules.calculators.orca import OrcaProfile
    kwargs['profile'] = OrcaProfile(command=kwargs['profile'])

mol.calc = {Code}(**kwargs)
if '{Code}' == 'Gaussian' and kwargs.get('chk'):
    mol.get_potential_energy()
    kwargs['guess'] = 'Read'
    mol.calc = {Code}(**kwargs)

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')
irc = IRC(mol,
          trajectory='{label}.traj',
          dx=0.1,
          eta=1e-4,
          gamma=0,
          logfile='{label}_sella.log',
          keep_going=True)
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
        db.write(mol, name='{label}', data={{'status': 'normal'}})
        success = True
    else:
        raise RuntimeError
except (RuntimeError, ValueError):
    if mol.positions is not None and mol.positions.any():
        db.write(mol, name='{label}', data={{'status': 'normal'}})
        success = True
    else:
        success = False
        db.write(mol, name='{label}', data={{'status': 'error'}})

if os.path.isfile('{label}.traj'):
    write('{label}.xyz', read('{label}.traj', index=':'), format='xyz')

with open('{label}_sella.log', 'a') as handle:
    handle.write('done\n')

if success:
    prod_kwargs = {prod_kwargs}
    if '{Code}' == 'ORCA':
        prod_kwargs['profile'] = kwargs['profile']

    mol.calc = {Code}(**prod_kwargs)
    if '{Code}' == 'Gaussian' and prod_kwargs.get('chk'):
        mol.get_potential_energy()
        prod_kwargs['guess'] = 'Read'
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
        if not converged_opt:
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
            if not converged_opt:
                raise RuntimeError
        except Exception:
            db.write(mol, name='{label}_prod', data={{'status': 'error'}})
        else:
            e = mol.get_potential_energy()
            freqs, zpe, hessian = profile_vibrations(mol, '{label}_prod', '{frequency_mode}')
            db.write(mol, name='{label}_prod', data={{
                'energy': e, 'status': 'normal', 'frequencies': freqs,
                'zpe': zpe, 'hess': hessian}})
    else:
        e = mol.get_potential_energy()
        freqs, zpe, hessian = profile_vibrations(mol, '{label}_prod', '{frequency_mode}')
        db.write(mol, name='{label}_prod', data={{
            'energy': e, 'status': 'normal', 'frequencies': freqs,
            'zpe': zpe, 'hess': hessian}})

    if os.path.isfile('{label}_prod.traj'):
        write('{label}_prod.xyz', read(
            '{label}_prod.traj', index=':'), format='xyz')
    with open('{label}_prod_sella.log', 'a') as handle:
        handle.write('done\n')
if os.path.isdir('{label}'):
    shutil.rmtree('{label}')
