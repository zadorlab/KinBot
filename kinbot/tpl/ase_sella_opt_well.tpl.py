import os
import sys
import shutil

import numpy as np
from ase import Atoms
from ase.db import connect
from ase.vibrations import Vibrations
from ase.optimize import BFGS
from ase.io import read, write
from sella import Sella

from kinbot.constants import EVtoHARTREE
from kinbot.ase_modules.calculators.{code} import {Code}
from kinbot.stationary_pt import StationaryPoint
from kinbot.frequencies import profile_vibrations, calc_vibrations, get_frequencies

scratch_dir = os.getcwd()

db = connect('{working_dir}/kinbot.db')
mol = Atoms(symbols={atom}, 
            positions={geom})

kwargs = {kwargs}
if '{Code}' == 'ORCA':
    from kinbot.ase_modules.calculators.orca import OrcaProfile
    kwargs['profile'] = OrcaProfile(command=kwargs['profile'])

mol.calc = {Code}(**kwargs)

basename = os.path.basename('{label}')
output_label = '{label}' if {profiled} else basename

if os.path.isfile(f'{{output_label}}_sella.log'):
    os.remove(f'{{output_label}}_sella.log')

# For monoatomic wells, just calculate the energy and exit. 
if len(mol) == 1:
    try:
        e = mol.get_potential_energy()
        data = {{'energy': e, 'frequencies': np.array([]), 'zpe': 0.0,
                 'hess': np.zeros([3, 3]), 'status': 'normal'}}
    except Exception as error:
        print(f'{label}: initial energy calculation failed: {{error}}', file=sys.stderr)
        data = {{'status': 'error'}}
    db.write(mol, name='{label}', data=data)

    if os.path.isdir(f'{{basename}}'):
        shutil.rmtree(f'{{basename}}')

    with open(f'{{output_label}}_sella.log', 'a') as f:
        f.write('Sella optimization is not needed for atoms.\ndone\n')
else:
    order = {order}
    sella_kwargs = {sella_kwargs}
    if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
        sella_kwargs['internal'] = False

    freqs = []
    mol.calc.label = '{label}'
    converged = False
    try:
        if '{Code}' == 'Gaussian':
            mol.get_potential_energy()
            kwargs['guess'] = 'Read'
            mol.calc = {Code}(**kwargs)
        if len(mol.symbols) > 2 and not ('{optimizer}' == 'bfgs' and order == 0):
            opt = Sella(mol,
                        order=order,
                        trajectory=f'{{output_label}}.traj',
                        logfile=f'{{output_label}}_sella.log',
                        **sella_kwargs)
        else:
            opt = BFGS(mol,
                       trajectory=f'{{output_label}}.traj',
                       logfile=f'{{output_label}}_sella.log')
        converged = opt.run(fmax={fmax}, steps={steps})
        traj = read(f'{{output_label}}.traj', index=':')
        write(f'{{output_label}}.xyz', traj, format='xyz')
    except Exception as error:
        print(f'{label}: energy calculation or optimization failed: {{error}}', file=sys.stderr)
        converged = False
    if converged:
        try:
            freqs, zpe, hessian = profile_vibrations(mol, '{label}', '{frequency_mode}')
            if freqs is None:
                converged = False
            elif order == 0 and (np.count_nonzero(np.array(freqs) < 0) > 1
                           or np.count_nonzero(np.array(freqs) < -50) >= 1):
                converged = False
            elif order == 1 and (np.count_nonzero(np.array(freqs) < 0) > 2  # More than two imag frequencies
                             or np.count_nonzero(np.array(freqs) < -50) >= 2  # More than one imag frequency larger than 50i
                             or np.count_nonzero(np.array(freqs) < 0) == 0):  # No imaginary frequencies
                converged = False
            else:
                e = mol.get_potential_energy()
                db.write(mol, name='{label}', 
                         data={{'energy': e, 'frequencies': freqs, 'zpe': zpe, 
                         'hess': hessian, 'status': 'normal'}})
        except:
            converged = False

    os.chdir(scratch_dir)

    if not converged:
        data = {{'status': 'error'}}
        db.write(mol, name='{label}', data=data)
    
    if os.path.isdir(f'{{basename}}'):
        shutil.rmtree(f'{{basename}}')

    if {profiled} and os.path.exists(f'{{basename}}_vib/vib.xyz'):
        shutil.copy2(f'{{basename}}_vib/vib.xyz', f'{{output_label}}_vib.xyz')
    if os.path.isdir(f'{{basename}}_vib'):
        shutil.rmtree(f'{{basename}}_vib')

    with open(f'{{output_label}}_sella.log', 'a') as f:
        f.write('done\n')
