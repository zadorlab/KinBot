import os
import sys
import pickle
import shutil

import numpy as np
from ase import Atoms
from ase.vibrations import Vibrations
from ase.optimize import BFGS
from ase.io import read, write
from sella import Sella

from kinbot.constants import EVtoHARTREE
from kinbot.ase_modules.calculators.{code} import {Code}
from kinbot.stationary_pt import StationaryPoint
from kinbot.frequencies import calc_vibrations, get_frequencies

scratch_dir = os.getcwd()

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

basename = os.path.basename('{label}')
frequency_mode = '{frequency_mode}'
pkl_file = '{label}.pkl'
if os.path.isfile(pkl_file):
    os.remove(pkl_file)


def publish(data):
    """Atomically hand one worker result to the KinBot driver."""
    payload = {{'sym': mol.symbols, 'pos': mol.positions,
               'calc': '{code}', 'name': '{label}', 'data': data}}
    temporary = f'{{pkl_file}}.tmp.{{os.getpid()}}'
    with open(temporary, 'wb') as stream:
        pickle.dump(payload, stream)
    os.replace(temporary, pkl_file)

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

# For monoatomic wells, just calculate the energy and exit. 
if len(mol) == 1:
    e = mol.get_potential_energy()
    data = {{'energy': e, 'frequencies': np.array([]), 'zpe': 0.0,
            'hess': np.zeros([3, 3]), 'status': 'normal'}}

    if os.path.isdir(f'{{basename}}'):
        shutil.rmtree(f'{{basename}}')

    with open('{label}_sella.log', 'a') as f:
        f.write('Sella optimization is not needed for atoms.\ndone\n')
else:
    data = {{'status': 'error'}}
    order = {order}
    sella_kwargs = {sella_kwargs}
    if sella_kwargs['internal'] == True and len(mol.symbols) < 5:
        sella_kwargs['internal'] = False

    if len(mol.symbols) > 2:
        opt = Sella(mol, 
                    order=order, 
                    trajectory='{label}.traj',
                    logfile='{label}_sella.log',
                    **sella_kwargs)
    else:
        opt = BFGS(mol,
                   trajectory='{label}.traj',
                   logfile='{label}_sella.log')
    freqs = []
    mol.calc.label = '{label}'
    converged = False
    try:
        converged = opt.run(fmax={fmax}, steps={steps})
        traj = read('{label}.traj', index=':')
        write('{label}.xyz', traj, format='xyz')
    except:
        converged = False
    if converged:
        try:
            if frequency_mode == 'native_hessian' and '{Code}' == 'Gaussian':
                # Sella owns the optimization. Run one fixed-geometry native
                # Gaussian frequency calculation afterward; no Gaussian
                # optimization keyword is present.
                from kinbot import reader_gauss
                from kinbot.utils import iowait

                frequency_kwargs = dict(kwargs)
                frequency_kwargs.pop('force', None)
                frequency_kwargs.pop('opt', None)
                frequency_kwargs['freq'] = ''
                frequency_kwargs['label'] = '{label}'
                # Rendered checkpoint option: {{'chk': '{checkpoint}'}}
                frequency_kwargs['chk'] = os.path.basename('{label}')
                mol.calc = Gaussian(**frequency_kwargs)
                e = mol.get_potential_energy()
                iowait('{label}.log', 'gauss')
                freqs = reader_gauss.read_freq('{label}.log', {atom})
                zpe = reader_gauss.read_zpe('{label}.log')
                hessian = None
            else:
                freqs, zpe, hessian = calc_vibrations(mol, f'{{basename}}')
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
                data = {{'energy': e, 'frequencies': freqs, 'zpe': zpe,
                        'status': 'normal'}}
                if hessian is not None:
                    data['hess'] = hessian
                pass
        except Exception as error:
            with open('{label}_sella.log', 'a') as f:
                f.write(f'Frequency evaluation failed: '
                        f'{{type(error).__name__}}: {{error}}\n')
            converged = False

    os.chdir(scratch_dir)

    if not converged:
        data = {{'status': 'error'}}
    
    if os.path.isdir(f'{{basename}}'):
        shutil.rmtree(f'{{basename}}')

    if os.path.isdir(f'{{basename}}_vib'):
        shutil.rmtree(f'{{basename}}_vib')

    with open('{label}_sella.log', 'a') as f:
        f.write('done\n')

publish(data)
