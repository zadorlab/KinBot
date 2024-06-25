import os
import shutil

import numpy as np
from ase import Atoms
from ase.db import connect
from ase.vibrations import Vibrations
from sella import Sella

from kinbot.constants import EVtoHARTREE
from kinbot.ase_modules.calculators.{code} import {Code}
from kinbot.stationary_pt import StationaryPoint
from kinbot.frequencies import get_frequencies


def calc_vibrations(mol):
        mol.calc.label = '{label}_vib'
        if 'chk' in mol.calc.parameters:
            del mol.calc.parameters['chk']
        # Compute frequencies in a separate temporary directory to avoid 
        # conflicts accessing the cache in parallel calculations.
        if not os.path.isdir('{label}_vib'):
            os.mkdir('{label}_vib')
        init_dir = os.getcwd()
        os.chdir('{label}_vib')
        if os.path.isdir('vib'):
            shutil.rmtree('vib')
        vib = Vibrations(mol)
        vib.run()
        # Use kinbot frequencies to avoid mixing low vib frequencies with 
        # the values associated with external rotations.
        _ = vib.get_frequencies()
        zpe = vib.get_zero_point_energy() * EVtoHARTREE
        hessian = vib.H / 97.17370087
        st_pt = StationaryPoint.from_ase_atoms(mol)
        st_pt.characterize()
        freqs, _ = get_frequencies(st_pt, hessian, st_pt.geom)
        os.chdir(init_dir)
        shutil.rmtree('{label}_vib')
        return freqs, zpe, hessian


db = connect('{working_dir}/kinbot.db')
mol = Atoms(symbols={atom},
            positions={geom})

kwargs = {kwargs}
mol.calc = {Code}(**kwargs)
if '{Code}' == 'Gaussian':
    mol.get_potential_energy()
    kwargs['guess'] = 'Read'
    mol.calc = {Code}(**kwargs)

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

sella_kwargs = {sella_kwargs}
opt = Sella(mol, order=1, 
            trajectory='{label}.traj',
            logfile='{label}_sella.log',
            **sella_kwargs)
freqs = []
try:
    converged = False
    fmax = 1e-4
    attempts = 1
    steps=500
    while not converged and attempts <= 3:
        mol.calc.label = '{label}'
        converged = opt.run(fmax=fmax, steps=steps)
        freqs, zpe, hessian = calc_vibrations(mol)
        if (np.count_nonzero(np.array(freqs) < 0) > 2  # More than two imag frequencies
                or np.count_nonzero(np.array(freqs) < -50) >= 2  # More than one frequency smaller than 50i
                or np.count_nonzero(np.array(freqs) < 0) == 0):  # No imaginary frequencies
            print(f'Wrong number of imaginary frequencies: {{freqs[6:]}}')
            converged = False
            mol.calc.label = '{label}'
            attempts += 1
            fmax *= 0.3
            if attempts <= 3:
                print(f'Retrying with a tighter criterion: fmax={{fmax}}.')
        else:
            converged = True
            e = mol.get_potential_energy()
            db.write(mol, name='{label}', 
                     data={{'energy': e, 'frequencies': freqs, 'zpe': zpe, 
                            'hess': hessian, 'status': 'normal'}})            
    if not converged:
        raise RuntimeError
except (RuntimeError, ValueError):
    data = {{'status': 'error'}}
    if freqs:
        data['frequencies'] = freqs
    db.write(mol, name='{label}', data=data)

with open('{label}.log', 'a') as f:
    f.write('done\n')
