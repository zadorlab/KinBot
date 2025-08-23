import os
import random
import sys
import shutil

import numpy as np
from ase import Atoms
from ase.db import connect
from ase.io import read, write
from ase.vibrations import Vibrations
from sella import Sella

from kinbot.constants import EVtoHARTREE
#from kinbot.ase_modules.calculators.{code} import {Code}
from fairchem.core import pretrained_mlip, FAIRChemCalculator
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
        vib.write_jmol()
        # Use kinbot frequencies to avoid mixing low vib frequencies with 
        # the values associated with external rotations.
        _ = vib.get_frequencies()
        zpe = vib.get_zero_point_energy() * EVtoHARTREE
        hessian = vib.H / 97.17370087
        st_pt = StationaryPoint.from_ase_atoms(mol)
        st_pt.characterize()
        freqs, _ = get_frequencies(st_pt, hessian, st_pt.geom)
        os.chdir(init_dir)
        #shutil.rmtree('{label}_vib')
        return freqs, zpe, hessian

with open('fairchem.log', 'a') as f:
    f.write('{label} | Beginning well optimization...\n')

db = connect('{working_dir}/kinbot.db')
mol = Atoms(symbols={atom}, 
            positions={geom})

kwargs = {kwargs}
mol.info.update({{"charge": kwargs['charge'], "spin": kwargs['mult']}})

mol.calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit("uma-s-1", device="cpu"), task_name="omol")
with open('fairchem.log', 'a') as f:
    f.write('{label} | Initial energy calculated\n')

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

# For monoatomic wells, just calculate the energy and exit. 
if len(mol) == 1:
    with open('fairchem.log', 'a') as f:
        f.write('{label} | Detected monatomic well\n')

    e = mol.get_potential_energy()
    forces = mol.calc.results['forces']
    del mol.calc.results['forces']
    random.seed()
    db.write(mol, name='{label}',
             data={{'energy': e, 'frequencies': np.array([]), 'zpe': 0.0,
                 'hess': np.zeros([3, 3]), 'status': 'normal'}})
    with open('{label}.log', 'a') as f:
        f.write('done\n')
    with open('fairchem.log', 'a') as f:
        f.write('{label} | Well optimization successful\n')
    sys.exit(0)

order = {order}
sella_kwargs = {sella_kwargs}
opt = Sella(mol, 
            order=order, 
            trajectory='{label}.traj', 
            logfile='{label}_sella.log',
            **sella_kwargs)

freqs = []
converged = False
fmax = 1e-4
attempts = 1
steps=500
while not converged and attempts <= 3:
    mol.calc.label = '{label}'
    with open('fairchem.log', 'a') as f:
        f.write(f'{label} | Optimizing well. Attempt {{attempts}}\n')
    try:
        converged = opt.run(fmax=fmax, steps=steps)
        traj = read('{label}.traj', index=':')
        write('{label}.xyz', traj, format='xyz')
    except ValueError:
        with open('fairchem.log', 'a') as f:
            f.write(f'{label} | Optimization failed. Perturbing coordinates\n')
        mol.set_positions(mol.get_positions() + np.random.normal(scale=0.05, size=(len(mol), 3)))
        opt = Sella(mol,
            order=order,
            trajectory='{label}.traj',
            logfile='{label}_sella.log',
            **sella_kwargs)
        converged = opt.run(fmax=fmax, steps=steps)
        traj = read('{label}.traj', index=':')
        write('{label}.xyz', traj, format='xyz')

    freqs, zpe, hessian = calc_vibrations(mol)
    if order == 0 and (np.count_nonzero(np.array(freqs) < 0) > 1
                        or np.count_nonzero(np.array(freqs) < -50) >= 1):
        print(f'Found one or more imaginary frequencies. {{freqs[1:6]}}')
        converged = False
        mol.calc.label = '{label}'
        attempts += 1
        fmax *= 0.3
        if attempts <= 3:
            print(f'Retrying with a tighter criterion: fmax={{fmax}}.')
    elif order == 1 and (np.count_nonzero(np.array(freqs) < 0) > 2  # More than two imag frequencies
                         or np.count_nonzero(np.array(freqs) < -50) >= 2  # More than one imag frequency larger than 50i
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
        forces = mol.calc.results['forces']
        del mol.calc.results['forces']
        random.seed()
        db.write(mol, name='{label}', 
                    data={{'energy': e, 'frequencies': freqs, 'zpe': zpe, 
                    'hess': hessian, 'forces': forces, 'status': 'normal'}})
        with open('fairchem.log', 'a') as f:
            f.write('{label} | Well optimization successful!\n')
if not converged:
    with open('fairchem.log', 'a') as f:
        f.write('{label} | Well did not converge\n')
    raise RuntimeError("Did not converge")

with open('{label}.log', 'a') as f:
    f.write('done\n')
