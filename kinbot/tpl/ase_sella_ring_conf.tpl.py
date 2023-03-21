import os
import sys
import shutil

import numpy as np
from ase import Atoms
from ase.db import connect
from ase.vibrations import Vibrations
from sella import Sella, Constraints

from kinbot.constants import EVtoHARTREE
from kinbot.modify_geom import modify_coordinates
from kinbot.ase_modules.calculators.{code} import {Code}
from kinbot.stationary_pt import StationaryPoint
from kinbot.frequencies import get_frequencies

db = connect('{working_dir}/kinbot.db')
mol = Atoms(symbols={atom}, 
            positions={geom})

kwargs = {kwargs}
mol.calc = {Code}(**kwargs)

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
        raise ValueError(f'¯\_(ツ)_/¯, Unexpected length of fix: {{fix}}')

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
        raise ValueError(f'¯\_(ツ)_/¯, Unexpected length of changes: {{change}}')

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')

# For monoatomic wells, just calculate the energy and exit. 
if len(mol) == 1:
    e = mol.get_potential_energy()
    db.write(mol, name='{label}',
             data={{'energy': e, 'frequencies': np.array([]), 'zpe': 0.0,
                    'hess': np.zeros([3, 3]), 'status': 'normal'}})
    with open('{label}.log', 'a') as f:
        f.write('done\n')
    sys.exit(0)

opt  = Sella(mol, order={order}, trajectory='{label}.traj', 
             logfile='{label}_sella.log')
try:
    cvgd = opt.run(fmax=0.0001, steps=300)
    if cvgd:
        e = mol.get_potential_energy()
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
        # Use kinbot frequencies to avoid problems with frequencies associated 
        # with external rotations.
        _ = vib.get_frequencies()
        zpe = vib.get_zero_point_energy() * EVtoHARTREE
        hessian = vib.H / 97.17370087
        st_pt = StationaryPoint.from_ase_atoms(mol)
        st_pt.characterize()
        freqs, __ = get_frequencies(st_pt, hessian, st_pt.geom)
        os.chdir(init_dir)
        shutil.rmtree('{label}_vib')
        db.write(mol, name='{label}', 
                 data={{'energy': e, 'frequencies': freqs, 'zpe': zpe, 
                        'hess': hessian, 'status': 'normal'}})
    else:  # TODO Eventually we might want to correct something in case it fails.
        raise RuntimeError
except RuntimeError:
    db.write(mol, name='{label}', data={{'status': 'error'}})

with open('{label}.log', 'a') as f:
    f.write('done\n')
