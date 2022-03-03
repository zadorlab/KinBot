from math import pi
import numpy as np
import ase
from ase import Atoms
from ase.calculators.gaussian import Gaussian
from ase.db import connect
from kinbot import reader_gauss

db = connect('{working_dir}/kinbot.db')

dummy = None
scan = {scan}
bimol = {bimol}
mol = Atoms(symbols={atom}, positions={geom})

kwargs = {kwargs}
Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.set_calculator(calc)

for tr in range({ntrial}):
    try:
        success = True
        e = mol.get_potential_energy() # use the Gaussian optimizer (task optimize)
        mol.positions = reader_gauss.read_geom('{label}.log', mol, dummy)
        db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
        break
    except RuntimeError: 
        success = False
        
if not success:
    if not bimol:
        try:
            mol.positions = reader_gauss.read_geom('{label}.log', mol, dummy)
            del kwargs['opt']  # this is when we give up optimization!!
            calc = Gaussian(**kwargs)
            e = mol.get_potential_energy() 
            mol.positions = reader_gauss.read_geom('{label}.log', mol, dummy)
            db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
        except: 
            if scan == 0:
                db.write(mol, name = '{label}', data = {{'status': 'error'}})
            elif scan == 1:
                # exception for scan-type calculations
                # write final geometry and energy even if all tries failed
                mol.positions = reader_gauss.read_geom('{label}.log', mol, dummy)
                e = reader_gauss.read_energy('{label}.log')
                mol.positions = reader_gauss.read_geom('{label}.log', mol, dummy)
                if mol.positions is not None and e is not None: 
                    db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
                else:
                    db.write(mol, name='{label}', data={{'status': 'error'}})
    else:
        try:
            mol.positions = reader_gauss.read_geom('{label}.log', mol, dummy)
            db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
        except: 
            db.write(mol, name = '{label}', data = {{'status': 'error'}})

with open(f'{label}.log','a') as f:
    f.write('done\n')
