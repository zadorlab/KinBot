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
mol = Atoms(symbols={atom}, positions={geom})

kwargs = {kwargs}
Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.set_calculator(calc)

try:
    e = mol.get_potential_energy() # use the Gaussian optimizer (task optimize)
    db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
except RuntimeError: 
    try:
        mol.positions = reader_gauss.read_geom('{label}.log', mol, dummy)
        e = mol.get_potential_energy() # use the Gaussian optimizer
        db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
    except:
        try:
            mol.positions = reader_gauss.read_geom('{label}.log', mol, dummy)
            e = mol.get_potential_energy() # use the Gaussian optimizer
            db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
        except:
            try:
                mol.positions = reader_gauss.read_geom('{label}.log', mol, dummy)
                del kwargs['opt']  # this is when we give up optimization!!
                calc = Gaussian(**kwargs)
                e = mol.get_potential_energy() 
                db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
            except: 
                if scan == 0:
                    db.write(mol, name = '{label}', data = {{'status': 'error'}})
                elif scan == 1:
                    # exception for scan-type calculations
                    # write final geometry and energy even if all tries failed
                    mol.positions = reader_gauss.read_geom('{label}.log', mol, dummy)
                    e = reader_gauss.read_energy('{label}.log')
                    if mol.positions is not None and e is not None: 
                        db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
                    else:
                        db.write(mol, name='{label}', data={{'status': 'error'}})

with open(f'{label}.log','a') as f:
    f.write('done\n')
