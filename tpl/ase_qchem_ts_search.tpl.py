from math import pi
import numpy as np
import ase
from ase import Atoms
from ase.calculators.qchem import QChem
from ase.db import connect
from kinbot import reader_qchem

db = connect('{working_dir}/kinbot.db')

dummy = None
scan = {scan}
mol = Atoms(symbols={atom}, positions={geom})

kwargs = {kwargs}
QChem.command = '{qc_command} PREFIX.in PREFIX.out'
calc = QChem(**kwargs)
mol.set_calculator(calc)

try:
    e = mol.get_potential_energy() # use the Gaussian optimizer (task optimize)
    db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
except RuntimeError: 
    try:
        mol.positions = reader_gauss.read_geom('{label}.out', mol, dummy)
        e = mol.get_potential_energy() # use the Gaussian optimizer
        db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
    except:
        try:
            mol.positions = reader_gauss.read_geom('{label}.out', mol, dummy)
            e = mol.get_potential_energy() # use the Gaussian optimizer
            db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
        except:
            try:
                mol.positions = reader_gauss.read_geom('{label}.out', mol, dummy)
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
                    mol.positions = reader_gauss.read_geom('{label}.out', mol, dummy)
                    e = reader_gauss.read_energy('{label}.out')
                    if mol.positions is not None and e is not None: 
                        db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
                    else:
                        db.write(mol, name='{label}', data={{'status': 'error'}})

with open(f'{label}.out','a') as f:
    f.write('done\n')
