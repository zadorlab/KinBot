import ase
from ase import Atoms
from ase.calculators.gaussian import Gaussian
from ase.db import connect
from kinbot import reader_gauss

db = connect('{working_dir}/kinbot.db')

dummy = {dummy}
mol = Atoms(symbols={atom}, positions={geom})

Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
kwargs = {kwargs}
calc = Gaussian(**kwargs)
mol.set_calculator(calc)

try:
    e = mol.get_potential_energy() # use the Gaussian optimizer
    #Positions (geom) updated in ase/ases/io/gaussian.py code
    for d in dummy:
        mol.pop()
    db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
except RuntimeError: 
    try:
        mol.positions = read_geom('{label}.log', mol, dummy)
        e = mol.get_potential_energy() # use the Gaussian optimizer
        for d in dummy:
            mol.pop()
        db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
    except:
        db.write(mol, name='{label}', data={{'status' : 'error'}})

with open('{label}.log', 'a')
    f.write('done\n')
