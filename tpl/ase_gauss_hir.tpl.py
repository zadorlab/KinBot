import ase
from ase import Atoms
from ase.calculators.gaussian import Gaussian
from ase.db import connect
from kinbot import reader_gauss

db = connect('{working_dir}/kinbot.db')

mol = Atoms(symbols={atom}, positions={geom})

Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
kwargs = {kwargs}
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy() # use the Gaussian optimizer
    mol.positions = reader_gauss.read_geom('{label}.log', mol)
    db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
except RuntimeError: 
    try:
        mol.positions = reader_gauss.read_geom('{label}.log', mol)
        e = mol.get_potential_energy() # use the Gaussian optimizer
        mol.positions = reader_gauss.read_geom('{label}.log', mol)
        db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
    except:
        db.write(mol, name='{label}', data={{'status' : 'error'}})

with open('{label}.log', 'a') as f:
    f.write('done\n')
