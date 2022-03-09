from ase import Atoms
from ase.calculators.qchem import QChem
from ase.db import connect
from kinbot import reader_qchem

db = connect('{working_dir}/kinbot.db')

dummy = {dummy}
mol = Atoms(symbols={atom}, positions={geom})

kwargs = {kwargs}
calc = QChem(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy() # use the QChem optimizer
    for d in dummy:
        mol.pop()
    db.write(mol, name='{label}', data={{'energy': e,'status': 'normal'}})
except RuntimeError: 
    try:
        mol.positions = reader_qchem.read_geom('{label}.out', mol, dummy)
        e = mol.get_potential_energy() # use the Gaussian optimizer
        for d in dummy:
            mol.pop()
        db.write(mol, name='{label}', data={{'energy': e, 'status': 'normal'}})
    except:
        db.write(mol, name='{label}', data={{'status': 'error'}})

with open('{label}.out', 'a') as f:
    f.write('done\n')
