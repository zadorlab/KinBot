from ase import Atoms
from ase.calculators.qchem import QChem
from ase.db import connect
from kinbot import reader_qchem

db = connect('{working_dir}/kinbot.db')

mol = Atoms(symbols={atom}, positions={geom})

QChem.command = '{qc_command} -nt {ppn} PREFIX.in PREFIX.out PREFIX.sv'
kwargs = {kwargs}
calc = QChem(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy()  # use the QChem optimizer
    mol.positions = reader_qchem.read_geom('{label}.out', mol)
    db.write(mol, name='{label}', data={{'energy': e, 'status': 'normal'}})
except RuntimeError: 
    try:
        mol.positions = reader_qchem.read_geom('{label}.out', mol)
        e = mol.get_potential_energy()  # use the QChem optimizer
        mol.positions = reader_qchem.read_geom('{label}.out', mol)
        db.write(mol, name='{label}', data={{'energy': e, 'status': 'normal'}})
    except:
        db.write(mol, name='{label}', data={{'status': 'error'}})

with open('{label}.out', 'a') as f:
    f.write('done\n')
