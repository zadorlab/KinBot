import numpy as np
from ase import Atoms
from ase.calculators.qchem import QChem
from ase.db import connect
from kinbot import reader_qchem

db = connect('{working_dir}/kinbot.db')

atom = {atom}
geom = {geom}
mol = Atoms(symbols=atom, positions=geom)

kwargs = {kwargs}
QChem.command = '{qc_command} -nt {ppn} PREFIX.in PREFIX.out PREFIX.sv'
calc = QChem(**kwargs)
mol.set_calculator(calc)

# Perform TS optimization.
try:
    e = mol.get_potential_energy()  # use the QChem optimizer
except RuntimeError:
    try:
        mol.positions = reader_qchem.read_geom('{label}.out', mol)
        e = mol.get_potential_energy()  # use the QChem optimizer
    except RuntimeError:
        db.write(mol, name='{label}', data={{'status': 'error'}})

# Compute Frequencies
mol.positions = reader_qchem.read_geom('{label}.out', mol)
kwargs['jobtype'] = 'freq'
calc = QChem(**kwargs)
mol.set_calculator(calc)
try:
    e = mol.get_potential_energy()  # use the QChem optimizer
except RuntimeError:
    try:
        mol.positions = reader_qchem.read_geom('{label}.out', mol)
        e = mol.get_potential_energy()  # use the QChem optimizer
    except RuntimeError:
        db.write(mol, name='{label}', data={{'status': 'error'}})

freq = reader_qchem.read_freq('{label}.out', {atom})
zpe = reader_qchem.read_zpe('{label}.out')
db.write(mol, name='{label}', data={{'energy': e, 'frequencies': np.asarray(freq),
                                     'zpe': zpe, 'status': 'normal'}})

with open(f'{label}.out', 'a') as f:
    f.write('done\n')
