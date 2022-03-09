import numpy as np
import ase
from ase import Atoms
from ase.calculators.qchem import QChem
from ase.db import connect
from kinbot import reader_qchem

db = connect('{working_dir}/kinbot.db')

atom = {atom}
geom = {geom}
dummy = None
mol = Atoms(symbols=atom, positions=geom)

kwargs = {kwargs}
QChem.command = '{qc_command} < PREFIX.in > PREFIX.out'
calc = QChem(**kwargs)
mol.set_calculator(calc)

try:
    e = mol.get_potential_energy() # use the Gaussian optimizer
    freq = reader_qchem.read_freq('{label}.out', {atom})
    zpe = reader_qchem.read_zpe('{label}.out')
    db.write(mol, name='{label}', data={{'energy': e,'frequencies': np.asarray(freq), 'zpe':zpe, 'status': 'normal'}})
except RuntimeError: 
    try:
        mol.positions = reader_qchem.read_geom('{label}.out', mol, dummy)
        e = mol.get_potential_energy() # use the Gaussian optimizer
        freq = reader_qchem.read_freq('{label}.out', {atom})
        zpe = reader_qchem.read_zpe('{label}.out')
        db.write(mol, name='{label}', data={{'energy': e,'frequencies': np.asarray(freq), 'zpe':zpe, 'status': 'normal'}})
    except:
        db.write(mol, name='{label}', data = {{'status': 'error'}})

with open(f'{label}.out','a') as f:
    f.write('done\n')
