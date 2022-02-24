import numpy as np

from ase import Atoms
from ase.calculators.qchem import QChem
from ase.db import connect
from kinbot import reader_qchem

db = connect('{working_dir}/kinbot.db')

dummy = {dummy}
mol = Atoms(symbols={atom}, positions={geom})

kwargs = {kwargs}
QChem.command = '{qc_command} -nt {ppn} PREFIX.in PREFIX.out PREFIX.sv'
calc = QChem(**kwargs)
mol.set_calculator(calc)

try:
    e = mol.get_potential_energy()  # use the QChem optimizer
    freq = reader_qchem.read_freq('{label}.log', {atom})
    zpe = reader_qchem.read_zpe('{label}.log')
    for d in dummy:
        mol.pop()
    db.write(mol, name='{label}', data={
        {'energy': e, 'frequencies': np.asarray(freq), 'zpe': zpe,
         'status': 'normal'}})

except RuntimeError:
    try:
        mol.positions = reader_qchem.read_geom('{label}.log', mol, dummy)
        e = mol.get_potential_energy()  # use the QChem optimizer
        freq = reader_qchem.read_freq('{label}.log', {atom})
        zpe = reader_qchem.read_zpe('{label}.log')
        for d in dummy:
            mol.pop()
        db.write(mol, name='{label}', data={
            {'energy': e, 'frequencies': np.asarray(freq), 'zpe': zpe,
             'status': 'normal'}})
    except:
        db.write(mol, name='{label}', data={{'status': 'error'}})

with open(f'{label}.log', 'a') as f:
    f.write('done\n')
