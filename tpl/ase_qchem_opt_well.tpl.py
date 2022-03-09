import numpy as np

from ase import Atoms
from ase.calculators.qchem import QChem
from ase.db import connect
from kinbot import reader_qchem

db = connect('{working_dir}/kinbot.db')

dummy = {dummy}
mol = Atoms(symbols={atom}, positions={geom})

# Perform Geometry optimization
kwargs = {kwargs}
QChem.command = '{qc_command} -nt {ppn} PREFIX.in PREFIX.out PREFIX.sv'
calc = QChem(**kwargs)
mol.calc = calc
try:
    e = mol.get_potential_energy()  # It's actually peforming an optimization
except RuntimeError:
    # Try again with the last geometry just in case.
    mol = reader_qchem.read_geom('{label}.out', dummy)
    mol.set_calculator(calc)
    e = mol.get_potential_energy()

# Perform vibrational analysis on the optimized geometry
opt_mol = reader_qchem.read_geom('{label}.out', dummy)
kwargs['jobtype'] = 'freq'
kwargs['label'] = kwargs['label'] + '_freq'
calc = QChem(**kwargs)
opt_mol.set_calculator(calc)

try:
    e = opt_mol.get_potential_energy()  # Compute frequencies

except RuntimeError:
    try:
        # Try again with the last geometry just in case.
        opt_mol = reader_qchem.read_geom('{label}.out', dummy)
        opt_mol.set_calculator(calc)
        e = opt_mol.get_potential_energy()  # Compute frequencies
    except RuntimeError:
        db.write(mol, name='{label}', data={{'status': 'error'}})

freq = reader_qchem.read_freq('{label}_freq.out', {atom})
zpe = reader_qchem.read_zpe('{label}_freq.out')
for d in dummy:
    mol.pop()

db.write(mol, name='{label}', data={{'energy': e,
                                     'frequencies': np.asarray(freq),
                                     'zpe': zpe, 'status': 'normal'}})

with open('{label}.out', 'a') as f:
    f.write('done\n')
