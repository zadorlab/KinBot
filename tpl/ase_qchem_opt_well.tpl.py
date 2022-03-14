import numpy as np

from ase import Atoms
from ase.calculators.qchem import QChem
from ase.db import connect
from kinbot import reader_qchem

db = connect('{working_dir}/kinbot.db')
label = '{label}'
outfile = '{label}.out'

mol = Atoms(symbols={atom}, positions={geom})

# Perform Geometry optimization
kwargs = {kwargs}
QChem.command = '{qc_command} -nt {ppn} PREFIX.in PREFIX.out PREFIX.sv'
calc = QChem(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy()  # It's actually peforming an optimization
except RuntimeError:
    try:
        # Try again with the last geometry just in case.
        mol.positions = reader_qchem.read_geom('{label}.out', mol)
        e = mol.get_potential_energy()
    except RuntimeError:
        db.write(mol, name=label, data={{'status': 'error'}})

# Perform vibrational analysis on the optimized geometry


try:
    opt_mol = reader_qchem.read_geom('{label}.out', mol)
    kwargs['jobtype'] = 'freq'
    kwargs['label'] = kwargs['label'] + '_freq'
    calc = QChem(**kwargs)
    opt_mol.set_calculator(calc)
    e = opt_mol.get_potential_energy()  # Compute frequencies
except RuntimeError:
    db.write(mol, name='{label}', data={{'status': 'error'}})

freq = reader_qchem.read_freq('{label}_freq.out', {atom})
zpe = reader_qchem.read_zpe('{label}_freq.out')

db.write(mol, name='{label}', data={{'energy': e,
                                     'frequencies': np.asarray(freq),
                                     'zpe': zpe, 'status': 'normal'}})

with open('{label}.out', 'a') as f:
    f.write('done\n')
