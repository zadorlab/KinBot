import numpy as np
from ase import Atoms
from ase.calculators.qchem import QChem
from ase.db import connect

from kinbot import reader_qchem
from kinbot.utils import iowait

db = connect('{working_dir}/kinbot.db')
label = '{label}'
logfile = '{label}.out'

mol = Atoms(symbols={atom}, positions={geom})

# Perform Geometry optimization
kwargs = {kwargs}
QChem.command = '{qc_command} -nt {ppn} PREFIX.in PREFIX.out PREFIX.sv'
calc = QChem(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy()  # It's actually peforming an optimization
    iowait(logfile, 'qchem')
except RuntimeError:
    try:
        # Try again with the last geometry just in case.
        mol.positions = reader_qchem.read_geom(logfile, mol)
        e = mol.get_potential_energy()
        iowait(logfile, 'qchem')
    except RuntimeError:
        db.write(mol, name=label, data={{'status': 'error'}})


# Perform a vibrational analysis on the optimized geometry
mol.positions = reader_qchem.read_geom(logfile, mol)
kwargs['jobtype'] = 'freq'
kwargs['xc_grid'] = '3'
kwargs['label'] = kwargs['label'] + '_freq'
calc = QChem(**kwargs)
mol.set_calculator(calc)

try:
    e = mol.get_potential_energy()  # Compute frequencies
    iowait(logfile, 'qchem')
except RuntimeError:
    db.write(mol, name=label, data={{'status': 'error'}})

freq = reader_qchem.read_freq('{label}_freq.out', {atom})
zpe = reader_qchem.read_zpe('{label}_freq.out')

db.write(mol, name=label, data={{'energy': e, 'frequencies': np.asarray(freq),
                                 'zpe': zpe, 'status': 'normal'}})

with open(logfile, 'a') as f:
    f.write('done\n')

with open('{label}_freq.out', 'a') as f:
    f.write('done\n')
