"""Frequency-only recovery at the selected geometry; no optimization."""
import numpy as np
from ase import Atoms
from ase.calculators.qchem import QChem
from ase.db import connect
from kinbot import reader_qchem
from kinbot.utils import iowait

db = connect('{working_dir}/kinbot.db')
label = '{label}'
mol = Atoms(symbols={atom}, positions={geom})
calc = QChem(**{kwargs})
calc.command = '{qc_command} -nt {ppn} PREFIX.inp PREFIX.out'
mol.calc = calc
try:
    energy = mol.get_potential_energy()
    iowait(label + '_freq.out', 'qchem')
    freq = reader_qchem.read_freq(label + '_freq.out', {atom})
    zpe = reader_qchem.read_zpe(label + '_freq.out')
    data = {{'energy': energy, 'frequencies': np.asarray(freq),
             'zpe': zpe, 'status': 'normal'}}
except RuntimeError:
    data = {{'status': 'error'}}
db.write(mol, name=label, data=data)
with open(label + '.out', 'a') as handle:
    handle.write('done\n')
