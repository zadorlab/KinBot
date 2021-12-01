import numpy as np
import ase
from ase import Atoms
from ase.calculators.gaussian import Gaussian
from ase.db import connect
from kinbot import reader_gauss

db = connect('{working_dir}/kinbot.db')

atom = {atom}
dummy = {dummy}
mol = Atoms(symbols=atom, positions={geom})

Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
kwargs = {kwargs}
calc = Gaussian(**kwargs)
mol.set_calculator(calc)

try:
    e = mol.get_potential_energy() # use the Gaussian optimizer (task optimize)
    freq = reader_gauss.read_freq('{label}.log', atom)
    zpe = reader_gauss.read_zpe('{label}.log')
    for d in dummy:
        mol.pop()
    db.write(mol, name='{label}',data={{'energy': e,'frequencies': np.asarray(freq), 'zpe':zpe, 'status': 'normal'}})
except RuntimeError: 
    db.write(mol, name='{label}', data={{'status': 'error'}})

with open('{label}.log', 'a')
    f.write('done\n')
