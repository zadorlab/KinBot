import re
import numpy as np
import ase
from ase import Atoms
from ase.calculators.gaussian import Gaussian
from ase.db import connect
from kinbot import reader_gauss

db = connect('{working_dir}/kinbot.db')
label = '{label}'
logfile = '{label}.log'

mol = Atoms(symbols={atom}, positions={geom})

kwargs = {kwargs}
Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.set_calculator(calc)

try:
    e = mol.get_potential_energy() # use the Gaussian optimizer
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, {atom})
    if freq[0] < 0. and freq[0] > -50.:  
        kwargs['opt'] = kwargs['opt'].replace('CalcFC', 'CalcAll')
        try:
            del kwargs['freq']
        except:
            pass
        e = mol.get_potential_energy() 
        mol.positions = reader_gauss.read_geom(logfile, mol)
        freq = reader_gauss.read_freq(logfile, {atom})
    zpe = reader_gauss.read_zpe(logfile)
    db.write(mol, name=label, data={{'energy': e,'frequencies': np.asarray(freq), 'zpe':zpe, 'status': 'normal'}})

except RuntimeError: 
    try:
        mol.positions = reader_gauss.read_geom(logfile, mol)
        e = mol.get_potential_energy() # use the Gaussian optimizer
        mol.positions = reader_gauss.read_geom(logfile, mol)
        freq = reader_gauss.read_freq(logfile, {atom})
        zpe = reader_gauss.read_zpe(logfile)
        db.write(mol, name=label, data={{'energy': e,'frequencies': np.asarray(freq), 'zpe':zpe, 'status': 'normal'}})
    except:
        db.write(mol, name=label, data={{'status': 'error'}})

with open(logfile,'a') as f:
    f.write('done\n')
