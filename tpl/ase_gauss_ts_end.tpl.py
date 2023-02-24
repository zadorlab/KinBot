import numpy as np
from ase import Atoms
#from ase.calculators.gaussian import Gaussian
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian  # New
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('{working_dir}/kinbot.db')
label = '{label}'
logfile = '{label}.log'

atom = {atom}
geom = {geom}
mol = Atoms(symbols=atom, positions=geom)

kwargs = {kwargs}
Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy() # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, {atom})
    zpe = reader_gauss.read_zpe(logfile)
    db.write(mol, name=label, data={{'energy': e,'frequencies': np.asarray(freq), 'zpe':zpe, 'status': 'normal'}})
except RuntimeError:
    try:
        iowait(logfile, 'gauss')
        mol.positions = reader_gauss.read_geom(logfile, mol)
        kwargs = reader_gauss.correct_kwargs(logfile, kwargs)
        mol.calc = Gaussian(**kwargs)
        e = mol.get_potential_energy()  # use the Gaussian optimizer
        iowait(logfile, 'gauss')
        mol.positions = reader_gauss.read_geom(logfile, mol)
        freq = reader_gauss.read_freq(logfile, {atom})
        zpe = reader_gauss.read_zpe(logfile)
        db.write(mol, name=label, data={{'energy': e,
                                         'frequencies': np.asarray(freq),
                                         'zpe': zpe, 'status': 'normal'}})
    except RuntimeError:
        db.write(mol, name=label, data={{'status': 'error'}})

with open(logfile,'a') as f:
    f.write('done\n')
