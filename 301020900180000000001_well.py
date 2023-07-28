import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('/Users/csoulie/Documents/Research/KinBot_Dev/KinBot/kinbot.db')
label = '301020900180000000001_well'
logfile = '301020900180000000001_well.log'

mol = Atoms(symbols=['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'], positions=[[0.27674, 0.60575, -0.0], [1.7959, 0.60575, -0.0], [-0.09983, -0.0876, -0.7813], [-0.09983, 1.62904, -0.20981], [-0.09983, 0.2758, 0.9911], [2.17247, 0.9357, -0.9911], [2.17247, -0.41754, 0.20981], [2.17247, 1.2991, 0.7813]])

kwargs = {'method': 'b3lyp', 'basis': '6-31+g', 'nprocshared': 8, 'mem': '700MW', 'chk': '301020900180000000001_well', 'label': '301020900180000000001_well', 'Symm': 'None', 'mult': 1, 'charge': 0, 'scf': 'xqc', 'freq': 'freq', 'opt': 'CalcFC, Tight'}
Gaussian.command = 'g16 < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy()  # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, ['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'])
    zpe = reader_gauss.read_zpe(logfile)
    db.write(mol, name=label, data={'energy': e, 'frequencies': np.asarray(freq),
                                     'zpe': zpe, 'status': 'normal'})

except RuntimeError:
    for i in range(3):
        try:
            iowait(logfile, 'gauss')
            mol.positions = reader_gauss.read_geom(logfile, mol)
            kwargs = reader_gauss.correct_kwargs(logfile, kwargs)
            mol.calc = Gaussian(**kwargs)
            e = mol.get_potential_energy()  # use the Gaussian optimizer
            iowait(logfile, 'gauss')
            mol.positions = reader_gauss.read_geom(logfile, mol)
            freq = reader_gauss.read_freq(logfile, ['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H'])
            zpe = reader_gauss.read_zpe(logfile)
            db.write(mol, name=label, data={'energy': e,
                                             'frequencies': np.asarray(freq),
                                             'zpe': zpe, 'status': 'normal'})
        except RuntimeError:
            if i == 2:
                db.write(mol, name=label, data={'status': 'error'})
            pass
        else:
            break

with open(logfile, 'a') as f:
    f.write('done\n')
