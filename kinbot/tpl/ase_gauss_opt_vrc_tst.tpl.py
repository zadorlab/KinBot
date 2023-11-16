import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

def correct_kwargs(kwargs, iteration):
    """Funtion that refine the optimization parameters depending on the number of trials"""
    match iteration:
        case 0:
            kwargs["iop"] = "1/8=2"
            kwargs['opt'] = kwargs['opt'].replace('CalcFC', 'ReCalcFC=5')
            kwargs['opt'] = kwargs['opt'].replace("MaxCycles=30", "MaxCycles=70")
        case 1:
            kwargs["iop"] = "1/8=2" #Use smaller max step
            kwargs['opt'] = kwargs['opt'].replace('ReCalcFC=5', 'CalcAll')
        case 2:
            kwargs["iop"] = "1/8=1,1/19=10" #Use different optimization procedure
            kwargs['opt'] = kwargs['opt'].replace("MaxCycles=70", "MaxCycles=50")
    return kwargs

db = connect('{working_dir}/kinbot.db')
label = '{label}'
logfile = '{label}.log'

mol = Atoms(symbols={atom}, positions={geom})

kwargs = {kwargs}
Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy()  # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    freq = reader_gauss.read_freq(logfile, {atom})
    zpe = reader_gauss.read_zpe(logfile)
    db.write(mol, name=label, data={{'energy': e, 'frequencies': np.asarray(freq),
                                     'zpe': zpe, 'status': 'normal'}})

except RuntimeError:
    for i in range(3):
        try:
            iowait(logfile, 'gauss')
            e, mol.positions = reader_gauss.read_lowest_geom_energy(logfile, mol)
            kwargs = correct_kwargs(kwargs, i)
            mol.calc = Gaussian(**kwargs)
            e = mol.get_potential_energy()  # use the Gaussian optimizer
            iowait(logfile, 'gauss')
            mol.positions = reader_gauss.read_geom(logfile, mol)
            freq = reader_gauss.read_freq(logfile, {atom})
            zpe = reader_gauss.read_zpe(logfile)
            db.write(mol, name=label, data={{'energy': e,
                                             'frequencies': np.asarray(freq),
                                             'zpe': zpe, 'status': 'normal'}})
            break
        except RuntimeError:
            iowait(logfile, 'gauss')
            #Save in db the lowest energy geometry if forces are converged
            if reader_gauss.read_convergence(logfile) != 0:
                e, mol.positions = reader_gauss.read_converged_geom_energy(logfile, mol)
                freq = reader_gauss.read_freq(logfile, {atom})
                zpe = reader_gauss.read_zpe(logfile)
                db.write(mol, name=label, data={{'energy': e,
                                            'frequencies': np.asarray(freq),
                                            'zpe': zpe, 'status': 'normal'}})
                break
            if i == 2:
                mol.positions = reader_gauss.read_geom(logfile, mol)
                freq = reader_gauss.read_freq(logfile, {atom})
                zpe = reader_gauss.read_zpe(logfile)
                db.write(mol, name=label, data={{'status': 'error'}})
        else:
            break

with open(logfile, 'a') as f:
    f.write('done\n')
