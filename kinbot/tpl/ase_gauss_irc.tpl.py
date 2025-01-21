import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('{working_dir}/kinbot.db')
label = '{label}'
logfile = '{label}.log'

mol = Atoms(symbols={atom}, positions={geom})

kwargs = {kwargs}
Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

success = True

try:
    e = mol.get_potential_energy()  # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    db.write(mol, name=label, data={{'energy': e, 'status': 'normal'}})
except RuntimeError:
    # Retry by correcting errors
    try:
        iowait(logfile, 'gauss')
        mol.positions = reader_gauss.read_geom(logfile, mol)
        kwargs = reader_gauss.correct_kwargs(logfile, kwargs)
        mol.calc = Gaussian(**kwargs)
        e = mol.get_potential_energy()  # use the Gaussian optimizer
        iowait(logfile, 'gauss')
        mol.positions = reader_gauss.read_geom(logfile, mol)
        db.write(mol, name=label, data={{'energy': e, 'status': 'normal'}})
    except RuntimeError:
        if mol.positions is not None:
            # although there is an error, continue from the final geometry
            db.write(mol, name=label, data={{'status': 'normal'}})
        else:
            db.write(mol, name=label, data={{'status': 'error'}})
            success = False

with open(logfile, 'a') as f:
    f.write('done\n')

if success:
    label = '{label}_prod'
    logfile = '{label}_prod.log'
    # start the product optimization
    prod_kwargs = {prod_kwargs}
    calc_prod = Gaussian(**prod_kwargs)
    mol_prod = Atoms(symbols={atom}, positions=mol.positions)
    mol_prod.calc = calc_prod
    try:
        e = mol_prod.get_potential_energy() # use the Gaussian optimizer
        iowait(logfile, 'gauss')
        mol_prod.positions = reader_gauss.read_geom(logfile, 
                                                    mol_prod, 
                                                    max2frag=True, 
                                                    charge=kwargs['charge'],
                                                    mult=kwargs['mult'])
        freq = reader_gauss.read_freq(logfile, {atom})
        zpe = reader_gauss.read_zpe(logfile)
        db.write(mol_prod, name=label, data={{'energy': e,
                                         'frequencies': np.asarray(freq),
                                         'zpe': zpe, 'status': 'normal'}})
    except RuntimeError: 
        for i in range(3):
            try:
                iowait(logfile, 'gauss')
                _, mol_prod.positions = reader_gauss.read_lowest_geom_energy(logfile, mol_prod)
                prod_kwargs = reader_gauss.correct_kwargs(logfile, prod_kwargs)
                mol_prod.calc = Gaussian(**prod_kwargs)
                e = mol_prod.get_potential_energy()  # use the Gaussian optimizer
                iowait(logfile, 'gauss')
                mol_prod.positions = reader_gauss.read_geom(logfile, mol_prod)
                freq = reader_gauss.read_freq(logfile, {atom})
                zpe = reader_gauss.read_zpe(logfile)
                db.write(mol_prod, name=label, data={{'energy': e,
                                                'frequencies': np.asarray(freq),
                                                'zpe': zpe, 'status': 'normal'}})
            except RuntimeError:
                if i == 2:
                    db.write(mol, name=label, data={{'status': 'error'}})
                pass
            else:
                break

    with open(logfile, 'a') as f:
        f.write('done\n')
