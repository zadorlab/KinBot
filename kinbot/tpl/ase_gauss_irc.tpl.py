import numpy as np
import logging
from pathlib import Path
from shutil import copyfile
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

if {bimol}:
    success = False
    initial_geometry = mol.positions.copy()
    for attempt in range(2):
        completed = False
        try:
            e = mol.get_potential_energy()
            completed = True
        except RuntimeError:
            pass
        iowait(logfile, 'gauss')
        _, accepted, _, _ = reader_gauss.read_all_irc_geoms(logfile)
        if accepted is not None and accepted.shape[1:] == initial_geometry.shape:
            mol.positions = accepted[-1]
            data = {{'status': 'normal', 'irc_completed': completed,
                     'irc_accepted_points': len(accepted)}}
            if completed:
                data['energy'] = e
            else:
                logging.warning('%s: using the last of %s accepted IRC points after early termination.',
                                label, len(accepted))
            db.write(mol, name=label, data=data)
            success = True
            break
        if attempt == 0:
            for suffix in ('.log', '.com', '.chk'):
                if Path(label + suffix).is_file():
                    copyfile(label + suffix, label + '.irc_attempt1' + suffix)
            kwargs = reader_gauss.irc_retry_kwargs(kwargs)
            mol.positions = initial_geometry.copy()
            mol.calc = Gaussian(**kwargs)
    if not success:
        logging.warning('%s: no accepted IRC point after retry; endpoint optimization skipped.', label)
        db.write(mol, name=label, data={{'status': 'error', 'irc_completed': False,
                                      'irc_accepted_points': 0}})

    # Preserve the actual accepted IRC points, including a usable partial path.
    # Save these before endpoint optimization can fail or overwrite its own log.
    reader_gauss.write_irc_trajectory(logfile, label + '.xyz')

else:
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
                                                    mult=kwargs['mult'],
                                                    keep_input={bimol})
        freq = reader_gauss.read_freq(logfile, {atom})
        zpe = reader_gauss.read_zpe(logfile)
        db.write(mol_prod, name=label, data={{'energy': e,
                                         'frequencies': np.asarray(freq),
                                         'zpe': zpe, 'status': 'normal'}})
    except RuntimeError: 
        for i in range(3):
            try:
                iowait(logfile, 'gauss')
                if {bimol}:
                    _, mol_prod.positions = reader_gauss.read_bimol_retry_geom_energy(logfile, mol_prod)
                else:
                    _, mol_prod.positions = reader_gauss.read_lowest_geom_energy(logfile, mol_prod)
                prod_kwargs = reader_gauss.correct_kwargs(logfile, prod_kwargs, bimol={bimol})
                mol_prod.calc = Gaussian(**prod_kwargs)
                e = mol_prod.get_potential_energy()  # use the Gaussian optimizer
                iowait(logfile, 'gauss')
                mol_prod.positions = reader_gauss.read_geom(logfile, mol_prod, keep_input={bimol})
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
