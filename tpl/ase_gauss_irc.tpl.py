from ase import Atoms
# from ase.calculators.gaussian import Gaussian
from ase.db import connect

from kinbot.ase_modules.calculators.gaussian import Gaussian  # New
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
        mol_prod.positions = reader_gauss.read_geom(logfile, mol_prod)
        db.write(mol_prod, name=label, data={{'energy': e, 'status': 'normal'}})
    except RuntimeError: 
        iowait(logfile, 'gauss')
        mol_prod.positions = reader_gauss.read_geom(logfile, mol_prod)
        if mol_prod.positions is not None:
            db.write(mol_prod, name=label, data={{'status': 'normal'}}) 
        else:
            db.write(mol_prod, name=label, data={{'status': 'error'}})

    with open(logfile, 'a') as f:
        f.write('done\n')
