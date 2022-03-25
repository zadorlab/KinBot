from ase import Atoms
from ase.calculators.qchem import QChem
from ase.db import connect
from kinbot import reader_qchem
from kinbot.utils import iowait

db = connect('{working_dir}/kinbot.db')
label = '{label}'
logfile = '{label}.out'

mol = Atoms(symbols={atom}, positions={geom})

kwargs = {kwargs}
calc = QChem(**kwargs)
mol.calc = calc
calc.command = '{qc_command} -nt {ppn} -save PREFIX.inp PREFIX.out PREFIX.sv'

# Compute the Hessian via a Frequencies calculation
e = mol.get_potential_energy()  # Compute Freq

# Carry out the IRC calculation with the Hessian from the last calculation
kwargs['jobtype'] = 'rpath'
calc = QChem(**kwargs)
calc.command = '{qc_command} -nt {ppn} PREFIX.inp PREFIX.out PREFIX.sv'
mol.set_calculator(calc)
success = True
try:
    e = mol.get_potential_energy()  # Perform IRC
    iowait(logfile, 'qchem')
    mol.positions = reader_qchem.read_geom(logfile, mol, irc=True)
    db.write(mol, name=label, data={{'energy': e, 'status': 'normal'}})
except RuntimeError:
    mol.positions = reader_qchem.read_geom(logfile, mol, irc=True)
    if mol.positions is not None:
        db.write(mol, name=label, data={{'status': 'normal'}})  # although there is an error, continue from the final
    else:
        db.write(mol, name=label, data={{'status': 'error'}})
        success = False

with open(logfile, 'a') as f:
    f.write('done\n')

if success:
    label = '{label}_prod'
    logfile = '{label}_prod.out'
    # start the product optimization
    prod_kwargs = {prod_kwargs}
    calc_prod = QChem(**prod_kwargs)
    calc_prod.command = '{qc_command} -nt {ppn} PREFIX.inp PREFIX.out'
    mol_prod = Atoms(symbols={atom}, positions=mol.positions)
    mol_prod.calc = calc_prod
    try:
        e = mol_prod.get_potential_energy()  # use the QChem optimizer
        iowait(logfile, 'qchem')
        mol_prod.positions = reader_qchem.read_geom(logfile, mol_prod)
        db.write(mol_prod, name=label, data={{'energy': e, 'status': 'normal'}})  # TODO mol_prod?
    except RuntimeError: 
        mol_prod.positions = reader_qchem.read_geom(logfile, mol_prod)
        if mol_prod.positions is not None:
            db.write(mol_prod, name=label, data={{'status': 'normal'}})
        else:
            db.write(mol_prod, name=label, data={{'status': 'error'}})

with open(logfile, 'a') as f:
    f.write('done\n')
