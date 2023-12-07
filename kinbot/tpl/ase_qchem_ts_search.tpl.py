from ase import Atoms
from ase.calculators.qchem import QChem
from ase.db import connect

from kinbot import reader_qchem
from kinbot.utils import iowait

db = connect('{working_dir}/kinbot.db')
label = '{label}'
logfile = '{label}.out'

scan = {scan}
bimol = {bimol}
mol = Atoms(symbols={atom}, positions={geom})

kwargs = {kwargs}
calc = QChem(**kwargs)
calc.command = '{qc_command} -nt {ppn} -save PREFIX.inp PREFIX.out PREFIX.sv'
mol.set_calculator(calc)

success = True
try:
    e = mol.get_potential_energy()  # use the QChem optimizer
    iowait(logfile, 'qchem')
    mol.positions = reader_qchem.read_geom(logfile, mol)
    db.write(mol, name=label, data={{'energy': e, 'status': 'normal'}})
except RuntimeError:
    success = False

if not success:
    if not bimol:
        try:
            mol.positions = reader_qchem.read_geom(logfile, mol)
            kwargs['jobtype'] = 'sp'  # this is when we give up optimization!!
            calc = QChem(**kwargs)
            e = mol.get_potential_energy()
            iowait(logfile, 'qchem')
            mol.positions = reader_qchem.read_geom(logfile, mol)
            db.write(mol, name=label, data={{'energy': e, 'status': 'normal'}})
        except RuntimeError:
            db.write(mol, name=label, data={{'status': 'error'}})
    else:
        try:
            mol.positions = reader_qchem.read_geom(logfile, mol)
            db.write(mol, name='{label}', data={{'energy': e, 'status': 'normal'}})
        except RuntimeError:
            db.write(mol, name='{label}', data={{'status': 'error'}})

with open(logfile, 'a') as f:
    f.write('done\n')
