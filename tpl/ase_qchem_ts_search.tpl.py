from ase import Atoms
from ase.calculators.qchem import QChem
from ase.db import connect
from kinbot import reader_qchem

db = connect('{working_dir}/kinbot.db')

scan = {scan}
bimol = {bimol}
mol = Atoms(symbols={atom}, positions={geom})

kwargs = {kwargs}
QChem.command = '{qc_command} -nt {ppn} PREFIX.in PREFIX.out PREFIX.sv'
calc = QChem(**kwargs)
mol.set_calculator(calc)

success = True
for tr in range({ntrial}):
    try:
        e = mol.get_potential_energy()  # use the QChem optimizer
        mol.positions = reader_qchem.read_geom('{label}.out', mol)
        db.write(mol, name='{label}', data={{'energy': e,
                                             'status': 'normal'}})
        break
    except RuntimeError:
        success = False

if not success:
    if not bimol:
        try:
            mol.positions = reader_qchem.read_geom('{label}.out', mol)
            kwargs['jobtype'] = 'sp'  # this is when we give up optimization!!
            calc = QChem(**kwargs)
            e = mol.get_potential_energy()
            db.write(mol, name='{label}', data={{'energy': e,
                                                 'status': 'normal'}})
        except RuntimeError:
            db.write(mol, name='{label}', data={{'status': 'error'}})
    else:
        try:  # TODO I don't understand this
            mol.positions = reader_qchem.read_geom('{label}.out', mol)
            db.write(mol, name='{label}', data={{'energy': e,
                                                 'status': 'normal'}})
        except RuntimeError:
            db.write(mol, name='{label}', data={{'status': 'error'}})

with open(f'{label}.out', 'a') as f:
    f.write('done\n')
