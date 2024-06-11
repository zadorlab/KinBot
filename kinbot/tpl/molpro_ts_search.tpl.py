from ase import Atoms
from ase.db import connect

from kinbot.molpro_calculator import Molpro_calc
from kinbot import reader_molpro
from kinbot.utils import iowait

db = connect('{working_dir}/kinbot.db')
label = '{label}'
xmlfile = '{label}.xml'

scan = {scan}
bimol = {bimol}
mol = Atoms(symbols={atom}, positions={geom})

kwargs = {kwargs}
Molpro_calc.command = '{qc_command} PREFIX.inp'
calc = Molpro_calc(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy()  # use the Molpro optimizer
except RuntimeError: 
    e = 0.
 
iowait(xmlfile, 'gauss')
mol.positions = reader_molpro.read_geom(xmlfile, mol)
if all([ci == 0 for mp in mol.positions for ci in mp]):
    mol.positions = {geom}  # reset to the original geometry
db.write(mol, name=label, data={{'energy': e, 'status': 'normal'}})

#for tr in range(ntrial):  # DELETED CURLY BRACKET
#    try:
#        success = True
#        e = mol.get_potential_energy() # use the Gaussian optimizer (task optimize)
#        iowait(logfile, 'gauss')
#        mol.positions = reader_gauss.read_geom(logfile, mol)
#        db.write(mol, name=label, data={{'energy': e,'status': 'normal'}})
#        break
#    except RuntimeError: 
#        success = False
#        
#if not success:
#    if not bimol:
#        try:
#            mol.positions = reader_gauss.read_geom(logfile, mol)
#            del kwargs['opt']  # this is when we give up optimization!!
#            calc = Gaussian(**kwargs)
#            e = mol.get_potential_energy() 
#            iowait(logfile, 'gauss')
#            mol.positions = reader_gauss.read_geom(logfile, mol)
#            db.write(mol, name=label, data={{'energy': e,'status': 'normal'}})
#        except: 
#            db.write(mol, name = label, data = {{'status': 'error'}})
#    else:
#        try:
#            mol.positions = reader_gauss.read_geom(logfile, mol)
#            db.write(mol, name=label, data={{'energy': e,'status': 'normal'}})
#        except: 
#            db.write(mol, name = label, data = {{'status': 'error'}})

with open(xmlfile,'a') as f:
    f.write('done\n')
