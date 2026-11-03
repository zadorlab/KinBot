from ase import Atoms
from ase.db import connect
import numpy as np

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait

db = connect('{working_dir}/kinbot.db')
label = '{label}'
logfile = '{label}.log'

scan = {scan}
bimol = {bimol}
mol = Atoms(symbols={atom}, positions={geom})

kwargs = {kwargs}
Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

try:
    e = mol.get_potential_energy()  # use the Gaussian optimizer
except RuntimeError:
    if {strict_scan}:
        raise
    # An unsuccessful path guess is not a normally completed energy point.
    db.write(mol, name=label, data={{'status': 'error'}})
    raise
 
iowait(logfile, 'gauss')
mol.positions = reader_gauss.read_geom(logfile, mol)
if all([ci == 0 for mp in mol.positions for ci in mp]):
    raise RuntimeError('Gaussian path calculation returned no geometry.')
if not np.isfinite(e) or not np.isfinite(mol.positions).all():
    raise RuntimeError('Gaussian path calculation returned nonfinite data.')
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

with open(logfile,'a') as f:
    f.write('done\n')
