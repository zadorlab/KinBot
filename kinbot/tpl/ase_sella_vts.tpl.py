import numpy as np
from ase import Atoms
from ase.db import connect
import time
import os
import copy
from sella import Sella, Constraints
from kinbot.ase_modules.calculators.{code} import {Code}

from kinbot.stationary_pt import StationaryPoint
from kinbot import geometry
from kinbot import zmatrix
import rmsd

db = connect('{working_dir}/kinbot.db')
label = '{label}'
logfile = '{label}.log'

mol = Atoms(symbols={atom}, positions={init_geom})
kwargs = {kwargs}
mol.calc = {Code}(**kwargs)
# RELAXED

scan_coo = {scan_coo}

if os.path.isfile('{label}_sella.log'):
    os.remove('{label}_sella.log')
sella_kwargs = {sella_kwargs}
cons = Constraints(mol)
cons.fix_bond((scan_coo[0], scan_coo[1]))
opt = Sella(mol, 
            order=0, 
            constraints=cons,
            internal=True,
            trajectory='{label}.traj', 
            logfile='{label}_sella.log',
            **sella_kwargs)

mol.calc.label = '{label}'
bonds = {bonds}
distances = np.array([np.linalg.norm(mol.positions[bond[0]] - mol.positions[bond[1]]) for bond in bonds])

last = True  # take the last geometry, otherwise the one before that
for i in opt.irun(fmax=1e-4, steps=100):
    if rmsd.kabsch_rmsd(np.array({init_geom}), mol.positions, translate=True) > {scan_deviation}:
        last = False
        print('rmsd is too large, optimization is stopped')
        break
    curr_distances = np.array([np.linalg.norm(mol.positions[bond[0]] - mol.positions[bond[1]]) for bond in bonds])
    ratio = curr_distances / distances
    if any([True if ri > 1.1 else False for ri in ratio]):
        last = False
        print('a bond is more than 10% strethed, optimization is stopped')
        break
    mol_prev = copy.deepcopy(mol)
    e = mol.get_potential_energy() 

if not last:
    mol = mol_prev

db.write(mol, name='{label}',
         data={{'energy': e, 'status': 'normal'}})

with open('{label}.log', 'a') as f:
    f.write('done\n')

# FROZEN

# create fragment geometries from the relaxed scan and save scan atoms' position and index
coo_A = np.array([mol.positions[i] for i in {frag_maps}[0]])
try:
    scan_atom_A = {frag_maps}[0].index(scan_coo[0])
except ValueError:
    scan_atom_A = {frag_maps}[0].index(scan_coo[1])
scan_pos_A = copy.deepcopy(coo_A[scan_atom_A])

coo_B = np.array([mol.positions[i] for i in {frag_maps}[1]])
try:
    scan_atom_B = {frag_maps}[1].index(scan_coo[0])
except ValueError:
    scan_atom_B = {frag_maps}[1].index(scan_coo[1])
scan_pos_B = copy.deepcopy(coo_B[scan_atom_B])

# Move fragments to own centroids (translation only)
cent_A = rmsd.centroid(coo_A)
cent_B = rmsd.centroid(coo_B)
coo_A -= cent_A
coo_B -= cent_B

# Move frozen fragments to own centroids (translation only)
coo_A_fr = np.array({froz_A_geom})
coo_B_fr = np.array({froz_B_geom})
cent_A_fr = rmsd.centroid(coo_A_fr)
cent_B_fr = rmsd.centroid(coo_B_fr)
coo_A_fr -= cent_A_fr
coo_B_fr -= cent_B_fr

# rotate each frozen fragment for maximum overlap with own relaxed version
coo_A_fr = rmsd.kabsch_rotate(coo_A_fr, coo_A)
coo_B_fr = rmsd.kabsch_rotate(coo_B_fr, coo_B)

# move fragments back to scan position
coo_A_fr -= coo_A_fr[scan_atom_A] - scan_pos_A
coo_B_fr -= coo_B_fr[scan_atom_B] - scan_pos_B

# reassemble full geometry from frozen fragments
frozen_geom = np.zeros((len({atom}), 3))
for fi, fa in enumerate({frag_maps}[0]):
    frozen_geom[fa] = coo_A_fr[fi]
for fi, fb in enumerate({frag_maps}[1]):
    frozen_geom[fb] = coo_B_fr[fi]
    
# here creating the run for the frozen version of this geometry based on the fully optimized one
label = '{label}_fr'
logfile = '{label}_fr.log'
kwargs.pop('opt', None)
kwargs.pop('addsec', None)
kwargs['label'] = label
calc = Gaussian(**kwargs)
mol = Atoms(symbols={atom}, positions=frozen_geom)
mol.calc = calc
try:
    e = mol.get_potential_energy()  # single point energy
    db.write(mol, name=label, data={{'energy': e, 'status': 'normal'}})
except:
    db.write(mol, name=label, data={{'energy': 10., 'status': 'normal'}})

time.sleep(1)  
with open(logfile, 'a') as f:
    f.write('done\n')
