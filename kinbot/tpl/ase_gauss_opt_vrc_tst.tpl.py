import numpy as np
from ase import Atoms
from ase.db import connect
import time
import os
import copy

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait
from kinbot.stationary_pt import StationaryPoint
from kinbot import geometry
from kinbot import zmatrix
import rmsd

db = connect('{working_dir}/kinbot.db')
label = '{label}'
logfile = '{label}.log'

mol = Atoms(symbols={atom}, positions={init_geom})

kwargs = {kwargs}
Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

# RELAXED

try:
    e = mol.get_potential_energy()  # use the Gaussian optimizer
except:
    pass  # it doesn't matter if it cannot converge

iowait(logfile, 'gauss')
# both geoms and energies are in reversed order
geoms = reader_gauss.read_all_geoms(logfile, mol)
energies = reader_gauss.read_all_energies(logfile)
if len(geoms) == 0:  # no optimization worked, assuming that we have at least one energy
    db.write(mol, name=label, data={{'energy': energies[-1], 'status': 'normal'}})
else:  # select the last geometry that was within the range of allowed change
    for gii, geom in enumerate(geoms): 
        if rmsd.kabsch_rmsd({init_geom}, geom, translate=True) < {scan_deviation}:
            mol.positions = geom  # update geometry
            db.write(mol, name=label, data={{'energy': energies(gii), 'status': 'normal'}})
            break

time.sleep(1)  
with open(logfile, 'a') as f:
    f.write('done\n')

# FROZEN

scan_coo = {scan_coo}

# create fragment geometries from the relaxed scan and save scan atoms' position and index
coo_A = [mol.positions[i] for i in {frag_maps}[0]]
try:
    scan_atom_A = {frag_maps}[0].index(scan_coo[0])
except ValueError:
    scan_atom_A = {frag_maps}[0].index(scan_coo[1])
scan_pos_A = coo_A[scan_atom_A]

coo_B = [mol.positions[i] for i in {frag_maps}[1]]
try:
    scan_atom_B = {frag_maps}[1].index(scan_coo[0])
except ValueError:
    scan_atom_B = {frag_maps}[1].index(scan_coo[1])
scan_pos_B = coo_B[scan_atom_B]

# Move fragments to own centroids (translation only)
cent_A = rmsd.centroid(coo_A)
cent_B = rmsd.centroid(coo_B)
coo_A -= cent_A
coo_B -= cent_B

# Move frozen fragments to own centroids (translation only)
cent_A_fr = rmsd.centroid({froz_A_geom})
cent_B_fr = rmsd.centroid({froz_B_geom})
coo_A_fr -= cent_A_fr
coo_B_fr -= cent_B_fr

# rotate each frozen fragment for maximum overlap with own relaxed version
coo_A_fr = rmsd.kabsch_rotate(coo_A_fr, coo_A)
coo_B_fr = rmsd.kabsch_rotate(coo_B_fr, coo_B)

# move fragments back to scan position
coo_A_fr -= coo_A_fr[scan_atom_A] - scan_pos_A
coo_B_fr -= coo_B_fr[scan_atom_B] - scan_pos_B

# reassemble full geometry from frozen fragments
frozen_geom = np.array(len({atom}))
for fi, fa in enumerate(frag_A_atoms):
    frozen_geom[fa] = coo_A_fr[fi]
for fi, fb in enumerate(frag_B_atoms):
    frozen_geom[fb] = coo_B_fr[fi]
    
# here creating the run for the frozen version of this geometry based on the fully optimized one
label = '{label}_fr'
logfile = '{label}_fr.log'
kwargs.pop('opt', None)
mol = Atoms(symbols={atom}, positions=frozen_geom)
try:
    e = mol.get_potential_energy()  # single point energy
    db.write(mol, name=label, data={{'energy': e, 'status': 'normal'}})
except:
    db.write(mol, name=label, data={{'energy': 10., 'status': 'normal'}})

time.sleep(1)  
with open(logfile, 'a') as f:
    f.write('done\n')
