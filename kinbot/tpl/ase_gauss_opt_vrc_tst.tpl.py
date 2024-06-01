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

def same_orientation(initial, final):
    is_true = True
    for initial_parameter, final_parameter in zip(initial, final):
        if len(initial_parameter) == 3: #Parameter is a bond
            if final_parameter[-1]/initial_parameter[-1] < 1 - {bond_deviation} or final_parameter[-1]/initial_parameter[-1] > 1 + {bond_deviation}:
                is_true = False
                break
        else:
            if final_parameter[-1] < initial_parameter[-1] - {angle_deviation} or final_parameter[-1] > initial_parameter[-1] + {angle_deviation}:
                is_true = False
                break
    return is_true
    
db = connect('{working_dir}/kinbot.db')
label = '{label}'
logfile = '{label}.log'

mol = Atoms(symbols={atom}, positions={geom})

# convert to z-matrix to get the important atoms and initial interfragmental parameters
spec = StationaryPoint.from_ase_atoms(mol)
zmat_atom, zmat_ref, zmat, zmatorder = make_zmat_from_cart(spec, {frozen_coo}, mol.positions, 2)
A = 
B = {frozen_coo}[0]
C = {frozen_coo}[1]

initial_inter_frag, key_dihedral = get_interfragments_param(mol, frozen_coo={frozen_coo})

kwargs = {kwargs}
Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc

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
        mol.positions = geom  # final mol is the updated one
        new_inter_frag, _ = get_interfragments_param(mol, instance={instance})
        if same_orientation(initial_inter_frag, new_inter_frag):
            db.write(mol, name=label, data={{'energy': energies(gii), 'status': 'normal'}})
            break


time.sleep(1)  # Avoid db errors
with open(logfile, 'a') as f:
    f.write('done\n')


# the original species
origmol = Atoms(symbols={atom}, positions={frozen_geom})
origspec = StationaryPoint.from_ase_atoms(origmol)
origspec.characterize()
# the frozen species (not yet at correct geometry)
frozenspec = StationaryPoint.from_ase_atoms(mol)
frozenspec.characterize()
# zmat_atom: element symbol in order
# zmat_ref: referencing atom for D, A, Dh
# zmat: values
# zmatorder: atom number in original species
zmat_atom, zmat_ref, zmat, zmatorder = make_zmat_from_cart(frozenspec, key_dihedral, mol.positions, 2)
# need to adjust the parameters in the zmat to the rigid values
# but not the values that are across the fragments
for zii, zaa in enumerate(zmat_atom):
    if zii == 0:  # nothing to change here
        continue
    if zii == 1:  # fixed interfrag distance
        continue
    A = zmatorder[zii]  # current atom's original index
    B, C, D = zmat_ref[zii]  # 0 indexed?
    zmat[zii][0] = origspec.dist[A][B]  # distance
    if zii == 2:  # fixed interfrag angle
        continue
    zmat[zii][1] = np.degree(geometry.calc_angle(origspec[A], origspec[B], origspec[C]))  # angle
    if zii == 3:  # fixed interfag dihed
        continue
    zmat[zii][2] = np.degree(geometry.calc_dihedral(origspec[A], origspec[B], origspec[C], origspec[D]))  # dihedral
# convert back to cartesian
cart = make_cart_from_zmat(zmat, zmat_atom, zmat_ref, frozenspec.natom, frozenspec.atom, zmatorder)
    
# here creating the run for the frozen version of this geometry based on the fully optimized one
label = '{frozen_label}'
logfile = '{frozen_label}.log'
kwargs.pop('opt', None)
mol = Atoms(symbols={atom}, positions=cart)
try:
    e = mol.get_potential_energy()  # single point energy
    db.write(mol, name=label, data={{'energy': e, 'status': 'normal'}})
except:
    db.write(mol, name=label, data={{'energy': 10., 'status': 'normal'}})

time.sleep(1)  # Avoid db errors
with open(logfile, 'a') as f:
    f.write('done\n')
