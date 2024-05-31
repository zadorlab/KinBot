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

def get_interfragments_param(atoms, instance):
    species = StationaryPoint.from_ase_atoms(atoms)
    species.characterize()
    """
    Function that selects the neighbors of the closest atom in each fragment,
    and returns the values of the interfragments angles and dihedrals.
    The level is only used as small orientation correction for planar fragments.
    """
    key_dihedral = None
    save_name = copy.copy(species.name)
    closest_to_RA = [[],[]] #each list contains the 2 (or less) closest atoms to the reactive atom in their respective fragment
    # cutting up structure into two fragments
    species.name = save_name
    species.bond[instance[0], instance[1]] = 0
    species.bond[instance[1], instance[0]] = 0
    species.bonds[0][instance[0], instance[1]] = 0
    species.bonds[0][instance[1], instance[0]] = 0
    fragments, maps = species.start_multi_molecular()
    for frag_number, ra in enumerate(instance):
        for neighbor_index, bond in zip(maps[frag_number], fragments[frag_number].bond[np.where(maps[frag_number] == ra)[0][0]]):
            if bond != 0:
                closest_to_RA[frag_number].append(neighbor_index)
                if len(closest_to_RA[frag_number]) == 2:
                    break
        #If no neighbor atoms were found, the fragment has to be mono-atomic
        if len(closest_to_RA[frag_number]) == 1 and fragments[frag_number].natom > 2:
            for neighbor_index, bond in zip(maps[frag_number], fragments[frag_number].bond[np.where(maps[frag_number] == closest_to_RA[frag_number][0])[0][0]]):
                if bond != 0 and neighbor_index != instance[frag_number]:
                    closest_to_RA[frag_number].append(neighbor_index)
                    if len(closest_to_RA[frag_number]) == 2:
                        break
    changes = []
    #Calculate the angle values:
    if len(closest_to_RA[0]) != 0:
        point_A = species.geom[closest_to_RA[0][0]]
        point_B = species.geom[instance[0]]
        point_C = species.geom[instance[1]]
        changes.append([closest_to_RA[0][0],
                                instance[0],
                                instance[1],
                                (np.degrees(geometry.calc_angle(point_A, point_B, point_C)))])
    if len(closest_to_RA[1]) != 0:
        point_A = species.geom[closest_to_RA[1][0]]
        point_B = species.geom[instance[1]]
        point_C = species.geom[instance[0]]
        changes.append([closest_to_RA[1][0],
                                instance[1],
                                instance[0],
                                (np.degrees(geometry.calc_angle(point_A, point_B, point_C)))])
    
    #Calculate the dihedral values:
    if len(closest_to_RA[0]) != 0 and len(closest_to_RA[1]) != 0:
        point_A = species.geom[closest_to_RA[0][0]]
        point_B = species.geom[instance[0]]
        point_C = species.geom[instance[1]]
        point_D = species.geom[closest_to_RA[1][0]]
        changes.append([closest_to_RA[0][0],
                        instance[0],
                        instance[1],
                        closest_to_RA[1][0],
                        geometry.calc_dihedral(point_A, point_B, point_C, point_D)[0]])
        key_dihedral = changes[-1][:4]
    if len(closest_to_RA[0]) == 2:
        point_A = species.geom[closest_to_RA[0][1]]
        point_B = species.geom[closest_to_RA[0][0]]
        point_C = species.geom[instance[0]]
        point_D = species.geom[instance[1]]
        changes.append([closest_to_RA[0][1],
                        closest_to_RA[0][0],
                        instance[0],
                        instance[1],
                        geometry.calc_dihedral(point_A, point_B, point_C, point_D)[0]])
    if len(closest_to_RA[1]) == 2:
        point_A = species.geom[closest_to_RA[1][1]]
        point_B = species.geom[closest_to_RA[1][0]]
        point_C = species.geom[instance[1]]
        point_D = species.geom[instance[0]]
        changes.append([closest_to_RA[1][1],
                        closest_to_RA[1][0],
                        instance[1],
                        instance[0],
                        geometry.calc_dihedral(point_A, point_B, point_C, point_D)[0]])
    
    return changes, key_dihedral

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

initial_inter_frag, key_dihedral = get_interfragments_param(mol, instance={instance})

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
