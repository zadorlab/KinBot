import numpy as np
from ase import Atoms
from ase.db import connect
import time
import os

from kinbot.ase_modules.calculators.gaussian import Gaussian
from kinbot import reader_gauss
from kinbot.utils import iowait
from kinbot.stationary_pt import StationaryPoint
from kinbot import geometry
import copy

def correct_kwargs(kwargs, iteration):
    """Funtion that refine the optimization parameters depending on the number of trials"""
    if 'opt' not in kwargs or iteration == -1:
        return kwargs
    match iteration:
        case 0:
            kwargs["iop"] = "1/19=4"
            kwargs['opt'] = kwargs['opt'].replace("MaxCycles=15", "MaxCycles=40")
        case 1:
            kwargs['opt'] = kwargs['opt'].replace('CalcFC', 'ReCalcFC=10')
        case 2: #Use different optimization procedure
            kwargs['opt'] = kwargs['opt'].replace('ReCalcFC=10', 'CalcAll')
            kwargs['opt'] = kwargs['opt'].replace("MaxCycles=40", "MaxCycles=30")
    return kwargs

def get_interfragments_param(atoms, instance):
    species = StationaryPoint.from_ase_atoms(atoms)
    species.characterize()
    """
    Function that selects the neighbourgs of the closest atom in each fragment,
    and returns the values of the interfragments angles and dihedrals.
    The level is only used as small orientation correction for planar fragments.
    """
    save_name = copy.copy(species.name)
    closest_to_RA = [[],[]] #each list contains the 2 (or less) closest atoms to the reactive atom in their respective fragment
    fragments, maps = species.start_multi_molecular()
    if len(fragments) == 1:
        species.name = save_name
        species.bond[instance[0], instance[1]] = 0
        species.bond[instance[1], instance[0]] = 0
        species.bonds[0][instance[0], instance[1]] = 0
        species.bonds[0][instance[1], instance[0]] = 0
        fragments, maps = species.start_multi_molecular()
    for frag_number, ra in enumerate(instance):
        for neighbourg_index, bond in zip(maps[frag_number], fragments[frag_number].bond[np.where(maps[frag_number] == ra)[0][0]]):
            if bond != 0:
                closest_to_RA[frag_number].append(neighbourg_index)
                if len(closest_to_RA[frag_number]) == 2:
                    break
        #If no neighbourgs atoms were found, the fragment has to be mono-atomic
        if len(closest_to_RA[frag_number]) == 1 and fragments[frag_number].natom > 2:
            for neighbourg_index, bond in zip(maps[frag_number], fragments[frag_number].bond[np.where(maps[frag_number] == closest_to_RA[frag_number][0])[0][0]]):
                if bond != 0 and neighbourg_index != instance[frag_number]:
                    closest_to_RA[frag_number].append(neighbourg_index)
                    if len(closest_to_RA[frag_number]) == 2:
                        break
    changes = []
    #Calculate the angles values:
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
    
    #Calculate the dihedrals values:
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
    
    return changes

def same_orientation(initial, final):
    is_true = True
    for initial_parameter, final_parameter in zip(initial, final):
        if len(initial_parameter) == 3: #Parameter is a bond
            if final_parameter[-1]/initial_parameter[-1] < 0.95 or final_parameter[-1]/initial_parameter[-1] > 1.05:
                is_true = False
                break
        else:
            if final_parameter[-1] < initial_parameter[-1] - 10 or final_parameter[-1] > initial_parameter[-1] + 10:
                is_true = False
                break
    return is_true
    
db = connect('{working_dir}/kinbot.db')
label = '{label}'
logfile = '{label}.log'

mol = Atoms(symbols={atom}, positions={geom})

initial_inter_frag = get_interfragments_param(mol, instance={instance})

kwargs = {kwargs}
Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)
mol.calc = calc
constrain_orientation = False

try:
    e = mol.get_potential_energy()  # use the Gaussian optimizer
    iowait(logfile, 'gauss')
    mol.positions = reader_gauss.read_geom(logfile, mol)
    new_inter_frag = get_interfragments_param(mol, instance={instance})
    if same_orientation(initial_inter_frag, new_inter_frag):
        freq = reader_gauss.read_freq(logfile, {atom})
        zpe = reader_gauss.read_zpe(logfile)
        db.write(mol, name=label, data={{'energy': e, 'frequencies': np.asarray(freq),
                                        'zpe': zpe, 'status': 'normal'}})
    else:
        constrain_orientation = True
        #Change kwargs
        if 'frozen' in label:
            kwargs.pop('opt', None)
        else:
            kwargs['addsec'] = ''
            for param in initial_inter_frag:
                for indx in param[:-1]:
                    kwargs['addsec'] += f'{{indx+1}} '
                kwargs['addsec'] += 'F\n'
        #reinitialize molecule
        mol = Atoms(symbols={atom}, positions={geom})
        #reinitialize kwargs
        calc = Gaussian(**kwargs)
        mol.calc = calc
        #recalculate
        e = mol.get_potential_energy()  # use the Gaussian optimizer
        mol.positions = reader_gauss.read_geom(logfile, mol)
        new_inter_frag = get_interfragments_param(mol, instance={instance})
        if same_orientation(initial_inter_frag, new_inter_frag):
            freq = reader_gauss.read_freq(logfile, {atom})
            zpe = reader_gauss.read_zpe(logfile)
            db.write(mol, name=label, data={{'energy': e, 'frequencies': np.asarray(freq),
                                            'zpe': zpe, 'status': 'normal'}})
        else:
            raise Exception("Has converged to a different reaction coordinate.")
except:
    i = 0
    while i < 3:
        try:
            iowait(logfile, 'gauss')
            e, mol.positions = reader_gauss.read_lowest_geom_energy(logfile, mol)
            new_inter_frag = get_interfragments_param(mol, instance={instance})
            if (not same_orientation(initial_inter_frag, new_inter_frag) and not constrain_orientation) or i == 3:
                constrain_orientation = True
                i = -1
                #reinitialize molecule
                mol = Atoms(symbols={atom}, positions={geom})
                #reinitialize kwargs
                kwargs = {kwargs}
                if 'frozen' in label:
                    kwargs.pop('opt', None)
                else:
                    kwargs['addsec'] = ''
                    for param in initial_inter_frag:
                        for indx in param[:-1]:
                            kwargs['addsec'] += f'{{indx+1}} '
                        kwargs['addsec'] += 'F\n'
            kwargs = correct_kwargs(kwargs, i)
            mol.calc = Gaussian(**kwargs)
            #recalculate
            e = mol.get_potential_energy()  # use the Gaussian optimizer
            mol.positions = reader_gauss.read_geom(logfile, mol)
            new_inter_frag = get_interfragments_param(mol, instance={instance})
            if same_orientation(initial_inter_frag, new_inter_frag):
                freq = reader_gauss.read_freq(logfile, {atom})
                zpe = reader_gauss.read_zpe(logfile)
                db.write(mol, name=label, data={{'energy': e, 'frequencies': np.asarray(freq),
                                                'zpe': zpe, 'status': 'normal'}})
                break
            elif constrain_orientation:
                db.write(mol, name=label, data={{'status': 'error'}})
                break
        except:
            iowait(logfile, 'gauss')
            #Save in db the lowest energy geometry if forces are converged
            if reader_gauss.read_convergence(logfile) != 0:
                try:
                    e, mol.positions = reader_gauss.read_converged_geom_energy(logfile, mol)
                    new_inter_frag = get_interfragments_param(mol, instance={instance})
                    if same_orientation(initial_inter_frag, new_inter_frag):
                        freq = reader_gauss.read_freq(logfile, {atom})
                        zpe = reader_gauss.read_zpe(logfile)
                        db.write(mol, name=label, data={{'energy': e, 'frequencies': np.asarray(freq),
                                                        'zpe': zpe, 'status': 'normal'}})
                        break
                    elif constrain_orientation:
                        db.write(mol, name=label, data={{'status': 'error'}})
                        break
                except:
                    pass
            else:
                pass
            if i == 2 and os.path.getsize(logfile) != 0 and constrain_orientation:
                e, mol.positions = reader_gauss.read_lowest_geom_energy(logfile, mol)
                freq = reader_gauss.read_freq(logfile, {atom})
                zpe = reader_gauss.read_zpe(logfile)
                db.write(mol, name=label, data={{'status': 'error'}})
                break
            elif i == 2 and constrain_orientation:
                db.write(mol, name=label, data={{'status': 'error'}})
                break
        else:
            break
    i += 1


time.sleep(1) #Avoid db errors
with open(logfile, 'a') as f:
    f.write('done\n')
