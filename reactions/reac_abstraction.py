import sys
import numpy as np
import copy
from kinbot import geometry
from kinbot.reac_General import GeneralReac
from kinbot import geometry
from kinbot import modify_geom


class Abstraction(GeneralReac):
    max_step = 2
    scan = 0
    skip = 0
    family_name = 'abstraction'

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        
        if step == 0:
            self.fix_bonds(fix)
            self.fix_bond_single(0, 2, fix)
            self.fix_angles(fix)
        # in step 1 the QST3 method is used

        self.clean_constraints(change, fix)
        return step, fix, change, release

def abstraction_align(startgeom, instance, atom, fragnatom):
    """
    Special function to align fragments. It gives a guess for the
    reactant and product geometries that are used in QST2.
    startgeom: starting geometry to be aligned, contains the two fragments, the first abstracts.
    instance: atoms in the instance
    atom: list of atoms
    fragnatom: the number of atoms in the first fragment as defined in the original stationary point(!)
    """

    # align O----H, O is at origin, H is on +Z axis
    g0 = copy.deepcopy(geometry.translate_and_rotate(startgeom, instance[0], instance[1]))
    # align H-C, H is at origin, C is on +Z axis at an angle
    g1 = copy.deepcopy(geometry.translate_and_rotate(startgeom, instance[1], instance[2]))
    for i in range(len(g1)):
        g1[i] = modify_geom.perform_rotation(g1[i], g1[instance[1]], np.array([0, 1, 0]), 10. / 180. * np.pi)
        #g1[i] = geometry.rotate_atom(g1[i], [0, 1, 0], 20. / 180. * np.pi)
    # C-H distance
    val0 = np.linalg.norm(startgeom[instance[1]] - startgeom[instance[2]])
    val = distance(atom, 0, 2) - val0
    g1[:, 2] += val  # shift the whole thing in z direction
    geom_reac = np.concatenate((g0[:fragnatom], g1[fragnatom:]))  # does not work for reverse
    geom_prod = copy.deepcopy(geom_reac)
    geom_prod[instance[1]] = [0, 0, 1.]
    geom_ts = copy.deepcopy(geom_reac)
    geom_ts[instance[1]] = [0, 0, 1.5]
    return geom_reac, geom_prod, geom_ts

def distance(atom, aa, cc):
    """
    Given atom types a and c, return distance.
    Thes are in an a--b--c abstraction, a is the abstractor
    Not used, can be deleted.
    """
    
    a = atom[aa]
    c = atom[cc]

    if a == 'O':
        if c == 'O': return 2.5
        elif c == 'C': return 3.0
        else: return 3.
    elif a == 'Cl':
        if c == 'O': return 3.0
        elif c == 'C': return 2.8
        else: return 3.
    elif a == 'H':
        if c == 'O': return 2.3
        elif c == 'C': return 2.3
        else: return 3.
    else: 
        return 3.
