import numpy as np
import copy
from kinbot import geometry
from kinbot.reac_General import GeneralReac
from kinbot import geometry


class Abstraction(GeneralReac):
    max_step = 14
    scan = 0
    skip = 0
    family_name = 'abstraction'

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        if step == 0:
            self.set_bond(1, 2, 1.2, change)
        if step < self.max_step:
            self.fix_bonds(fix)
        if step < self.max_step - 1:
            self.fix_angle_single(0, 1, 2, fix)
            if self.species.atom[self.instance[2]] == 'O': fval = 1.35
            if self.species.atom[self.instance[2]] == 'C': fval = 1.0
            self.set_bond(0, 1, -999, change, step=step+1, stmax=self.max_step-1, findist=fval, geom=geom)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release

def abstraction_align(startgeom, instance, fragnatom):
    """
    Special function to align the fragments before any optimization.
    Maybe later this can be generalized for other bimolecular reactions.
    Then this has to be moved to reac_General or perhaps modify_geom.
    startgeom: starting geometry to be aligned
    instance: atoms in the instance
    fragnatom: the number of atoms in the first fragment as defined in the original stationary point(!)
    The first fragment is assumed to be the abstractor for now. Need to make sure we can do the reverse as well.
    """

    # align O----H, O is at origin, H is on +Z axis
    g0 = copy.deepcopy(geometry.translate_and_rotate(startgeom, instance[0], instance[1]))
    # align H-C, H is at origin, C is on +Z axis at an angle
    g1 = copy.deepcopy(geometry.translate_and_rotate(startgeom, instance[1], instance[2]))
    for i in range(len(g1)):
        g1[i] = geometry.rotate_atom(g1[i], [0, 1, 0], 20. / 180. * np.pi)
    g1[:, 2] += 2.  # shift the whole thing in z direction
    geom = np.concatenate((g0[:fragnatom], g1[fragnatom:]))  # does not work for reverse
    return geom
