from ase.data import atomic_numbers, covalent_radii

from kinbot import constants
from kinbot.reac_General import GeneralReac

class S12ShiftF(GeneralReac):
    max_step = 12
    scan = 0
    skip = 0
    family_name = 's12shiftf'
        

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)
           
        if step < self.max_step:
            try:
                final_dist1 = constants.st_bond[''.join(sorted(self.species.atom[self.instance[1]] 
                                                               + self.species.atom[self.instance[2]]))]
            except KeyError:
                final_dist1 = 1.2 * (covalent_radii[atomic_numbers[self.species.atom[self.instance[1]]]] 
                                     + covalent_radii[atomic_numbers[self.species.atom[self.instance[2]]]])                
            self.set_bond(1, 2, -999, change, step=step, stmax=self.max_step, findist=final_dist1, geom=geom)
            try:       
                final_dist2 = constants.st_bond[''.join(sorted(self.species.atom[self.instance[0]] 
                                                           + self.species.atom[self.instance[2]]))]
            except KeyError:
                final_dist2 = 1.2 * (covalent_radii[atomic_numbers[self.species.atom[self.instance[0]]]] 
                                     + covalent_radii[atomic_numbers[self.species.atom[self.instance[2]]]])
            self.set_bond(0, 2, -999, change, step=step, stmax=self.max_step, findist=final_dist2, geom=geom)

        self.clean_constraints(change, fix)        
        
        return step, fix, change, release
