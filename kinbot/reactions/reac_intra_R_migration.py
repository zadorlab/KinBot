from ase.data import atomic_numbers, covalent_radii

from kinbot.reac_General import GeneralReac
from kinbot import constants


class IntraRMigration(GeneralReac):
    max_step = 14
    scan = 0
    skip = 1
    dihstep = max_step - 2
    family_name = 'intrarmigration'
   

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change, step)

        elif step == self.dihstep:
            if len(self.instance) > 3:
                self.fix_dihedrals(fix)
                self.set_angles(change)
            else:
                fval = 1.35
                if self.species.atom[self.instance[0]] == 'O': fval = 1.2
                self.set_bond(0, -1, fval, change)
                
                fval = 1.35
                if self.species.atom[self.instance[-2]] == 'O': fval = 1.2
                self.set_bond(-2, -1, fval, change)
                #  step = 13  TODO this was a typo I think
        elif step == 13:
            self.release_angles(release)
            self.release_dihedrals(release)
            
            try:
                fval = constants.st_bond[''.join(sorted(self.species.atom[self.instance[0]]+self.species.atom[self.instance[-1]]))]
            except KeyError:
                fval = 1.2 * (covalent_radii[atomic_numbers[self.species.atom[self.instance[0]]]] 
                              + covalent_radii[atomic_numbers[self.species.atom[self.instance[-1]]]])

            self.set_bond(0, -1, fval, change)
            
            try:
                fval = constants.st_bond[''.join(sorted(self.species.atom[self.instance[-2]]+self.species.atom[self.instance[-1]]))]
            except KeyError:
                fval = (covalent_radii[atomic_numbers[self.species.atom[self.instance[-2]]]] 
                        + covalent_radii[atomic_numbers[self.species.atom[self.instance[-1]]]])

            self.set_bond(-2, -1, fval, change)
        
        self.clean_constraints(change, fix)
        
        return step, fix, change, release
