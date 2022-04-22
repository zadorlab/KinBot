import numpy as np
from kinbot.reac_General import GeneralReac
from kinbot import geometry

class RAdditionMultipleBond(GeneralReac):
    max_step = 12
    scan = 0
    skip = 0
    family_name = 'radditionmultiplebond'
    

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []

        self.fix_bonds(fix)
        self.fix_angles(fix)
        self.fix_dihedrals(fix)

        if step < self.max_step:
#            # verify if anchor atoms have more than one neighbor and 
#            # change the dihedral to 90 degrees in that case
            neigh = [i for i, ni in enumerate(self.species.bond[self.instance[0]]) if ni > 0]
            if len(neigh) > 1:
                for ni in neigh:
                    if not ni in self.instance:
                        change.append([ni + 1, self.instance[0] + 1, self.instance[1] + 1, self.instance[2] + 1, 90.])
                        break
            # anchors for middle atom, improper dihedral
            neigh = [i for i, ni in enumerate(self.species.bond[self.instance[1]]) if ni > 0]
            if len(neigh) > 2:
                for ni in neigh:
                    if not ni in self.instance:
                        change.append([ni + 1, self.instance[1] + 1, self.instance[0] + 1, self.instance[2] + 1, 90.])
                        self.set_angle_single(ni, self.instance[1], self.instance[2], 90., change)
                        break

            self.set_angle_single(self.instance[0], self.instance[1], self.instance[2], 90., change)

            final_dist = 1.42
            if self.species.atom[self.instance[0]] == 'C' and self.species.atom[self.instance[1]] == 'C' and self.species.atom[self.instance[2]] == 'C':
                final_dist = 2.20
            if self.species.atom[self.instance[0]] == 'C' and self.species.atom[self.instance[1]] == 'C' and self.species.atom[self.instance[2]] == 'H':
                final_dist = 1.79
            if self.species.atom[self.instance[0]] == 'C' and self.species.atom[self.instance[1]] == 'C' and self.species.atom[self.instance[2]] == 'O':
                final_dist = 2.04
            if self.species.atom[self.instance[0]] == 'O' and self.species.atom[self.instance[1]] == 'C' and self.species.atom[self.instance[2]] == 'C':
                final_dist = 2.12
            if self.species.atom[self.instance[0]] == 'O' and self.species.atom[self.instance[1]] == 'C' and self.species.atom[self.instance[2]] == 'H':
                final_dist = 1.84
            if self.species.atom[self.instance[0]] == 'O' and self.species.atom[self.instance[1]] == 'C' and self.species.atom[self.instance[2]] == 'O':
                final_dist = 2.04 #TODO: verify if this value is OK
            if self.species.atom[self.instance[0]] == 'C' and self.species.atom[self.instance[1]] == 'O' and self.species.atom[self.instance[2]] == 'C':
                final_dist = 2.04 #TODO: verify if this value is OK
            if self.species.atom[self.instance[0]] == 'C' and self.species.atom[self.instance[1]] == 'O' and self.species.atom[self.instance[2]] == 'H':
                final_dist = 1.42
            if self.species.atom[self.instance[0]] == 'C' and self.species.atom[self.instance[1]] == 'O' and self.species.atom[self.instance[2]] == 'O':
                final_dist = 2.04 #TODO: verify if this value is OK
            if self.species.atom[self.instance[0]] == 'O' and self.species.atom[self.instance[1]] == 'O' and self.species.atom[self.instance[2]] == 'C':
                final_dist = 2.04 #TODO: verify if this value is OK
            
            self.set_bond(1, 2, -999, change, step=step, stmax=self.max_step, findist=final_dist, geom=geom)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
