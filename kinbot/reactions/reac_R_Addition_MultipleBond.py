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

# instance is a*-b-c, and the product is a=b + c
        if step < self.max_step:
#            # verify if anchor atoms have more than one neighbor and 
#            # change the dihedral to 90 degrees in that case
#            neigh = [i for i, ni in enumerate(self.species.bond[self.instance[0]]) if ni > 0]
#            if len(neigh) > 1:
#                for ni in neigh:
#                    if not ni in self.instance:
#                        change.append([ni + 1, self.instance[0] + 1, self.instance[1] + 1, self.instance[2] + 1, 90.])
#                        break
#            # anchors for middle atom, improper dihedral
#            neigh = [i for i, ni in enumerate(self.species.bond[self.instance[1]]) if ni > 0]
#            if len(neigh) > 2:
#                for ni in neigh:
#                    if not ni in self.instance:
#                        change.append([ni + 1, self.instance[1] + 1, self.instance[0] + 1, self.instance[2] + 1, 90.])
#                        self.set_angle_single(ni, self.instance[1], self.instance[2], 90., change)
#                        break

            self.set_angle_single(self.instance[0], self.instance[1], self.instance[2], 90., change)

            atom_pattern = f'{self.species.atom[self.instance[0]]}{self.species.atom[self.instance[1]]}{self.species.atom[self.instance[2]]}'

            match atom_pattern:
                case 'CCC':
                    final_dist = 2.20
                case 'CCH':
                    final_dist = 1.79
                case 'CCO':
                    final_dist = 2.04
                case 'OCC':
                    final_dist = 2.12
                case 'OCH':
                    final_dist = 1.84
                case 'OCO':
                    final_dist = 2.04 
                case 'COC':
                    final_dist = 2.04 
                case 'COH':
                    final_dist = 1.42
                case 'COO': 
                    final_dist = 2.04 
                case 'OOC':
                    final_dist = 2.04 
                case 'NOO':
                    final_dist = 1.8
                case 'NCH':
                    final_dist = 1.8
                case 'OON':
                    final_dist = 2.1 
                case 'CNC':
                    final_dist = 2.1 
                case 'CNH':
                    final_dist = 1.8
                case 'CCS':
                    final_dist = 2.5 
                case 'COS':
                    final_dist = 2.5 
                case 'CNS':
                    final_dist = 2.5 
                case 'CSC':
                    final_dist = 2.5 
                case 'CSH':
                    final_dist = 2.5 
                case _:
                    final_dist = 1.42

#                final_dist = 2.20
#                final_dist = 1.79
#                final_dist = 2.04
#                final_dist = 2.12
#                final_dist = 1.84
#                final_dist = 2.04 
#                final_dist = 2.04 
#                final_dist = 1.42
#                final_dist = 2.04 
#                final_dist = 2.04 
#                final_dist = 1.8
#                final_dist = 1.8
#                final_dist = 2.1 
#                final_dist = 2.5 
# 
            
            self.set_bond(1, 2, -999, change, step=step, stmax=self.max_step, findist=final_dist, geom=geom)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
