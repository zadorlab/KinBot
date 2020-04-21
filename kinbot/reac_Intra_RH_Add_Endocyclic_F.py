from kinbot.reac_General import GeneralReac
from kinbot import geometry

class IntraRHAddEndoF(GeneralReac):
    max_step = 22
    scan = 0
    skip = 0
    dihstep = 12

    
    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change, step, 1)

        elif step < self.max_step:
            self.realease_dihedrals(release)

            fvals = [2.0, 1.4, 1.3, 1.8, 1.3]
            val1 = geometry.new_bond_length(self.species, self.instance[0], self.instance[-2], step-11, 10, fvals[0], geom)
            self.set_bond(0, -2, val1, change)
            val2 = geometry.new_bond_length(self.species, self.instance[0], self.instance[1], step-11, 10, fvals[1], geom)
            self.set_bond(0, 1, val2, change)
            val3 = geometry.new_bond_length(self.species, self.instance[1], self.instance[-1], step-11, 10, fvals[2], geom)
            self.set_bond(1, -1, val3, change)
            val4 = geometry.new_bond_length(self.species, self.instance[0], self.instance[-1], step-11, 10, fvals[3], geom)
            self.set_bond(0, -1, val4, change)
            val5 = geometry.new_bond_length(self.species, self.instance[-1], self.instance[-2], step-11, 10, fvals[4], geom)
            self.set_bond(-1, -2, val5, change)
        
        self.clean_constraints(change, fix)

        return step, fix, change, release
