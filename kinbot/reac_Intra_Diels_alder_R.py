from kinbot import geometry
from kinbot.reac_General import GeneralReac

class IntraDielsAlder(GeneralReac):
    max_step = 22
    scan = 0
    skip = 0
    dihstep = 12
    

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)
        if step < self.dihstep:
            self.set_dihedrals(change)
        elif step < self.max_step:
            self.release_dihedrals(release)

            first_dih = [self.instance[i] + 1 for i in range(4)]
            if step < 18: # make sure that the forming double bond stays in trans instead of moving to cis
                fix.append(first_dih)
            else:
                release.append(first_dih)

            fval = 2.2
            val = geometry.new_bond_length(self.species, self.instance[0], self.instance[-1], step-11, 10, fval,geom)
            self.set_bond(0, -1, val, change)
            
            fval = 1.8
            val = geometry.new_bond_length(self.species, self.instance[-2], self.instance[3], step-11, 10, fval, geom)
            self.set_bond(-2, 3, val, change)

        self.clean_constraints(change, fix)

        return step, fix, change, release
