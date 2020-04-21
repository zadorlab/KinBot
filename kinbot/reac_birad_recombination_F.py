from kinbot import geometry
from kinbot.reac_General import GeneralReac

class BiradRecombinationF(GeneralReac):
    max_step = 22
    scan = 0
    skip = 0
    dihstep = 12
    

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        if step < self.self.max_step:
            self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change, step)

        elif step < max_step:
            self.release_dihedrals(release)
                
            fval = 2.0
            val = geometry.new_bond_length(self.species, self.instance[0], self.instance[-1], step-11, 10, fval, geom)
            self.set_bond(0, -1, val, change)
        
        self.clean_constraints(change, fix)
        
        return step, fix, change, release
