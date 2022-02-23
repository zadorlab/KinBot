from kinbot.reac_General import GeneralReac

class KetoEnol(GeneralReac):
    max_step = 14
    scan = 0
    skip = 0
    dihstep = max_step - 2
    family_name = 'ketoenol'
    

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []

        if step < self.max_step:
            self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change, step)

        elif step == self.dihstep:
            self.fix_dihedrals(fix)
            self.set_angles(change)

        elif step == self.dihstep + 1:
            self.release_angles(release)
            self.release_dihedrals(release)
            
            fval = 1.35
            if self.species.atom[self.instance[-1]] != 'H': fval = 1.9
            self.set_bond(0, -1, fval, change) 
            self.set_bond(-2, -1, fval, change) 
            
        self.clean_constraints(change, fix)

        return step, fix, change, release
