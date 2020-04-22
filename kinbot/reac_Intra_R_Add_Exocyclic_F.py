from kinbot.reac_General import GeneralReac


class IntraRAddExocyclicF(GeneralReac):
    max_step = 14
    scan = 0
    skip = 1
    dihstep = max_step - 2


    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change, fix, cut=1)

        elif step == self.dihstep:
            self.fix_dihedrals(fix)
            self.set_angles(change)

        elif step == self.dihstep + 1:
            self.release_angles(release)
            self.release_dihedrals(release)
                
            fval = 2.2
            self.set_bond(0, -2, fval, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
