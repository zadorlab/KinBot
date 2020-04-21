class HO2Elimination(GeneralReac):
    max_step = 14
    scan = 0
    skip = 0
    dihstep = max_step - 2


    def get_constraints(self,step, geom):
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
                
            fval = 1.3
            self.set_bond(0, -1, fval)
            
            fval = 1.3
            self.set_bond(0, 1, fval)
            
            fval = 2.0
            self.set_bond(2, 1, fval)
        
        self.clean_constraints(change, fix)
        
        return step, fix, change, release
