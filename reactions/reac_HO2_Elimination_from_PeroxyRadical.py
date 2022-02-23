from kinbot.reac_General import GeneralReac


class HO2Elimination(GeneralReac):
    max_step = 14
    scan = 0
    skip = 0
    dihstep = max_step - 2
    family_name = 'ho2elimination'


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
            self.set_bond(0, -1, fval, change)
            
            fval = 1.3
            self.set_bond(0, 1, fval, change)
            
            fval = 2.0
            #self.set_bont(2, 1, fval)
            self.set_bond(2, 3, fval, change)
        
        self.clean_constraints(change, fix)
        
        return step, fix, change, release
