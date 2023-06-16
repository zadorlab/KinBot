from kinbot.reac_General import GeneralReac


class IntraOHMigration(GeneralReac):
    max_step = 14
    scan = 0
    skip = 1
    dihstep = max_step - 2
    family_name = 'intraohmigration'
    

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
            
            if self.species.atom[self.instance[0]] == 'C':
                fval1 = 1.5
                fval2 = 2.1
            else:
                fval1 = 1.7 
                fval2 = 2.0

            self.set_bond(0, -1, fval1, change)
            self.set_bond(-2, -1, fval2, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
