from kinbot.reac_General import GeneralReac


class IntraHMigrationSuprafacial(GeneralReac):
    max_step = 14
    scan = 0
    skip = 1
    dihstep = 12
    family_name = 'intrahmigrationsuprafacial'
    

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change, step)
            if step == 0 and len(self.instance) > 6:
                self.set_angles(change)
            else:
                self.fix_angles(fix)

        elif step == self.dihstep:
            if len(self.instance) > 3:
                self.fix_dihedrals(fix)
                self.set_angles(change)
            else:
                fval = 1.35
                if self.species.atom[self.instance[0]] == 'O': fval = 1.2
                self.set_bond(0, -1, fval, change)
                
                fval = 1.35
                if self.species.atom[self.instance[-2]] == 'O': fval = 1.2
                self.set_bond(-2, -1, fval, change)
                step += 1

        elif step == self.dihstep + 1:
            self.release_angles(release)
            self.release_dihedrals(release)
                
            fval = 1.35
            if self.species.atom[self.instance[0]] == 'O': fval = 1.2
            self.set_bond(0, -1, fval, change)
            
            fval = 1.0
            self.set_bond(1, -1, fval, change)
            
            fval = 1.35
            if self.species.atom[self.instance[-2]] == 'O': fval = 1.2
            self.set_bond(-2, -1, fval, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
