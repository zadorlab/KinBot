from kinbot.reac_General import GeneralReac


class IntraRAddExocyclicF(GeneralReac):
    max_step = 22
    scan = 0
    skip = 1
    dihstep = 12
    family_name = 'intraraddexocyclicf'


    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

        if step == 0: 
            self.instance = self.filter_cycle(self.species, self.instance, cut=1)
        if step < self.dihstep:
            #self.set_dihedrals(change, step, cut=1)
            self.set_dihedrals(change, step)
            self.fix_angles(fix)

        elif step < self.max_step:
            self.release_angles(release)
            self.release_dihedrals(release)
                
            fval = 2.2
            if self.species.atom[self.instance[0]] == 'O' and self.species.atom[self.instance[1]] == 'O' and self.species.atom[self.instance[-2]] == 'C':
                fval = 1.9
            self.set_bond(0, -1, -999, change, step=step-11, stmax=10, findist=fval, geom=geom)

        self.clean_constraints(change, fix)
        return step, fix, change, release
