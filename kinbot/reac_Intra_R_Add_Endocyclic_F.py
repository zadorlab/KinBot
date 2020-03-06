from kinbot import geometry
from kinbot.reac_General import GeneralReac

class IntraRAddEndocyclicF(GeneralReac):
    max_step = 14
    scan = 0
    skip = 1
    dihstep = 12


    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change)
        elif step == self.dihstep:
            if len(self.instance) > 3:
                self.fix_dihedrals(fix)
                self.set_angles(change)
        elif step == self.dishtep + 1:
            self.release_angles(release)
            self.releave_dihedrals(release)

            if len(self.instance) > 4:
                fval = 2.2
            else:
                fval = 2.0
            self.set_bond(0, -1, fval, change)

        self.clean_constraints(change, fix)

        return step, fix, change, release
