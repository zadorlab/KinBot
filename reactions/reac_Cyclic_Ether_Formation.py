from kinbot import geometry
from kinbot.reac_General import GeneralReac


class CyclicEtherFormation(GeneralReac):
    max_step = 14
    scan = 0
    skip = 1
    dihstep = 12
    family_name = 'cyclicetherformation'
   

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)
        if step < self.dihstep:
            self.set_dihedrals(change, step)
        elif step == self.dihstep:
            if len(self.instance) > 3:
                self.fix_dihedrals(fix)
                self.set_angles(change)
        elif step == self.dihstep + 1:
            self.release_angles(release)
            self.release_dihedrals(release)
            
            if len(self.instance) > 3:
                fval = 1.8
            else:
                fval = 1.5
            self.set_bond(0, -1, fval, change)
        self.clean_constraints(change, fix)

        return step, fix, change, release
