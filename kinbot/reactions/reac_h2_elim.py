from kinbot.reac_General import GeneralReac
from kinbot import geometry

class H2Elim(GeneralReac):
    max_step = 22
    scan = 0
    skip = 1
    dihstep = 12
    family_name = 'h2elim'
    
    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        self.fix_bonds(fix)

        if step < self.max_step:
            self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change, step)
        elif step == self.dihstep:
            self.fix_dihedrals(fix)
            self.set_angles(change) 
        elif step < self.max_step:
            self.release_angles(release)
            self.release_dihedrals(release)
            self.set_bond(0, 1, -999, change, step=step-13, stmax=9, findist=1.35, geom=geom)
            self.set_bond(-2, -1, -999, change, step=step-13, stmax=9, findist=1.35, geom=geom)
            self.set_bond(0, -1, -999, change, step=step-13, stmax=9, findist=0.9, geom=geom)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
