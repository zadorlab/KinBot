from kinbot.reac_General import GeneralReac
from kinbot import geometry

class R14CyclicBiradScission(GeneralReac):
    max_step = 22
    scan = 0
    skip = 0
    dihstep = 12
    family_name = 'r14cyclicbiradscission'
    

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []

        self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change, step)

        elif step == self.dihstep:  # originally was < 12, must be a typo
            self.release_dihedrals(release)
                
            fval = 1.8
            self.set_bond(0, -1, -999, change, step=step-11, stmax=10, findist=fval, geom=geom)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
