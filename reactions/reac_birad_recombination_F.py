from kinbot import geometry
from kinbot.reac_General import GeneralReac

class BiradRecombinationF(GeneralReac):
    max_step = 22
    scan = 0
    skip = 0
    dihstep = 12
    family_name = 'biradrecombinationf'
    

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        if step < self.self.max_step:
            self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change, step)

        elif step < self.max_step:
            self.release_dihedrals(release)
                
            fval = 2.0
            self.set_bond(0, -1, -999, change, step=step-11, stmax=10, findist=fval, geom=geom)
        
        self.clean_constraints(change, fix)
        
        return step, fix, change, release
