from kinbot.reac_General import GeneralReac
from kinbot import geometry

class R13InsertionCO2(GeneralReac):
    max_step = 22
    scan = 0
    skip = 0
    dihstep = 12
    family_name = 'r13insertionco2'

    
    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change, step)

        elif step < self.max_step:
            self.release_dihedrals(release)
                
            fval = [2.0, 1.3, 2.0, 2.0]
            if self.species.atom[self.instance[-1]] == 'H':
                fval[2] = 1.35
                fval[3] = 1.35
            
            self.set_bond(0, 1, -999, change, step=step-11, stmax=10, findist=fval[0], geom=geom)
            self.set_bond(1, 2, -999, change, step=step-11, stmax=10, findist=fval[1], geom=geom)
            self.set_bond(2, 3, -999, change, step=step-11, stmax=10, findist=fval[2], geom=geom)
            self.set_bond(3, 0, -999, change, step=step-11, stmax=10, findist=fval[3], geom=geom)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
