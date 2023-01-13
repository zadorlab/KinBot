from kinbot.reac_General import GeneralReac
from kinbot import geometry

class H2Elim(GeneralReac):
    max_step = 12
    scan = 0
    skip = 0
    family_name = 'h2elim'
    
    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        self.fix_bonds(fix)

        if step <= self.max_step:
            self.set_bond(0, 1, -999, change, step=step, stmax=self.max_step, findist=1.35, geom=geom)
            self.set_bond(2, 3, -999, change, step=step, stmax=self.max_step, findist=1.35, geom=geom)
            self.set_bond(0, 3, -999, change, step=step, stmax=self.max_step, findist=0.9, geom=geom)
            #self.set_dihedrals(change, step)

#        elif step == self.max_step:
#            self.release_dihedrals(release)
#                
#            self.set_bond(0, 1, 1.35, change)
#            self.set_bond(2, 3, 1.35, change)
#            self.set_bond(0, 3, 0.9, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
