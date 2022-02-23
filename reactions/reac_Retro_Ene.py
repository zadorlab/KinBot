from kinbot.reac_General import GeneralReac
from kinbot import geometry

class RetroEne(GeneralReac):
    max_step = 22
    scan = 0
    skip = 0
    dihstep = 12
    family_name = 'retorene'
    

    def get_constraints(self,step,geom):
        fix = []
        change = []
        release = []

        self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change, step)

        elif step < self.max_step:
            self.release_dihedrals(release)
            
            self.set_bond(0, -1, -999, change, step=step-11, stmax=10, findist=1.35, geom=geom)
            self.set_bond(-1, -2, -999, change, step=step-11, stmax=10, findist=1.35, geom=geom)
            self.set_bond(2, 3, -999, change, step=step-11, stmax=10, findist=2.0, geom=geom)
            
        self.clean_constraints(change, fix)
        
        return step, fix, change, release
