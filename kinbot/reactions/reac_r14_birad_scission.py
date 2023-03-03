from kinbot.reac_General import GeneralReac
import numpy as np

class R14BiradScission(GeneralReac):
    scan = 1
    skip = 0
    mp2 = 1
    family_name = 'r14iradscission'
    
    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []

        self.fix_bonds(fix)

        if step < self.max_step:
            val = np.linalg.norm(geom[self.instance[1]] - geom[self.instance[2]]) + 0.1
            self.set_bond(1, 2, val, change)
        
        self.clean_constraints(change, fix)
        
        return step, fix, change, release
