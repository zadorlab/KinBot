from kinbot.reac_General import GeneralReac
import numpy as np
from kinbot.parameters import Parameters

class R12Cycloaddition(GeneralReac):
    scan = 1
    skip = 0
    par = Parameters()
    max_step = par.par['scan_step']    
    mp2 = 1

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        self.fix_bonds(fix)
        if step < self.max_step:
            val = np.linalg.norm(geom[self.instance[0]] - geom[self.instance[2]]) + 0.02
            self.set_bond(0, 2, val, change)
            
            val = np.linalg.norm(geom[self.instance[1]] - geom[self.instance[2]]) + 0.02
            self.set_bond(1, 2, val, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
