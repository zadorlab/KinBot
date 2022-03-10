from kinbot.reac_General import GeneralReac
from kinbot import geometry
import numpy as np

class BarrierlessSaddle(GeneralReac):
    scan = 1
    skip = 0
    family_name = 'barrierlesssaddle'
    
    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []

        if step < self.max_step:
            if step == 0:
                delta = self.par['barrierless_saddle_start'] - np.linalg.norm(geom[self.instance[0]] - geom[self.instance[1]])
            else:
                delta = self.par['barrierless_saddle_step']
            val = np.linalg.norm(geom[self.instance[0]] - geom[self.instance[1]]) + delta
            self.set_bond(0, 1, val, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
