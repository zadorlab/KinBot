from kinbot.reac_General import GeneralReac
from kinbot import geometry
import numpy as np

class BarrierlessSaddle(GeneralReac):
    max_step = 20
    scan = 1
    skip = 0
    

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []

        if step < self.max_step:
            if step == 0:
                delta = 0.5
            else:
                delta = 0.2
            val = np.linalg.norm(geom[self.instance[0]] - geom[self.instance[1]]) + delta
            self.set_bond(0, 1, val, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
