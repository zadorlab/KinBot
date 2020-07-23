from kinbot.reac_General import GeneralReac
from kinbot import geometry
from kinbot.parameters import Parameters
import numpy as np

class BarrierlessSaddle(GeneralReac):
    par = Parameters()
    max_step = par.par['scan_step']
    scan = 1
    skip = 0
    

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []

        if step < self.max_step:
            if step == 0:
                delta = par.par['barrierless_saddle_start'] - np.linalg.norm(geom[self.instance[0]] - geom[self.instance[1]])
            else:
                delta = par.par['barrierless_saddle_step']
            val = np.linalg.norm(geom[self.instance[0]] - geom[self.instance[1]]) + delta
            self.set_bond(0, 1, val, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
