from kinbot.reac_General import GeneralReac
from kinbot import geometry
import numpy as np

class HS(GeneralReac):
    scan = 0
    skip = 0
    max_step = 0
    family_name = 'hs'
    
    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
