from kinbot.reac_General import GeneralReac
from kinbot import geometry

class HS(GeneralReac):
    max_step = 1
    scan = 0
    skip = 0
    
    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        self.fix_bonds(fix)
        self.fix_angles(fix)
        self.fix_dihedrals(fix)

        if step <= self.max_step:
            self.set_bond(0, 1, 10., change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
