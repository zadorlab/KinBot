from kinbot.reac_General import GeneralReac

class R22Cycloaddition(GeneralReac):
    max_step = 1
    scan = 0
    skip = 0
    family_name = 'r22cycloaddition'
    

    def get_constraints(self,step,geom):
        fix = []
        change = []
        release = []

        self.fix_bonds(fix)

        if step == 0:
            fval = 2.2
            self.set_bond(0, 1, fval, change)
            self.set_bond(2, 3, fval, change)
            self.fix_angle_single(0, 1, 2, fix) 
            
        self.clean_constraints(change, fix) 
        
        return step, fix, change, release
