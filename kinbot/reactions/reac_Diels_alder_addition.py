from kinbot.reac_General import GeneralReac


class DielsAlder(GeneralReac):
    max_step = 1
    scan = 0
    skip = 0
    family_name = 'dielsalder'


    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

        if step ==0 :
            fval = 2.2
            self.set_bond(2, 3, fval, change)
            self.set_bond(4, 5, fval, change)
        
        self.clean_constraints(change, fix)
        
        return step, fix, change, release
