from kinbot.reac_General import GeneralReac
from kinbot import geometry

class R12Insertion(GeneralReac):
    max_step = 12
    scan = 0
    skip = 0
    family_name = 'r12insertion'
    

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        self.fix_bonds(fix)

        if step < self.max_step:
            
            fval = [1.67, 2.2, 1.9]
            if self.species.atom[self.instance[0]] == 'H':
                fval = [1.7, 1.09, 2.2]
            if self.species.atom[self.instance[2]] == 'H': 
                fval = [2.2, 1.09, 1.7]
            if self.species.atom[self.instance[0]] == 'O' or self.species.atom[self.instance[2]] == 'O': 
                fval[1] = 1.8

            self.set_bond(0, 1, -999, change, step=step+1, stmax=self.max_step, findist=fval[0], geom=geom)
            self.set_bond(1, 2, -999, change, step=step+1, stmax=self.max_step, findist=fval[1], geom=geom)
            self.set_bond(2, 0, -999, change, step=step+1, stmax=self.max_step, findist=fval[2], geom=geom)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
