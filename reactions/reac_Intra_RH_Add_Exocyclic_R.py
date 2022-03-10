from kinbot.reac_General import GeneralReac
from kinbot import geometry

class IntraRHAddExoR(GeneralReac):
    max_step = 12
    scan = 0
    skip = 0
    family_name = 'intrarhaddexor'
    

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

        if step < self.max_step:
            final_dist = [2.0, 1.45, 1.4, 1.9]
            for i in range(len(self.instance)):
                if i == len(self.instance) - 1:
                    j = 0
                else:
                    j = i + 1
                self.set_bond(i, j, -999, change, step=step+1, stmax=self.max_step, findist=final_dist[i], geom=geom)

            f_dist = 1.3
            self.set_bond(1, 3, -999, change, step=step+1, stmax=self.max_step, findist=f_dist, geom=geom)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
