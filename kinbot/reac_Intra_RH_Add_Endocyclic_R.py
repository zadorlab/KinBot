from kinbot.reac_General import GeneralReac
from kinbot import geometry

class IntraRHAddEndoR(GeneralReac):
    max_step = 12
    scan = 0
    skip = 0
    dihstep = 12
    

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

        if step < self.dihstep:
            final_dist = [2.0, 1.3, 1.35, 1.3]
            for i in range(len(self.instance)):
                if i == len(self.instance) - 1:
                    j = 0
                else:
                    j = i + 1
                val = geometry.new_bond_length(self.species,self.instance[i],self.instance[j],step,12,final_dist[i],geom)
                self.set_bond(i, j, val, change)

        self.clean_constraints(change, fix)

        return step, fix, change, release
