from kinbot import geometry
from kinbot.reac_General import GeneralReac

class BetaDelta(GeneralReac):
    max_step = 12
    scan = 0
    skip = 0
    family_name = 'betadelta'
   

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)
        if step == 0:
            # verify if the radical atom has more than one neighbor and 
            # change the dihedral to 90 degrees in that case
            neigh = [i for i,ni in enumerate(self.species.bond[self.instance[0]]) if ni > 0]
            if len(neigh) > 1:
                for ni in neigh:
                    if not ni in self.instance:
                        change.append([ni + 1, self.instance[0] + 1, self.instance[1] + 1, self.instance[2] + 1, 90.])
                        break

        if step < self.max_step:
            final_dist = 1.42
            if self.species.atom[self.instance[0]] == 'C' and self.species.atom[self.instance[1]] == 'C' and self.species.atom[self.instance[2]] == 'C':
                final_dist = 1.6
            if self.species.atom[self.instance[0]] == 'C' and self.species.atom[self.instance[1]] == 'C' and self.species.atom[self.instance[2]] == 'H':
                final_dist = 1.6
            if self.species.atom[self.instance[0]] == 'C' and self.species.atom[self.instance[1]] == 'C' and self.species.atom[self.instance[2]] == 'O':
                final_dist = 1.6
            if self.species.atom[self.instance[0]] == 'O' and self.species.atom[self.instance[1]] == 'C' and self.species.atom[self.instance[2]] == 'C':
                final_dist = 1.6
            if self.species.atom[self.instance[0]] == 'O' and self.species.atom[self.instance[1]] == 'C' and self.species.atom[self.instance[2]] == 'H':
                final_dist = 1.6
            if self.species.atom[self.instance[0]] == 'O' and self.species.atom[self.instance[1]] == 'C' and self.species.atom[self.instance[2]] == 'O':
                final_dist = 1.6 
            if self.species.atom[self.instance[0]] == 'C' and self.species.atom[self.instance[1]] == 'O' and self.species.atom[self.instance[2]] == 'C':
                final_dist = 1.6 
            if self.species.atom[self.instance[0]] == 'C' and self.species.atom[self.instance[1]] == 'O' and self.species.atom[self.instance[2]] == 'H':
                final_dist = 1.6
            if self.species.atom[self.instance[0]] == 'C' and self.species.atom[self.instance[1]] == 'O' and self.species.atom[self.instance[2]] == 'O':
                final_dist = 1.6 
            if self.species.atom[self.instance[0]] == 'O' and self.species.atom[self.instance[1]] == 'O' and self.species.atom[self.instance[2]] == 'C':
                final_dist = 1.6 
            
            self.set_bond(1, 2, -999, change, step=step, stmax=self.max_step, findist=final_dist, geom=geom)

            # breaking the delta bond
            final_dist = 1.7 # this might need to be refined
            self.set_bond(3, 4, -999, change, step=step, stmax=self.max_step, findist=final_dist, geom=geom)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
