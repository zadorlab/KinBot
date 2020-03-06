from kinbot import geometry
from kinbot.reac_General import GeneralReac

class betadelta(GeneralReac):
    max_step = 12
    scan = 0
    skip = 0
   

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bond(fix)
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
            
            val = geometry.new_bond_length(self.species,self.instance[1], self.instance[2], step, self.max_step, final_dist, geom)
            self.set_bond(1, 2, val, change)

            # breaking the delta bond
            final_dist = 1.7 # this might need to be refined
            val = geometry.new_bond_length(self.species,self.instance[3], self.instance[4], step, self.max_step, final_dist, geom)
            self.set_bond(3, 4, val, change)

        self.clean_contraints(change, fix)
        
        return step, fix, change, release
