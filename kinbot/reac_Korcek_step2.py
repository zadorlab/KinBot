from kinbot.reac_General import GeneralReac
from kinbot import geometry

class KorcekStep2(GeneralReac):
    max_step = 12
    scan = 0
    skip = 0
    

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        self.fix_bonds(fix)
        if step < 12:
            
            if step < 10:
                dih = geometry.calc_dihedral(geom[self.instance[-4]], geom[self.instance[-3]], geom[self.instance[-2]], geom[self.instance[-1]])[0]
                if np.abs(dih) < 160:
                    #move the dihedral to 160 degrees in 10 steps
                    frac = 1. / (10 - step + 0.)
                    new_dih = dih + frac * (160. - dih)
                    constraint = [self.instance[-4] + 1,self.instance[-3] + 1,self.instance[-2] + 1,self.instance[-1] + 1,new_dih]
                    change.append(constraint)
            
            fval = [2.0,2.0,1.8,1.8]
            if self.species.atom[self.instance[-1]] == 'H':
                fval[2] = 1.35
                fval[3] = 1.35
            val = geometry.new_bond_length(self.species,self.instance[0],self.instance[1],step+1,12,fval[0],geom)
            self.set_bond(0, 1, val, change)
            
            val = geometry.new_bond_length(self.species,self.instance[2],self.instance[3],step+1,12,fval[1],geom)
            self.set_bond(2, 3, val, change)
            
            if self.species.bond[self.instance[-1]][self.instance[-2]] == 1:
                val = geometry.new_bond_length(self.species,self.instance[-1],self.instance[-2],step+1,12,fval[2],geom)
                self.set_bond(-1, -2, val, change)
                #else do not change this bond length, the bond needs to stay and just change in order
            
            val = geometry.new_bond_length(self.species,self.instance[-1],self.instance[3],step+1,12,fval[3],geom)
            self.set_bond(-1, 3, val, change)  #todo: larger rings, this only work for 5 membered rings

        self.clean_constraints(change, fix) 
        
        return step, fix, change, release
