from kinbot.reac_General import GeneralReac
from kinbot import geometry

class KorcekStep2(GeneralReac):
    max_step = 12
    scan = 0
    skip = 0
    family_name = 'korcekstep2'
    

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)
            
            if step < 10:
                dih = geometry.calc_dihedral(geom[self.instance[-4]], geom[self.instance[-3]], geom[self.instance[-2]], geom[self.instance[-1]])[0]
                if abs(dih) < 160:
                    #move the dihedral to 160 degrees in 10 steps
                    frac = 1. / (10 - step + 0.)
                    new_dih = dih + frac * (160. - dih)
                    constraint = [self.instance[-4] + 1, self.instance[-3] + 1, self.instance[-2] + 1, self.instance[-1] + 1, new_dih]
                    change.append(constraint)
            
            fval = [2.0, 2.0, 1.8, 1.8]
            if self.species.atom[self.instance[-1]] == 'H':
                fval[2] = 1.35
                fval[3] = 1.35

            self.set_bond(0, 1, -999, change, step=step+1, stmax=self.max_step, findist=fval[0], geom=geom)
            self.set_bond(2, 3, -999, change, step=step+1, stmax=self.max_step, findist=fval[1], geom=geom)
            
            if self.species.bond[self.instance[-1]][self.instance[-2]] == 1:
                self.set_bond(-1, -2, -999, change, step=step+1, stmax=self.max_step, findist=fval[2], geom=geom)
                #else do not change this bond length, the bond needs to stay and just change in order
            
            self.set_bond(-1, 3, -999, change, step=step+1, stmax=self.max_step, findist=fval[3], geom=geom)
            # TODO larger rings, this only work for 5 membered rings

        self.clean_constraints(change, fix) 
        
        return step, fix, change, release
