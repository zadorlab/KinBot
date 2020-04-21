from kinbot.reac_General import GeneralReac
from kinbot import geometry

class R13InsertionRSR(GeneralReac):
    max_step = 22
    scan = 0
    skip = 0
    dihstep = 12
    

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change)

        elif step < self.max_step:
            self.release_dihedrals(release)
            
            fval = [2.0,1.45,2.0,2.0]
            if self.species.atom[self.instance[0]] == 'H':
                fval[0] = 1.3
                fval[3] = 1.3
            
            val = geometry.new_bond_length(self.species,self.instance[0],self.instance[1],step-11,10,fval[0],geom)
            self.set_bond(0, 1, val, change)
            
            val = geometry.new_bond_length(self.species,self.instance[1],self.instance[2],step-11,10,fval[1],geom)
            self.set_bond(1, 2, val, change)

            val = geometry.new_bond_length(self.species,self.instance[2],self.instance[3],step-11,10,fval[2],geom)
            self.set_bond(2, 3, val, change)

            val = geometry.new_bond_length(self.species,self.instance[3],self.instance[0],step-11,10,fval[3],geom)
            self.set_bond(3, 0, val, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
