from kinbot import geometry

class IntraDisproportionationF(GeneralReac):
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
            self.set_dihedrals(change, step)

        elif step < self.max_step:
            self.release_dihedrals(release)
                
            fval = 1.35

            val = geometry.new_bond_length(self.species, self.instance[0], self.instance[-1], step-11, 10, fval, geom)
            self.set_bond(0, -1, val)
            
            val = geometry.new_bond_length(self.species, self.instance[-2], self.instance[-1], step-11, 10, fval, geom)
            self.set_bond(-2, -1, val)
        
        self.clean_constraints(change, fix)
        
        return step, fix, change, release
