from kinbot import geometry

class IntraRAddExoTetCyclicF(GeneralReac):
    max_step = 22
    scan = 0
    skip = 1
    dihstep = 12
    

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change, step, cut)

            ldih = [] # constraint for the last dihedral, which needs to be 180 degrees    
            for i in range(4):
                ldih.append(self.instance[len(self.instance)-4+i] + 1)
            dih = geometry.calc_dihedral(geom[ldih[0] - 1], geom[ldih[1] - 1], geom[ldih[2] - 1], geom[ldih[3] - 1])[0]
            frac = 1./(12. - step)
            if dih < 0:
                new_dih = dih - frac * (180. + dih) 
                ldih.append(new_dih)
            else:
                new_dih = dih + frac * (180. - dih) 
                ldih.append(new_dih)
            change.append(ldih)

        elif step < self.max_step:
            self.release_dihedrals(release)

            fdist1 = constants.st_bond[''.join(sorted(self.species.atom[self.instance[0]] + self.species.atom[self.instance[-2]]))] * 1.0
            if ''.join(sorted(self.species.atom[self.instance[0]] + self.species.atom[self.instance[-2]])) == 'CO':
                fdist1 = 1.68
            ndist1 = geometry.new_bond_length(self.species, self.instance[0], self.instance[-2], step-11, 10, fdist1, geom)
            self.set_bond(0, -2, ndist1, change)
            
            fdist2 = constants.st_bond[''.join(sorted(self.species.atom[self.instance[-1]] + self.species.atom[self.instance[-2]]))] * 1.0
            if ''.join(sorted(self.species.atom[self.instance[-1]] + self.species.atom[self.instance[-2]])) == 'CO':
                fdist2 = 1.68
            ndist2 = geometry.new_bond_length(self.species, self.instance[-1], self.instance[-2], step-11, 10, fdist2, geom)
            self.set_bond(-1, -2, ndist2, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
