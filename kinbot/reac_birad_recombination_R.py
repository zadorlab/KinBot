class BiradRecombinationR(GeneralReac):
    max_step = self.par.par['scan_step']
    scan = 1
    skip = 0
    
    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

        if step < self.max_step:
            val = np.linalg.norm(geom[self.instance[0]] - geom[self.instance[1]]) + 0.1
            self.set_bond(0, 1, val, change)

        self.clean_constraints(change, fix)

        return step, fix, change, release
