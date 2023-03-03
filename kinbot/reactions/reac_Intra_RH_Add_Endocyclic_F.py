from kinbot.reac_General import GeneralReac
from kinbot import geometry

class IntraRHAddEndoF(GeneralReac):
    max_step = 22
    scan = 0
    skip = 0
    dihstep = 12
    family_name = 'intrarhaddendof'

    
    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

        if step < self.dihstep:
            # this might be necessary
            #self.instance = self.filter_cycle(self.species, self.instance)
            self.set_dihedrals(change, step, 1)
            self.fix_angles(fix)

        elif step < self.max_step:
            self.release_dihedrals(release)

            fvals = [2.0, 1.4, 1.3, 1.8, 1.3]
            self.set_bond(0, -2, -999, change, step=step-11, stmax=10, findist=fvals[0], geom=geom)
            self.set_bond(0, 1, -999, change, step=step-11, stmax=10, findist=fvals[1], geom=geom)
            self.set_bond(1, -1, -999, change, step=step-11, stmax=10, findist=fvals[2], geom=geom)
            self.set_bond(0, -1, -999, change, step=step-11, stmax=10, findist=fvals[3], geom=geom)
            self.set_bond(-1, -2, -999, change, step=step-11, stmax=10, findist=fvals[4], geom=geom)
        
        self.clean_constraints(change, fix)

        return step, fix, change, release
