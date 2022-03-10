from kinbot.reac_General import GeneralReac
from kinbot.geometry import calc_dihedral


class IntraOHMigrationExocyclicF(GeneralReac):
    max_step = 14
    scan = 0
    skip = 1
    dihstep = max_step - 2
    family_name = 'intraohmigrationexocyclicf'

    def __init__(self, species, qc, par, instance, instance_name):  # special
        GeneralReac.__init__(self, species, qc, par, instance, instance_name)
        self.cistrans = self.instance[-1]  # (nominal) cis -1, trans -2
        self.double = self.instance[0]  # the dangling atom with double bond
        self.instance = self.instance[1:-1]  # cut off both
        self.dih0, _ = calc_dihedral(self.species.geom[self.double],  # the original dihedral
                                     self.species.geom[self.instance[0]], 
                                     self.species.geom[self.instance[1]],
                                     self.species.geom[self.instance[2]])

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)
            #if self.cistrans == -1:
            #    self.fix_dihedral_single(self.double, self.instance[0], self.instance[1], self.instance[2], fix)
            #else:
            if self.cistrans == -2:
                val = self.dih0 + 180. / self.max_step * step
                self.set_dihedral_single(self.double, self.instance[0], self.instance[1], self.instance[2], val, change)

        if step < self.dihstep:
            self.set_dihedrals(change, step)

        elif step == self.dihstep:
            self.fix_dihedrals(fix)
            #self.set_dihedral_single(self.instance[0], self.instance[-1], self.instance[-2], self.instance[-3], 0., change)
            self.set_angles(change) 

        elif step == self.dihstep + 1:
            self.release_angles(release)
            self.release_dihedrals(release)
            
            if self.species.atom[self.instance[0]] == 'C':
                fval1 = 1.5
                fval2 = 1.9
            else:
                fval1 = 1.7 
                fval2 = 2.0

            self.set_bond(0, -1, fval1, change)
            self.set_bond(-2, -1, fval2, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
