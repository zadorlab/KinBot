from kinbot import geometry


class IntraHMigration:
    max_step = 14
    scan = 0
    skip = 1
    dihstep = 12
  

    def __init__(self, species, qc,par, instance, instance_name):
        self.species = species
        self.ts = None
        self.products = []
        self.product_bonds = [] 
        self.broken_bonds = []
        self.formed_bonds = []
        
        self.ts_opt = None
        self.prod_opt = []
        
        self.qc = qc
        self.par = par
        
        self.instance = instance
        self.instance_name = instance_name


    def clean_constraints(self, change, fix):

        for c in change:
            if len(c) == 3:
                index = -1
                for i, fi in enumerate(fix):
                    if len(fi) == 2:
                        if sorted(fi) == sorted(c[:2]):
                            index = i
                if index > -1:
                    del fix[index]

        return change, fix


    def fix_bonds(self, fix):
        for i in range(self.species.natom - 1):
            for j in range(i + 1, self.species.natom):
                if self.species.bond[i][j] > 0:
                    fix.append([i + 1, j + 1])

        return fix


    def fix_angles(self, fix):
        for angle in range(len(self.instance) - 2):
            constraint = []
            for i in range(3):
                constraint.append(self.instance[angle+i] + 1)
            fix.append(constraint)

        return fix


    def fix_dihedrals(self, fix):
        for dih in range(len(self.instance)-3):
            f = []
            for i in range(4):
                f.append(self.instance[dih+i] + 1)
            fix.append(f)

        return fix


    def set_bond(self, atom1, atom2, fval, change):
        constraint = [self.instance[atom1] + 1,self.instance[atom2] + 1,fval]
        change.append(constraint)


    def set_angles(self, change):
        for angle in range(len(self.instance) - 2):
            constraint = []
            for i in range(3):
                constraint.append(self.instance[angle+i] + 1)
            constraint.append(180. * (len(self.instance) - 2.) / len(self.instance))
            change.append(constraint)

        return change


    def set_dihedrals(self, change, step):
        new_dihs = geometry.new_ring_dihedrals(self.species, self.instance, step, self.dihstep)
        for dih in range(len(self.instance) - 3):
            constraint = []
            for i in range(4):
                constraint.append(self.instance[dih+i] + 1)
            constraint.append(new_dihs[dih])
            change.append(constraint)

        return change


    def release_angles(self, release):
        for angle in range(len(self.instance)-2):
            constraint = []
            for i in range(3):
                constraint.append(self.instance[angle+i] + 1)
            release.append(constraint)

        return release

    
    def release_dihedrals(self, release):
        for dih in range(len(self.instance)-3):  
            constraint = []
            for i in range(4):
                constraint.append(self.instance[dih+i] + 1)
            release.append(constraint)

        return release


    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

        if step < self.dihstep:
           self.set_dihedrals(change, step)

           if step == 0 and len(self.instance) > 6:
               self.set_angles(change)
           else:
               self.fix_angles(fix)

        elif step == self.dihstep:
            if len(self.instance) > 3:
                self.fix_dihedrals(fix)
                self.set_angles(change)
            else:
                fval = 1.35
                if self.species.atom[self.instance[0]] == 'O': fval = 1.2
                self.set_bond(0, -1, fval, change)
               
                fval = 1.35
                if self.species.atom[self.instance[-2]] == 'O': fval = 1.2
                self.set_bond(-2, -1, fval, change)
                step += 1
        elif step == self.dihstep + 1:
            self.release_angles(release)
            self.release_dihedrals(release)
               
            fval = 1.35
            if self.species.atom[self.instance[0]] == 'O': fval = 1.2
            self.set_bond(0, -1, fval, change)

            fval = 1.35
            if self.species.atom[self.instance[-2]] == 'O': fval = 1.2
            self.set_bond(-2, -1, fval, change)


        self.clean_constraints(change, fix)
        
        return step, fix, change, release
