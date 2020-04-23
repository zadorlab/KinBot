from kinbot import geometry


class GeneralReac:
    max_step = -1
    scan = -1
    skip = -1
    dihstep = -1


    def __init__(self, species, qc, par, instance, instance_name):
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

        if self.scan:
            max_step = self.par.par['scan_step']


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

    def filter_cycle(self, species, instances):
        cycle_atoms = self.species.cycle_chain
        instances = instances
        filtered_instances = instances
        for c, cycle in enumerate(cycle_atoms):
            for a, atoms in enumerate(instances):
                j = a - 1
                k = a + 1
                if k < len(instances) and j >= 0:
                    ca = instances[j]
                    cb = instances[a]
                    cc = instances[k]
                    if ca in cycle_atoms[c] and cb in cycle_atoms[c] and cc in cycle_atoms[c]:
                        filtered_instances.pop(a)

        return filtered_instances

    def fix_bonds(self, fix):
        for i in range(self.species.natom - 1):
            for j in range(i + 1, self.species.natom):
                if self.species.bond[i][j] > 0:
                    fix.append([i + 1, j + 1])


    def fix_angles(self, fix):
        for angle in range(len(self.instance) - 2):
            constraint = []
            for i in range(3):
                constraint.append(self.instance[angle + i] + 1)
            fix.append(constraint)


    def fix_angle_single(self, a, b, c, fix):
        constraint = [self.instance[a] + 1,self.instance[b] + 1, self.instance[c] + 1]
        fix.append(constraint)


    def fix_dihedrals(self, fix):
        for dih in range(len(self.instance) - 3):
            f = []
            for i in range(4):
                f.append(self.instance[dih + i] + 1)
            fix.append(f)


    def set_bond(self, a, b, val, change, step=None, stmax=None, findist=None, geom=None):
        if step is not None:
            val = geometry.new_bond_length(self.species, self.instance[a], self.instance[b],
                    step, stmax, findist, geom)
        constraint = [self.instance[a] + 1, self.instance[b] + 1, val]
        change.append(constraint)


    def set_angles(self, change):
        for angle in range(len(self.instance) - 2):
            constraint = []
            for i in range(3):
                constraint.append(self.instance[angle + i] + 1)
            constraint.append(180. * (len(self.instance) - 2.) / len(self.instance))
            change.append(constraint)


    def set_angle_single(self, a, b, c, val, change):
        constraint = [self.instance[a] + 1, self.instance[b] + 1, self.instance[c] + 1]
        constraint.append(val)
        change.append(constraint)


    def set_dihedrals(self, change, step, cut=0):
        new_dihs = geometry.new_ring_dihedrals(self.species, self.instance, step, self.dihstep)
        for dih in range(len(self.instance) - 3 - cut):
            constraint = []
            for i in range(4):
                constraint.append(self.instance[dih + i] + 1)
            constraint.append(new_dihs[dih])
            change.append(constraint)


    def release_angles(self, release):
        for angle in range(len(self.instance) - 2):
            constraint = []
            for i in range(3):
                constraint.append(self.instance[angle + i] + 1)
            release.append(constraint)

    
    def release_dihedrals(self, release, start=0):
        for dih in range(start, len(self.instance) - 3):  
            constraint = []
            for i in range(4):
                constraint.append(self.instance[dih + i] + 1)
            release.append(constraint)
