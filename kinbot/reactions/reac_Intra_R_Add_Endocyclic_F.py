from kinbot import geometry
from kinbot.reac_General import GeneralReac

class IntraRAddEndocyclicF(GeneralReac):
    max_step = 23
    scan = 0
    skip = 1
    dihstep = 13
    family_name = 'intraraddendocyclicf'
    invert = None


    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []

        if step < self.max_step:
            self.fix_bonds(fix)

        if step == 0:
            hatom = -999  # the atom that is in the way
            head = -999  # the head atom
            axis = -999  # the other axis atom
            lever = -999  # the one on the chain before the head
            # target atom is in a cycle
            if self.species.cycle[self.instance[-1]] == 1:
                # attacking atom is not in a cycle
                if self.species.cycle[self.instance[0]] == 0:
                    # has to be at least 2 removed from (the) cycle
                    if self.species.cycle[self.instance[1]] == 0:
                        # need to find head atom
                        for i, h in enumerate(self.instance):
                            if self.species.cycle[h] == 1:
                                head = h
                                break
                        if head != -999:
                            axis = self.instance[i-1]
                            lever = self.instance[i-2]
                        # neighbors of head atom
                        neigh = []
                        for at in range(self.species.natom):
                            if self.species.bond[head][at]:
                                neigh.append(at)  # index list of neighbors
                        if len(neigh) == 4:  # need 4 neighbors
                            incyc = 0
                            for n in neigh:  # need two in cycles
                                if self.species.cycle[n] == 1:
                                    incyc += 1  
                            if incyc == 2:
                                for n in neigh:
                                    if self.species.cycle[n] == 0 and (n not in self.instance):
                                        hatom = n
            if hatom != -999 and head != -999 and axis != -999 and lever != -999:
                self.set_dihedral_single(hatom, head, axis, lever, 180., change)
                self.invert = [hatom + 1, head + 1, axis + 1, lever + 1]  # one-indexing
        elif step < self.dihstep:
            self.instance = self.filter_cycle(self.species, self.instance)
            self.set_dihedrals(change, step)
            if self.invert is not None:
                newchange = []
                for ch in change:
                    if ch[1] in self.invert[1:3] and ch[2] in self.invert[1:3]:
                        pass
                    else:
                        newchange.append(ch)
                change = newchange
            self.fix_angles(fix)
        elif step < self.max_step:
            self.release_angles(release)
            self.release_dihedrals(release)

            fval = 1.9
            if len(self.instance) > 4:
                fval = 2.2
            if self.species.atom[self.instance[0]] == 'O' and self.species.atom[self.instance[1]] == 'O' and self.species.atom[self.instance[-1]] == 'C':
                fval = 1.9
            self.set_bond(0, -1, -999, change, step=step-12, stmax=10, findist=fval, geom=geom)

        self.clean_constraints(change, fix)

        return step, fix, change, release
