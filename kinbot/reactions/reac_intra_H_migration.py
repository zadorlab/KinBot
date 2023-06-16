from kinbot.reac_General import GeneralReac


class IntraHMigration(GeneralReac):
    max_step = 15
    scan = 0
    skip = 1
    dihstep = 13
    family_name = 'intrahmigration' 

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []

        if step < self.max_step:
            self.fix_bonds(fix)

        if step == 0:
            # ins[0]: attacker, ins[-1]: H atom
            hatom = -999  # the atom that is in the way
            head = -999  # the head atom
            axis = -999  # the other axis atom
            lever = -999  # the one on the chain before the head
            # target atom is directly attached to a cycle, -1 is H, -2 is the next atom
            if self.species.cycle[self.instance[-2]] == 1:
                # attacking atom, e.g., O*, is not in a cycle
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
            if step == 0 and len(self.instance) > 6:
                self.set_angles(change)
            else:
                self.fix_angles(fix)

        elif step == self.dihstep:
            if len(self.instance) > 3:
                self.fix_dihedrals(fix)
                self.set_angles(change)
            else:
                if self.species.atom[self.instance[0]] == 'O': 
                    fval = 1.2
                elif self.species.atom[self.instance[0]] == 'S': 
                    fval = 1.5
                else: 
                    fval = 1.35
                self.set_bond(0, -1, fval, change)
               
                if self.species.atom[self.instance[-2]] == 'O': 
                    fval = 1.2
                elif self.species.atom[self.instance[-2]] == 'S': 
                    fval = 1.5
                else: 
                    fval = 1.35
                self.set_bond(-2, -1, fval, change)
                step += 1

        elif step == self.dihstep + 1:
            self.release_angles(release)
            self.release_dihedrals(release)
               
            if self.species.atom[self.instance[0]] == 'O': 
                fval = 1.2
            elif self.species.atom[self.instance[0]] == 'S': 
                fval = 1.5
            else: 
                fval = 1.35
            self.set_bond(0, -1, fval, change)

            if self.species.atom[self.instance[-2]] == 'O': 
                fval = 1.2
            elif self.species.atom[self.instance[-2]] == 'S': 
                fval = 1.5
            else: 
                fval = 1.35
            self.set_bond(-2, -1, fval, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
