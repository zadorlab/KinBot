from kinbot import geometry
from kinbot.reac_General import GeneralReac

class IntraDielsAlder(GeneralReac):
    max_step = 22
    scan = 0
    skip = 0
    dihstep = 12
    family_name = 'intradielsalder'
    

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)
        if step < self.dihstep:
            self.set_dihedrals(change, step)
        elif step < self.max_step:
            self.release_dihedrals(release, start=1)

            first_dih = [self.instance[i] + 1 for i in range(4)]
            if step < 18: # make sure that the forming double bond stays in trans instead of moving to cis
                fix.append(first_dih)
            else:
                release.append(first_dih)

            fval = 2.2
            self.set_bond(0, -1, -999, change, step=step-11, stmax=10, findist=fval, geom=geom)
            
            fval = 1.8
            self.set_bond(-2, 3, -999, change, step=step-11, stmax=10, findist=fval, geom=geom)

        self.clean_constraints(change, fix)

        return step, fix, change, release
