from kinbot.reac_General import GeneralReac
from kinbot import geometry
import numpy as np

class KorcekStep2Odd(GeneralReac):
    max_step = 12
    scan = 0
    skip = 0
    family_name = 'korcekstep2odd'


    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []

        if step < self.max_step:
            self.fix_bonds(fix)

            bondbreak = np.array(self.instance[1:-3])
            bondbreak = np.delete(bondbreak, np.where(bondbreak == self.instance[-3])[0][0])
            bondbreak = np.append(bondbreak, self.instance[0])  # the O-O bond

            for ii in range(int(len(bondbreak) / 2)):
                fval = 1.8
                a = np.where(self.instance == bondbreak[2 * ii])[0][0]
                b = np.where(self.instance == bondbreak[2 * ii + 1])[0][0]
                self.set_bond(a, b, -999, change, step=step, stmax=self.max_step, findist=fval, geom=geom)

            fval = 1.35
            self.set_bond(-3, -1, -999, change, step=step, stmax=self.max_step, findist=fval, geom=geom)
            self.set_bond(-2, -1, -999, change, step=step, stmax=self.max_step, findist=fval, geom=geom)

        self.clean_constraints(change, fix)

        return step, fix, change, release
