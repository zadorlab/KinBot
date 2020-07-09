from kinbot.reac_General import GeneralReac
from kinbot import geometry
import numpy as np

class KorcekStep2Odd(GeneralReac):
    max_step = 12
    scan = 0
    skip = 0


    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []

        if step < self.max_step:
            self.fix_bonds(fix)

            bondbreak = np.array(self.instance[1:-3])
            bondbreak = np.delete(bondbreak(np.where(bondbreak == self.instance[-3])[0][0]))
            bondbreak = np.append(bondbreak, self.instance[0])  # the O-O bond

            for ati in bondbreak[1::2]:
                for atj in np.roll(bondbreak[::2], int(len(bondbreak) / 2 - 1):
                    ii = np.where(self.instance, ati)[0][0]
                    jj = np.where(self.instance, atj)[0][0]
                    fval = 2.0  # to be refined based on atom types
                    self.set_bond(ii, jj, -999, change, step=step, stmax=self.max_step, findist=fval, geom=geom)

            fval = 1.35
            self.set_bond(-3, -1, -999, change, step=step, stmax=self.max_step, findist=fval, geom=geom)

        if step == self.max_step:
            set_angle_single(-1, -3, -2, 80., change):

        self.clean_constraints(change, fix)

        return step, fix, change, release
