from kinbot.reac_General import GeneralReac
from kinbot import geometry
import numpy as np

class KorcekStep2Even(GeneralReac):
    max_step = 12
    scan = 0
    skip = 0
    family_name = 'korcekstep2even'


    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

            bondbreak = np.array(self.instance[1:])
            bondbreak = np.append(bondbreak, self.instance[0])  # the O-O bond

            for ii in range(int(len(bondbreak) / 2)):
                fval = 1.8
                a = np.where(self.instance == bondbreak[2 * ii])[0][0]
                b = np.where(self.instance == bondbreak[2 * ii + 1])[0][0]
                finalAtomIndex = len(self.instance) - 1
                if a == finalAtomIndex and b == 0:
                    fval = 1.9  # Should be the O-O bond length at the TS
                self.set_bond(a, b, -999, change, step=step, stmax=self.max_step, findist=fval, geom=geom)

        self.clean_constraints(change, fix)

        return step, fix, change, release
