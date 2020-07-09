from kinbot.reac_General import GeneralReac
from kinbot import geometry

class KorcekStep2Even(GeneralReac):
    max_step = 12
    scan = 0
    skip = 0


    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

            for ati in self.instance[1::2]:
                for atj in np.roll(self.instance[::2], int(len(self.instance) / 2 - 1)):
                    fval = 2.0  # to be refined based on atom types
                    ii = np.where(self.instance, ati)[0][0]
                    jj = np.where(self.instance, atj)[0][0]
                    self.set_bond(ii, jj, -999, change, step=step, stmax=self.max_step, findist=fval, geom=geom)

        self.clean_constraints(change, fix)

        return step, fix, change, release
