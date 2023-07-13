from ase.data import atomic_numbers, covalent_radii

from kinbot import geometry
from kinbot.reac_General import GeneralReac
from kinbot import constants


class IntraRAddExoTetCyclicF(GeneralReac):
    max_step = 22
    scan = 0
    skip = 1
    dihstep = 12
    family_name = 'intraraddexotetcyclicf'
    

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)
        if step < self.dihstep:
            self.set_dihedrals(change, step, cut=1)

            ldih = [] # constraint for the last dihedral, which needs to be 180 degrees    
            for i in range(4):
                ldih.append(self.instance[len(self.instance)-4+i] + 1)
            dih = geometry.calc_dihedral(geom[ldih[0] - 1], geom[ldih[1] - 1], geom[ldih[2] - 1], geom[ldih[3] - 1])[0]
            frac = 1./(12. - step)
            if dih < 0:
                new_dih = dih - frac * (180. + dih) 
                ldih.append(new_dih)
            else:
                new_dih = dih + frac * (180. - dih) 
                ldih.append(new_dih)
            change.append(ldih)

        elif step < self.max_step:
            self.release_dihedrals(release)
            try:
                fdist1 = constants.st_bond[''.join(sorted(self.species.atom[self.instance[0]] + self.species.atom[self.instance[-2]]))]
            except KeyError:
                fdist1 = 1.2 * (covalent_radii[atomic_numbers[self.species.atom[self.instance[0]]]] 
                                + covalent_radii[atomic_numbers[self.species.atom[self.instance[-2]]]])
            if ''.join(sorted(self.species.atom[self.instance[0]] + self.species.atom[self.instance[-2]])) == 'CO':
                if ''.join(sorted(self.species.atom[self.instance[-1]] + self.species.atom[self.instance[-2]])) == 'OO':
                    fdist1 = 1.96
                else:
                    fdist1 = 1.68
            self.set_bond(0, -2, -999, change, step=step-11, stmax=10, findist=fdist1, geom=geom)
            try:
                fdist2 = constants.st_bond[''.join(sorted(self.species.atom[self.instance[-1]] + self.species.atom[self.instance[-2]]))]
            except KeyError:
                fdist2 = 1.2 * (covalent_radii[atomic_numbers[self.species.atom[self.instance[-1]]]] 
                                + covalent_radii[atomic_numbers[self.species.atom[self.instance[-2]]]])
            
            if ''.join(sorted(self.species.atom[self.instance[-1]] + self.species.atom[self.instance[-2]])) == 'CO':
                fdist2 = 1.68
            self.set_bond(-1, -2, -999, change, step=step-11, stmax=10, findist=fdist2, geom=geom)
        self.clean_constraints(change, fix)
        
        return step, fix, change, release
