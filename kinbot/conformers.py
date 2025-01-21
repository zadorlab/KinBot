import os
import random
import time
import copy
import logging
from shutil import copyfile

import numpy as np
from ase.db import connect
from ase import Atoms
from ase.units import invcm, Hartree, kcal, mol
from ase.thermochemistry import IdealGasThermo
import rmsd

from kinbot import geometry
from kinbot import zmatrix
from kinbot.stationary_pt import StationaryPoint
from kinbot import constants

logger = logging.getLogger('KinBot')


class Conformers:
    """
    Class that does all the steps for the conformers of one species
    """
    def __init__(self, species, par, qc, semi_emp=0):
        """
        species: instance of StationaryPoint
        qc: instance of QuantumChemistry
        semi_emp: is this search at low level (e.g. am1) or at the L1 level. The latter is the default.
        """
        self.species = species
        self.qc = qc
        # status of the conformational analysis
        # -1: not yet started
        # 0: running
        # 1: finished
        # -999:failed
        self.scycconf = -1
        self.sconf = -1

        self.grid = par['conf_grid']
        self.diffthrs = par['difference_threshold']

        # final geometries of the cyclic search
        self.cyc_conf_geoms = []

        # for each conformer, check the index of its progress
        self.cyc_conf_index = []
        # dihedral atoms of the cyclic conformers
        self.cyc_dih_atoms = []
        # dihedral values of the cyclic conformers
        self.cyc_dih_values = []
        # number of cyclic conformers generated
        self.cyc_conf = 0

        # number of open chain conformers generated
        self.conf = 0

        # -1 (not finished), 0 (successful) or
        # 1 (failed) for each cyclic conformer
        self.cyc_conf_status = []
        # -1 (not finished), 0 (successful) or
        # 1 (failed) for each open chain conformer
        self.conf_status = []
        self.zf = par['zf']

        # do semi empirical conformer search?
        self.semi_emp = semi_emp

        # Maximum number of diherals for which exhaustive
        # conformation searches are done
        self.max_dihed = par['max_dihed']
        # Number of random conformers in case no
        # exhaustive search is done
        self.nconfs = par['random_conf']

        self.info = True
        
        if semi_emp:
            # Maximum number of diherals for which exhaustive
            # conformation searches are done
            self.max_dihed = par['max_dihed_semi_emp']
            # Number of random conformers in case no
            # exhaustive search is done
            self.nconfs = par['random_conf_semi_emp']

        # db to be used for skipping conf generation
        self.db = connect('kinbot.db')

        self.imagfreq_threshold = par['imagfreq_threshold']
        self.flat_ring_dih_angle = par['flat_ring_dih_angle']
        self.print_warning = True

    def generate_ring_conformers(self, cart):
        """
        Generate the conformers of a cyclic structure
        by randomly sampling the dihedrals of the ring
        """
        # iterate the different rings in the species
        ncyc = len(self.species.cycle_chain)
        for cyc in self.species.cycle_chain:
            if len(cyc) > 3:  # three membered rings don't have conformers
                dihs = []  # list of the ring dihedrals
                for i, at in enumerate(cyc):
                    dihs.append([cyc[i-3], cyc[i-2], cyc[i-1], cyc[i]])

                # define the flatness of the ring by the sum of the
                # absolute values of the dihedrals along the ring
                # divided by the number of atoms in the ring
                cycdih = 0.
                # list of flat sections of the ring
                flat_ring_dih = []  
                for dih in dihs:
                    val = geometry.calc_dihedral(cart[dih[0]], cart[dih[1]],
                                                 cart[dih[2]], cart[dih[3]])[0]
                    if abs(val) < self.flat_ring_dih_angle:
                        flat_ring_dih.append(True)
                    else:
                        flat_ring_dih.append(False)
                        cycdih += np.abs(val)
                # only care about flatness for non-flat parts
                n_dih_nonflat = len(cyc) - np.sum(flat_ring_dih)
                if n_dih_nonflat > 0:
                    cycdih /= float(n_dih_nonflat) 
                else:  # for instance benzene
                    self.cyc_conf_geoms.append(copy.deepcopy(cart))
                    continue

                # randomly select at most N-3 dihedrals,
                # with N the number of non-flat dihedrals in the ring
                # number of independent dihedrals
                if n_dih_nonflat < 4:
                    nd = 1
                else:
                    nd = n_dih_nonflat - 3
                random_dihs = list(random.sample(list(np.array(dihs)[[not f for f in flat_ring_dih]]), nd))
                
                # number of conformers (nc) per ring conformer:
                # 4, 5, 6 member rings nc = 3 ^ nd
                # 7+ member rings = 27
                if len(cyc) < 7:
                    nc = np.power(3, nd)
                else:
                    nc = 27  # 3 ^ (6-3)

                # new thing for structures with multiple rings (whether fused or not)
                nc = int(nc / ncyc)

                for i in range(nc):
                    self.cyc_dih_atoms.append(random_dihs)
                    # values the dihedrals will be modified to
                    values = []
                    for j in range(nd):
                        values.append(cycdih*(np.mod(i // np.power(3, j), 3) - 1))
                    self.cyc_dih_values.append(values)
                    self.cyc_conf_index.append(-1)
                    self.cyc_conf += 1
            else:
                self.cyc_conf_geoms.append(copy.deepcopy(cart))

        for ci in range(self.cyc_conf):
            self.start_ring_conformer_search(ci, copy.deepcopy(self.species.geom))

    def start_ring_conformer_search(self, index, cart):
        """
        index: number of the conformer
        In each iteration a given dihedral is changed, and then in the 
        next one it's fixed and another one is changed, and then two are
        fixed and the next one is changed, until all are at their desired values
        """
        if self.cyc_conf_index[index] == len(self.cyc_dih_atoms[index]) - 1:
            # this conformer has finished
            return 0
        else:
            self.cyc_conf_index[index] += 1
            fix = []
            change = []
            for j, da in enumerate(self.cyc_dih_atoms[index]):
                if j == self.cyc_conf_index[index]:
                    new_dih = self.cyc_dih_values[index][j]
                    change.append([da[0] + 1, da[1] + 1, da[2] + 1, da[3] + 1, new_dih])
                    break
                else:
                    fix.append([da[0] + 1, da[1] + 1, da[2] + 1, da[3] + 1])
            for i in range(self.species.natom - 1):
                for j in range(i+1, self.species.natom):
                    if self.species.bond[i][j] > 0:
                        fix.append([i + 1, j + 1])
            self.qc.qc_ring_conf(self.species, cart, fix, change, index, self.cyc_conf_index[index])
        return 1

    def test_ring_conformer(self, index):
        """
        Test whether a conformer has the same bond matrix as the original structure.
        Returns the conformer object and -1 if not yet finished, 0 if same, and 1 if not.
        """
        job = self.get_job_name(index, cyc=1)

        status, geom = self.qc.get_qc_geom(job, self.species.natom)
        if status == 1:  # still running
            return np.zeros((self.species.natom, 3)), -1
        elif status == -1:  # conformer search failed
            logger.debug('Conformer search failed for scan point {}'.format(job))
            return np.zeros((self.species.natom, 3)), 1
        else:
            if self.start_ring_conformer_search(index, geom):
                logger.debug('Running the next dihedral for conformer {}'.format(job))
                return geom, -1
            else:
                # check if all the bond lenghts are withing 10% or the original bond lengths
                dummy = StationaryPoint('dummy',
                                       self.species.charge,
                                       self.species.mult,
                                       atom=self.species.atom,
                                       geom=geom)
                dummy.bond_mx()
                dummy.calc_chemid()
                if geometry.equal_geom(self.species, dummy, 0.10):
                    logger.debug('Successfully finished conformer {}'.format(job))
                    return geom, 0
                else:
                    logger.debug('Conformer too far from original structure {}'.format(job))
                    return np.zeros((self.species.natom, 3)), 1

    def check_ring_conformers(self, wait=0):
        """
        Check if the conformer optimizations finished.
        Test them, and submit frequency calculations.
        Then select the lowest energy one.
        returns:
        *status: 0 if still running, 1 if finished
        *geometries of all the conformers
        wait: wait for all the conformer calculations to finish before returning anything
        """
        if len(self.cyc_conf_status) < self.cyc_conf:
            for i in range(self.cyc_conf):
                self.cyc_conf_status.append(-1)
        while 1:
            # check if conformational search is finished
            for i, si in enumerate(self.cyc_conf_status):
                if si == -1:
                    self.cyc_conf_status[i] = self.test_ring_conformer(i)[1]
            if all([si >= 0 for si in self.cyc_conf_status]):
                geoms = [self.species.geom]  # list used for intermediate ring conformer generation
                final_geoms = []
                for ci in range(self.cyc_conf):
                    si = self.cyc_conf_status[ci]
                    if si == 0:  # this is a valid confomer
                        job = self.get_job_name(ci, cyc=1)
                        err, geom = self.qc.get_qc_geom(job, self.species.natom)
                        geoms.append(geom)
                        final_geoms.append(geom)
                    else:
                        final_geoms.append(np.zeros((self.species.natom, 3)))
                self.write_profile(self.cyc_conf_status, final_geoms, [0 for gi in final_geoms], ring=1)
                return 1, geoms
            else:
                if wait:
                    time.sleep(1)
                else:
                    return 0, np.zeros((self.species.natom, 3))

    def generate_conformers(self, rotor, cart, print_warning=False):
        """
        Generate guesses for all of the canonical conformers.
        This is a recursive routine to generate them.
        rotor: the rotor number in the order it was discovered
        if -999, then just do a single calculation at the given geometry
        """
        if self.cyc_conf == 0:
            cycles = 1
        else:
            cycles = self.cyc_conf

        name = self.get_name()
 
        # what is the value of cycles
        # what is value of all things associated w/ conf generation
        # what is length of conf_dihed?
        theoretical_confs = np.power(self.grid, len(self.species.conf_dihed)) * cycles

        if rotor != -999:
            if len(self.species.conf_dihed) > self.max_dihed or theoretical_confs > self.nconfs:
                if rotor == 0:
                    if self.info: 
                        logger.info('\tRandom conformer search is carried out for {}.'.format(name))
                        self.info = False

                    # skipping generation if done
                    if self.cyc_conf > 1:
                        nrandconf = int(round(self.nconfs / self.cyc_conf) + 2)
                    else:
                        nrandconf = self.nconfs
                    if os.path.exists('{}.log'.format(self.get_job_name(nrandconf - 1))) and os.path.exists('conf/{}_low.log'.format(name)): 
                        rows = self.db.select(name=self.get_job_name(nrandconf - 1))
                        for row in rows:
                            self.conf = nrandconf
                            logger.debug('\tLast conformer was found in kinbot.db, generation is skipped for {}.'.format(self.get_job_name(nrandconf)))
                            return 1

                self.generate_conformers_random_sampling(cart)
                return 0

        # retraction from the recursion
        if rotor == len(self.species.conf_dihed) or rotor == -999:
            self.qc.qc_conf(self.species, cart, self.conf, semi_emp=self.semi_emp)
            if self.conf == 0:
                logger.debug('Theoretical number of conformers is {} for {}.'.format(theoretical_confs, name))
            self.conf += 1
            return 0

        # skipping generation if done
        if os.path.exists('{}.log'.format(self.get_job_name(theoretical_confs - 1))) and os.path.exists('conf/{}_low.log'.format(name)): 
            rows = self.db.select(name=self.get_job_name(theoretical_confs - 1))
            for row in rows:
                self.conf = theoretical_confs
                if print_warning:
                    logger.debug('Theoretical number of conformers is {} for {}.'.format(theoretical_confs, name))
                    logger.info('\tLast conformer was found in kinbot.db, generation is skipped for {}.'.format(name))
                return 1

        cart = np.asarray(cart)
        zmat_atom, zmat_ref, zmat, zmatorder = zmatrix.make_zmat_from_cart(self.species, rotor, cart, 1)

        rotor += 1

        for gr in range(self.grid):
            zmat[3][2] += 360. / self.grid
            for i in range(4, self.species.natom):
                if zmat_ref[i][2] == 4:
                    zmat[i][2] += 360. / self.grid
                if zmat_ref[i][2] == 1:
                    zmat[i][2] += 360. / self.grid
            cartmod = zmatrix.make_cart_from_zmat(zmat, zmat_atom, zmat_ref, self.species.natom, self.species.atom, zmatorder)
            self.generate_conformers(rotor, cartmod)

        return 0

    def generate_conformers_random_sampling(self, ini_cart):
        """
        Generate a random sampling of each dihedral for a number nconfs of conformers
        """
        if self.cyc_conf > 1:
            nrandconf = int(round(self.nconfs/self.cyc_conf) + 2)
        else:
            nrandconf = self.nconfs

        for ni in range(nrandconf):
            cart = copy.deepcopy(ini_cart)
            if ni == 0:
                sample = [0. for di in self.species.conf_dihed]
            else:
                sample = [random.choice([0., 120., 240.]) for di in self.species.conf_dihed]
            for rotor in range(len(self.species.conf_dihed)):
                zmat_atom, zmat_ref, zmat, zmatorder = zmatrix.make_zmat_from_cart(self.species, rotor, cart, 1)
                zmat[3][2] += sample[rotor]
                for i in range(4, self.species.natom):
                    if zmat_ref[i][2] == 4:
                        zmat[i][2] += sample[rotor]
                    if zmat_ref[i][2] == 1:
                        zmat[i][2] += sample[rotor]
                cart = zmatrix.make_cart_from_zmat(zmat, zmat_atom, zmat_ref, self.species.natom, self.species.atom, zmatorder)
            self.qc.qc_conf(self.species, cart, self.conf, semi_emp=self.semi_emp)
            self.conf += 1
        return 0

    def test_conformer(self, conf):
        """
        Test whether a conformer has the same bond matrix as the original structure.
        Returns the conformer object and -1 if not yet finished, 0 if same, and 1 if not.
        """
        add = ''
        if self.semi_emp:
            add = 'semi_emp_'
        job = self.get_job_name(conf, add=add)

        status, geom = self.qc.get_qc_geom(job, self.species.natom)
        if status == 1:  # still running
            return np.zeros((self.species.natom, 3)), -1
        elif status == -1:  # conformer search failed
            return np.zeros((self.species.natom, 3)), 1
        else:
            # check if all the bond lenghts are withing 10% of the original bond lengths
            dummy = StationaryPoint('dummy',
                                    self.species.charge,
                                    self.species.mult,
                                    atom=self.species.atom,
                                    geom=geom)
            dummy.bond_mx()
            dummy.calc_chemid()
            if geometry.equal_geom(self.species, dummy, 0.10):
                return geom, 0
            else:
                return np.zeros((self.species.natom, 3)), 1

    def check_conformers(self, wait=0):
        """
        Check if the conformer optimizations finished.
        Test connectivity and frequencies.
        Then select the lowest energy one (incl zpe).
        returns:
        *status: 0 if still running, 1 if finished
        *geometry of lowest energy conformer
        wait: wait for all the conformer calculations to finish before returning anything
        """
        if len(self.conf_status) < self.conf:
            for i in range(len(self.conf_status), self.conf):
                self.conf_status.append(-1)
        status = self.conf_status

        lowest_conf = str(0).zfill(self.zf)  # the index of the lowest conf, to be updated as we go

        while 1:
            # check if conformational search is finished
            name = self.get_name()
            for i, si in enumerate(status):
                if si == -1:
                    status[i] = self.test_conformer(i)[1]
            # problem: if the first sample has 2 imaginary frequencies, what to do then?
            if all([si >= 0 for si in status]):
                if isinstance(name, int):
                    lowest_job = f'{name}_well'
                else:
                    lowest_job = name
                *_, last_row = self.db.select(name=f'{lowest_job}')
                # The following refers to the conformer with lowest E + ZPE, 
                # not the individual lowest.
                lowest_energy = np.inf
                lowest_zpe = np.inf
                lowest_e_geom = last_row.positions
                final_geoms = []  # list of all final conformer geometries
                totenergies = []
                frequencies = []

                if all(status):  # if all conformers are invalid, 1 (different) or fail (-1)
                    if self.qc.qc == 'gauss':
                        ext = 'log'
                    elif self.qc.qc == 'qchem':
                        ext = 'out'
                    elif self.qc.qc == 'nn_pes':
                        ext = 'log'
                    else:
                        raise NotImplementedError(f'Code {self.qc.qc} not available.')
                    copyfile(f'{lowest_job}.{ext}', f'conf/{name}_low.{ext}')
                    mol = Atoms(symbols=last_row.symbols, positions=last_row.positions)
                    data = {'energy': last_row.data.get('energy'),
                            'frequencies': last_row.data.get('frequencies'),
                            'zpe': last_row.data.get('zpe'),
                            'status': last_row.data.get('status')}
                    self.db.write(mol, name='conf/{}_low'.format(name), 
                                  data=data)
                    #logger.warning(f'All conformer optimizations failed for {name}.')

                    return 1, lowest_conf, lowest_e_geom, last_row.data.get('energy'),\
                           final_geoms, totenergies, frequencies, status

                for ci in range(self.conf):
                    if status[ci] == 0:  # this is a valid confomer
                        add = ''
                        if self.semi_emp:
                            add = 'semi_emp_'
                        job = self.get_job_name(ci, add=add)
                        err, energy = self.qc.get_qc_energy(job)
                        err, zpe = self.qc.get_qc_zpe(job)
                        err, geom = self.qc.get_qc_geom(job, self.species.natom)
                        err, freq = self.qc.get_qc_freq(job, self.species.natom)
                        final_geoms.append(geom)
                        totenergies.append(energy + zpe)
                        if freq != []:
                            frequencies.append(freq)
                        else:
                            frequencies.append(None)
                        if lowest_energy is np.inf:
                            if self.species.natom > 1:
                                # job fails if conformer freq array is empty
                                if len(freq) > 0:
                                    if self.species.wellorts:
                                        if freq[0] >= 0.:
                                            err = -1
                                        if self.species.natom > 2 and freq[1] <= -1 * self.imagfreq_threshold:
                                            err = -1
                                    else:
                                        if freq[0] <= -1 * self.imagfreq_threshold:
                                            err = -1
                                else:
                                    logger.warning("Conformer {} failed due to empty freq array".format(ci))
                                    err = -1
                            if err == 0:
                                lowest_energy = energy
                                lowest_zpe = zpe
                        if energy + zpe <= lowest_energy + lowest_zpe:
                            err, freq = self.qc.get_qc_freq(job, self.species.natom)
                            ratio = 0.8
                            # job fails if conformers freq array is empty
                            if len(freq) > 0:
                                if self.species.wellorts:
                                    if freq[0] / self.species.freq[0] < ratio:
                                        err = -1 
                                    if freq[0] / self.species.freq[0] > 1. / ratio:
                                        err = -1 
                                    if self.species.natom > 2 and freq[1] <= 0.:
                                        err = -1
                                else:
                                    if freq[0] <= -1. * self.imagfreq_threshold:  # note that now we allow negative frequencies here as well
                                        err = -1
                            elif self.species.natom > 1:
                                logger.warning("Conformer {} failed due to empty freq array".format(ci))
                                err = -1
                            if err == 0:
                                lowest_job = job
                                lowest_conf = str(ci).zfill(self.zf)
                                lowest_energy = energy
                                lowest_zpe = zpe
                                lowest_e_geom = geom
                        if err == -1:
                            status[ci] = 1  # make it invalid
                    else:
                        totenergies.append(0.)
                        final_geoms.append(np.zeros((self.species.natom, 3)))
                        if self.species.natom == 1:
                            frequencies.append(None)
                        elif self.species.natom == 2:
                            frequencies.append(np.zeros(self.species.natom * 3 - 5)) 
                        else:
                            frequencies.append(np.zeros(self.species.natom * 3 - 6))

                self.write_profile(status, final_geoms, totenergies)

                # Check if at least one conformer has the same enrgy as the L1 parent structure
                if self.species.wellorts:
                    *_, l1_last_row = self.db.select(name=self.species.name)
                else:
                    *_, l1_last_row = self.db.select(name=f'{self.species.name}_well')
                l1energy = l1_last_row.data.get('energy') * constants.EVtoHARTREE
                l1energy += l1_last_row.data.get('zpe')
                if not any([abs(en - l1energy) < self.diffthrs * constants.KCALtoHARTREE for en in totenergies]): # 0.1 kcal/mol
                    if self.print_warning:
                        logger.warning(f'\tNone of {self.species.name} '
                                       'conformers has the same energy as its '
                                       'parent structure.')
                        self.print_warning = False
                    lowest_job = l1_last_row.name
                    lowest_conf = 'low'
                    lowest_e_geom = l1_last_row.positions
                    lowest_energy = l1energy
                
                low_row = None
                low_rows = self.db.select(name='conf/{}_low'.format(name))
                for lrow in low_rows:
                    low_row = lrow
                try:
                    if self.qc.qc == 'gauss':
                        copyfile(f'{lowest_job}.log', f'conf/{name}_low.log')
                    elif self.qc.qc == 'qchem':
                        copyfile(f'{lowest_job}.out', f'conf/{name}_low.out')
                    elif self.qc.qc == 'nn_pes':
                        pass
                    else:
                        raise NotImplementedError(f'Code {self.qc.qc} not available.')
                    rows = self.db.select(name='{}'.format(lowest_job))
                    for row in rows:
                        row_last = row
                    mol = Atoms(symbols=row_last.symbols, positions=row_last.positions)
                    data = {'energy': row_last.data.get('energy'),
                            'frequencies': row_last.data.get('frequencies'),
                            'zpe': row_last.data.get('zpe'),
                            'status': row_last.data.get('status')}
                    if low_row and hasattr(low_row, 'data') \
                            and all((np.all(data.get(k) == v) 
                                     for k, v in low_row.data.items())):
                        pass
                    else:
                        self.db.write(mol, name='conf/{}_low'.format(name), 
                                      data=data)
                except UnboundLocalError:
                    pass

                return 1, lowest_conf, lowest_e_geom, lowest_energy,\
                       final_geoms, totenergies, frequencies, status

            else:
                if wait:
                    time.sleep(1)
                else:
                    return 0, lowest_conf, np.zeros((self.species.natom, 3)), \
                           self.species.energy, np.zeros((self.species.natom, 3)), \
                           np.zeros(1), np.zeros(1), np.zeros(1)


    def add_new_conf_from_hir(self, geom) -> None:
        """Generate new conformer found during hir calculations.

        Args:
            geom (_type_): cartesian geometry of the lowest energy point along all hir
        """
        self.qc.qc_conf(self.species, geom, self.conf, semi_emp=self.semi_emp)
        if os.path.exists(f'conf/{self.get_name()}_low.log'):
            logger.debug(f'Removing file conf/{self.get_name()}_low.log')
            os.remove(f'conf/{self.get_name()}_low.log')
        self.conf += 1


    def lowest_conf_info(self):
        """
        in case conformer search was successfully skipped on restart, 
        this reads the energy and geometry of the lowest energy conf directly
        """
        name = self.get_name()

        job = 'conf/{}_low'.format(name)
        try:
            err, energy = self.qc.get_qc_energy(job)
            err, zpe = self.qc.get_qc_zpe(job)
            err, geom = self.qc.get_qc_geom(job, self.species.natom)
        except ValueError:
            _ = self.check_conformers()
            err, energy = self.qc.get_qc_energy(job)
            err, zpe = self.qc.get_qc_zpe(job)
            err, geom = self.qc.get_qc_geom(job, self.species.natom)
                
        return geom, energy, zpe 

    def find_unique(self, conformers, energies, frequencies, valid, temp=None, boltz=None):
        """
        Given a set of conformers, finds the set of unique ones.
        Algorithm:
        If valid:
            If energy is close to test energy:
                If moments of inertia are close to test moi:
                    If rmsd is small:
                        They are the same
        Otherwise unique

        test all previous structures.

        temp is temperature, and only exp(-G/RT) > boltz conformers are considered if defined.
        returns the geometries, total energies, frequencies, and indices (as in the /conf directory)
        """

        conformers_unq = []
        energies_unq = []
        zeroenergies_unq = []
        frequencies_unq = []
        indices_unq = []

        if temp is not None:
            # calculate the Gibbs free energy for all conformers
            # at T = temp, P = 101325 Pa
            gibbs = []
            geo_type = 'nonlinear'
            for vi, val in enumerate(valid):
                if frequencies == [None]:
                    vib_energies = [0]
                    geo_type = 'monatomic'
                else:
                    vib_energies = [ff * invcm for ff in frequencies[vi] if ff > 0]  # convert to eV
                    if np.shape(frequencies)[1] == 3 * len(self.species.atom) - 5:
                        geo_type = 'linear'
                potentialenergy = energies[vi] * Hartree  # convert to eV
                atoms = Atoms(symbols=self.species.atom, positions=conformers[vi])
                thermo = IdealGasThermo(vib_energies=vib_energies,
                                        potentialenergy=potentialenergy,
                                        atoms=atoms,
                                        geometry=geo_type,
                                        symmetrynumber=1, spin=(self.species.mult-1)/2)
                gibbs.append(thermo.get_gibbs_energy(temperature=temp, pressure=101325., verbose=False))

        for vi, val in enumerate(valid):
            unique = True
            if val == 0:
                for ei, en in enumerate(energies_unq):
                    if abs(energies[vi] - en) * constants.AUtoKCAL < 0.2:
                        moi_test, _ = geometry.get_moments_of_inertia(conformers[vi], self.species.atom)
                        moi_unq, _ = geometry.get_moments_of_inertia(conformers_unq[ei], self.species.atom)
                        if all(moi_test / moi_unq) < 1.1 and all(moi_test / moi_unq) > 0.9:
                            p_coord = copy.deepcopy(conformers[vi])
                            q_coord = copy.deepcopy(conformers_unq[ei])
                            p_atoms = copy.deepcopy(self.species.atom)
                            q_atoms = copy.deepcopy(self.species.atom)
                            p_cent = rmsd.centroid(p_coord)
                            q_cent = rmsd.centroid(q_coord)
                            p_coord -= p_cent
                            q_coord -= q_cent
                            rotation_method = rmsd.kabsch_rmsd
                            reorder_method = rmsd.reorder_hungarian
                            #q_review = reorder_method(p_atoms, q_atoms, p_coord, q_coord)
                            #q_coord = q_coord[q_review]
                            #q_atoms = q_atoms[q_review]
                            result_rmsd = rotation_method(p_coord, q_coord)
                            if result_rmsd < 0.05:
                                unique = False
                                break
                if unique:
                    if temp is None:
                        conformers_unq.append(conformers[vi])
                        energies_unq.append(energies[vi])
                        frequencies_unq.append(frequencies[vi])
                        indices_unq.append(vi)
                    else:
                        if np.exp(-1000. * (gibbs[vi] - min(gibbs)) / (kcal / mol) /\
                                  (constants.R / constants.CALtoJ * temp)) > boltz:
                            conformers_unq.append(conformers[vi])
                            energies_unq.append(energies[vi])
                            frequencies_unq.append(frequencies[vi])
                            indices_unq.append(vi)

        zeroenergies_unq = [0.] * len(energies_unq)
        return conformers_unq, energies_unq, zeroenergies_unq, frequencies_unq, indices_unq

    def write_profile(self, status, final_geoms, energies, ring=0):
        """
        Write a molden-readable file with the CONF analysis (geometries and total energies)
        """
        r = ''
        if ring:
            r = 'r'
        if self.species.wellorts:
            ff = open('conf/' + self.species.name + r + '.xyz', 'w')
        else:
            ff = open('conf/' + str(self.species.chemid) + r + '.xyz', 'w')
        for i, st in enumerate(status):
            s = str(self.species.natom) + '\n'
            s += 'energy = ' + str(energies[i]) + '\n'
            for j, at in enumerate(self.species.atom):
                x, y, z = final_geoms[i][j]
                s += '{} {:.8f} {:.8f} {:.8f}\n'.format(at, x, y, z)
            if st == 0:  # valid conformer:
                ff.write(s)
        ff.close()

    def get_name(self):
        if self.species.wellorts:
            name = self.species.name
        else:
            name = self.species.chemid
        return name

    def get_job_name(self, idx, cyc=0, add=''):
        name = str(self.get_name())

        if cyc == 0:
            job = 'conf/' + name + '_' + add + str(idx).zfill(self.zf)
        else:
            job = 'conf/' + name + '_r' + str(idx).zfill(self.zf) + '_' + str(self.cyc_conf_index[idx]).zfill(self.zf)
        return job
