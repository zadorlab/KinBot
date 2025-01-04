import os
import copy
import logging
import time

import rmsd
import numpy as np

from kinbot import frequencies
from kinbot import geometry
from kinbot import symmetry
from kinbot.conformers import Conformers
from kinbot.hindered_rotors import HIR
from kinbot.molpro import Molpro
from kinbot.orca import Orca
from kinbot import reader_gauss, reader_qchem
from kinbot.stationary_pt import StationaryPoint
from kinbot import constants

logger = logging.getLogger('KinBot')


class Optimize:
    """
    This class does the following:

    1. Conformational search of the species
    2. High level optimization and freq calc of the species
    3. Hindered rotor scans
    4. Repeat steps 2-3 as long as lower energy structures are found
    """

    def __init__(self, species, par, qc, wait=0, just_high=False):
        """
        just_high: ony do high_level calculation
        """
        self.species = species

        try:
            delattr(self.species, 'cycle_chain')
        except AttributeError:
            logger.debug(f"{self.species.chemid} has no cycle_chain attribute "
                         f"to delete")
        if self.species.wellorts or par['cluster']:
            self.species.characterize(bond_mx=self.species.bond)
            self.name = str(self.species.name)
        else:
            self.species.characterize()
            self.name = self.species.name
        self.par = par
        self.qc = qc
        # wait for all calculations to finish before returning
        self.wait = wait

        self.shigh = -1
        self.just_high = just_high

        if not self.just_high:
            # status of the various parts
            # -1: not yet started
            #  0: running
            #  1: finished
            # -999:failed
            self.scycconf = -1
            self.sconf = -1
            self.shir = -1
        else:
            # Do not perform conformer search for vdW wells
            self.scycconf = 1
            self.sconf = 1
            self.shir = 1

        # restart counter: number of times the high-level and hir calculations
        # has been restarted in case a lower energy structure has been found
        self.restart = 0

        self.skip_conf_check = 0  # initialize

    def do_optimization(self):
        while 1:
            # do the conformational search
            if self.par['conformer_search'] == 1:
                if self.scycconf == -1 and self.sconf == -1:
                    logger.info('\tStarting conformational search of '
                                f'{self.name}')
                    self.species.confs = Conformers(self.species, self.par, self.qc)

                # first do the cyclic part of the molecule
                if self.scycconf == -1:
                    # start the ring conf search
                    if len(self.species.cycle_chain) > 0:
                        # there are rings in the molecule, do a search
                        self.species.confs.generate_ring_conformers(copy.deepcopy(self.species.geom))
                        # set the cyclic conf status to running
                        self.scycconf = 0
                    else:
                        # there are no rings in the molecule, continue from the current geometry
                        self.species.confs.cyc_conf_geoms = [copy.deepcopy(self.species.geom)]
                        # no ring conf search has to be done, set status to finished
                        self.scycconf = 1
                if self.scycconf == 0:
                    # ring conf search is running, check if finished
                    status, self.species.confs.cyc_conf_geoms = self.species.confs.check_ring_conformers()
                    if status:
                        # ring conf search is finished
                        self.scycconf = 1
                # first do an semi empirical optimization if requested by the user
                if self.par['semi_emp_conformer_search'] == 1:
                    logger.info('\tSemi-empirical conformer search is starting '
                                f'for {self.name}')
                    if self.ssemi_empconf == -1:
                        # semi empirical part has not started yet
                        self.species.semi_emp_confs = Conformers(self.species, self.par, self.qc, semi_emp=1)
                        for geom in self.species.confs.cyc_conf_geoms:
                            # take all the geometries from the cyclic part
                            # generate the conformers for the current geometry
                            self.species.semi_emp_confs.generate_conformers(0, geom)
                        # set conf status to running
                        self.ssemi_empconf = 0
                        if self.ssemi_empconf == 0:
                            # semi empirical conformational search is running
                            # check if the conformational search is done
                            status, lowest_conf, geom, self.semi_emp_low_energy, self.semi_emp_conformers, self.semi_emp_energies = self.species.semi_emp_confs.check_conformers(wait=self.wait)
                            if status == 1:
                                logger.info("\tSemi- empirical lowest energy "
                                            f"conformer for species {self.name}"
                                            f" is number {lowest_conf}")
                                # set conf status to finished
                                self.ssemi_empconf = 1
                else:
                    self.ssemi_empconf = 1
                if self.ssemi_empconf == 1 and self.scycconf == 1:
                    # do open chain part if cyclic (if there were any) and semi empirical (if requested) parts are done
                    if self.sconf == -1:
                        # open chain part has not started yet
                        # if semi empirical conformer were searched for, start from those,
                        # else start from cyclic conformers
                        if self.par['semi_emp_conformer_search'] == 1:
                            self.species.confs.nconfs = 1
                            for i, geom in enumerate(self.semi_emp_conformers):
                                if (self.semi_emp_energies[i] - self.semi_emp_low_energy) * constants.AUtoKCAL < self.par['semi_emp_confomer_threshold']:
                                    self.species.confs.generate_conformers(-999, geom)
                            logger.info("\tThere are {} structures below the {} kcal/mol threshold for species {} in the semiempirical search.". \
                                         format(i, self.par['semi_emp_confomer_threshold'], self.name))
                        else:
                            print_warning = True
                            for geom in self.species.confs.cyc_conf_geoms:
                                # take all the geometries from the cyclic part
                                # generate the conformers for the current geometry
                                self.skip_conf_check = self.species.confs.generate_conformers(0, geom, print_warning=print_warning)
                                print_warning = False
                        # set conf status to running
                        self.sconf = 0
                    if self.sconf == 0:
                        # conformational search is running
                        # check if the conformational search is done
                        if self.skip_conf_check == 0 \
                                or self.par['multi_conf_tst'] \
                                or self.par['print_conf'] \
                                or self.par['calc_aie']:
                            status, lowest_conf, geom, self.species.low_energy, conformers, energies, frequency_vals, valid =\
                                    self.species.confs.check_conformers(wait=self.wait)
                            if status == 1:
                                self.species.conformer_geom, self.species.conformer_energy, \
                                        self.species.conformer_zeroenergy, \
                                        self.species.conformer_freq, self.species.conformer_index = \
                                        self.species.confs.find_unique(conformers, 
                                        energies, 
                                        frequency_vals, 
                                        valid,
                                        self.par['multi_conf_tst_temp'],
                                        self.par['multi_conf_tst_boltz'])
                                logger.info(f'\tLowest energy conformer for species {self.name} is number {lowest_conf}')
                                if self.par['multi_conf_tst'] or self.par['print_conf'] or self.par['calc_aie']:
                                    logger.info(f'\tAt {self.par["multi_conf_tst_temp"]} K '
                                                f'with {100*self.par["multi_conf_tst_boltz"]}% '
                                                f'cutoff the unique conformers for species {self.name} '
                                                f'are {self.species.conformer_index}')
                                # save lowest energy conformer as species geometry
                                self.species.geom = geom
                                # save lowest energy conformer energy
                                self.species.energy = self.species.low_energy
                                # set conf status to finished
                                self.sconf = 1
                        elif self.skip_conf_check == 1 and self.restart == 0:
                            self.species.geom, self.species.energy, self.species.zpe = self.species.confs.lowest_conf_info()
                            logger.info('\tEnergy and geometry updated based on conf/{}_low file.'.format(self.name))
                            self.sconf = 1

            else:
                # no conf search necessary, set status to finished
                self.sconf = 1
            if self.sconf == 1:  # conf search is finished
                # if the conformers were already done in a previous run
                if self.par['conformer_search'] == 1 and not self.just_high and self.restart == 0:
                    status, lowest_conf, self.species.geom, self.species.low_energy, conformers, energies, frequency_vals, valid = \
                        self.species.confs.check_conformers(wait=self.wait)
                        
                while self.restart <= self.par['rotation_restart']:
                    # do the high level calculations
                    if self.par['high_level'] == 1:
                        if self.shigh == -1:
                            if self.species.wellorts:
                                name = self.species.name
                                self.qc.qc_opt_ts(self.species, self.species.geom, high_level=1)
                                if self.par['multi_conf_tst']:
                                    for ci, conindx in enumerate(self.species.conformer_index):
                                        self.qc.qc_opt_ts(self.species, 
                                                          self.species.conformer_geom[ci], 
                                                          high_level=1,
                                                          ext=f'_{str(conindx).zfill(4)}_high',
                                                          )
                            else:
                                if self.just_high:
                                    name = self.species.name
                                    self.qc.qc_opt(self.species, self.species.geom, high_level=1, do_vdW=True)
                                    #Set conformer_zeroenergy to the species zero-energy if no conformer analysis has been done
                                    self.species.conformer_zeroenergy = [self.species.energy + self.qc.get_qc_zpe(f"{self.species.name}_high", wait=0)[1]]
                                else:
                                    name = self.species.chemid
                                    self.qc.qc_opt(self.species, self.species.geom, high_level=1)
                                if self.par['multi_conf_tst']:
                                    for ci, conindx in enumerate(self.species.conformer_index):
                                        self.qc.qc_opt(self.species, 
                                                       self.species.conformer_geom[ci], 
                                                       high_level=1,
                                                       ext=f'_{str(conindx).zfill(4)}_high',
                                                       )
                                if self.par['calc_aie']:
                                     for ci, conindx in enumerate(self.species.conformer_index):
                                        self.qc.qc_aie(self.species, 
                                                       self.species.conformer_geom[ci], 
                                                       ext=f'{str(conindx).zfill(4)}',
                                                       )
                                     logger.info('\tStarting AIE optimization(s) of {}'.format(name))
                            logger.info('\tStarting high level optimization(s) of {}'.format(name))
                            self.shigh = 0  # set the high status to running
                        if self.shigh == 0:
                            # high level calculation is running
                            # check if finished
                            status = self.qc.check_qc(self.log_name(1))
                            if status == 'error':
                                # found an error
                                logger.warning('High level optimization failed for {}'.format(self.name))
                                self.shigh = -999
                            elif status == 'normal':
                                self.compare_structures()  # this switches shigh to 0.5 or 1 and updates the geometry
                        if self.shigh == 0.5:  # the top one was tested already and was ok
                            stati = [0] * len(self.species.conformer_index)
                            for ci, conindx in enumerate(self.species.conformer_index):
                                status = self.qc.check_qc(self.log_name(1, conf=conindx))
                                if status == 'error':
                                    stati[ci] = 1
                                    self.species.conformer_index[ci] = -999
                                elif status == 'normal':
                                    stati[ci] = 1
                                    self.compare_structures(conf=conindx)
                            if sum(stati) == len(self.species.conformer_index):
                                self.shigh = 1
                    else:
                        # no high-level calculations necessary, set status to finished
                        self.shigh = 1
                    if self.shigh == 1:
                        # do the HIR calculation
                        if self.par['rotor_scan'] == 1 and self.par['multi_conf_tst'] != 1:
                            if self.shir == -1:
                                # hir not stated yet
                                logger.info('\tStarting hindered rotor calculations of {}'.format(self.name))
                                self.species.hir = HIR(self.species, self.qc, self.par)
                                self.species.hir.generate_hir_geoms(copy.deepcopy(self.species.geom), self.par['rigid_hir'])
                                self.shir = 0
                            if self.shir == 0:
                                # hir is running
                                # check if it is done:
                                status = self.species.hir.check_hir(wait=self.wait)
                                if status:  # Finished correctly
                                    if len(self.species.hir.hir_energies) > 0:
                                        # check if along the hir potential a structure was found with a lower energy
                                        min_en = self.species.hir.hir_energies[0][0]
                                        min_rotor = -1
                                        min_ai = -1
                                        for rotor in range(len(self.species.dihed)):
                                            for ai in range(self.species.hir.nrotation):
                                                # use a 0.1 kcal/mol cutoff for numerical noise
                                                if self.species.hir.hir_status[rotor][ai] == 0:  # do not test for fails
                                                    if self.species.hir.hir_energies[rotor][ai] < min_en - 1.6E-4:
                                                        min_en = self.species.hir.hir_energies[rotor][ai]
                                                        min_rotor = rotor
                                                        min_ai = ai
                                        if min_rotor > -1:
                                            self.restart += 1
                                            job = self.log_name(1, hir=1, r=min_rotor, s=min_ai)
                                            diff = (self.species.hir.hir_energies[min_rotor][0] - 
                                                    self.species.hir.hir_energies[min_rotor][min_ai]) * constants.AUtoKCAL 
                                            if self.restart < self.par['rotation_restart']:
                                                # lower energy structure found
                                                logger.warning(f'Lower energy conformer during HIR: {job}, diff = {np.round(diff, 2)} kcal/mol. Restart #{self.restart}')

                                                err, self.species.geom = self.qc.get_qc_geom(job, self.species.natom)
                                                # err, geom = self.qc.get_qc_geom(job, self.species.natom)
                                                # self.species.confs.add_new_conf_from_hir(geom)
                                                # delete the high_level log file and the hir log files
                                                if os.path.exists(self.log_name(1) + '.log'):
                                                    logger.debug(f'Removing file {self.log_name(1)}.log')
                                                    os.remove(self.log_name(1) + '.log')
                                                for rotor in range(len(self.species.dihed)):
                                                    for ai in range(self.species.hir.nrotation):
                                                        if os.path.exists(self.log_name(1, hir=1, r=rotor, s=ai) + '.log'):
                                                            logger.debug('Removing file ' + self.log_name(1, hir=1, r=rotor, s=ai) + '.log')
                                                            os.remove(self.log_name(1, hir=1, r=rotor, s=ai)  + '.log')
                                                # set the status of high and hir back to not started
                                                self.shigh = -1
                                                self.shir = -1
                                            else:
                                                logger.warning(f'Lower energy during HIR: {job}, diff = {np.round(diff, 2)} kcal/mol. Reached max restart for this species.')
                                                self.shir = 1
                                        else:
                                            self.shir = 1
                                    else:
                                        self.shir = 1
                        else:
                            # no hir calculations necessary, set status to finished
                            self.shir = 1
                            if self.par['multi_conf_tst'] == 1:
                                logger.debug('No rotor scans performed because multi conformer TST is requested.')
                    if not self.wait or self.shir == 1 or self.shigh == -999:
                        # break the loop if no waiting is required or
                        # if the hir calcs are done or
                        # if the high level calc failed
                        break
                    else:
                        time.sleep(1)
            if self.shir == 1:
                # finalize if no waiting is required or if the hir calcs are done
                # calculate the symmetry numbers
                symmetry.calculate_symmetry(self.species)

                # calculate the new frequencies with the internal rotations projected out
                if self.par['multi_conf_tst'] == 0 and self.par['rotor_scan'] \
                        and not self.just_high:
                    fr_file = self.log_name(self.par['high_level'])
                    hess = self.qc.read_qc_hess(fr_file, self.species.natom)
                    if self.qc.qc == 'qchem':
                        massweighted = True
                    else:
                        massweighted = False
                    self.species.kinbot_freqs, self.species.reduced_freqs = frequencies.get_frequencies(self.species, 
                                                                                                        hess, 
                                                                                                        self.species.geom, 
                                                                                                        massweighted=massweighted)
                else:
                    self.species.kinbot_freqs = self.species.freq
                    self.species.reduced_freqs = self.species.freq

                # write the L3 input and read the L3 energy, if available
                if self.par['L3_calc'] == 1:
                    if self.par['single_point_qc'] == 'molpro':
                        molp = Molpro(self.species, self.par)
                        if 'barrierless_saddle' in self.name:
                            key = self.par['barrierless_saddle_single_point_key']
                            molp.create_molpro_input(bls=1)
                        else:
                            key = self.par['single_point_key']
                            molp.create_molpro_input(do_vdW=self.just_high)
                        if self.par['queuing'] != 'local':
                            molp.create_molpro_submit(do_vdW=self.just_high)
                        status, molpro_energy = molp.get_molpro_energy(key, do_vdW=self.just_high)
                        if status:
                            self.species.l3energy = molpro_energy

                    elif self.par['single_point_qc'] == 'orca':
                        orca = Orca(self.species, self.par)
                        key = self.par['single_point_key']
                        orca.create_orca_input()
                        orca.create_orca_submit()
                        status, orca_energy = orca.get_orca_energy(key)
                        if status:
                            self.species.l3energy = orca_energy

                    self.delete_files()

            if self.wait:
                if self.shir == 1 or self.shigh == -999:
                    return 0
                time.sleep(1)
            else:
                return 0
            

    def delete_files(self):
        if self.par['delete_intermediate_files'] == 0:
            return 0

        extensions = ['chk', 'py', 'sbatch', 'ase', 'pbs', 'com']

        # job names
        names = []
        zf = self.par['zf']
        if self.species.wellorts:
            names.append(self.name)
            names.append(self.name + '_high')
            names.append(self.name + '_IRC_F')
            names.append(self.name + '_IRC_R')
            names.append(self.name + '_IRC_F_prod')
            names.append(self.name + '_IRC_R_prod')

            if self.par['high_level'] == 1:
                for count in range(self.species.hir.nrotation):
                    for rot_num in range(self.par['nrotation']):
                        names.append('hir/' + self.name + '_hir_' + str(count) + '_' + str(rot_num).zfill(2))
            if self.par['conformer_search'] == 1:
                for count in range(self.species.confs.conf):
                    names.append('conf/' + self.name + '_' + str(count).zfill(zf))
                for count in range(self.species.confs.cyc_conf):
                    for num in range(self.species.confs.cyc_conf_index[count]):
                        names.append('conf/' + self.name + '_r' + str(count).zfill(zf) + '_' + str(num).zfill(zf))
        else:
            names.append(self.name + '_well')
            names.append(self.name + '_well_high')
            if self.par['high_level'] == 1 and \
               'prod' not in self.species.name:
                for count in range(self.species.hir.nrotation):
                    for rot_num in range(self.par['nrotation']):
                        names.append('hir/' + self.name + '_hir_' + str(count) + '_' + str(rot_num).zfill(2))
            if self.par['conformer_search'] == 1 and \
               'prod' not in self.species.name:
                for count in range(self.species.confs.conf):
                    names.append('conf/' + self.name + '_' + str(count).zfill(zf))
                for count in range(self.species.confs.cyc_conf):
                    for num in range(self.species.confs.cyc_conf_index[count]):
                        names.append('conf/' + self.name + '_r' + str(count).zfill(zf) + '_' + str(num).zfill(zf))

        for name in names:
            for ext in extensions:
                # delete file
                file = '.'.join([name, ext])
                try:
                    os.remove(file)
                # except FileNotFoundError:
                except:
                    pass
        return 0


    def log_name(self, high, conf=-1, hir=0, r=None, s=None):
        """
        Provide base file names.
        """
        if hir == 1:
            return f'hir/{self.name}_hir_{str(r)}_{str(s).zfill(2)}'

        if conf >= 0 and high:
            return f'{self.name}_{str(conf).zfill(4)}_high'  # running in the main dir
        if conf >= 0:
            return f'conf/{self.name}_{str(conf).zfill(4)}'

        name = str(self.name)
        if not self.species.wellorts and not self.just_high:
            name += '_well'
        if high:
            name += '_high'

        return(name)


    def compare_structures(self, conf=-1):
        """
        Function to compare L1 and L2 strctures and decide is high is successful or not.
        If conf is >= 0, then we are testing for conformer number conf in the conf/ directory.
        """

        # creating a species for the L2
        err, new_geom = self.qc.get_qc_geom(self.log_name(1, conf=conf), self.species.natom, wait=self.wait)
        dummy = StationaryPoint('dummy',
                                self.species.charge,
                                self.species.mult,
                                atom=self.species.atom,
                                geom=new_geom)
        dummy.bond_mx()
        dummy.calc_chemid()

        # comparing L1 and L2 geometries and imaginary mode if TS
        if self.species.wellorts:  # for TS we need reasonable geometry agreement and normal mode correlation
            if self.par['conformer_search'] == 0:
                l1_file = self.log_name(0)  # name of the original L1 TS file
            elif self.par['multi_conf_tst'] or self.skip_conf_check == 0:
                l1_file = self.log_name(0, conf=conf)
            else:
                l1_file = 'conf/{}_low'.format(self.log_name(0))
            l2_file = self.log_name(1, conf=conf)
            if self.qc.qc == 'gauss':
                imagmode = reader_gauss.read_imag_mode(l1_file, self.species.natom)
                imagmode_high = reader_gauss.read_imag_mode(l2_file, self.species.natom)
            elif self.qc.qc == 'qchem':
                imagmode = reader_qchem.read_imag_mode(l1_file, self.species.natom)
                imagmode_high = reader_qchem.read_imag_mode(l2_file, self.species.natom)
            # either geom is roughly same with closely matching imaginary modes, or geometry is very close
            # maybe we need to do IRC at the high level as well...
            same_geom = ((geometry.matrix_corr(imagmode, imagmode_high) > 0.9) and \
                    (geometry.equal_geom(self.species, dummy, 0.3))) \
                    or (geometry.equal_geom(self.species, dummy, 0.2))
            if self.par['multi_conf_tst'] != 1:  # for now skipping this
                p_coord = copy.deepcopy(self.species.geom)
                q_coord = copy.deepcopy(dummy.geom)
                p_atoms = self.species.atom
                q_atoms = self.species.atom
                p_cent = rmsd.centroid(p_coord)
                q_cent = rmsd.centroid(q_coord)
                p_coord -= p_cent
                q_coord -= q_cent
                rotation_method = rmsd.kabsch_rmsd
                reorder_method = rmsd.reorder_brute
                #q_review = reorder_method(p_atoms, q_atoms, p_coord, q_coord)
                #q_coord = q_coord[q_review]
                #q_atoms = q_atoms[q_review]
                result_rmsd = rotation_method(p_coord, q_coord)
                #if result_rmsd > 0.15:
                #    same_geom = 0
            else:
                result_rmsd = 'not done'
            logger.info(f'\t{self.name} high level rmsd: {result_rmsd}, '\
                        f'same(0.15): {geometry.equal_geom(self.species, dummy, 0.15)}, '\
                        f'corr: {geometry.matrix_corr(imagmode, imagmode_high):.2f}, '\
                        f'same: {same_geom}')
        else:
            same_geom = geometry.equal_geom(self.species, dummy, 0.1)

        # checking if L2 frequencies are okay
        err, freq = self.qc.get_qc_freq(self.log_name(1, conf), self.species.natom)
        if self.species.natom == 1:
            freq_ok = 1
        elif len(freq) == 1 and freq[0] == 0:
            freq_ok = 0
        elif self.species.wellorts == 0 and freq[0] > -1. * self.par['imagfreq_threshold']:
            freq_ok = 1
            if freq[0] < 0.:
                logger.warning(f'Negative frequency {freq[0]} cm-1 detected in '
                            f'{self.name}. Flipped to {-freq[0]}.')
                freq[0] *= -1.
        elif self.species.wellorts == 1 \
                and not (np.count_nonzero(np.array(freq) < 0) >= 3  # Three or more imag frequencies
                        or np.count_nonzero(np.array(freq) < -1 * self.par['imagfreq_threshold']) >= 2  # More than one imaginary frequency beyond the threshold
                        or np.count_nonzero(np.array(freq) < 0) == 0):  # No imaginary frequencies:
            freq_ok = 1
        else:
            freq_ok = 0

        if conf == -1:
            # update properties for base structure
            if same_geom and freq_ok:
                err, self.species.geom = self.qc.get_qc_geom(self.log_name(1), self.species.natom)
                err, self.species.energy = self.qc.get_qc_energy(self.log_name(1))
                err, self.species.freq = self.qc.get_qc_freq(self.log_name(1), self.species.natom)   # TODO use fr variable
                err, self.species.zpe = self.qc.get_qc_zpe(self.log_name(1))
                if self.par['multi_conf_tst'] == 0:
                    self.shigh = 1
                else:
                    self.shigh = 0.5
            else:
                # geometry diverged to other structure
                if not same_geom:
                    logger.warning('High level optimization converged to different structure for {}, related channels are deleted.'.format(self.name))
                if not freq_ok:
                    logger.warning('Wrong number of imaginary frequencies for {}, related channels are deleted.'.format(self.name))
                self.shigh = -999
        else:
            # update property of conformers
            inx = self.species.conformer_index.index(conf) 
            if same_geom and freq_ok:
                err, self.species.conformer_geom[inx] = self.qc.get_qc_geom(self.log_name(1, conf=conf), self.species.natom)
                err, self.species.conformer_energy[inx] = self.qc.get_qc_energy(self.log_name(1, conf=conf))
                err, self.species.conformer_freq[inx] = self.qc.get_qc_freq(self.log_name(1, conf=conf), self.species.natom)   # TODO use fr variable
                err, zpe = self.qc.get_qc_zpe(self.log_name(1, conf=conf))
                self.species.conformer_zeroenergy[inx] = self.species.conformer_energy[inx] + zpe
            else:
                self.species.conformer_index[inx] = -999
                logger.warning(f'High level optimization failed for {self.log_name(1, conf=conf)}')

        return
