from __future__ import print_function
import os
import copy
import logging
import time
import numpy as np

from kinbot import frequencies
from kinbot import geometry
from kinbot import symmetry
from kinbot.conformers import Conformers
from kinbot.hindered_rotors import HIR
from kinbot.molpro import Molpro
from kinbot import reader_gauss
from kinbot.stationary_pt import StationaryPoint

class Optimize:
    """
    This class does the following:

    1. Conformational search of the species
    2. High level optimization and freq calc of the species
    3. Hindered rotor scans
    4. Repeat steps 2-3 as long as lower energy structures are found

    TODO: find better name for this module and class
    """

    def __init__(self, species, par, qc, wait=0):
        self.species = species
        try:
            delattr(self.species, 'cycle_chain')
        except AttributeError:
            logging.info("{} has no cycle_chain attribute to delete".format(self.species.chemid))
        if self.species.wellorts:
            self.species.characterize(bond_mx=self.species.bond)
        else:
            self.species.characterize()
        self.par = par
        self.qc = qc
        # wait for all calculations to finish before returning
        self.wait = wait
        # high level job name
        if self.species.wellorts:
            self.job_high = self.species.name + '_high'
            self.job_hir = 'hir/' + self.species.name + '_hir_'
        else:
            self.job_high = str(self.species.chemid) + '_well_high'
            self.job_hir = 'hir/' + str(self.species.chemid) + '_hir_'
        # status of the various parts
        # -1: not yet started
        #  0: running
        #  1: finished
        # -999:failed
        self.scycconf = -1
        self.sconf = -1
        self.shigh = -1
        self.shir = -1

        # restart counter: number of times the high-level and hir calculations
        # has been restarted in case a lower energy structure has been found
        self.restart = 0
        # maximum restart count
        self.max_restart = par.par['rotation_restart']

    def do_optimization(self):
        #print("start do_opt {}".format(self.species.chemid))
        while 1:
            # do the conformational search
            if self.par.par['conformer_search'] == 1:
                if self.scycconf == -1 and self.sconf == -1:
                    # conformational analysis has to be started
                    logging.info('\tStarting conformational search of {}'.format(self.species.chemid))
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
                # do the open chain part of the molecule
                if self.scycconf == 1:
                    # do open chain part if cyclic part is done
                    if self.sconf == -1:
                        # open chain part has not started yet
                        for geom in self.species.confs.cyc_conf_geoms:
                            # take all the geometries from the cyclic part
                            # generate the conformers for the current geometry
                            self.species.confs.generate_conformers(0, geom)
                        # set conf status to running
                        self.sconf = 0
                    if self.sconf == 0:
                        # conformational search is running
                        # check if the conformational search is done
                        status, lowest_conf, geom, low_energy, conformers, energies = self.species.confs.check_conformers(wait=self.wait)
                        if status == 1:
     
                            #create stationary points for chirality check
                            well0_stpt = StationaryPoint(name='well0', charge=self.par.par['charge'], mult=self.par.par['mult'],                                                                                natom=self.species.natom, atom=self.species.atom, geom=self.species.geom)
                            lowConf_stpt = StationaryPoint(name='conf', charge=self.par.par['charge'], mult=self.par.par['mult'],                                                                                 natom=self.species.natom, atom=self.species.atom, geom=geom)

                            #characterize & generate chirality label for well0 & conf
                            well0_stpt.characterize()
                            well0_stpt.bond = self.species.bond
                            well0_chiral = well0_stpt.calc_chiral()
            
                            lowConf_stpt.characterize()
                            lowConf_stpt.bond = self.species.bond
                            lowConfChiral = lowConf_stpt.calc_chiral()
                            well0Chiral_str = ' '.join(str(val) for val in well0_chiral)
                            lowConfChiral_str = ' '.join(str(val) for val in lowConfChiral)

                            print("{}\nwell0: {}\nlow conf: {}".format(self.species.chemid, well0Chiral_str, lowConfChiral_str))
                            npe = np.array(energies)
                            badconfs = []
                            #if well0Chiral_str != lowConfChiral_str:
			    #logging.info("Low energy conformer chirality differs from well chirality")
			    #logging.info("\tChecking conformation of other conformers")
			    print("Checking conformation of other conformers")
                            print("conformers: {}".format(len(conformers)))
			    
                            for i, conf in enumerate(conformers):
			        print(i)
                                #conf_stpt = StationaryPoint(name='conf', charge=self.par.par['charge'], mult=self.par.par['mult'],                                                                                                                   natom=self.species.natom, atom=self.species.atom, geom=conf)
			        #conf_stpt.characterize()
			        #conf_stpt.bond = self.species.bond
			        #conf_stptChiral = conf_stpt.calc_chiral()
			        #conf_stptChiralStr = ' '.join(str(val) for val in conf_stptChiral)
			        #if conf_stptChiralStr == well0Chiral_str:
				#    print("conf {} matches chirality".format(i))
			        #else:
				#    print("conf {} does NOT match chirality".format(i))
				#    badconfs.append(i)
                                    
                                #print("done with confs")
                                #print(badconfs, i)
                                
                            print("done with confs") 
                            """
                            #Implement the following
                            # 1. check next conf chirality
                            # 2. if chirality changes remove conf & energy + create log
                            # 3. if lowest E conf changes chirality check array for lowest energy & report lowest E as new conf
                            #if well0Chiral_str != lowConfChiral_str:
			    print("reading through conformers")
			    for i, conf in enumerate(conformers):
			        conf_stpt.characterize()
			        conf_stpt.bond = self.species.bond
			        conf_stptChiral = conf_stpt.calc_chiral()
			        conf_stptChiralStr = ' '.join(str(val) for val in conf_stptChiral)
			        print("{} conf {}, energy: {} chiral: {}".format(self.species.chemid, i, energies[i], conf_stptChiral))
			    while len(energies) > 0:
			        lowe_index = energies.index(np.amin(npe))
			        print("min energy at index {}, e = {}".format(lowe_index, energies[lowe_index]))
			        conformers.pop(lowe_index)
			        energies.pop(lowe_index)
			        npe = np.delete(npe, lowe_index)
			        print(len(energies), len(npe), len(conformers))
                            """
                            print("past conf check")
                            # conf search is done
                            if self.species.wellorts:
                                self.name = self.species.name
                            else:
                                self.name = self.species.chemid
                            logging.info("lowest energy conformer for species: {} is number {}".format(self.name, lowest_conf))
                            # save lowest energy conformer as species geometry
                            self.species.geom = geom
                            # save lowest energy conformer energy
                            self.species.energy = low_energy
                            # set conf status to finished
                            self.sconf = 1

            else:
                # no conf search necessary, set status to finished
                self.sconf = 1
            if self.sconf == 1:  # conf search is finished
                # if the conformers were already done in a previous run
                if self.par.par['conformer_search'] == 1:
                    status, lowest_conf, geom, low_energy, conformers, energies = self.species.confs.check_conformers(wait=self.wait)
                    # perform conformer check at this point
                    filteredConf = [] 
                    #for conf in conformers:
                        #print(conf)
                        #confStPt = StationaryPoint(self.name, self.species.charge, self.species.mult, self.species.natom, self.species.atom, geom, self.species.wellorts)
                        
                while self.restart < self.max_restart:
                    # do the high level calculations
                    if self.par.par['high_level'] == 1:
                        if self.shigh == -1:
                            if self.species.wellorts:
                                name = self.species.name
                            else:
                                name = self.species.chemid
                            # high level calculation did not start yet
                            logging.info('\tStarting high level optimization of {}'.format(name))
                            if self.species.wellorts:
                                # do the high level optimization of a ts
                                self.qc.qc_opt_ts(self.species, self.species.geom, high_level=1)
                            else:
                                # do the high level optimization of a well
                                self.qc.qc_opt(self.species, self.species.geom, high_level=1)
                            self.shigh = 0  # set the high status to running
                        if self.shigh == 0:
                            # high level calculation is running
                            # check if it is finished
                            status = self.qc.check_qc(self.job_high)
                            if status == 'error':
                                # found an error
                                logging.info('\tHigh level optimization failed for {}'.format(self.species.name))
                                self.shigh = -999
                            if status == 'normal':
                                # finished successfully
                                err, new_geom = self.qc.get_qc_geom(self.job_high, self.species.natom, wait=self.wait)

                                if self.species.wellorts:  # for TS we need reasonable geometry agreement and normal mode correlation
                                    if self.par.par['conformer_search'] == 0:
                                        fr_file = self.fr_file_name(0)  # name of the original TS file

                                    else:
                                        fr_file = 'conf/{}_{}'.format(self.fr_file_name(0), lowest_conf)
                                    if self.qc.qc == 'gauss':
                                        imagmode = reader_gauss.read_imag_mode(fr_file, self.species.natom)
                                    fr_file = self.fr_file_name(1)
                                    if self.qc.qc == 'gauss':
                                        imagmode_high = reader_gauss.read_imag_mode(fr_file, self.species.natom)
                                    # either geom is roughly same with closely matching imaginary modes, or geometry is very close
                                    # maybe we need to do IRC at the high level as well...
                                    same_geom = ((geometry.matrix_corr(imagmode, imagmode_high) > 0.9) and \
                                            (geometry.equal_geom(self.species.bond, self.species.geom, new_geom, 0.3))) \
                                            or (geometry.equal_geom(self.species.bond, self.species.geom, new_geom, 0.15))
                                else:
                                    same_geom = geometry.equal_geom(self.species.bond, self.species.geom, new_geom, 0.1)

                                if same_geom:
                                    # geometry is as expected and normal modes are the same for TS
                                    err, self.species.geom = self.qc.get_qc_geom(self.job_high, self.species.natom)
                                    err, self.species.energy = self.qc.get_qc_energy(self.job_high)
                                    err, self.species.freq = self.qc.get_qc_freq(self.job_high, self.species.natom)
                                    err, self.species.zpe = self.qc.get_qc_zpe(self.job_high)
                                    self.shigh = 1
                                else:
                                    # geometry diverged to other structure
                                    logging.info('\tHigh level optimization converged to different structure for {}, related channels are deleted.'.format(self.species.name))
                                    self.shigh = -999
                              
                    else:
                        # no high-level calculations necessary, set status to finished
                        self.shigh = 1
                    logging.info("done with conformer search for {}".format(self.species.name))
                    if self.shigh == 1:
                        # do the HIR calculation
                        if self.par.par['rotor_scan'] == 1:
                            if self.shir == -1:
                                # hir not stated yet
                                logging.info('\tStarting hindered rotor calculations of {}'.format(self.species.name))
                                self.species.hir = HIR(self.species, self.qc, self.par)
                                self.species.hir.generate_hir_geoms(copy.deepcopy(self.species.geom))
                                self.shir = 0
                            if self.shir == 0:
                                # hir is running
                                # check if it is done:
                                status = self.species.hir.check_hir(wait=self.wait)
                                if status:
                                    if len(self.species.hir.hir_energies) > 0:
                                        # check if along the hir potential a structure was found with a lower energy
                                        min = self.species.hir.hir_energies[0][0]
                                        min_rotor = -1
                                        min_ai = -1
                                        for rotor in range(len(self.species.dihed)):
                                            for ai in range(self.species.hir.nrotation):
                                                # use a 0.1kcal/mol cutoff for numerical noise
                                                if self.species.hir.hir_energies[rotor][ai] < min - 1.6E-4:
                                                    min = self.species.hir.hir_energies[rotor][ai]
                                                    min_rotor = rotor
                                                    min_ai = ai
                                        if min_rotor > -1:
                                            self.restart += 1
                                            if self.restart < self.max_restart:
                                                # lower energy structure found
                                                logging.info('\t\tLower energy found during hindered rotor scan for {}'.format(self.species.name))
                                                logging.info('\t\tRestart number: ' + str(self.restart))
                                                logging.info('\t\tRotor: ' + str(min_rotor))
                                                logging.info('\t\tScan point: ' + str(min_ai))
                                                job = self.job_hir + str(min_rotor) + '_' + str(min_ai).zfill(2)

                                                err, self.species.geom = self.qc.get_qc_geom(job, self.species.natom)
                                                # delete the high_level log file and the hir log files
                                                if os.path.exists(self.job_high + '.log'):
                                                    # logging.info("\t\t\tRemoving file " + self.job_high + '.log')
                                                    os.remove(self.job_high + '.log')
                                                for rotor in range(len(self.species.dihed)):
                                                    for ai in range(self.species.hir.nrotation):
                                                        if os.path.exists(self.job_hir + str(rotor) + '_' + str(ai).zfill(2) + '.log'):
                                                            # logging.info("\t\t\tRemoving file " + self.job_hir + str(rotor) + '_' + str(ai).zfill(2) + '.log')
                                                            os.remove(self.job_hir + str(rotor) + '_' + str(ai).zfill(2) + '.log')
                                                # set the status of high and hir back to not started
                                                self.shigh = -1
                                                self.shir = -1
                                            else:
                                                logging.info('\t\tLower energy found, but readched max restart for {}'.format(self.species.name))
                                                self.shir = 1
                                        else:
                                            self.shir = 1
                                    else:
                                        self.shir = 1
                        else:
                            # no hir calculations necessary, set status to finished
                            self.shir = 1
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
                fr_file = self.species.name
                if not self.species.wellorts:
                    fr_file = str(self.species.chemid)
                    fr_file += '_well'
                if self.par.par['high_level']:
                    fr_file += '_high'
                fr_file = self.fr_file_name(self.par.par['high_level'])
                hess = self.qc.read_qc_hess(fr_file, self.species.natom)
                self.species.kinbot_freqs, self.species.reduced_freqs = frequencies.get_frequencies(self.species, hess, self.species.geom)

                # write the molpro input and read the molpro energy, if available
                if self.par.par['single_point_qc'] == 'molpro':
                    molp = Molpro(self.species, self.par)
                    molp.create_molpro_input()
                    molp.create_molpro_submit()
                    print(self.par.par['single_point_key'])
                    status, molpro_energy = molp.get_molpro_energy(self.par.par['single_point_key'])
                    if status:
                        self.species.energy = molpro_energy

                # delete unnecessary files
                if self.par.par['delete_intermediate_files'] == 1:
                    self.delete_files()

            if self.wait:
                if self.shir == 1 or self.shigh == -999:
                    return 0
                time.sleep(1)
            else:
                return 0
            

    def delete_files(self):
        # job names
        names = []
        zf = self.par.par['zf']
        if self.species.wellorts:
            names.append(self.species.name)
            names.append(self.species.name + '_high')
            names.append(self.species.name + '_IRC_F')
            names.append(self.species.name + '_IRC_R')
            names.append(self.species.name + '_IRC_F_prod')
            names.append(self.species.name + '_IRC_R_prod')

            if self.par.par['high_level'] == 1:
                for count in range(self.species.hir.nrotation):
                    for rot_num in range(self.par.par['nrotation']):
                        names.append('hir/' + self.species.name + '_hir_' + str(count) + '_' + str(rot_num).zfill(2))
            if self.par.par['conformer_search'] == 1:
                for count in range(self.species.confs.conf):
                    names.append('conf/' + self.species.name + '_' + str(count).zfill(zf))
                for count in range(self.species.confs.cyc_conf):
                    for num in range(self.species.confs.cyc_conf_index[count]):
                        names.append('conf/' + self.species.name + '_r' + str(count).zfill(zf) + '_' + str(num).zfill(zf))
        else:
            names.append(str(self.species.chemid) + '_well')
            names.append(str(self.species.chemid) + '_well_high')
            if self.par.par['high_level'] == 1:
                for count in range(self.species.hir.nrotation):
                    for rot_num in range(self.par.par['nrotation']):
                        names.append('hir/' + str(self.species.chemid) + '_hir_' + str(count) + '_' + str(rot_num).zfill(2))
            if self.par.par['conformer_search'] == 1:
                for count in range(self.species.confs.conf):
                    names.append('conf/' + str(self.species.chemid) + '_' + str(count).zfill(zf))
                for count in range(self.species.confs.cyc_conf):
                    for num in range(self.species.confs.cyc_conf_index[count]):
                        names.append('conf/' + self.species.name + '_r' + str(count).zfill(zf) + '_' + str(num).zfill(zf))

        extensions = ['chk', 'py', 'sbatch']

        for name in names:
            for ext in extensions:
                # delete file
                file = '.'.join([name, ext])
                try:
                    os.remove(file)
                # except FileNotFoundError:
                except:
                    pass

    def fr_file_name(self, high):
        fr_file = self.species.name
        if not self.species.wellorts:
            fr_file += '_well'
        # if self.par.par['high_level']:
        if high:
            fr_file += '_high'

        return(fr_file)
