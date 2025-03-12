import os, sys
import shutil
import time
import logging
import copy
import itertools

import numpy as np

from kinbot import constants
#from kinbot import filecopying
from kinbot import pes
from kinbot import postprocess
from kinbot import reac_family
from kinbot import cheminfo
from kinbot.irc import IRC
from kinbot.optimize import Optimize
from kinbot.reac_General import GeneralReac
from kinbot.stationary_pt import StationaryPoint
from kinbot.molpro import Molpro
from ase.db import connect
from ase import Atoms
from kinbot.utils import reorder_coord


logger = logging.getLogger('KinBot')

class ReactionGenerator:
    '''
    This class generates the reactions using the qc codes
    It builds an initial guess for the ts, runs that ts towards a saddle point
    and does IRC calculations
    '''

    def __init__(self, species, par, qc, input_file):
        self.species = species
        self.par = par
        self.inp = input_file
        self.qc = qc

    def generate(self):
        '''
        Creates the input for each reaction, runs them, and tests for success.
        If successful, it creates the barrier and product objects.
        It also then does the conformational search, and finally, the hindered rotor scans.
        To make the code the most efficient, all of these happen in parallel, in a sense that
        the jobs are not waiting for each other. E.g., one reaction can still be in the stage
        of TS search, while the other can be already at the hindered rotor scan. This way,
        all cores are occupied efficiently.

        The switching between the various stages are done via the reac_ts_done variable.
        0: initiate the TS search
        1: check barrier height and errors in TS, and initiates normal mode displacement test, start the irc calculations
        2: submit product optimization
        3: submit the frequency calculation
        4: do the optimization of the ts and the products
        5: follow up on the optimizations
        6: finalize calculations, check for wrong number of negative frequencies

        If at any times the calculation fails, reac_ts_done is set to -999.
        If all steps are successful, reac_ts_done is set to -1.
        '''
        deleted = []
        if len(self.species.reac_inst) > 0:
            alldone = 1
        else:
            alldone = 0

        for ro in self.species.reac_obj:
            ro.prod_done = 0
            ro.valid_prod = []

        frag_unique = []

        while alldone:
            for index, instance in enumerate(self.species.reac_inst):
                obj = self.species.reac_obj[index]
                if obj.instance_name in self.par['skip_reactions'] \
                        and self.species.reac_ts_done[index] != -999:
                    logger.info(f'\tRemoving reaction {obj.instance_name}.')
                    self.species.reac_ts_done[index] = -999
                # START REACTION SEARCH
                if self.species.reac_ts_done[index] == 0 and self.species.reac_step[index] == 0:
                    # verify after restart if search has failed in previous kinbot run
                    status = self.qc.check_qc(obj.instance_name)
                    if status == 'error':
                        logger.info('\tRxn search failed for {} - had error status on previous round.'
                                     .format(obj.instance_name))
                        self.species.reac_ts_done[index] = -999
                if self.species.reac_type[index] == 'hom_sci' and self.species.reac_ts_done[index] == 0:  # no matter what, set to 2
                    # somewhat messy manipulation to force the new bond matrix for hom_sci
                    obj.irc_prod = StationaryPoint(f'{instance}_prod', self.species.charge, self.species.mult,
                                                   atom=self.species.atom, geom=self.species.geom, wellorts=0)
                    obj.irc_prod.characterize()
                    obj.irc_prod.bonds[0][obj.instance[0]][obj.instance[1]] = 0  # delete bond
                    obj.irc_prod.bonds[0][obj.instance[1]][obj.instance[0]] = 0  # delete bond
                    obj.irc_prod.bond[obj.instance[0]][obj.instance[1]] = 0  # delete bond
                    obj.irc_prod.bond[obj.instance[1]][obj.instance[0]] = 0  # delete bond
                    self.species.reac_ts_done[index] = 2

                if self.species.reac_ts_done[index] == 0:  # ts search is ongoing
                    if obj.scan == 0:  # don't do a scan of a bond
                        if self.species.reac_step[index] == obj.max_step + 1:
                            status, freq = self.qc.get_qc_freq(obj.instance_name, self.species.natom)
                            if status == 0 and (np.count_nonzero(np.array(freq) < 0) >= 3  # Three or more imag frequencies
                                                or np.count_nonzero(np.array(freq) < -1 * self.par['imagfreq_threshold']) >= 2  # More than one imaginary frequency beyond the threshold
                                                or np.count_nonzero(np.array(freq) < 0) == 0):  # No imaginary frequencies
                                logger.info(f'\tReaction search failed for {obj.instance_name}: '
                                            'Wrong number of imaginary frequencies.')
                                self.species.reac_ts_done[index] = -999
                            elif status == 0:
                                self.species.reac_ts_done[index] = 1
                            elif status == -1: 
                                logger.info(f'\tReaction search failed for {obj.instance_name}.')
                                self.species.reac_ts_done[index] = -999
                        else:
                            self.species.reac_step[index] = reac_family.carry_out_reaction(
                                                            obj, self.species.reac_step[index], self.par['qc_command'],
                                                            bimol=self.par['bimol'])
                            if self.species.reac_step[index] == -1:
                                self.species.reac_ts_done[index] = -999
                                logger.info(f'\tReaction search failed for {obj.instance_name}: '
                                            'Invalid geometry.')

                    elif obj.scan == 1:  # do a bond scan 

                        if (self.species.reac_step[index] == self.par['scan_step'] + 1):
                            status, freq = self.qc.get_qc_freq(obj.instance_name, self.species.natom)
                            if status == 0 and (np.count_nonzero(np.array(freq) < 0) >= 3  # More than two imag frequencies
                                                or np.count_nonzero(np.array(freq) < -1 * self.par['imagfreq_threshold']) >= 2  # More than one imaginary frequency beyond the threshold
                                                or np.count_nonzero(np.array(freq) < 0) == 0):  # No imaginary frequencies
                                logger.info(f'\tReaction search failed for {obj.instance_name}: '
                                            'wrong number of imaginary frequencies.')
                                self.species.reac_ts_done[index] = -999
                            elif status == 0:
                                self.species.reac_ts_done[index] = 1
                            elif status == -1:
                                logger.info(f'\tRxn search using scan failed for {obj.instance_name} '
                                            'in TS optimization stage.')
                                self.species.reac_ts_done[index] = -999
                        else:  
                            # First point of the scan
                            if self.species.reac_step[index] == 0:
                                command = self.par['qc_command']
                                self.species.reac_step[index] = reac_family.carry_out_reaction(
                                                                obj, self.species.reac_step[index], command,
                                                                bimol=self.par['bimol'])
                                if self.species.reac_step[index] == -1:
                                    self.species.reac_ts_done[index] = -999
                                    logger.info('\tRxn search failed for {} because of 0 0 0 geometry.'
                                                 .format(obj.instance_name))
                            # Middle points of the scan
                            if (self.species.reac_step[index] < self.par['scan_step'] and obj.family_name):
                                status = self.qc.check_qc(obj.instance_name)
                                if status == 'error':
                                    logger.info('\tRxn search using scan failed for {} in step {}'
                                                 .format(obj.instance_name, self.species.reac_step[index]))
                                    self.species.reac_ts_done[index] = -999
                                else:
                                    err, energy = self.qc.get_qc_energy(obj.instance_name)
                                    if err == 0:
                                        self.species.reac_scan_energy[index].append(energy)
                                        logger.debug(f'Scan energy for {obj.instance_name} in step {self.species.reac_step[index]}:')
                                        logger.debug(f'{self.species.reac_scan_energy[index][-1]} Hartree.')
                                        # need at least 3 points for a maximum
                                        if len(self.species.reac_scan_energy[index]) >= 3:
                                            ediff = np.diff(self.species.reac_scan_energy[index])
                                            if ediff[-1] < 0 and ediff[-2] > 0:  # max
                                                logger.info(f'\tMaximum found for {obj.instance_name}.')
                                                e_in_kcal = [constants.AUtoKCAL * (self.species.reac_scan_energy[index][ii] - 
                                                               self.species.reac_scan_energy[index][0])
                                                               for ii in range(len(self.species.reac_scan_energy[index]))]
                                                e_in_kcal = np.round(e_in_kcal, 2)
                                                logger.info(f'\tEnergies: {e_in_kcal}')
                                                logger.debug(f'Derivatives: {ediff}')
                                                self.species.reac_step[index] = self.par['scan_step']  # ending the scan
                                            if len(ediff) >= 3:
                                                if any([edf == 0 for edf in ediff[-2:]]):
                                                    logger.warning(f'Calculation failed in the bond scan of {obj.instance_name}.')
                                                    self.species.reac_ts_done[index] = -999
                                                    continue
                                                if 10. * (ediff[-3] / ediff[-2]) < (ediff[-2] / ediff[-1]):  # sudden change in slope
                                                    logger.info(f'\tSudden change in slope for {obj.instance_name}.')
                                                    logger.info(f'\tRelative energies (kcal/mol): {self.species.reac_scan_energy[index]}')
                                                    logger.debug(f'Derivatives: {ediff}')
                                                    self.species.reac_step[index] = self.par['scan_step']  # ending the scan
                                     
                                        # scan continues, and if reached scan_step, then goes for full optimization
                                        self.species.reac_step[index] = reac_family.carry_out_reaction(
                                                                        obj, self.species.reac_step[index], command,
                                                                        bimol=self.par['bimol'])
                                        if self.species.reac_step[index] == -1:
                                            self.species.reac_ts_done[index] = -999
                                            logger.info('\tRxn search failed for {} because of 0 0 0 geometry.'
                                                         .format(obj.instance_name))
                            else:  # the last step was reached, and no max or inflection was found
                                status = self.qc.check_qc(obj.instance_name)
                                if status == 'running':
                                    continue
                                logger.info('\tRxn search using scan failed for {}, no saddle guess found.'
                                            .format(obj.instance_name))
                                db = connect('{}/kinbot.db'.format(os.getcwd()))
                                # error line, H atom is just placeholder
                                db.write(Atoms('H'), name=obj.instance_name, data={'status': 'error'})
                                # this is copied here so that a non-AM1 file is in place
                                shutil.copy(f'{os.getcwd()}/{self.species.chemid}_well.log', f'{os.getcwd()}/{obj.instance_name}.log')
                                self.species.reac_ts_done[index] = -999

                elif self.species.reac_ts_done[index] == 1:
                    status = self.qc.check_qc(obj.instance_name)
                    if status == 'running':
                        continue
                    elif status == 'error':
                        logger.info('\tRxn search failed (gaussian error) for {}'
                                     .format(obj.instance_name))
                        self.species.reac_ts_done[index] = -999
                    else:
                        # check the barrier height:
                        ts_energy = self.qc.get_qc_energy(obj.instance_name)[1]
                        ts_zpe = self.qc.get_qc_zpe(obj.instance_name)[1]
                        if self.species.reac_type[index] == 'R_Addition_MultipleBond' \
                                and self.qc.qc != 'nn_pes':
                            ending = 'well_mp2'
                            thresh = self.par['barrier_threshold']  # need to fix for mp2 specific
                        elif self.species.reac_type[index] == 'barrierless_saddle' \
                                and self.qc.qc != 'nn_pes':
                            ending = 'well_bls'
                            thresh = self.par['barrier_threshold']
                        else:
                            ending = 'well'
                            thresh = self.par['barrier_threshold']
                        sp_energy = self.qc.get_qc_energy('{}_{}'.format(str(self.species.chemid), ending))[1]
                        sp_zpe = self.qc.get_qc_zpe('{}_{}'.format(str(self.species.chemid), ending))[1]
                        try:
                            barrier = (ts_energy + ts_zpe - sp_energy - sp_zpe) * constants.AUtoKCAL
                        except TypeError:
                            logger.error(f'Faulty calculations, check or delete files for {obj.instance_name}.')
                            sys.exit(-1)
                        if barrier > thresh:
                            logger.info('\tRxn barrier too high ({0:.2f} kcal/mol) at L1 for {1}'
                                         .format(barrier, obj.instance_name))
                            self.species.reac_ts_done[index] = -999
                        else:
                            obj.irc = IRC(obj, self.par)  
                            irc_status = obj.irc.check_irc()
                            if 0 in irc_status:
                                logger.info('\tRxn barrier is ({0:.2f} kcal/mol) at L1 for {1}'
                                             .format(barrier, obj.instance_name))
                                # No IRC started yet, start the IRC now
                                logger.info('\tStarting IRC calculations for {}'
                                             .format(obj.instance_name))
                                obj.irc.do_irc_calculations()
                            elif irc_status[0] == 'running' or irc_status[1] == 'running':
                                continue
                            else:
                                # IRC's have successfully finished, have an error, in any case
                                # read the geometries and try to make products out of them
                                # verify which of the ircs leads back to the reactant, if any
                                logger.info('\tRxn barrier is {0:.2f} kcal/mol for {1}'
                                             .format(barrier, obj.instance_name))
                                obj.irc_prod = obj.irc.irc2stationary_pt()
                                if obj.irc_prod == 0:
                                    logger.info('\tNo product found for {}'.format(obj.instance_name))
                                    self.species.reac_ts_done[index] = -999
                                else:
                                    self.species.reac_ts_done[index] = 2

                elif self.species.reac_ts_done[index] == 2:
                    # obj.valid_prod: list marking fragments for deletion (if breaks apart or changes)
                    # frag_unique: list of unique fragments across all reactions for this well
                    # obj.products: list of products for given reaction, which includes changes and further dissociation 
                    if obj.prod_done == 0:  # not started optimization yet
                        # identify bimolecular products and wells from IRC - do it once
                        obj.products, _ = obj.irc_prod.start_multi_molecular(vary_charge=True)
                        if self.species.charge == 0:
                            logger.info(f'\tBased on the end of IRC, reaction {obj.instance_name} leads to products '
                                        f'{[fr.chemid for fr in obj.products]} {[fr.atom for fr in obj.products]}')
                        else:
                            logger.info(f'\tBased on the end of IRC, reaction {obj.instance_name} leads to products '
                                        f'{[fr.chemid for fr in obj.products]} {[fr.atom for fr in obj.products]} '
                                        '(including all possible charge distributions)')

                        self.equate_identical(obj.products)
                        obj.valid_prod = len(obj.products) * [True]

                        # make the geom of products in frag_unique the one from the multi_molecular (not optimized)
                        self.equate_unique(obj.products, frag_unique)
                        obj.prod_done = 1

                    for frag in obj.products:
                        self.qc.qc_opt(frag, frag.geom)
                        e, _ = self.qc.get_qc_geom(str(frag.chemid) + '_well', frag.natom) # check if finished without updating geom
                        if e == 1:  # it's running
                            continue

                    # initial fragment calculations finished, reading results...
                    hom_sci_energy = 0
                    products_orig = [copy.copy(opr) for opr in obj.products] 
                    ndone = 0
                    for fragii, frag in enumerate(products_orig):
                        if frag.chemid == self.species.chemid:
                            logger.info(f'Product in {obj.instance_name} is identical to the reactant. Reaction deleted.')
                            self.species.reac_ts_done[index] = -999 
                            break
                        # do not look at already invalid fragments again
                        elif not obj.valid_prod[fragii]:
                            ndone += 1
                            continue
                        chemid_orig = frag.chemid
                        e, frag.geom, frag.atom = self.qc.get_qc_geom(
                            str(frag.chemid) + '_well',
                            frag.natom,
                            reorder=True)
                        if e < 0:
                            logger.info(f'\tProduct optimization failed for {obj.instance_name}, product {frag.chemid}')
                            self.species.reac_ts_done[index] = -999
                            ndone += 1
                        elif e == 1:
                            break
                        else:
                            ndone += 1
                            _, frag.energy = self.qc.get_qc_energy(str(frag.chemid) + '_well')
                            _, frag.zpe = self.qc.get_qc_zpe(str(frag.chemid) + '_well')
                            if self.species.reac_type[index] == 'hom_sci': # TODO energy is the sum of all possible fragments  
                                hom_sci_energy += frag.energy + frag.zpe
                            # Reinitialize rads and bonds
                            frag.reset_order()
                            # connectivity changed
                            if chemid_orig != frag.chemid:
                                for fri, fr in enumerate(obj.products):
                                    if fr.chemid == chemid_orig:
                                        obj.valid_prod[fri] = False
                                newfrags, _ = frag.start_multi_molecular(vary_charge=True)  
                                self.equate_identical(newfrags)
                                self.equate_unique(newfrags, frag_unique)
                                logger.warning(f'Product {chemid_orig} optimized to {[nf.chemid for nf in newfrags]} '
                                               f'in reaction {obj.instance_name}')
                                for nf in newfrags:
                                    obj.products.append(nf)
                                    obj.valid_prod.append(True)
                            else:
                                for fri, fr in enumerate(obj.products):
                                    if fr.chemid == chemid_orig:
                                        obj.products[fri].energy = frag.energy
                                        obj.products[fri].zpe = frag.zpe
                                        # Reorder the coordinates of frag in case the atom order is different
                                        diff = [p1 for p1, p2 in zip(obj.products[fri].atom, frag.atom) if p1 != p2]
                                        if diff != []:
                                            reorder_coord(mol_A=obj.products[fri],
                                                          mol_B=frag)
                                        obj.products[fri].geom = frag.geom

                    if ndone == len(obj.products) and self.species.reac_ts_done[index] != -999:  # all currently recognized fragments are done
                        # delete invalid ones
                        obj.products = list(np.array(obj.products)[obj.valid_prod])
                        if self.species.charge != 0:  # select the lower energy combination
                            # brute force all combinations for ions
                            combs = np.array([np.array(i) for i in itertools.product([0, 1], repeat = len(obj.products))])
                            val = 1000.
                            ens = np.array([i.energy for i in obj.products])
                            zpes = np.array([i.zpe for i in obj.products])
                            masses = np.array([i.mass for i in obj.products])
                            charges = np.array([i.charge for i in obj.products])
                            frag_symbols = [p.atom for p in obj.products]
                            low_e_comb = combs[0]
                            for comb in combs:
                                comb_symbols = []
                                for i, frag_symb in enumerate(frag_symbols):
                                    if comb[i]:
                                        comb_symbols.extend(frag_symb)
                                if sum(comb * charges) != self.species.charge \
                                        or sum(comb * masses) != self.species.mass \
                                        or sorted(comb_symbols) != sorted(self.species.atom):
                                    continue
                                rel_energy = constants.AUtoKCAL * (sum(comb * (ens + zpes)) - \
                                             (self.species.start_energy + self.species.start_zpe))
                                relev_comb = [p.chemid for i, p in enumerate(obj.products) if comb[i]]
                                if len(obj.products) > 2:
                                    logger.info(f'\tPossible ion combination and energy for {obj.instance_name} products: '
                                                f'{", ".join([str(c) for c in relev_comb])} '
                                                f'at {rel_energy:.1f} kcal/mol.')
                                if sum(comb * (ens + zpes)) < val:
                                    val = sum(comb * (ens + zpes))
                                    low_e_comb = comb
                            obj.products = list(np.array(obj.products)[low_e_comb.astype(bool)])
                        obj.irc_fragments = [copy.copy(this_frag) for this_frag in obj.products]
                        prods_energy = sum([p.energy + p.zpe for p in obj.products])

                        # Select reactions potentially leading to vdW wells
                        if len(obj.products) == 2 and 'hom_sci' not in obj.instance_name:
                            logger.info('\tChecking vdW well for {}.'.format(obj.instance_name))
                            obj.irc_prod.characterize()
                            e, obj.irc_prod.energy = self.qc.get_qc_energy(f'{obj.irc_prod.name}')  # e is the error code: should be 0 (success) at this point.
                            e, obj.irc_prod.zpe = self.qc.get_qc_zpe(f'{obj.irc_prod.name}')
                            e, obj.irc_prod.geom = self.qc.get_qc_geom(obj.irc_prod.name, obj.irc_prod.natom)
                            e, obj.irc_prod.freq = self.qc.get_qc_freq(obj.irc_prod.name, obj.irc_prod.natom) 
                            for this_frag in obj.irc_fragments:
                                e, this_frag.freq = self.qc.get_qc_freq(f'{this_frag.name}_well', this_frag.natom) 
                            fragments_energies = sum([(this_frag.energy + this_frag.zpe) for this_frag in obj.irc_fragments ])
                            obj.vdW_depth = (fragments_energies - (obj.irc_prod.energy + obj.irc_prod.zpe)) * constants.AUtoKCAL
                            if obj.vdW_depth > self.par['vdW_detection']:
                                logger.info(f'\tvdW well detected for {obj.irc_prod.name}. Depth: {np.round(obj.vdW_depth, 2)} kcal/mol.')
                                obj.do_vdW = True
                            else:
                                logger.info(f'\t{obj.irc_prod.name} was not succesfully identified as a vdW well. Well depth is {np.round(obj.vdW_depth, 2)} kcal/mol.')
                            
                        if self.species.reac_type[index] == 'hom_sci': # TODO energy is the sum of all possible fragments
                            hom_sci_energy = (prods_energy - self.species.start_energy - self.species.start_zpe) * constants.AUtoKCAL
                            if hom_sci_energy < self.par['barrier_threshold'] + self.par['hom_sci_threshold_add']:
                                self.species.reac_ts_done[index] = 3
                            else:
                                logger.info(f'\thom_sci energy is too high at {np.round(hom_sci_energy, 2)} kcal/mol for {obj.instance_name}')
                                self.species.reac_ts_done[index] = -999
                        else:
                            self.species.reac_ts_done[index] = 3

                elif self.species.reac_ts_done[index] == 3:
                    for frag in obj.products:
                        # # Reordering in case different fragment
                        e, frag.geom, frag.atom = self.qc.get_qc_geom(
                            str(frag.chemid) + '_well',
                            frag.natom,
                            reorder=True)
                        # Create a new stp instead of updating geom and atom
                        # because rads and bonds need to be changed
                        frag.reset_order()
                        # e, frag.geom = self.qc.get_qc_geom(str(frag.chemid) + '_well', frag.natom)

                    # Do the TS and product optimization
                    # make a stationary point object of the ts
                    bond_mx = np.zeros((self.species.natom, self.species.natom), dtype=int)
                    for i in range(self.species.natom):
                        for j in range(self.species.natom):
                            bond_mx[i][j] = max(self.species.bond[i][j], obj.irc_prod.bonds[0][i][j])

                    if self.species.reac_type[index] != 'hom_sci':
                        err, geom = self.qc.get_qc_geom(obj.instance_name, self.species.natom)
                        ts = StationaryPoint(obj.instance_name, self.species.charge, self.species.mult,
                                             atom=self.species.atom, geom=geom, wellorts=1)
                        err, ts.energy = self.qc.get_qc_energy(obj.instance_name)
                        err, ts.zpe = self.qc.get_qc_zpe(obj.instance_name)
                        err, ts.freq = self.qc.get_qc_freq(obj.instance_name, self.species.natom)
                        ts.distance_mx()
                        ts.bond_mx()
                        ts.bond = np.maximum(ts.bond, bond_mx)
                        # -1: broken, 0: no change, +1: formed
                        ts.reac_bond = np.array(obj.irc_prod.bond01) - np.array(self.species.bond01) 
                        ts.find_cycle()
                        ts.find_conf_dihedral()
                        obj.ts = ts
                        # do the ts optimization
                        obj.ts_opt = Optimize(obj.ts, self.par, self.qc)
                        obj.ts_opt.do_optimization()
                    else:
                        obj.ts = copy.copy(obj.species)  # the TS will be for now the species itself
                        obj.ts.wellorts = 1

                    # do the vdW optimizations
                    if obj.do_vdW:
                        # do the high level optimization if requested
                        obj.irc_prod_opt = Optimize(obj.irc_prod, self.par, self.qc, just_high=True)
                        obj.irc_prod_opt.do_optimization()
                        if obj.irc_prod_opt.shigh == -999: 
                            logger.info('\tVdW search failed for {}, prod_opt shigh fail for {}.'
                                            .format(obj.irc_prod.name, obj.irc_prod.chemid))
                            obj.do_vdW = False
                        # non-conformationally optimized fragments high-level qc to be added here if we want to
                                                  
                    # do the products optimizations
                    temp_prod_opt = []  # holding the optimization objects temporarily
                    for st_pt in obj.products:
                        new = 1
                        # do the products optimizations
                        # check for products of other reactions that are the same as this product
                        # in the case such products are found, use the same Optimize object for both
                        for i, inst_i in enumerate(self.species.reac_inst):
                            if not new:
                                break
                            if i != index:
                                obj_i = self.species.reac_obj[i]
                                if self.species.reac_ts_done[i] > 2:
                                    for j, st_pt_i in enumerate(obj_i.products):
                                        if st_pt_i.chemid == st_pt.chemid:
                                            if len(obj_i.prod_opt) > j:
                                                prod_opt = obj_i.prod_opt[j]
                                                new = 0
                                                break
                        if new:
                            prod_opt = Optimize(st_pt, self.par, self.qc)
                            prod_opt.do_optimization()
                            if prod_opt.shigh == -999:
                                logger.info('\tRxn search failed for {}, prod_opt shigh fail for {}.'
                                             .format(obj.instance_name, prod_opt.species.chemid))
                                self.species.reac_ts_done[index] = -999
                                #break  # breaks so that other species is not looked at
                        temp_prod_opt.append(prod_opt)
                    if self.species.reac_ts_done[index] != -999:
                        for tpo in temp_prod_opt:
                            obj.prod_opt.append(tpo)

                    if self.species.reac_ts_done[index] != -999:  # so we don't reset faulty calculation
                        for st_pt in obj.products:
                            # section where comparing products in same reaction occurs
                            if len(obj.prod_opt) > 0:
                                for j, st_pt_opt in enumerate(obj.prod_opt):
                                    if st_pt.chemid == st_pt_opt.species.chemid:
                                        if len(obj.prod_opt) > j:
                                            prod_opt = obj.prod_opt[j]
                                            break

                        self.species.reac_ts_done[index] = 4

                elif self.species.reac_ts_done[index] == 4:
                    # check up on the TS and product optimizations
                    opts_done = 1
                    fails = 0
                    # check if ts and vdW are done
                    if self.species.reac_type[index] != 'hom_sci':
                        if not obj.ts_opt.shir == 1:  # last stage in optimize
                            opts_done = 0
                            obj.ts_opt.do_optimization()
                        if obj.ts_opt.shigh == -999:
                            logger.warning('Reaction {} ts_opt_shigh failure'.format(obj.instance_name))
                            fails = 1
                        if obj.do_vdW:
                            if obj.irc_prod_opt.shigh == -999:
                                logger.warning('High level vdW well optimization of {} failed'.format(obj.irc_prod.name))
                                obj.do_vdW = False
                            if not obj.irc_prod_opt.shigh == 1:
                                opts_done = 0
                                obj.irc_prod_opt.do_optimization()
                    for pr_opt in obj.prod_opt:
                        if not pr_opt.shir == 1:
                            opts_done = 0
                            pr_opt.do_optimization()
                        if pr_opt.shigh == -999:
                            logger.warning('Reaction {} pr_opt_shigh failure'.format(obj.instance_name))
                            fails = 1
                        continue
                            
                    if fails:
                        self.species.reac_ts_done[index] = -999
                    elif opts_done:
                        self.species.reac_ts_done[index] = 5

                elif self.species.reac_ts_done[index] == 5:
                    # Finilize the calculations
                    st_pt = obj.prod_opt[0].species
                    # kill reaction if higher than L2 threshold
                    if self.par['barrier_threshold_L2'] and self.par['high_level'] and 'hom_sci' not in obj.instance_name:
                        # check the barrier height again at L2 if requested
                        ts_energy = self.qc.get_qc_energy(f'{obj.instance_name}_high')[1]
                        ts_zpe = self.qc.get_qc_zpe(f'{obj.instance_name}_high')[1]
                        valid = (ts_energy + ts_zpe - self.species.energy - self.species.zpe) * constants.AUtoKCAL - self.par['barrier_threshold_L2']
                        if  valid > 0. :
                            logger.info(f'\t{obj.instance_name} is higher than the L2 threshold by {np.round(valid, 2)} kcal/mol, reaction is deleted.')
                            self.species.reac_ts_done[index] = -999
                            continue
                    # continue to PES search in case a new well was found
                    if self.par['pes']:
                        # verify if product is monomolecular, and if it is new
                        if len(obj.products) == 1:
                            st_pt = obj.prod_opt[0].species
                            chemid = st_pt.chemid
                            # if high level was requested, it is L2, otherwise L1
                            rel_en = (st_pt.energy + st_pt.zpe - self.species.energy - self.species.zpe) * constants.AUtoKCAL 
                            logger.info(f'\tProduct {obj.instance_name} energy is {np.round(rel_en, 2)} kcal/mol.')
                            if self.par['barrier_threshold_L2'] and self.par['high_level']:
                                new_barrier_threshold = None
                                new_barrier_threshold_L2 = self.par['barrier_threshold_L2'] - rel_en 
                            else:
                                new_barrier_threshold = self.par['barrier_threshold'] - rel_en 
                                new_barrier_threshold_L2 = None
                            dirwell = os.path.dirname(os.getcwd())
                            jobs = open(dirwell + '/chemids', 'r').read().split('\n')
                            jobs = [ji for ji in jobs]
                            if not str(chemid) in jobs:
                                # this well is new, add it to the jobs
                                while 1:
                                    try:
                                        # try to open the file and write to it
                                        logger.info(f'\tLaunching new KinBot as {chemid}')
                                        pes.write_input(self.inp, obj.products[0], new_barrier_threshold, new_barrier_threshold_L2, dirwell, self.par['me'])
                                        with open(dirwell + '/chemids', 'a') as f:
                                            f.write('{}\n'.format(chemid))
                                        break
                                    except IOError:
                                        # wait a second and try again
                                        time.sleep(1)
                                        pass

                    # check for wrong number of negative frequencies
                    neg_freq = 0
                    for st_pt in obj.products:
                        if len(st_pt.reduced_freqs):
                            if -1 * self.par['imagfreq_threshold'] <= st_pt.reduced_freqs[0] <= 0.:
                                logger.warning(f'Found negative frequency {st_pt.reduced_freqs[0]} cm-1 for a product of {obj.instance_name}. Flipped.')
                                st_pt.reduced_freqs[0] *= -1.
                            elif st_pt.reduced_freqs[0] < -1 * self.par['imagfreq_threshold']:
                                logger.warning(f'Found negative frequency {st_pt.reduced_freqs[0]} cm-1 for a product of {obj.instance_name}.')
                                self.species.reac_ts_done[index] = -999
                                neg_freq = 1
                    if any([fi < 0. for fi in obj.ts.reduced_freqs[1:]]):
                        logger.warning('Found more than one negative frequency for ' + obj.instance_name)
                        logger.warning(obj.ts.reduced_freqs)
                        self.species.reac_ts_done[index] = -999
                        neg_freq = 1
                        
                    if not neg_freq:
                        # the reaction search is finished
                        self.species.reac_ts_done[index] = -1  # this is the success code

                        # write a temporary pes input file
                        # remove old xval and im_extent files
                        if os.path.exists('{}_xval.txt'.format(self.species.chemid)):
                            os.remove('{}_xval.txt'.format(self.species.chemid))
                        if os.path.exists('{}_im_extent.txt'.format(self.species.chemid)):
                            os.remove('{}_im_extent.txt'.format(self.species.chemid))
                        postprocess.createPESViewerInput(self.species, self.qc, self.par)

                elif self.species.reac_ts_done[index] == -999:
                    if self.par['delete_intermediate_files'] == 1:
                        if not self.species.reac_obj[index].instance_name in deleted:
                            self.delete_files(self.species.reac_obj[index].instance_name)
                            deleted.append(self.species.reac_obj[index].instance_name)

            alldone = 1
            for index, instance in enumerate(self.species.reac_inst):
                if any(self.species.reac_ts_done[index] >= 0
                       for i in range(len(self.species.reac_inst))):
                    alldone = 1
                    break
                else:
                    alldone = 0

            # write a small summary while running
            with open('kinbot_monitor.out', 'w') as f_out:
                for index, instance in enumerate(self.species.reac_inst):
                    if self.species.reac_ts_done[index] == -1:
                        prodstring = []
                        for pp in self.species.reac_obj[index].products:
                            prodstring.append(str(pp.chemid))
                        f_out.write('{}\t{}\t{}\t{}\n'.format(self.species.reac_ts_done[index], 
                                                              self.species.reac_step[index], 
                                                              self.species.reac_obj[index].instance_name,
                                                              ' '.join(prodstring)))
                    else:
                        f_out.write('{}\t{}\t{}\n'.format(self.species.reac_ts_done[index], 
                                                          self.species.reac_step[index], 
                                                          self.species.reac_obj[index].instance_name))
            time.sleep(1)

        # Create molpro file for the BLS products
        for index, instance in enumerate(self.species.reac_inst):
            if self.species.reac_type[index] == 'barrierless_saddle':
                obj = self.species.reac_obj[index]
                if len(obj.products) == 2:
                    blsprodatom = np.append(obj.products[0].atom, obj.products[1].atom)
                    # vector connecting the centers
                    # assumed: the two structures are naturally separated at the end of the IRC
                    x0 = 0
                    y0 = 0
                    z0 = 0
                    for g in obj.products[0].geom:
                        x0 += g[0]
                        y0 += g[1]
                        z0 += g[1]
                    x0 /= obj.products[0].natom
                    y0 /= obj.products[0].natom
                    z0 /= obj.products[0].natom

                    x1 = 0
                    y1 = 0
                    z1 = 0
                    for g in obj.products[1].geom:
                        x1 += g[0]
                        y1 += g[1]
                        z1 += g[1]
                    x1 /= obj.products[1].natom
                    y1 /= obj.products[1].natom
                    z1 /= obj.products[1].natom

                    center_vec = [x1 - x0, y1 - y0, z1- z0]
                    blsprodgeom = np.append(obj.products[0].geom, obj.products[1].geom)
                    blsprodgeom = np.reshape(blsprodgeom, (-1, 3))
                    bps = list(zip(blsprodatom, blsprodgeom))
                    bps1 = [item for sublist in bps for item in sublist]
                    blsprodstruct = [item for sublist in bps1 for item in sublist]
                    blsprod = StationaryPoint('blsprod',
                                              self.par['charge'],
                                              self.par['mult'],
                                              structure=blsprodstruct)
                    blsprod.characterize()
                    molp = Molpro(blsprod, self.par)
                    molpname = '{}_prod'.format(obj.instance_name)
                    molp.create_molpro_input(bls=1, name=molpname, shift_vec=center_vec, natom1=len(obj.products[1].geom))
                    molp.create_molpro_submit(name=molpname)

        s = []
        for index, instance in enumerate(self.species.reac_inst):
            obj = self.species.reac_obj[index]
            # Write a summary on the combinatorial exploration
            if 'combinatorial' in obj.instance_name:
                s.append('NAME\t' + obj.instance_name)

                # Write the bonds that were broken and formed
                s.append('BROKEN_BONDS\t' + '\t'.join('[{}, {}]'.format(re[0], re[1]) for re in obj.reac))
                s.append('FORMED_BONDS\t' + '\t'.join('[{}, {}]'.format(pr[0], pr[1]) for pr in obj.prod))

                # Populate the ts_bond_lengths dict with the values
                # of this reaction

                if self.species.reac_ts_done[index] == -1:
                    for i in range(self.species.natom - 1):
                        for j in range(i + 1, self.species.natom):
                            if self.species.bond[i][j] != obj.irc_prod.bonds[0][i][j]:
                                if (self.species.bond[i][j] == 0 or obj.irc_prod.bonds[0][i][j] == 0):
                                    syms = []
                                    syms.append(self.species.atom[i])
                                    syms.append(self.species.atom[j])
                                    syms = ''.join(sorted(syms))
                                    dist = np.linalg.norm(obj.ts.geom[i] - obj.ts.geom[j])
                                    s.append('TS_BOND_LENGTHS\t{}\t{}'.format(syms, dist))
                # write the expected inchis

                s.append('EXPECTED_INCHIS\t' + '\t'.join(inchi for inchi in obj.prod_inchi))
                # get the inchis the reaction found
                if self.species.reac_ts_done[index] == -1:
                    inchis = obj.get_final_inchis()
                    s.append('FOUND_INCHIS\t' + '\t'.join(inchis))
                s.append('\n')
            with open('combinatorial.txt', 'w') as f:
                f.write('\n'.join(s) + '\n')

        logger.info('Reaction generation done!')

    def delete_files(self, name):
        # job names
        names = []

        names.append(name)
        names.append(name + '_high')
        names.append(name + '_IRC_F')
        names.append(name + '_IRC_R')
        names.append(name + '_IRC_F_prod')
        names.append(name + '_IRC_R_prod')
        extensions = ['chk', 'py', 'sbatch']

        for name in names:
            for ext in extensions:
                # delete file
                file = '.'.join([name, ext])
                try:
                    os.remove(file)
                except OSError:
                    pass

    def get_stereochemistry(self, atom, geom):
        inchi = cheminfo.create_inchi_from_geom(atom, geom)
        if '/t' in str(inchi):
            stereochem = inchi.split('/t')[1].split('/')[0]
        else:
            stereochem = ''
        return stereochem

    def equate_identical(self, frag):
        ''' Make identical fragments for a given reaction be exactly the same
        '''
        for ii in range(len(frag)):
            for jj in range(ii + 1, len(frag)):
                if frag[ii].chemid == frag[jj].chemid:
                    frag[jj] = frag[ii]
        return

    def equate_unique(self, fragments, frag_unique):
        '''Make identical fragments across reactions be exatly the same
        '''
        for fi, frag in enumerate(fragments):
            if len(frag_unique) == 0:
                frag_unique.append(frag)
            elif len(frag_unique) > 0:
                new = 1
                for fragu in frag_unique:  
                    if frag.chemid == fragu.chemid:
                        fragments[fi] = fragu
                        new = 0
                        break
                if new:
                    frag_unique.append(frag)
        return
