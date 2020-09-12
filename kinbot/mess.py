from __future__ import print_function
import logging
import os
import subprocess
import time
import pkg_resources

from kinbot import constants
from kinbot import frequencies
from kinbot.uncertaintyAnalysis import UQ


class MESS:
    """
    Class that read and writes MESS files
    UQ analysis parameter (uq) can be used to generate 'n' number of mess input files
    with the following parameters randomized within the alloted UQ tolerance.
    UQ tolerance is set to the default values as follows
       1. Stationary point energy (E+ZPE, +/- 0.5 kcal/mol)
       2. Barrier (E+ZPE, +/- 1.0 kcal/mol)
       3. Frequencies (cm-1 +/- 20%)
    Default parameters were chosen/based on the following paper:  Goldsmith, C. F. PCI, 2013, 177-185
    New parameters can be set within the input json file with the following keywords
    See parameters.py file for more information.
    """

    def __init__(self, par, species):
        self.par = par
        self.species = species
        self.well_names = {}
        self.bimolec_names = {}
        self.fragment_names = {}
        self.ts_names = {}
        self.termolec_names = {}
        self. barrierless_names = {}
        # read all templates to create mess input
        with open(pkg_resources.resource_filename('tpl', 'mess_header.tpl')) as f:
            self.headertpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_dummy.tpl')) as f:
            self.dummytpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_termol.tpl')) as f:
            self.termoltpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_fragment.tpl')) as f:
            self.fragmenttpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_hinderedrotor.tpl')) as f:
            self.hinderedrotortpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_atom.tpl')) as f:
            self.atomtpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_tunneling.tpl')) as f:
            self.tunneltpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_well.tpl')) as f:
            self.welltpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_bimol.tpl')) as f:
            self.bimoltpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_barrierless.tpl')) as f:
            self.blbimoltpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_barrier.tpl')) as f:
            self.barriertpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_rrho.tpl')) as f:
            self.rrhotpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_core_rr.tpl')) as f:
            self.corerrtpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_pst.tpl')) as f:
            self.psttpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_variational.tpl')) as f:
            self.variationaltpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_2tst.tpl')) as f:
            self.twotstpl = f.read()

    def write_header(self):
        """
        Create the header block for MESS
        """
        # Read the header template
        header = self.headertpl.format(TemperatureList=' '.join([str(ti) for ti in self.par['TemperatureList']]),
                                       PressureList=' '.join([str(pi) for pi in self.par['PressureList']]),
                                       EnergyStepOverTemperature=self.par['EnergyStepOverTemperature'],
                                       ExcessEnergyOverTemperature=self.par['ExcessEnergyOverTemperature'],
                                       ModelEnergyLimit=self.par['ModelEnergyLimit'],
                                       CalculationMethod=self.par['CalculationMethod'],
                                       ChemicalEigenvalueMax=self.par['ChemicalEigenvalueMax'],
                                       Reactant=self.well_names[self.species.chemid],
                                       EnergyRelaxationFactor=self.par['EnergyRelaxationFactor'],
                                       EnergyRelaxationPower=self.par['EnergyRelaxationPower'],
                                       EnergyRelaxationExponentCutoff=self.par['EnergyRelaxationExponentCutoff'],
                                       e_coll=constants.epsilon[self.par['collider']],
                                       s_coll=constants.sigma[self.par['collider']],
                                       m_coll=constants.mass[self.par['collider']],
                                       e_well=self.par['epsilon'],
                                       s_well=self.par['sigma'],
                                       m_well=self.species.mass,
                                       )
        return header

    def create_short_names(self):
        """
        Create a short name for all the wells, all the bimolecular products and all the transition states
        """
        # add the initial well to the well names:
        self.well_names[self.species.chemid] = 'w_1'

        for index, reaction in enumerate(self.species.reac_obj):
            if self.species.reac_ts_done[index] == -1:
                self.ts_names[reaction.instance_name] = 'ts_' + str(len(self.ts_names) + 1)
                if len(reaction.products) == 1:
                    st_pt = reaction.products[0]
                    if st_pt.chemid not in self.well_names:
                        self.well_names[st_pt.chemid] = 'w_' + str(len(self.well_names) + 1)
                elif len(reaction.products) == 2:
                    for st_pt in reaction.products:
                        if st_pt.chemid not in self.fragment_names:
                            self.fragment_names[st_pt.chemid] = 'fr_' + str(len(self.fragment_names) + 1)
                    bimol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                    if bimol_name not in self.bimolec_names:
                        self.bimolec_names[bimol_name] = 'b_' + str(len(self.bimolec_names) + 1)
                else:
                    # TER MOLECULAR
                    for st_pt in reaction.products:
                        if st_pt.chemid not in self.fragment_names:
                            self.fragment_names[st_pt.chemid] = 'fr_' + str(len(self.fragment_names) + 1)
                    termol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                    if termol_name not in self.termolec_names:
                        self.termolec_names[termol_name] = 't_' + str(len(self.termolec_names) + 1)

        # Barrierless short names
        try:
            for hs in self.species.homolytic_scissions.hss:
                if hs.status == -1:
                    if len(hs.products) == 1:
                        st_pt = hs.products[0]
                        if st_pt.chemid not in self.well_names:
                            self.well_names[st_pt.chemid] = 'w_' + str(len(self.well_names) + 1)
                    elif len(hs.products) == 2:
                        for st_pt in hs.products:
                            if st_pt.chemid not in self.fragment_names:
                                self.fragment_names[st_pt.chemid] = 'fr_' + str(len(self.fragment_names) + 1)
                    bimol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in hs.products]))
                    if bimol_name not in self.bimolec_names:
                        self.bimolec_names[bimol_name] = 'b_' + str(len(self.bimolec_names) + 1)
                    else:
                        # TER MOLECULAR
                        for st_pt in hs.products:
                            if st_pt.chemid not in self.fragment_names:
                                self.fragment_names[st_pt.chemid] = 'fr_' + str(len(self.fragment_names) + 1)
                        termol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in hs.products]))
                        if termol_name not in self.termolec_names:
                            self.termolec_names[termol_name] = 't_' + str(len(self.termolec_names) + 1)

        except:
            logging.info("No Homolytic Scission Reactions")

    def write_input(self, uq, uq_n, qc):
        """
        write the input for all the wells, bimolecular products and barriers
        both in a separate file, as well as in one large ME file
        """

        uq_obj = UQ()
        well_uqVal = float(self.par['well_uq'])
        barrier_uqVal = float(self.par['barrier_uq'])
        freq_uqVal = float(self.par['freq_uq'])
        imagfreq_uqVal = float(self.par['imagfreq_uq'])

        if uq == 1:
            logging.info('Uncertainty Analysis is turned on,'
                         'number of mess files being generated = {}'.format(uq_n))
        else:
            logging.info("Uncertainty Analysis is turned off.")

        # create short names for all the species, bimolecular products and barriers
        self.create_short_names()
        header = self.write_header()

        # filter ts's with the same reactants and products:
        ts_unique = {}  # key: ts name, value: [prod_name, energy]
        ts_all = {}
        for index, reaction in enumerate(self.species.reac_obj):
            if self.species.reac_ts_done[index] == -1:
                prod_list = reaction.products
                rxnProds = []
                for x in prod_list:
                    rxnProds.append(x.chemid)
                rxnProds.sort()
                prod_name = '_'.join([str(pi) for pi in rxnProds])
                energy = reaction.ts.energy
                zpe = reaction.ts.zpe
                new = 1
                remove = []
                ts_all[reaction.instance_name] = [prod_name, energy + zpe]
                for ts in ts_unique:
                    if ts_unique[ts][0] == prod_name:
                        # check for the barrier with the lowest energy
                        if ts_unique[ts][1] > energy + zpe:
                            # remove the current barrier
                            remove.append(ts)
                        else:
                            new = 0
                for ts in remove:
                    ts_unique.pop(ts, None)
                if new:
                    ts_unique[reaction.instance_name] = [prod_name, energy + zpe]

        # write the mess input for the different blocks
        uq_iter = 0
        # list of energy variations through UQ
        ts_e_iter = []
        well_e_iter = []
        prod_e_iter = []
        bimol_e_iter = []
        termol_e_iter = []
        barrierless_e_iter = []
        
        for uq_iter in range(0, uq_n):
            with open('uq.log', 'a') as f:
                f.write("uq iteration {}".format(uq_iter))
            # set UQ factors for each i run
            if uq_iter == 0:
                well_add = barrier_add = 0.0
                freqFactor = imagfreqFactor = 1.0

            well_blocks = {}
            ts_blocks = {}
            bimolec_blocks = {}
            termolec_blocks = {}
            termolec_ts_blocks = {}
            barrierless_blocks = {}
            allTS = {}

            # energy and freq for each iteration of UQ
            # arrays reset at the start of each iteration to hold new values
            ts_imagFreq_iter = []
            ts_freq_iter = []
            ts_rxnName = []

            well_fr_iter = []
            well_name = []
       
            prod_fr_iter = []
            prod_names = []
            bimol_fr_iter = []
            bimol_names = []

            termol_fr_iter = []
            termolec_names = []

            barrierless_fr_iter = []
            barrierless_name = []

            well_energyAdd = uq_obj.calc_energyUQ(well_uqVal)
            well_freqFactor = uq_obj.calc_freqUQ(freq_uqVal)
            well_blocks[self.species.chemid], well_e, well_fr = self.write_well(self.species,
                                                                                uq,
                                                                                uq_n,
                                                                                well_energyAdd,
                                                                                well_freqFactor,
                                                                                uq_iter)
            well_e_iter.append(well_e)
            well_fr_iter.append(well_fr)
            well_name.append(self.species.chemid)
            count = 0
            for index, reaction in enumerate(self.species.reac_obj):
                barrier_adds = []
                for rxn in ts_all:
                    barrier_add = uq_obj.calc_energyUQ(barrier_uqVal)
                    barrier_adds.append(barrier_add)
                ts_freqFactor = uq_obj.calc_freqUQ(freq_uqVal)
                imagfreqFactor = uq_obj.calc_freqUQ(imagfreq_uqVal)
                if reaction.instance_name in ts_all:
                    allTS[reaction.instance_name], ts_e, ts_freq, ts_imagFreq = self.write_barrier(reaction,
                                                                                                   index,
                                                                                                   uq,
                                                                                                   uq_n,
                                                                                                   barrier_adds,
                                                                                                   ts_freqFactor,
                                                                                                   imagfreqFactor,
                                                                                                   uq_iter,
                                                                                                   qc,
                                                                                                   count)

                if reaction.instance_name in ts_unique:
                    lenProd = len(reaction.products)
                    ts_blocks[reaction.instance_name] = allTS[reaction.instance_name]
                    freqFactor = uq_obj.calc_freqUQ(freq_uqVal)
                    energyAdd = uq_obj.calc_energyUQ(well_uqVal)
                    if uq_iter == (uq_n - 1):
                        ts_e_iter.append(ts_e)
                        ts_imagFreq_iter.append(ts_imagFreq)
                        ts_freq_iter.append(ts_freq)
                        count = count + 1
                        ts_rxnName.append(reaction.instance_name)
                    if len(reaction.products) == 1:
                        st_pt = reaction.prod_opt[0].species
                        well_blocks[st_pt.chemid], prod_e, prod_fr = self.write_well(st_pt,
                                                                                     uq,
                                                                                     uq_n,
                                                                                     energyAdd,
                                                                                     freqFactor,
                                                                                     uq_iter)
                        prod_e_iter.append(prod_e)
                        prod_fr_iter.append(prod_fr)
                        prod_names.append(st_pt.chemid)
                    elif len(reaction.products) == 2:
                        bimol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                        bimolec_blocks[bimol_name], bimol_e, bimol_fr = self.write_bimol([opt.species for opt in reaction.prod_opt],
                                                                                         0,
                                                                                         uq,
                                                                                         uq_n,
                                                                                         lenProd,
                                                                                         energyAdd,
                                                                                         freqFactor,
                                                                                         uq_iter)
                        bimol_e_iter.append(bimol_e)
                        bimol_fr_iter.append(bimol_fr)
                        bimol_names.append(bimol_name)
                    else:
                        # termol
                        termolec_ts_blocks[reaction.instance_name] = allTS[reaction.instance_name]
                        termol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                        termolec_blocks[termol_name] = self.write_termol([opt.species for opt in reaction.prod_opt], reaction, uq, uq_n, energyAdd, freqFactor, 0, uq_iter)
                        termolec_names.append(termol_name)
            # Homolytic scission - barrierless reactions
            barrierless = {}
            if self.species.homolytic_scissions is not None:
                for hs in self.species.homolytic_scissions.hss:
                    bar = 1
                    barrierless_freqFactor = uq_obj.calc_freqUQ(freq_uqVal)
                    barrierless_energyAdd = uq_obj.calc_energyUQ(well_uqVal)
                    new = 1
                    if hs.status == -1:
                        hs_prod_name = '_'.join(sorted([str(prod.chemid) for prod in hs.products]))
                        if hs_prod_name not in self.bimolec_names and hs_prod_name not in self.termolec_names:
                            if hs_prod_name in self.barrierless_names:
                                new = 0
                            if new:
                                self.barrierless_names[hs_prod_name] = hs_prod_name
                            if new == 1 or uq_iter >= 1:
                                barrierless_blocks[hs_prod_name], barrierless_e, barrierless_fr = self.write_barrierless([opt.species for opt in hs.prod_opt],
                                                                                                                          hs,
                                                                                                                          uq,
                                                                                                                          uq_n,
                                                                                                                          barrierless_energyAdd,
                                                                                                                          barrierless_freqFactor,
                                                                                                                          bar,
                                                                                                                          uq_iter)
                            barrierless_e_iter.append(barrierless_e)
                            barrierless_fr_iter.append(barrierless_fr)
                            barrierless_name.append(hs_prod_name)
       
            wells = ''
            divider = '\n!****************************************\n'
            for well in well_blocks:
                wells += well_blocks[well] + divider
            bimols = ''
            for bimol in bimolec_blocks:
                bimols += bimolec_blocks[bimol] + divider
            termols = ''
            for termol in termolec_blocks:
                termols += termolec_blocks[termol] + divider
            tss = ''
            for ts in ts_blocks:
                tss += ts_blocks[ts] + divider
            barrierless = ''
            for rxn in barrierless_blocks:
                barrierless += barrierless_blocks[rxn] + divider

            dum = self.dummytpl.format(barrier='tsd', reactant=self.well_names[self.species.chemid], dummy='d1')

            mess_iter = "{0:04d}".format(uq_iter)
            with open('me/mess_%s.inp' % mess_iter, 'w') as f_out:
                f_out.write(header + divider + wells + bimols + tss + termols + barrierless + divider + 'End ! end kinetics\n')

        return 0

    def write_barrierless(self, species_list, reaction, uq, uq_n, energyAdd, freqFactor, bar, uq_iter):

        if len(reaction.products) == 2:
            lenProd = len(reaction.products)
            barrierless, barrierless_e, barrierless_fr = self.write_bimol(species_list,
                                                                          bar,
                                                                          uq,
                                                                          uq_n,
                                                                          lenProd,
                                                                          energyAdd,
                                                                          freqFactor,
                                                                          uq_iter)
        else:
            barrierless = self.write_termol(species_list, reaction, uq, uq_n, energyAdd, freqFactor, bar, uq_iter)
            barrierless_e, barrierless_fr = 0.0

        return barrierless, barrierless_e, barrierless_fr

    def write_termol(self, species_list, reaction, uq, uq_n, energyAdd, freqFactor, bar, uq_iter):
        # Create the dummy MESS block for ter-molecular products.
        fragments = ''
        termol = ''

        terPr_name = '_'.join(sorted([str(species.chemid) for species in species_list]))
        for species in species_list:
            if species.natom > 1:
                freq = ''
                termolArrayFreq = []
                termolArrayEnergy = []
                logFile = open('uq.log', 'a')

                logFile.write("Bimol species: {}\n".format(species.chemid))
                termolArrayFreq.append(species.chemid)
                termolArrayEnergy.append(species.chemid)
                for i, fr in enumerate(species.reduced_freqs):
                    if uq_iter == 0:
                        freqFactor = 1.0
                        logFile.write("\ttermol posFreq factor: {}\n".format(freqFactor))
                        logFile.write("\t\tOriginal first frequency: {}\n".format(fr))
                        termolArrayFreq.append(fr)
                        if i == 0:
                            logFile.write("\t\tUpdated first frequency: {}.\n\t\t\tFrequency should be unchanged.\n".format(fr))
                            freq += '! {:.4f}'.format(fr)
                        elif i > 0 and i % 3 == 0:
                            freq += '\n !           {:.4f}'.format(fr)
                        else:
                            freq += '!    {:.4f}'.format(fr)
                    elif uq_iter > 0:
                        if i == 0:
                            logFile.write("\ttermol posFreq factor: {}\n".format(freqFactor))
                            logFile.write("\t\tOriginal first frequency: {}\n".format(fr))
                        fr = fr * freqFactor
                        termolArrayFreq.append(fr)
                        if i == 0:
                            logFile.write("\t\tUpdated first frequency: {}\n".format(fr))
                            freq += '! {:.4f}'.format(fr)
                        elif i > 0 and i % 3 == 0:
                            freq += ' \n !            {:.4f}'.format(fr)
                        else:
                            freq += '!    {:.4f}'.format(fr)
                    else:
                        logging.error('uq_n is negative')

                    allFreqs = ",".join(str(termolFreq) for termolFreq in termolArrayFreq)

                energy = ((species.energy + species.zpe) - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL

                if self.par['pes']:
                    name = 'fr_name_{}'.format(species.chemid)
                    name = '{' + name + '}'
                    energy = '{ground_energy}'
                else:
                    name = self.fragment_names[species.chemid] + ' ! ' + str(species.chemid)
                # molecule template
                fragments += self.fragmenttpl.format(chemid=name,
                                                     natom=species.natom,
                                                     geom=self.make_geom(species),
                                                     symm=float(species.sigma_ext) / float(species.nopt),
                                                     nfreq=len(species.reduced_freqs),
                                                     freq=freq,
                                                     hinderedrotor=self.make_rotors(species),
                                                     nelec=1,
                                                     mult=species.mult)
            else:
                if self.par['pes']:
                    name = 'fr_name_{}'.format(species.chemid)
                    name = '{' + name + '}'
                else:
                    name = self.fragment_names[species.chemid] + ' ! ' + str(species.chemid)

                # atom template
                fragments += self.atomtpl.format(chemid=name,
                                                 element=species.atom[0],
                                                 nelec=1,
                                                 mult=species.mult)

        prod_name = self.termolec_names[terPr_name]

        if bar == 0:
            rxn_name = self.ts_names[reaction.instance_name]
            for species in species_list:
                if self.par['pes']:
                    name = 'fr_name_{}'.format(species.chemid)
                    name = 't_' + name
                else:
                    name = self.termolec_names[terPr_name] + ' ! ' + terPr_name
        else:
            rxn_name = 'termol_nobar_' + str(bar)
        # full termol template in progress need to figure out how to show data without code reading data
        # termol += tpl.format(product=prod_name,dummy=terPr_name, fragments=fragments, ground_energy=energy)
        termol += self.termoltpl.format(name=prod_name, product=terPr_name)
        if uq == 0:
            f = open(terPr_name + '.mess', 'w')
        else:
            mess_iter = "{0:04d}".format(uq_iter)
            f = open(terPr_name + '_' + mess_iter + '.mess', 'w')

        f.write(termol)
        f.close()

        return termol

    def write_bimol(self, species_list, bar, uq, uq_n, lenProd, well_add, freqFactor, uq_iter):
        """
        Create the block for MESS for a bimolecular product.
        well0: reactant on this PES (zero-energy reference)
        uq_n = number of uncertainty runs
        """
        # open the templates
        logFile = open('uq.log', 'a')

        fragments = ''
        for species in species_list:
            if species.natom > 1:
                freq = ''

                bimolArrayFreq = []
                bimolArrayEnergy = []
                logFile.write("Bimol species: {}\n".format(species.chemid))
                bimolArrayFreq.append(species.chemid)
                bimolArrayEnergy.append(species.chemid)
                for i, fr in enumerate(species.reduced_freqs):
                    if uq_iter == 0:
                        freqFactor = 1.0
                        logFile.write("\tBimol posFreq factor: {}\n".format(freqFactor))
                        logFile.write("\t\tOriginal first frequency: {}\n".format(fr))
                        bimolArrayFreq.append(fr)
                        if i == 0:
                            logFile.write("\t\tUpdated first frequency: {}.\n\t\t\tFrequency should be unchanged.\n".format(fr))
                            freq += '{:.4f}'.format(fr)
                        elif i > 0 and i % 3 == 0:
                            freq += '\n            {:.4f}'.format(fr)
                        else:
                            freq += '    {:.4f}'.format(fr)
                    elif uq_iter > 0:
                        if i == 0:
                            logFile.write("\tBimol posFreq factor: {}\n".format(freqFactor))
                            logFile.write("\t\tOriginal first frequency: {}\n".format(fr))
                        fr = fr * freqFactor
                        bimolArrayFreq.append(fr)
                        if i == 0:
                            logFile.write("\t\tUpdated first frequency: {}\n".format(fr))
                            freq += '{:.4f}'.format(fr)
                        elif i > 0 and i % 3 == 0:
                            freq += '\n            {:.4f}'.format(fr)
                        else:
                            freq += '    {:.4f}'.format(fr)
                    else:
                        logging.error('uq_n is negative')

                    allFreqs = ",".join(str(bimolFreq) for bimolFreq in bimolArrayFreq)

                energy = ((species.energy + species.zpe) - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL

                if self.par['pes']:
                    name = 'fr_name_{}'.format(species.chemid)
                    name = '{' + name + '}'
                    energy = '{ground_energy}'
                else:
                    name = self.fragment_names[species.chemid] + ' ! ' + str(species.chemid)
                # molecule template
                fragments += self.fragmenttpl.format(chemid=name,
                                                    natom=species.natom,
                                                    geom=self.make_geom(species),
                                                    symm=float(species.sigma_ext) / float(species.nopt),
                                                    nfreq=len(species.reduced_freqs),
                                                    freq=freq,
                                                    hinderedrotor=self.make_rotors(species),
                                                    nelec=1,
                                                    mult=species.mult)
            else:
                if self.par['pes']:
                    name = 'fr_name_{}'.format(species.chemid)
                    name = '{' + name + '}'
                else:
                    name = self.fragment_names[species.chemid] + ' ! ' + str(species.chemid)

                fragments += self.atomtpl.format(chemid=name,
                                                 element=species.atom[0],
                                                 nelec=1,
                                                 mult=species.mult)

        pr_name = '_'.join(sorted([str(species.chemid) for species in species_list]))
        if self.par['pes']:
            name = '{name}'
            energy = '{ground_energy}'
        else:
            if bar == 0:
                name = self.bimolec_names[pr_name] + ' ! ' + pr_name
            elif bar == 1:
                name = self.barrierless_names[pr_name] + '! barrierless'

            energy = (sum([sp.energy for sp in species_list]) + sum([sp.zpe for sp in species_list]) - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL

            logFile.write("Total energy for Bimol/Termol Product: {}\n".format(name))
            if uq_iter == 0:
                well_add = 0.0
                logFile.write("\tBimol/Termol well_add: {}\n".format(well_add))
                logFile.write("\t\tOriginal Energy = {}\n".format(energy))
            if uq_iter > 0:
                logFile.write("\tBimol/Termol well_add: {}\n".format(well_add))
                logFile.write("\t\tOriginal Energy = {}\n".format(energy))
                energy = energy + well_add
            logFile.write("\t\tUpdated energy = {}\n\n".format(energy))

        logFile.close()
        bimolArrayEnergy.append(energy)

        if bar == 0:
            bimol = self.bimoltpl.format(chemids=name,
                                         fragments=fragments,
                                         ground_energy=energy)

        elif bar == 1:
            index = self.barrierless_names.keys().index(pr_name)
            bimol = self.blbimoltpl.format(barrier='nobar_{}'.format(index),
                                           reactant=self.well_names[self.species.chemid],
                                           prod=name,
                                           dummy='',
                                           fragments=fragments,
                                           ground_energy=energy)

        if self.par['uq'] == 0:
            f = open(pr_name + '.mess', 'w')
        else:
            mess_iter = "{0:04d}".format(uq_iter)
            f = open(pr_name + '_' + str(mess_iter) + '.mess', 'w')
        f.write(bimol)
        f.close()

        strBimolArray = [str(i) for i in bimolArrayEnergy]
        joiner = ","
        strBimols = joiner.join(strBimolArray)

        energyVals = bimolArrayEnergy
        freqVals = bimolArrayFreq
        e = energyVals
        freqVals.pop(0)
        energyVals.pop(0)
        e = energyVals
        fr = freqVals

        return bimol, e, fr

    def write_well(self, species, uq, uq_n, well_add, freqFactor, uq_iter):
        """
        Create the block for MESS for a well.
        well0: reactant on this PES (zero-energy reference)
        """
        logFile = open('uq.log', 'a')

        freq = ''

        wellsArrayEnergy = []
        wellsArrayFreq = []
        wellsArrayEnergy.append(species.chemid)
        wellsArrayFreq.append(species.chemid)
        if uq == 1:
            logFile.write("Species: {}\n".format(species.chemid))
            logFile.write("\tWell freqFactor: {}\n".format(freqFactor))
        for i, fr in enumerate(species.reduced_freqs):
            if i == 0:
                logFile.write("\t\tOriginal first frequency: {}\n".format(fr))
            if uq_iter == 0:
                freqFactor = 1.0
                wellsArrayFreq.append(fr)
                if i == 0:
                    logFile.write("\t\tUpdated first frequency: {}\n".format(fr))
                    freq += '{:.4f}'.format(fr)
                elif i > 0 and i % 3 == 0:
                    freq += '\n            {:.4f}'.format(fr)
                else:
                    freq += '    {:.4f}'.format(fr)
            elif uq_iter > 0:
                fr = fr * freqFactor
                wellsArrayFreq.append(fr)
                if i == 0:
                    logFile.write("\t\tUpdated first frequency: {}\n".format(fr))
                    freq += '{:.4f}'.format(fr)
                elif i > 0 and i % 3 == 0:
                    freq += '\n            {:.4f}'.format(fr)
                else:
                    freq += '    {:.4f}'.format(fr)
            else:
                logging.error('uq_n is negative')

            allWellFreqs = ",".join(str(wellFreq) for wellFreq in wellsArrayFreq)

        if self.par['pes']:
            name = '{name}'
            energy = '{zeroenergy}'
        else:
            name = self.well_names[species.chemid] + ' ! ' + str(species.chemid)
            energy = ((species.energy + species.zpe) - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL
            if uq_iter == 0:
                well_add = 0.0
                logFile.write("\tWell_add: {}\n".format(well_add))
                logFile.write("\t\tOriginal energy = {}\n".format(energy))
            if uq_iter > 0:
                logFile.write("\tWell_factor: {}\n".format(well_add))
                logFile.write("\t\tOriginal energy = {}\n".format(energy))
                energy = energy + well_add
            logFile.write("\t\tUpdated energy = {}\n\n".format(energy))
        logFile.close()
        wellsArrayEnergy.append(energy)
        mess_well = self.welltpl.format(chemid=name,
                                        natom=species.natom,
                                        geom=self.make_geom(species),
                                        symm=float(species.sigma_ext) / float(species.nopt),
                                        nfreq=len(species.reduced_freqs),
                                        freq=freq,
                                        hinderedrotor=self.make_rotors(species),
                                        nelec=1,
                                        mult=species.mult,
                                        zeroenergy=energy)

        if self.par['uq'] == 0:
            with open(str(species.chemid) + '.mess', 'w') as f:
                f.write(mess_well)
        else:
            mess_iter = "{0:04d}".format(uq_iter)
            with open(str(species.chemid) + '_' + str(mess_iter) + '.mess', 'w') as f:
                f.write(mess_well)

        wellsArrayFreq.pop(0)
        wellsArrayEnergy.pop(0)
        e = wellsArrayEnergy
        fr = wellsArrayFreq

        return mess_well, e, fr

    def write_barrier(self, reaction, index, uq, uq_n, barrier_adds, freqFactor, imagfreqFactor, uq_iter, qc, count):
        """
        Create the block for a MESS barrier.
        """
        logFile = open('uq.log', 'a')

        freq = ''

        barrierArrayEnergy = []
        barrierArrayPosFreq = []
        barrierArrayImagfreq = []
        barrierArrayEnergy.append(reaction.instance_name)
        barrierArrayPosFreq.append(reaction.instance_name)
        barrierArrayImagfreq.append(reaction.instance_name)
        logFile.write("Reaction: {}\n".format(reaction.instance_name))
        logFile.write("\tBarrier freq factor: {}\n".format(freqFactor))

        for i, fr in enumerate(reaction.ts.reduced_freqs[1:]):
            if i == 0:
                logFile.write("\t\tOriginal first frequency: {}\n".format(fr))
            if uq_iter == 0:
                barrierArrayPosFreq.append(fr)
                if i == 0:
                    logFile.write("\t\tUpdated first frequency: {}.\n\t\t\t\tFrequency should be unchanged\n".format(fr))
                    freq += '{:.4f}'.format(fr)
                elif i > 0 and i % 3 == 0:
                    freq += '\n            {:.4f}'.format(fr)
                else:
                    freq += '    {:.4f}'.format(fr)
            elif uq_iter > 0:
                fr = fr * freqFactor
                if i == 0:
                    logFile.write("\t\tUpdated first frequency: {}\n".format(fr))
                barrierArrayPosFreq.append(fr)
                if i == 0:
                    freq += '{:.4f}'.format(fr)
                elif i > 0 and i % 3 == 0:
                    freq += '\n            {:.4f}'.format(fr)
                else:
                    freq += '    {:.4f}'.format(fr)
            else:
                logging.error('uq_iter not recognized')

            allBarrierFreqs = ",".join(str(barrierFreq) for barrierFreq in barrierArrayPosFreq)

        # get left-right barrier
        species_zeroenergy = (self.species.energy + self.species.zpe) * constants.AUtoKCAL
        if self.species.reac_ts_done[index] == -1:
            ts_zeroenergy = (reaction.ts.energy + reaction.ts.zpe) * constants.AUtoKCAL

            if not self.par['high_level'] and reaction.mp2 == 1:
                jobname = '{}_well_mp2'.format(str(self.species.chemid))
                well_zeroenergy = self.get_zeroenergy(jobname, qc)
            elif not self.par['high_level'] and self.species.reac_type[index] == 'barrierless_saddle':
                jobname = '{}_well_bls'.format(str(self.species.chemid))
                well_zeroenergy = self.get_zeroenergy(jobname, qc)
            else:
                well_zeroenergy = species_zeroenergy 
            left_zeroenergy = ts_zeroenergy - well_zeroenergy 

            prod_zeroenergy = 0
            for opt in reaction.prod_opt:
                prod_zeroenergy += (opt.species.energy + opt.species.zpe) * constants.AUtoKCAL
            right_zeroenergy = left_zeroenergy - (prod_zeroenergy - well_zeroenergy)

            barriers = [left_zeroenergy, right_zeroenergy]
    
        # write tunneling block
        if any([bi < 0 for bi in barriers]):
            tun = ''
        else:
            imagfreq = -reaction.ts.reduced_freqs[0]
            imagfreqFactor = imagfreqFactor
            logFile.write("\tTS imagfreqFactor: {}\n".format(imagfreqFactor))
            logFile.write("\t\tOriginal imaginary frequency: {}\n".format(imagfreq))
            if uq_iter == 0:
                logFile.write("\t\tUpdated imaginary frequency: {}.\n\t\t\tFrequency should be unchanged.\n".format(imagfreq))
                barrierArrayImagfreq.append(imagfreq)
            elif uq_iter > 0:
                imagfreq = imagfreq * imagfreqFactor
                logFile.write("\t\tUpdated imaginary frequency: {}\n".format(imagfreq))
                barrierArrayImagfreq.append(imagfreq)
            tun = self.tunneltpl.format(cutoff=min(barriers),
                                        imfreq=imagfreq,
                                        welldepth1=barriers[0],
                                        welldepth2=barriers[1])

        # name the product
        if len(reaction.products) == 1:
            prod_name = self.well_names[reaction.products[0].chemid]
        elif len(reaction.products) == 2:
            long_name = '_'.join(sorted([str(pi.chemid) for pi in reaction.products]))
            prod_name = self.bimolec_names[long_name]
        else:
            long_name = '_'.join(sorted([str(pi.chemid) for pi in reaction.products]))
            prod_name = self.termolec_names[long_name]

        if self.par['pes']:
            name = '{name}'
            chemid_reac = ''
            chemid_prod = ''
            long_rxn_name = ''
            zeroenergy = '{zeroenergy}'
        else:
            name = self.ts_names[reaction.instance_name]
            chemid_reac = self.well_names[self.species.chemid]
            chemid_prod = prod_name
            long_rxn_name = reaction.instance_name
            if uq_iter == 0:
                barrier_add = 0.0
            else:
                barrier_add = barrier_adds[count]
            logFile.write("\tBarrier factor: {}\n".format(barrier_adds[count]))
            logFile.write("\t\tOriginal barrier energy: {}\n".format(energies[count]))
            energy = left_barrier + barrier_add
            logFile.write("\t\tUpdated barrier energy: {}\n\n".format(energy))
        logFile.close()
    
        # TODO working here
        if self.species.reac_type[index] == 'barrierless_saddle':

            for opt in reaction.prod_opt:
                geom1=opt.species.geom
            outerts = self.psttpl.format(natom1=reaction.prod_opt[0].species.natom,
                                         geom1=self.make_geom(reaction.prod_opt[0].species),
                                         natom2=reaction.prod_opt[1].species.natom,
                                         geom2=self.make_geom(reaction.prod_opt[1].species),
                                         symm='xxx',
                                         prefact='prefactor',
                                         exponent=3.)
            twotst = self.twotstpl.format(outerts=outerts)
            corerr = self.corerrtpl.format(symm=float(reaction.ts.sigma_ext) / float(reaction.ts.nopt))
            rrho = self.rrhotpl.format(natom=reaction.ts.natom,
                                       geom=self.make_geom(reaction.ts),
                                       core=corerr,
                                       nfreq=len(reaction.ts.reduced_freqs)-1,
                                       freq=freq,
                                       rotor=self.make_rotors(self.species.reac_obj[index].ts, norot=self.ts_names[reaction.instance_name]),
                                       tunneling='',
                                       nelec=1,
                                       mult=reaction.ts.mult,
                                       zeroenergy=energy)
            variational = self.variationaltpl.format(twotst=twotst,
                                                     variationalmodel=rrho,
                                                     tunneling=tun)
            mess_barrier = self.barriertpl.format(rxn_name=name,
                                                  chemid_reac=chemid_reac,
                                                  chemid_prod=chemid_prod,
                                                  long_rxn_name=long_rxn_name,
                                                  model=variational)
        else:
            corerr = self.corerrtpl.format(symm=float(reaction.ts.sigma_ext) / float(reaction.ts.nopt))
            rrho = self.rrhotpl.format(natom=reaction.ts.natom,
                                       geom=self.make_geom(reaction.ts),
                                       core=corerr,
                                       nfreq=len(reaction.ts.reduced_freqs)-1,
                                       freq=freq,
                                       rotor=self.make_rotors(self.species.reac_obj[index].ts, norot=self.ts_names[reaction.instance_name]),
                                       tunneling=tun,
                                       nelec=1,
                                       mult=reaction.ts.mult,
                                       zeroenergy=energy)
            mess_barrier = self.barriertpl.format(rxn_name=name,
                                                  chemid_reac=chemid_reac,
                                                  chemid_prod=chemid_prod,
                                                  long_rxn_name=long_rxn_name,
                                                  model=rrho)

        if self.par['uq'] == 0:
            with open(reaction.instance_name + '.mess', 'w') as f:
                f.write(mess_barrier)
        else:
            with open('{}_{:04d}.mess'.format(reaction.instance_name, uq_iter), 'w') as f:
                f.write(mess_barrier)

        barrierArrayImagfreq.pop(0)
        barrierArrayPosFreq.pop(0)
        barrierArrayEnergy.pop(0)
        e = barrierArrayEnergy
        frPos = barrierArrayPosFreq
        frNeg = barrierArrayImagfreq

        return mess_barrier, energy, frPos, frNeg

    def run(self, uq_n):
        """
        write a pbs or slurm file for the me/all.inp mess input file
        submit the pbs/slurm file to the queue
        wait for the mess run to finish
        """

        # open the the header and the specific templates

        if self.par['queue_template'] == '':
            q_file = pkg_resources.resource_filename('tpl', self.par['queuing'] + '.tpl')
        else:
            q_file = self.par['queue_template']
        with open(q_file) as f:
            tpl_head = f.read()

        q_file = pkg_resources.resource_filename('tpl', self.par['queuing'] + '_mess_uq.tpl')
        with open(q_file) as f:
            tpl = f.read()
        # queue_name = self.par['queuing']
        submitscript = 'run_mess' + constants.qext[self.par['queuing']]

        '''
           copy submit script so that new runs can be done from run_mess_n files
           need to check for jobs running
           if < x jobs running submit next batch
           runs through counter adding up to n mess jobs total
           need to optimize the number that can be submitted
        '''

        uq_iter = 0  # counter for jobs
        maxRunning = 1
        if self.par['uq'] == 1:
            # change to same as the other max jobs parameter?
            maxRunning = self.par['uq_max_runs']  # max jobs running at once
        job_counter = maxRunning
        previousLoop = 0  # number of running jobs on previous run so that jobs are not double counted as finishing
        pids = []  # list of job pids
        while(uq_iter < uq_n):
            if uq_n < job_counter:
                maxRunning = job_counter = uq_n
            while(uq_iter < job_counter):
                with open(submitscript, 'w') as f:
                    mess_iter = "{0:04d}".format(uq_iter)
                    if self.par['queue_template'] == '':
                        if self.par['queuing'] == 'pbs':
                            f.write((tpl_head).format(name='mess', ppn=self.par['ppn'], queue_name=self.par['queue_name'], dir='me'))
                            f.write((tpl).format(n=mess_iter))
                        elif self.par['queuing'] == 'slurm':
                            f.write((tpl_head).format(name='mess', ppn=self.par['ppn'], queue_name=self.par['queue_name'], dir='me'), slurm_feature='')
                            f.write((tpl).format(n=mess_iter))
                    else:
                        f.write(tpl_head)
                        f.write((tpl).format(n=mess_iter))

                command = [constants.qsubmit[self.par['queuing']], submitscript]
                process = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = process.communicate()
                out = out.decode()
                if self.par['queuing'] == 'pbs':
                    pid = out.split('\n')[0].split('.')[0]
                elif self.par['queuing'] == 'slurm':
                    pid = out.split('\n')[0].split()[-1]
                pids.append(pid)
                time.sleep(5)
                uq_iter = uq_iter + 1

            currentJob = 0
            runningJobs = 0
            exit = 0  # bool to exit loop
            while(exit == 0):
                devnull = open(os.devnull, 'w')
                currentLoop = 0  # number of jobs finished in current loop through pids[]
                while(currentJob < len(pids)):
                    pid_currentJob = pids[currentJob]
                    if self.par['queuing'] == 'pbs':
                        command = 'qstat -f | grep ' + '"Job Id: ' + pid_currentJob + '"' + ' > /dev/null'
                    elif self.par['queuing'] == 'slurm':
                        command = 'scontrol show job ' + pid_currentJob + ' | grep "JobId=' + pid_currentJob + '"' + ' > /dev/null'

                    stat = int(subprocess.call(command, shell=True, stdout=devnull, stderr=devnull))
                    if stat == 0:
                        # time.sleep(1)
                        currentJob = currentJob + 1
                        runningJobs = runningJobs + 1

                    if stat == 1:
                        currentLoop = currentLoop + 1
                        currentJob = currentJob + 1
                    if currentJob == len(pids) and runningJobs == maxRunning:
                        if runningJobs == (uq_n - 1):
                            exit = 1
                        currentJob = 0
                        runningJobs = 0
                        currentLoop = 0

                    if((currentLoop > previousLoop) and (currentJob == len(pids))):
                        job_counter = job_counter + currentLoop - previousLoop
                        previousLoop = currentLoop
                        currentLoop = 0
                        runningJobs = 0
                        currentJob = 0
                        exit = 1

            if job_counter > uq_n:
                job_counter = uq_n

        return 0

    def make_geom(self, species):
        geom = ''
        for i, at in enumerate(species.atom):
            if i > 0:
                geom += '            '
            x, y, z = species.geom[i]
            geom += '! {} {:.6f} {:.6f} {:.6f}\n'.format(at, x, y, z)
        return geom

    def make_rotorpot(self, species, i, rot):
        rotorsymm = self.rotorsymm(species, rot)
        ens = species.hir.hir_energies[i]
        rotorpot = [(ei - ens[0]) * constants.AUtoKCAL for ei in ens]
        rotorpot = ' '.join(['{:.2f}'.format(ei) for ei in rotorpot[:species.hir.nrotation // rotorsymm]])
        return rotorpot

    def rotorsymm(self, species, rot):
        return species.sigma_int[rot[1]][rot[2]]

    def nrotorpot(self, species, rot): 
        rotorsymm = self.rotorsymm(species, rot)
        return species.hir.nrotation // rotorsymm

    def make_rotors(self, species, norot=None):
        rotors = []
        if self.par['rotor_scan']:
            for i, rot in enumerate(species.dihed):
                if norot is not None:
                    if frequencies.skip_rotor(norot, rot) == 1:
                        continue
                rotors.append(self.hinderedrotortpl.format(group=' '.join([str(pi + 1) for pi in frequencies.partition(species, rot, species.natom)[0][1:]]),
                                                           axis='{} {}'.format(str(rot[1] + 1), str(rot[2] + 1)),
                                                           rotorsymm=self.rotorsymm(species, rot),
                                                           nrotorpot=self.nrotorpot(species, rot),
                                                           rotorpot=self.make_rotorpot(species, i, rot)))
        rotors = '\n'.join(rotors)
        return rotors

    def get_zeroenergy(self, jobname, qc):
        energy = qc.get_qc_energy(jobname)[1]
        zpe = qc.get_qc_zpe(jobname)[1]
        return (energy + zpe) * constants.AUtoKCAL
