from __future__ import print_function
from shutil import copyfile
import sys
import json
import logging
import os
import datetime
import subprocess
import time
import pkg_resources
import math
import random
import numpy as np

from kinbot import constants
from kinbot import frequencies
from kinbot import license_message
from kinbot.uncertaintyAnalysis import UQ
from kinbot.parameters import Parameters


class MESS:
    """
    Class that read and writes MESS files for UQ analysis
    Code is the same as mess.py exceppt that it loops to create 'n' files
    and it randomizes the following parameters within the alloted UQ tolerance
       1. Stationary point energy (E+ZPE, +/- 0.5 kcal/mol)
       2. Barrier (E+ZPE, +/- 1.0 kcal/mol)
       3. Frequencies (cm-1 +/-20%)
    Ranges were chosen/based on the following paper:  Goldsmith, C. F. PCI, 2013, 177-185
    'n' = number of input files generated, default = X (will be implemented once code works)
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

    def write_header(self):
        """
        Create the header block for MESS
        """
        # Read the header template
        header_file = pkg_resources.resource_filename('tpl', 'mess_header.tpl')
        with open(header_file) as f:
            tpl = f.read()
        header = tpl.format(TemperatureList=' '.join([str(ti) for ti in self.par.par['TemperatureList']]),
                            PressureList=' '.join([str(pi) for pi in self.par.par['PressureList']]),
                            EnergyStepOverTemperature=self.par.par['EnergyStepOverTemperature'],
                            ExcessEnergyOverTemperature=self.par.par['ExcessEnergyOverTemperature'],
                            ModelEnergyLimit=self.par.par['ModelEnergyLimit'],
                            CalculationMethod=self.par.par['CalculationMethod'],
                            ChemicalEigenvalueMax=self.par.par['ChemicalEigenvalueMax'],
                            Reactant=self.well_names[self.species.chemid],
                            EnergyRelaxationFactor=self.par.par['EnergyRelaxationFactor'],
                            EnergyRelaxationPower=self.par.par['EnergyRelaxationPower'],
                            EnergyRelaxationExponentCutoff=self.par.par['EnergyRelaxationExponentCutoff'],
                            Epsilons=' '.join([str(ei) for ei in self.par.par['Epsilons']]),
                            Sigmas=' '.join([str(si) for si in self.par.par['Sigmas']]),
                            Masses=' '.join([str(mi) for mi in self.par.par['Masses']]))
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
                    #TER MOLECULAR
                    for st_pt in reaction.products:
                        if st_pt.chemid not in self.fragment_names:
                            self.fragment_names[st_pt.chemid] = 'fr_' + str(len(self.fragment_names)+1)
                    termol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                    if termol_name not in self.termolec_names:
                        self.termolec_names[termol_name] = 't_' + str(len(self.termolec_names)+1)

        #Barrierless short names
        try:
            for hs in self.species.homolytic_scissions.hss:
                if hs.status == -1:
                    if len(hs.products) == 1:
                        st_pt = hs.products[0]
                        if st_pt.chemid not in self.well_names:
                            self.well_names[st_pt.chemid] = 'w_' + str(len(self.well_names)+1)
                    elif len(hs.products) == 2:
                        for st_pt in hs.products:
                            if st_pt.chemid not in self.fragment_names:
                                self.fragment_names[st_pt.chemid] = 'fr_' + str(len(self.fragment_names)+1)
                    bimol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in hs.products]))
                    if bimol_name not in self.bimolec_names:
                        self.bimolec_names[bimol_name] = 'b_' + str(len(self.bimolec_names) + 1)
                    else:
                        #TER MOLECULAR
                        for st_pt in hs.products:
                            if st_pt.chemid not in self.fragment_names:
                                self.fragment_names[st_pt.chemid] = 'fr_' + str(len(self.fragment_names)+1)
                        termol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in hs.products]))
                        if termol_name not in self.termolec_names:
                            self.termolec_names[termol_name] = 't_' + str(len(self.termolec_names)+1)
        
        except:
            logging.info("No Homolytic Scission Reactions")

    def write_input(self, uq, uq_n, qc):
        """
        write the input for all the wells, bimolecular products and barriers
        both in a separate file, as well as in one large ME file
        """
        
        uq_obj = UQ()
        well_uq = self.par.par['well_uq']
        barrier_uq = self.par.par['barrier_uq']
        freq_uq = self.par.par['freq_uq']
        imagfreq_uq = self.par.par['imagfreq_uq']
        qc = qc

        if uq == 1:
            logging.info("Uncertainty Analysis is turned on, number of mess files being generated = {}".format(uq_n))
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
                for x in sorted(prod_list):
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
        while ( uq_iter < uq_n ):
            #set UQ factors for each i run
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

            if uq_iter >= 0:
                # energy and freq for each iteration of UQ
                # arrays reset at the start of each iteration to hold new values
                ts_e_iter = []
                ts_imagFreq_iter = []
                ts_freq_iter = []
                ts_rxnName_iter = []

                well_e_iter = []
                well_fr_iter = []

                prod_e_iter = []
                prod_fr_iter = []
                
                bimol_e_iter = []
                bimol_fr_iter = []

                termol_e_iter = []
                termol_fr_iter = []
                
                barrierless_e_iter = []
                barrierless_fr_iter = []

                well_uqVal = float(well_uq)
                freq_uqVal = float(freq_uq)
                imagfreq_uqVal = float(imagfreq_uq)
                barrier_uqVal = float(barrier_uq)
        
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
                for index, reaction in enumerate(self.species.reac_obj):
                    barrier_add = uq_obj.calc_energyUQ(barrier_uqVal)
                    ts_freqFactor = uq_obj.calc_freqUQ(freq_uqVal)
                    imagfreqFactor = uq_obj.calc_freqUQ(imagfreq_uqVal)
                    if reaction.instance_name in ts_all:
                        allTS[reaction.instance_name], ts_e, ts_freq, ts_imagFreq = self.write_barrier(reaction,
                                                                                                       uq, 
                                                                                                       uq_n,
                                                                                                       barrier_add,
                                                                                                       ts_freqFactor,
                                                                                                       imagfreqFactor,
                                                                                                       uq_iter,
                                                                                                       qc)
                        ts_e_iter.append(ts_e)
                        ts_imagFreq_iter.append(ts_imagFreq)
                        ts_freq_iter.append(ts_freq)

                    if reaction.instance_name in ts_unique:
                        lenProd = len(reaction.products)
                        ts_blocks[reaction.instance_name] = allTS[reaction.instance_name]
                        freqFactor = uq_obj.calc_freqUQ(freq_uqVal)
                        energyAdd = uq_obj.calc_energyUQ(well_uqVal)
                        if uq_iter == (uq_n-1):
                            ts_rxnName_iter.append(reaction.instance_name)
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
                        else:
                            #termol
                            termolec_ts_blocks[reaction.instance_name] = allTS[reaction.instance_name]
                            termol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                            termolec_blocks[termol_name] = self.write_termol([opt.species for opt in reaction.prod_opt], reaction, uq, uq_n, energyAdd, freqFactor, 0, uq_iter)
                            
                #Homolytic scission - barrierless reactions
                barrierless = {}
                if self.species.homolytic_scissions is not None:
                    for hs in self.species.homolytic_scissions.hss:
                        bar = 1
                        barrierless_freqFactor = uq_obj.calc_freqUQ(freq_uqVal)
                        barrierless_energyAdd = uq_obj.calc_energyUQ(well_uqVal)
                        new=1
                        if hs.status == -1:
                            prod_name = '_'.join(sorted([str(prod.chemid) for prod in hs.products]))
                            if prod_name in self.barrierless_names:
                                new=0
                            if new:
                                self.barrierless_names[prod_name] = prod_name
                            if new == 1 or uq_iter >= 1:
                                barrierless_blocks[prod_name], barrierless_e, barrierless_fr = self.write_barrierless([opt.species for opt in hs.prod_opt],
                                                                                                                      hs,
                                                                                                                      uq,
                                                                                                                      uq_n,
                                                                                                                      barrierless_energyAdd,
                                                                                                                      barrierless_freqFactor,
                                                                                                                      bar,
                                                                                                                      uq_iter)

            #uq_obj.norm_energy(all_wellE, 'well energy', all_tsRxnName, n) 
                            
            wells = ''
            for well in well_blocks:
                wells += well_blocks[well] + '\n!****************************************\n'
            bimols = ''
            for bimol in bimolec_blocks:
                bimols += bimolec_blocks[bimol] + '\n!****************************************\n'
            termols = ''
            for termol in termolec_blocks:
                termols += termolec_blocks[termol] + '\n!****************************************\n'
            tss = ''
            for ts in ts_blocks:
                tss += ts_blocks[ts] + '\n!****************************************\n'
            barrierless = ''
            for rxn in barrierless_blocks:
                barrierless += barrierless_blocks[rxn] + '\n!****************************************\n'

            dummy_template = pkg_resources.resource_filename('tpl', 'mess_dummy.tpl')
            with open(dummy_template) as f:
                dummy = f.read()
            dum = dummy.format(barrier='tsd', reactant=self.well_names[self.species.chemid], dummy='d1')


            mess_iter = "{0:04d}".format(uq_iter)
            f_out = open('me/mess_%s.inp' %mess_iter, 'w')
            f_out.write(header + '\n!****************************************\n')
            f_out.write(wells)
            f_out.write(bimols)
            f_out.write(tss)
            f_out.write(termols)
            f_out.write(barrierless)
            f_out.write('\n!****************************************\nEnd ! end kinetics\n')
            f_out.close()
       
            uq_iter = uq_iter + 1 
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
        #Create the dummy MESS block for ter-molecular products.
        #open the dummy template
        termol_file = pkg_resources.resource_filename('tpl', 'mess_termol.tpl')
        with open(termol_file) as f:
            tpl = f.read()
        fragment_file = pkg_resources.resource_filename('tpl', 'mess_fragment.tpl')
        with open(fragment_file) as f:
            fragment_tpl = f.read()
        hir_file = pkg_resources.resource_filename('tpl', 'mess_hinderedrotor.tpl')
        with open(hir_file) as f:
            rotor_tpl = f.read()
        atom_file = pkg_resources.resource_filename('tpl', 'mess_atom.tpl')
        with open(atom_file) as f:
            atom_tpl = f.read()

        fragments = ''
        termol = ''

        terPr_name = '_'.join(sorted([str(species.chemid) for species in species_list]))
        for species in species_list:
            if species.natom > 1:
                rotors = []
                if self.par.par['rotor_scan']:
                    for i, rot in enumerate(species.dihed):
                        group = ' '.join([str(pi+1) for pi in frequencies.partition(species, rot, species.natom)[0][1:]])
                        axis = '{} {}'.format(str(rot[1]+1), str(rot[2]+1))
                        rotorsymm = species.sigma_int[rot[1]][rot[2]]
                        nrotorpot = species.hir.nrotation // rotorsymm
                        ens = species.hir.hir_energies[i]
                        rotorpot = [(ei - ens[0])*constants.AUtoKCAL for ei in ens]
                        rotorpot = ' '.join(['{:.2f}'.format(ei) for ei in rotorpot[:species.hir.nrotation // rotorsymm]])
                        rotors.append(rotor_tpl.format(group=group,
                                                       axis=axis,
                                                       rotorsymm=rotorsymm,
                                                       nrotorpot=nrotorpot,
                                                       rotorpot=rotorpot))
                rotors = '\n'.join(rotors)

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

                    allFreqs=",".join(str(termolFreq) for termolFreq in termolArrayFreq)

                geom = ''
                for i, at in enumerate(species.atom):
                    if i > 0:
                        geom += '            '
                    x, y, z = species.geom[i]
                    geom += '! {} {:.6f} {:.6f} {:.6f}\n'.format(at, x, y, z)

                energy = ((species.energy + species.zpe) - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL

                if self.par.par['pes']:
                    name = 'fr_name_{}'.format(species.chemid)
                    name = '{' + name + '}'
                    energy = '{ground_energy}'
                else:
                    name = self.fragment_names[species.chemid] + ' ! ' + str(species.chemid)
                #molecule template
                fragments += fragment_tpl.format(chemid=name,
                                                 natom=species.natom,
                                                 geom=geom,
                                                 symm=float(species.sigma_ext) / float(species.nopt),
                                                 nfreq=len(species.reduced_freqs),
                                                 freq=freq,
                                                 hinderedrotor=rotors,
                                                 nelec=1,
                                                 charge=species.charge,
                                                 mult=species.mult)
                fragments += '\n'
            else:
                if self.par.par['pes']:
                    name = 'fr_name_{}'.format(species.chemid)
                    name = '{' + name + '}'
                else:
                    name = self.fragment_names[species.chemid] + ' ! ' + str(species.chemid)

                # atom template
                fragments += atom_tpl.format(chemid=name,
                                             element=species.atom[0],
                                             nelec=1,
                                             charge=species.charge,
                                             mult=species.mult)
                fragments += '\n'

        prod_name = self.termolec_names[terPr_name]

        if bar == 0:
            rxn_name = self.ts_names[reaction.instance_name]
            for species in species_list:
                if self.par.par['pes']:
                    name = 'fr_name_{}'.format(species.chemid)
                    name = 't_' + name
                else:
                    name = self.termolec_names[terPr_name] + ' ! ' + terPr_name
        else:
            rxn_name = 'termol_nobar_' + str(bar)
        # full termol template in progress need to figure out how to show data without code reading data
        #termol += tpl.format(product=prod_name,dummy=terPr_name, fragments=fragments, ground_energy=energy)
        termol += tpl.format(name=prod_name, product=terPr_name)
        f=open(terPr_name + '.mess', 'w')
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
 
        if bar == 0:
            bimol_file = pkg_resources.resource_filename('tpl', 'mess_bimol.tpl')
        elif bar == 1: 
            bimol_file = pkg_resources.resource_filename('tpl', 'mess_barrierless.tpl')

        with open(bimol_file) as f:
            tpl = f.read()
        fragment_file = pkg_resources.resource_filename('tpl', 'mess_fragment.tpl')
        with open(fragment_file) as f:
            fragment_tpl = f.read()
        hir_file = pkg_resources.resource_filename('tpl', 'mess_hinderedrotor.tpl')
        with open(hir_file) as f:
            rotor_tpl = f.read()
        atom_file = pkg_resources.resource_filename('tpl', 'mess_atom.tpl')
        with open(atom_file) as f:
            atom_tpl = f.read()
        
        fragments = ''
        for species in species_list:
            if species.natom > 1:
                rotors = []
                if self.par.par['rotor_scan']:
                    for i, rot in enumerate(species.dihed):
                        group = ' '.join([str(pi+1) for pi in frequencies.partition(species, rot, species.natom)[0][1:]])
                        axis = '{} {}'.format(str(rot[1]+1), str(rot[2]+1))
                        rotorsymm = species.sigma_int[rot[1]][rot[2]]
                        nrotorpot = species.hir.nrotation // rotorsymm
                        ens = species.hir.hir_energies[i]
                        rotorpot = [(ei - ens[0])*constants.AUtoKCAL for ei in ens]
                        rotorpot = ' '.join(['{:.2f}'.format(ei) for ei in rotorpot[:species.hir.nrotation // rotorsymm]])
                        rotors.append(rotor_tpl.format(group=group,
                                                       axis=axis,
                                                       rotorsymm=rotorsymm,
                                                       nrotorpot=nrotorpot,
                                                       rotorpot=rotorpot))
                rotors = '\n'.join(rotors)   
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
                     
                    allFreqs=",".join(str(bimolFreq) for bimolFreq in bimolArrayFreq)

                geom = ''
                for i, at in enumerate(species.atom):
                    if i > 0:
                        geom += '            '
                    x, y, z = species.geom[i]
                    geom += '{} {:.6f} {:.6f} {:.6f}\n'.format(at, x, y, z)

                energy = ((species.energy + species.zpe) - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL
                
                if self.par.par['pes']:
                    name = 'fr_name_{}'.format(species.chemid)
                    name = '{' + name + '}'
                    energy = '{ground_energy}'
                else:
                    name = self.fragment_names[species.chemid] + ' ! ' + str(species.chemid)
                #molecule template
                fragments += fragment_tpl.format(chemid=name,
                                                 natom=species.natom,
                                                 geom=geom,
                                                 symm=float(species.sigma_ext) / float(species.nopt),
                                                 nfreq=len(species.reduced_freqs),
                                                 freq=freq,
                                                 hinderedrotor=rotors,
                                                 nelec=1,
                                                 charge=species.charge,
                                                 mult=species.mult)
                fragments += '\n'
            else:
                if self.par.par['pes']:
                    name = 'fr_name_{}'.format(species.chemid)
                    name = '{' + name + '}'
                else:
                    name = self.fragment_names[species.chemid] + ' ! ' + str(species.chemid)
            
                # atom template
                fragments += atom_tpl.format(chemid=name,
                                             element=species.atom[0],
                                             nelec=1,
                                             charge=species.charge,
                                             mult=species.mult)
                fragments += '\n'

        pr_name = '_'.join(sorted([str(species.chemid) for species in species_list]))
        if self.par.par['pes']:
            name = '{name}'
            energy = '{ground_energy}'
        else:
            if bar == 0:
                name = self.bimolec_names[pr_name] + ' ! ' + pr_name
                 
            elif bar == 1:
                name = self.barrierless_names[pr_name] + '! barrierless'

            energy = (sum([sp.energy for sp in species_list]) + sum([sp.zpe for sp in species_list]) -
                      (self.species.energy + self.species.zpe)) * constants.AUtoKCAL

            logFile.write("Total energy for Bimol/Termol Product: {}\n".format(name))
            if uq_iter == 0:
                well_add = 0.0
                logFile.write("\tBimol/Termol well_add: {}\n".format(well_add))
                logFile.write("\t\tOriginal Energy = {}\n".format(energy))
            if uq_iter > 0:
                logFile.write("\tBimol/Termol well_add: {}\n".format(well_add))
                logFile.write("\t\tOriginal Energy = {}\n".format(energy))
                energy=energy + well_add
            logFile.write("\t\tUpdated energy = {}\n\n".format(energy))

        logFile.close()
        bimolArrayEnergy.append(energy)

        if bar == 0:
            bimol = tpl.format(chemids=name,
                               fragments=fragments,
                               ground_energy=energy)
            
        elif bar == 1:
            reac = self.well_names[self.species.chemid]
            index = self.barrierless_names.keys().index(pr_name)
            bimol = tpl.format(barrier='nobar_{}'.format(index),
                               reactant=reac,
                               prod=name,
                               dummy='',
                               fragments=fragments,
                               ground_energy=energy)
            
        if self.par.par['uq'] == 0:
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
        # open the templates
        well_file = pkg_resources.resource_filename('tpl', 'mess_well.tpl')
        with open(well_file) as f:
            tpl = f.read()
        hir_file = pkg_resources.resource_filename('tpl', 'mess_hinderedrotor.tpl')
        with open(hir_file) as f:
            rotor_tpl = f.read()

        rotors = []
        if self.par.par['rotor_scan']:
            for i, rot in enumerate(species.dihed):
                group = ' '.join([str(pi+1) for pi in frequencies.partition(species, rot, species.natom)[0][1:]])
                axis = '{} {}'.format(str(rot[1]+1), str(rot[2]+1))
                rotorsymm = species.sigma_int[rot[1]][rot[2]]
                nrotorpot = species.hir.nrotation // rotorsymm
                ens = species.hir.hir_energies[i]
                rotorpot = [(ei - ens[0]) * constants.AUtoKCAL for ei in ens]
                rotorpot = ' '.join(['{:.2f}'.format(ei) for ei in rotorpot[:species.hir.nrotation // rotorsymm]])
                rotors.append(rotor_tpl.format(group=group,
                                               axis=axis,
                                               rotorsymm=rotorsymm,
                                               nrotorpot=nrotorpot,
                                               rotorpot=rotorpot))
        rotors = '\n'.join(rotors)

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

        geom = ''
        for i, at in enumerate(species.atom):
            if i > 0:
                geom += '            '
            x, y, z = species.geom[i]
            geom += '{} {:.6f} {:.6f} {:.6f}\n'.format(at, x, y, z)

        if self.par.par['pes']:
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
        mess_well = tpl.format(chemid=name,
                               natom=species.natom,
                               geom=geom,
                               symm=float(species.sigma_ext) / float(species.nopt),
                               nfreq=len(species.reduced_freqs),
                               freq=freq,
                               hinderedrotor=rotors,
                               nelec=1,
                               charge=species.charge,
                               mult=species.mult,
                               zeroenergy=energy)

        if self.par.par['uq'] == 0:
            f = open(str(species.chemid) + '.mess', 'w')
        else:
            mess_iter = "{0:04d}".format(uq_iter)
            f = open(str(species.chemid) + '_' + str(mess_iter) + '.mess', 'w')
        f.write(mess_well)
        f.close()

        wellsArrayFreq.pop(0)
        wellsArrayEnergy.pop(0)
        e = wellsArrayEnergy
        fr = wellsArrayFreq       

        return mess_well, e, fr

    def write_barrier(self, reaction, uq, uq_n, barrier_add, freqFactor, imagfreqFactor, uq_iter, qc):
        """
        Create the block for a MESS barrier.
        """
        logFile = open('uq.log', 'a')

        # open the templates
        ts_file = pkg_resources.resource_filename('tpl', 'mess_ts.tpl')
        with open(ts_file) as f:
            tpl = f.read()
        hir_file = pkg_resources.resource_filename('tpl', 'mess_hinderedrotor.tpl')
        with open(hir_file) as f:
            rotor_tpl = f.read()
        tunn_file = pkg_resources.resource_filename('tpl', 'mess_tunneling.tpl')
        with open(tunn_file) as f:
            tun_tpl = f.read()
        
        rotors = []
        if self.par.par['rotor_scan']:
            for i, rot in enumerate(reaction.ts.dihed):
                group = ' '.join([str(pi+1) for pi in frequencies.partition(reaction.ts, rot, reaction.ts.natom)[0][1:]])
                axis = '{} {}'.format(str(rot[1]+1), str(rot[2]+1))
                rotorsymm = reaction.ts.sigma_int[rot[1]][rot[2]]
                nrotorpot = reaction.ts.hir.nrotation // rotorsymm
                ens = reaction.ts.hir.hir_energies[i]
                rotorpot = [(ei - ens[0]) * constants.AUtoKCAL for ei in ens]
                rotorpot = ' '.join(['{:.2f}'.format(ei) for ei in rotorpot[:reaction.ts.hir.nrotation // rotorsymm]])
                rotors.append(rotor_tpl.format(group=group,
                                               axis=axis,
                                               rotorsymm=rotorsymm,
                                               nrotorpot=nrotorpot,
                                               rotorpot=rotorpot))
        rotors = '\n'.join(rotors)

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

        geom = ''
        for i, at in enumerate(reaction.ts.atom):
            if i > 0:
                geom += '            '
            x, y, z = reaction.ts.geom[i]
            geom += '{} {:.6f} {:.6f} {:.6f}\n'.format(at, x, y, z)

        for index in range(len(self.species.reac_inst)):
            if self.species.reac_ts_done[index] == -1:
                ts = self.species.reac_obj[index].ts
                if self.species.reac_type[index] == 'R_Addition_MultipleBond' and not self.par.par['high_level']:
                    mp2_energy = qc.get_qc_energy(str(self.species.chemid) + '_well_mp2')[1]
                    mp2_zpe = qc.get_qc_zpe(str(self.species.chemid) + '_well_mp2')[1]
                    energy = (ts.energy + ts.zpe - mp2_energy - mp2_zpe) * constants.AUtoKCAL
                    prod_mp2_energy = 0
                    prod_mp2_zpe = 0
                    for opt in reaction.prod_opt:
                        energy = qc.get_qc_energy(str(opt.species.chemid) + 'well_mp2')[1]
                        zpe = qc.get_qc_zpe(str(opt.species.chemid) + 'well_mp2')[1]
                        prod_mp2_energy = prod_mp2_energy + energy
                        prod_zpe_energy = prod_mp2_zpe + zpe
                    energy2 = (ts.energy + ts.zpe - prod_mp2_energy - prod_zpe_energy) * constants.AUtoKCAL
                else:
                    energy = (ts.energy + ts.zpe - self.species.energy - self.species.zpe) * constants.AUtoKCAL
                    energy2 = (reaction.ts.energy + reaction.ts.zpe) - sum([(opt.species.energy + opt.species.zpe) for opt in reaction.prod_opt]) * constants.AUtoKCAL 

                barriers = [ energy, energy2, ]

                """
                barriers = [
                    ((reaction.ts.energy + reaction.ts.zpe) - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL,
                    ((reaction.ts.energy + reaction.ts.zpe) - sum([(opt.species.energy + opt.species.zpe) for opt in reaction.prod_opt])) * constants.AUtoKCAL,
                ]
                """
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
            tun = tun_tpl.format(cutoff=min(barriers),
                                 imfreq=imagfreq,
                                 welldepth1=barriers[0],
                                 welldepth2=barriers[1])


        if len(reaction.products) == 1:
            prod_name = self.well_names[reaction.products[0].chemid]
        elif len(reaction.products) == 2:
            long_name = '_'.join(sorted([str(pi.chemid) for pi in reaction.products]))
            prod_name = self.bimolec_names[long_name]
        else:
            long_name = '_'.join(sorted([str(pi.chemid) for pi in reaction.products]))
            prod_name = self.termolec_names[long_name]


        if self.par.par['pes']:
            name = '{name}'    
            chemid_reac = ''
            chemid_prod = ''
            long_rxn_name = ''
            energy = '{zeroenergy}'
        else:
            name = self.ts_names[reaction.instance_name]
            chemid_reac = self.well_names[self.species.chemid]
            chemid_prod = prod_name
            long_rxn_name = reaction.instance_name
            for index in range(len(self.species.reac_inst)):
                if self.species.reac_ts_done[index] == -1:
                    ts = self.species.reac_obj[index].ts
                    if self.species.reac_type[index] == 'R_Addition_MultipleBond' and not self.par.par['high_level']:
                        mp2_energy = qc.get_qc_energy(str(self.species.chemid) + '_well_mp2')[1]
                        mp2_zpe = qc.get_qc_zpe(str(self.species.chemid) + '_well_mp2')[1]
                        energy = (ts.energy + ts.zpe - mp2_energy - mp2_zpe) * constants.AUtoKCAL
                    else:
                        energy = (ts.energy + ts.zpe - self.species.energy - self.species.zpe) * constants.AUtoKCAL
            
            if uq_iter == 0:
                barrier_add = 0.0
                logFile.write("\tBarrier_add: {}\n".format(barrier_add))
                logFile.write("\t\tOriginal barrier energy: {}\n".format(energy))
            if uq_iter > 0:
                logFile.write("\tBarrier factor: {}\n".format(barrier_add))
                logFile.write("\t\tOriginal barrier energy: {}\n".format(energy))
                energy = energy + barrier_add
            logFile.write("\t\tUpdated barrier energy: {}\n\n".format(energy))
        logFile.close()

        barrierArrayEnergy.append(energy)

        mess_ts = tpl.format(rxn_name=name,
                             chemid_reac=chemid_reac,
                             chemid_prod=chemid_prod,
                             long_rxn_name=long_rxn_name,
                             natom=reaction.ts.natom,
                             geom=geom,
                             symm=float(reaction.ts.sigma_ext) / float(reaction.ts.nopt),
                             nfreq=len(reaction.ts.reduced_freqs) - 1,
                             freq=freq,
                             hinderedrotor=rotors,
                             tunneling=tun,
                             nelec=1,
                             charge=reaction.ts.charge,
                             mult=reaction.ts.mult,
                             zeroenergy=energy)

        if self.par.par['uq'] == 0:
            f = open(reaction.instance_name + '.mess', 'w')
        else:
            mess_iter = "{0:04d}".format(uq_iter)
            f = open(reaction.instance_name + '_' + str(mess_iter) + '.mess', 'w')
        f.write(mess_ts)
        f.close()

        barrierArrayImagfreq.pop(0)
        barrierArrayPosFreq.pop(0)
        barrierArrayEnergy.pop(0)
        e = barrierArrayEnergy
        frPos = barrierArrayPosFreq
        frNeg = barrierArrayImagfreq

        return mess_ts, e, frPos, frNeg

    def write_high_p_rates(all_rates, well_short_names, bimol_short_names):
        f_out = open('rates.out', 'w')
        if len(all_rates) > 0:
            f_out.write('temperatures\t\t')
            f_out.write('\t'.join([str(ti) for ti in all_rates[0].temperatures]))
            f_out.write('\n')

            for rates in all_rates:
                true_rxn = 1
                if rates.reactant in well_short_names.values():
                    for name in well_short_names:
                        if well_short_names[name] == rates.reactant:
                            reac_name = name
                elif rates.reactant in bimol_short_names.values():
                    for name in bimol_short_names:
                        if bimol_short_names[name] == rates.reactant:
                            reac_name = name
                else:
                    true_rxn = 0
                if rates.product in well_short_names.values():
                    for name in well_short_names:
                        if well_short_names[name] == rates.product:
                            prod_name = name
                elif rates.product in bimol_short_names.values():
                    for name in bimol_short_names:
                        if bimol_short_names[name] == rates.product:
                            prod_name = name
                else:
                    true_rxn = 0

                if true_rxn:
                    s = reac_name + '\t'
                    s += prod_name + '\t'
                    s += '\t'.join([str(ri) for ri in rates.rates[-1]])
                    s += '\n'
                    f_out.write(s)

        f_out.close()

    def run(self, uq_n):
        """
        write a pbs or slurm file for the me/all.inp mess input file
        submit the pbs/slurm file to the queue
        wait for the mess run to finish
        """
        
        # open the the header and the specific templates
        
        if self.par.par['queue_template'] == '':
            q_file = pkg_resources.resource_filename('tpl', self.par.par['queuing'] + '.tpl')
        else:
            q_file = self.par.par['queue_template']
        with open(q_file) as f:
            tpl_head = f.read()
        
        q_file = pkg_resources.resource_filename('tpl', self.par.par['queuing'] + '_mess_uq.tpl')
        with open(q_file) as f:
            tpl = f.read()
        queue_name = self.par.par['queuing']
        submitscript = 'run_mess' + constants.qext[self.par.par['queuing']]

        '''
           copy submit script so that new runs can be done from run_mess_n files
           need to check for jobs running
           if < x jobs running submit next batch
           runs through counter adding up to n mess jobs total
           need to optimize the number that can be submitted
        '''

        uq_iter = 0 #counter for jobs
        maxRunning = 1
        if self.par.par['uq'] == 1:
            maxRunning = self.par.par['uq_max_runs'] #max jobs running at once, can make this an input parameter at somepoint if neccessary
        job_counter = maxRunning
        previousLoop = 0 #number of running jobs on previous run so that jobs are not double counted as finishing
        pids = [] #list of job pids
        while(uq_iter < uq_n):
            if uq_n < job_counter:
                maxRunning = job_counter = uq_n
            while(uq_iter < job_counter):
                with open(submitscript, 'w') as f:
                    mess_iter = "{0:04d}".format(uq_iter)
                    if self.par.par['queue_template'] == '':
                        if self.par.par['queuing'] == 'pbs':
                            f.write((tpl_head).format(name='mess', ppn=self.par.par['ppn'], queue_name=self.par.par['queue_name'], dir='me'))
                            f.write((tpl).format(n=mess_iter))
                        elif self.par.par['queuing'] == 'slurm':
                            f.write((tpl_head).format(name='mess', ppn=self.par.par['ppn'], queue_name=self.par.par['queue_name'], dir='me'), slurm_feature='')
                            f.write((tpl).format(n=mess_iter))
                    else:
                        f.write(tpl_head)
                        f.write((tpl).format(n=mess_iter))
                         
                command = [constants.qsubmit[self.par.par['queuing']], submitscript ]
                process = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = process.communicate()
                out = out.decode()
                if self.par.par['queuing'] == 'pbs':
                    pid = out.split('\n')[0].split('.')[0]
                elif self.par.par['queuing'] == 'slurm':
                    pid = out.split('\n')[0].split()[-1]
                pids.append(pid)
                time.sleep(5)
                uq_iter = uq_iter + 1 

            currentJob = 0 
            runningJobs = 0
            exit = 0 #bool to exit loop
            while(exit == 0):
                devnull = open(os.devnull, 'w')
                currentLoop = 0 #number of jobs finished in current loop through pids[]
                while(currentJob < len(pids)):
                    pid_currentJob = pids[currentJob]
                    if self.par.par['queuing'] == 'pbs':
                        command = 'qstat -f | grep ' + '"Job Id: ' + pid_currentJob + '"' + ' > /dev/null'
                    elif self.par.par['queuing'] == 'slurm':
                        command = 'scontrol show job ' + pid_currentJob + ' | grep "JobId=' + pid_currentJob + '"' + ' > /dev/null'

                    stat = int(subprocess.call(command, shell=True, stdout=devnull, stderr=devnull))
                    if stat == 0:
                        #time.sleep(1)
                        currentJob = currentJob + 1
                        runningJobs = runningJobs + 1

                    if stat == 1:
                        currentLoop = currentLoop + 1
                        currentJob = currentJob + 1
                    if currentJob == len(pids) and runningJobs == maxRunning:
                        if runningJobs == (uq_n-1):
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
