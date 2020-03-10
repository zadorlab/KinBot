##################################################
##                                               ##
## This file is part of the KinBot code v2.0     ##
##                                               ##
## The contents are covered by the terms of the  ##
## BSD 3-clause license included in the LICENSE  ##
## file, found at the root.                      ##
##                                               ##
## Copyright 2018 National Technology &          ##
## Engineering Solutions of Sandia, LLC (NTESS). ##
## Under the terms of Contract DE-NA0003525 with ##
## NTESS, the U.S. Government retains certain    ##
## rights to this software.                      ##
##                                               ##
## Authors:                                      ##
##   Judit Zador                                 ##
##   Ruben Van de Vijver                         ##
##                                               ##
###################################################
from __future__ import division
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
                self.ts_names[reaction.instance_name] = 'ts_' + str(len(self.ts_names)+1)
                if len(reaction.products) == 1:
                    st_pt = reaction.products[0]
                    if st_pt.chemid not in self.well_names:
                        self.well_names[st_pt.chemid] = 'w_' + str(len(self.well_names)+1)
                else:
                    for st_pt in reaction.products:
                        if st_pt.chemid not in self.fragment_names:
                            self.fragment_names[st_pt.chemid] = 'fr_' + str(len(self.fragment_names)+1)
                    bimol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                    if bimol_name not in self.bimolec_names:
                        self.bimolec_names[bimol_name] = 'b_' + str(len(self.bimolec_names)+1)

    def write_input(self,uq,n):
        """
        write the input for all the wells, bimolecular products and barriers
        both in a separate file, as well as in one large ME file
        """

        uq=uq
        n=n
        stPt_uq=self.par.par['stPt_uq']
        barrier_uq=self.par.par['barrier_uq']
        posFreq_uq=self.par.par['posFreq_uq']
        negFreq_uq=self.par.par['negFreq_uq']

        # create short names for all the species, bimolecular products and barriers
        self.create_short_names()
        header = self.write_header()
        # filter ts's with the same reactants and products:
        ts_unique = {}  # key: ts name, value: [prod_name, energy]
        ts_all = {}
        for index, reaction in enumerate(self.species.reac_obj):
            if self.species.reac_ts_done[index] == -1:
                prod_name = '_'.join([str(pi.chemid) for pi in reaction.products])
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
        logFile=open('uq.log','w')
        if uq == 1: 
            logFile.write("UQ analysis on, n = {}".format(n))
        elif uq == 0:
            logFile.write("UQ analysis off, n = {}".format(n))
        else:
            logFile.write("UQ analysis parameter is invalid")
        i=0
        if i == 0:
            # actual E/Fr values
            all_tsE=[]
            all_tsFrPos=[]
            all_tsFrNeg=[]
            all_tsRxnName=[]
            all_bimolE=[]
            all_bimolFr=[] 
            all_wellE=[]
            all_wellFr=[]
            all_prodE=[]
            all_prodFr=[]
        while ( i<n ):
            #set UQ factors for each i run
            if i == 0:
                stPt_factor=barrier_factor=posFreq_factor=negFreq_factor=0

            logFile=open('uq.log','a')
            logFile.write("\nUQ round: {}\n".format(i))
            if i == 0:
                logFile.write("\tInitial Run, all values should remain unchanged\n") 
            if i > 0 :
                logFile.write("\tStPt_factor:\t{}\n".format(stPt_factor))
                logFile.write("\tBarrier_factor:\t{}\n".format(barrier_factor))
                logFile.write("\tPosFreq_factor:\t{}\n".format(posFreq_factor))
                logFile.write("\tNegFreq_factor:\t{}\n\n".format(negFreq_factor))
            logFile.close()
    
            well_blocks = {}
            ts_blocks = {}
            bimolec_blocks = {}
            allTS= {}
            if i >= 0: 
                ts_e_i=[]
                ts_frNeg_i=[]
                ts_frPos_i=[]
                ts_rxnName_i=[]

                well_e_i=[]
                well_frPos_i=[]

                prod_e_i=[]
                prod_frPos_i=[]
                
                bimol_e_i=[]
                bimol_frPos_i=[]
                
                stPt_uqVal=float(stPt_uq)
                posFreq_uqVal=float(posFreq_uq)
                negFreq_uqVal=float(negFreq_uq)
                barrier_uqVal=float(barrier_uq)
     
                well_stPt_factor=random.uniform(-stPt_uqVal,stPt_uqVal)
                well_posFreqPercent=random.uniform(-posFreq_uqVal,posFreq_uqVal)
                well_posFreq_factor=(100+well_posFreqPercent)/100
                well_blocks[self.species.chemid], well_e, well_fr = self.write_well(self.species, i, uq, well_stPt_factor, well_posFreq_factor,n)
                well_e_i.append(well_e)
                well_frPos_i.append(well_fr)
                for index, reaction in enumerate(self.species.reac_obj):
                    barrier_factor=random.uniform(-barrier_uqVal,barrier_uqVal)
                    negFreqPercent=random.uniform(-negFreq_uqVal,negFreq_uqVal)
                    negFreq_factor=(100+negFreqPercent)/100
                    ts_posFreqPercent=random.uniform(-posFreq_uqVal,posFreq_uqVal)
                    ts_posFreq_factor=(100+ts_posFreqPercent)/100
             
                    if reaction.instance_name in ts_all:
                        allTS[reaction.instance_name], ts_e, ts_frPos, ts_frNeg = self.write_barrier(reaction, i, uq, barrier_factor, ts_posFreq_factor, negFreq_factor,n)
                        ts_e_i.append(ts_e)
                        ts_frNeg_i.append(ts_frNeg)
                        ts_frPos_i.append(ts_frPos)
                    if reaction.instance_name in ts_unique:
                        ts_blocks[reaction.instance_name] = allTS[reaction.instance_name]
                        if i == (n-1):
                            ts_rxnName_i.append(reaction.instance_name)
                        if len(reaction.products) == 1:
                            prod_stPt_factor=random.uniform(-stPt_uqVal,stPt_uqVal)
                            prod_posFreqPercent=random.uniform(-posFreq_uqVal,posFreq_uqVal)
                            prod_posFreq_factor=(100+prod_posFreqPercent)/100
                            st_pt = reaction.prod_opt[0].species
                            well_blocks[st_pt.chemid], prod_e, prod_fr = self.write_well(st_pt, i, uq, prod_stPt_factor, prod_posFreq_factor,n)
                            prod_e_i.append(prod_e)
                            prod_frPos_i.append(prod_fr)
                        else:
                            bimol_stPt_factor=random.uniform(-stPt_uqVal,stPt_uqVal)
                            bimol_posFreqPercent=random.uniform(-posFreq_uqVal,posFreq_uqVal)
                            bimol_posFreq_factor=(100+bimol_posFreqPercent)/100
                            bimol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                            bimolec_blocks[bimol_name], bimol_e, bimol_fr = self.write_bimol([opt.species for opt in reaction.prod_opt], i, uq, bimol_stPt_factor, bimol_posFreq_factor,n)
                            bimol_e_i.append(bimol_e)
                            bimol_frPos_i.append(bimol_fr)

                # Arrays for all UQ values
                # Normalization of values for postprocessing
                all_tsE.append(ts_e_i)
                all_tsFrPos.append(ts_frPos_i)
                all_tsFrNeg.append(ts_frNeg_i)
                all_tsRxnName.append(ts_rxnName_i)
                all_wellE.append(well_e_i)
                all_wellFr.append(well_frPos_i)
                all_prodE.append(prod_e_i)
                all_prodFr.append(prod_frPos_i)
                all_bimolE.append(bimol_e_i)
                all_bimolFr.append(bimol_frPos_i)
                
                allTSPosFr=[] 
                allWellFr=[]
                allProdFr=[] 
                allBimolFr=[] 
                
                # TO DO - ADD ALL VALUES TO ARRAY - only final itteration is being added currently
                # Array above has all of the iterations
                """
                if i == (n-1):
                    sets=len(all_tsE)
                    nts=len(all_tsE[0])
                    nwell=len(all_wellE[0])
                    nprods=len(all_prodE[0])
                    nbimols=len(all_bimolE[0])

                    all_tsE_t=[] #transposed array
                    all_wellE_t=[]
                    all_prodE_t=[]
                    all_bimolE_t=[]
                    i=0
                    while i < nts:
                        j=0
                        tsE_t=[]
                        while j < sets:
                            tsE_t.append(all_tsE[j][i])
                            j=j+1
                        all_tsE_t.append(tsE_t)
                        i=i+1
                    i=0
                    while i < nwell:
                        j=0
                        wellE_t=[]
                        while j < sets:
                            wellE_t.append(all_wellE[j][i])
                            j=j+1
                        all_wellE_t.append(wellE_t)
                        i=i+1
                    i=0
                    while i < nprods:
                        j=0
                        prodE_t=[]
                        while j < sets:
                            prodE_t.append(all_prodE[j][i])
                            j=j+1
                        all_prodE_t.append(prodE_t)
                        i=i+1
                    i=0
                    while i < nbimols:
                        j=0
                        bimolE_t=[]
                        while j < sets:
                            bimolE_t.append(all_bimolE[j][i])
                            j=j+1
                        all_bimolE_t.append(bimolE_t)
                        i=i+1

                    #self.norm(all_tsE_t, "tsE", 0, all_tsRxnName, len(all_tsE_t))
                    #self.norm(all_wellE_t, "wellE", 0, all_tsRxnName, len(all_wellE_t))
                    #self.norm(all_bimolE_t, "bimolE", 0, all_tsRxnName, len(all_bimolE_t))
                    #self.norm(all_prodE_t, "prodE", 0, all_tsRxnName, len(all_prodE_t))
                    
                    #negative freq for TSs
                    all_tsFrNeg_t=[]
                    i=0
                    while i < nts:
                        j=0
                        tsFrNeg_t=[]
                        while j < sets:
                            tsFrNeg_t.append(all_tsFrNeg[j][i])
                            j=j+1
                        all_tsFrNeg_t.append(tsFrNeg_t)
                        i=i+1
                    self.norm(all_tsFrNeg_t, "tsFrNeg", 1, all_tsRxnName, len(all_tsFrNeg_t))
             
                    #positive freqs
                    k=len(all_tsFrPos[0])
                    tsPosFr=[]
                    for i, x in enumerate(all_tsFrPos):
                        j=0
                        while j < k:
                            posfr=all_tsFrPos[i][j]
                            tsPosFr.append(posfr)
                            j=j+1
                        allTSPosFr.append(tsPosFr)
                        tsPosFr=[]
                
                    fr_i=[]
                    frSet_i=[]
                    m=len(all_wellFr[0])
                    #all_wellFr = [set][species][fr]
                    #goal = [species][fr][set]
                    k=0
                    for i, x in enumerate(all_wellFr): #sets
                        for j, y in enumerate(all_wellFr[i]): #species
                            while k < m: #freqs
                                posfr=all_wellFr[k][i][j]
                                fr_i.append(posfr) #all fr for set
                                k=k+1
                            frSet_i.append(fr_i)
                            fr_i=[]
                        all_wellFr.append(frSet_i)
                        frSet_i=[]

                    
                    for i, well in enumerate(all_wellFr):
                        for j, fr in enumerate(all_wellFr[i]):
 
                    k=len(all_bimolFr[0])
                    bimolFr=[]
                    for i, x in enumerate(all_bimolFr):
                        j=0
                        while j < k:
                            posfr=all_bimolFr[i][j]
                            bimolFr.append(posfr)
                            j=j+1
                        allBimolFr.append(bimolFr)
                        bimolFr=[]

                    k=len(all_prodFr[0])
                    prodPosFr=[]
                    for i, x in enumerate(all_prodFr):
                        j=0
                        while j < k:
                            posfr=all_prodFr[i][j]
                            prodPosFr.append(posfr)
                            j=j+1
                        allProdFr.append(prodPosFr)
                        prodPosFr=[]

                    #self.norm(allTSPosFr, "TS Pos Fr", 1, all_tsRxnName, len(allTSPosFr))
                    #self.norm(allWellFr, "Well Pos Fr", 2, all_tsRxnName, len(allWellFr))
                    #self.norm(allBimolFr, "Bimol Pos Fr", 1, all_tsRxnName, len(allBimolFr))
                    #self.norm(allProdFr, "Prod Pos Fr", 1, all_tsRxnName, len(allProdFr))
            """        
            wells = ''
            for well in well_blocks:
                wells += well_blocks[well] + '\n!****************************************\n'
            bimols = ''
            for bimol in bimolec_blocks:
                bimols += bimolec_blocks[bimol] + '\n!****************************************\n'
            tss = ''
            for ts in ts_blocks:
                tss += ts_blocks[ts] + '\n!****************************************\n'

            dummy_template = pkg_resources.resource_filename('tpl', 'mess_dummy.tpl')
            with open(dummy_template) as f:
                dummy = f.read()
            dum = dummy.format(barrier='tsd', reactant=self.well_names[self.species.chemid], dummy='d1')
            num=str(i)
            f_out = open('me/mess_%s.inp' %num, 'w')
            f_out.write(header + '\n!****************************************\n')
            f_out.write(wells)
            f_out.write(bimols)
            f_out.write(tss)
            f_out.write(dum)
            f_out.write('\n!****************************************\nEnd ! end kinetics\n')
            f_out.close()
       
            i=i+1 
        return 0

    def write_bimol(self, species_list, n_uq, uq, stPt_factor, posFreq_factor,n):
        """
        Create the block for MESS for a bimolecular product.
        well0: reactant on this PES (zero-energy reference)
        n_uq = number of uncertainty runs
        """
        # open the templates
        logFile=open('uq.log', 'a')

        bimol_file = pkg_resources.resource_filename('tpl', 'mess_bimol.tpl')
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
        
        n=n
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

                #calculate error (+/- 20%)
                # ! TO DO ! CHANGE +/- 20% to mult/div 20% !
                bimolArrayFreq=[]
                bimolArrayEnergy=[]
                posFreq_factor=posFreq_factor
                logFile.write("Bimol species: {}\n".format(species.chemid))
                bimolArrayFreq.append(species.chemid)
                bimolArrayEnergy.append(species.chemid)
                for i, fr in enumerate(species.reduced_freqs):
                    if n_uq == 0:
                        posFreq_factor=1
                        if i == 0:
                            logFile.write("\tBimol posFreq factor: {}\n".format(posFreq_factor))
                            logFile.write("\t\tOriginal first frequency: {}\n".format(fr)) 
                        bimolArrayFreq.append(fr)
                        if i == 0:
                            logFile.write("\t\tUpdated first frequency: {}.\n\t\t\tFrequency should be unchanged.\n".format(fr))
                            freq += '{:.4f}'.format(fr)
                        elif i > 0 and i % 3 == 0:
                            freq += '\n            {:.4f}'.format(fr)
                        else:
                            freq += '    {:.4f}'.format(fr)
                    elif n_uq > 0:
                        if i == 0:
                            logFile.write("\tBimol posFreq factor: {}\n".format(posFreq_factor))
                            logFile.write("\t\tOriginal first frequency: {}\n".format(fr))
                        fr = fr*posFreq_factor
                        bimolArrayFreq.append(fr)
                        if i == 0:
                            logFile.write("\t\tUpdated first frequency: {}\n".format(fr))
                            freq += '{:.4f}'.format(fr)
                        elif i > 0 and i % 3 == 0:
                            freq += '\n            {:.4f}'.format(fr)
                        else:
                            freq += '    {:.4f}'.format(fr)
                    else:
                        logging.error('n_uq is negative')
                     
                    allFreqs=",".join(str(bimolFreq) for bimolFreq in bimolArrayFreq )

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
            name = self.bimolec_names[pr_name] + ' ! ' + pr_name
            logFile.write("Total energy for Bimol Product: {}\n".format(name))
            energy = (sum([sp.energy for sp in species_list]) + sum([sp.zpe for sp in species_list]) -
                      (self.species.energy + self.species.zpe)) * constants.AUtoKCAL
            if n_uq == 0:
                stPt_factor=0
                logFile.write("\tBimol stPt_factor: {}\n".format(stPt_factor))
                logFile.write("\t\tOriginal Energy = {}\n".format(energy))
            if n_uq > 0:
                stPt_factor=stPt_factor
                logFile.write("\tBimol stPt_factor: {}\n".format(stPt_factor))
                logFile.write("\t\tOriginal Energy = {}\n".format(energy))
                energy=energy+stPt_factor
            logFile.write("\t\tUpdated energy = {}\n\n".format(energy))

        logFile.close()
        bimolArrayEnergy.append(energy) 
        bimol = tpl.format(chemids=name,
                           fragments=fragments,
                           ground_energy=energy)
        if self.par.par['uq'] == 0:
            f = open(pr_name + '.mess', 'w')
        else:
            f = open(pr_name + '_uq_' + str(n_uq) + '.mess', 'w')
        f.write(bimol)
        f.close()
       
        strBimolArray = [str(i) for i in bimolArrayEnergy]
        joiner=","
        strBimols=joiner.join(strBimolArray)
       

        energyVals=bimolArrayEnergy
        freqVals=bimolArrayFreq
        e=energyVals
        freqVals.pop(0)
        energyVals.pop(0)
        e=energyVals
        fr=freqVals        

        return bimol, e, fr

    def write_well(self, species, n_uq, uq, stPt_factor, posFreq_factor,n):
        """
        Create the block for MESS for a well.
        well0: reactant on this PES (zero-energy reference)
        """
        logFile=open('uq.log', 'a')
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
                rotorpot = [(ei - ens[0])*constants.AUtoKCAL for ei in ens]
                rotorpot = ' '.join(['{:.2f}'.format(ei) for ei in rotorpot[:species.hir.nrotation // rotorsymm]])
                rotors.append(rotor_tpl.format(group=group,
                                               axis=axis,
                                               rotorsymm=rotorsymm,
                                               nrotorpot=nrotorpot,
                                               rotorpot=rotorpot))
        rotors = '\n'.join(rotors)

        freq = ''
 
        wellsArrayEnergy=[]
        wellsArrayFreq=[]
        wellsArrayEnergy.append(species.chemid)
        wellsArrayFreq.append(species.chemid)
        if uq == 1:
            posFreq_factor=posFreq_factor
            logFile.write("Species: {}\n".format(species.chemid))
            logFile.write("\tWell posFreq_factor: {}\n".format(posFreq_factor))
        for i, fr in enumerate(species.reduced_freqs):
            if i == 0:
                logFile.write("\t\tOriginal first frequency: {}\n".format(fr))
            if n_uq == 0:
                posFreq_factor=1
                wellsArrayFreq.append(fr)
                if i == 0:
                    logFile.write("\t\tUpdated first frequency: {}\n".format(fr))
                    freq += '{:.4f}'.format(fr)
                elif i > 0 and i % 3 == 0:
                    freq += '\n            {:.4f}'.format(fr)
                else:
                    freq += '    {:.4f}'.format(fr)
            elif n_uq > 0:
                fr = fr*posFreq_factor
                wellsArrayFreq.append(fr)
                if i == 0:
                    logFile.write("\t\tUpdated first frequency: {}\n".format(fr))
                    freq += '{:.4f}'.format(fr)
                elif i > 0 and i % 3 == 0:
                    freq += '\n            {:.4f}'.format(fr)
                else:
                    freq += '    {:.4f}'.format(fr)
            else:
                logging.error('n_uq is negative')
            
            allWellFreqs=",".join(str(wellFreq) for wellFreq in wellsArrayFreq)

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
            if n_uq == 0:
                stPt_factor=0
                logFile.write("\tWell stPt factor: {}\n".format(stPt_factor))
                logFile.write("\t\tOriginal energy = {}\n".format(energy))
            if n_uq > 0:
                stPt_factor=stPt_factor
                logFile.write("\tWell stPt factor: {}\n".format(stPt_factor))
                logFile.write("\t\tOriginal energy = {}\n".format(energy))
                energy=energy+stPt_factor
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
            f = open(str(species.chemid) + '_uq_' + str(n_uq) + '.mess', 'w')
        f.write(mess_well)
        f.close()

        wellsArrayFreq.pop(0)
        wellsArrayEnergy.pop(0)
        e=wellsArrayEnergy
        fr=wellsArrayFreq       

        return mess_well, e, fr

    def write_barrier(self, reaction, n_uq, uq, barrier_factor, posFreq_factor, negFreq_factor,n):
        """
        Create the block for a MESS barrier.
        """
        logFile=open('uq.log', 'a')

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
                rotorpot = [(ei - ens[0])*constants.AUtoKCAL for ei in ens]
                rotorpot = ' '.join(['{:.2f}'.format(ei) for ei in rotorpot[:reaction.ts.hir.nrotation // rotorsymm]])
                rotors.append(rotor_tpl.format(group=group,
                                               axis=axis,
                                               rotorsymm=rotorsymm,
                                               nrotorpot=nrotorpot,
                                               rotorpot=rotorpot))
        rotors = '\n'.join(rotors)

        freq = ''

        barrierArrayEnergy=[]
        barrierArrayPosFreq=[]
        barrierArrayNegFreq=[]
        barrierArrayEnergy.append(reaction.instance_name)
        barrierArrayPosFreq.append(reaction.instance_name)
        barrierArrayNegFreq.append(reaction.instance_name)

        posFreq_factor=posFreq_factor
        logFile.write("Reaction: {}\n".format(reaction.instance_name))
        logFile.write("\tBarrier posFreq factor: {}\n".format(posFreq_factor))

        for i, fr in enumerate(reaction.ts.reduced_freqs[1:]):
            if i == 0:
                logFile.write("\t\tOriginal first frequency: {}\n".format(fr))
            if n_uq == 0:
                barrierArrayPosFreq.append(fr)
                if i == 0:
                    logFile.write("\t\tUpdated first frequency: {}.\n\t\t\t\tFrequency should be unchanged\n".format(fr))
                    freq += '{:.4f}'.format(fr)
                elif i > 0 and i % 3 == 0:
                    freq += '\n            {:.4f}'.format(fr)
                else:
                    freq += '    {:.4f}'.format(fr)
            elif n_uq > 0:
                fr = fr*posFreq_factor
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
                logging.error('n_uq is negative')

            allBarrierFreqs=",".join(str(barrierFreq) for barrierFreq in barrierArrayPosFreq)

        geom = ''
        for i, at in enumerate(reaction.ts.atom):
            if i > 0:
                geom += '            '
            x, y, z = reaction.ts.geom[i]
            geom += '{} {:.6f} {:.6f} {:.6f}\n'.format(at, x, y, z)

        barriers = [
            ((reaction.ts.energy + reaction.ts.zpe) - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL,
            ((reaction.ts.energy + reaction.ts.zpe) - sum([(opt.species.energy + opt.species.zpe) for opt in reaction.prod_opt])) * constants.AUtoKCAL,
        ]
        if any([bi < 0 for bi in barriers]):
            tun = ''
        else:
            imgfreq=-reaction.ts.reduced_freqs[0]
            negFreq_factor=negFreq_factor
            logFile.write("\tTS negFreq factor: {}\n".format(negFreq_factor))
            logFile.write("\t\tOriginal imaginary frequency: {}\n".format(imgfreq))
            if n_uq == 0:
                imgfreq=-reaction.ts.reduced_freqs[0]
                logFile.write("\t\tUpdated imaginary frequency: {}.\n\t\t\tFrequency should be unchanged.\n".format(imgfreq))
                barrierArrayNegFreq.append(imgfreq)
            elif n_uq > 0:
                imgfreq=-reaction.ts.reduced_freqs[0]
                imgfreq=imgfreq*negFreq_factor
                logFile.write("\t\tUpdated imaginary frequency: {}\n".format(imgfreq))
                barrierArrayNegFreq.append(imgfreq)
            tun = tun_tpl.format(cutoff=min(barriers),
                                 imfreq=imgfreq,
                                 welldepth1=barriers[0],
                                 welldepth2=barriers[1])



        if len(reaction.products) == 1:
            prod_name = self.well_names[reaction.products[0].chemid]
        else:
            long_name = '_'.join(sorted([str(pi.chemid) for pi in reaction.products]))
            prod_name = self.bimolec_names[long_name]

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
            energy = ((reaction.ts.energy + reaction.ts.zpe) - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL
            barrier_factor=barrier_factor 
            if n_uq == 0:
                barrier_factor=0
                logFile.write("\tBarrier factor: {}\n".format(barrier_factor))
                logFile.write("\t\tOriginal barrier energy: {}\n".format(energy))
            if n_uq > 0:
                logFile.write("\tBarrier factor: {}\n".format(barrier_factor))
                logFile.write("\t\tOriginal barrier energy: {}\n".format(energy))
                energy=energy+barrier_factor
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
            f = open(reaction.instance_name + '_uq_' + str(n_uq) + '.mess', 'w')
        f.write(mess_ts)
        f.close()

        barrierArrayNegFreq.pop(0)
        barrierArrayPosFreq.pop(0)
        barrierArrayEnergy.pop(0)
        e=barrierArrayEnergy
        frPos=barrierArrayPosFreq
        frNeg=barrierArrayNegFreq

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

    def run(self, n):
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
        queue_name=self.par.par['queuing']
        submitscript = 'run_mess' + constants.qext[self.par.par['queuing']]
        with open(submitscript, 'a') as qu: 
            if self.par.par['queue_template'] == '':
                if self.par.par['queuing'] == 'pbs':
                    qu.write((tpl_head).format(name='mess', ppn=self.par.par['ppn'], queue_name=self.par.par['queue_name'], dir='me'))
                elif self.par.par['queuing'] == 'slurm':
                    qu.write((tpl_head).format(name='mess', ppn=self.par.par['ppn'], queue_name=self.par.par['queue_name'], dir='me'), slurm_feature='')
            else:
                qu.write(tpl_head)

        '''
           copy submit script so that new runs can be done from run_mess_n files
           need to check for jobs running
           if < x jobs running submit next batch
           runs through counter adding up to n mess jobs total
           need to optimize the number that can be submitted
        '''

        #if uq=1 then run section below
        #else just run mess command 

        subcp= 'run_mess_uq' + constants.qext[self.par.par['queuing']]

        i=0 #counter for jobs
        x=self.par.par['uq_max_runs'] #max jobs running at once, can make this an input parameter at somepoint if neccessary
        n=n #total number of mess jobs to run
        xcounter=x #total jobs that have been submitted throughout process
        previousLoop=0 #number of running jobs on previous run so that jobs are not double counted as finishing
        pids=[] #list of job pids
        while(i<n):
                while(i<xcounter):
                    copyfile(submitscript, subcp)
                    with open(subcp, 'a') as f:
                        if self.par.par['queue_template'] == '':
                            if self.par.par['queuing'] == 'pbs':
                                f.write((tpl).format(n=i))
                            elif self.par.par['queuing'] == 'slurm':
                                f.write((tpl).format(n=i))
                        else:
                            f.write((tpl).format(n=i))
                         
                    command = [constants.qsubmit[self.par.par['queuing']], subcp ]
                    #command = [constants.qsubmit[self.par.par['queuing']], submitscript ]
                    
                    process = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
                    out, err = process.communicate()
                    out = out.decode()
                    if self.par.par['queuing'] == 'pbs':
                        pid = out.split('\n')[0].split('.')[0]
                    elif self.par.par['queuing'] == 'slurm':
                        pid = out.split('\n')[0].split()[-1]
                    pids.append(pid)
                    time.sleep(5)
                    i=i+1 

                a=0 #job being evaluated
                runningJobs=0 #number of running jobs
                exit=0 #bool to exit loop
                while(exit == 0):
                    devnull = open(os.devnull, 'w')
                    currentLoop=0 #number of jobs finished in current loop through pids[]
                    while(a<len(pids)):
                        pida=pids[a]
                        if self.par.par['queuing'] == 'pbs':
                            command = 'qstat -f | grep ' + '"Job Id: ' + pida + '"' + ' > /dev/null'
                        elif self.par.par['queuing'] == 'slurm':
                            command = 'scontrol show job ' + pida + ' | grep "JobId=' + pida + '"' + ' > /dev/null'

                        stat = int(subprocess.call(command, shell=True, stdout=devnull, stderr=devnull))

                        if stat == 0:
                            #time.sleep(1)
                            a=a+1
                            runningJobs=runningJobs+1

                        if stat == 1:
                            currentLoop=currentLoop+1
                            a=a+1
                        
                        if a==len(pids) and runningJobs==x:
                            a=0
                            runningJobs=0
                            currentLoop=0

                        if((currentLoop > previousLoop) and (a==len(pids))):
                            xcounter=xcounter+currentLoop-previousLoop
                            previousLoop=currentLoop
                            currentLoop=0
                            runningJobs=0
                            a=0
                            exit=1


                if xcounter > n:
                    xcounter=n
                 
        return 0
    """
    def norm(self, values, variable, varType, names,n):
    #function to normalize UQ values
        fi='{}_uq.log'.format(variable)
        fio=open(fi,'w')
        fio.write("{} UQ Normalization\nEnergy\tNormalized Value\n".format(variable))
        if varType == 0: #energy/barrier
            rxns=names
            data=values
            allNormVals=[]
            allMin=[]
            allMax=[]
            allAvg=[]
            allStDev=[]
            ogs=data[0]
            for i, x in enumerate(data):
                npData=np.array(data[i])
                low=np.min(npData)
                allMin.append(low)
                high=np.max(npData)
                allMax.append(high)
                avg=np.mean(npData)
                allAvg.append(avg)
                stDev=np.std(npData)
                allStDev.append(stDev)
                og=data[i][0]
                normVal=[]
                for j in data[i]:
                    val=2*((j-low)/(high-low))-1
                    normVal.append(val)
                allNormVals.append(normVal)
                npAllNormVals=np.array(allNormVals)
            if len(npAllNormVals) == n:
                for k, x in enumerate(npAllNormVals):
                    for j, y in enumerate(npAllNormVals[k]):
                        en=npAllNormVals[k][j]
                        en=float(en)
                        enf=('{:.4}'.format(en))
                        #str_j = str(enf)[1 : -1]
                        datap=data[k][j]
                        str_k = str(datap)[1 : -1]
                        en1=float(str_k)
                        en1f=('{:.4}'.format(en1))
                        fio.write("\n{}\t{}".format(en1f,str(enf)))
                    fio.write("\n")
    
        elif varType == 1:
            name=variable
            rxn=names
            data=values
            allnormVal=[]
            allNormVals=[]
            allMin=[]
            allMax=[]
            allAvg=[]
            allStDev=[]
            ogs=data[0]
            for i, x in enumerate(data):
                npData=np.array(data[i])
                low=np.min(npData)
                allMin.append(low)
                high=np.max(npData)
                allMax.append(high)
                avg=np.mean(npData)
                allAvg.append(avg)
                stDev=np.std(npData)
                allStDev.append(stDev)
                og=data[i][0]
                normVal=[]
                for j in data[i]:
                    lnx=np.log(j)
                    lnmin=np.log(low)
                    lnmax=np.log(high)       
                    val=2*((lnx-lnmin)/(lnmax-lnmin))-1
                    normVal.append(val)
                allNormVals.append(normVal)
                npAllNormVals=np.array(allNormVals)
            if len(npAllNormVals) == n:
                for k, x in enumerate(npAllNormVals):
                    for j, y in enumerate(npAllNormVals[k]):
                        en=npAllNormVals[k][j]
                        en=float(en)
                        enf=('{:.4}'.format(en))
                        #str_j = str(enf)[1 : -1]
                        datap=data[k][j]
                        str_k = str(datap)[1 : -1]
                        en1=float(str_k)
                        en1f=('{:.8}'.format(en1))
                        fio.write("\n{}\t{}".format(en1f,str(enf)))
                    fio.write("\n")
        elif varType == 2:

        return 0
        """
