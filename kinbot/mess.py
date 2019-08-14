###################################################
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
import os
import subprocess
import time
import pkg_resources
import logging
from kinbot import constants
from kinbot import frequencies


class MESS:
    """
    Class that reads and write all MESS files and runs MESS
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

    def write_input(self):
        """
        write the input for all the wells, bimolecular products and barriers
        both in a separate file, as well as in one large ME file
        """
        # create short names for all the species, bimolecular products and barriers
        self.create_short_names()
        header = self.write_header()

        # filter ts's with the same reactants and products:
        ts_unique = {}  # key: ts name, value: [prod_name, energy]
        for index, reaction in enumerate(self.species.reac_obj):
            if self.species.reac_ts_done[index] == -1:
                prod_name = '_'.join([str(pi.chemid) for pi in reaction.products])
                energy = reaction.ts.energy
                zpe = reaction.ts.zpe
                new = 1
                remove = []
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
        well_blocks = {}
        ts_blocks = {}
        bimolec_blocks = {}
        well_blocks[self.species.chemid] = self.write_well(self.species)
        for index, reaction in enumerate(self.species.reac_obj):
            ts=open("ts.log", 'a')
            ts.write('{0} {1}'.format(reaction, "\n"))
            ts.close()
            if reaction.instance_name in ts_unique:
                ts_blocks[reaction.instance_name] = self.write_barrier(reaction)
                if len(reaction.products) == 1:
                    st_pt = reaction.prod_opt[0].species
                    well_blocks[st_pt.chemid] = self.write_well(st_pt)
                else:
                    bimol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                    bimolec_blocks[bimol_name] = self.write_bimol([opt.species for opt in reaction.prod_opt])

        # write the mess input file
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

        f_out = open('me/mess.inp', 'w')
        f_out.write(header + '\n!****************************************\n')
        f_out.write(wells)
        f_out.write(bimols)
        f_out.write(tss)
        f_out.write(dum)
        f_out.write('\n!****************************************\nEnd ! end kinetics\n')
        f_out.close()

        return 0

    def write_bimol(self, species_list):
        """
        Create the block for MESS for a bimolecular product.
        well0: reactant on this PES (zero-energy reference)
        """
        # open the templates
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
                #reduced freqs used always bc HIR creates more accurate mess input
                for i, fr in enumerate(species.reduced_freqs):
                    if i == 0:
                        freq += '{:.4f}'.format(fr)
                    elif i > 0 and i % 3 == 0:
                        freq += '\n            {:.4f}'.format(fr)
                    else:
                        freq += '    {:.4f}'.format(fr)
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
            energy = (sum([sp.energy for sp in species_list]) + sum([sp.zpe for sp in species_list]) -
                      (self.species.energy + self.species.zpe)) * constants.AUtoKCAL

        
        bimol = tpl.format(chemids=name,
                           fragments=fragments,
                           ground_energy=energy)

        f = open(pr_name + '.mess', 'w')
        f.write(bimol)
        f.close()

        return bimol

    def write_well(self, species):
        """
        Create the block for MESS for a well.
        well0: reactant on this PES (zero-energy reference)
        """
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
        #reduced freqs used for mess input to make more accurate with HIR corrections
        for i, fr in enumerate(species.reduced_freqs):
            if i == 0:
                freq += '{:.4f}'.format(fr)
            elif i > 0 and i % 3 == 0:
                freq += '\n            {:.4f}'.format(fr)
            else:
                freq += '    {:.4f}'.format(fr)

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

        f = open(str(species.chemid) + '.mess', 'w')
        f.write(mess_well)
        f.close()

        return mess_well

    def write_barrier(self, reaction):
        """
        Create the block for a MESS barrier.
        """
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
        #reduced freqs used for better accuracy of mess input files
        for i, fr in enumerate(reaction.ts.reduced_freqs[1:]):
            if i == 0:
                freq += '{:.4f}'.format(fr)
            elif i > 0 and i % 3 == 0:
                freq += '\n            {:.4f}'.format(fr)
            else:
                freq += '    {:.4f}'.format(fr)

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
            tun = tun_tpl.format(cutoff=min(barriers),
                                 imfreq=-reaction.ts.reduced_freqs[0],
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

        f = open(reaction.instance_name + '.mess', 'w')
        f.write(mess_ts)
        f.close()

        return mess_ts

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

    def run(self):
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
            tpl_head = f.read() #file for main body of submission script
        
        q_file = pkg_resources.resource_filename('tpl', self.par.par['queuing'] + '_mess.tpl')
        with open(q_file) as f:
            tpl = f.read() #file for submitting mess.inp 'mess mess.inp'
        queue_name=self.par.par['queuing']
        submitscript = 'run_mess' + constants.qext[self.par.par['queuing']]
        with open(submitscript, 'a') as qu:
            if self.par.par['queue_template'] == '':
                if self.par.par['queuing'] == 'pbs':
                    #pbs.tpl has no analogue to slurm_feature in slurm.tpl
                    qu.write((tpl_head + tpl).format(name='mess', ppn=self.par.par['ppn'], queue_name=self.par.par['queue_name'], dir='me'))
                elif self.par.par['queuing'] == 'slurm':
                    qu.write((tpl_head + tpl).format(name='mess', ppn=self.par.par['ppn'], queue_name=self.par.par['queue_name'], dir='me', slurm_feature=''))
            else:
                #run_mess file - generate with tpl files
                q_file = pkg_resources.resource_filename('tpl', self.par.par['queuing'] + '.tpl')
                with open(q_file) as f:
                    tpl_head=f.read()
                qu.write(tpl_head.format(name='mess', pn=self.par.par['ppn'], queue_name=par.par['queue_name'], dir='me', slurm_feature=''))
                qu.write(tpl)

        command = [constants.qsubmit[self.par.par['queuing']], submitscript ]
        process = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        out = out.decode()
        if self.par.par['queuing'] == 'pbs':
            pid = out.split('\n')[0].split('.')[0]
        elif self.par.par['queuing'] == 'slurm':
            pid = out.split('\n')[0].split()[-1]

        while 1:  
            devnull = open(os.devnull, 'w')
            if self.par.par['queuing'] == 'pbs':
                command = 'qstat -f | grep ' + '"Job Id: ' + pid + '"' + ' > /dev/null'
            elif self.par.par['queuing'] == 'slurm':
                command = 'scontrol show job ' + pid + ' | grep "JobId=' + pid + '"' + ' > /dev/null'
            if int(subprocess.call(command, shell=True, stdout=devnull, stderr=devnull)) == 0:
                time.sleep(1)
            else:
                break
        return 0



