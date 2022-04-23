import os
import numpy as np
import subprocess
import time
import pkg_resources
from collections import Counter

from kinbot import constants
from kinbot import frequencies
from kinbot.uncertaintyAnalysis import UQ


class MESS:
    """
    Class that reads and writes MESS files
    UQ analysis parameter (uq) can be used to generate 'n' number of mess input files
    with the following parameters randomized within the alloted UQ tolerance.
    By default UQ ranges are set to these values
       1. Stationary point energy (E+ZPE, +/- 0.5 kcal/mol)
       2. Barrier (E+ZPE, +/- 1.0 kcal/mol)
       3. Frequencies (cm-1 */ / 1.2)
    Default parameters were chosen/based on the following paper:  Goldsmith, C. F. PCI, 2013, 177-185
    New parameters can be set within the input 
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
        with open(pkg_resources.resource_filename('tpl', 'mess_fragment_OH.tpl')) as f:
            self.fragmenttplOH = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_pstfragment.tpl')) as f:
            self.pstfragmenttpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_hinderedrotor.tpl')) as f:
            self.hinderedrotortpl = f.read()
        with open(pkg_resources.resource_filename('tpl', 'mess_freerotor.tpl')) as f:
            self.freerotortpl = f.read()
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
                                       e_coll=round(constants.epsilon[self.par['collider']], 2),
                                       s_coll=constants.sigma[self.par['collider']],
                                       m_coll=constants.mass[self.par['collider']],
                                       e_well=round(self.par['epsilon']),
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
                self.ts_names[reaction.instance_name] = 'ts_{}'.format(len(self.ts_names) + 1)
                if len(reaction.products) == 1:
                    st_pt = reaction.products[0]
                    if st_pt.chemid not in self.well_names:
                        self.well_names[st_pt.chemid] = 'w_{}'.format(len(self.well_names) + 1)
                elif len(reaction.products) == 2:
                    for st_pt in reaction.products:
                        if st_pt.chemid not in self.fragment_names:
                            self.fragment_names[st_pt.chemid] = 'fr_{}'.format(len(self.fragment_names) + 1)
                    bimol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                    if bimol_name not in self.bimolec_names:
                        self.bimolec_names[bimol_name] = 'b_{}'.format(len(self.bimolec_names) + 1)
                else:
                    # TER MOLECULAR
                    for st_pt in reaction.products:
                        if st_pt.chemid not in self.fragment_names:
                            self.fragment_names[st_pt.chemid] = 'fr_{}'.format(len(self.fragment_names) + 1)
                    termol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                    if termol_name not in self.termolec_names:
                        self.termolec_names[termol_name] = 't_{}'.format(len(self.termolec_names) + 1)


    def write_input(self, qc):
        """
        write the input for all the wells, bimolecular products and barriers
        both in a separate file, as well as in one large ME file
        """
        uq = UQ(self.par)

        # create short names for all the species, bimolecular products and barriers
        self.create_short_names()
        header = self.write_header()

        # filter ts's with the same reactants and products:
        ts_unique = {}  # key: ts name, value: [prod_name, energy]
        ts_all = {}
        for index, reaction in enumerate(self.species.reac_obj):
            if self.species.reac_ts_done[index] == -1:
                rxnProds = []
                for x in reaction.products:
                    rxnProds.append(x.chemid)
                rxnProds.sort()
                prod_name = '_'.join([str(pi) for pi in rxnProds])
                new = 1
                remove = []
                ts_all[reaction.instance_name] = [prod_name, reaction.ts.energy  + reaction.ts.zpe]
                for ts in ts_unique:
                    if ts_unique[ts][0] == prod_name:
                        # check for the barrier with the lowest energy
                        if ts_unique[ts][1] > reaction.ts.energy + reaction.ts.zpe:
                            # remove the current barrier
                            remove.append(ts)
                        else:
                            new = 0
                for ts in remove:
                    ts_unique.pop(ts, None)
                if new:
                    ts_unique[reaction.instance_name] = [prod_name, reaction.ts.energy + reaction.ts.zpe]

        # write the mess input for the different blocks
        for uq_iter in range(self.par['uq_n']):
            well_blocks = {}
            ts_blocks = {}
            bimolec_blocks = {}
            termolec_blocks = {}
            termolec_ts_blocks = {}
            barrierless_blocks = {}
            allTS = {}
            # arrays reset at the start of each iteration to hold new values
            written_bimol_names = []
            written_termolec_names = []

            well_energy_add = uq.calc_factor('energy', self.species.chemid, uq_iter)
            well_freq_factor = uq.calc_factor('freq', self.species.chemid, uq_iter)
            well_blocks[self.species.chemid] = self.write_well(self.species,
                                                               well_energy_add,
                                                               well_freq_factor,
                                                               uq_iter)
            
            for index, reaction in enumerate(self.species.reac_obj):
                if reaction.instance_name in ts_all:
                    barrier_add = uq.calc_factor('barrier', reaction.instance_name, uq_iter)
                    freq_factor = uq.calc_factor('freq', reaction.instance_name, uq_iter)
                    imagfreq_factor = uq.calc_factor('imagfreq', reaction.instance_name, uq_iter)
        
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

                    allTS[reaction.instance_name], zeroenergy = self.write_barrier(reaction,
                                                                                   index,
                                                                                   left_zeroenergy,
                                                                                   right_zeroenergy,
                                                                                   barrier_add,
                                                                                   freq_factor,
                                                                                   imagfreq_factor,
                                                                                   uq_iter)

                # Only write products once, stops duplicate product writing
                if reaction.instance_name in ts_unique:
                    ts_blocks[reaction.instance_name] = allTS[reaction.instance_name]
                    if len(reaction.products) == 1:
                        st_pt = reaction.prod_opt[0].species
                        energy_add = uq.calc_factor('energy', st_pt.chemid, uq_iter)
                        freq_factor = uq.calc_factor('freq', st_pt.chemid, uq_iter)
                        well_blocks[st_pt.chemid] = self.write_well(st_pt,
                                                                    energy_add,
                                                                    freq_factor,
                                                                    uq_iter)
                    elif len(reaction.products) == 2:
                        bimol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                        energy_add = uq.calc_factor('energy', bimol_name, uq_iter)
                        freq_factor = uq.calc_factor('freq', bimol_name, uq_iter)
                        bimolec_blocks[bimol_name] = self.write_bimol([opt.species for opt in reaction.prod_opt],
                                                                      energy_add,
                                                                      freq_factor,
                                                                      uq_iter)
                        written_bimol_names.append(bimol_name)
                    else:
                        # termol
                        termolec_ts_blocks[reaction.instance_name] = allTS[reaction.instance_name]
                        termol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                        termolec_blocks[termol_name] = self.write_termol([opt.species for opt in reaction.prod_opt], 
                                                                         reaction,
                                                                         uq_iter)
                        written_termolec_names.append(termol_name)

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

        uq.format_uqtk_data() 

        return 0


    def write_barrierless(self, species_list, reaction, energy_add, freq_factor, uq_iter):

        if len(reaction.products) == 2:
            barrierless = self.write_bimol(species_list,
                                           energy_add,
                                           freq_factor,
                                           uq_iter,
                                           bless=1)
        else:
            barrierless = self.write_termol(species_list, reaction, uq_iter, bless=1)

        return barrierless

    def write_termol(self, species_list, reaction, uq_iter, bless=0):
        # Create the dummy MESS block for ter-molecular products.
        termol = ''
        terPr_name = '_'.join(sorted([str(species.chemid) for species in species_list]))
        prod_name = self.termolec_names[terPr_name]
        termol += self.termoltpl.format(name=prod_name, product=terPr_name)
        mess_iter = "{0:04d}".format(uq_iter)
        with open(terPr_name + '_' + mess_iter + '.mess', 'w') as f:
            f.write(termol)

        return termol


    def write_bimol(self, prod_list, well_add, freq_factor, uq_iter, bless=0):
        """
        Create the block for MESS for a bimolecular product.
        well0: reactant on this PES (zeroenergy reference)
        uq_n = number of uncertainty runs
        """

        fragments = ''
        if bless == 1:
            tot_nfreq = 0
            combined_freq = ''
            combined_hir = ''
        smi = []
        for nsp, species in enumerate(prod_list):
            smi.append(species.smiles)
            if species.natom > 1:

                if self.par['pes']:
                    name = '{{fr_name_{}}}'.format(species.chemid)
                else:
                    name = self.fragment_names[species.chemid] + ' ! ' + str(species.chemid)
                # molecule template
                if species.chemid == 170170000000000000002:  # exception for OH
                    fragments += self.fragmenttplOH.format(chemid=name,
                                                         smi=species.smiles,
                                                         natom=species.natom,
                                                         geom=self.make_geom(species),
                                                         symm=float(species.sigma_ext) / float(species.nopt),
                                                         freq=self.make_freq(species, freq_factor, 0))
                else:
                    fragments += self.fragmenttpl.format(chemid=name,
                                                         smi=species.smiles,
                                                         natom=species.natom,
                                                         geom=self.make_geom(species),
                                                         symm=float(species.sigma_ext) / float(species.nopt),
                                                         nfreq=len(species.reduced_freqs),
                                                         freq=self.make_freq(species, freq_factor, 0),
                                                         hinderedrotor=self.make_rotors(species, freq_factor),
                                                         nelec=1,
                                                         mult=species.mult)
                if bless == 1:
                    tot_nfreq += len(species.reduced_freqs)
                    combined_freq += self.make_freq(species, freq_factor, 0)
                    combined_hir += self.make_rotors(species, freq_factor)

                    if nsp == 0: 
                        combined_mult = species.mult
                        frag1 = self.pstfragmenttpl.format(chemid=name,
                                                           smi=species.smiles,
                                                           natom=species.natom,
                                                           geom=self.make_geom(species))
                    if nsp == 1: 
                        if combined_mult == 1 and species.mult == 1:
                            combined_mult = 1
                        elif combined_mult == 2 and species.mult == 1:
                            combined_mult = 2
                        elif combined_mult == 1 and species.mult == 2:
                            combined_mult = 2
                        elif combined_mult == 2 and species.mult == 2:
                            combined_mult = 1
                        frag2 = self.pstfragmenttpl.format(chemid=name,
                                                           smi=species.smiles,
                                                           natom=species.natom,
                                                           geom=self.make_geom(species))
            else:
                if self.par['pes']:
                    name = '{{fr_name_{}}}'.format(species.chemid)
                else:
                    name = self.fragment_names[species.chemid] + ' ! ' + str(species.chemid)

                fragments += self.atomtpl.format(chemid=name,
                                                 element=species.atom[0],
                                                 nelec=1,
                                                 mult=species.mult)
                if bless == 1:
                    if nsp == 0: 
                        frag1 = self.pstfragmenttpl.format(chemid=name,
                                                           smi=species.smiles,
                                                           natom=species.natom,
                                                           geom=self.make_geom(species))
                    if nsp == 1: 
                        frag2 = self.pstfragmenttpl.format(chemid=name,
                                                           smi=species.smiles,
                                                           natom=species.natom,
                                                           geom=self.make_geom(species))
 

        pr_name = '_'.join(sorted([str(species.chemid) for species in prod_list]))
        if self.par['pes']:
            name = '{{name}} ! {} {}'.format(smi[0], smi[1])
            energy = '{ground_energy}'
        else:
            if bless == 0:
                name = '{} ! {}'.format(self.bimolec_names[pr_name], pr_name)
            elif bless == 1:
                name = '{} ! barrierless'.format(self.barrierless_names[pr_name])

            energy = (sum([sp.energy for sp in prod_list]) + sum([sp.zpe for sp in prod_list]) 
                      - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL
            energy += well_add
            energy = round(energy, 2)
        
        if bless == 0:
            bimol = self.bimoltpl.format(chemids=name,
                                         smi=species.smiles,
                                         fragments=fragments,
                                         ground_energy=energy)

        elif bless == 1:
            values = list(self.barrierless_names.values())
            index = values.index(pr_name)
            stoich = ''
            el_counter = Counter(self.species.atom)
            for el in constants.elements:
                if el_counter[el]:
                    stoich += '{}{}'.format(el, el_counter[el])
            bimol = self.blbimoltpl.format(barrier='nobar_{}'.format(index),
                                           reactant=self.well_names[self.species.chemid],
                                           prod=name,
                                           stoich=stoich,
                                           frag1=frag1,
                                           frag2=frag2,
                                           nfreq=tot_nfreq,
                                           freq=combined_freq,
                                           mult=combined_mult,
                                           hinderedrotor=combined_hir, 
                                           fragments=fragments,
                                           ground_energy=energy)

        with open('{}_{:04d}.mess'.format(pr_name, uq_iter), 'w') as f:
            f.write(bimol)

        return bimol


    def write_well(self, species, well_add, freq_factor, uq_iter):
        """
        Create the block for MESS for a well.
        well0: reactant on this PES (zeroenergy reference)
        """

        if self.par['pes']:
            name = '{name}'
            zeroenergy = '{zeroenergy}'
        else:
            name = self.well_names[species.chemid] + ' ! ' + str(species.chemid)
            zeroenergy = ((species.energy + species.zpe) - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL
            zeroenergy += well_add
            zeroenergy = round(zeroenergy, 2)

        nunq_confs = 0  # number of unique conformers
        for co in self.species.conformer_index:
            if co >= 0:
                nunq_confs += 1

        if not self.par['multi_conf_tst'] or nunq_confs == 1: 
            mess_well = self.welltpl.format(chemid=name,
                                            smi=species.smiles,
                                            natom=species.natom,
                                            geom=self.make_geom(species),
                                            symm=float(species.sigma_ext) / float(species.nopt),
                                            nfreq=len(species.reduced_freqs),
                                            freq=self.make_freq(species, freq_factor, 0),
                                            hinderedrotor=self.make_rotors(species, freq_factor),
                                            nelec=1,
                                            mult=species.mult,
                                            zeroenergy=zeroenergy)
        else:
            rrho = ''
            for ci, co in enumerate(self.species.conformer_index):
                corerr = self.corerrtpl.format(symm=float(species.sigma_ext) / float(species.nopt))
                rrho += self.rrhotpl.format(natom=species.natom,
                                            geom=self.make_geom(species),
                                            core=corerr,
                                            nfreq=len(species.kinbot_freqs),
                                            freq=self.make_freq(species, freq_factor, 0),
                                            rotors=self.make_rotors(species, freq_factor),
                                            tunneling='',
                                            nelec=1,
                                            mult=species.mult,
                                            zeroenergy=zeroenergy)
 

        with open('{}_{:04d}.mess'.format(species.chemid, uq_iter), 'w') as f:
            f.write(mess_well)

        return mess_well

    def write_barrier(self, reaction, index, left_zeroenergy, right_zeroenergy, barrier_add, freq_factor, imagfreq_factor, uq_iter):
        """
        Create the block for a MESS barrier.
        """

        left_zeroenergy += barrier_add
        right_zeroenergy += barrier_add

        # write tunneling block
        if left_zeroenergy < 0 or right_zeroenergy < 0:
            tun = f'! barrier is submerged {left_zeroenergy} {right_zeroenergy}'
        elif self.par['pes'] == 0:
            tun = self.tunneltpl.format(cutoff=round(min(left_zeroenergy, right_zeroenergy), 2),
                                        imfreq=round(-reaction.ts.reduced_freqs[0] * imagfreq_factor, 2),
                                        welldepth1=round(left_zeroenergy, 2),
                                        welldepth2=round(right_zeroenergy, 2))
        else: 
            tun = self.tunneltpl.format(cutoff='{cutoff}',
                                        imfreq=round(-reaction.ts.reduced_freqs[0] * imagfreq_factor, 2),
                                        welldepth1='{welldepth1}',
                                        welldepth2='{welldepth2}')

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
            zeroenergy = round(left_zeroenergy, 2)
    
        # TODO working here
        if self.species.reac_type[index] == 'barrierless_saddle':
            freq = self.make_freq(reaction.prod_opt[0].species, freq_factor, 0) + \
                   self.make_freq(reaction.prod_opt[1].species, freq_factor, 0) 
            rotors = self.make_rotors(reaction.prod_opt[0].species, freq_factor) + \
                     self.make_rotors(reaction.prod_opt[1].species, freq_factor) 
            nfreq = len(reaction.prod_opt[0].species.reduced_freqs) + \
                    len(reaction.prod_opt[1].species.reduced_freqs)
            if self.par['pes']:
                prodzeroenergy = '{prodzeroenergy}'
            else:
                prodzeroenergy = ((reaction.prod_opt[0].species.energy + reaction.prod_opt[0].species.zpe + \
                                 reaction.prod_opt[1].species.energy + reaction.prod_opt[1].species.zpe) - \
                                 (self.species.energy + self.species.zpe)) * constants.AUtoKCAL

            outerts = self.psttpl.format(natom1=reaction.prod_opt[0].species.natom,
                                         geom1=self.make_geom(reaction.prod_opt[0].species),
                                         natom2=reaction.prod_opt[1].species.natom,
                                         geom2=self.make_geom(reaction.prod_opt[1].species),
                                         symm=float(reaction.ts.sigma_ext) / float(reaction.ts.nopt),
                                         prefact='prefactor',
                                         exponent=3.,
                                         nfreq=nfreq,
                                         freq=freq,
                                         hinderedrotor=rotors,
                                         nelec=1,
                                         mult=reaction.ts.mult,
                                         prodzeroenergy=prodzeroenergy
                                         )
            twotst = self.twotstpl.format(outerts=outerts)
            corerr = self.corerrtpl.format(symm=float(reaction.ts.sigma_ext) / float(reaction.ts.nopt))
            rrho = self.rrhotpl.format(natom=reaction.ts.natom,
                                       geom=self.make_geom(reaction.ts),
                                       core=corerr,
                                       nfreq=len(reaction.ts.reduced_freqs)-1,
                                       freq=self.make_freq(reaction.ts, freq_factor, 1),
                                       rotors=self.make_rotors(reaction.ts, freq_factor, norot=self.ts_names[reaction.instance_name]),
                                       tunneling='',
                                       nelec=1,
                                       mult=reaction.ts.mult,
                                       zeroenergy=zeroenergy)
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
                                       freq=self.make_freq(reaction.ts, freq_factor, 1),
                                       rotors=self.make_rotors(reaction.ts, freq_factor, norot=self.ts_names[reaction.instance_name]),
                                       tunneling=tun,
                                       nelec=1,
                                       mult=reaction.ts.mult,
                                       zeroenergy=zeroenergy)
            mess_barrier = self.barriertpl.format(rxn_name=name,
                                                  chemid_reac=chemid_reac,
                                                  chemid_prod=chemid_prod,
                                                  long_rxn_name=long_rxn_name,
                                                  model=rrho)

        with open('{}_{:04d}.mess'.format(reaction.instance_name, uq_iter), 'w') as f:
            f.write(mess_barrier)
    
        return mess_barrier, zeroenergy 

    def run(self):
        """
        Submit the pbs/slurm file to the queue
        wait for the mess run to finish
        """

        submitscript = 'run_mess' + constants.qext[self.par['queuing']]

        pids = []  # list of job pids

        uq_iter = 0
        pid_stats = []
        while uq_iter < self.par['uq_n']:
            self.write_submitscript(submitscript, uq_iter)
            while len(pids) > self.par['uq_max_runs']:
                time.sleep(5)
                for pid in pids:
                    stat = self.check_running(pid)
                    if stat == 0:
                        pids.remove(pid)
                        pid_stats.append(stat)
            
            pid = self.submit(submitscript)
            pids.append(pid)

            if self.par['uq_n'] < self.par['uq_max_runs']:
                stat = 1
                while stat != 0:
                    stat = self.check_running(pid)
                    time.sleep(5)
                pid_stats.append(stat)
            uq_iter += 1
  
        if all(stat == 0 for s in pid_stats):
            return 0

    def write_submitscript(self, submitscript, uq_iter):
        """
        write a pbs or slurm file for the me/all.inp mess input file
        """

        if self.par['queue_template'] == '':
            q_file = pkg_resources.resource_filename('tpl', self.par['queuing'] + '.tpl')
        else:
            q_file = self.par['queue_template']
        with open(q_file) as f:
            tpl_head = f.read()

        q_file = pkg_resources.resource_filename('tpl', self.par['queuing'] + '_mess_uq.tpl')
        with open(q_file) as f:
            tpl = f.read()

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
        return 0

    def submit(self, submitscript):
        command = [constants.qsubmit[self.par['queuing']], submitscript]
        process = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        out = out.decode()
        if self.par['queuing'] == 'pbs':
            pid = out.split('\n')[0].split('.')[0]
        elif self.par['queuing'] == 'slurm':
            pid = out.split('\n')[0].split()[-1]
        return pid

    def check_running(self, pid):
        devnull = open(os.devnull, 'w')
        if self.par['queuing'] == 'pbs':
            command = 'qstat -f | grep ' + '"Job Id: ' + pid + '"' + ' > /dev/null'
        elif self.par['queuing'] == 'slurm':
            command = 'scontrol show job ' + pid + ' | grep "JobId=' + pid + '"' + ' > /dev/null'
 
        stat = int(subprocess.call(command, shell=True, stdout=devnull, stderr=devnull))
        return stat

    def make_geom(self, species):
        geom = ''
        for i, at in enumerate(species.atom):
            x, y, z = species.geom[i]
            geom += '        {} {:.6f} {:.6f} {:.6f}\n'.format(at, x, y, z)
        return geom[:-1]

    def make_freq(self, species, factor, wellorts):
        """
        Frequencies are scaled with factor in UQ.
        At 100 cm-1 the scaling is applied as is.
        For lower frequencies the scaling is amplified.
        For higher frequencies the scaling is dampened.
        """
        freq = '        '
        #wellorts: 0 for wells and 1 for saddle points
        if wellorts == 0:
            frequencies = species.reduced_freqs
        else:
            frequencies = species.reduced_freqs[1:]
        for i, fr in enumerate(frequencies):
            if factor >= 1:
                fr = fr * (1 / fr * (factor - 1 ) * 100 + 1)
            else:
                fr = fr / ( 1 / fr * (1 - factor) * 100 + 1)
            freq += '{:.1f} '.format(fr)
            if i % 3 == 2:
                freq += '\n        '
        return(freq[:-1])

    def make_rotorpot(self, species, i, rot, freq_factor):
        rotortype = 'hindered'
        rotorsymm = self.rotorsymm(species, rot)
        ens = species.hir.hir_energies[i]
        rotorpot_num = [(ei - ens[0]) * constants.AUtoKCAL for ei in ens]
        maxen = max(rotorpot_num)
        # solution for 6-fold symmetry, not general enough
        if species.hir.nrotation // rotorsymm == 2:  # MESS needs at least 3 potential points
            fit_angle = 15. * 2. * np.pi / 360. 
            fit_energy = species.hir.get_fit_value(fit_angle)  # kcal/mol
            rotorpot_num.insert(1, fit_energy)
            rotorpot_num = [freq_factor * rpn for rpn in rotorpot_num]
            rotorpot = ' '.join(['{:.2f}'.format(ei) for ei in rotorpot_num[:species.hir.nrotation // rotorsymm + 1]])
        else:
            rotorpot_num = [freq_factor * rpn for rpn in rotorpot_num]
            rotorpot = ' '.join(['{:.2f}'.format(ei) for ei in rotorpot_num[:species.hir.nrotation // rotorsymm]])
        rotorpot = '        {}'.format(rotorpot)
        if maxen < self.par['free_rotor_thrs']:
            rotortype = 'free'
        return rotorpot, rotortype

    def rotorsymm(self, species, rot):
        return species.sigma_int[rot[1]][rot[2]]

    def nrotorpot(self, species, rot): 
        rotorsymm = self.rotorsymm(species, rot)
        if species.hir.nrotation // rotorsymm > 2:
            return species.hir.nrotation // rotorsymm
        else:
            return species.hir.nrotation // rotorsymm + 1

    def make_rotors(self, species, freq_factor, norot=None):
        rotors = []
        if self.par['rotor_scan']:
            for i, rot in enumerate(species.dihed):
                if norot is not None:
                    if frequencies.skip_rotor(norot, rot) == 1:
                        continue
                rotorpot, rotortype = self.make_rotorpot(species, i, rot, freq_factor)
                if rotortype == 'hindered':
                    rotors.append(self.hinderedrotortpl.format(group=' '.join([str(pi + 1) for pi in frequencies.partition(species, rot, species.natom)[0][1:]]),
                                                               axis='{} {}'.format(str(rot[1] + 1), str(rot[2] + 1)),
                                                               rotorsymm=self.rotorsymm(species, rot),
                                                               nrotorpot=self.nrotorpot(species, rot),
                                                               rotorpot=rotorpot))
                elif rotortype == 'free':
                    rotors.append(self.freerotortpl.format(geom=self.make_geom(species),
                                                           natom=species.natom,
                                                           group=' '.join([str(pi + 1) for pi in frequencies.partition(species, rot, species.natom)[0][1:]]),
                                                           axis='{} {}'.format(str(rot[1] + 1), str(rot[2] + 1)),
                                                           ))

        rotors = '\n'.join(rotors)
        return rotors

    def get_zeroenergy(self, jobname, qc):
        energy = qc.get_qc_energy(jobname)[1]
        zpe = qc.get_qc_zpe(jobname)[1]
        return (energy + zpe) * constants.AUtoKCAL

    def check_running(self, pid):
        devnull = open(os.devnull, 'w')
        if self.par['queuing'] == 'pbs':
            command = 'qstat -f | grep ' + '"Job Id: ' + pid + '"' + ' > /dev/null'
        elif self.par['queuing'] == 'slurm':
            command = 'scontrol show job ' + pid + ' | grep "JobId=' + pid + '"' + ' > /dev/null'

        stat = int(subprocess.call(command, shell=True, stdout=devnull, stderr=devnull))
        return stat
