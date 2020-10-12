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
        self.pes = self.par['pes']
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

        # Barrierless short names
        if self.species.homolytic_scissions is not None:
            for hs in self.species.homolytic_scissions.hss:
                if hs.status == -1:
                    if len(hs.products) == 1:
                        st_pt = hs.products[0]
                        if st_pt.chemid not in self.well_names:
                            self.well_names[st_pt.chemid] = 'w_{}'.format(len(self.well_names) + 1)
                    elif len(hs.products) == 2:
                        for st_pt in hs.products:
                            if st_pt.chemid not in self.fragment_names:
                                self.fragment_names[st_pt.chemid] = 'fr_{}'.format(len(self.fragment_names) + 1)
                    bimol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in hs.products]))
                    if bimol_name not in self.bimolec_names:
                        self.bimolec_names[bimol_name] = 'b_{}'.format(len(self.bimolec_names) + 1)
                    else:
                        # TER MOLECULAR
                        for st_pt in hs.products:
                            if st_pt.chemid not in self.fragment_names:
                                self.fragment_names[st_pt.chemid] = 'fr_{}'.format(len(self.fragment_names) + 1)
                        termol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in hs.products]))
                        if termol_name not in self.termolec_names:
                            self.termolec_names[termol_name] = 't_{}'.format(len(self.termolec_names) + 1)


    def write_input(self, qc):
        """
        write the input for all the wells, bimolecular products and barriers
        both in a separate file, as well as in one large ME file
        """
        uq_obj = UQ(self.par)

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

            well_energyAdd = uq_obj.calc_factor('energy', self.species.chemid, uq_iter, self.pes)
            well_freqFactor = uq_obj.calc_factor('freq', self.species.chemid, uq_iter, self.pes)
            well_blocks[self.species.chemid] = self.write_well(self.species,
                                                                                well_energyAdd,
                                                                                well_freqFactor,
                                                                                uq_iter)
            
            for index, reaction in enumerate(self.species.reac_obj):
                if reaction.instance_name in ts_all:
                    barrierAdd = uq_obj.calc_factor('barrier', reaction.instance_name, uq_iter, self.pes)
                    freqFactor = uq_obj.calc_factor('freq', reaction.instance_name, uq_iter, self.pes)
                    imagfreqFactor = uq_obj.calc_factor('imagfreq', reaction.instance_name, uq_iter, self.pes)
        

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
                                                                       barrierAdd,
                                                                       freqFactor,
                                                                       imagfreqFactor,
                                                                       uq_iter)

                # Only write products once, stops duplicate product writing
                if reaction.instance_name in ts_unique:
                    ts_blocks[reaction.instance_name] = allTS[reaction.instance_name]
                    if len(reaction.products) == 1:
                        st_pt = reaction.prod_opt[0].species
                        energyAdd = uq_obj.calc_factor('energy', st_pt.chemid, uq_iter, self.pes)
                        freqFactor = uq_obj.calc_factor('freq', st_pt.chemid, uq_iter, self.pes)

                        well_blocks[st_pt.chemid] = self.write_well(st_pt,
                                                                    energyAdd,
                                                                    freqFactor,
                                                                    uq_iter)
                    elif len(reaction.products) == 2:
                        bimol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                        energyAdd = uq_obj.calc_factor('energy', bimol_name, uq_iter, self.pes)
                        freqFactor = uq_obj.calc_factor('freq', bimol_name, uq_iter, self.pes)
                        bimolec_blocks[bimol_name] = self.write_bimol([opt.species for opt in reaction.prod_opt],
                                                                      energyAdd,
                                                                      freqFactor,
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

            # Homolytic scission - barrierless reactions with no TS search
            if self.species.homolytic_scissions is not None:
                for hs in self.species.homolytic_scissions.hss:
                    new = 1
                    if hs.status == -1:
                        hs_prod_name = '_'.join(sorted([str(prod.chemid) for prod in hs.products]))
                        if hs_prod_name not in written_bimol_names and hs_prod_name not in written_termolec_names:
                            if hs_prod_name in self.barrierless_names:
                                new = 0
                            if new:
                                self.barrierless_names[hs_prod_name] = hs_prod_name
                            if new == 1 or uq_iter >= 0:
                                barrierless_energyAdd = uq_obj.calc_factor('energy', hs_prod_name, uq_iter, self.pes)
                                barrierless_freqFactor = uq_obj.calc_factor('freq', hs_prod_name, uq_iter, self.pes)
                                barrierless_blocks[hs_prod_name] = self.write_barrierless([opt.species for opt in hs.prod_opt],
                                                                                          hs,
                                                                                          barrierless_energyAdd,
                                                                                          barrierless_freqFactor,
                                                                                          uq_iter)
       
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

        uq_obj.format_uqtk_data() 

        return 0


    def write_barrierless(self, species_list, reaction, energyAdd, freqFactor, uq_iter):

        if len(reaction.products) == 2:
            barrierless = self.write_bimol(species_list,
                                           energyAdd,
                                           freqFactor,
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


    def write_bimol(self, prod_list, well_add, freqFactor, uq_iter, bless=0):
        """
        Create the block for MESS for a bimolecular product.
        well0: reactant on this PES (zero-energy reference)
        uq_n = number of uncertainty runs
        """

        fragments = ''
        for species in prod_list:
            if species.natom > 1:

                if self.par['pes']:
                    name = '{{fr_name_{}}}'.format(species.chemid)
                    freq = '{freq}'
                else:
                    name = self.fragment_names[species.chemid] + ' ! ' + str(species.chemid)
                    freq = self.make_freq(species, freqFactor, 0)
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
                    name = '{{fr_name_{}}}'.format(species.chemid)
                else:
                    name = self.fragment_names[species.chemid] + ' ! ' + str(species.chemid)

                fragments += self.atomtpl.format(chemid=name,
                                                 element=species.atom[0],
                                                 nelec=1,
                                                 mult=species.mult)

        pr_name = '_'.join(sorted([str(species.chemid) for species in prod_list]))
        if self.par['pes']:
            name = '{name}'
            energy = '{ground_energy}'
        else:
            if bless == 0:
                name = '{} ! {}'.format(self.bimolec_names[pr_name], pr_name)
            elif bless == 1:
                name = '{} ! barrierless'.format(self.barrierless_names[pr_name])

            energy = (sum([sp.energy for sp in prod_list]) + sum([sp.zpe for sp in prod_list]) - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL

            energy += well_add
        

        if bless == 0:
            bimol = self.bimoltpl.format(chemids=name,
                                         fragments=fragments,
                                         ground_energy=energy)

        elif bless == 1:
            values = list(self.barrierless_names.values())
            index = values.index(pr_name)
            bimol = self.blbimoltpl.format(barrier='nobar_{}'.format(index),
                                           reactant=self.well_names[self.species.chemid],
                                           prod=name,
                                           dummy='',
                                           fragments=fragments,
                                           ground_energy=energy)

            with open('{}_{:04d}.mess'.format(pr_name, uq_iter), 'w') as f:
                f.write(bimol)

        return bimol


    def write_well(self, species, well_add, freqFactor, uq_iter):
        """
        Create the block for MESS for a well.
        well0: reactant on this PES (zero-energy reference)
        """

        if self.par['pes']:
            name = '{name}'
            energy = '{zeroenergy}'
            freq = '{freq}'
        else:
            name = self.well_names[species.chemid] + ' ! ' + str(species.chemid)
            energy = ((species.energy + species.zpe) - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL
            energy += well_add
            freq = self.make_freq(species, freqFactor, 0)

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

        with open('{}_{:04d}.mess'.format(species.chemid, uq_iter), 'w') as f:
            f.write(mess_well)

        return mess_well

    def write_barrier(self, reaction, index, left_zeroenergy, right_zeroenergy, barrier_add, freqFactor, imagfreqFactor, uq_iter):
        """
        Create the block for a MESS barrier.
        """

        freq = ''

        if self.par['pes'] == 1:
            imfreq = {'imfreq'}
            cutoff = {'cutoff'}
            welldepth1 = {'welldepth1'}
            welldepth2 = {'welldepth2'}
        else:
            left_zeroenergy += barrier_add
            right_zeroenergy += barrier_add
            cutoff = min(left_zeroenergy, right_zeroenergy)
            imfreq = -reaction.ts.reduced_freqs[0] * imagfreqFactor
            welldepth1 = left_zeroenergy
            welldepth2 = right_zeroenergy
 
            # write tunneling block
        if left_zeroenergy < 0 or right_zeroenergy < 0:
            tun = ''
        else:
            tun = self.tunneltpl.format(cutoff=cutoff,
                                        imfreq=imfreq,
                                        welldepth1=welldepth1,
                                        welldepth2=welldepth2)

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
            freq = '{freq}'
        else:
            name = self.ts_names[reaction.instance_name]
            chemid_reac = self.well_names[self.species.chemid]
            chemid_prod = prod_name
            long_rxn_name = reaction.instance_name
            zeroenergy = left_zeroenergy + barrier_add
            freq = self.make_freq(reaction.ts, freqFactor, 0)
    
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
                                       rotors=self.make_rotors(reaction.ts, norot=self.ts_names[reaction.instance_name]),
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
                                       freq=freq,
                                       rotors=self.make_rotors(reaction.ts, norot=self.ts_names[reaction.instance_name]),
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


        return 0


    def make_geom(self, species):
        geom = ''
        for i, at in enumerate(species.atom):
            geom += '            '
            x, y, z = species.geom[i]
            geom += '{} {:.6f} {:.6f} {:.6f}\n'.format(at, x, y, z)
        return geom


    def make_freq(self, species, factor, wellorts):
        freq = ''
        freqarray = []
        #wellorts: 0 for wells and 1 for saddle points
        if wellorts == 0:
            frequencies = species.reduced_freqs
        else:
            frequencies = species.reduced_freqs[1:]
        for i, fr in enumerate(frequencies):
            fr = fr * factor
            freqarray.append(fr)
            freq += '{:.1f} '.format(fr)
            if i % 3 == 2:
                freq += '\n         '
        return(freq)


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


    def check_running(self, pid):
        devnull = open(os.devnull, 'w')
        if self.par['queuing'] == 'pbs':
            command = 'qstat -f | grep ' + '"Job Id: ' + pid + '"' + ' > /dev/null'
        elif self.par['queuing'] == 'slurm':
            command = 'scontrol show job ' + pid + ' | grep "JobId=' + pid + '"' + ' > /dev/null'

        stat = int(subprocess.call(command, shell=True, stdout=devnull, stderr=devnull))
        return stat
