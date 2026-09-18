import os
import stat
import re
import logging
import json
import numpy as np
import subprocess
import time
import logging
from collections import Counter

from kinbot import kb_path
from kinbot import constants
from kinbot import frequencies
from kinbot.energy import species_zero_k_hartree
from kinbot.uncertaintyAnalysis import UQ

logger = logging.getLogger('KinBot')


def apply_conformer_shifts(contents, ground_min=None):
    """Resolve deferred conformer offsets after PES energy placeholders.

    Older intermediate files also defer imaginary frequencies in comments.
    Consume those comments once. If requested, apply the existing submerged
    barrier correction to each resolved member, before adjusting its depths.
    """
    output = []
    correction = None
    for line in contents.splitlines():
        words = line.split()
        changed = False
        if words and words[0].startswith('ZeroEnergy['):
            correction = None
            shift = float(words[3]) if len(words) == 4 and words[2] == '!' else 0.
            original = float(words[1])
            energy = original + shift
            if ground_min is not None:
                energy = max(energy, ground_min)
            if shift or ground_min is not None or len(words) == 4:
                correction = energy - original
                words = [words[0], str(round(energy, 2))]
                changed = True
        elif 'End ! RRHO' in line:
            correction = None
        elif words and words[0].startswith(('CutoffEnergy[', 'WellDepth[')):
            if correction is not None:
                words[1] = str(round(float(words[1]) + correction, 2))
                changed = True
        elif words and words[0].startswith('ImaginaryFrequency['):
            if len(words) == 4 and words[2] == '!':
                words = [words[0], str(abs(float(words[3])))]
                changed = True
        if changed:
            line = line[:len(line) - len(line.lstrip())] + ' '.join(words)
        output.append(line)
    return '\n'.join(output) + ('\n' if contents.endswith('\n') else '')


def finalize_mc_mess(contents, correct_submerged=False):
    """Resolve MC offsets, then enforce bounds from the serialized endpoints.

    MESS uses the lowest rendered member as a well's ground. This bound is
    distinct from the selected-parent approximation used for Eckart depths.
    Keep that shape convention, adjusting depths only when the barrier itself
    is raised, and omit tunneling at/below the actual connected-well ground.
    """
    contents = apply_conformer_shifts(contents)
    starts = list(re.finditer(r'^\s*(Well|Bimolecular|Barrier)\s+(\S+)[^\n]*',
                              contents, re.MULTILINE))
    blocks = [contents[match.start():starts[i+1].start() if i+1 < len(starts) else len(contents)]
              for i, match in enumerate(starts)]
    grounds, wells = {}, set()
    for match, block in zip(starts, blocks):
        kind, name = match.group(1, 2)
        if kind == 'Barrier':
            continue
        keyword = 'ZeroEnergy' if kind == 'Well' else 'GroundEnergy'
        values = [float(line.split()[1]) for line in block.splitlines()
                  if line.strip().startswith(keyword + '[kcal/mol]')]
        if values:
            grounds[name] = min(values)
        if kind == 'Well':
            wells.add(name)
    for i, (match, block) in enumerate(zip(starts, blocks)):
        if match.group(1) != 'Barrier':
            continue
        endpoints = match.group().split()[2:4]
        floor = max((grounds[name] for name in endpoints if name in grounds), default=None)
        well_floor = max((grounds[name] for name in endpoints if name in wells and name in grounds),
                         default=float('-inf'))
        if correct_submerged and floor is not None:
            block = apply_conformer_shifts(block, ground_min=floor)
        lines = block.splitlines(keepends=True)
        output, energy, j = [], None, 0
        while j < len(lines):
            words = lines[j].split()
            if words and words[0].startswith('ZeroEnergy[kcal/mol]'):
                energy = float(words[1])
            if words[:2] == ['Tunneling', 'Eckart']:
                end = j + 1
                while end < len(lines) and lines[end].split()[:1] != ['End']:
                    end += 1
                parameters = [float(line.split()[1]) for line in lines[j+1:end]
                              if line.strip().startswith(('CutoffEnergy[', 'WellDepth['))]
                if ((energy is not None and energy <= well_floor)
                        or any(value <= 0. for value in parameters)):
                    output.append('! barrier is submerged or has zero serialized tunneling depth\n')
                    j = end + 1
                    continue
            output.append(lines[j])
            j += 1
        blocks[i] = ''.join(output)
    return (contents[:starts[0].start()] + ''.join(blocks)) if starts else contents


class MESS:
    """
    Class that reads and writes MESS files
    UQ analysis parameter (uq) can be used to generate 'n' number of mess input files
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
        with open(f'{kb_path}/tpl/mess_header.tpl') as f:
            self.headertpl = f.read()
        with open(f'{kb_path}/tpl/mess_dummy.tpl') as f:
            self.dummytpl = f.read()
        with open(f'{kb_path}/tpl/mess_termol.tpl') as f:
            self.termoltpl = f.read()
        with open(f'{kb_path}/tpl/mess_fragment.tpl') as f:
            self.fragmenttpl = f.read()
        with open(f'{kb_path}/tpl/mess_fragment_OH.tpl') as f:
            self.fragmenttplOH = f.read()
        with open(f'{kb_path}/tpl/mess_pstfragment.tpl') as f:
            self.pstfragmenttpl = f.read()
        with open(f'{kb_path}/tpl/mess_hinderedrotor.tpl') as f:
            self.hinderedrotortpl = f.read()
        with open(f'{kb_path}/tpl/mess_hinderedrotorgeom.tpl') as f:  # to include geometry when needed
            self.hinderedrotorgeomtpl = f.read()
        with open(f'{kb_path}/tpl/mess_freerotor.tpl') as f:
            self.freerotortpl = f.read()
        with open(f'{kb_path}/tpl/mess_atom.tpl') as f:
            self.atomtpl = f.read()
        with open(f'{kb_path}/tpl/mess_tunneling.tpl') as f:
            self.tunneltpl = f.read()
        with open(f'{kb_path}/tpl/mess_well.tpl') as f:
            self.welltpl = f.read()
        with open(f'{kb_path}/tpl/mess_well_union.tpl') as f:
            self.welluniontpl = f.read()
        with open(f'{kb_path}/tpl/mess_bimol.tpl') as f:
            self.bimoltpl = f.read()
        with open(f'{kb_path}/tpl/mess_barrierless.tpl') as f:
            self.blbimoltpl = f.read()
        with open(f'{kb_path}/tpl/mess_barrier.tpl') as f:
            self.barriertpl = f.read()
        with open(f'{kb_path}/tpl/mess_barrier_union.tpl') as f:
            self.barrieruniontpl = f.read()
        with open(f'{kb_path}/tpl/mess_rrho.tpl') as f:
            self.rrhotpl = f.read()
        with open(f'{kb_path}/tpl/mess_core_rr.tpl') as f:
            self.corerrtpl = f.read()
        with open(f'{kb_path}/tpl/mess_pst.tpl') as f:
            self.psttpl = f.read()
        with open(f'{kb_path}/tpl/mess_variational.tpl') as f:
            self.variationaltpl = f.read()
        with open(f'{kb_path}/tpl/mess_2tst.tpl') as f:
            self.twotstpl = f.read()


    def write_header(self, lot):
        """
        Create the header block for MESS
        """
        if 'prod' in self.species.name:
            reactant = self.species.name
        else:
            reactant = self.species.chemid
        header = self.headertpl.format(LevelOfTheory=lot,
                                       TemperatureList=' '.join([str(ti) for ti in self.par['TemperatureList']]),
                                       PressureList=' '.join([str(pi) for pi in self.par['PressureList']]),
                                       EnergyStepOverTemperature=self.par['EnergyStepOverTemperature'],
                                       ExcessEnergyOverTemperature=self.par['ExcessEnergyOverTemperature'],
                                       ModelEnergyLimit=self.par['ModelEnergyLimit'],
                                       CalculationMethod=self.par['CalculationMethod'],
                                       ChemicalEigenvalueMax=self.par['ChemicalEigenvalueMax'],
                                       Reactant=self.well_names[reactant],
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
                    if reaction.do_vdW:
                        vdW_well = reaction.irc_prod
                        if vdW_well.name not in self.well_names:
                            self.well_names[vdW_well.name] = 'w_{}'.format(len(self.well_names) + 1)
                    for st_pt in reaction.products:
                        if st_pt.chemid not in self.fragment_names:
                            self.fragment_names[st_pt.chemid] = 'fr_{}'.format(len(self.fragment_names) + 1)
                    bimol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                    if bimol_name not in self.bimolec_names:
                        self.bimolec_names[bimol_name] = 'b_{}'.format(len(self.bimolec_names) + 1)
                else:
                    # TERMOLECULAR
                    for st_pt in reaction.products:
                        if st_pt.chemid not in self.fragment_names:
                            self.fragment_names[st_pt.chemid] = 'fr_{}'.format(len(self.fragment_names) + 1)
                    termol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                    if termol_name not in self.termolec_names:
                        self.termolec_names[termol_name] = 't_{}'.format(len(self.termolec_names) + 1)

    def _validate_composite_coverage(self):
        """A direct MESS network needs a complete accepted E0 for each species."""
        species = [self.species]
        for index, reaction in enumerate(self.species.reac_obj):
            if self.species.reac_ts_done[index] != -1:
                continue
            species.append(reaction.ts)
            if reaction.do_vdW:
                species.append(reaction.irc_prod_opt.species)
            species.extend(opt.species for opt in reaction.prod_opt)
        finals = [getattr(item, 'final_zero_k_energy', None) for item in species]
        if not any(final is not None for final in finals):
            return
        if any(final is None for final in finals):
            raise ValueError('MESS needs accepted ANL energies for every well, '
                             'transition state, and product in this network.')
        for item in species:
            species_zero_k_hartree(item)

    def _formation_metadata(self):
        """Keep absolute Hf(0) next to MESS's relative zero-energy input."""
        species = [self.species]
        for index, reaction in enumerate(self.species.reac_obj):
            if self.species.reac_ts_done[index] != -1:
                continue
            if reaction.do_vdW:
                species.append(reaction.irc_prod_opt.species)
            species.extend(opt.species for opt in reaction.prod_opt)
        records = {}
        for item in species:
            formation = getattr(item, 'formation_enthalpy_0k', None)
            if formation is None:
                continue
            final = getattr(item, 'final_zero_k_energy', None)
            from kinbot.anl.atct import _canonical_smiles
            if (final is None
                    or _canonical_smiles(formation.target_smiles) !=
                       _canonical_smiles(final.smiles)
                    or formation.method != final.method
                    or formation.energy_sources.get(formation.target_smiles) != final.source):
                raise ValueError(f'{item.name}: CBH formation provenance does not match E0.')
            key = str(item.chemid)
            record = {'smiles': final.smiles, 'method': formation.method,
                      'zero_k_energy_hartree': final.hartree,
                      'formation_0k_kj_mol': formation.formation_0k_kj_mol,
                      'cbh_rung': formation.rung,
                      'atct_version': formation.atct_version,
                      'atct_source_sha256': formation.atct_source_sha256,
                      'energy_sources': dict(formation.energy_sources),
                      'reference_ids': {smiles: reference.atct_id for smiles, reference
                                        in formation.references.items()}}
            if key in records and records[key] != record:
                raise ValueError(f'{key}: conflicting 0 K formation enthalpies.')
            records[key] = record
        return records

    def write_input(self, qc):
        """
        write the input for all the wells, bimolecular products and barriers
        both in a separate file, as well as in one large ME file
        """
        uq = UQ(self.par)

        self._validate_composite_coverage()
        formation_metadata = self._formation_metadata()

        # create short names for all the species, bimolecular products and barriers
        self.create_short_names()
        final = getattr(self.species, 'final_zero_k_energy', None)
        header = self.write_header(f'{self.par["high_level_method"]}/'
                                   f'{self.par["high_level_basis"]}'
                                   if final is None else 'accepted ANL ladder')

        # filter ts's with the same reactants and products:
        ts_unique = {}  # key: ts name, value: [prod_name, energy]
        ts_all = {}

        for index, reaction in enumerate(self.species.reac_obj):
            if self.species.reac_ts_done[index] == -1:
                rxnProds = []
                if reaction.do_vdW:
                    rxnProds.append(reaction.irc_prod.name)
                else:
                    for x in reaction.products:
                        rxnProds.append(x.chemid)
                rxnProds.sort()
                prod_name = '_'.join([str(pi) for pi in rxnProds])
                new = 1
                remove = []
                ts_all[reaction.instance_name] = [
                    prod_name, species_zero_k_hartree(reaction.ts)]
                for ts in ts_unique:
                    if ts_unique[ts][0] == prod_name:
                        # check for the barrier with the lowest energy  # check if hom_sci is first
                        if (ts_unique[ts][1] > species_zero_k_hartree(reaction.ts)
                                or 'hom_sci' in ts) \
                                and 'hom_sci' not in reaction.instance_name:
                            # remove the current barrier
                            remove.append(ts)
                        else:
                            new = 0
                for ts in remove:
                    ts_unique.pop(ts, None)
                if new:
                    ts_unique[reaction.instance_name] = [
                        prod_name, species_zero_k_hartree(reaction.ts)]

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

            well_energy_add = uq.calc_factor('energy', uq_iter)
            well_freq_factor = uq.calc_factor('freq', uq_iter)
            well_blocks[self.species.chemid] = self.write_well(self.species,
                                                               well_energy_add,
                                                               well_freq_factor,
                                                               uq_iter)
            
            for index, reaction in enumerate(self.species.reac_obj):
                if reaction.instance_name in ts_all: #should be TS unique?
                    barrier_add = uq.calc_factor('barrier', uq_iter)
                    freq_factor = uq.calc_factor('freq', uq_iter)
                    imagfreq_factor = uq.calc_factor('imagfreq', uq_iter)
        
        # get left-right barrier
                    species_zeroenergy = species_zero_k_hartree(self.species) * constants.AUtoKCAL
                    if self.species.reac_ts_done[index] == -1:
                        ts_zeroenergy = species_zero_k_hartree(reaction.ts) * constants.AUtoKCAL
                        if (final is None and not self.par['high_level']
                                and reaction.mp2 == 1 and self.par['qc'] != 'nn_pes'):
                            jobname = '{}_well_mp2'.format(str(self.species.chemid))
                            well_zeroenergy = self.get_zeroenergy(jobname, qc)
                        elif (final is None and not self.par['high_level']
                              and self.species.reac_type[index] == 'barrierless_saddle'):
                            jobname = '{}_well_bls'.format(str(self.species.chemid))
                            well_zeroenergy = self.get_zeroenergy(jobname, qc)
                        else:
                            well_zeroenergy = species_zeroenergy
                        left_zeroenergy = ts_zeroenergy - well_zeroenergy

                        prod_zeroenergy = 0
                        if reaction.do_vdW:
                            prod_zeroenergy += species_zero_k_hartree(
                                reaction.irc_prod_opt.species) * constants.AUtoKCAL
                        else:
                            for opt in reaction.prod_opt:
                                prod_zeroenergy += species_zero_k_hartree(
                                    opt.species) * constants.AUtoKCAL
                        right_zeroenergy = ts_zeroenergy - prod_zeroenergy

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
                    if reaction.do_vdW:
                        st_pt = reaction.irc_prod_opt.species
                        energy_add = uq.calc_factor('energy', uq_iter)
                        freq_factor = uq.calc_factor('freq', uq_iter)
                        well_blocks[st_pt.name] = self.write_well(st_pt,
                                                                    energy_add,
                                                                    freq_factor,
                                                                    uq_iter)
                        bimol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                        energy_add = uq.calc_factor('energy', uq_iter)
                        freq_factor = uq.calc_factor('freq', uq_iter)
                        bless = 1
                        pstsymm_factor = uq.calc_factor('pstsymm', uq_iter)
                        bimolec_blocks[bimol_name] = self.write_bimol([opt.species for opt in reaction.prod_opt],
                                                                      energy_add,
                                                                      freq_factor,
                                                                      pstsymm_factor,
                                                                      uq_iter,
                                                                      bless=bless,
                                                                      vdW=True)
                        written_bimol_names.append(bimol_name)
                    elif len(reaction.products) == 1:
                        st_pt = reaction.prod_opt[0].species
                        energy_add = uq.calc_factor('energy', uq_iter)
                        freq_factor = uq.calc_factor('freq', uq_iter)
                        well_blocks[st_pt.chemid] = self.write_well(st_pt,
                                                                    energy_add,
                                                                    freq_factor,
                                                                    uq_iter)
                    elif len(reaction.products) == 2:
                        bimol_name = '_'.join(sorted([str(st_pt.chemid) for st_pt in reaction.products]))
                        energy_add = uq.calc_factor('energy', uq_iter)
                        freq_factor = uq.calc_factor('freq', uq_iter)
                        if 'hom_sci' not in reaction.instance_name:
                            bless = 0
                            pstsymm_factor = 1
                        else:
                            bless = 1
                            pstsymm_factor = uq.calc_factor('pstsymm', uq_iter)

                        bimolec_blocks[bimol_name] = self.write_bimol([opt.species for opt in reaction.prod_opt],
                                                                      energy_add,
                                                                      freq_factor,
                                                                      pstsymm_factor,
                                                                      uq_iter,
                                                                      bless=bless)
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

            if 'prod' in self.species.name:
                dum = self.dummytpl.format(barrier='tsd', reactant=self.well_names[self.species.name], dummy='d1')
            else:
                dum = self.dummytpl.format(barrier='tsd', reactant=self.well_names[self.species.chemid], dummy='d1')

            mess_iter = "{0:04d}".format(uq_iter)

            with open('me/mess_%s.inp' % mess_iter, 'w') as f_out:
                contents = header + divider + wells + bimols + tss + termols + barrierless + divider + 'End ! end kinetics\n'
                if self.par['multi_conf_tst']:
                    contents = finalize_mc_mess(contents, self.par.get('correct_submerged', 0))
                f_out.write(contents)

        if formation_metadata:
            with open('me/formation_0k.json', 'w') as f_out:
                json.dump({'schema': 1, 'units': 'kJ/mol',
                           'species': formation_metadata}, f_out,
                          indent=2, sort_keys=True)
                f_out.write('\n')

        return 0


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


    def write_bimol(self, prod_list, well_add, freq_factor, pstsymm_factor, uq_iter, bless, vdW=False):
        """
        Create the block for MESS for a bimolecular product.
        In case of a barrierless reaction (bless=1) also add a phase-space theory barrier.
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
                                                           geom=self.rotor_geom(species),
                                                           symm=float(species.sigma_ext) / float(species.nopt),
                                                           rotconst=self.rotor_core_line(species),
                                                           freq=self.make_freq(species.reduced_freqs, freq_factor, 0))
                else:
                    fragments += self.fragmenttpl.format(chemid=name,
                                                         smi=species.smiles,
                                                         natom=species.natom,
                                                         geom=self.rotor_geom(species),
                                                         symm=float(species.sigma_ext) / float(species.nopt),
                                                         rotconst=self.rotor_core_line(species),
                                                         nfreq=len(species.reduced_freqs),
                                                         freq=self.make_freq(species.reduced_freqs, freq_factor, 0),
                                                         hinderedrotor=self.make_rotors(species, freq_factor),
                                                         nelec=1,
                                                         mult=species.mult)
                if bless == 1:
                    tot_nfreq += len(species.reduced_freqs)
                    combined_freq += self.make_freq(species.reduced_freqs, freq_factor, 0)
                    combined_hir += self.make_rotors(species, freq_factor, bless=True)

                    if nsp == 0: 
                        frag1 = self.pstfragmenttpl.format(chemid=name,
                                                           smi=species.smiles,
                                                           natom=species.natom,
                                                           geom=self.rotor_geom(species))
                    if nsp == 1: 
                        frag2 = self.pstfragmenttpl.format(chemid=name,
                                                           smi=species.smiles,
                                                           natom=species.natom,
                                                           geom=self.rotor_geom(species))
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
                                                           geom=self.rotor_geom(species))
                    if nsp == 1: 
                        frag2 = self.pstfragmenttpl.format(chemid=name,
                                                           smi=species.smiles,
                                                           natom=species.natom,
                                                           geom=self.rotor_geom(species))
 

        pr_name = '_'.join(sorted([str(species.chemid) for species in prod_list]))
        if self.par['pes']:
            name = '{{name}} ! {} {}'.format(smi[0], smi[1])
            energy = '{ground_energy}'
        else:
            name = '{} ! {}'.format(self.bimolec_names[pr_name], pr_name)
            energy = (sum(species_zero_k_hartree(sp) for sp in prod_list)
                      - species_zero_k_hartree(self.species)) * constants.AUtoKCAL
            energy += well_add
            energy = round(energy, 2)
        
        if bless == 0:
            bimol = self.bimoltpl.format(chemids=name,
                                         smi=species.smiles,
                                         fragments=fragments,
                                         ground_energy=energy)

        elif bless == 1:
            stoich = ''
            el_counter = Counter(self.species.atom)
            for el in constants.elements:
                if el_counter[el]:
                    stoich += '{}{}'.format(el, el_counter[el])
            bimol = self.blbimoltpl.format(barrier='{blessname}',
                                           reactant='{wellname}',
                                           prod='{prodname}',
                                           chemids=name,
                                           pstsymm=pstsymm_factor,
                                           stoich=stoich,
                                           frag1=frag1,
                                           frag2=frag2,
                                           nfreq=tot_nfreq,
                                           freq=combined_freq,
                                           mult=self.species.mult,
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
            if 'prod' not in species.name:
                norot = None
            else:
                norot = str(species.name)
        else:
            if 'prod' not in species.name:
                name = self.well_names[species.chemid] + ' ! ' + str(species.chemid)
                norot = None
            else:
                name = self.well_names[species.name] + ' ! ' + str(species.name)
                norot = str(species.name)
            zeroenergy = (species_zero_k_hartree(species) -
                          species_zero_k_hartree(self.species)) * constants.AUtoKCAL
            zeroenergy += well_add
            zeroenergy = round(zeroenergy, 2)

        valid_conformers = [ci for ci, co in enumerate(species.conformer_index)
                            if co >= 0]
        nunq_confs = len(valid_conformers)

        if not self.par['multi_conf_tst'] or not valid_conformers:
            mess_well = self.welltpl.format(chemid=name,
                                            smi=species.smiles,
                                            natom=species.natom,
                                            geom=self.rotor_geom(species),
                                            symm=float(species.sigma_ext) / float(species.nopt),
                                            rotconst=self.rotor_core_line(species),
                                            nfreq=len(species.reduced_freqs),
                                            freq=self.make_freq(species.reduced_freqs, freq_factor, 0),
                                            hinderedrotor=self.make_rotors(species, freq_factor, norot=norot),
                                            nelec=1,
                                            mult=species.mult,
                                            zeroenergy=zeroenergy)
        else:
            rrho = '      '
            base_zeroen = species.energy + species.zpe
            for ci in valid_conformers:
                shift = constants.AUtoKCAL * (species.conformer_zeroenergy[ci] - base_zeroen)
                conformer_freq = frequencies.thermochemical_frequencies(
                    species.conformer_freq[ci], 0, self.par.get('imagfreq_threshold', 50.))
                conformer_zeroenergy = (zeroenergy if self.par['pes'] else
                                        round(zeroenergy + shift, 2))
                corerr = self.corerrtpl.format(symm=float(species.sigma_ext) / float(species.nopt),
                                               rotconst=self.rotor_core_line(species, species.conformer_geom[ci], species.conformer_freq[ci]))
                rrho += self.rrhotpl.format(natom=species.natom,
                                            geom=self.rotor_geom(species, species.conformer_geom[ci], species.conformer_freq[ci]),
                                            core=corerr,
                                            nfreq=len(species.conformer_freq[ci]),
                                            freq=self.make_freq(conformer_freq, freq_factor, 0),
                                            rotors='',
                                            tunneling='',
                                            nelec=1,
                                            mult=species.mult,
                                            zeroenergy=conformer_zeroenergy,
                                            shift=shift if self.par['pes'] else '',
                                           )
            rrho = '      '.join(rrho.splitlines(True))  # indent
            mess_well = self.welluniontpl.format(chemid=name,
                                                 smi=species.smiles,
                                                 nunion=nunq_confs,
                                                 rrho=rrho)
        if 'prod' not in species.name:
            with open('{}_{:04d}.mess'.format(species.chemid, uq_iter), 'w') as f:
                f.write(mess_well)
        else:
            with open('{}_{:04d}.mess'.format(species.name, uq_iter), 'w') as f:
                f.write(mess_well)

        return mess_well

    def write_barrier(self, reaction, index, left_zeroenergy, right_zeroenergy, barrier_add, freq_factor, imagfreq_factor, uq_iter):
        """
        Create the block for a MESS barrier.
        """

        left_zeroenergy += barrier_add
        right_zeroenergy += barrier_add

        valid_conformers = [ci for ci, co in enumerate(reaction.ts.conformer_index)
                            if co >= 0]
        nunq_confs = len(valid_conformers)

        # write tunneling block
        if (not self.par['pes']
                and min(round(left_zeroenergy, 2), round(right_zeroenergy, 2)) <= 0):
            tun = f'! barrier is submerged {left_zeroenergy} {right_zeroenergy}'
        elif self.par['pes'] == 0:
            tun = self.tunneltpl.format(cutoff=round(min(left_zeroenergy, right_zeroenergy), 2),
                                        imfreq=round(-reaction.ts.reduced_freqs[0] * imagfreq_factor, 2),
                                        welldepth1=round(left_zeroenergy, 2),
                                        welldepth2=round(right_zeroenergy, 2))
        else: 
            logger.info(f'Barrier {reaction.instance_name}')
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
    
        if self.species.reac_type[index] == 'barrierless_saddle':
            freq = self.make_freq(reaction.prod_opt[0].species.reduced_freqs, freq_factor, 0) + \
                   self.make_freq(reaction.prod_opt[1].species.reduced_freqs, freq_factor, 0) 
            rotors = self.make_rotors(reaction.prod_opt[0].species, freq_factor) + \
                     self.make_rotors(reaction.prod_opt[1].species, freq_factor) 
            nfreq = len(reaction.prod_opt[0].species.reduced_freqs) + \
                    len(reaction.prod_opt[1].species.reduced_freqs)
            if self.par['pes']:
                prodzeroenergy = '{prodzeroenergy}'
            else:
                prodzeroenergy = (
                    species_zero_k_hartree(reaction.prod_opt[0].species) +
                    species_zero_k_hartree(reaction.prod_opt[1].species) -
                    species_zero_k_hartree(self.species)) * constants.AUtoKCAL

            outerts = self.psttpl.format(natom1=reaction.prod_opt[0].species.natom,
                                         geom1=self.rotor_geom(reaction.prod_opt[0].species),
                                         natom2=reaction.prod_opt[1].species.natom,
                                         geom2=self.rotor_geom(reaction.prod_opt[1].species),
                                         symm=float(reaction.ts.sigma_ext) / float(reaction.ts.nopt),
                                         prefact='prefactor',
                                         exponent=6,
                                         nfreq=nfreq,
                                         freq=freq,
                                         hinderedrotor=rotors,
                                         nelec=1,
                                         mult=reaction.ts.mult,
                                         prodzeroenergy=prodzeroenergy
                                         )
            twotst = self.twotstpl.format(outerts=outerts)
            corerr = self.corerrtpl.format(symm=float(reaction.ts.sigma_ext) / float(reaction.ts.nopt),
                                           rotconst=self.rotor_core_line(reaction.ts))
            rrho = self.rrhotpl.format(natom=reaction.ts.natom,
                                       geom=self.rotor_geom(reaction.ts),
                                       core=corerr,
                                       nfreq=len(reaction.ts.reduced_freqs)-1,
                                       freq=self.make_freq(reaction.ts.reduced_freqs, freq_factor, 1),
                                       rotors=self.make_rotors(reaction.ts, freq_factor, norot=self.ts_names[reaction.instance_name]),
                                       tunneling='',
                                       nelec=1,
                                       mult=reaction.ts.mult,
                                       zeroenergy=zeroenergy,
                                       shift='',
                                      )
            variational = self.variationaltpl.format(twotst=twotst,
                                                     variationalmodel=rrho,
                                                     tunneling=tun)
            mess_barrier = self.barriertpl.format(rxn_name=name,
                                                  chemid_reac=chemid_reac,
                                                  chemid_prod=chemid_prod,
                                                  long_rxn_name=long_rxn_name,
                                                  model=variational)
        elif not self.par['multi_conf_tst'] or not valid_conformers:
            corerr = self.corerrtpl.format(symm=float(reaction.ts.sigma_ext) / float(reaction.ts.nopt),
                                           rotconst=self.rotor_core_line(reaction.ts))
            rrho = self.rrhotpl.format(natom=reaction.ts.natom,
                                       geom=self.rotor_geom(reaction.ts),
                                       core=corerr,
                                       nfreq=len(reaction.ts.reduced_freqs)-1,
                                       freq=self.make_freq(reaction.ts.reduced_freqs, freq_factor, 1),
                                       rotors=self.make_rotors(reaction.ts, freq_factor, norot=self.ts_names[reaction.instance_name]),
                                       tunneling=tun,
                                       nelec=1,
                                       mult=reaction.ts.mult,
                                       zeroenergy=zeroenergy,
                                       shift='',
                                      )
            mess_barrier = self.barriertpl.format(rxn_name=name,
                                                  chemid_reac=chemid_reac,
                                                  chemid_prod=chemid_prod,
                                                  long_rxn_name=long_rxn_name,
                                                  model=rrho)
        else:
            rrho = '      '
            base_zeroen = reaction.ts.energy + reaction.ts.zpe
            for ci in valid_conformers:
                shift = constants.AUtoKCAL * (reaction.ts.conformer_zeroenergy[ci] - base_zeroen)
                conformer_freq = frequencies.thermochemical_frequencies(
                    reaction.ts.conformer_freq[ci], 1, self.par.get('imagfreq_threshold', 50.))
                conformer_zeroenergy = (zeroenergy if self.par['pes'] else
                                        round(zeroenergy + shift, 2))
                corerr = self.corerrtpl.format(symm=float(reaction.ts.sigma_ext) / float(reaction.ts.nopt),
                                               rotconst=self.rotor_core_line(reaction.ts, reaction.ts.conformer_geom[ci], reaction.ts.conformer_freq[ci]))
                imfreq = round(-reaction.ts.conformer_freq[ci][0] * imagfreq_factor, 2)
                if self.par['pes']:
                    tun_conf = self.tunneltpl.format(
                        cutoff='{cutoff}', imfreq=imfreq,
                        welldepth1='{welldepth1}', welldepth2='{welldepth2}')
                else:
                    left, right = left_zeroenergy + shift, right_zeroenergy + shift
                    if min(round(left, 2), round(right, 2)) <= 0:
                        tun_conf = f'! barrier is submerged {left} {right}'
                    else:
                        tun_conf = self.tunneltpl.format(
                            cutoff=round(min(left, right), 2), imfreq=imfreq,
                            welldepth1=round(left, 2), welldepth2=round(right, 2))
                rrho += self.rrhotpl.format(natom=reaction.ts.natom,
                                            geom=self.rotor_geom(reaction.ts, reaction.ts.conformer_geom[ci], reaction.ts.conformer_freq[ci]),
                                            core=corerr,
                                            nfreq=len(reaction.ts.conformer_freq[ci])-1,
                                            freq=self.make_freq(conformer_freq, freq_factor, 1),
                                            rotors='',
                                            tunneling=tun_conf,
                                            nelec=1,
                                            mult=reaction.ts.mult,
                                            zeroenergy=conformer_zeroenergy,
                                            shift=shift if self.par['pes'] else '',
                                           )  
            rrho = '      '.join(rrho.splitlines(True))  # indent
            mess_barrier = self.barrieruniontpl.format(rxn_name=name,
                                                       chemid_reac=chemid_reac,
                                                       chemid_prod=chemid_prod,
                                                       long_rxn_name=long_rxn_name,
                                                       nunion=nunq_confs,
                                                       model=rrho)


        with open('{}_{:04d}.mess'.format(reaction.instance_name, uq_iter), 'w') as f:
            f.write(mess_barrier)

    
        return mess_barrier, zeroenergy 

    def run(self):
        """
        Submit the pbs/slurm file to the queue
        wait for the mess run to finish
        """

        pids = []  # list of job pids

        batch_list = ''
        uq_iter = 0
        pid_stats = []
        for uq_iter in range(self.par['uq_n']):
            if self.par["queuing"] == 'local':
                return 0
            submitscript = f'me/run_mess_{str(uq_iter).zfill(4)}{constants.qext[self.par["queuing"]]}'
            self.write_submitscript(submitscript, uq_iter)
            batch_list += f'sbatch {submitscript}\n'
            if not self.par['run_me']:
                continue
            while len(pids) > self.par['uq_max_runs']:
                time.sleep(5)
                for pid in pids:
                    stati = self.check_running(pid)
                    if stati == 0:
                        pids.remove(pid)
                        pid_stats.append(stati)
            
            pid = self.submit(submitscript)
            pids.append(pid)

            if self.par['uq_n'] < self.par['uq_max_runs']:
                stati = 1
                while stati != 0:
                    stati = self.check_running(pid)
                    time.sleep(5)
                pid_stats.append(stati)
        
        # for manual submission
        batch_me = 'batch_me.sub'
        with open(batch_me, 'w') as f:
            f.write(batch_list)
        os.chmod(batch_me, stat.S_IRWXU)

        if all(stati == 0 for s in pid_stats):
            return 0

    def write_submitscript(self, submitscript, uq_iter):
        """
        write a pbs or slurm file for the me/all.inp mess input file
        """

        if self.par['queue_template'] == '':
            q_file = f'{kb_path}/tpl/{self.par["queuing"]}.tpl'
        else:
            q_file = self.par['queue_template']
        with open(q_file) as f:
            tpl_head = f.read()

        q_file = f'{kb_path}/tpl/{self.par["queuing"]}_mess_uq.tpl'
        with open(q_file) as f:
            tpl = f.read()

        with open(submitscript, 'w') as f:
            mess_iter = "{0:04d}".format(uq_iter)
            if self.par['queue_template'] == '':
                if self.par['queuing'] == 'pbs':
                    f.write((tpl_head).format(name='mess', ppn=self.par['ppn'], queue_name=self.par['queue_name'], errdir='me'))
                    f.write((tpl).format(n=mess_iter))
                elif self.par['queuing'] == 'slurm':
                    f.write((tpl_head).format(name='mess', ppn=self.par['ppn'], queue_name=self.par['queue_name'], errdir='me', slurm_feature=self.par['slurm_feature']))
                    f.write((tpl).format(n=mess_iter))
            else:
                f.write((tpl_head).format(name='mess', ppn=self.par['ppn'], queue_name=self.par['queue_name'], errdir='me', slurm_feature=self.par['slurm_feature']))
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

    def rotor_geom(self, species, geom=None, freqs=None):
        """Geometry block for a rigid-rotor species.

        A species that carries 3N-5 frequencies was judged linear when its
        frequencies were derived (see frequencies.assess_linearity). Write a
        linearised copy so that MESS's own test on the moments of inertia
        agrees, and say so in the input and in the log. The stored geometry is
        never modified, and other species are written as calculated.
        """
        geom = species.geom if geom is None else geom
        if freqs is None:
            freqs = getattr(species, 'reduced_freqs', None) or species.freq
        natom = len(species.atom)
        block = self.make_geom(geom, species.atom)
        if natom < 3 or freqs is None or len(freqs) != 3 * natom - 5:
            return block
        deviation = frequencies.max_bend_deviation(geom, species.bond)
        if deviation > frequencies.LINEAR_ANGLE_TOLERANCE:
            logger.warning(f'{species.name}: carries 3N-5 frequencies but its geometry deviates '
                           f'{deviation:.1f} deg from linear; written to MESS as calculated.')
            return block + (f'\n        ! WARNING: 3N-5 frequencies but the geometry deviates '
                            f'{deviation:.1f} deg from linear; written as calculated')
        logger.info(f'{species.name}: written to MESS as an exactly linear rotor '
                    f'(max bend deviation {deviation:.2f} deg).')
        return (self.make_geom(frequencies.linearize(geom, species.atom), species.atom)
                + f'\n        ! geometry linearised for the rigid-rotor model '
                  f'(max bend deviation {deviation:.2f} deg, 3N-5 frequencies)')

    def rotor_core_line(self, species, geom=None, freqs=None):
        """Extra line for the RigidRotor core, normally empty.

        The reverse mismatch to rotor_geom: KinBot judged the species
        non-linear (3N-6 frequencies, a genuinely bent minimum) but its
        geometry is so close to linear that MESS's own test (I_min/I_mid <
        1e-5) would make it a linear rotor, losing one degree of freedom.
        Give MESS the three rotational constants explicitly so it keeps the
        3D rotor KinBot's frequency count assumes, and say so.
        """
        geom = species.geom if geom is None else geom
        if freqs is None:
            freqs = getattr(species, 'reduced_freqs', None) or species.freq
        natom = len(species.atom)
        if natom < 3 or freqs is None or len(freqs) != 3 * natom - 6:
            return ''
        ratio = frequencies.moment_ratio(geom, species.atom)
        if ratio >= frequencies.MESS_LINEAR_MOMENT_RATIO:
            return ''
        constants_cm = frequencies.rotational_constants(geom, species.atom)
        logger.warning(f'{species.name}: quasi-linear species (I_min/I_mid = {ratio:.1e}) with '
                       '3N-6 frequencies; explicit rotational constants written so MESS keeps a '
                       'non-linear rotor. RRHO is unreliable for this species.')
        return ('          RotationalConstants[1/cm]   '
                + ' '.join(f'{b:.6f}' for b in constants_cm)
                + f'\n          ! quasi-linear species (I_min/I_mid = {ratio:.1e}): constants '
                  'given explicitly to keep a non-linear rotor consistent with 3N-6 frequencies')

    def make_geom(self, g, a):
        geom = ''
        for i, at in enumerate(a):
            x, y, z = g[i]
            geom += '        {} {:.6f} {:.6f} {:.6f}\n'.format(at, x, y, z)
        return geom[:-1]

    def scale_freq(self, fr, factor):
        """
        Apply the UQ factor to a single frequency, amplified at low
        frequencies and dampened at high ones. See make_freq.
        """
        exponent = min(self.par['freq_uq_ref'] / fr,
                       self.par['freq_uq_max_exp'])
        return fr * factor ** exponent

    def make_freq(self, fr, factor, wellorts):
        """
        Frequencies are scaled with factor in UQ.
        At freq_uq_ref the scaling is applied as is.
        For lower frequencies the scaling is amplified.
        For higher frequencies the scaling is dampened.

        The scaling is the power law v' = v * factor**(freq_uq_ref / v). The
        UQ factor is sampled log-uniformly, so factor and 1/factor are equally
        likely and have to cancel exactly; the power law guarantees that at
        every frequency, whereas a shift or a division does not. The exponent
        is capped at freq_uq_max_exp so that the amplification stays bounded
        as v approaches zero.
        """
        freq = '        '
        #wellorts: 0 for wells and 1 for saddle points
        if wellorts == 0:
            frequencies = fr
        else:
            frequencies = fr[1:]
        for i, fr in enumerate(frequencies):
            freq += '{:.1f} '.format(self.scale_freq(fr, factor))
            if i % 3 == 2:
                freq += '\n        '
        return(freq[:-1])

    def make_rotorpot(self, species, i, rot, freq_factor):
        rotortype = 'hindered'
        rotorsymm = self.rotorsymm(species, rot)
        ens = species.hir.hir_energies[i]
        rotorpot_num = [(ei - ens[0]) * constants.AUtoKCAL for ei in ens]
        maxen = max(rotorpot_num)
        count = self.nrotorpot(species, rot)
        if species.hir.nrotation % rotorsymm or species.hir.nrotation // rotorsymm < 3:
            # MESS interprets these as equally spaced points in one symmetry
            # period, excluding its repeated endpoint. Interpolate on that grid.
            rotorpot_num = [species.hir.get_fit_value(
                2 * np.pi * point / (rotorsymm * count), rotor=i)
                for point in range(count)]
        else:
            rotorpot_num = rotorpot_num[:count]
        rotorpot = '        ' + ' '.join(f'{freq_factor * energy:.2f}' for energy in rotorpot_num)
        if maxen < self.par['free_rotor_thrs']:
            rotortype = 'free'
        return rotorpot, rotortype

    def rotorsymm(self, species, rot):
        return species.sigma_int[rot[1]][rot[2]]

    def nrotorpot(self, species, rot): 
        rotorsymm = self.rotorsymm(species, rot)
        return max(3, int(np.ceil(species.hir.nrotation / rotorsymm)))

    def make_rotors(self, species, freq_factor, norot=None, bless=False):
        rotors = []
        if self.par['rotor_scan']:
            for i, rot in enumerate(species.dihed):
                if norot is not None:
                    if frequencies.skip_rotor(norot, rot) == 1:
                        continue
                if species.hir is None:
                    if (getattr(species, 'rotor_projection', None) or {}).get('method') != 'harmonic_fallback':
                        # Species optimized with just_high never ran a scan.
                        continue
                    why = 'no usable selected-geometry Hessian; hindered rotors omitted'
                else:
                    why = species.hir.invalid_rotor_reason(i)
                if why is not None:
                    # Leave a trace in the MESS input: this torsion stays a
                    # harmonic oscillator in the frequency list above.
                    rotors.append(f'      ! Rotor about atoms {rot[1] + 1}-{rot[2] + 1} '
                                  f'kept as a harmonic oscillator: {why}')
                    continue
                rotorpot, rotortype = self.make_rotorpot(species, i, rot, freq_factor)
                if rotortype == 'hindered' and bless == False:
                    rotors.append(self.hinderedrotortpl.format(group=' '.join([str(pi + 1) for pi in frequencies.partition(species, rot, species.natom)[0][1:]]),
                                                               axis='{} {}'.format(str(rot[1] + 1), str(rot[2] + 1)),
                                                               rotorsymm=self.rotorsymm(species, rot),
                                                               nrotorpot=self.nrotorpot(species, rot),
                                                               rotorpot=rotorpot))
                if rotortype == 'hindered' and bless == True:
                    rotors.append(self.hinderedrotorgeomtpl.format(geom=self.make_geom(species.geom, species.atom),
                                                                   natom=species.natom,
                                                                   group=' '.join([str(pi + 1) for pi in frequencies.partition(species, rot, species.natom)[0][1:]]),
                                                                   axis='{} {}'.format(str(rot[1] + 1), str(rot[2] + 1)),
                                                                   rotorsymm=self.rotorsymm(species, rot),
                                                                   nrotorpot=self.nrotorpot(species, rot),
                                                                   rotorpot=rotorpot))
                elif rotortype == 'free':
                    rotors.append(self.freerotortpl.format(geom=self.make_geom(species.geom, species.atom),
                                                           natom=species.natom,
                                                           rotorsymm=self.rotorsymm(species, rot),
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
