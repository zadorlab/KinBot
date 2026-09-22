from kinbot.species_routing import routing_key, routing_name, mess_filename
import os
import re
import logging
import numpy as np
import logging
from collections import Counter
from itertools import product

from kinbot import kb_path
from kinbot import constants
from kinbot import frequencies
from kinbot.reaction_path import reaction_path_id, compare_pathways, reject_invalid_pathway
from kinbot.product_complex import reassess_product_complex
from kinbot.uncertaintyAnalysis import UQ

from kinbot.mess_mirrors import annotate_population, annotate_endpoints, complete_mirror_channels

logger = logging.getLogger('KinBot')


def union_stereochemical_barriers(blocks):
    """Sum distinct routes while emitting one MESS Barrier per endpoint pair.

    Each supplied block already contains its own resolved energy, rotors and
    tunneling model. A nested Union preserves a route's MC ensemble, if any.
    This does not enable MC-TST or manufacture degeneracy factors.
    """
    groups = {}
    for block in blocks:
        match = re.search(r'^[ \t]*Barrier[ \t]+[^\n]+\n', block, re.MULTILINE)
        if match is None or len(match.group().split('!')[0].split()) != 4:
            raise ValueError('Expected one resolved MESS barrier with two endpoints.')
        endpoints = tuple(sorted(match.group().split()[2:4]))
        groups.setdefault(endpoints, []).append((block, match))
    result = []
    for entries in groups.values():
        if len(entries) == 1:
            result.append(entries[0][0])
            continue
        keys = []
        for block, _ in entries:
            labels = re.findall(r'^! kinbot_stereopath (\S+)$', block, re.MULTILINE)
            if len(labels) != 1 or labels[0] in keys:
                raise ValueError('Parallel MESS barriers require distinct classified '
                                 'stereochemical routes; do not sum duplicate or unclassified routes.')
            keys.append(labels[0])
        body = ''.join('! pathway source: ' + match.group().strip() + '\n'
                       + block[:match.start()] + block[match.end():] + '\n'
                       for block, match in entries)
        result.append(entries[0][1].group()
                      + f'    Union ! {len(entries)} stereochemical pathways\n'
                      + body + '    End ! stereochemical pathways\n')
    return result


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
        elif words and words[0].startswith('GroundEnergy['):
            correction = None
            if len(words) == 4 and words[2] == '!':
                words = [words[0], str(round(float(words[1]) + float(words[3]), 2))]
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


def validate_mess_populations(contents):
    """Do not represent a full racemate twice under separately named mirrors."""
    # Separately named R/S wells must not each represent the same full racemate.
    # These comments survive deferred PES assembly as well as direct rendering.
    populations = {}
    for key, family in re.findall(r'^! kinbot_racemic_population (\S+) (\S+)$',
                                  contents, re.MULTILINE):
        if family in populations and populations[family] != key:
            raise ValueError('Overlapping racemic populations in MESS: '
                             f'{populations[family]} and {key}. Use one declared '
                             'racemic population or separate specified stereoisomers.')
        populations[family] = key


def finalize_mc_mess(contents, correct_submerged=False):
    """Resolve MC offsets, then enforce bounds from the serialized endpoints.

    MESS uses the lowest rendered member as a well's ground. This bound is
    distinct from the selected-parent approximation used for Eckart depths.
    Keep that shape convention, adjusting depths only when the barrier itself
    is raised, and omit tunneling at/below the actual connected-well ground.
    """
    contents = apply_conformer_shifts(contents)
    validate_mess_populations(contents)
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
from kinbot.conformer_counting import writer_members, representative_record
from kinbot.stereo_identity import canonical_identity, optical_scope
from kinbot.stereo_routing import refuse_routing


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
        with open(f'{kb_path}/tpl/mess_pst_rrho.tpl') as f:
            self.pstrrhotpl = f.read()
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
            reactant = routing_key(self.species)
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
        self.well_names[routing_key(self.species)] = 'w_1'
        for index, reaction in enumerate(self.species.reac_obj):
            if self.species.reac_ts_done[index] == -1:
                self.ts_names[reaction.instance_name] = 'ts_{}'.format(len(self.ts_names) + 1)
                if len(reaction.products) == 1:
                    st_pt = reaction.products[0]
                    if routing_key(st_pt) not in self.well_names:
                        self.well_names[routing_key(st_pt)] = 'w_{}'.format(len(self.well_names) + 1)
                elif len(reaction.products) == 2:
                    if self.par['pes'] and reaction.do_vdW:
                        vdW_well = reaction.irc_prod
                        if vdW_well.name not in self.well_names:
                            self.well_names[vdW_well.name] = 'w_{}'.format(len(self.well_names) + 1)
                    for st_pt in reaction.products:
                        if routing_key(st_pt) not in self.fragment_names:
                            self.fragment_names[routing_key(st_pt)] = 'fr_{}'.format(len(self.fragment_names) + 1)
                    bimol_name = '_'.join(sorted([routing_name(st_pt) for st_pt in reaction.products]))
                    if bimol_name not in self.bimolec_names:
                        self.bimolec_names[bimol_name] = 'b_{}'.format(len(self.bimolec_names) + 1)
                else:
                    # TERMOLECULAR
                    for st_pt in reaction.products:
                        if routing_key(st_pt) not in self.fragment_names:
                            self.fragment_names[routing_key(st_pt)] = 'fr_{}'.format(len(self.fragment_names) + 1)
                    termol_name = '_'.join(sorted([routing_name(st_pt) for st_pt in reaction.products]))
                    if termol_name not in self.termolec_names:
                        self.termolec_names[termol_name] = 't_{}'.format(len(self.termolec_names) + 1)

    def calculation_label(self):
        """Describe the calculation level used by direct and PES MESS output."""
        if self.par['qc'] == 'fc':
            return f"FairChem {self.par['fc_model_path']} ({self.par['fc_task_name']})"
        if self.par['qc'] == 'nn_pes':
            return f"nn_pes {self.par['nn_model']}"
        prefix = 'high_level_' if self.par['high_level'] else ''
        return f"{self.par[prefix + 'method']}/{self.par[prefix + 'basis']}"

    def write_input(self, qc):
        """
        write the input for all the wells, bimolecular products and barriers
        both in a separate file, as well as in one large ME file
        """
        uq = UQ(self.par)
        self.mess_jobs = []

        for index, reaction in enumerate(self.species.reac_obj):
            reject_invalid_pathway(self.species, index, self.par)
            if self.par['pes'] and self.species.reac_ts_done[index] == -1:
                reassess_product_complex(reaction, self.par)
        # create short names for all the species, bimolecular products and barriers
        self.create_short_names()
        header = self.write_header(self.calculation_label())
        for reaction in self.species.reac_obj:
            reason = getattr(reaction, 'stereochemical_rejection', None)
            if reason:
                header += (f'! WARNING: omitted channel {reaction.instance_name}: {reason}. '
                           'Reaction network is incomplete.\n')

        # filter ts's with the same reactants and products:
        ts_unique = {}  # key: ts name, value: [prod_name, energy, stereochemical path]
        ts_all = {}

        for index, reaction in enumerate(self.species.reac_obj):
            if self.species.reac_ts_done[index] == -1:
                # IRC filenames label observations, not distinct products.
                rxnProds = [routing_name(x) for x in reaction.products]
                rxnProds.sort()
                prod_name = '_'.join([str(pi) for pi in rxnProds])
                path_id = reaction_path_id(reaction)
                new = 1
                remove = []
                ts_all[reaction.instance_name] = [prod_name, reaction.ts.energy
                                                + reaction.ts.zpe]
                for ts in ts_unique:
                    if ts_unique[ts][0] == prod_name:
                        decision = compare_pathways(ts, ts_unique[ts][2], ts_unique[ts][1],
                            reaction.instance_name, path_id, reaction.ts.energy + reaction.ts.zpe)
                        if decision == 'replace':
                            remove.append(ts)
                        elif decision == 'keep':
                            new = 0
                for ts in remove:
                    ts_unique.pop(ts, None)
                if new:
                    ts_unique[reaction.instance_name] = [prod_name, reaction.ts.energy + reaction.ts.zpe, path_id]

        if not self.par['multi_conf_tst']:
            from kinbot.hindered_rotors import recover_hir_model
            states = [self.species]
            for reaction in self.species.reac_obj:
                if reaction.instance_name in ts_all:
                    states.extend(opt.species for opt in getattr(reaction, 'prod_opt', ()))
                    if 'hom_sci' not in reaction.instance_name:
                        states.append(reaction.ts)
                    if self.par['pes'] and getattr(reaction, 'do_vdW', False):
                        states.append(reaction.irc_prod_opt.species)
            checked = set()
            for state in states:
                if id(state) not in checked and getattr(state, 'hir', None) is not None:
                    checked.add(id(state))
                    recover_hir_model(state, qc, self.par)
        self._check_optical_models(ts_all, ts_unique)
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
            well_blocks[routing_key(self.species)] = self.write_well(self.species,
                                                               well_energy_add,
                                                               well_freq_factor,
                                                               uq_iter)
            
            for index, reaction in enumerate(self.species.reac_obj):
                # A homolytic scission has no optimized saddle. Its barrier is
                # the existing product-based phase-space model below; never
                # interpret the copied parent as a TS (including its optics).
                if (reaction.instance_name in ts_all
                        and (self.par['pes'] or reaction.instance_name in ts_unique)
                        and self.species.reac_type[index] != 'hom_sci'):
                    barrier_add = uq.calc_factor('barrier', uq_iter)
                    freq_factor = uq.calc_factor('freq', uq_iter)
                    imagfreq_factor = uq.calc_factor('imagfreq', uq_iter)
        
        # get left-right barrier
                    species_zeroenergy = (self.species.energy + self.species.zpe) * constants.AUtoKCAL
                    if self.species.reac_ts_done[index] == -1:
                        ts_zeroenergy = (reaction.ts.energy + reaction.ts.zpe) * constants.AUtoKCAL
                        left_reference_job = getattr(self.species, 'source_job', None)
                        if not self.par['high_level'] and reaction.mp2 == 1 and self.par['qc'] != 'nn_pes':
                            jobname = '{}_well_mp2'.format(routing_name(self.species))
                            well_zeroenergy = self.get_zeroenergy(jobname, qc)
                            left_reference_job = jobname
                        elif not self.par['high_level'] and self.species.reac_type[index] == 'barrierless_saddle':
                            jobname = '{}_well_bls'.format(routing_name(self.species))
                            well_zeroenergy = self.get_zeroenergy(jobname, qc)
                            left_reference_job = jobname
                        else:
                            well_zeroenergy = species_zeroenergy
                        left_zeroenergy = ts_zeroenergy - well_zeroenergy

                        prod_zeroenergy = 0
                        # Retain master's inner-barrier tunneling reference,
                        # even when direct kinetics omits the fast complex exit.
                        if reaction.do_vdW:
                            complex_species = reaction.irc_prod_opt.species
                            prod_zeroenergy += (complex_species.energy + complex_species.zpe) * constants.AUtoKCAL
                        else:
                            for opt in reaction.prod_opt:
                                prod_zeroenergy += (opt.species.energy + opt.species.zpe) * constants.AUtoKCAL
                        right_zeroenergy = ts_zeroenergy - prod_zeroenergy
                        reaction.mess_left_endpoint_source = left_reference_job

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
                    if self.species.reac_type[index] != 'hom_sci':
                        ts_blocks[reaction.instance_name] = allTS[reaction.instance_name]
                    if self.par['pes'] and reaction.do_vdW:
                        st_pt = reaction.irc_prod_opt.species
                        energy_add = uq.calc_factor('energy', uq_iter)
                        freq_factor = uq.calc_factor('freq', uq_iter)
                        well_blocks[st_pt.name] = self.write_well(st_pt,
                                                                    energy_add,
                                                                    freq_factor,
                                                                    uq_iter)
                        bimol_name = '_'.join(sorted([routing_name(st_pt) for st_pt in reaction.products]))
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
                        well_blocks[routing_key(st_pt)] = self.write_well(st_pt,
                                                                    energy_add,
                                                                    freq_factor,
                                                                    uq_iter)
                    elif len(reaction.products) == 2:
                        bimol_name = '_'.join(sorted([routing_name(st_pt) for st_pt in reaction.products]))
                        if bimol_name in written_bimol_names:
                            continue
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
                        termol_name = '_'.join(sorted([routing_name(st_pt) for st_pt in reaction.products]))
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
            for block in ts_blocks.values():
                tss += block + divider
            barrierless = ''
            for rxn in barrierless_blocks:
                barrierless += barrierless_blocks[rxn] + divider

            if 'prod' in self.species.name:
                dum = self.dummytpl.format(barrier='tsd', reactant=self.well_names[self.species.name], dummy='d1')
            else:
                dum = self.dummytpl.format(barrier='tsd', reactant=self.well_names[routing_key(self.species)], dummy='d1')

            mess_iter = "{0:04d}".format(uq_iter)

            contents = header + divider + wells + bimols + tss + termols + barrierless + divider + 'End ! end kinetics\n'
            if not self.par['pes']:
                contents = complete_mirror_channels(contents)
                validate_mess_populations(contents)
                if self.par['multi_conf_tst']:
                    contents = finalize_mc_mess(contents, self.par.get('correct_submerged', 0))
            # PES worker files still contain deferred energy/name placeholders.
            # Fold and validate their populations in final PES assembly instead.
            if self.par['pes']:
                with open('me/mess_%s.inp' % mess_iter, 'w') as f_out:
                    f_out.write(contents)
            else:
                from kinbot.mess_networks import write_network_inputs
                write_network_inputs(self, contents, uq_iter)

        return 0


    def write_termol(self, species_list, reaction, uq_iter, bless=0):
        # Create the dummy MESS block for ter-molecular products.
        if self.par.get('optical_population', 'specified') == 'racemic':
            self._mc_product_identities(species_list)
        termol = ''
        terPr_name = '_'.join(sorted([routing_name(species) for species in species_list]))
        prod_name = '{name}' if self.par['pes'] else self.termolec_names[terPr_name]
        termol += self.termoltpl.format(name=prod_name, product=terPr_name)
        termol = annotate_population(termol, species_list, self.par.get('optical_population', 'specified'))
        mess_iter = "{0:04d}".format(uq_iter)
        with open(mess_filename(terPr_name, uq_iter), 'w') as f:
            f.write(termol)

        return termol


    def _check_optical_models(self, ts_all, ts_unique):
        """Reject inconsistent data; geometric uncertainty uses a warned weight one."""
        states = [self.species]
        representative_only = set()
        for index, reaction in enumerate(self.species.reac_obj):
            if (reaction.instance_name in ts_all
                    and (self.par['pes'] or reaction.instance_name in ts_unique)
                    and self.species.reac_type[index] != 'hom_sci'):
                self._set_barrier_population(reaction)
                states.append(reaction.ts)
                if self.species.reac_type[index] == 'barrierless_saddle':
                    representative_only.add(id(reaction.ts))
            if reaction.instance_name in ts_unique:
                states.extend(opt.species for opt in reaction.prod_opt)
                if self.par['pes'] and reaction.do_vdW:
                    states.append(reaction.irc_prod_opt.species)
        unresolved, checked = [], set()
        for species in states:
            if id(species) in checked or species.natom == 1:
                continue
            checked.add(id(species))
            try:
                if self.par['multi_conf_tst'] and id(species) not in representative_only:
                    members = writer_members(species, self.par.get('optical_population', 'specified'),
                                             preserve_errors=False)
                    if not members:
                        representative_record(species, self.par.get('optical_population', 'specified'),
                                              preserve_errors=False)
                else:
                    self._parent_symmetry(species, preserve_errors=False,
                                          single_structure=id(species) in representative_only)
            except ValueError as error:
                unresolved.append((species, str(error)))
        if unresolved:
            refuse_routing('Unresolved optical decisions:\n' + '\n'.join(
                f'{species.name}: {reason}' for species, reason in unresolved),
                [species for species, _ in unresolved])

    def _set_barrier_population(self, reaction):
        from kinbot.reaction_path import set_endpoint_populations
        products = [opt.species for opt in getattr(reaction, 'prod_opt', [])] or reaction.products
        set_endpoint_populations(reaction.ts, [self.species], products)
        if not hasattr(reaction.ts, 'optical_reference'):
            reaction.ts.optical_reference = canonical_identity(self.species)
        return products

    def _parent_symmetry(self, species, *, preserve_errors=True, single_structure=False):
        species.optical_population = self.par.get('optical_population', 'specified')
        if self.par['multi_conf_tst'] and not single_structure:
            record = representative_record(species, self.par.get('optical_population', 'specified'))
            species.mess_optical_counting = record.optical_evidence
            return record.sigma_ext / record.remaining_optical_weight
        from copy import copy
        from kinbot.counting_contract import optical_counting
        from kinbot.thermochemistry import hir_evidence
        from kinbot.optical import bind_assumption
        bind_assumption(species, self.par)
        view = copy(species)
        view.conformer_representation = 'single structure'
        counting = optical_counting(view, hir_evidence(view))
        if (counting['status'] == 'unresolved'
                and counting.get('remaining_multiplier') is None
                and getattr(species, 'hir', None) is not None):
            # write_input already attempted recovery with QC access. Standalone
            # block writers must also omit an unusable HIR model consistently.
            from kinbot.hindered_rotors import use_harmonic_model
            use_harmonic_model(species, self.par, 'Unusable HIR model: ' + counting['reason'])
            view = copy(species)
            view.conformer_representation = 'single structure'
            counting = optical_counting(view, hir_evidence(view))
        species.mess_optical_counting = counting
        if (counting['status'] not in ('resolved', 'assumed', 'legacy_unverified')
                and counting.get('fallback') != 'unresolved_symmetry'):
            reason = 'Unresolved optical coverage: ' + counting['reason']
            if preserve_errors:
                refuse_routing(reason, [species])
            raise ValueError(reason)
        return float(species.sigma_ext) / counting['remaining_multiplier']


    def _member_rrho(self, species, record, freq_factor, zeroenergy, *,
                     saddle=False, tunneling='', shift=''):
        """Render one MC member from its associated geometry and properties."""
        modes = frequencies.thermochemical_frequencies(
            record.frequencies_cm1, saddle, self.par.get('imagfreq_threshold', 50.))
        return self._counting_comment(record.optical_evidence, record.member_id) + self.rrhotpl.format(
            natom=species.natom, geom=self.rotor_geom(species, record.geometry, record.frequencies_cm1),
            core=self.corerrtpl.format(symm=record.sigma_ext / record.remaining_optical_weight,
                rotconst=self.rotor_core_line(species, record.geometry, record.frequencies_cm1)),
            nfreq=len(modes) - int(saddle), freq=self.make_freq(modes, freq_factor, int(saddle)),
            rotors='', tunneling=tunneling, nelec=1, mult=species.mult,
            zeroenergy=zeroenergy, shift=shift)

    def _population_comment(self, species):
        comment = self._optical_comment(species)
        if self.par.get('optical_population', 'specified') == 'racemic':
            identity = optical_scope(species, 'racemic')['identity']
            if identity.get('status') == 'assigned' and identity['is_chiral_configuration']:
                return comment + (f"! kinbot_racemic_population {routing_name(species)} "
                        f"{identity['mirror_family_id']}\n")
        return comment

    @staticmethod
    def _optical_comment(species):
        count = getattr(species, 'mess_optical_counting', {}) or {}
        return MESS._counting_comment(count, species.name)

    @staticmethod
    def _counting_comment(count, name):
        count = count or {}
        notes = ''.join('! WARNING: ' + str(note) + '\n' for note in count.get('warnings', ()))
        warning = MESS._optical_warning(count, name)
        if warning:
            return notes + warning
        pairs = count.get('explicit_mirror_comparisons', {})
        comments = []
        for partner in count.get('explicit_mirror_ids', ()):
            pair = pairs.get(partner, {})
            if pair.get('status') == 'match':
                energies = [pair[key]['stable_midpoint_energy_kcal_mol'] for key in ('forward', 'reverse')]
                comments.append(f'! Optical factor 1: explicit mirror {partner}; harmonic midpoint '
                    f'approximation {energies[0]:.6g}/{energies[1]:.6g} kcal/mol, '
                    f'cutoff {pair["cutoff_kcal_mol"]:g} kcal/mol.\n')
        if comments:
            return notes + ''.join(comments)
        if count.get('heuristic') == 'harmonic_midpoint':
            message = count['reason']
            logger.warning('%s: %s', name, message)
            return notes + '! ' + message + '\n'
        if count.get('status') == 'assumed':
            return notes + (f"! optical factor {count['remaining_multiplier']}: explicit assumption; "
                    f"mirror coverage undetermined. {count['reason']}\n")
        if count.get('status') == 'legacy_unverified':
            return notes + '! ' + count['reason'] + '\n'
        return notes

    @staticmethod
    def _optical_warning(count, name):
        if not count or count.get('fallback') != 'unresolved_symmetry':
            return ''
        message = (f'WARNING: unresolved symmetry number for {name}; using optical factor 1 '
                   f'(external rotational symmetry unchanged). {count["reason"]}')
        midpoint = count.get('harmonic_midpoint', {})
        if midpoint.get('status') in ('complete', 'limited'):
            message += (f" Harmonic midpoint: {midpoint['stable_midpoint_energy_kcal_mol']:.4g} kcal/mol "
                        f"in stable modes; negative-mode share of squared mass-weighted displacement "
                        f"{midpoint['negative_mode_displacement_fraction']:.1%}. "
                        "Diagnostic only, not an inversion barrier; no energy cutoff applied.")
            if midpoint['status'] == 'limited':
                message += ' Atom-mapping search incomplete.'
        logger.warning(message)
        return '! ' + message + '\n'

    def _mc_product_identities(self, products):
        identities = [canonical_identity(product) for product in products]
        if any(identity['status'] != 'assigned' for identity in identities):
            from kinbot.stereo_identity import legacy_stereo_warning
            for product, identity in zip(products, identities):
                if identity['status'] != 'assigned':
                    legacy_stereo_warning(product, identity.get('reason'))
            return identities
        if (self.par.get('optical_population', 'specified') == 'racemic'
                and sum(identity['is_chiral_configuration'] for identity in identities) > 1):
            refuse_routing('A global racemic pair is not independent racemates of multiple fragments', products)
        return identities

    def _validate_mc_endpoints(self, reaction, members):
        products = ([opt.species for opt in getattr(reaction, 'prod_opt', [])]
                    or reaction.products)
        identities = self._mc_product_identities(products)
        if any(identity['status'] != 'assigned' for identity in identities):
            return  # The products explicitly use the warned legacy treatment.
        scope = optical_scope(reaction.ts, self.par.get('optical_population', 'specified'))
        if (not getattr(reaction.ts, 'ts_endpoint_identities', ())
                and scope['population'] == 'specified' and scope['mirror_allowed']
                and any(record.mirror_states == 2 for record in members)
                and Counter(identity['id'] for identity in identities)
                    != Counter(identity['mirror_id'] for identity in identities)):
            refuse_routing('The full TS mirror orbit is incompatible with the specified product population',
                           [self.species, reaction.ts, *products])

    def write_bimol(self, prod_list, well_add, freq_factor, pstsymm_factor, uq_iter, bless, vdW=False):
        """
        Create the block for MESS for a bimolecular product.
        In case of a barrierless reaction (bless=1) also add a phase-space theory barrier.
        well0: reactant on this PES (zeroenergy reference)
        uq_n = number of uncertainty runs
        """

        if (self.par['multi_conf_tst']
                or self.par.get('optical_population', 'specified') == 'racemic'):
            self._mc_product_identities(prod_list)

        fragments = ''
        fragment_reference_shift = 0.
        smi = []
        for nsp, species in enumerate(prod_list):
            smi.append(species.smiles)
            if species.natom > 1:

                if self.par['pes']:
                    name = '{{fr_name_{}}}'.format(routing_key(species))
                else:
                    name = self.fragment_names[routing_key(species)] + ' ! ' + routing_name(species)
                members = (writer_members(species, self.par.get('optical_population', 'specified'))
                           if self.par['multi_conf_tst'] and species.chemid != 170170000000000000002
                           else {})
                # molecule template
                if species.chemid == 170170000000000000002:  # exception for OH
                    fragment = self.fragmenttplOH.format(chemid=name,
                                                           smi=species.smiles,
                                                           natom=species.natom,
                                                           geom=self.rotor_geom(species),
                                                           symm=float(species.sigma_ext) / float(species.nopt),
                                                           rotconst=self.rotor_core_line(species),
                                                           freq=self.make_freq(species.reduced_freqs, freq_factor, 0))
                elif members:
                    base = min(record.zero_energy_hartree for record in members.values())
                    fragment_reference_shift += (base - species.energy - species.zpe) * constants.AUtoKCAL
                    fragment = f'  Fragment {name} ! {species.smiles}\n    Union\n'
                    for record in members.values():
                        fragment += self._member_rrho(species, record, freq_factor,
                            round((record.zero_energy_hartree - base) * constants.AUtoKCAL, 2))
                    fragment += '    End ! Union\n'
                else:
                    fragment = self.fragmenttpl.format(chemid=name,
                                                         smi=species.smiles,
                                                         natom=species.natom,
                                                         geom=self.rotor_geom(species),
                                                         symm=self._parent_symmetry(species),
                                                         rotconst=self.rotor_core_line(species),
                                                         nfreq=len(species.reduced_freqs),
                                                         freq=self.make_freq(species.reduced_freqs, freq_factor, 0),
                                                         hinderedrotor=self.make_rotors(species, freq_factor),
                                                         nelec=1,
                                                         mult=species.mult)
            else:
                if self.par['pes']:
                    name = '{{fr_name_{}}}'.format(routing_key(species))
                else:
                    name = self.fragment_names[routing_key(species)] + ' ! ' + routing_name(species)

                fragment = self.atomtpl.format(chemid=name,
                                                 element=species.atom[0],
                                                 nelec=1,
                                                 mult=species.mult)

            # Compute comments after rendering (which sets optical warnings),
            # but put the identity marker before the fragment it describes.
            fragments += self._population_comment(species) + fragment


        pr_name = '_'.join(sorted([routing_name(species) for species in prod_list]))
        if self.par['pes']:
            name = '{{name}} ! {} {}'.format(smi[0], smi[1])
            energy_reference = '{ground_energy}'
            energy = '{ground_energy}' + (f' ! {fragment_reference_shift}' if fragment_reference_shift else '')
        else:
            name = '{} ! {}'.format(self.bimolec_names[pr_name], pr_name)
            energy = (sum([sp.energy for sp in prod_list]) + sum([sp.zpe for sp in prod_list]) 
                      - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL
            energy_reference = round(energy + well_add, 2)
            energy += well_add + fragment_reference_shift
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
            well_name = self.well_names[routing_key(self.species)] if not self.par['pes'] else None
            bimol = self.blbimoltpl.format(barrier=('{blessname}' if self.par['pes'] else
                                                   f'bl_{well_name}_{self.bimolec_names[pr_name]}'),
                                           reactant=('{wellname}' if self.par['pes'] else well_name),
                                           prod=('{prodname}' if self.par['pes'] else self.bimolec_names[pr_name]),
                                           chemids=name,
                                           model=self._phase_space_models(prod_list, stoich, energy,
                                               fragment_reference_shift, freq_factor, pstsymm_factor),
                                           fragments=fragments,
                                           ground_energy=energy)

        bimol = annotate_population(bimol, prod_list, self.par.get('optical_population', 'specified'),
            energy_reference=energy_reference if self.par['multi_conf_tst'] else None)
        with open(mess_filename(pr_name, uq_iter), 'w') as f:
            f.write(bimol)

        return bimol

    def _phase_space_models(self, products, stoich, ground_energy, reference_shift,
                            freq_factor, symmetry_factor):
        """Sum the same product conformers used by the separated fragments.

        The capture potential, total electronic spin and parent-product
        symmetry normalization remain KinBot's existing phase-space model.
        MC terms carry relative rotational/optical divisors with respect to
        those parent products; this is not an absolute PST symmetry correction.
        No saddle geometry is involved.
        """
        choices = []
        for species in products:
            # OH uses the existing single-geometry spin-orbit fragment model;
            # atoms do not have a vibrational conformer ensemble.
            members = (writer_members(species, self.par.get('optical_population', 'specified'))
                       if self.par['multi_conf_tst'] and species.natom > 1
                       and species.chemid != 170170000000000000002 else {})
            if self.par['multi_conf_tst'] and species.natom > 1 and not members \
                    and species.chemid != 170170000000000000002:
                choices.append([representative_record(species, self.par.get('optical_population', 'specified'))])
            else:
                choices.append(list(members.values()) or [None])
        models = []
        for combination in product(*choices):
            geometries, modes, rotors = [], [], []
            offset = 0.
            divisor = symmetry_factor
            for species, record in zip(products, combination):
                if record is None:
                    geom, freq = species.geom, species.reduced_freqs
                    rotors.append(self.make_rotors(species, freq_factor, bless=True))
                else:
                    geom = record.geometry
                    freq = frequencies.thermochemical_frequencies(
                        record.frequencies_cm1, False, self.par.get('imagfreq_threshold', 50.))
                    sigma = record.sigma_ext / record.remaining_optical_weight
                    offset += (record.zero_energy_hartree - species.energy - species.zpe) * constants.AUtoKCAL
                if record is not None:
                    divisor *= sigma / self._parent_symmetry(species)
                modes.extend(freq)
                geometries.append(self.pstfragmenttpl.format(chemid=routing_name(species),
                    smi=species.smiles, natom=species.natom,
                    geom=self.rotor_geom(species, geom, freq)))
            zero = ('{ground_energy}' + (f' ! {offset}' if offset else '') if self.par['pes']
                    else round(float(ground_energy) - reference_shift + offset, 2))
            models.append(self.pstrrhotpl.format(stoich=stoich,
                frag1=geometries[0], frag2=geometries[1], pstsymm=divisor,
                nfreq=len(modes), freq=self.make_freq(modes, freq_factor, 0),
                hinderedrotor='\n'.join(rotors), mult=self.species.mult, zeroenergy=zero))
        if len(models) == 1:
            model = models[0]
        else:
            model = ('    Union ! product conformer combinations\n' + ''.join(models)
                     + '    End ! product conformer combinations\n')
        if self.par['multi_conf_tst']:
            model = ('! MC phase-space weights are relative to the selected product structures.\n'
                     '! The existing absolute capture normalization is retained.\n' + model)
        return model


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
                name = self.well_names[routing_key(species)] + ' ! ' + routing_name(species)
                norot = None
            else:
                name = self.well_names[species.name] + ' ! ' + str(species.name)
                norot = str(species.name)
            zeroenergy = ((species.energy + species.zpe) - (self.species.energy + self.species.zpe)) * constants.AUtoKCAL
            zeroenergy += well_add
            zeroenergy = round(zeroenergy, 2)

        member_records = (writer_members(species, self.par.get('optical_population', 'specified'))
                          if self.par['multi_conf_tst'] else {})
        valid_conformers = list(member_records)
        nunq_confs = len(valid_conformers)

        if not self.par['multi_conf_tst'] or not valid_conformers:
            mess_well = self.welltpl.format(chemid=name,
                                            smi=species.smiles,
                                            natom=species.natom,
                                            geom=self.rotor_geom(species),
                                            symm=self._parent_symmetry(species),
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
            for record in member_records.values():
                shift = constants.AUtoKCAL * (record.zero_energy_hartree - base_zeroen)
                conformer_zeroenergy = (zeroenergy if self.par['pes'] else
                                        round(zeroenergy + shift, 2))
                rrho += self._member_rrho(species, record, freq_factor,
                    conformer_zeroenergy, shift=shift if self.par['pes'] else '')
            rrho = '      '.join(rrho.splitlines(True))  # indent
            mess_well = self.welluniontpl.format(chemid=name,
                                                 smi=species.smiles,
                                                 nunion=nunq_confs,
                                                 rrho=rrho)
        mess_well = annotate_population(mess_well, [species], self.par.get('optical_population', 'specified'),
            energy_reference=zeroenergy if self.par['multi_conf_tst'] else None)
        mess_well = self._population_comment(species) + mess_well
        if 'prod' not in species.name:
            with open('{}_{:04d}.mess'.format(routing_key(species), uq_iter), 'w') as f:
                f.write(mess_well)
        else:
            with open('{}_{:04d}.mess'.format(species.name, uq_iter), 'w') as f:
                f.write(mess_well)

        return mess_well

    def write_barrier(self, reaction, index, left_zeroenergy, right_zeroenergy, barrier_add, freq_factor, imagfreq_factor, uq_iter):
        """Create the block for a MESS barrier."""
        variational = self.species.reac_type[index] == 'barrierless_saddle'
        use_ensemble = self.par['multi_conf_tst'] and not variational
        variational_warning = ''
        if self.par['multi_conf_tst'] and variational:
            message = (f'{reaction.instance_name}: variational barrier uses the selected '
                       'representative TS and product structures, as in the legacy writer; '
                       'additional conformers are not included in this barrier model.')
            logger.warning(message)
            variational_warning = '! WARNING: ' + message + '\n'

        left_zeroenergy += barrier_add
        right_zeroenergy += barrier_add

        products = self._set_barrier_population(reaction)
        member_records = (writer_members(reaction.ts, self.par.get('optical_population', 'specified'))
                          if use_ensemble else {})
        valid_conformers = list(member_records)
        nunq_confs = len(valid_conformers)

        if use_ensemble:
            counted = list(member_records.values()) or [representative_record(
                reaction.ts, self.par.get('optical_population', 'specified'))]
            self._validate_mc_endpoints(reaction, counted)
            tunneling_products = ([reaction.irc_prod_opt.species] if getattr(reaction, 'do_vdW', False)
                        else [opt.species for opt in getattr(reaction, 'prod_opt', [])]
                        or reaction.products)
            def parent_observation(point):
                return {'source_job': getattr(point, 'source_job', None),
                        'electronic_energy_hartree': float(point.energy),
                        'zpe_hartree': float(point.zpe)}
            reaction.ts.mess_tunneling_reference = {
                'convention': 'selected-parent endpoint approximation',
                'member_connectivity': 'not individually revalidated',
                'left_parent_observation': parent_observation(self.species),
                'right_parent_observations': [parent_observation(p) for p in tunneling_products],
                'left_endpoint_source_override': getattr(reaction, 'mess_left_endpoint_source', None),
                'input_depths_kcal_mol': (None if self.par['pes'] else
                                         [float(left_zeroenergy), float(right_zeroenergy)]),
                'energy_resolution': ('final depths resolved by PES at its selected energy level'
                                      if self.par['pes'] else 'caller-selected endpoint depths'),
                'ensemble_minima_redefine_depths': False,
            }

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
        if self.par['pes'] and getattr(reaction, 'do_vdW', False):
            prod_name = self.well_names[reaction.irc_prod_opt.species.name]
        elif len(reaction.products) == 1:
            prod_name = self.well_names[routing_key(reaction.products[0])]
        elif len(reaction.products) == 2:
            long_name = '_'.join(sorted([routing_name(pi) for pi in reaction.products]))
            prod_name = self.bimolec_names[long_name]
        else:
            long_name = '_'.join(sorted([routing_name(pi) for pi in reaction.products]))
            prod_name = self.termolec_names[long_name]

        if self.par['pes']:
            name = '{name}'
            chemid_reac = ''
            chemid_prod = ''
            long_rxn_name = ''
            zeroenergy = '{zeroenergy}'
        else:
            name = self.ts_names[reaction.instance_name]
            chemid_reac = self.well_names[routing_key(self.species)]
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
                prodzeroenergy = ((reaction.prod_opt[0].species.energy + reaction.prod_opt[0].species.zpe + \
                                 reaction.prod_opt[1].species.energy + reaction.prod_opt[1].species.zpe) - \
                                 (self.species.energy + self.species.zpe)) * constants.AUtoKCAL

            outerts = self.psttpl.format(natom1=reaction.prod_opt[0].species.natom,
                                         geom1=self.rotor_geom(reaction.prod_opt[0].species),
                                         natom2=reaction.prod_opt[1].species.natom,
                                         geom2=self.rotor_geom(reaction.prod_opt[1].species),
                                         symm=self._parent_symmetry(reaction.ts, single_structure=True),
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
            corerr = self.corerrtpl.format(symm=self._parent_symmetry(reaction.ts, single_structure=True),
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
            corerr = self.corerrtpl.format(symm=self._parent_symmetry(reaction.ts),
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
            for record in member_records.values():
                shift = constants.AUtoKCAL * (record.zero_energy_hartree - base_zeroen)
                conformer_zeroenergy = (zeroenergy if self.par['pes'] else
                                        round(zeroenergy + shift, 2))
                imfreq = round(-record.frequencies_cm1[0] * imagfreq_factor, 2)
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
                rrho += self._member_rrho(reaction.ts, record, freq_factor,
                    conformer_zeroenergy, saddle=True, tunneling=tun_conf,
                    shift=shift if self.par['pes'] else '')
            rrho = '      '.join(rrho.splitlines(True))  # indent
            mess_barrier = self.barrieruniontpl.format(rxn_name=name,
                                                       chemid_reac=chemid_reac,
                                                       chemid_prod=chemid_prod,
                                                       long_rxn_name=long_rxn_name,
                                                       nunion=nunq_confs,
                                                       model=rrho)


        if self.par['multi_conf_tst']:
            mess_barrier = ('! MC Eckart convention: selected-parent endpoint approximation.\n'
                            '! MC ensemble minima do not redefine these WellDepth values.\n'
                            + mess_barrier)
        if not self.par['pes'] and getattr(reaction, 'do_vdW', False):
            mess_barrier = ('! Direct model: separated products; optional IRC complex and its fast exit omitted.\n'
                            '! Inner-TS tunneling depths retain the original IRC complex reference.\n'
                            + mess_barrier)
        mess_barrier = variational_warning + self._optical_comment(reaction.ts) + mess_barrier
        connected_products = ([reaction.irc_prod_opt.species]
                              if self.par['pes'] and getattr(reaction, 'do_vdW', False)
                              else products)
        mess_barrier = annotate_endpoints(mess_barrier, [self.species], connected_products,
                                          self.par.get('optical_population', 'specified'))
        path_id = reaction_path_id(reaction)
        if path_id is not None:
            mess_barrier = f'! kinbot_stereopath {path_id}\n' + mess_barrier
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
