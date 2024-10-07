import numpy as np
import sys
import copy
import logging

from kinbot import bond_combinations
from kinbot import find_motif
from kinbot.reactions.reac_Cyclic_Ether_Formation import CyclicEtherFormation
from kinbot.reactions.reac_Diels_alder_addition import DielsAlder
from kinbot.reactions.reac_Intra_Diels_alder_R import IntraDielsAlder
from kinbot.reactions.reac_12_shift_S_F import S12ShiftF
from kinbot.reactions.reac_12_shift_S_R import S12ShiftR
from kinbot.reactions.reac_cpd_H_migration import CpdHMigration
from kinbot.reactions.reac_intra_H_migration import IntraHMigration
from kinbot.reactions.reac_intra_H_migration_suprafacial import IntraHMigrationSuprafacial
from kinbot.reactions.reac_intra_OH_migration import IntraOHMigration
from kinbot.reactions.reac_intra_OH_migration_Exocyclic_F import IntraOHMigrationExocyclicF
from kinbot.reactions.reac_Intra_R_Add_Endocyclic_F import IntraRAddEndocyclicF
from kinbot.reactions.reac_Intra_R_Add_Exocyclic_F import IntraRAddExocyclicF
from kinbot.reactions.reac_Intra_R_Add_ExoTetCyclic_F import IntraRAddExoTetCyclicF
from kinbot.reactions.reac_intra_R_migration import IntraRMigration
from kinbot.reactions.reac_Retro_Ene import RetroEne
from kinbot.reactions.reac_r22_cycloaddition import R22Cycloaddition
from kinbot.reactions.reac_r12_insertion_R import R12Insertion
from kinbot.reactions.reac_r13_insertion_RSR import R13InsertionRSR
from kinbot.reactions.reac_r13_insertion_ROR import R13InsertionROR
from kinbot.reactions.reac_r13_insertion_CO2 import R13InsertionCO2
from kinbot.reactions.reac_r12_cycloaddition import R12Cycloaddition
from kinbot.reactions.reac_R_Addition_MultipleBond import RAdditionMultipleBond
from kinbot.reactions.reac_R_Addition_CSm_R import RAdditionCS
from kinbot.reactions.reac_R_Addition_COm3_R import RAdditionCO
from kinbot.reactions.reac_Korcek_step2_odd import KorcekStep2Odd
from kinbot.reactions.reac_Korcek_step2_even import KorcekStep2Even
from kinbot.reactions.reac_Korcek_step2 import KorcekStep2
from kinbot.reactions.reac_ketoenol import KetoEnol
from kinbot.reactions.reac_Intra_RH_Add_Exocyclic_R import IntraRHAddExoR
from kinbot.reactions.reac_Intra_RH_Add_Exocyclic_F import IntraRHAddExoF
from kinbot.reactions.reac_Intra_RH_Add_Endocyclic_R import IntraRHAddEndoR
from kinbot.reactions.reac_Intra_RH_Add_Endocyclic_F import IntraRHAddEndoF
from kinbot.reactions.reac_HO2_Elimination_from_PeroxyRadical import HO2Elimination
from kinbot.reactions.reac_beta_delta import BetaDelta
from kinbot.reactions.reac_birad_recombination_F import BiradRecombinationF
from kinbot.reactions.reac_birad_recombination_R import BiradRecombinationR
from kinbot.reactions.reac_Intra_disproportionation_R import IntraDisproportionationR
from kinbot.reactions.reac_Intra_disproportionation_F import IntraDisproportionationF
from kinbot.reactions.reac_r14_birad_scission import R14BiradScission
from kinbot.reactions.reac_r14_cyclic_birad_scission_R import R14CyclicBiradScission
from kinbot.reactions.reac_barrierless_saddle import BarrierlessSaddle
from kinbot.reactions.reac_h2_elim import H2Elim
from kinbot.reactions.reac_bimol_disproportionation_R import BimolDisproportionationR
from kinbot.reactions.reac_homolytic_scission import HS
from kinbot.reactions.reac_combinatorial import Combinatorial

logger = logging.getLogger('KinBot')


class ReactionFinder:
    '''
    Class to find all the potential reactions starting from a well.
    '''
    
    def __init__(self, species, par, qc):
        self.species = species
        self.qc = qc
        self.par = par
        self.families = par['families']
        self.skip_families = par['skip_families']
        self.specific_reaction = par['specific_reaction']
        self.break_bond = par['break_bonds']
        self.form_bond = par['form_bonds']
        self.ringrange = range(self.par['ringrange'][0], self.par['ringrange'][1])
        self.one_reaction_comb = par['one_reaction_comb']
        self.one_reaction_fam = par['one_reaction_fam']
        # make a set of frozen sets from the breaking and forming bond lists
        self.reac_bonds = set()
        for i, bond in enumerate(par['break_bonds']):
            self.reac_bonds.add(frozenset(par['break_bonds'][i]))
        self.prod_bonds = set()
        for i, bond in enumerate(par['form_bonds']):
            self.prod_bonds.add(frozenset(par['form_bonds'][i]))
        try:
            self.barrierless_saddle = par['barrierless_saddle'][str(self.species.chemid)]
        except KeyError:
            self.barrierless_saddle = None
            
        # keys: names of the families
        # values: list of instances
        # this dict is used to keep track of the unique reactions found,
        # and to verify whether a new reaction is indeed unique 
        self.reactions = {}

    def find_reactions(self):
        '''
        List all reaction types available, and find the key atoms for them.
        '''

        reaction_names = {'intra_H_migration': self.search_intra_H_migration,
                          'intra_H_migration_suprafacial': self.search_intra_H_migration_suprafacial,
                          'intra_R_migration': self.search_intra_R_migration,
                          'intra_OH_migration': self.search_intra_OH_migration,
                          'intra_OH_migration_Exocyclic_F': self.search_intra_OH_migration_Exocyclic_F,
                          'cpd_H_migration': self.search_cpd_H_migration, 
                          'Intra_RH_Add_Endocyclic_F': self.search_Intra_RH_Add_Endocyclic_F,
                          'Intra_RH_Add_Endocyclic_R': self.search_Intra_RH_Add_Endocyclic_R,
                          'Cyclic_Ether_Formation': self.search_Cyclic_Ether_Formation,
                          'Intra_RH_Add_Exocyclic_F': self.search_Intra_RH_Add_Exocyclic_F,
                          'Intra_RH_Add_Exocyclic_R': self.search_Intra_RH_Add_Exocyclic_R,
                          'Retro_Ene': self.search_Retro_Ene,
                          'Intra_R_Add_Endocyclic_F': self.search_Intra_R_Add_Endocyclic_F, 
                          'Intra_R_Add_ExoTetCyclic_F': self.search_Intra_R_Add_ExoTetCyclic_F,
                          'Intra_R_Add_Exocyclic_F': self.search_Intra_R_Add_Exocyclic_F, 
                          'Korcek_step2_odd': self.search_Korcek_step2_odd, 
                          'Korcek_step2_even': self.search_Korcek_step2_even, 
                          'Korcek_step2': self.search_Korcek_step2, 
                          'r22_cycloaddition': self.search_r22_cycloaddition, 
                          'r12_cycloaddition': self.search_r12_cycloaddition, 
                          'r12_insertion_R': self.search_r12_insertion_R, 
                          'r13_insertion_CO2': self.search_r13_insertion_CO2, 
                          'r13_insertion_ROR': self.search_r13_insertion_ROR, 
                          'Diels_alder_addition': self.search_Diels_alder_addition, 
                          'Intra_Diels_alder_R': self.search_Intra_Diels_alder_R, 
                          'ketoenol': self.search_ketoenol, 
                          'HO2_Elimination_from_PeroxyRadical': self.search_HO2_Elimination_from_PeroxyRadical, 
                          'R_Addition_COm3_R': self.search_R_Addition_COm3_R, 
                          'R_Addition_MultipleBond': self.search_R_Addition_MultipleBond, 
                          '12_shift_S_F': self.search_12_shift_S_F, 
                          '12_shift_S_R': self.search_12_shift_S_R, 
                          'R_Addition_CSm_R': self.search_R_Addition_CSm_R, 
                          'r13_insertion_RSR': self.search_r13_insertion_RSR, 
                          'beta_delta': self.search_beta_delta, 
                          'h2_elim': self.search_h2_elim,
                          'hom_sci': self.search_hom_sci,
                          'barrierless_saddle': self.search_barrierless_saddle,
                          'Intra_disproportionation_F': self.search_Intra_disproportionation_F, 
                          'Intra_disproportionation_R': self.search_Intra_disproportionation_R, 
                          'bimol_disproportionation_R': self.search_bimol_disproportionation_R, 
                          }

        if 'combinatorial' in self.families:
            reaction_names['combinatorial'] = self.search_combinatorial        
 
        atom = self.species.atom
        natom = self.species.natom
        
        for i, bond in enumerate(self.species.bonds):  # for all resonance structures
            rad = self.species.rads[i]
            
            if self.one_reaction_comb:
                # search for just one reaction, given by the list of bonds to be 
                # broken or formed
                
                # based on the combinatorial reaction family, because they are also
                # defined by the list of bonds to be broken or formed
                name = 'combinatorial'
                self.reactions[name] = []
                
                self.reac_bonds = self.par['break_bonds']
                self.prod_bonds = self.par['form_bonds']
                ts = bond_combinations.generate_ts(self.reac_bonds, self.prod_bonds, self.species.bond)
                self.reactions[name].append([self.reac_bonds, self.prod_bonds, ts, 1])
                
            else:
                for rn in reaction_names:
                    if 'all' in self.families or rn in self.families:
                        if not rn in self.skip_families:
                            reaction_names[rn](natom, atom, bond, rad)

        for name in self.reactions:
            self.reaction_matrix(self.reactions[name], name) 
        
        for index in range(len(self.species.reac_name)-1):
            if self.species.reac_name[index] in self.species.reac_name[index + 1:]:
                logger.error(f'Found reaction name {self.species.reac_name[index]} more than once')
                logger.error('Exiting')
                sys.exit()

        logger.info('\tFound the following reactions:')
        for rxn in self.species.reac_name:
            logger.info('\t\t{}'.format(rxn))
        
        return 0  
   

    def search_combinatorial(self, natom, atom, bond, rad):
        ''' 
        This is a method to create all possible combinations of maximum 3 bond breakings 
        and maximum 3 bond formations.
        TODO: allow bond breaking without the atoms forming new bond (only for radicals)
        '''
        
        name = 'combinatorial'

        if not name in self.reactions:
            self.reactions[name] = []

        instances = bond_combinations.generate_all_product_bond_matrices(self.species, self.par)
        for inst in instances:
            self.reactions[name].append(inst)
        #~ self.reactions[name] = []
        #~ reac = [[0, 5], [1, 2], [3, 4]]
        #~ prod = [[0, 1], [2, 3], [4, 5]]
        #~ ts = self.species.bond
        #~ self.reactions[name].append([reac, prod, ts])
        return 0


    def search_intra_H_migration(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        H-R~~~~~~~R* <==> R*~~~~~~~R-H

        Works in both directions.
        H is moved to
        * radical site
        * multiple bond
        * lone pair
        '''
        
        name = 'intra_H_migration'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        # if np.sum(rad) == 0:
        # find H-migrations over double bonds and to lone pairs

        for ringsize in self.ringrange:
            motif = ['X' for i in range(ringsize)]
            motif[-1] = 'H'
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)

            # double bonds
            for instance in instances:
                if any([bi > 1 for bi in bond[instance[0]]]):
                    rxns += [instance]

            # lone pairs
            for instance in instances:
                if (self.species.atom[instance[0]] == 'O' or
                        self.species.atom[instance[0]] == 'S' or
                        self.species.atom[instance[0]] == 'N'):
                    rxns += [instance]

                # carbene
                if self.species.atom[instance[0]] == 'C' and rad[instance[0]] == 2:
                    rxns += [instance]

            # cations (try all)
            if self.species.charge == 1:
                for instance in instances:
                    if self.species.atom[instance[0]] != 'H':
                        rxns += [instance]

        if np.sum(rad) != 0:
            #        else:
            instances = []
            for ringsize in self.ringrange:
                motif = ['X' for i in range(ringsize)]
                motif[-1] = 'H'
                for rad_site in np.nonzero(rad)[0]:
                    instances += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)
            for instance in instances:
                rxns.append(instance)
        rxns = self.clean_rigid(name, rxns, 0, -1)

        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
#                    new = 0
        
        return 0

    def search_intra_H_migration_suprafacial(self, natom, atom, bond, rad):
        ''' 
        This is a special case of H migration reactions over a double bond 
        (keto-enol type) that proceeds through a suprafacial instead of the
        common antrafacial TS.
        '''
        
        name = 'intra_H_migration_suprafacial'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        # search for keto-enol type reactions
        motif = ['X', 'X', 'X', 'H']
        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        
        # filter for the double bond
        for instance in instances:
            if bond[instance[0]][instance[1]] == 2:
                rxns += [instance]

        self.new_reaction(rxns, name, a=0, b=-1)

#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
#                    new = 0
        
        return 0

    
    def search_intra_R_migration(self, natom, atom, bond, rad):
        ''' 
        This is an class that covers several RMG classes.
        
        R cannot be an H, this is already taken care of in the intra_H_migration
        
        TODO: merge this with intra H migration families?
        yes, because it is the same rule
        no, because then it's hard to search for just one of the types
        TODO: this should also include migration to lone pair electrons?
        currently it moves atoms to radical sites only
        '''
        
        if np.sum(rad) != 1: return 

        name = 'intra_R_migration'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        instances = []
        for ringsize in self.ringrange:
            motif = ['X' for i in range(ringsize)]
            for rad_site in np.nonzero(rad)[0]:
                instances += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)

        for instance in instances: 
            if not atom[instance[-1]] == 'H':
                rxns.append(instance)
        
        rxns = self.clean_rigid(name, rxns, 0, -1)

        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
#                    new = 0
        
        return 0

    def search_cpd_H_migration(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        H-C1-C=C-C=C-1 <==> C1=C-C=C-C(-H)-1

        '''
        
        if not any([len(ci) == 5 for ci in self.species.cycle_chain]) : return
        
        name = 'cpd_H_migration'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        bondsum = 0
        
        for cycle in self.species.cycle_chain:
            if len(cycle) == 5:
                for index, atomi in enumerate(cycle):
                    if index < 4:
                        atomj = cycle[index + 1]
                    else:
                        atomj = cycle[0]
                    if index == 0:
                        atomk = cycle[-1]
                    else:
                        atomk = cycle[index - 1]
                    bondsum += bond[atomi][atomj]
                    if bond[atomi][atomj] == 1 and bond[atomi][atomk] == 1:
                        start = atomi
                        startindex = index
                if bondsum != 7: return # exactly two double bonds
                ring_forw = np.ndarray.tolist(np.roll(cycle, 5 - startindex))
                ring_rev = ring_forw[::-1] # look at the ring in the reverse direction for an H-shift to the other side
                ring_rev = np.ndarray.tolist(np.roll(ring_rev, 1))
                rings = [ring_forw,ring_rev]
                
                Hatomi = -1
                for atomi in range(natom):
                    if atom[atomi] == 'H':
                        if bond[atomi][start] == 1:
                            Hatomi = atomi
                if Hatomi > -1:
                    for ring in rings:
                        instance = ring[:]
                        instance.append(Hatomi)
                        rxns += [instance]

        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                # TODO need to check if this is correct
#                if self.reac_bonds != {frozenset({inst[0], inst[1]})} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
#                    new = 0

        return 0
        

    def search_intra_OH_migration(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class extended.

        R*~~~~~~~O-OH <==> HOR~~~~~~~O*

        The H atom is not counted in the cycle size but has to be there.
        OH transfer to:
        radical sites
        double bonds on closed shell (just forward)
        '''
        
        name = 'intra_OH_migration'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        if np.sum(rad) == 0: 
        #find OH-migrations over double bonds and to lone pairs
            for ringsize in self.ringrange:
                # double bonds 
                motif = ['X' for i in range(ringsize)]
                motif[-1] = 'H'
                motif[-2] = 'O'
                motif[-3] = 'O'
                instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
           
                for instance in instances:
                    if any([bi > 1 for bi in bond[instance[0]]]):
                        rxns += [instance]

        else: 
            for ringsize in self.ringrange:
                instances = []
                # forward direction
                motif = ['X' for i in range(ringsize+1)]
                motif[-1] = 'H'
                motif[-2] = 'O'
                motif[-3] = 'O'
                for rad_site in np.nonzero(rad)[0]:
                    instances += find_motif.start_motif(motif, natom, bond, atom, 
                                                        rad_site, self.species.atom_eqv)
                # reverse direction
                motif = ['X' for i in range(ringsize+1)]
                motif[-1] = 'H'
                motif[-2] = 'O'
                motif[0] = 'O'
                for rad_site in np.nonzero(rad)[0]:
                    instances += find_motif.start_motif(motif, natom, bond, atom, 
                                                        rad_site, self.species.atom_eqv)
                for ins in instances:
                    rxns.append(ins)

        for case in range(len(rxns)):
            rxns[case] = rxns[case][:-1] #cut off H
            
        rxns = self.clean_rigid(name, rxns, 0, -1)

        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[-3], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-2]})}:
#                    new = 0

        return 0


    def search_intra_OH_migration_Exocyclic_F(self, natom, atom, bond, rad):
        ''' 
        This is the same as search_intra_OH_migration but for double bonds only

          0 .....-2-1
        R=R~~~~~~~O-OH <==> R=R~~?? + ??~~=O
                              |
                              OH

        The H atom is not counted in the cycle size but has to be there.
        OH transfer to double bonds on closed shell
        This is just the forward step as the expectation is that
        the product will fall apart while at least one of the framents 
        also rearrange, yielding two closed shell products.
        A special feature is to test both cis and trans transfer, therefore,
        in addition to testing for the extra H atom (which is deleted from the motif)
        the double bond is also registered and kept in the motif. 
        '''
        
        name = 'intra_OH_migration_Exocyclic_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        if np.sum(rad) == 0: 
            for ringsize in self.ringrange:
                # double bonds 
                motif = ['X' for i in range(ringsize)]
                motif[-1] = 'H'
                motif[-2] = 'O'
                motif[-3] = 'O'
                instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
           
                for instance in instances:
                    if bond[instance[0]][instance[1]] == 2:
                        rxns += [instance[:-1]]  # cut off H and add a -1 for nominal cis
                        rxns[-1].append(-1)
                        rxns += [instance[:-1]]  # cut off H and add a -2 for nominal trans
                        rxns[-1].append(-2)

        self.new_reaction(rxns, name, a=0, b=-1, c=-2)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[-3], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-2]})}:
#                    new = 0

        return 0


    def search_Intra_RH_Add_Endocyclic_F(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

                                  H
                                  | 
        H-R~~~~~~~R=R ==> R~~~~~~~R-R
                          |         |
                           ---------
        This is for the forward direction.
        '''
        
        if np.sum(rad) != 0: return
        if len(self.species.cycle_chain) > 0: return
        
        name = 'Intra_RH_Add_Endocyclic_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in range(5, 9):
            motif = ['X' for i in range(ringsize + 1)]
            motif[-1] = 'H'
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)

            bondpattern = ['X' for i in range(ringsize)]
            bondpattern[0] = 2
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance] 
            
        self.new_reaction(rxns, name, a=0, b=-2, length=True)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-2]}), frozenset({inst[-1], inst[1]})}:
#                    new = 0
                
        return 0


    def search_Intra_RH_Add_Endocyclic_R(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

                H
                | 
        R~~~~~~~R-R ==> H-R~~~~~~~R=R
        |         |
         ---------
        This is for the reverse direction.
        '''
        
        if len(self.species.cycle_chain) == 0: return
        if np.sum(rad) != 0: return
        
        name = 'Intra_RH_Add_Endocyclic_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ci in self.species.cycle_chain:
            motif = ['X' for i in range(len(ci) + 1)]
            motif[-1] = 'H'
            
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            
            # check if there is a bond between the first and second to last atom
            for instance in instances:
                if bond[instance[0]][instance[-2]] > 0:
                    rxns += [instance[-4:]]
        
        self.new_reaction(rxns, name, a=0, b=1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[-3], inst[-4]}), frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[-1], inst[-4]})}:
#                    new = 0
        
        self.reactions[name]

        return 0


    def search_Cyclic_Ether_Formation(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        R*~~~~~~~O-OR ==> R~~~~~~~O + OR
                          |_______|

        The OR groups are not counted in the cycle size but have to be there.
        Only the forward direction is included.
        '''
        
        
        if np.sum(rad) == 0: return
        
        name = 'Cyclic_Ether_Formation'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(4, 10):
            motif = ['X' for i in range(ringsize)]
            motif[-2] = 'O'
            motif[-3] = 'O'
            motif[0] = 'C'
            for rad_site in np.nonzero(rad)[0]:
                rxns += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)

        for instance in range(len(rxns)):
            rxns[instance] = rxns[instance][:-2] #cut off OR
            
        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[-2], inst[-3]})} or self.prod_bonds != {frozenset({inst[0], inst[-3]})}:
#                    new = 0
        
        return 0


    def search_Intra_R_Add_Endocyclic_F(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        *R~~~~~~~~R=R ==> R~~~~~~~~R*-R
                          |___________|

        '''
        
        if np.sum(rad) == 0: return

        name = 'Intra_R_Add_Endocyclic_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in self.ringrange:
            motif = ['X' for i in range(ringsize)]
            instances = []
            for rad_site in np.nonzero(rad)[0]:
                instances += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)
            bondpattern = ['X' for i in range(ringsize-1)]
            bondpattern[-1] = 2
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance]
            
            bondpattern[-1] = 3
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance]

        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset()} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
#                    new = 0
                
        return 0


    def search_Intra_R_Add_ExoTetCyclic_F(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        *R~~~~~~~~R-R ==> R~~~~~~~~R + R*
                          |________|

        '''

        if np.sum(rad) == 0: return

        name = 'Intra_R_Add_ExoTetCyclic_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in self.ringrange:
            motif = ['X' for i in range(ringsize + 1)]
            for rad_site in np.nonzero(rad)[0]:
                rxns += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)

        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-2]})}:
#                    new = 0

        return 0


    def search_Intra_R_Add_Exocyclic_F(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        *R~~~~~~~~R=R ==> R~~~~~~~~R-R*
                          |________|

        '''
        
        if np.sum(rad) == 0: return

        name = 'Intra_R_Add_Exocyclic_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in self.ringrange:
            motif = ['X' for i in range(ringsize + 1)]
            instances = []
            for rad_site in np.nonzero(rad)[0]:
                instances += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)
            bondpattern = ['X' for i in range(ringsize)]
            bondpattern[-1] = 2
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance]
                    
            bondpattern[-1] = 3
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance]

        self.new_reaction(rxns, name, a=0, b=-2)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset()} or self.prod_bonds != {frozenset({inst[0], inst[-2]})}:
#                    new = 0

        return 0


    def search_Intra_RH_Add_Exocyclic_F(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        The general scheme is:

        H-R~~~~~~R=R ==> R~~~~~~R-R-H
                         |      |
                          ------

        The special case of this reaction is Korcel_step1:
              
                R        R   OH 
        R      /          \ /
         \ / \C=O          C
          |       ==>     / \
          O   H          |   O
           \ /          / \ /
            O          R   O

        Implemented as:

                                   --O--O--
                                  |        |
        O=C~~~~~~~~C-O-O-H ==> HO-C~~~~~~~~C
          |                       |
          R                       R

        The carbonyl dangling R and the
        tail H are included, but are not counted as the ring size, but these two atoms are kept
        because they are needed in the geometry manipulation step.
        '''
        
        if len(self.species.cycle_chain) > 0: return
        
        name = 'Intra_RH_Add_Exocyclic_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in self.ringrange:
            motif = ['X' for i in range(ringsize+2)]
            motif[-1] = 'H'

            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            bondpattern = ['X' for i in range(ringsize+1)]
            bondpattern[0] = 2
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance]

        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[1], inst[-2]}), frozenset({inst[-1], inst[0]})}:
#                    new = 0
        
        return 0


    def search_Intra_RH_Add_Exocyclic_R(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

                                    H -1
                          1       2 | 
        H-R~~~~~~~R=R <== R~~~~~~~R-R
                          |_______|

        This is for the reverse direction.
        '''
        
        name = 'Intra_RH_Add_Exocyclic_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer


        if len(self.species.cycle_chain) == 0: return
        if np.sum(rad) != 0: return

        for ci in self.species.cycle_chain:
            motif = ['X' for i in range(len(ci) + 2)]
            motif[-1] = 'H'
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)

            # check if there is a bond between the first and second to last atom
            for instance in instances:
                if bond[instance[0]][instance[-3]] > 0:
                    rxns += [instance[-4:]]

        self.new_reaction(rxns, name, a=0, b=1, c=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[-1], inst[-2]}), frozenset({inst[-3], inst[0]})} or self.prod_bonds != {frozenset({inst[-1], inst[0]})}:
#                    new = 0

        return 0


    def search_Retro_Ene(self, natom, atom, bond, rad):
        ''' 
        This is not an RMG class.

        R-R-R-R=R ==> R=R + R=R-R

        '''
        
        
        if np.sum(rad) != 0: return
        
        name = 'Retro_Ene'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X' for i in range(6)]
        motif[-1] = 'H'
        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)

        bondpattern = ['X' for i in range(5)]
        bondpattern[0] = 2
        for instance in instances:
            if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                rxns += [instance] 
        
        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[2], inst[3]})} or self.prod_bonds != {frozenset()}:
#                    new = 0
        
        return 0


    def search_Korcek_step2_odd(self, natom, atom, bond, rad):
        ''' 
        Korcek step 2 for cyclic peroxides originating with odd number of atoms in the cycle.
        Ring breaks at O-O and then forms 1 three-ringatom and (ringsize-3)/2 two-ringatom 
        fragments. 
        The three-ringatom fragment needs to have a H atom transfer. 
        Numbering:
        0 1 2 3 4 5 6 7 8
        O-X-X-X-X-X-X-X-O
        if the fragment is set to 2 2 3 2, then it breaks like
        O-X X-X X-X-X X-X
        and the X-X-X will have 0, 1, or 2 possible H transfers. 
        If the H transfer is not possible, the path is abandoned.
        The instance name is:
            all atoms in the chain, the middle atom in the triplet,
            and the atom to which the H is migrates in the triplet
            and the hydrogen itself
        '''

        name = 'Korcek_step2_odd'

        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in range(5, 15, 2):  # odd number of atoms in the ring
            motif = ['X' for i in range(ringsize)]
            motif[-1] = 'O'
            motif[0] = 'O'
            korcek_chain = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            # filter clockwise and anti clockwise hits
            korcek_chain_filt = []
            for kch in korcek_chain:
                k = copy.deepcopy(kch)  # need in order to prevent changes to korcek_chain with reverse()
                l = copy.deepcopy(kch)
                l.reverse()
                if k not in korcek_chain_filt and l not in korcek_chain_filt:
                    korcek_chain_filt.append(kch)

            for ins in korcek_chain_filt:
                if bond[ins[0]][ins[-1]] == 1:  # it is a ring
                    fragment = [2] * int((ringsize - 3) / 2 + 1)
                    fragment[0] = 3  # [3, 2, 2, ...]
                    for ii in range(len(fragment)):  # loop over all possible 2/3 fragmentation
                        threefrag = ins[ii * 2 : ii * 2 + 3]  # atoms in the 3-long fragment
                        for at in range(natom):
                            if bond[threefrag[1]][at] == 1 and atom[at] == 'H':  # there is H on the middle atom
                                # if there are 2 hydrogens, they are treated separately, as they are not 
                                # in general equivalent due to the ring
                                ins_full = ins + [threefrag[1]] + [threefrag[0]] + [at] # H adds to the first atom of this fragment
                                rxns += [ins_full]
                                ins_full = ins + [threefrag[1]] + [threefrag[2]] + [at] # H adds to the second atom of this fragment
                                rxns += [ins_full]

        self.new_reaction(rxns, name, full=True)
#        for n, inst in enumerate(rxns):
#            new = 1
#            #filter for the same reactions
#            for instance in self.reactions[name]:
#                if inst == instance:
#                    new = 0
#            # filter for specific reaction after this # TODO
#            #if self.one_reaction_fam and new:
#            #    if ring_var[n] == 7: 
#            #        if (not {frozenset({inst[-2], inst[-3]}), frozenset({inst[0], inst[1]})}.issubset(self.reac_bonds)) or self.prod_bonds != {frozenset()}:
#            #            new = 0
#            #    if ring_var[n] == 8: 
#            #        #  TODO this is an incomplete check
#            #        if self.reac_bonds != {frozenset({inst[-2], inst[-3]}), frozenset({inst[-4], inst[-5]}), frozenset({inst[0], inst[1]})}:
#            #            new = 0
#            if new:
#                self.reactions[name].append(inst)

        return 0


    def search_Korcek_step2_even(self, natom, atom, bond, rad):
        ''' 
        Korcek step 2 for cyclic peroxides with even number of atoms in the ring.
        Still, the 4 membered ring equals a 2,2 cycloaddition and is not considered here.
        Ring breaks at O-O and then at every second bond, no H shift is needed.
        '''

        name = 'Korcek_step2_even'

        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in range(6, 14, 2):  # even number of atoms in the ring
            motif = ['X' for i in range(ringsize)]
            motif[-1] = 'O'
            motif[0] = 'O'
            korcek_chain =  find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            # filter clockwise and anti clockwise hits
            korcek_chain_filt = []
            for kch in korcek_chain:
                k = copy.deepcopy(kch)  # need in order to prevent changes to korcek_chain with reverse()
                l = copy.deepcopy(kch)
                l.reverse()
                if k not in korcek_chain_filt and l not in korcek_chain_filt:
                    korcek_chain_filt.append(kch)

            for ins in korcek_chain_filt:
                if bond[ins[0]][ins[-1]] == 1:  # it is a ring
                    rxns += [ins]

        self.new_reaction(rxns, name, full=True)
#        for n, inst in enumerate(rxns):
#            new = 1
#            #filter for the same reactions
#            for instance in self.reactions[name]:
#                if inst == instance:
#                    new = 0
#            # filter for specific reaction after this
#            #if self.one_reaction_fam and new:
#            #    if ring_var[n] == 7: 
#            #        if (not {frozenset({inst[-2], inst[-3]}), frozenset({inst[0], inst[1]})}.issubset(self.reac_bonds)) or self.prod_bonds != {frozenset()}:
#            #            new = 0
#            #    if ring_var[n] == 8: 
#            #        #  TODO this is an incomplete check
#            #        if self.reac_bonds != {frozenset({inst[-2], inst[-3]}), frozenset({inst[-4], inst[-5]}), frozenset({inst[0], inst[1]})}:
#            #            new = 0
#            if new:
#                self.reactions[name].append(inst)

        return 0


    def search_Korcek_step2(self, natom, atom, bond, rad):
        ''' 
        Generalized Korcek step 
        
        The 4 membered ring equals a 2,2 cycloaddition and is not considered here (no H shift involved)
        
        The 5 membered ring proceeds through a 6 membered transition state (including a 1,2 H migration):

            --O--O--
           |        |
        HO-C---C----C-R  ==> RCOOH + R3CC(R)O
           |  / \   |
           R R   R  R

        6-membered ring: TODO

        Only the forward direction is included.

        '''
        
        
        name = 'Korcek_step2'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        ring_var = [] #  a helper variable to temporarily mark the ring size within this function

        for ringsize in range(5, 6):
            motif = ['X' for i in range(ringsize + 1)]
            #motif[-1] = 'H'  #  deleted because atom types are no longer checked
            korcek_chain =  find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            for ins in korcek_chain:
                if bond[ins[0]][ins[-2]] == 1:
                    rxns += [ins]
                    ring_var.append(ringsize)

        self.new_reaction(rxns, name, a=0, b=-1)
#        for n, inst in enumerate(rxns):
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if ring_var[n] == 7: 
#                    if (not {frozenset({inst[-2], inst[-3]}), frozenset({inst[0], inst[1]})}.issubset(self.reac_bonds)) or self.prod_bonds != {frozenset()}:
#                        new = 0
#                if ring_var[n] == 8: 
#                    #  TODO this is an incomplete check
#                    if self.reac_bonds != {frozenset({inst[-2], inst[-3]}), frozenset({inst[-4], inst[-5]}), frozenset({inst[0], inst[1]})}:
#                        new = 0
        
        return 0


    def search_r22_cycloaddition(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        R      R         R---R
        ||  +  ||  <==   |   |
        R      R         R---R

        N.B.: only the reverse direction is available. Also, the 3 related RMG classes are treated as one.

        '''
        
        
        name = 'r22_cycloaddition'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        if not any([len(ci) == 4 for ci in self.species.cycle_chain]): return 
        
        for ci in self.species.cycle_chain:
            if len(ci) == 4:
                # there are two ways to slice a 4-mem ring
                ring1 = ci
                ring2 = np.ndarray.tolist(np.roll(ring1, 1))

                # FIXME only works for 1 cycle
                rxns += [ring1]
                rxns += [ring2]

        self.new_reaction(rxns, name, a=0, b=1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                # TODO need to make sure that these are the bonds that are broken, see the reaction details
#                if self.reac_bonds != {frozenset({inst[0], inst[1]}), frozenset({inst[2], inst[3]})} or self.prod_bonds != {frozenset()}:
#                    new = 0
        
        return 0



    def search_r12_cycloaddition(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

                       R--R
        R=R + R: <==   \  /
                        R 

        N.B.: only the reverse direction is available. 

        '''
        
        
        name = 'r12_cycloaddition'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        if not any([len(ci) == 3 for ci in self.species.cycle_chain]):
            return
        
        for ci in self.species.cycle_chain:
            if len(ci) == 3:
                # there are three ways to slice a 3-mem ring
                ring1 = ci
                ring2 = np.ndarray.tolist(np.roll(ring1, 1))
                ring3 = np.ndarray.tolist(np.roll(ring1, 2))

                rxns += [ring1]
                rxns += [ring2]
                rxns += [ring3]

        self.new_reaction(rxns, name, a=0, b=1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                # TODO need to make sure that these are the bonds that are broken, see the reaction details
#                if self.reac_bonds != {frozenset({inst[0], inst[2]}), frozenset({inst[1], inst[2]})} or self.prod_bonds != {frozenset()}:
#                    new = 0
        
        return 0



    def search_r12_insertion_R(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

                          X
                          |
        X-P + R-R <==   R-P-R

        '''
        
        #if np.sum(rad) != 0: return
        
        name = 'r12_insertion_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        motif = ['X','X','X']
        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        
        for instance in instances:
            #if all([atom[atomi] != 'H' for atomi in instance]):
            rxns += [instance]

        self.new_reaction(rxns, name, a=0, b=1, c=2)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != set({frozenset({inst[0], inst[1]}), frozenset({inst[1], inst[2]})}) or self.prod_bonds != {frozenset({inst[0], inst[2]})}:
#                    new = 0
        
        return 0


    def search_r13_insertion_CO2(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

                          O
                          ||
        O=C=O + R-R <== R-C-O-R


        '''
        
        #if np.sum(rad) != 0: return
        
        name = 'r13_insertion_CO2'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        motif = ['X','C','O','X']
        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        for instance in instances:
            for atomi in range(natom):
                if not atomi in instance:
                    if atom[atomi] == 'O':
                        if bond[atomi][instance[1]] == 2:
                            rxns += [instance]

        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[0], inst[1]}), frozenset({inst[2], inst[3]})} or self.prod_bonds != {frozenset({inst[0], inst[3]})}:
#                    new = 0
        
        return 0


    def search_r13_insertion_ROR(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        R1-O-R2 + R=R <== R1-R-R-O-R2


        '''
        
        #if np.sum(rad) != 0: return
        name = 'r13_insertion_ROR'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X','X','X','O']
        rxns = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        
        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != set({frozenset({inst[0], inst[1]}), frozenset({inst[2], inst[3]})}) or self.prod_bonds != {frozenset({inst[0], inst[3]})}:
#                    new = 0
        
        return 0


    def search_Diels_alder_addition(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

          R                  R
        //                 /   \
        R       R         R     R
        |  +    ||  <==   ||    |
        R       R         R     R
         \\                \   /
           R                 R

        N.B.: only the reverse direction is available. 

        '''
        
        
        name = 'Diels_alder_addition'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        if not any([len(ci) == 6 for ci in self.species.cycle_chain]): return 

        for ci in self.species.cycle_chain:
            if len(ci) == 6:
                bondsum = 0
                for index, atomi in enumerate(ci):
                    if index < 5:
                        atomj = ci[index + 1]
                    else:
                        atomj = ci[0]
                    bondsum += bond[atomi][atomj]
                    if bond[atomi][atomj] == 2:
                        start = atomi
                        startindex = index
                if bondsum != 7: return # exactly one double bond
                ring = np.ndarray.tolist(np.roll(ci, 6 - startindex))

                rxns += [ring] # FIXME only works for 1 cycle

        self.new_reaction(rxns, name, a=0, b=1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != set({frozenset({inst[2], inst[3]}), frozenset({inst[4], inst[5]})}) or self.prod_bonds != {frozenset()}:
#                    new = 0
        
        return 0



    def search_Intra_Diels_alder_R(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.
        The reaction is a lot more general now, simply requires multiple bonds that can add to each other.
        No middle double bond is required.
        If there is a middle double bond, there will be a just a rearrangement of bonds.
        If there is no middle double bond, that bond will break up.
        At the moment the family is only meaningful for 6-mem rings, 
        behavior for other ring sizes is unknown. Although this is the DA reaction.

                             C
                            / \\
        1 2 3 4  -2-1      C   C
        C=C-C=C~~~C=C ==>  |   |
                           C   C
                            \ //
                              C

        '''
        
        name = 'Intra_Diels_alder_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in self.ringrange:  # TODO what is the meaning of these larger rings?
            motif = ['X' for i in range(ringsize)]
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)

            bondpattern = ['X' for i in range(ringsize - 1)]
            bondpattern[0] = 2
            #bondpattern[2] = 2  # by removing this, this family becomes a lot more general
            bondpattern[-1] = 2

            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern, atleast=True) == 0:
                    #inst = instance[:4] + instance[-2:]
                    rxns += [instance] 
     

        self.new_reaction(rxns, name, a=0, b=-1, cross=True)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset()} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
#                    new = 0
                
        return 0


    def search_ketoenol(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        R=R-O-R1 <==> R1-R-R=O
        '''
        
        name = 'ketoenol'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        # enol to keto
        motif = ['C', 'C', 'O', 'X']
        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)

        # keto to enol
        motif = ['O', 'C', 'C', 'X']
        instances += find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        bondpattern = [2, 'X', 'X', 'X']
        for instance in instances:
            if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                rxns += [instance]
            
        self.new_reaction(rxns, name, a=0, b=1, c=2, d=3)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[2], inst[3]})} or self.prod_bonds != {frozenset({inst[0], inst[1]})}:
#                    new = 0
        
        return 0
 


    def search_HO2_Elimination_from_PeroxyRadical(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        H-R-R-O-O* ==> R=R + HO2

        N.B.: only the forward direction is available.
        '''
        
        
        if np.sum(rad) == 0: return
        
        name = 'HO2_Elimination_from_PeroxyRadical'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['H', 'X', 'X', 'O', 'O']
        rxns += find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            
        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != set({frozenset({inst[0], inst[1]}), frozenset({inst[2], inst[3]})}) or self.prod_bonds != {frozenset({inst[0], inst[4]})}:
#                    new = 0
        
        return 0
        

    def search_R_Addition_COm3_R(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        C#O + R* <== R-C*=O

        N.B.: only the reverse direction is available. 
        '''
        
        
        if np.sum(rad) == 0: return
        
        name = 'R_Addition_COm3_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X', 'C', 'O']

        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)

        for instance in instances:
            bondpattern = [1, 2]
            if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                if rad[instance[1]] == 1:
                    rxns += [instance]
        
        self.new_reaction(rxns, name, a=0, b=1, c=2)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[0], inst[1]})} or self.prod_bonds != {frozenset()}:
#                    new = 0

        return 0


        
    def search_R_Addition_MultipleBond(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        R=R + r* <== R*-R-r 
                     0  1 2

        - only the reverse direction is available. 
        - this is not a scan class
        - we are also allowing to form anticipated resonance stabilized species:
            R=R-R-r ==> [R=R-R <--> R-R=R] + r 
        '''
        
        
        name = 'R_Addition_MultipleBond'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        # simple beta scission for radicals
        motif = ['X', 'X', 'X']
        for rad_site in np.nonzero(rad)[0]:
            rxns += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)

        # anticipated resonance stabilized radical
        motif = ['X', 'X', 'X', 'X']
        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        bondpattern = [2, 'X', 'X', 'X']
        for instance in instances:
            if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                rxns += [instance[1:]]  # cutting off the first atom
        bondpattern = [3, 'X', 'X', 'X']
        for instance in instances:
            if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                rxns += [instance[1:]]
    
        self.new_reaction(rxns, name, a=0, b=1, c=2)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[1], inst[2]})} or self.prod_bonds != {frozenset()}:
#                    new = 0
        
        return 0


    def search_12_shift_S_F(self, natom, atom, bond, rad):
        '''
        This is an RMG class.
        '''

        if np.sum(rad) != 1: return
        
        name = '12_shift_S_F'
        
        if not name in self.reactions:
            self.reactions[name] = []
        
        motif = ['X','S','X']
        rxns = []
        for rad_site in np.nonzero(rad)[0]:
            rxns += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)

        #filter for identical reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1] and inst[2] == instance[2]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[1], inst[2]})} or self.prod_bonds != {frozenset()}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        return 0


    def search_12_shift_S_R(self, natom, atom, bond, rad):
        '''
        This is an RMG class.

        C-S-R* <== *S-R-C

        TODO: why not forward??

        '''
        
        if np.sum(rad) != 1: return
        
        name = '12_shift_S_R'
        
        if not name in self.reactions:
            self.reactions[name] = []
        
        motif = ['S','X','X']
        rxns = []
        for rad_site in np.nonzero(rad)[0]:
            rxns += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)
        
        for inst in rxns:
            new = 1
            # filter for identical reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1] and inst[2] == instance[2]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[0], inst[1]})} or self.prod_bonds != {frozenset({inst[0], inst[2]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        return 0


    def search_r13_insertion_RSR(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        R-S-R + R1=R2 <== R-R1-R2-S-R


        '''
        
        #if np.sum(rad) != 0: return
        name = 'r13_insertion_RSR'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X','X','X','S']
        rxns = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        

        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != set({frozenset({inst[0], inst[1]}), frozenset({inst[0], inst[1]})}) or self.prod_bonds != {frozenset({inst[0], inst[3]})}:
#                    new = 0
        
        return 0


    def search_R_Addition_CSm_R(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        C#S + R* <== R-C*=S

        N.B.: only the reverse direction is available. 
        '''
        
        
        if np.sum(rad) == 0: return
        
        name = 'R_Addition_CSm_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X', 'C', 'S']

        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)

        for instance in instances:
            bondpattern = [1, 2]
            if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                if rad[instance[1]] == 1:
                    rxns += [instance]
        
        self.new_reaction(rxns, name, a=0, b=1, c=2)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[0], inst[1]})} or self.prod_bonds != {frozenset()}:
#                    new = 0

        return 0


    def search_r14_birad_scission(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        It is now renamed to 1,4_Linear_birad_scission on the RMG website,

        *R-R-R-R* ==> R=R + R=R

        Problematic reaction because of the biradical character.

        '''

        if np.sum(rad) != 2: return
        
        
        name = 'r14_birad_scission'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        
        motif = ['X','X','X','X']
        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        for instance in instances: 
            if rad[instance[0]] == 1 and rad[instance[-1]] == 1:
                rxns += [instance]

        self.new_reaction(rxns, name, a=1, b=2)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[1], inst[2]})} or self.prod_bonds != {frozenset()}:
#                    new = 0
        
        return 0


    def search_r14_cyclic_birad_scission_R(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        R1-R*~~~~~~R*-R2   <==  R1=R~~~~~~R=R2
        |______________|
        (this one bond)

        TODO forward?

        '''

        if np.sum(rad) != 0: return
        
        name = 'r14_cyclic_birad_scission_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in range(5, 9):
            motif = ['X' for i in range(ringsize)]
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            
            bondpattern = ['X' for i in range(ringsize - 1)]
            bondpattern[0] = 2
            bondpattern[-1] = 2
            
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance] 
                    
        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset()} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
#                    new = 0
        
        return 0


    def search_birad_recombination_F(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        *R~~~~~~~~R* ==> R~~~~~~~~R
                         |________|

        '''

        if np.sum(rad) != 2: return
        
        name = 'birad_recombination_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in self.ringrange:
            motif = ['X' for i in range(ringsize)]
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
           
            for instance in instances: 
                if rad[instance[0]] == 1 and rad[instance[-1]] == 1:
                    rxns += [instance]

        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset()} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
#                    new = 0
        
        return 0


    def search_birad_recombination_R(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        *R~~~~~~~~R* <== R~~~~~~~~R
                         |________|

        '''

        if np.sum(rad) != 0: return
        if len(self.species.cycle_chain) == 0: return
        
        name = 'birad_recombination_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X','X']
        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        for instance in instances: 
            if instance[0] in self.cycle and instance[1] in self.cycle :
                rxns += [instance]

        self.new_reaction(rxns, name, a=0, b=1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[0], inst[1]})} or self.prod_bonds != {frozenset()}:
#                    new = 0
        
        return 0


    def search_Intra_disproportionation_F(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        *R~~~~~R*-R-H ==> H-R~~~~~R=R

        '''

        if np.sum(rad) != 2: return
        
        name = 'Intra_disproportionation_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in range(5, 9):
            motif = ['X' for i in range(ringsize)]
            motif[-1] = 'H'
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
           
            for instance in instances: 
                if rad[instance[0]] == 1 and rad[instance[-3]] == 1:
                    rxns += [instance]

        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
#                    new = 0
        
        return 0


    def search_Intra_disproportionation_R(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.

        *R~~~~~R*-R-H <== H-R~~~~~R=R

        '''

        if np.sum(rad) != 0: return
        
        name = 'Intra_disproportionation_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(5, 9):
            motif = ['X' for i in range(ringsize)]
            motif[-1] = 'H'
            
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            
            bondpattern = ['X' for i in range(ringsize - 1)]
            bondpattern[0] = 2
            
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance] 

        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
#                    new = 0
        
        return 0


    def search_bimol_disproportionation_R(self, natom, atom, bond, rad):
        ''' 
        This is an RMG class.
          X                  X
          |                  |
        R=R  *R*-R-H <== H-R-R-R=R
          |   |              | | 
          Y   Z              Y Z
        '''

        if np.sum(rad) != 0: return
        
        name = 'bimol_disproportionation_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(5, 6):
            motif = ['X' for i in range(ringsize)]
            motif[-1] = 'H'
            
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            
            bondpattern = ['X' for i in range(ringsize - 1)]
            bondpattern[0] = 2
            
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance] 

        self.new_reaction(rxns, name, a=0, b=-1)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
#                    new = 0
        
        return 0


    def search_beta_delta(self, natom, atom, bond, rad):
        '''
        This is not an RMG class.

        A*-B-C-D-E ==> A=B + C=D + E* 

        It is the parallel breaking of not just the beta but also of the gamma bond, resulting in two unsaturated bonds and a radical.
        '''


        if np.sum(rad) == 0: return

        name = 'beta_delta'

        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        motif = ['X', 'X', 'X', 'X', 'X']
        for rad_site in np.nonzero(rad)[0]:
            rxns += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)

        self.new_reaction(rxns, name, a=0, b=1, c=2, d=3, e=4)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[1], inst[2]}), frozenset({inst[3], inst[4]})} or self.prod_bonds != {frozenset()}:
#                    new = 0

        return 0


    def search_h2_elim(self, natom, atom, bond, rad):
        ''' 
        This is not an RMG class. Now extended for general H2 elimination, not just 1,2


        H   H
        |   |
        X - X ==> X=X  + H2 

        '''

        name = 'h2_elim'

        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

#        motif = ['H','X','X','H']
#        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
#        for instance in instances: 
#            rxns += [instance]

#        self.new_reaction(rxns, name, a=0, b=-1, cross=True)

        for ringsize in self.ringrange:
            motif = ['X' for i in range(ringsize)]
            motif[0] = 'H'
            motif[-1] = 'H'
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            for instance in instances: 
                rxns += [instance]

        rxns = self.clean_rigid(name, rxns, 0, -1)
        self.new_reaction(rxns, name, a=0, b=-1, cross=True)

#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[0], inst[1]}), frozenset({inst[2], inst[3]})} or self.prod_bonds != {frozenset({inst[0], inst[3]})}:
#                    new = 0

        return 0


    def search_hom_sci(self, natom, atom, bond, rad):
        ''' 
        This is not an RMG class.

        R-R ==> R + R
        
        We exclude ring bonds.
        '''

        name = 'hom_sci'

        if not name in self.reactions:
            self.reactions[name] = []

        rxns = []  # reactions found with the current resonance isomer

        if self.par['homolytic_bonds'] == {}:
            motif = ['X','X']
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            for instance in instances: 
                if not self.species.cycle[instance[0]] or not self.species.cycle[instance[1]]:
                    rxns += [instance]
        else: 
            try:
                rxns = self.par['homolytic_bonds'][str(self.species.chemid)]
            except KeyError:
                pass
                

        self.new_reaction(rxns, name, a=0, b=1, cross=True, aid=True)
#            # filter for specific reaction after this
#            if self.one_reaction_fam and new:
#                if self.reac_bonds != {frozenset({inst[0], inst[1]}), frozenset({inst[2], inst[3]})} or self.prod_bonds != {frozenset({inst[0], inst[3]})}:
#                    new = 0

        return 0


    def search_barrierless_saddle(self, natom, atom, bond, rad):
        ''' 
        This is not an RMG class.

        R - R ==> R + R

        Attempts to find a saddle point for a nominally barrierless reaction.
        '''

        name = 'barrierless_saddle'

        if not name in self.reactions:
            self.reactions[name] = []
        if self.barrierless_saddle is not None:
            rxns = self.barrierless_saddle  # defined by the user
        else:
            return 0

        self.new_reaction(rxns, name, a=0, b=-1, cross=True)

        return 0


    def reaction_matrix(self, reac_list, reac_id):
        ''' 
        Create arrays to store all reactions for species.
        input: 
        reac_list: atom motifs from individual searches
        reac_id: reaction name (e.g., HO2_Elimination_from_PeroxyRadical) from individual search functions
        Every reaction type just makes the below arrays longer, generated as reactions are found.

        generated:
        reac_type: reaction class identifier
        reac_inst: reaction instance defined by the important atoms
        reac_step: the step at which the search is at
        reac_scan_energy: for each reaction the energy as a function of steps, only used for scanning type searches, e.g. R_Addition_MultipleBond
        rec_ts_done: the last calculations is submitted in the sequence
        reac_ts_geom: the geometry of the TS
        reac_ts_freq: the freqencies of the TS
        reac_name: the base name of the file to run - created for each reaction later
        '''
        
        self.species.reac_type += [reac_id for i in range(len(reac_list))]
        self.species.reac_inst += reac_list
        self.species.reac_step += [0 for i in range(len(reac_list))]
        self.species.reac_scan_energy += [[] for i in range(len(reac_list))]
        self.species.reac_ts_done += [0 for i in range(len(reac_list))] 
        self.species.reac_ts_geom += [0 for i in range(len(reac_list))]
        self.species.reac_ts_freq += [0 for i in range(len(reac_list))]
        
        for i in range(len(reac_list)):
            if reac_id == 'intra_H_migration':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraHMigration(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'intra_H_migration_suprafacial':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraHMigrationSuprafacial(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'intra_R_migration':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraRMigration(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'intra_OH_migration':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraOHMigration(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'intra_OH_migration_Exocyclic_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][-2] + 1) + '_' + str(reac_list[i][-1])  # last element is cis/trans (-1, -2)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraOHMigrationExocyclicF(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'cpd_H_migration':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) + '_' + str(reac_list[i][-2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(CpdHMigration(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'Intra_RH_Add_Endocyclic_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(len(reac_list[i])) + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraRHAddEndoF(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'Intra_RH_Add_Endocyclic_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraRHAddEndoR(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'Cyclic_Ether_Formation':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(CyclicEtherFormation(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'Intra_RH_Add_Exocyclic_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraRHAddExoF(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'Intra_RH_Add_Exocyclic_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraRHAddExoR(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'Retro_Ene':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(RetroEne(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'Intra_R_Add_Endocyclic_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraRAddEndocyclicF(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'Intra_R_Add_ExoTetCyclic_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-2] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraRAddExoTetCyclicF(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'Intra_R_Add_Exocyclic_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraRAddExocyclicF(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'Korcek_step2_odd':
                name = str(self.species.chemid) + '_' + reac_id
                for j in range(len(reac_list[i])):
                    name += '_' + str(reac_list[i][j] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(KorcekStep2Odd(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'Korcek_step2_even':
                name = str(self.species.chemid) + '_' + reac_id
                for j in range(len(reac_list[i])):
                    name += '_' + str(reac_list[i][j] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(KorcekStep2Even(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'Korcek_step2':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(KorcekStep2(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'r22_cycloaddition':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(R22Cycloaddition(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'r12_cycloaddition':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(R12Cycloaddition(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'r12_insertion_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(R12Insertion(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'r13_insertion_CO2':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(R13InsertionCO2(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'r13_insertion_ROR':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1) + '_' + str(reac_list[i][3] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(R13InsertionROR(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'r14_birad_scission':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(R14BiradScission(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'r14_cyclic_birad_scission_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(R14CyclicBiradScission(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'birad_recombination_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(BiradRecombinationF(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'birad_recombination_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(BiradRecombinationR(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'Intra_disproportionation_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraDisproportionationF(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'Intra_disproportionation_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraDisproportionationR(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'Diels_alder_addition':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(DielsAlder(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'Intra_Diels_alder_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraDielsAlder(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'ketoenol':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1)  + '_' + str(reac_list[i][2] + 1)  + '_' + str(reac_list[i][3] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(KetoEnol(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'HO2_Elimination_from_PeroxyRadical':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(HO2Elimination(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'R_Addition_COm3_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(RAdditionCO(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'R_Addition_MultipleBond':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(RAdditionMultipleBond(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == '12_shift_S_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(S12ShiftF(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == '12_shift_S_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(S12ShiftR(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'R_Addition_CSm_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(RAdditionCS(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'r13_insertion_RSR':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1) + '_' + str(reac_list[i][3] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(R13InsertionRSR(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'beta_delta':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1) + '_' + str(reac_list[i][3] + 1) + '_' + str(reac_list[i][4] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(BetaDelta(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'h2_elim':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(H2Elim(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'hom_sci':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(HS(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'bimol_disproportionation_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(BimolDisproportionationR(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'barrierless_saddle':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(BarrierlessSaddle(self.species, self.qc, self.par, reac_list[i], name))
            elif reac_id == 'combinatorial':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(i)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(Combinatorial(self.species, self.qc, self.par, reac_list[i], name))
            else:
                self.species.reac_name.append(0)
        return 0


    def clean_rigid(self, name, instances, pivot1, pivot2):
        '''
        Getting rid of instances where the rigid structure would not allow the 
        transfer of atoms, e.g., H transfer across a large rigid ring structure.
        It is based on the presence of (partial) double bonds along the motif.
        If the structure is rigid, and the selected pivot atoms are further than a cutoff
        then the instance will be deleted from the list.
        Pivots requires manual determination for each family, where this is important.
        Not applied to all families.
        '''

        #cutoff = 3.  # Angstrom
        mask = [True] * len(instances)
        try:
            self.species.maxbond
        except AttributeError:
            self.species.calc_maxbond()
        for inst, instance in enumerate(instances):
            if all(self.species.maxbond[instance[ii]][instance[ii + 1]] > 1 for ii in range(len(instance) - 2)):
                if np.linalg.norm(self.species.geom[instance[pivot1]] - self.species.geom[instance[pivot2]]) > self.par['rigid_reaction_cutoff']:
                    mask[inst] = False
                    numbers = [ii + 1 for ii in instance]
                    logger.info(f'{name} reaction {numbers} over rigid backbone with cutoff {self.par["rigid_reaction_cutoff"]} A is removed.')
        return list(np.array(instances, dtype=object)[mask])


    def new_reaction(self, rxns, name, a=None, b=None, c=None, d=None, e=None, 
                     length=None, full=False, cross=False, aid=False):
        '''
        Checks a variable number of identical elements
        Also can check full equivalency (full=True), same lenght (length=True), and 
        equivalency between elements that are interchangeable (cross=True)
        if aid is True, then it will throw away reactions where there is already one
           with the same atom IDs involved - at least important for hom_sci
        now also filters non-solute reactions in cluster mode
        '''

        for inst in rxns:
            new = True
            for instance in self.reactions[name]:
                if aid == True:
                    if (self.species.atomid[inst[a]] == self.species.atomid[instance[a]] and
                            self.species.atomid[inst[b]] == self.species.atomid[instance[b]]):
                        new = False
                        break
                    if (self.species.atomid[inst[b]] == self.species.atomid[instance[a]] and
                            self.species.atomid[inst[a]] == self.species.atomid[instance[b]]):
                        new = False
                        break
                if cross == True:
                    if (inst[a] == instance[a] and inst[b] == instance[b]):
                        new = False
                        break 
                    if (inst[a] == instance[b] and inst[b] == instance[a]):
                        new = False
                        break 
                if a is not None:
                    if inst[a] != instance[a]:
                        continue
                if b is not None:
                    if inst[b] != instance[b]:
                        continue
                if c is not None:
                    if inst[c] != instance[c]:
                        continue
                if d is not None:
                    if inst[d] != instance[d]:
                        continue
                if e is not None:
                    if inst[e] != instance[e]:
                        continue
                if length is not None:
                    if len(inst) != len(instance):
                        continue
                if full == True:
                    if len(inst) != len(instance):
                        continue
                    if any([inst[i] != instance[i] for i, _ in enumerate(inst)]):
                        continue
                new = False
                continue
            if new:
                if self.par['cluster']:  # only append is solute is part of it
                    labels = [a, b, c, d, e]
                    testing = []
                    for ll in labels:
                        if ll != None:
                            testing.append(inst[ll]) 
                    if len([tt for tt in testing if tt in self.par['solute']]) > 0:
                        self.reactions[name].append(inst)
                else:
                    self.reactions[name].append(inst)
        return 0

def main():
    '''
    Find reaction patterns
    '''

    if __name__ == '__main__':
        main()
