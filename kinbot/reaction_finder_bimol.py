import numpy as np
import sys
import logging

from kinbot.reactions.reac_abstraction import Abstraction

logger = logging.getLogger('KinBot')


class ReactionFinderBimol:
    """
    Class to find all reactions for bimolecular starting species.
    """
    
    def __init__(self, species, par, qc):
        self.species = species
        self.qc = qc
        self.par = par
        self.families = par['bimol_families']
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
        """
        List all reaction types available, and find the key atoms for them.
        """

        reaction_names = {'abstraction': self.search_abstraction,
                          }

        if 'combinatorial' in self.families:
            reaction_names['combinatorial'] = self.search_combinatorial        
 
        atomA = self.species.fragA.atom
        atomB = self.species.fragB.atom
        natomA = self.species.fragA.natom
        natomB = self.species.fragB.natom
        uniqA = self.species.fragA.atom_uniq
        uniqB = self.species.fragB.atom_uniq
        
        for i, bondA in enumerate(self.species.fragA.bonds):  # for all resonance structures
            radA = self.species.fragA.rads[i]
            for j, bondB in enumerate(self.species.fragB.bonds):  # for all resonance structures
                radB = self.species.fragB.rads[j]
                for rn in reaction_names:
                    if rn in self.families or 'all' in self.families:
                        if not rn in self.skip_families:
                            reaction_names[rn](natomA, atomA, bondA, radA, uniqA, natomB, atomB, bondB, radB, uniqB, natomA)  # last natomA is a shift
                            reaction_names[rn](natomB, atomB, bondB, radB, uniqB, natomA, atomA, bondA, radA, uniqA, natomA)  # swapping roles for all families

        for name in self.reactions:
            self.reaction_matrix(self.reactions[name], name) 
        
        for index in range(len(self.species.reac_name)-1):
            if self.species.reac_name[index] in self.species.reac_name[index + 1:]:
                logger.error('Found reaction name "{}" more than once'
                               .format(self.species.reac_name[index]))
                logger.error('Exiting')
                sys.exit()

        logger.info('\tFound the following bimolecular reactions:')
        for rxn in self.species.reac_name:
            logger.info('\t\t{}'.format(rxn))
        
        return 0  
   

    def search_abstraction(self, natomA, atomA, bondA, radA, uniqA, natomB, atomB, bondB, radB, uniqB, shift):
        """ 
        This is a general atom abstraction family.

        A-a + b-B --> A-a-b + B
          0   1 2

        b is  
            - leaf atom 
            - with a single bond 
            - and no free valence, 
            - e.g., -H or -Cl, but not =O
        a is an abstractor site
            - radical location
            - ?
        """
        
        name = 'abstraction'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        instances = []
        
        for radA_site in list(set(np.nonzero(radA)[0]) & set(uniqA)):  # unique radical centers
            for bb in set(uniqB):
                if atomB[bb] in ['H']:  # TODO expand to other atoms
                    for ee, bond in enumerate(bondB[bb]):
                        if bond == 1:
                            instances.append([radA_site, bb + shift, ee + shift])
                            break

        for instance in instances: 
            rxns.append(instance)
        
        self.new_reaction(rxns, name, a=0, b=1, c=2)
        
        return 0

    def reaction_matrix(self, reac_list, reac_id):
        """ 
        Create arrays to store all reactions for species.
        input: 
        reac_list: atom motifs from individual searches
        reac_id: reaction name (e.g., HO2_Elimination_from_PeroxyRadical) from individual searc functions
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
        """
        
        self.species.reac_type += [reac_id for i in range(len(reac_list))]
        self.species.reac_inst += reac_list
        self.species.reac_step += [0 for i in range(len(reac_list))]
        self.species.reac_scan_energy += [[] for i in range(len(reac_list))]
        self.species.reac_ts_done += [0 for i in range(len(reac_list))] 
        self.species.reac_ts_geom += [0 for i in range(len(reac_list))]
        self.species.reac_ts_freq += [0 for i in range(len(reac_list))]
        
        for i in range(len(reac_list)):
            if reac_id == 'abstraction':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(Abstraction(self.species, self.qc, self.par, reac_list[i], name))
            else:
                self.species.reac_name.append(0)
        return 0

    def new_reaction(self, rxns, name, a=None, b=None, c=None, d=None, e=None, length=None, full=False, cross=False):
        """
        Returns 1 if new, and 0 if not new
        Checks a variable number of identical elements
        Also can check full equivalency (full=True), same lenght (length=True), and
        equivalency between elements that are interchangeable (cross=True)
        """

        for inst in rxns:
            new = True
            for instance in self.reactions[name]:
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
                    if any([inst[i] != instance[i] for i, _ in enumerate(inst)]):
                        continue
                new = False
                continue
            if new:
                self.reactions[name].append(inst)

        return 0

def main():
    """
    Find reaction patterns
    """

    if __name__ == "__main__":
        main()
