import numpy as np
import sys
import logging

from kinbot.reactions.reac_abstraction import Abstraction
from kinbot.species_routing import routing_name
from kinbot.stereochemistry import (refine_equivalence_group, stereotopic_classes,
                                     stereotopic_site_signature)

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
        # keys: names of the families
        # values: list of instances
        # this dict is used to keep track of the unique reactions found,
        # and to verify whether a new reaction is indeed unique 
        self.reactions = {}
        self.site_metadata = {}
        self._stereotopic_cache = {}
        self._site_signature_cache = {}

    def _site_classes(self, species, atoms):
        key = (id(species), tuple(sorted(int(atom) for atom in atoms)))
        if key not in self._stereotopic_cache:
            self._stereotopic_cache[key] = stereotopic_classes(
                species, key[1])
        return self._stereotopic_cache[key]

    def _site_signature(self, atom):
        if atom not in self._site_signature_cache:
            self._site_signature_cache[atom] = stereotopic_site_signature(
                self.species, atom)
        return self._site_signature_cache[atom]

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
                            reaction_names[rn](
                                natomA, atomA, bondA, radA, uniqA,
                                natomB, atomB, bondB, radB, uniqB, 0, natomA)
                            reaction_names[rn](
                                natomB, atomB, bondB, radB, uniqB,
                                natomA, atomA, bondA, radA, uniqA, natomA, 0)

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
   

    def search_abstraction(self, natomA, atomA, bondA, radA, uniqA,
                           natomB, atomB, bondB, radB, uniqB,
                           shiftA, shiftB):
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
            - radical location or either atom of neutral O2

        shiftA and shiftB map fragment-local indices to the combined order.
        """
        
        name = 'abstraction'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        instances = []
        radical_sites = [int(atom + shiftA)
                         for atom in np.nonzero(radA)[0]]
        fragment = self.species.fragA if shiftA == 0 else self.species.fragB
        oxygen = (natomA == 2 and all(atom == 'O' for atom in atomA)
                  and bondA[0][1] > 0 and fragment.charge == 0)
        if oxygen:
            # Singlet and triplet O2 have an O=O graph without radical sites.
            # Either oxygen can accept H; keep the graph and spin unchanged.
            radical_sites = [shiftA, shiftA + 1]
        hydrogen_donors = {}
        for atom in range(natomB):
            if atomB[atom] != 'H':  # TODO expand to other atoms
                continue
            donor = next((index for index, order in enumerate(bondB[atom])
                          if order == 1), None)
            if donor is not None:
                hydrogen_donors[int(atom + shiftB)] = int(donor + shiftB)

        # Classify on the disconnected combined graph so a fixed chiral
        # reacting partner turns an otherwise enantiotopic pair into two
        # diastereotopic reaction paths.
        if oxygen:
            # O2's atom_eqv groups can leave its equivalent O atoms separate.
            radical_classes = refine_equivalence_group(self.species, radical_sites)
        else:
            radical_classes = self._site_classes(self.species, radical_sites)
        hydrogen_classes = self._site_classes(
            self.species, hydrogen_donors)
        for radical_class in radical_classes:
            for hydrogen_class in hydrogen_classes:
                radical = radical_class['representative']
                hydrogen = hydrogen_class['representative']
                donor = hydrogen_donors[hydrogen]
                instance = [int(radical), int(hydrogen), int(donor)]
                instances.append(instance)
                stereopath_id = (
                    'abstraction:'
                    f'{self._site_signature(radical)}:'
                    f'{self._site_signature(hydrogen)}')
                self.site_metadata.setdefault(tuple(instance), {
                    'stereopath_id': stereopath_id,
                    'site_degeneracy': (
                        len(radical_class['members'])
                        * len(hydrogen_class['members'])),
                    'site_relations': {
                        'abstractor': radical_class['relation'],
                        'transferred_atom': hydrogen_class['relation'],
                    },
                    'site_members': {
                        'abstractor': [int(atom) for atom in
                                       radical_class['members']],
                        'transferred_atom': [int(atom) for atom in
                                             hydrogen_class['members']],
                    },
                    'allow_reflection': not np.any(self.species.chiral),
                })

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
                name = routing_name(self.species) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1)
                self.species.reac_name.append(name)
                reaction = Abstraction(
                    self.species, self.qc, self.par, reac_list[i], name)
                for key, value in self.site_metadata.get(
                        tuple(reac_list[i]), {}).items():
                    setattr(reaction, key, value)
                self.species.reac_obj.append(reaction)
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
