from kinbot.fragments import Fragment
from ase.db import connect
from kinbot.stationary_pt import StationaryPoint

class VRC_TST_surfaces()
'''
Class that generates a list of pivot points and a list of pivot points
and their associated distances depending on the intermolecular distance.
'''

    def __init__(self, par, fragments, distances=None, surfaces=None, nfrag=None, reactive_atoms=None)
        self.distances = []
        self.distances.append(i for i in arange(par['vrc_tst_dist_start'],
                                                par['vrc_tst_dist_stop'],
                                                par['vrc_tst_dist_step']))
        self.lr_fragments = fragments #Array of Fragment objects
        self.par = par
        self.surfaces = []
        self.nfrag = len(fragments)
        self.reactive_atoms = [] #Contains a list of integers. These are the order number of the atoms involved in the reaction.
        self.setup_fragments()
        self.set_parent()
        self.set_reactive_atoms()
        self.set_order()

    def setup_fragments(self):
        for this_frag in self.lr_fragments:
            this_frag.characterize()


    #Create a stationary point object for the parent
    def set_parent(self):
        parent_chemid = self.fragment[0].parent_chemid
        if par['high_level']:
            #Will read info from L2 structure
            basename = '{parent_chemid}_well_high'
        else
            #Will read info from L1 structure
            basename = '{parent_chemid}_well'

        db = connect('{parent_chemid}/kinbot.db')
        for row in db.select(name='{basename}'):
            tmp = row.toatoms() #This is an ase.atoms object
            self.parent = StationaryPoint.from_ase_atoms(tmp)
            self.parent.characterize()

    def set_reactive_atoms(self, fragments)
        #Recover the index of the reactive atoms from the reaction name.
        summary = open(fragments[0].parent_chemid + '/summary_' + fragments[0].parent_chemid + '.out', 'r').readlines()
        #Check if the line corresponds to the reaction before taking the indices
        for line in summary:
            if line.startswith('SUCCESS') and 'hom_sci' in line:
                pieces = line.split()
                current_parent = pieces[2].split('_')[0]

                corresponds = 0
                #If this reaction has a number of product different than the number of fragments,
                #then skip this iteration.
                if len(pieces[3:]) == len(fragments):
                    pass
                else
                    continue

                for current_product, this_frag in zip(pieces[3:], fragments):
                    if current_parent == this_frag.parent_chemid and current_product == this_frag.chemid:
                        corresponds = 1
                    else
                        corresponds = 0
                        break
                #When the reaction is found, create the list of reactive atoms        
                if corresponds:
                    for i in pieces[2].split('_')[3:]:
                        self.reactive_atoms.append(int(i))
                    break

    def set_order(self):
    '''
    Funtion that map the ids of the long-range fragments' atoms
    from the ids of the parent.
    '''
        #Cut the bond in parent
        self.parent.bond[self.reactive_atoms[0]][self.reactive_atoms[1]] = 0
        self.parent.bond[self.reactive_atoms[1]][self.reactive_atoms[0]] = 0

        #Find the short-range fragment in the reactant
        self.sr_fragments, maps = self.parent.start_multi_molecular()
        
        #Map the sr fragments into the individually optimized lr fragments
        for sr_fragment, lr_fragment, sr_map in zip(self.sr_fragments, self.lr_fragments, maps):
            lr_map = []
            for lr_atom in range(lr_fragment.natom):
                #Find atoms with matching chemid in lr_fragment
                sr_match = np.where(sr_fragment.atomid == lr_fragment.atomid[lr_atom])[0] #list of matching indexes
                if sr_match.size == 1:
                    lr_map.append(sr_map[sr_match[0]])
                else
                #add something if array empty
                    for index in sr_match:
                        if sr_map[index] in lr_map:
                            pass
                        else
                            lr_map.append(sr_map[index])

            lr_fragment.set_map(lr_map)


    def get_surfaces(self):
        for frag, atom in zip(self.lr_fragments, self.reactive_atoms):
            frag.set_pivot_points(self.distances, atom) #Create the attribute frag.pivot_points: list of atom objects

        #Create a "surface" dictionary as expected by rotd_py
        #TODO Must return the coordinates of the pivot points from the fragment object,
        #and a pivot point distances matrix.
        self.set_surfaces()

    def set_surfaces(self):
        pass

