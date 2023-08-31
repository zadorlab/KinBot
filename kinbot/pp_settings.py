import numpy as np
import logging

logger = logging.getLogger('KinBot')

def get_ra(reac):
    reactive_atoms = []
    #Recover the index of the reactive atoms from the reaction name.
    with open(reac[0] + '/summary_' + reac[0] + '.out', 'r') as summary:
        for line in summary.readlines()[4:]:
            if line.startswith('SUCCESS') and 'hom_sci' in line: #TO change to accept different type of barrierless reactions
                pieces = line.split()

                if sorted(pieces[3:]) == sorted(reac[2]) :
                    for ra in pieces[2].split('_')[3:]:
                        reactive_atoms.append(int(ra)-1)
                    break
                else:
                    logger.warning("Identification of reactive atoms is not implemented yet for this type of reactions.")
                    return "Not implemented"
                    
    return reactive_atoms

def reset_reactive_atoms(reactive_atoms, maps):
    '''
    Function that reorder the list of the reactive atoms to be sure they are in the same
    order as the maps of the fragments.
    '''
    new_list = [None] * len(reactive_atoms)
    maps_array = np.full((len(maps), max([len(mmap) for mmap in maps])), None)

    for index, mmap in enumerate(maps):
        maps_array[index,0:len(mmap)]= mmap

    for ra in reactive_atoms:
        new_list[np.array(list(zip(*np.where(maps_array == ra))))[0,0]] = ra
    return new_list

def set_order(parent, reactive_atoms, fragments):
    '''
    Funtion that map the ids of the long-range fragments' atoms
    from the ids of the parent.
    '''

    #Cut the bond in parent
    parent.bond[reactive_atoms[0],reactive_atoms[1]] = 0
    parent.bond[reactive_atoms[1],reactive_atoms[0]] = 0

    #Find the short-range fragment in the reactant
    #rearange the sr_frags and their map to have same order as lr_frags, because chemids are sorted when creating barrierless
    sr_fragments, maps = zip(*sorted(zip(*parent.start_multi_molecular()), key = lambda fm : fm[0].chemid ))

    reactive_atoms = reset_reactive_atoms(reactive_atoms, maps)

    #Map the sr fragments into the individually optimized lr fragments
    for sr_fragment, lr_fragment, sr_map in zip(sr_fragments, fragments, maps):
        lr_map = []
        #sr_map: list of indexes in the parent. The position corresponds to the index in the short-range fragment.
        for lr_atom in range(lr_fragment.natom):
            #Find atoms with matching atomid in lr_fragment
            sr_match = np.where(np.array(sr_fragment.atomid) == np.array(lr_fragment.atomid[lr_atom]))[0] #list of matching indexes
            if sr_match.size == 1:
                lr_map.append(sr_map[sr_match[0]])
            else:
            #add something if array empty
                for index in sr_match:
                    if sr_map[index] in lr_map:
                        continue
                    else:
                        lr_map.append(sr_map[index])
        lr_fragment.set_map(lr_map)

def set_surface(self):
    pp_dict = {}
    for index, frag in enumerate(self.lr_fragments):
        pp_dict['{index}'] = "np.array({self.lr_fragments.pivot_points})"
    distances = "distances = np.array({self.pp_dist})"

    self.divid_surf_block = self.divid_surf_block +\
                            "Surface(\
                                    {pp_dict},\n\
                                    {distances}),"
