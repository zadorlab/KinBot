from types import NoneType
import numpy as np
import logging
from ase.db import connect
from kinbot.stationary_pt import StationaryPoint
from kinbot import constants
from kinbot.vrc_tst_surfaces import VRC_TST_Surface

logger = logging.getLogger('KinBot')

def get_ra(reac, parent_chemid):
    reactive_atoms = []
    #Recover the index of the reactive atoms from the reaction name.
    with open(reac[0] + '/summary_' + reac[0] + '.out', 'r') as summary:
        lines = summary.readlines()
    for line in lines[4:]:
        if line.startswith('SUCCESS') and ("vdW" in line or "hom_sci" in line):
            if 'vdW' not in line :
                success, ts_energy, reaction_name, *products = line.split()
                if reaction_name == reac[1]:
                    break
            elif 'vdW' in line : #Unpack differently when a vdW well is in line
                success, ts_energy, reaction_name, *products, vdW_energy, vdW_direction = line.split()
                if reaction_name + vdW_direction[3:] == reac[1]:
                    break

    if sorted(products) == sorted(reac[2]):
        if 'hom_sci' in reaction_name:
            for ra in reaction_name.split('_')[3:]:
                reactive_atoms.append([int(ra)-1])
        elif 'vdW' in line:
            vdW_well = reaction_name + vdW_direction[3:]
            db = connect(f"{parent_chemid}/kinbot.db")
            *_, last_row = db.select(name=f"{vdW_well}")
            sp_vdW_well = StationaryPoint.from_ase_atoms(last_row.toatoms())
            sp_vdW_well.characterize()
            frags, maps = sp_vdW_well.start_multi_molecular()
            min_dist = np.inf
            for ra1 in maps[0]:
                for ra2 in maps[1]:
                    if ra1 != ra2 and sp_vdW_well.dist[ra1, ra2] < min_dist:
                        min_dist = sp_vdW_well.dist[ra1, ra2]
                        closest = [[ra1], [ra2]]
            reactive_atoms = closest

        else:
            logger.warning("Identification of reactive atoms is not implemented yet for this type of reactions.")

    #return list of reactive atoms with indexes in parent
    #reactive_atoms = [[Reactive atoms for frag=reac[2][0]], [Reactive atoms for frag=reac[2][1]]]
    return reactive_atoms

def create_all_surf_for_dist(dist, fragments, par, reactive_atoms):
    surfaces = []
    pps_lengths = [[],[]]
    pps_lengths_indexes = []
    if dist <= max(par['pp_on_atom']):
        for findex, frag in enumerate(fragments):
            frag.set_ra(reactive_atoms[findex]) #detect equivalent atoms
    
    #Create 'next_to_atom' surfaces
    if dist <= max(par['pp_oriented']):
        for findex, frag in enumerate(fragments):
            for ra in frag.ra:
                #Creates the list of pp_length to use for each ra
                #Here, values are in Angstroms
                pps_lengths[findex].append(par['pp_length'][frag.atom[ra]])
            pps_lengths_indexes.append([0 for i in range(len(frag.ra))])
        
        next = 0
        stop = False
        #Fragment number
        fn = 0
        while not stop:

            while pps_lengths_indexes[fn][next] <= len(pps_lengths[fn][next])-1:
                #List of 1 pp_length per ra on frag1
                pp_length_f1 = []
                for ra_index_f1, length_index_f1 in enumerate(pps_lengths_indexes[0]):
                    pp_length_f1.append(pps_lengths[0][ra_index_f1][length_index_f1])
                #List of 1 pp_length per ra on frag2
                pp_length_f2 = []
                for ra_index_f2, length_index_f2 in enumerate(pps_lengths_indexes[1]):
                    pp_length_f2.append(pps_lengths[1][ra_index_f2][length_index_f2])

                #Create surface for a combination of pp_length for all RA
                surfaces.append(create_surface(dist, fragments, [pp_length_f1, pp_length_f2]))

                if pps_lengths_indexes[fn][next] < len(pps_lengths[fn][next])-1:
                    pps_lengths_indexes[fn][next] +=1
                else:
                    break
            if next + 1 < len(pps_lengths[fn]):
                next += 1
            
            #Find next indice
            if pps_lengths_indexes[fn][next] == len(pps_lengths[fn][next])-1:
                if next + 1 < len(pps_lengths[fn]) and \
                pps_lengths_indexes[fn][next+1] != len(pps_lengths[fn][next+1])-1:
                    next +=1
                elif next + 2 < len(pps_lengths[fn]) and \
                len(pps_lengths[fn][next+1])-1 == 0:
                    next +=2
                else:
                    fn += 1
                    next = 0
                    if pps_lengths_indexes[fn][next] == len(pps_lengths[fn][next])-1:
                        if next + 1 < len(pps_lengths[fn]) and \
                        pps_lengths_indexes[fn][next+1] != len(pps_lengths[fn][next+1])-1:
                            next +=1
                        #Add more elif, or a while loop if multiple ra on same frag only have one pp_length
                        elif next + 2 < len(pps_lengths[fn]) and \
                        len(pps_lengths[fn][next+1])-1 == 0:
                            next +=2
            #Exit loop when finished
            for fidx in range(len(pps_lengths)):
                for ra_idx, ra_lidx in enumerate(pps_lengths_indexes[fidx]):
                    if ra_lidx != len(pps_lengths[fidx][ra_idx])-1:
                        stop = False
                        break
                    else:
                        stop = True
                if stop == False:
                    break
                    
            #Increase next indice
            pps_lengths_indexes[fn][next] +=1

            #reset previous indexes after increase
            for rst in range(next):
                pps_lengths_indexes[fn][rst] = 0
            while fn != 0:
                fn -= 1
                for rst in range(len(pps_lengths_indexes[fn])):
                    pps_lengths_indexes[fn][rst] = 0

            next = 0
            fn = 0

    #Create 'on_atom' surfaces
    if dist <= max(par['pp_on_atom']) and\
       dist >= min(par['pp_on_atom']):
        surfaces.append(create_surface(dist, fragments=fragments))

    if dist >= par['pp_on_COM']:
        surfaces.append(create_surface(dist))

    return surfaces

def create_surface(dist, fragments=None, pps_dists=None):
    """Collect the coord of all pivot points and
    return a list containing VRC_TST surfaces."""
    pps_coords = [[],[]]
    info = [None, None, f'Distance: {dist} Angstrom']
    if fragments != None:
        for findex, frag in enumerate(fragments):
            info[findex] = f'Fragment {findex}: '
            if frag.natom == 1:
                pps_coords[findex].append(frag.get_pp_on_com())
                info[findex] += 'COM'
                continue
            elif pps_dists is None:
                info[findex] += 'on atom'
                for ra in frag.ra:
                    pps_coords[findex].append(frag.get_pp_on_atom(ra))
                    info[findex] += f' {ra}'
            else:
                info[findex] += 'next to atom'
                for ra_index, ra in enumerate(frag.ra):
                    coord = frag.get_pp_next_to_ra(ra, pps_dists[findex][ra_index])
                    info[findex] += f' {ra} ({pps_dists[findex][ra_index]/ constants.BOHRtoANGSTROM} bohr)'
                    if isinstance(coord[0], list):
                        pps_coords[findex].extend(coord)
                    else:
                        pps_coords[findex].append(coord)
    else:
        pps_coords = [[[0.0, 0.0, 0.0]],[[0.0, 0.0, 0.0]]]
        info = ['Fragment 0: COM','Fragment 1: COM', f'Distance: {dist} Angstrom']

    #Create dist matrix
    dist_dim = (len(pps_coords[1]),len(pps_coords[0]))
    dist_matrix = np.full(dist_dim, dist)


    return VRC_TST_Surface(pps_coords, dist_matrix, info)


def reset_reactive_atoms(reactive_atoms, maps):
    '''
    Function that reorder the list of the reactive atoms to be sure they are in the same
    order as the maps of the fragments.
    '''
    new_list = [[None] * len(reactive_atoms[0]), [None] * len(reactive_atoms[1])]
    maps_array = np.full((len(maps), max([len(mmap) for mmap in maps])), None)

    for index, mmap in enumerate(maps):
        maps_array[index,0:len(mmap)]= mmap

    for ra in reactive_atoms:
        new_list[np.array(list(zip(*np.where(maps_array == ra))))[0,0]] = ra
    return new_list

def set_parent_map(parent, reactive_atoms, fragments):
    '''
    Funtion that map the ids of the long-range fragments' atoms
    from the ids of the parent.
    '''
    for ra_f1 in reactive_atoms[0]:
        for ra_f2 in reactive_atoms[1]:
            #Cut the bonds in parent between frag1 and frag2
            parent.bond[ra_f1,ra_f2] = 0
            parent.bond[ra_f2,ra_f1] = 0

    #Find the short-range fragment in the reactant
    #rearange the sr_frags and their map to have same order as lr_frags, because chemids are sorted when creating barrierless
    sr_fragments, maps = zip(*sorted(zip(*parent.start_multi_molecular()), key = lambda fm : fm[0].chemid ))

    if len(sr_fragments) != 2:
        #This error can happen in some cases where
        #two equivalent atoms break their bond with the other fragment
        #during the fragmentation. Ex: cycle divided in two parts/ double homolytic scission
        raise NotImplementedError('Cutting bonds in parent from identified reactive atoms did not lead to two fragment.')

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
