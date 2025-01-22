import logging
from typing import Any

import numpy as np
from numpy.typing import NDArray
from ase import Atoms

from kinbot.fragments import Fragment
from kinbot.vrc_tst_surfaces import VRC_TST_Surface

logger = logging.getLogger('KinBot')


def create_all_surf_for_dist(dist: float,
                             equiv_ra: list[list[list[int]]],
                             fragments: list[Fragment],
                             par: dict
                             ) -> tuple[list[list[int]],
                                        list[list[int]],
                                        list[VRC_TST_Surface]]:
    """Create all the surfaces, and detect which surface should
    be sampled, and with wich weight.

    Args:
        dist (float): distance along the reaction coordinate (Angstrom)
        equiv_ra (list[list[list[int]]]):
            [0] frag1, [1] frag2
            for each, contain multiple lists of indexes of atoms chemically
            equivalent.
        fragments (list[Fragment]): list of Fragments objects
        par (dict): User input parameters.

    Returns:
        tuple[list[list[int]],      :[[faces weights for a surf],[...]...]
              list[list[int]],      :[[indexes of faces to sample for a surf]]
              list[VRC_TST_Surface]]: Contains pivot points and a dist matrix
    """
    surfaces: list[VRC_TST_Surface] = []
    pps_l: list[list[Any]] = [[], []]
    pp_l_idxs = []
    fw: list[int]
    sf: list[int]
    surf: VRC_TST_Surface
    selected_faces: list[list[int]] = []
    faces_weights: list[list[int]] = []

    # Create 'oriented_pp' surfaces
    if dist <= max(par['pp_oriented']):
        for findex, frag in enumerate(fragments):
            for ra in frag.ra:
                # Creates the list of pp_length to use for each ra
                # Here, values are in Angstroms
                if frag.atom[ra] in par['pp_length']:
                    pps_l[findex].append(par['pp_length'][frag.atom[ra]])
                else:
                    pps_l[findex].append(par['pp_length']['X'])
            pp_l_idxs.append([0 for i in range(len(frag.ra))])

        if 'X' not in par['pp_length']:
            (surfs, fws, sfs) = \
                create_combinatorial_surfaces(
                    pp_l_idxs=pp_l_idxs,
                    pps_l=pps_l,
                    dist=dist,
                    equiv_ra=equiv_ra,
                    fragments=fragments)
            surfaces.extend(surfs)
            faces_weights.extend(fws)
            selected_faces.extend(sfs)
        else:
            for i in range(len(par['pp_length']['X'])):
                pp_l_f1 = np.array(pps_l[0])[:, i].tolist()
                pp_l_f2 = np.array(pps_l[1])[:, i].tolist()
                (fw, sf, surf) = create_surface(
                    dist=dist,
                    equiv_ra=equiv_ra,
                    fragments=fragments,
                    pps_dists=[pp_l_f1, pp_l_f2])
                surfaces.append(surf)
                faces_weights.append(fw)
                selected_faces.append(sf)
    # Create 'on_atom' surfaces
    if dist <= max(par['pp_on_atom']) and\
       dist >= min(par['pp_on_atom']):
        (fw, sf, surf) = create_surface(
            dist=dist,
            equiv_ra=equiv_ra,
            fragments=fragments)
        surfaces.append(surf)
        faces_weights.append(fw)
        selected_faces.append(sf)

    if dist >= par['pp_on_COM']:
        (fw, sf, surf) = create_surface(
            dist=dist,
            equiv_ra=equiv_ra)
        surfaces.append(surf)
        faces_weights.append(fw)
        selected_faces.append(sf)

    return (faces_weights,
            selected_faces,
            surfaces)


def create_combinatorial_surfaces(pp_l_idxs,
                                  pps_l,
                                  dist,
                                  equiv_ra,
                                  fragments):
        surfaces = []
        faces_weights = []
        selected_faces = []
        next = 0
        stop = False
        # Fragment number
        fn = 0
        while not stop:
            while pp_l_idxs[fn][next] <= len(pps_l[fn][next])-1:
                # List of 1 pp_length per ra on frag1
                pp_l_f1 = []
                for ra_idx_f1, l_idx_f1 in enumerate(pp_l_idxs[0]):
                    pp_l_f1.append(pps_l[0][ra_idx_f1][l_idx_f1])
                # List of 1 pp_length per ra on frag2
                pp_l_f2 = []
                for ra_idx_f2, l_idx_f2 in enumerate(pp_l_idxs[1]):
                    pp_l_f2.append(pps_l[1][ra_idx_f2][l_idx_f2])

                # Create surface for a combination of pp_length for all RA
                (fw, sf, surf) = create_surface(
                    dist=dist,
                    equiv_ra=equiv_ra,
                    fragments=fragments,
                    pps_dists=[pp_l_f1, pp_l_f2])
                surfaces.append(surf)
                faces_weights.append(fw)
                selected_faces.append(sf)

                if pp_l_idxs[fn][next] < len(pps_l[fn][next])-1:
                    pp_l_idxs[fn][next] += 1
                else:
                    break
            if next + 1 < len(pps_l[fn]):
                next += 1

            # Find next indice
            if pp_l_idxs[fn][next] == len(pps_l[fn][next])-1:
                if next + 1 < len(pps_l[fn]) and \
                   pp_l_idxs[fn][next+1] != len(pps_l[fn][next+1])-1:
                    next += 1
                elif next + 2 < len(pps_l[fn]) and\
                     len(pps_l[fn][next+1])-1 == 0:
                    next += 2
                else:
                    fn += 1
                    next = 0
                    if pp_l_idxs[fn][next] == len(pps_l[fn][next])-1:
                        if next + 1 < len(pps_l[fn]) and \
                           pp_l_idxs[fn][next+1] != len(pps_l[fn][next+1])-1:
                            next += 1
                        # Add more elif, or a while loop if multiple ra
                        # on same frag only have one pp_length
                        elif next + 2 < len(pps_l[fn]) and \
                        len(pps_l[fn][next+1])-1 == 0:
                            next += 2
            # Exit loop when finished
            for fidx in range(len(pps_l)):
                for ra_idx, ra_lidx in enumerate(pp_l_idxs[fidx]):
                    if ra_lidx != len(pps_l[fidx][ra_idx])-1:
                        stop = False
                        break
                    else:
                        stop = True
                if stop is False:
                    break

            # Increase next indice
            pp_l_idxs[fn][next] += 1

            # reset previous indexes after increase
            for rst in range(next):
                pp_l_idxs[fn][rst] = 0
            while fn != 0:
                fn -= 1
                for rst in range(len(pp_l_idxs[fn])):
                    pp_l_idxs[fn][rst] = 0

            next = 0
            fn = 0
        return (surfaces, faces_weights, selected_faces)


def create_surface(dist,
                   equiv_ra: list[list[list[int]]],
                   fragments: list[Fragment] = [],
                   pps_dists: list[list[float]] | None = None
                   ) -> tuple[list[int], list[int], VRC_TST_Surface]:
    """Collect the coord of all pivot points and
    return a list containing VRC_TST surfaces."""
    pps_coords = [[], []]
    all_pp_dists = [[], []]
    # Contains the index of the faces to sample for a surface.
    selected_faces: list[int] = []
    # Contains the weight that multiplies the flux for the face
    # at a given index.
    weights: list[list[int]] = [[], []]
    # Used to detect faces belonging to same reaction,
    # and avoid sampling the other
    reac_weights: list[list[int]] = [[], []]
    faces_weights: list[int] = []
    info: list[str] = ['',
                       '',
                       f'Distance: {dist} Angstrom']
    if len(fragments) != 0:
        for findex, frag in enumerate(fragments):
            info[findex] = f'Fragment {findex}: '
            if frag.natom == 1:
                pps_coords[findex].append(frag.get_pp_on_com())
                info[findex] += 'COM'
                weights[findex].append(1)
                reac_weights[findex].append(1)
                all_pp_dists[findex].append(0.0)
                continue
            elif pps_dists is None:
                info[findex] += 'on atom'
                for ra in frag.ra:
                    all_pp_dists[findex].append(0.0)
                    pps_coords[findex].append(frag.get_pp_on_atom(ra))
                    info[findex] += f' {ra}'
                    for skip_face, equiv in enumerate(equiv_ra[findex]):
                        if ra == equiv[0]:
                            weights[findex].append(len(equiv))
                            if not skip_face:
                                reac_weights[findex].append(1)
                            else:
                                reac_weights[findex].append(0)
                            break
                        elif ra in equiv:
                            weights[findex].append(0)
                            reac_weights[findex].append(0)
                            break
            else:
                info[findex] += 'pp oriented'
                for ra_index, ra in enumerate(frag.ra):
                    coord: list[list[float]]
                    pp_dist: list[float]
                    coord, pp_dist = \
                        frag.get_pp_next_to_ra(
                            index=ra,
                            dist_from_ra=pps_dists[findex][ra_index])
                    all_pp_dists[findex].extend(pp_dist)
                    for skip_face, equiv in enumerate(equiv_ra[findex]):
                        if ra == equiv[0]:
                            weights[findex].append(len(coord)*len(equiv))
                            if not skip_face:
                                reac_weights[findex].append(1)
                            else:
                                reac_weights[findex].append(0)
                            for i in range(len(coord)-1):
                                weights[findex].append(0)
                                reac_weights[findex].append(0)
                            break
                        elif ra in equiv:
                            weights[findex].append(0)
                            reac_weights[findex].append(0)
                            if len(coord) == 2:
                                weights[findex].append(0)
                                reac_weights[findex].append(0)
                            break

                    info[findex] += ' {} ({} bohr)'.format(
                        ra,
                        pps_dists[findex][ra_index])
                    pps_coords[findex].extend(coord)
        for pp2 in range(len(pps_coords[1])):
            for pp1 in range(len(pps_coords[0])):
                faces_weights.append(
                    weights[0][pp1] * reac_weights[0][pp1] * weights[1][pp2] * reac_weights[1][pp2]
                    )
                # Only sample the faces that have a non-zero weight
                if faces_weights[-1] != 0:
                    selected_faces.append(len(faces_weights)-1)
    else:
        pps_coords: list[list[list[float]]] = [[[0.0, 0.0, 0.0]],
                                               [[0.0, 0.0, 0.0]]]
        all_pp_dists: list[list[float]] = [[0.0], [0.0]]
        info = ['Fragment 0: COM',
                'Fragment 1: COM',
                f'Distance: {dist} Angstrom']
        faces_weights: list[int] = [1]
        selected_faces: list[int] = [0]

    # Create distance matrix
    dist_dim: tuple[int, int] = (len(pps_coords[0]), len(pps_coords[1]))
    dist_matrix: NDArray[Any] = np.zeros(dist_dim)
    # Adjust the distance matrix from the distances used for each pivot point
    for ppf1 in range(len(pps_coords[0])):
        for ppf2 in range(len(pps_coords[1])):
            dist_matrix[ppf1, ppf2] = dist-(all_pp_dists[0][ppf1]+all_pp_dists[1][ppf2])
    for fnum, frag in enumerate(fragments):
        pps = ''
        for pp in range(len(pps_coords[fnum])):
            pps += 'X'
        pos = frag.geom.tolist()
        pos.extend(pps_coords[fnum])
        atm = Atoms(symbols=f'{frag.atom}{pps}',
                    positions=pos)
        atm.write(f'rotdPy/S_{VRC_TST_Surface.__id__}_F_{frag.chemid}.xyz')

    return (faces_weights,
            selected_faces,
            VRC_TST_Surface(pp_coords=pps_coords,
                            dist_matrix=dist_matrix,
                            info=info)
            )


# def reset_reactive_atoms(reactive_atoms, maps):
#     '''
#     Function that reorders the list of the reactive atoms to make sure they are in the same
#     order as the maps of the fragments.
#     '''
#     new_list = [[None] * len(reactive_atoms[0]), [None] * len(reactive_atoms[1])]
#     maps_array = np.full((len(maps), max([len(mmap) for mmap in maps])), None)

#     for index, mmap in enumerate(maps):
#         maps_array[index,0:len(mmap)]= mmap

#     for ra in reactive_atoms:
#         new_list[np.array(list(zip(*np.where(maps_array == ra))))[0,0]] = ra
#     return new_list

# def set_parent_map(parent, reactive_atoms, fragments):
#     '''
#     Funtion that maps the ids of the long-range fragments' atoms
#     from the ids of the parent.
#     '''
#     for ra_f1 in reactive_atoms[0]:
#         for ra_f2 in reactive_atoms[1]:
#             #Cut the bonds in parent between frag1 and frag2
#             parent.bond[ra_f1, ra_f2] = 0
#             parent.bond[ra_f2, ra_f1] = 0

#     # Find the short-range fragment in the reactant
#     # Rearange the sr_frags and their map to have same order as lr_frags, because chemids are sorted when creating barrierless
#     sr_fragments, maps = zip(*sorted(zip(*parent.start_multi_molecular()), key = lambda fm : fm[0].chemid ))

#     if len(sr_fragments) != 2:
#         #This error can happen in some cases where
#         #two equivalent atoms break their bond with the other fragment
#         #during the fragmentation. Ex: cycle divided in two parts/ double homolytic scission
#         raise NotImplementedError('Cutting bonds in parent from identified reactive atoms did not lead to two fragment.')

#     #Map the sr fragments into the individually optimized lr fragments
#     for sr_fragment, lr_fragment, sr_map in zip(sr_fragments, fragments, maps):
#         lr_map = []
#         #sr_map: list of indexes in the parent. The position corresponds to the index in the short-range fragment.
#         for lr_atom in range(lr_fragment.natom):
#             #Find atoms with matching atomid in lr_fragment
#             sr_match = np.where(np.array(sr_fragment.atomid) == np.array(lr_fragment.atomid[lr_atom]))[0] #list of matching indexes
#             if sr_match.size == 1:
#                 lr_map.append(sr_map[sr_match[0]])
#             else:
#             #add something if array empty
#                 for index in sr_match:
#                     if sr_map[index] in lr_map:
#                         continue
#                     else:
#                         lr_map.append(sr_map[index])
#         lr_fragment.set_map(lr_map)

# def set_surface(self):
#     pp_dict = {}
#     for index, frag in enumerate(self.lr_fragments):
#         pp_dict['{index}'] = "np.array({self.lr_fragments.pivot_points})"
#     distances = "distances = np.array({self.pp_dist})"

#     self.divid_surf_block = self.divid_surf_block +\
#                             "Surface(\
#                                     {pp_dict},\n\
#                                     {distances}),"
