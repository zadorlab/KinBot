from kinbot.reactions.reac_vrc_tst_scan import VrcTstScan
from kinbot import geometry
import numpy as np
import logging
import copy
from ase import Atoms
import ase

logger = logging.getLogger('KinBot')

class VrcTstScanFrozen(VrcTstScan):
    def __init__(self, species, qc, par, instance, instance_name):
        super().__init__(species, qc, par, instance, instance_name, scan_type="sample")
        """
        Frozen fragments are updated at each step of the scan.
        The internal geometry shouldn't change, but 
        """
        self.species = self.get_frozen_species()

    def get_frozen_species(self, level="L1", distance=None):
        """
        Function that returns a stationary_point object 
        with a given distance between the fragments.
        The fragments internal coordinates are optimized for the desired level.
        The orientation between the fragments is kept the same as in self.species.
        """
        #distance between the two fragments
        if distance == None or not isinstance(distance, float):
            distance = self.scan_list[0]

        interfragments_param = self.get_interfragments_param(level=level)

        #Recover the geometry of the frozen fragments for desired level
        geometries = [geom for geom in self.long_range["geom"][level].values()]

        #Set inter-fragment distance
        frag_number = [[frag_index for frag_index, fmap in enumerate(self.long_range["maps"]) if atom_index in fmap][0] for atom_index in self.instance]
        point_A = geometries[frag_number[0]][np.where(self.long_range["maps"][frag_number[0]] == self.instance[0])[0]][0]
        point_B = geometries[frag_number[1]][np.where(self.long_range["maps"][frag_number[1]] == self.instance[1])[0]][0]

        current_distance = point_B - point_A + 0.0000001 #avoids superposition of identical fragments
        shift = geometry.unit_vector(current_distance)*distance - current_distance
        geometries[frag_number[1]] += shift

        #modify relative orientation of frozen fragments (dihedrals and angles)
        for modif in reversed(interfragments_param):
            match len(modif):
                case 5:
                    frag_number = [[frag_index for frag_index, fmap in enumerate(self.long_range["maps"]) if atom_index in fmap][0] for atom_index in modif[:-1]]

                    #Center the fragment coordinates on the atom to rotate around
                    
                    point_A = geometries[frag_number[0]][np.where(self.long_range["maps"][frag_number[0]] == modif[0])[0]][0]
                    point_B = geometries[frag_number[1]][np.where(self.long_range["maps"][frag_number[1]] == modif[1])[0]][0]
                    point_C = geometries[frag_number[2]][np.where(self.long_range["maps"][frag_number[2]] == modif[2])[0]][0]
                    point_D = geometries[frag_number[3]][np.where(self.long_range["maps"][frag_number[3]] == modif[3])[0]][0]
                    coord = geometries[frag_number[1]] - point_B

                    axis = geometry.unit_vector(point_C - point_B)

                    current_dihedral = geometry.calc_dihedral(point_A, point_B, point_C, point_D)[0]
                    rotation_angle = np.radians(modif[-1] - current_dihedral)

                    geometries[frag_number[0]] = np.dot(coord, geometry.rotation_matrix(axis, rotation_angle)) + point_B
                case 4:
                    frag_number = [[frag_index for frag_index, fmap in enumerate(self.long_range["maps"]) if atom_index in fmap][0] for atom_index in modif[:-1]]

                    #Center the fragment coordinates on the atom to rotate around

                    point_A = geometries[frag_number[0]][np.where(self.long_range["maps"][frag_number[0]] == modif[0])[0]][0]
                    point_B = geometries[frag_number[1]][np.where(self.long_range["maps"][frag_number[1]] == modif[1])[0]][0]
                    point_C = geometries[frag_number[2]][np.where(self.long_range["maps"][frag_number[2]] == modif[2])[0]][0]
                    coord = geometries[frag_number[1]] - point_B

                    axis = geometry.unit_vector(np.cross(point_A - point_B, point_C - point_B))

                    current_angle = np.degrees(geometry.calc_angle(point_A, point_B, point_C))
                    rotation_angle = np.radians(modif[-1] - current_angle)

                    geometries[frag_number[0]] = np.dot(coord, geometry.rotation_matrix(axis, rotation_angle)) + point_B

        #Set inter-fragment distance
        frag_number = [[frag_index for frag_index, fmap in enumerate(self.long_range["maps"]) if atom_index in fmap][0] for atom_index in self.instance]
        point_A = geometries[frag_number[0]][np.where(self.long_range["maps"][frag_number[0]] == self.instance[0])[0]][0]
        point_B = geometries[frag_number[1]][np.where(self.long_range["maps"][frag_number[1]] == self.instance[1])[0]][0]

        current_distance = point_B - point_A
        shift = geometry.unit_vector(current_distance)*distance - current_distance
        geometries[frag_number[1]] += shift

        #recreate bimolecular system with frag geometries in same order as parent
        bimolecular_geom = np.empty((self.species.natom,3))
        for line_number in range(self.species.natom):
            for frag_number, frag_map in enumerate(self.long_range["maps"]):
                if line_number in frag_map:
                    geom = geometries[frag_number]
                    bimolecular_geom[line_number] = geom[np.where(frag_map == line_number)[0][0]]
        tmp_species = copy.deepcopy(self.species)
        tmp_species.geom = bimolecular_geom
        tmp_species.characterize()

        # self.print_dif(tmp_species, interfragments_param, distance)
        return tmp_species

    def print_dif(self, geometries, interfragments_param, distance): #debugging funtion
        #recreate bimolecular system with frag geometries in same order as parent
        
        bimolecular_geom = np.empty((self.species.natom,3))
        for line_number in range(self.species.natom):
            for frag_number, frag_map in enumerate(self.long_range["maps"]):
                if line_number in frag_map:
                    geom = geometries[frag_number]
                    bimolecular_geom[line_number] = geom[np.where(frag_map == line_number)[0][0]]

        myatom = Atoms(self.species.atom, positions=bimolecular_geom)

        changes = []
        for frozen_coords in self.long_range["coord"]["L1"].values():
            changes.extend(frozen_coords)

        changes.extend(interfragments_param)
        changes.append(copy.copy(self.instance))
        changes[-1].append(distance)

        for index, frozen_coord in enumerate(changes):
            target = frozen_coord[-1]
            match len(frozen_coord[:-1]):
                case 2:
                    diff = myatom.get_distance(frozen_coord[0], frozen_coord[1]) - target
                    coord_type = 'bond'
                case 3:
                    diff = myatom.get_angle(frozen_coord[0], frozen_coord[1], frozen_coord[2]) - target
                    coord_type = 'angle'
                case 4:
                    diff = myatom.get_dihedral(frozen_coord[0], frozen_coord[1], frozen_coord[2], frozen_coord[3]) - target
                    coord_type = 'dihedral'
            message = "Diff for {coord_type:9s} {indexes}".format(coord_type=coord_type, indexes=frozen_coord[:-1])
            print(f"{message:35s}: {diff:11.6f}")
        print("=======================================")
        ase.io.write(filename="print_dif.xyz", images=myatom, format="xyz")

                
    def get_interfragments_param(self, level="L1"):
        """
        Function that selects the neighbourgs of the closest atom in each fragment,
        and returns the values of the interfragments angles and dihedrals.
        The level is only used as small orientation correction for planar fragments.
        """
        save_name = copy.copy(self.species.name)
        closest_to_RA = [[],[]] #each list contains the 2 (or less) closest atoms to the reactive atom in their respective fragment
        fragments, maps = self.species.start_multi_molecular()
        if len(fragments) == 1:
            self.species.name = save_name
            self.species.bond[self.instance[0], self.instance[1]] = 0
            self.species.bond[self.instance[1], self.instance[0]] = 0
            self.species.bonds[0][self.instance[0], self.instance[1]] = 0
            self.species.bonds[0][self.instance[1], self.instance[0]] = 0
            fragments, maps = self.species.start_multi_molecular()
        for frag_number, ra in enumerate(self.instance):
            for neighbourg_index, bond in zip(maps[frag_number], fragments[frag_number].bond[np.where(maps[frag_number] == ra)[0][0]]):
                if bond != 0:
                    closest_to_RA[frag_number].append(neighbourg_index)
                    if len(closest_to_RA[frag_number]) == 2:
                        break
            #If no neighbourgs atoms were found, the fragment has to be mono-atomic
            if len(closest_to_RA[frag_number]) == 1 and fragments[frag_number].natom > 2:
                for neighbourg_index, bond in zip(maps[frag_number], fragments[frag_number].bond[np.where(maps[frag_number] == closest_to_RA[frag_number][0])[0][0]]):
                    if bond != 0 and neighbourg_index != self.instance[frag_number]:
                        closest_to_RA[frag_number].append(neighbourg_index)
                        if len(closest_to_RA[frag_number]) == 2:
                            break
        changes = []
        #Calculate the angles values:
        modif_angle = self.check_deformation_into_planar(fragments[0], self.long_range["geom"][level][f"{0}"], maps[0], self.instance[0])
        if len(closest_to_RA[0]) != 0:
            point_A = self.species.geom[closest_to_RA[0][0]]
            point_B = self.species.geom[self.instance[0]]
            point_C = self.species.geom[self.instance[1]]
            changes.append([closest_to_RA[0][0],
                                   self.instance[0],
                                   self.instance[1],
                                   (np.degrees(geometry.calc_angle(point_A, point_B, point_C)))- modif_angle])
        modif_angle = self.check_deformation_into_planar(fragments[1], self.long_range["geom"][level][f"{1}"], maps[1], self.instance[1])
        if len(closest_to_RA[1]) != 0:
            point_A = self.species.geom[closest_to_RA[1][0]]
            point_B = self.species.geom[self.instance[1]]
            point_C = self.species.geom[self.instance[0]]
            changes.append([closest_to_RA[1][0],
                                   self.instance[1],
                                   self.instance[0],
                                   (np.degrees(geometry.calc_angle(point_A, point_B, point_C)))- modif_angle])
        
        #Calculate the dihedrals values:
        if len(closest_to_RA[0]) != 0 and len(closest_to_RA[1]) != 0:
            point_A = self.species.geom[closest_to_RA[0][0]]
            point_B = self.species.geom[self.instance[0]]
            point_C = self.species.geom[self.instance[1]]
            point_D = self.species.geom[closest_to_RA[1][0]]
            changes.append([closest_to_RA[0][0],
                            self.instance[0],
                            self.instance[1],
                            closest_to_RA[1][0],
                            geometry.calc_dihedral(point_A, point_B, point_C, point_D)[0]])
        modif_angle = self.check_deformation_into_planar(fragments[0], self.long_range["geom"][level][f"{0}"], maps[0], self.instance[0])
        if len(closest_to_RA[0]) == 2:
            point_A = self.species.geom[closest_to_RA[0][1]]
            point_B = self.species.geom[closest_to_RA[0][0]]
            point_C = self.species.geom[self.instance[0]]
            point_D = self.species.geom[self.instance[1]]
            changes.append([closest_to_RA[0][1],
                            closest_to_RA[0][0],
                            self.instance[0],
                            self.instance[1],
                            self.calc_corrected_dihedral(point_A, point_B, point_C, point_D, modif_angle)[0]])
        modif_angle = self.check_deformation_into_planar(fragments[1], self.long_range["geom"][level][f"{1}"], maps[1], self.instance[1])
        if len(closest_to_RA[1]) == 2:
            point_A = self.species.geom[closest_to_RA[1][1]]
            point_B = self.species.geom[closest_to_RA[1][0]]
            point_C = self.species.geom[self.instance[1]]
            point_D = self.species.geom[self.instance[0]]
            changes.append([closest_to_RA[1][1],
                            closest_to_RA[1][0],
                            self.instance[1],
                            self.instance[0],
                            self.calc_corrected_dihedral(point_A, point_B, point_C, point_D, modif_angle)[0]])
        
        return changes
    
    def calc_corrected_dihedral(self, a, b, c, d, alpha):
        """
        Calculate the A - B - C - D dihedral angle in radians.
        For collinear or close to collinear structures return a warning.
        alpha: Deformation angle (from out-of-plane to planar, in degree) 
        between long-range and short range geometry (ex: ch3)

        Returns the value in degrees
        np.arctan2 returns a value on the [-Pi,+Pi] interval
        """
        collinear_cutoff = 175./180.
        collinear = 0
        if (abs(geometry.calc_angle(a, b, c)) > np.pi * collinear_cutoff or
                abs(geometry.calc_angle(b, c, d)) > np.pi * collinear_cutoff):
            collinear = 1

        b0 = a - b  # reversed on purpose
        b1 = c - b
        b2 = d - c

        # normalize b1 so that it does not influence magnitude of vector
        b1 /= np.linalg.norm(b1)

        # v = projection of b0 onto plane perpendicular to b1
        #   = b0 minus component that aligns with b1
        # w = projection of b2 onto plane perpendicular to b1
        #   = b2 minus component that aligns with b1
        v = b0 - np.dot(b0, b1)*b1
        w = b2 - np.dot(b2, b1)*b1

        # angle between v and w in a plane is the torsion angle
        # v and w may not be normalized but that's fine since tan is y/x
        x = np.dot(v, w)
        y = np.dot(np.cross(b1, v), (w+ np.sin(np.deg2rad(alpha))*np.linalg.norm(c - b)*geometry.unit_vector(b0))/np.cos(np.deg2rad(alpha))) 
        return np.degrees(np.arctan2(y, x)), collinear

    def check_deformation_into_planar(self, stationary_point, geom, stp_map, reactive_atom):
        """Function that checks the planarity of individualy optimised fragments around the reactive atom,
        and returns the out-of-plane angle of the short-range optimized fragment"""
        neighbourgs = 0
        neighbourgs_index = []
        for index, bond in enumerate(stationary_point.bond[np.where(stp_map == reactive_atom)[0][0]]):
            if bond > 0:
                neighbourgs += 1
                neighbourgs_index.append(index)
        neighbourgs_index.append(np.where(stp_map == reactive_atom)[0][0])

        if neighbourgs == 3:
            point_A = stationary_point.geom[neighbourgs_index[0]]
            point_B = stationary_point.geom[neighbourgs_index[1]]
            point_C = stationary_point.geom[neighbourgs_index[2]]
            point_D = stationary_point.geom[neighbourgs_index[3]]
            original_angle, colinear = np.degrees(geometry.calc_out_of_plane_angle(point_A, point_B, point_C, point_D))

            point_A = geom[neighbourgs_index[0]]
            point_B = geom[neighbourgs_index[1]]
            point_C = geom[neighbourgs_index[2]]
            point_D = geom[neighbourgs_index[3]]
            new_angle, colinear = np.degrees(geometry.calc_out_of_plane_angle(point_A, point_B, point_C, point_D))
            angle = original_angle - new_angle
            return angle
        else:
            return 0.

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []

        val = self.scan_list[step]
        self.set_bond(0, 1, val, change)

        self.clean_constraints(change, fix)
        for frozen_coords in self.long_range["coord"]["L1"].values():
            for c in frozen_coords:
                c_new = [ci + 1 for ci in c[:-1]]
                c_new.append(c[-1])
                change.append(c_new)

        return step, fix, change, release
    