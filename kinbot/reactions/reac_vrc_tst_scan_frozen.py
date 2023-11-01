from kinbot.reactions.reac_vrc_tst_scan import VrcTstScan
from kinbot import geometry
from sella import Internals
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

        interfragments_param = self.get_interfragments_param()

        geometries = [geom for geom in self.long_range["geom"][level].values()]

        #modify dihedrals and angles
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

                    current_dihedral = geometry.calc_dihedral(point_A, point_B, point_C, point_D)[0]%360
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

                    current_angle = np.degrees(geometry.calc_angle(point_A, point_B, point_C))%360
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

        #self.print_dif(tmp_species, interfragments_param, distance)
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
                    diff = myatom.get_angle(frozen_coord[0], frozen_coord[1], frozen_coord[2])%360 - target
                    coord_type = 'angle'
                case 4:
                    diff = myatom.get_dihedral(frozen_coord[0], frozen_coord[1], frozen_coord[2], frozen_coord[3])%360 - target
                    coord_type = 'dihedral'
            message = "Diff for {coord_type:9s} {indexes}".format(coord_type=coord_type, indexes=frozen_coord[:-1])
            print(f"{message:35s}: {diff:11.6f}")
        print("=======================================")
        ase.io.write(filename="print_dif.xyz", images=myatom, format="xyz")

                
    def get_interfragments_param(self):
        """
        Function that selects the neighbourgs of the closest atom in each fragment,
        and returns the values of the interfragments angles and dihedrals.
        """
        closest_to_RA = [[],[]] #each list contains the 2 (or less) closest atoms to the reactive atom in their respective fragment
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
        if len(closest_to_RA[0]) != 0:
            point_A = self.species.geom[closest_to_RA[0][0]]
            point_B = self.species.geom[self.instance[0]]
            point_C = self.species.geom[self.instance[1]]
            changes.append([closest_to_RA[0][0],
                                   self.instance[0],
                                   self.instance[1],
                                   np.degrees(geometry.calc_angle(point_A, point_B, point_C))%360])
        if len(closest_to_RA[1]) != 0:
            point_A = self.species.geom[closest_to_RA[1][0]]
            point_B = self.species.geom[self.instance[1]]
            point_C = self.species.geom[self.instance[0]]
            changes.append([closest_to_RA[1][0],
                                   self.instance[1],
                                   self.instance[0],
                                   np.degrees(geometry.calc_angle(point_A, point_B, point_C))%360])
        
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
                            geometry.calc_dihedral(point_A, point_B, point_C, point_D)[0]%360])
        if len(closest_to_RA[0]) == 2:
            point_A = self.species.geom[closest_to_RA[0][1]]
            point_B = self.species.geom[closest_to_RA[0][0]]
            point_C = self.species.geom[self.instance[0]]
            point_D = self.species.geom[self.instance[1]]
            changes.append([closest_to_RA[0][1],
                            closest_to_RA[0][0],
                            self.instance[0],
                            self.instance[1],
                            geometry.calc_dihedral(point_A, point_B, point_C, point_D)[0]%360])
        if len(closest_to_RA[1]) == 2:
            point_A = self.species.geom[closest_to_RA[1][1]]
            point_B = self.species.geom[closest_to_RA[1][0]]
            point_C = self.species.geom[self.instance[1]]
            point_D = self.species.geom[self.instance[0]]
            changes.append([closest_to_RA[1][1],
                            closest_to_RA[1][0],
                            self.instance[1],
                            self.instance[0],
                            geometry.calc_dihedral(point_A, point_B, point_C, point_D)[0]%360])
        
        return changes

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
    