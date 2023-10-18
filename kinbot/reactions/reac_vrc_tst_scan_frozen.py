from kinbot.reactions.reac_vrc_tst_scan import VrcTstScan
from kinbot import geometry
from kinbot.molpro import Molpro
from sella import Internals
from kinbot import geometry
import numpy as np
import logging
import copy
from ase import Atoms

logger = logging.getLogger('KinBot')

class VrcTstScanFrozen(VrcTstScan):
    def __init__(self, species, qc, par, instance, instance_name):
        super().__init__(species, qc, par, instance, instance_name)
        """
        Frozen fragments are updated at each step of the scan.
        The internal geometry shouldn't change, but 
        """

        self.frozen = {
            "coord":{
                "L1": {},
                "L2": {}
                },
            "geom":{
                "L1": {},
                "L2": {}
                }
                }#maps is added just below
        self.prepare_frozen_scan()
        
    def prepare_frozen_scan(self):
        fragments, self.frozen["maps"] = self.species.start_multi_molecular()
        fragments_optimized = 0
        hl_opt = 0
        #Recover the LR frozen fragments geometries
        while fragments_optimized < len(fragments) or hl_opt < len(fragments):
            fragments_optimized = 0
            hl_opt = 0
            for frag_number, frag in enumerate(fragments):
                frag.characterize()
                frag.name = f"{frag.chemid}_well_VTS"
                err, geom = self.qc.get_qc_geom(frag.name, frag.natom)
                if err == 0:
                    frag.geom = geom
                    frag.characterize()
                    broken_frag, broken_maps = frag.start_multi_molecular()
                    if len(broken_frag) != 1:
                        logger.warning(f"Optimization of {frag.name} has failed. Cancelling scan {self.name}.")
                        fragments_optimized = len(fragments)
                        self.scan = 0
                        break
                    frag_atom = Atoms(frag.atom, positions=frag.geom)

                    #Only define internals here. Use same set for L2.
                    internals = Internals(frag_atom)
                    internals.find_all_bonds()
                    internals.find_all_angles()
                    internals.find_all_dihedrals()
                    self.set_frozen_coord(level="L1", frag=frag, internals=internals, frag_number=frag_number)#Must only be called once before start of the scan
                    fragments_optimized += 1
                    #Save the cartesian coordinates
                    self.frozen["geom"]["L1"][f"{frag_number}"] = frag.geom
                else:
                    self.qc.qc_opt(frag, frag.geom, ext="_well_VTS")
                if self.par["high_level"]:
                    err, geom = self.qc.get_qc_geom(f"{frag.name}_well_VTS_high", frag.natom)
                    if err == 0:
                        frag.geom = geom
                        frag.characterize()
                        broken_frag, broken_maps = frag.start_multi_molecular()
                        if len(broken_frag) != 1:
                            logger.warning(f"Optimization of {frag.name} has failed. Cancelling scan {self.name}.")
                            self.scan = 0
                            hl_opt = len(fragments)
                            break
                        self.set_frozen_coord(level="L2", frag=frag, internals=internals, frag_number=frag_number) #Must only be called once before start of the scan
                        hl_opt += 1
                        #Save the cartesian coordinates
                        self.frozen["geom"]["L2"][f"{frag_number}"] = frag.geom
                    else:
                        self.qc.qc_opt(frag, frag.geom, ext="_well_VTS_high", high_level=1)

            if not self.par["high_level"]:
                hl_opt = fragments_optimized

        if self.scan == 1:
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

        geometries = [geom for geom in self.frozen["geom"][level].values()]

        #modify dihedrals and angles
        for modif in reversed(interfragments_param):
            match len(modif):
                case 5:
                    frag_number = [[frag_index for frag_index, fmap in enumerate(self.frozen["maps"]) if atom_index in fmap][0] for atom_index in modif[:-1]]

                    #Center the fragment coordinates on the atom to rotate around
                    origin = geometries[frag_number[0]][np.asarray(self.frozen["maps"][frag_number[0]] == self.instance[frag_number[0]]).nonzero()][0]
                    coord = geometries[frag_number[0]] - origin

                    point_A = geometries[frag_number[0]][np.asarray(self.frozen["maps"][frag_number[0]] == modif[0]).nonzero()][0]
                    point_B = geometries[frag_number[1]][np.asarray(self.frozen["maps"][frag_number[1]] == modif[1]).nonzero()][0]
                    point_C = geometries[frag_number[2]][np.asarray(self.frozen["maps"][frag_number[2]] == modif[2]).nonzero()][0]
                    point_D = geometries[frag_number[3]][np.asarray(self.frozen["maps"][frag_number[3]] == modif[3]).nonzero()][0]

                    axis = geometry.unit_vector(point_C - point_B)

                    current_dihedral = geometry.calc_dihedral(point_A, point_B, point_C, point_D)[0]%360
                    rotation_angle = np.radians(modif[-1] - current_dihedral)

                    geometries[frag_number[0]] = np.dot(coord, geometry.rotation_matrix(axis, rotation_angle)) + origin
                case 4:
                    frag_number = [[frag_index for frag_index, fmap in enumerate(self.frozen["maps"]) if atom_index in fmap][0] for atom_index in modif[:-1]]

                    #Center the fragment coordinates on the atom to rotate around
                    origin = geometries[frag_number[0]][np.asarray(self.frozen["maps"][frag_number[0]] == self.instance[frag_number[0]]).nonzero()][0]
                    coord = geometries[frag_number[0]] - origin

                    point_A = geometries[frag_number[0]][np.asarray(self.frozen["maps"][frag_number[0]] == modif[0]).nonzero()][0] - origin
                    point_B = geometries[frag_number[1]][np.asarray(self.frozen["maps"][frag_number[1]] == modif[1]).nonzero()][0] - origin
                    point_C = geometries[frag_number[2]][np.asarray(self.frozen["maps"][frag_number[2]] == modif[2]).nonzero()][0] - origin

                    axis = geometry.unit_vector(np.cross(point_B - point_A, point_B - point_C))

                    current_angle = np.degrees(geometry.calc_angle(point_A, point_B, point_C))%360
                    rotation_angle = np.radians(modif[-1] - current_angle)

                    geometries[frag_number[0]] = np.dot(coord, geometry.rotation_matrix(axis, rotation_angle)) + origin
        #Set inter-fragment distance
        frag_number = [[frag_index for frag_index, fmap in enumerate(self.frozen["maps"]) if atom_index in fmap][0] for atom_index in self.instance]
        point_A = geometries[frag_number[0]][np.asarray(self.frozen["maps"][frag_number[0]] == self.instance[0]).nonzero()][0]
        point_B = geometries[frag_number[1]][np.asarray(self.frozen["maps"][frag_number[1]] == self.instance[1]).nonzero()][0]

        current_distance = point_B - point_A
        shift = geometry.unit_vector(current_distance)*distance - current_distance
        geometries[frag_number[1]] += shift

        #recreate bimolecular system with frag geometries in same order as parent
        bimolecular_geom = np.empty((self.species.natom,3))
        for line_number in range(self.species.natom):
            for frag_number, frag_map in enumerate(self.frozen["maps"]):
                if line_number in frag_map:
                    geom = geometries[frag_number]
                    bimolecular_geom[line_number] = geom[np.where(frag_map == line_number)[0][0]]
        tmp_species = copy.deepcopy(self.species)
        tmp_species.geom = bimolecular_geom
        tmp_species.characterize()

        return tmp_species

    def set_frozen_coord(self, level="L1", frag=None, internals=None, frag_number=None):
        """
        Function which saves the internal parameters of the individually optimized fragments.
        The internal coordinates are defined by Sella.
        self.frozen["coord"] is specific to frozen scans
        self.frozen_param exists both in frozen and relaxed scans and is used reaction generator.
        """

        #Save the bonds
        self.frozen["coord"][level][f"{frag_number}"] = []
        for bond in internals.internals["bonds"]:
            point_A = frag.geom[bond.indices[0]]
            point_B = frag.geom[bond.indices[1]]
            bond_length = np.linalg.norm(point_B - point_A)
            self.frozen["coord"][level][f"{frag_number}"].append([self.frozen["maps"][frag_number][bond.indices[0]],
                                                                  self.frozen["maps"][frag_number][bond.indices[1]],
                                                                  bond_length])
            if level == "L1":
                self.frozen_param.append([self.frozen["maps"][frag_number][bond.indices[0]]+1,
                                          self.frozen["maps"][frag_number][bond.indices[1]]+1])

        #Save the angles
        for angle in internals.internals["angles"]:
            point_A = frag.geom[angle.indices[0]]
            point_B = frag.geom[angle.indices[1]]
            point_C = frag.geom[angle.indices[2]]
            angle_value = np.degrees(geometry.calc_angle(point_A, point_B, point_C))
            self.frozen["coord"][level][f"{frag_number}"].append([self.frozen["maps"][frag_number][angle.indices[0]],
                                                                  self.frozen["maps"][frag_number][angle.indices[1]],
                                                                  self.frozen["maps"][frag_number][angle.indices[2]],
                                                                  angle_value])
            if level == "L1":
                self.frozen_param.append([self.frozen["maps"][frag_number][angle.indices[0]]+1,
                                          self.frozen["maps"][frag_number][angle.indices[1]]+1,
                                          self.frozen["maps"][frag_number][angle.indices[2]]+1])

        #Save the dihedrals
        for dihedral in internals.internals["dihedrals"]:
            point_A = frag.geom[dihedral.indices[0]]
            point_B = frag.geom[dihedral.indices[1]]
            point_C = frag.geom[dihedral.indices[2]]
            point_D = frag.geom[dihedral.indices[3]]
            dihedral_value = geometry.calc_dihedral(point_A, point_B, point_C, point_D)[0]
            self.frozen["coord"][level][f"{frag_number}"].append([self.frozen["maps"][frag_number][dihedral.indices[0]],
                                                                  self.frozen["maps"][frag_number][dihedral.indices[1]],
                                                                  self.frozen["maps"][frag_number][dihedral.indices[2]],
                                                                  self.frozen["maps"][frag_number][dihedral.indices[3]],
                                                                  dihedral_value])
            if level == "L1":
                self.frozen_param.append([self.frozen["maps"][frag_number][dihedral.indices[0]]+1,
                                          self.frozen["maps"][frag_number][dihedral.indices[1]]+1,
                                          self.frozen["maps"][frag_number][dihedral.indices[2]]+1,
                                          self.frozen["maps"][frag_number][dihedral.indices[3]]+1])
                
    def get_interfragments_param(self):
        """
        Function that selects the neighbourgs of the closest atom in each fragment,
        and returns the values of the interfragments angles and dihedrals.
        """
        closest_to_RA = [[],[]] #each list contains the 2 (or less) closest atoms to the reactive atom in their respective fragment
        fragments, maps = self.species.start_multi_molecular()
        for frag_number, ra in enumerate(self.instance):
            for neighbourg_index, bond in zip(maps[frag_number], fragments[frag_number].bond[ra]):
                if bond != 0:
                    closest_to_RA[frag_number].append(neighbourg_index)
                    if len(closest_to_RA[frag_number]) == 2:
                        break
            #If no neighbourgs atoms were found, the fragment has to be mono-atomic
            if len(closest_to_RA[frag_number]) == 1 and fragments[frag_number].natom > 2:
                for neighbourg_index, bond in zip(maps[frag_number], fragments[frag_number].bond[closest_to_RA[frag_number][0]]):
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
        if len(closest_to_RA[0]) != 0 and len(closest_to_RA[1]) != 0:
            point_A = self.species.geom[closest_to_RA[0][0]]
            point_B = self.species.geom[self.instance[0]]
            point_C = self.species.geom[self.instance[1]]
            point_D = self.species.geom[closest_to_RA[1][0]]
            changes.append([closest_to_RA[0][0],
                            self.instance[0],
                            self.instance[1],
                            closest_to_RA[1][0],
                            np.degrees(geometry.calc_dihedral(point_A, point_B, point_C, point_D)[0])%360])
        if len(closest_to_RA[0]) == 2:
            point_A = self.species.geom[closest_to_RA[0][1]]
            point_B = self.species.geom[closest_to_RA[0][0]]
            point_C = self.species.geom[self.instance[0]]
            point_D = self.species.geom[self.instance[1]]
            changes.append([closest_to_RA[0][1],
                            closest_to_RA[0][0],
                            self.instance[0],
                            self.instance[1],
                            np.degrees(geometry.calc_dihedral(point_A, point_B, point_C, point_D)[0])%360])
        if len(closest_to_RA[1]) == 2:
            point_A = self.species.geom[closest_to_RA[1][1]]
            point_B = self.species.geom[closest_to_RA[1][0]]
            point_C = self.species.geom[self.instance[1]]
            point_D = self.species.geom[self.instance[0]]
            changes.append([closest_to_RA[1][1],
                            closest_to_RA[1][0],
                            self.instance[1],
                            self.instance[0],
                            np.degrees(geometry.calc_dihedral(point_A, point_B, point_C, point_D)[0])%360])
        
        return changes

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []

        val = self.scan_list[step]
        self.set_bond(0, 1, val, change)

        self.clean_constraints(change, fix)
        for frozen_coords in self.frozen["coord"]["L1"].values():
            for c in frozen_coords:
                c_new = [ci + 1 for ci in c[:-1]]
                c_new.append(c[-1])
                change.append(c_new)

        return step, fix, change, release
    