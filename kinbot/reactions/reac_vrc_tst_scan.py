from kinbot.reac_General import GeneralReac
from kinbot import utils
from kinbot import constants
from kinbot.stationary_pt import StationaryPoint
from kinbot.molpro import Molpro
from kinbot import geometry
from sella import Internals
from ase import Atoms
import numpy as np
import copy
import logging

logger = logging.getLogger('KinBot')

class VrcTstScan(GeneralReac):
    def __init__(self, species, qc, par, instance, instance_name):
        super().__init__(species, qc, par, instance, instance_name)
        self.scan = 1
        self.skip = 0
        self.family_name = 'VrcTstScan'
        self.instance_basename = self.instance_name
        #self.scanned is a dict such as each key is a point of the scan.
        #Each point value is a dict itself such as:
        #self.scanned = {
        #   "0": {
        #       "energy": {
        #           "L1": energyL1
        #           "L2": energyL2
        #           "L3": energyL3
        #            },
        #       "stationary_point": StationaryPoint object used for this point
        #       "opt": Optimization object only doing high level
        #       "molp": Molpro object
        #       },
        #   "1": {...}
        #   }
        self.scanned = {}
        self.long_range = {
            "coord":{
                "L1": {},
                "L2": {}
                },
            "geom":{
                "L1": {},
                "L2": {}
                },
            "energies":{
                "L1": {},
                "L2": {},
                "L3": {}
                }
                }
        self.e_in_kcal = {}
        #Reset the species to have the broken bond between fragments
        self.species, self.instance, self.shortest = self.find_bond_to_scan()
        self.set_scan_list() #Internally, the list is in bhor
        self.max_step = len(self.scan_list) - 1 #Step starts at 0
        self.removed = []
        self.points_to_remove = []
        self.frozen_param = [[ index+1 for index in self.instance]]
        self.prepare_scan()

    def prepare_scan(self):
        fragments, self.long_range["maps"] = self.species.start_multi_molecular()
        fragments_optimized = 0
        hl_opt = 0
        #Recover the LR frozen fragments geometries
        while fragments_optimized < len(fragments) or hl_opt < len(fragments):
            fragments_optimized = 0
            hl_opt = 0
            for frag_number, frag in enumerate(fragments):
                frag.characterize()
                err, geom = self.qc.get_qc_geom(f"{frag.chemid}_well_VTS", frag.natom)
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
                    self.long_range["geom"]["L1"][f"{frag_number}"] = frag.geom
                    err, self.long_range["energies"]["L1"][f"{frag_number}"] = self.qc.get_qc_energy(f"{frag.chemid}_well_VTS")
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
                        self.long_range["geom"]["L2"][f"{frag_number}"] = frag.geom
                        err, self.long_range["energies"]["L2"][f"{frag_number}"] = self.qc.get_qc_energy(f"{frag.chemid}_well_VTS_high")
                    else:
                        self.qc.qc_opt(frag, frag.geom, ext="_well_VTS_high", high_level=1)

                    # write the L3 input and read the L3 energy, if available
                    if self.par['L3_calc'] == 1:
                        key = self.par['vrc_tst_scan_parameters']["molpro_key"].upper()
                        molp = Molpro(self.species, self.par)
                        molp.create_molpro_input(name=f"{frag.chemid}_well_VTS", VTS=True)
                        status, molpro_energy = molp.get_molpro_energy(key, name=f"{frag.chemid}_well_VTS")
                        if status:
                            self.long_range["energies"]["L3"][f"{frag_number}"] = molpro_energy
                            self.try_l3 = True
                        else:
                            self.try_l3 = False

            if not self.par["high_level"]:
                hl_opt = fragments_optimized

    def set_frozen_coord(self, level="L1", frag=None, internals=None, frag_number=None):
        """
        Function which saves the internal parameters of the individually optimized fragments.
        The internal coordinates are defined by Sella.
        self.long_range["coord"] is specific to frozen scans
        self.frozen_param exists both in frozen and relaxed scans and is used reaction generator.
        """

        #Save the bonds
        self.long_range["coord"][level][f"{frag_number}"] = []
        for bond in internals.internals["bonds"]:
            point_A = frag.geom[bond.indices[0]]
            point_B = frag.geom[bond.indices[1]]
            bond_length = np.linalg.norm(point_B - point_A)
            self.long_range["coord"][level][f"{frag_number}"].append([self.long_range["maps"][frag_number][bond.indices[0]],
                                                                      self.long_range["maps"][frag_number][bond.indices[1]],
                                                                      bond_length])
            if level == "L1":
                self.frozen_param.append([self.long_range["maps"][frag_number][bond.indices[0]]+1,
                                          self.long_range["maps"][frag_number][bond.indices[1]]+1])

        #Save the angles
        for angle in internals.internals["angles"]:
            point_A = frag.geom[angle.indices[0]]
            point_B = frag.geom[angle.indices[1]]
            point_C = frag.geom[angle.indices[2]]
            angle_value = np.degrees(geometry.calc_angle(point_A, point_B, point_C))
            self.long_range["coord"][level][f"{frag_number}"].append([self.long_range["maps"][frag_number][angle.indices[0]],
                                                                  self.long_range["maps"][frag_number][angle.indices[1]],
                                                                  self.long_range["maps"][frag_number][angle.indices[2]],
                                                                  angle_value])
            if level == "L1":
                self.frozen_param.append([self.long_range["maps"][frag_number][angle.indices[0]]+1,
                                          self.long_range["maps"][frag_number][angle.indices[1]]+1,
                                          self.long_range["maps"][frag_number][angle.indices[2]]+1])

        #Save the dihedrals
        for dihedral in internals.internals["dihedrals"]:
            point_A = frag.geom[dihedral.indices[0]]
            point_B = frag.geom[dihedral.indices[1]]
            point_C = frag.geom[dihedral.indices[2]]
            point_D = frag.geom[dihedral.indices[3]]
            dihedral_value = geometry.calc_dihedral(point_A, point_B, point_C, point_D)[0]
            self.long_range["coord"][level][f"{frag_number}"].append([self.long_range["maps"][frag_number][dihedral.indices[0]],
                                                                  self.long_range["maps"][frag_number][dihedral.indices[1]],
                                                                  self.long_range["maps"][frag_number][dihedral.indices[2]],
                                                                  self.long_range["maps"][frag_number][dihedral.indices[3]],
                                                                  dihedral_value])
            if level == "L1":
                self.frozen_param.append([self.long_range["maps"][frag_number][dihedral.indices[0]]+1,
                                          self.long_range["maps"][frag_number][dihedral.indices[1]]+1,
                                          self.long_range["maps"][frag_number][dihedral.indices[2]]+1,
                                          self.long_range["maps"][frag_number][dihedral.indices[3]]+1])

    def set_scan_list(self):
        if self.par["vrc_tst_scan_parameters"]["distances"] == None or\
           not isinstance(self.par["vrc_tst_scan_parameters"]["distances"], list):
            if self.shortest < self.par["vrc_tst_scan_parameters"]["start"]:
                self.scan_list = np.arange(self.shortest,\
                                           self.par["vrc_tst_scan_parameters"]["start"],\
                                           self.par["vrc_tst_scan_parameters"]["step"])
                self.scan_list = np.append(self.scan_list,\
                                np.arange(self.par["vrc_tst_scan_parameters"]["start"],\
                                        self.par["vrc_tst_scan_parameters"]["stop"],\
                                        self.par["vrc_tst_scan_parameters"]["step"]))
            else:
                self.scan_list = np.arange(self.par["vrc_tst_scan_parameters"]["start"],\
                                        self.par["vrc_tst_scan_parameters"]["stop"],\
                                        self.par["vrc_tst_scan_parameters"]["step"])
        else:
            if self.shortest < self.par["vrc_tst_scan_parameters"]["distances"][0]:
                self.scan_list = np.arange(self.shortest,\
                                           self.par["vrc_tst_scan_parameters"]["distances"][0],\
                                           self.par["vrc_tst_scan_parameters"]["step"])
                self.scan_list = np.append(self.scan_list,\
                                           self.par["vrc_tst_scan_parameters"]["distances"])
            else:
                self.scan_list = self.par["vrc_tst_scan_parameters"]["distances"]
        for index, dist in enumerate(self.scan_list):
            if dist < self.shortest:
                self.scan_list.pop(index)
        self.scan_list = self.scan_list

    def find_bond_to_scan(self):
        reaction_name, products = self.instance
        if "hom_sci" in reaction_name: #The two numbers at the end of a hom_sci reaction should be where the bond was broken
            atoms_index = list(reaction_name.split("_")[-2:])
            initial_well = copy.deepcopy(self.species)
            initial_well.bond[atoms_index[0], atoms_index[1]] = 0
            initial_well.bond[atoms_index[1], atoms_index[0]] = 0
            fragments, maps = initial_well.start_multi_molecular()
        else: #For vdW well, find the direction, and find the closest atoms between fragments
            for direction in ["F","R"]:
                status, geom = self.qc.get_qc_geom(f"{reaction_name}_IRC_{direction}_prod", self.species.natom)
                if status == 0:
                    vdW_name = f"{reaction_name}_IRC_{direction}_prod"
                    initial_well = StationaryPoint(vdW_name, self.species.charge, self.species.mult, natom=self.species.natom, atom=self.species.atom, geom=geom)
                    initial_well.characterize()
                    fragments, maps = initial_well.start_multi_molecular()
                    if len(fragments) == 2:
                        break
            if status != 0 : #Could not find the vdW_well, don't do the scan
                self.scan = 0
        #Additional check to see if this is the correct vdW well
        frag_chemids = [str(frag.chemid) for frag in fragments]
        prod_chemids = list(products.split("_"))
        for frag_chemid, prod_chemid in zip(sorted(frag_chemids), sorted(prod_chemids)):
            if frag_chemid != prod_chemid:
                self.scan = 0 #products from db and from input are different, don't do the scan
                break
        shortest = np.inf
        for atom1 in maps[0]:
            for atom2 in maps[1]:
                if initial_well.dist[atom1, atom2] < shortest:
                    shortest = initial_well.dist[atom1, atom2]
                    atoms_index = [atom1, atom2]
        return [initial_well, atoms_index, shortest]
    
    def finish_vrc_tst_scan(self, level):

        self.filter_points()
        if level != "L1":
            self.filter_points(level=level)
            if self.try_l3:
                key = self.par['vrc_tst_scan_parameters']["molpro_key"].upper()
                L3 = 1
                for point in self.scanned:
                    self.scanned[point]["molp"] = Molpro(self.scanned[point]["stationary_point"], self.par)
                    status, molpro_energy = self.scanned[point]["molp"].get_molpro_energy(key, name=self.scanned[point]["stationary_point"].name)
                    if not status:
                        logger.info(f"Missing {self.scanned[point]['stationary_point'].name} at L3")
                        L3 = 0
                        break
                    else:
                        self.scanned[point]["energy"]["L3"] = molpro_energy
                if L3:
                    level = "L3"
        self.removed.sort()
        for point in reversed(self.removed):
            self.scan_list = np.delete(self.scan_list, int(point))
        e_in_kcal = self.get_e_in_kcal(level)
        logger.info('\tSuccessful scan for {}.'.format(self.instance_name))
        logger.info(f"\tEnergies: {e_in_kcal}")
        logger.info(f"Points removed: {self.removed}")
        self.print_scan_results(level=level)
    
    def filter_points(self, level="L1"):
        for point in self.scanned:
            if self.scanned[f"{point}"]["energy"][level] == 0.0 and point not in self.points_to_remove and point not in self.removed:
                self.points_to_remove.append(point)
        for point in self.points_to_remove:
            self.scanned.pop(point)
            if point not in self.removed:
                self.removed.append(int(point))
        self.points_to_remove = []

    def print_scan_results(self, level="L1"):
        x_label = f"{self.species.atom[self.instance[0]]}{self.instance[0]} - {self.species.atom[self.instance[1]]}{self.instance[1]} distance ($\AA$)"
        y_label = f"Energy $(kcal/mol)$"
        y = list([self.get_e_in_kcal(level="L1")])
        x = list(np.array(self.scan_list))
        data_legends = []
        data_legends.append(f"{self.qc.VTS_methods['L1']}/{self.qc.VTS_basis['L1']}")
        if level == "L2":
            y.append(list(self.get_e_in_kcal(level="L2")))
            data_legends.append(f"{self.qc.VTS_methods['L2']}/{self.qc.VTS_basis['L2']}")
            if level == "L3":
                y.append(list(self.get_e_in_kcal(level="L3")))
                if "frozen" in self.instance_basename:
                    data_legends.append(f"{self.qc.VTS_methods['L3']}/{self.qc.VTS_basis['L3'][0]}")
                else:
                    data_legends.append(f"{self.qc.VTS_methods['L3']}/{self.qc.VTS_basis['L3'][1]}")
        utils.create_matplotlib_graph(x = x, data = y, name=f"{self.instance_basename}", x_label=x_label, y_label=y_label, data_legends=data_legends)

    def get_e_in_kcal(self, level="L1"):
        assymptote = self.long_range["energies"][level]["0"] + self.long_range["energies"][level]["1"]
        e_in_kcal = [constants.AUtoKCAL * (self.scanned[f"{point}"]["energy"][level] - assymptote)
                    for point in self.scanned]
        e_in_kcal = list(np.round(e_in_kcal, 2))
        return e_in_kcal

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []

        val = self.scan_list[step]
        self.set_bond(0, 1, val, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
