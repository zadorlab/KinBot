from kinbot.reac_General import GeneralReac
from kinbot import utils
from kinbot import constants
from kinbot.stationary_pt import StationaryPoint
from kinbot.molpro import Molpro
from kinbot import geometry
from kinbot import modify_geom
from kinbot import zmatrix
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
        #self.scan is a dict such as each key is a point of the scan.
        #Each point value is a dict itself such as:
        #self.scan = {
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
        self.relaxed = {}
        self.frozen = {
            "fragments": {
                "L1": {},
                "L2": {},
                "L3":{}
                }
                }
        self.frozen_coord = {"L1": [], "L2": []}
        self.e_in_kcal = {}
        #Reset the species to have the broken bond between fragments
        self.species, self.instance, self.shortest = self.find_bond_to_scan()
        self.set_scan_list() #Internally, the list is in bhor
        self.max_step = len(self.scan_list) - 1 #Step starts at 0
        self.removed = []
        self.points_to_remove = []
    
    def set_scan_list(self):
        if self.par["vrc_tst_scan_parameters"]["distances"] == None or\
           not isinstance(self.par["vrc_tst_scan_parameters"]["distances"], list):
            if self.shortest*constants.BOHRtoANGSTROM < self.par["vrc_tst_scan_parameters"]["start"]:
                self.scan_list = np.arange(self.shortest*constants.BOHRtoANGSTROM,\
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
            if self.shortest*constants.BOHRtoANGSTROM < self.par["vrc_tst_scan_parameters"]["distances"][0]:
                self.scan_list = np.arange(self.shortest*constants.BOHRtoANGSTROM,\
                                           self.par["vrc_tst_scan_parameters"]["distances"][0],\
                                           self.par["vrc_tst_scan_parameters"]["step"])
                self.scan_list = np.append(self.scan_list,\
                                           self.par["vrc_tst_scan_parameters"]["distances"])
            else:
                self.scan_list = self.par["vrc_tst_scan_parameters"]["distances"]
        for index, dist in enumerate(self.scan_list):
            if dist < self.shortest * constants.BOHRtoANGSTROM:
                self.scan_list.pop(index)
        self.scan_list = self.scan_list / constants.BOHRtoANGSTROM

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
        e_in_kcal = self.get_e_in_kcal(level)
        logger.info('\tSuccessful scan for {}.'.format(self.instance_name))
        logger.info(f"\tEnergies: {e_in_kcal}")
        logger.info(f"Points removed: {self.removed}")
        self.print_scan_results(level=level)
        if level == "L2":
            for point in self.relaxed:
                self.relaxed[point]["molp"] = Molpro(self.relaxed[point]["stationary_point"])
                self.relaxed[point]["molp"].create_molpro_input(VTS=True,\
                                                                name=self.relaxed[point]["stationary_point"])
    
    def filter_points(self, level="L1"):
        for point in self .relaxed:
            if self.relaxed[f"{point}"]["energy"][level] == 0.0:
                self.points_to_remove.append(point)
        for point in self.points_to_remove:
            self.relaxed.pop(str(point))
            self.scan_list = np.delete(self.scan_list, int(point))
            if point not in self.removed:
                self.removed.append(point)
        self.points_to_remove = []

    def print_scan_results(self, level="L1"):
        x_label = f"{self.species.atom[self.instance[0]]}{self.instance[0]} - {self.species.atom[self.instance[1]]}{self.instance[1]} distance ($\AA$)"
        y_label = f"Energy $(kcal/mol)$"
        y = list([self.get_e_in_kcal(level="L1")])
        x = list(np.array(self.scan_list) * constants.BOHRtoANGSTROM)
        data_legends = []
        data_legends.append(f"{self.qc.VTS_methods['L1']}/{self.qc.VTS_basis['L1']}")
        if level == "L2":
            y.append(list(self.get_e_in_kcal(level="L2")))
            data_legends.append(f"{self.qc.VTS_methods['L2']}/{self.qc.VTS_basis['L2']}")
        utils.create_matplotlib_graph(x = x, data = y, name=f"{self.instance_name}", x_label=x_label, y_label=y_label, data_legends=data_legends)

    def get_e_in_kcal(self, level="L1"):
        e_in_kcal = [constants.AUtoKCAL * (self.relaxed[f"{point}"]["energy"][level] - 
                                self.relaxed["0"]["energy"][level]) 
                    for point in self.relaxed]
        e_in_kcal = list(np.round(e_in_kcal, 2))
        return e_in_kcal
    
    def prepare_frozen_scan(self):
        fragments, self.frozen_maps = self.species.start_multimolecular()
        #Recover the LR frozen fragments geometries
        for index, frag in enumerate(fragments):
            frag.characterize()
            frag.name = f"{frag.chemid}_well_VTS"
            err, geom = self.qc.get_qc_geom(frag.name, frag.natom)
            if err == 0:
                frag.geom = geom
                frag.characterize()
                broken_frag, broken_maps = frag.start_multi_molecular()
                if len(broken_frag) != 1:
                    logger.warning(f"Optimization of {frag.name} has failed. Cancelling scan {self.name}.")
                    self.scan = 0
                    break
                self.frozen["fragments"]["L1"][f"{index}"] = frag
            else:
                self.qc.qc_opt(self.species, self.species.geom, ext="_well_VTS")
            if self.par["high_level"]:
                err, geom = self.qc.get_qc_geom(f"{frag.name}_high", frag.natom)
                if err == 0:
                    self.frag.geom = geom
                    self.frag.characterize()
                    broken_frag, broken_maps = frag.start_multi_molecular()
                    if len(broken_frag) != 1:
                        logger.warning(f"Optimization of {frag.name} has failed. Cancelling scan {self.name}.")
                        self.scan = 0
                        break
                    self.frozen["fragments"]["L2"][f"{index}"] = frag
                    self.frozen["fragments"]["L3"][f"{index}"] = frag
                else:
                    self.qc.qc_opt(self.species, self.species.geom, ext="_well_VTS_high", high_level=1)

        #Create initial geometry
        if self.scan == 1:
            #Check that both fragments are present and optimized
            fragment_missing = 0
            for index in range(2):
                if str(index) not in self.frozen["fragments"]["L1"]:
                    fragment_missing = 1

            #Create starting bimolecular stationary_point from optimized fragments
            if not fragment_missing:
                self.set_frozen_coord("L1")
                if self.par["high_level"]:
                    self.set_frozen_coord("L2")
                #Create the list of constrain for frozen fragments
                changes = [self.instance.append(self.shortest)] #This is the initial distance between the fragments
                changes.append(self.frozen_coord["L1"])

                axis = self.species.geom[self.instance[1]] - self.species.geom[self.instance[0]]
                shift = [0, 10] #Distance in bhor to move each fragment along the bond axis
                x_shift, y_shift, z_shift = np.array(shift) * geometry.unit_vector(axis)
                bimolecular_geom = np.empty(0,3)
                bimolecular_atom = []
                for line_number in range(self.species.natom):#index, frag_map in enumerate(maps):
                    for frag_number, frag_map in enumerate(self.frozen_maps):
                        if line_number in frag_map:
                            fragment = self.frozen["fragments"]["L1"][f"{frag_number}"]
                            x, y, z = fragment.geom[frag_map.index(line_number)]
                            np.append(bimolecular_geom, [x + x_shift[frag_number], y + y_shift[frag_number], z + z_shift[frag_number]])
                            bimolecular_atom += fragment.atom[frag_map.index(line_number)]
                self.frozen_species = StationaryPoint(f"{self.instance_name}_frozen",\
                                                    self.species.charge,\
                                                    self.species.mult,\
                                                    natom=self.species.natom,\
                                                    atom=bimolecular_atom,\
                                                    geom=bimolecular_geom)
                self.frozen_species.characterize()
                success, geom = modify_geom.modify_coordinates(self.frozen_species,\
                                                       f"{self.instance_name}_frozen",\
                                                       self.frozen_species.geom,\
                                                       changes,\
                                                       self.frozen_species.bond)
                if success:
                    self.frozen_species.geom = geom
                    self.frozen_species.characterize() #Frozen species is the starting geom of L1 frozen_scan

    def set_frozen_coord(self, level="L1"):
        parameters = [] #[[[dist_frag0][ang_frag0][dihed_frag0]][[dist_frag1][ang_frag1][dihed_frag1]]]
        for index in range(2):
            parameters.append([geometry.make_simple_zmat_from_cart(self.frozen["fragments"][level][f"{index}"])])
        
        for index, frag in self.frozen["fragments"][level].items():
                    for atom1 in range(0, frag.natom-2):
                        for atom2 in range(atom1+1, frag.natom-1):
                            #Constrain fragments bonds
                            self.frozen_coord[level].append([self.frozen_maps[int(index)][atom1], self.frozen_maps[int(index)][atom2], parameters[int(index)][0][atom1]])
                            if atom2 < frag.natom-1:
                                for atom3 in range(atom2+1, frag.natom-1):
                                    #Constrain fragments angles
                                    self.frozen_coord[level].append([self.frozen_maps[int(index)][atom1], self.frozen_maps[int(index)][atom2], self.frozen_maps[int(index)][atom3], parameters[int(index)][1][atom1]])
                                    if atom3 < frag.natom-1:
                                        for atom4 in range(atom3+1, frag.natom-1):
                                            #Constrain fragments dihedrals
                                            self.frozen_coord[level].append([self.frozen_maps[int(index)][atom1], self.frozen_maps[int(index)][atom2], self.frozen_maps[int(index)][atom3], self.frozen_maps[int(index)][atom4], parameters[int(index)][1][atom1]])

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []

        if step < self.max_step:
            if step == 0:
                delta = 0
            else:
                delta = self.scan_list[step] - self.scan_list[step-1]

            val = np.linalg.norm(geom[self.instance[0]] - geom[self.instance[1]]) + delta
            self.set_bond(0, 1, val, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release
