from kinbot.reac_General import GeneralReac
from kinbot import utils
from kinbot import constants
from kinbot.stationary_pt import StationaryPoint
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
        #       },
        #   "1": {...}
        #   }
        self.scanned = {}
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
    
    def filter_points(self, level="L1"):
        for point in self .scanned:
            if self.scanned[f"{point}"]["energy"][level] == 0.0:
                self.points_to_remove.append(point)
        for point in self.points_to_remove:
            self.scanned.pop(str(point))
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
        e_in_kcal = [constants.AUtoKCAL * (self.scanned[f"{point}"]["energy"][level] - 
                                self.scanned["0"]["energy"][level]) 
                    for point in self.scanned]
        e_in_kcal = list(np.round(e_in_kcal, 2))
        return e_in_kcal
    
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
