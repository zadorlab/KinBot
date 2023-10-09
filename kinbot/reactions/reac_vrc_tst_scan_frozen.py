from kinbot.reactions.reac_vrc_tst_scan import VrcTstScan
from kinbot.stationary_pt import StationaryPoint
from kinbot import zmatrix
from kinbot.molpro import Molpro
from sella import Sella, Constraints
from kinbot import geometry
from kinbot import modify_geom
import numpy as np
import logging
import copy

logger = logging.getLogger('KinBot')

class VrcTstScanFrozen(VrcTstScan):
    def __init__(self, species, qc, par, instance, instance_name):
        super().__init__(species, qc, par, instance, instance_name)
        """
        Frozen fragments are updated at each step of the scan.
        The internal geometry shouldn't change, but 
        """

        self.frozen = {
            "fragments": {
                "L1": {},
                "L2": {}
                },
            "coord":{
                "L1": {},
                "L2": {}
                }
                }#maps is added just below
        self.prepare_frozen_scan()
        
        #self.family_name = 'VrcTstScan'
        
    def prepare_frozen_scan(self):
        fragments, self.frozen["maps"] = self.species.start_multi_molecular()
        fragments_optimized = 0
        hl_opt = 0
        #Recover the LR frozen fragments geometries
        while fragments_optimized < len(fragments) and hl_opt < len(fragments):
            fragments_optimized = 0
            hl_opt = 0
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
                        fragments_optimized = len(fragments)
                        self.scan = 0
                        break
                    self.frozen["fragments"]["L1"][f"{index}"] = frag
                    fragments_optimized += 1
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
                        self.frozen["fragments"]["L2"][f"{index}"] = frag
                        hl_opt += 1
                    else:
                        self.qc.qc_opt(frag, frag.geom, ext="_well_VTS_high", high_level=1)

            if not self.par["high_level"]:
                hl_opt = fragments_optimized

        if self.scan == 1:
            #Check that both fragments are present and optimized
            self.set_frozen_coord(level="L1") #Must only be called once before start of the scan
            if self.par["high_level"]:
                self.set_frozen_coord(level="L2") #Must only be called once before start of the scan
                    
            self.species = self.get_frozen_species()

    def set_frozen_fragments(self, level="L1"):
        """
        Take the fragments of the current point(any level)
        and modify their geometry to fit the frozen parameters for the level asked.
        """
        fragments, maps = self.species.start_multi_molecular()
        for frag_number, frag in enumerate(fragments):
            frag.characterize()
            changes = []
            for coord in self.frozen["coord"][level][f"{frag_number}"]:
                parent_atom_indexes = coord[:-1]
                frag_atom_indexes = []
                for atm_idx in parent_atom_indexes:
                    if atm_idx in self.frozen["maps"][frag_number]:
                        frag_atom_indexes.append(np.where( self.frozen["maps"][frag_number] == atm_idx)[0][0])
                frag_atom_indexes.append(coord[-1])
                changes.append(frag_atom_indexes)
        retry = 0
        success = 0
        while not success or retry < 3:
            success, geom = modify_geom.modify_coordinates(frag,\
                                                       self.instance_name,\
                                                       frag.geom,\
                                                       changes,\
                                                       frag.bond,\
                                                       write_files=0)
            frag.geom = geom
            frag.characterize()
            if success:
                self.frozen["fragments"][level][f"{frag_number}"] = frag
            else:
                logger.warning("The frozen fragment {} wasn't set correctly for {}".format(frag.name, self.instance_name))
                if retry < 3:
                    retry += 1
                    logger.warning("Restarting BFGS optimization from last geometry.")
                else:
                    self.scan = 0
                    return self.species


    def get_frozen_species(self, level="L1", distance=None):
        #distance between the two fragments
        if distance == None or not isinstance(distance, float):
            distance = self.shortest

        #Create the list of constrain for frozen fragments
        changes = [copy.copy(self.instance)]
        changes[0].append(distance)
        for frozen_coords in self.frozen["coord"][level].values():
            changes.extend(frozen_coords)

        axis = self.species.geom[self.instance[1]] - self.species.geom[self.instance[0]] #Use last orientation
        shift = [0, 1] #Distance in bhor to move each fragment along the bond axis
        x_shift, y_shift, z_shift = np.multiply(np.reshape(geometry.unit_vector(axis), (3, 1)), np.array(shift))
        bimolecular_geom = np.empty((self.species.natom,3))
        bimolecular_atom = []
        self.set_frozen_fragments(level=level)
        #recreate bimolecular system with frag geometries in same order as parent
        for line_number in range(self.species.natom):
            for frag_number, frag_map in enumerate(self.frozen["maps"]):
                if line_number in frag_map:
                    fragment = self.frozen["fragments"][level][f"{frag_number}"]
                    x, y, z = fragment.geom[np.where(frag_map == line_number)[0][0]]
                    bimolecular_geom[line_number] = np.array([x + x_shift[frag_number], y + y_shift[frag_number], z + z_shift[frag_number]])
                    bimolecular_atom += fragment.atom[np.where(frag_map == line_number)[0][0]]
        tmp_frozen_species = StationaryPoint(self.species.name,\
                                            self.species.charge,\
                                            self.species.mult,\
                                            natom=self.species.natom,\
                                            atom=bimolecular_atom,\
                                            geom=bimolecular_geom)
        tmp_frozen_species.characterize()
        retry = 0
        success = 0
        while not success or retry < 3:
            geom = tmp_frozen_species.geom
            success, geom = modify_geom.modify_coordinates(tmp_frozen_species,\
                                                    self.instance_name,\
                                                    geom,\
                                                    changes,\
                                                    tmp_frozen_species.bond,\
                                                    write_files=0) #Turn 1 for debugging
            if success:
                tmp_frozen_species.geom = geom
                tmp_frozen_species.characterize() #Frozen species is the starting geom of L1 frozen_scan
                return tmp_frozen_species
            else:
                logger.warning("The frozen species wasn't set correctly for {}".format(self.instance_name))
                if retry < 3:
                    retry += 1
                    logger.warning("Restarting BFGS optimization from last geometry.")
                else:
                    self.scan = 0
                    return self.species

    def set_frozen_coord(self, level="L1"):
        parameters = {} 
        for frag_number in range(2):
            parameters[f"{frag_number}"]= {}
            parameters[f"{frag_number}"]["values"], parameters[f"{frag_number}"]["indexes"] = zmatrix.make_simple_zmat_from_cart(self.frozen["fragments"][level][f"{frag_number}"])
            for angle in parameters[f"{frag_number}"]["values"][1]:
                angle %= 360
            for angle in parameters[f"{frag_number}"]["values"][2]:
                angle %= 360

        for index, frag in self.frozen["fragments"][level].items(): #index is the fragment number  
            self.frozen["coord"][level][index] = []
            for atom_index in range(1, frag.natom-1):
                atom1 = atom_index
                atom2 = atom_index-1
                #Constrain fragments bonds
                self.frozen["coord"][level][index].append([self.frozen["maps"][int(index)][atom2], self.frozen["maps"][int(index)][atom1], parameters[index]["values"][0][atom2]])
                if level == "L1":
                    self.frozen_param.append([self.frozen["maps"][int(index)][atom2]+1, self.frozen["maps"][int(index)][atom1]+1])
                if atom_index > 1:
                    atom3 = atom_index-2
                    #Constrain fragments angles
                    self.frozen["coord"][level][index].append([self.frozen["maps"][int(index)][atom3], self.frozen["maps"][int(index)][atom2], self.frozen["maps"][int(index)][atom1], parameters[index]["values"][1][atom3]])
                    if level == "L1":
                        self.frozen_param.append([self.frozen["maps"][int(index)][atom3]+1, self.frozen["maps"][int(index)][atom2]+1, self.frozen["maps"][int(index)][atom1]+1])
                    if atom_index > 2:
                        atom4 = atom_index-3
                        #Constrain fragments dihedrals
                        self.frozen["coord"][level][index].append([self.frozen["maps"][int(index)][atom4], self.frozen["maps"][int(index)][atom3], self.frozen["maps"][int(index)][atom2], self.frozen["maps"][int(index)][atom1], parameters[index]["values"][2][atom4]])
                        if level == "L1":
                            self.frozen_param.append([self.frozen["maps"][int(index)][atom4]+1, self.frozen["maps"][int(index)][atom3]+1, self.frozen["maps"][int(index)][atom2]+1, self.frozen["maps"][int(index)][atom1]+1])

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