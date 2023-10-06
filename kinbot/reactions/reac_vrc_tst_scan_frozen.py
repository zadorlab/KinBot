from kinbot.reactions.reac_vrc_tst_scan import VrcTstScan
from kinbot.stationary_pt import StationaryPoint
from kinbot import zmatrix
from kinbot.molpro import Molpro
from kinbot import geometry
from kinbot import modify_geom
import numpy as np
import logging
import copy

logger = logging.getLogger('KinBot')

class VrcTstScanFrozen(VrcTstScan):
    def __init__(self, species, qc, par, instance, instance_name):
        super().__init__(species, qc, par, instance, instance_name)

        self.frozen = {
            "fragments": {
                "L1": {},
                "L2": {},
                "L3": {}
                }
                }
        self.frozen_coord = {"L1": {}, "L2": {}}
        self.prepare_frozen_scan()
        
        #self.family_name = 'VrcTstScan'
        
    def prepare_frozen_scan(self):
        fragments, self.frozen_maps = self.species.start_multi_molecular()
        fragments_optimized = 0
        #Recover the LR frozen fragments geometries
        while fragments_optimized < len(fragments):
            fragments_optimized = 0
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
                    if not self.par["high_level"]:
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
                            fragments_optimized = len(fragments)
                            break
                        self.frozen["fragments"]["L2"][f"{index}"] = frag
                        self.frozen["fragments"]["L3"][f"{index}"] = frag
                        fragments_optimized += 1
                    else:
                        self.qc.qc_opt(frag, frag.geom, ext="_well_VTS_high", high_level=1)


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
                changes = [copy.copy(self.instance)]
                changes[0].append(self.shortest) #This is the initial distance between the fragments
                for frozen_coords in self.frozen_coord["L1"].values():
                    changes.extend(frozen_coords)

                axis = self.species.geom[self.instance[1]] - self.species.geom[self.instance[0]]
                shift = [0, 10] #Distance in bhor to move each fragment along the bond axis
                x_shift, y_shift, z_shift = np.multiply(np.reshape(geometry.unit_vector(axis), (3, 1)), np.array(shift))
                bimolecular_geom = np.empty((self.species.natom,3))
                bimolecular_atom = []
                #recreate bimolecular system with frag geometries in same order as parent
                for line_number in range(self.species.natom):
                    for frag_number, frag_map in enumerate(self.frozen_maps):
                        if line_number in frag_map:
                            fragment = self.frozen["fragments"]["L1"][f"{frag_number}"]
                            x, y, z = fragment.geom[np.where(frag_map == line_number)[0][0]]
                            bimolecular_geom[line_number] = np.array([x + x_shift[frag_number], y + y_shift[frag_number], z + z_shift[frag_number]])
                            bimolecular_atom += fragment.atom[np.where(frag_map == line_number)[0][0]]
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
                    self.species = self.frozen_species

    def set_frozen_coord(self, level="L1"):
        parameters = {} #[[[dist_frag0][ang_frag0][dihed_frag0]][[dist_frag1][ang_frag1][dihed_frag1]]]
        for index in range(2):
            parameters[f"{index}"]= {}
            parameters[f"{index}"]["values"], parameters[f"{index}"]["indexes"] = zmatrix.make_simple_zmat_from_cart(self.frozen["fragments"][level][f"{index}"])
        
        for index, frag in self.frozen["fragments"][level].items(): #index is the fragment number  
            self.frozen_coord[level][index] = []
            for atom_index in range(1, frag.natom-1):
                atom1 = atom_index
                atom2 = atom_index-1
                #Constrain fragments bonds
                self.frozen_coord[level][index].append([self.frozen_maps[int(index)][atom2], self.frozen_maps[int(index)][atom1], parameters[index]["values"][0][atom2]])
                if atom_index > 1:
                    atom3 = atom_index-2
                    #Constrain fragments angles
                    self.frozen_coord[level][index].append([self.frozen_maps[int(index)][atom3], self.frozen_maps[int(index)][atom2], self.frozen_maps[int(index)][atom1], parameters[index]["values"][1][atom3]])
                    if atom_index > 2:
                        atom4 = atom_index-3
                        #Constrain fragments dihedrals
                        self.frozen_coord[level][index].append([self.frozen_maps[int(index)][atom4], self.frozen_maps[int(index)][atom3], self.frozen_maps[int(index)][atom2], self.frozen_maps[int(index)][atom1], parameters[index]["values"][2][atom4]])

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
        change.append(self.frozen_coord["L1"])
        
        return step, fix, change, release