import os
from typing import Any
from kinbot import kb_path
from kinbot import constants
from kinbot.stationary_pt import StationaryPoint
from ase.atoms import Atoms
import numpy as np
from numpy.typing import NDArray
from numpy import float32, float64, int64, int16, bool_, pi, floating
from kinbot import geometry
import logging
import copy
import math
from kinbot.constants import BOHRtoANGSTROM
from kinbot.utils import reorder_coord

logger = logging.getLogger('KinBot')


class Fragment(StationaryPoint):
    """
    Class that creates stationary points
    with specific methods to setup VRC TST calculations.
    """
    _instances = []

    def __init__(self,
                 frag_num: int,
                 max_frag: int,
                 symbols: str,
                 geom: list[list[float]],
                 ra: list[int],
                 par: dict[str, Any],
                 parent: str,
                 equiv: list[list[int]],
                 charge: int | None = None,
                 mult: int | None = None,
                 atoms: Atoms | None = None
                 ) -> None:
        """
        Class generator.
        The fragment is part of an emsemble of fragments
        used to generate pivot points in VRC TST.
        """
        Fragment._instances.append(self)

        self.frag_number: int = frag_num
        self.max_frag: int = max_frag
        self.pivot_points: list[list[float]] = []
        self.par: dict[str, Any] = par
        self.ra: NDArray[int64] = np.array(ra, dtype=int64)
        self.atoms = Atoms(symbols=symbols,
                           positions=geom
                           )
        self.equiv = equiv
        self.orient: dict[str, NDArray[float32]] = {}
        self.parent: str = parent
        if charge is None:
            self.charge: int = sum(self.atoms.get_initial_charges())
        else:
            self.charge: int = charge
        if mult is None:
            self.mult: int = 1
        elif mult is not None:
            self.mult: int = mult
        self.formula: list[str] = self.atoms.get_chemical_symbols()
        self.frag_name: str = self.atoms.get_chemical_formula()
        self.geom = np.array(geom)

        self.recentre()
        self.frag_name: str
        if len(self.geom) == 1:
            self.frag_type = 'Monoatomic'
        elif abs(np.prod(self.atoms.get_moments_of_inertia())) < 1e-6:
            self.frag_type = 'Linear'
        else:
            self.frag_type = 'Nonlinear'

        Fragment.set_fragnames(self)

        super(Fragment, self).__init__(name=self.frag_name,
                                       charge=self.charge,
                                       mult=self.mult,
                                       atom=self.atoms.symbols,
                                       geom=self.geom)
        self.characterize()
        self.com: NDArray[float64] = np.zeros(3)

    def __repr__(self) -> str:
        """Return the string representation of the fragment used in rotdPy.

        Returns:
            str: string to rebuild the fragment in rotdPy input.
        """
        with open(f'{kb_path}/tpl/rotdPy_frag.tpl', 'r') as f:
            tpl: str = f.read()
        rpr: str = tpl.format(
            frag_name=self.frag_name,
            frag_type=self.frag_type,
            formula=self.formula,
            positions=np.round(np.array(self.geom),
                               decimals=4).tolist()
        )
        rpr = rpr.replace('], [', '],\n                     [')

        return rpr

    @classmethod
    def set_fragnames(cls, self) -> None:
        """Create an intuitive name for the fragments to use in rotdPy input.
        """
        cls._fragnames: list[str] = [inst.frag_name for inst in cls._instances[:-1]]
        if f"{self.frag_name}_0" in cls._fragnames:
            self.frag_name = f"{self.frag_name}_{self.frag_number}"
        elif self.frag_name in cls._fragnames:
            index: int = cls._fragnames.index(self.frag_name)
            cls._instances[:-1][index].frag_name = f"{self.frag_name}_{index}"
            cls._fragnames = [inst.frag_name for inst in cls._instances[:-1]]
            self.frag_name = f"{self.frag_name}_{self.frag_number}"
        cls._fragnames.append(self.frag_name)

    @classmethod
    def get_fragnames(cls) -> list[str]:
        return cls._fragnames

    def recentre(self):
        self.atoms = Atoms(symbols=self.formula,
                           positions=self.geom)
        self.geom: NDArray[Any] = np.subtract(self.geom,
                                              self.atoms.get_center_of_mass())

    def get_pp_on_com(self) -> list[float]:
        """Return the coordinate of the center of mass
        of the fragment

        Returns:
            list[float]: COM coordinates (origin)
        """

        return np.round(copy.copy(self.com),
                        decimals=5).tolist()

    def get_pp_on_atom(self, index) -> list[float]:
        """Return the coordinate of the desired atom

        Args:
            index (int): Index of the atom

        Returns:
            list[float]: coordinates
        """
        return np.round(copy.copy(self.geom[index]),
                        decimals=5).tolist()

    def get_pp_next_to_ra(self,
                          index: int,
                          dist_from_ra: float = 0.0
                          )  -> tuple[Any, Any]:
        """Get the atom type and use it to get
        the coordinates of the pivot point.

        Args:
            index (int): index of the reactive atom
            dist_from_ra (float, optional): pivot point distance from
                                            the reactive atom.
                                            Defaults to 0.0.

        Returns:
            _type_: coordinates
        """
        orientation_vect: NDArray[float32] = self.get_pp_orientation(
            index=index)
        coord = (
            orientation_vect *
            dist_from_ra*constants.BOHRtoANGSTROM +
            self.geom[index]).tolist()
        pp_dist = np.full(
            shape=len(coord),
            fill_value=dist_from_ra*constants.BOHRtoANGSTROM
            ).tolist()
        return coord, pp_dist

    def get_pp_orientation(self,
                           index: int) -> NDArray[float32]:
        """Create a list of orientation vectors for all
        pivot points.

        Args:
            index (int): index of atom in fragment's geometry
                         from which the orientation starts.

        Raises:
            NotImplementedError: _description_
            NotImplementedError: _description_

        Returns:
            NDArray[float32]: _description_
        """
        if str(index) in self.orient:
            return self.orient[str(index)]
        orient: NDArray[float32] = np.array([[0, 0, 0]], dtype=float32)
        # if self.par['pp_orient'] == 'geometric':
        # element = self.atom[index]
        # DETECT THE CONNECTIVITY
        nconnect = 0
        ndouble = 0
        ntriple = 0
        for this_bond in self.bond[index]:
            if this_bond == 0:
                continue
            nconnect += this_bond
            if this_bond == 2:
                ndouble += 1
            elif this_bond == 3:
                ntriple += 1
        if nconnect == 1:
            orient = self.pp_aligned_with_bond(index=index)
        elif nconnect == 2:
            if ndouble == 1:
                orient = self.pp_aligned_with_bond(index=index)
            else:
                orient = self.pp_from_homo(index=index)
        elif nconnect == 3:
            if ntriple == 1:
                orient = self.pp_aligned_with_bond(index=index)
            elif ndouble == 1:
                # To adapt in case linear
                orient = self.pp_from_homo(index=index)
            else:
                orient = self.pp_bipyramide_triangle_base(index)
        elif nconnect == 4:
            if ntriple == 1:
                # To adapt in case linear
                orient = self.pp_from_homo(index=index)
            elif ndouble == 2:
                # To adapt in case linear
                orient = self.pp_from_homo(index=index)
            elif ndouble == 1:
                orient = self.pp_bipyramide_triangle_base(index)
            else:
                raise NotImplementedError('Automatic pivot points orientation from connectivity not implemented yet for this configuration.')
        else:
            raise NotImplementedError('Automatic pivot points orientation from connectivity not implemented yet for this configuration.')
        self.orient[str(index)] = orient
        # Find the orientation of the pivot point from the orbital analysis
        # Requires the generation of a gaussian cubefile
        # elif self.par['pp_orient'] == 'homo':
        #     if os.path.isfile(f'{self.parent}/vrctst/{self.chemid}_vts.cube'):
        #         with open(f'{self.parent}/vrctst/{self.chemid}_vts.cube', 'r') as f:
        #             cubefile: list[str] = f.readlines()
        #     else:
        #         raise TypeError('HOMO mode selected for pivot point orientation, but cubefile is not available.')
        #     step: NDArray[float32] = np.zeros(3, dtype=float32)
        #     dim: NDArray[int16] = np.zeros(3, dtype=int16)
        #     # Saves info from the header: step size, origin and dimension
        #     for ln, line in enumerate(cubefile):
        #         if ln == 2:
        #             natom, x, y, z = line.split()[:-1]
        #             origin: NDArray[float32] = np.array([x, y, z],
        #                                                 dtype=float32)
        #             natm = abs(int(natom))
        #             cube_geom: NDArray[float32] = np.empty((natm,3), dtype=float32)
        #             continue
        #         elif ln == 3:
        #             dim[0] = int16(line.split()[0])
        #             step[0] = float32(line.split()[1])
        #             continue
        #         elif ln == 4:
        #             dim[1] = int16(line.split()[0])
        #             step[1] = float32(line.split()[2])
        #             continue
        #         elif ln == 5:
        #             dim[2] = int16(line.split()[0])
        #             step[2] = float32(line.split()[3])
        #             continue
        #         # Save the cube geometry
        #         elif ln > 5 and ln <= 5 + natm:
        #             x, y, z = line.split()[2:]
        #             cube_geom[ln-6] = np.array([x, y, z], dtype=float32)
        #         elif ln > 5 and len(line.split()) == 2:
        #             start: int = ln + 1
        #             break

        #     cube_geom *= BOHRtoANGSTROM
        #     origin *= BOHRtoANGSTROM
        #     step *= BOHRtoANGSTROM

        #     # Create the box

        #     # Size of a z-axis block in number of lines
        #     ysize = int(math.ceil(dim[2]/6))
        #     xsize = int(ysize*dim[1])

        #     # Change the geom for the one in the cube to ensure same orientation
        #     self.geom = cube_geom
        #     self.recentre()

        #     for eq in self.equiv:
        #         if index in eq:
        #             to_integrate: list[int] = eq
        #             break

        #     # max_integral = -1
        #     max_var: float = np.inf
        #     all_orient = []
        #     for atm_order, atm in enumerate(to_integrate):

        #         # define a searchbox for max/min value around atom of interest
        #         searchbox: NDArray[Any] = np.array(
        #             [cube_geom[atm]-1, cube_geom[atm]+1])
        #         search_origin: NDArray[int16] = np.array(np.trunc(
        #             np.absolute(origin - searchbox[0]) / step), dtype=int16)
        #         search_end: NDArray[int16] = np.array(np.trunc(
        #             np.absolute(origin - searchbox[1]) / step)+1, dtype=int16)

        #         box_dim = (
        #                 search_end[0] - search_origin[0],
        #                 search_end[1] - search_origin[1],
        #                 search_end[2] - search_origin[2]
        #                 )

        #         start_line: int = \
        #             start + \
        #             search_origin[0] * xsize + \
        #             search_origin[1] * ysize
        #         stop_line: int = \
        #             start + \
        #             search_end[0] * xsize + \
        #             search_end[1] * ysize + 1

        #         box_origin = search_origin * step + origin

        #         # Defines a sphere in the box
        #         radius = 0.8
        #         sphere = np.empty(box_dim, dtype=bool)
        #         for x in range(box_dim[0]):
        #             for y in range(box_dim[1]):
        #                 for z in range(box_dim[2]):
        #                     coord = np.array([x, y, z]) * step + box_origin
        #                     if np.linalg.norm(coord-cube_geom[atm]) < radius:
        #                         sphere[x, y, z] = True
        #                     else:
        #                         sphere[x, y, z] = False

        #         # Saves the part of interest of the orbital in numpy array
        #         orb_box: NDArray[float32] = np.empty(
        #             shape=box_dim,
        #             dtype=float32
        #         )

        #         # Read data for the box only
        #         save_pos = np.arange(0, box_dim[2])
        #         for ln, line in enumerate(cubefile[start_line:stop_line]):
        #             nx = int(np.trunc((ln+start_line-start) / xsize))
        #             if np.in1d(nx, np.arange(search_origin[0], search_end[0])):
        #                 rest_y = (ln+start_line-start) % xsize
        #                 ny = int(np.trunc(rest_y / ysize))
        #                 if np.in1d(ny, np.arange(search_origin[1], search_end[1])):
        #                     # rest_z goes from 0 to ysize-1
        #                     rest_z = rest_y % ysize
        #                     indexes_in_line: NDArray[Any] = np.arange(rest_z*6, rest_z*6+6)
        #                     is_in_box: NDArray[bool_] = np.in1d(
        #                         indexes_in_line,
        #                         np.arange(search_origin[2], search_end[2]))
        #                     zidx_in_box = (
        #                         indexes_in_line[is_in_box] - search_origin[2]
        #                         ).tolist()
        #                     # Boolean mask
        #                     where2save = np.in1d(save_pos, zidx_in_box)
        #                     which2save = np.where(is_in_box)[0]
        #                     if not any(where2save):
        #                         continue
        #                     orb_box[nx-search_origin[0],
        #                             ny-search_origin[1],
        #                             where2save] = np.take(line.split(),
        #                                                   indices=which2save)

        #         # find indexes of max value in sphere
        #         orb_max: float32 = np.max(orb_box[sphere])

        #         # Search indexes of the value in the box
        #         tmp_max_idx: NDArray[Any] = np.transpose(np.where(orb_box == orb_max))
        #         # Take closest value in case several values are minimum
        #         if len(tmp_max_idx) > 1:
        #             dist: list[floating[Any]] = [
        #                 np.linalg.norm(i*step+origin-cube_geom[atm])
        #                 for i in tmp_max_idx]
        #             max_idx: NDArray[Any] = tmp_max_idx[np.where(dist == np.min(dist))][0]
        #         else:
        #             max_idx: NDArray[Any] = tmp_max_idx[0]

        #         max_coord = step * max_idx + box_origin

        #         max_orient: NDArray[float32] = np.array(
        #             geometry.unit_vector(max_coord - cube_geom[atm]),
        #             dtype=float32)

        #         # Define a sphere around highest density in the box
        #         max_radius = np.linalg.norm(max_coord - cube_geom[atm])
        #         min_radius = 0.8 * max_radius
        #         max_sphere = np.empty(box_dim, dtype=bool)
        #         for x in range(box_dim[0]):
        #             for y in range(box_dim[1]):
        #                 for z in range(box_dim[2]):
        #                     coord = np.array([x, y, z]) * step + box_origin
        #                     if np.linalg.norm(coord-max_coord) < max_radius and\
        #                     np.linalg.norm(coord-max_coord) > min_radius:
        #                         max_sphere[x, y, z] = True
        #                     else:
        #                         max_sphere[x, y, z] = False

        #         # Variance in the valence of the sphere, centered on the electron             
        #         tot_var = np.var(orb_box[max_sphere])

        #         # Check if minimum should be a pivot point
        #         orb_min: float32 = np.min(orb_box[sphere])
        #         if abs(orb_min) > 0.8*orb_max:
        #             # Search indexes of the value in the box
        #             tmp_min_idx = np.transpose(np.where(orb_box == orb_min))
        #             # Take closest value in case several values are minimum
        #             if len(tmp_min_idx) > 1:
        #                 dist = [np.linalg.norm(i*step+origin-cube_geom[atm]) for i in tmp_min_idx]
        #                 min_idx = tmp_min_idx[np.where(dist == np.min(dist))][0]
        #             else:
        #                 min_idx = tmp_min_idx[0]

        #             min_coord = step * min_idx + box_origin

        #             min_orient: NDArray[float32] = np.array(
        #                 geometry.unit_vector(min_coord - cube_geom[atm]),
        #                 dtype=float32)

        #             all_orient.append([max_orient, min_orient])

        #             # Define a sphere around lowest density in the box
        #             max_radius = np.linalg.norm(min_coord - cube_geom[atm])
        #             min_radius = 0.8 * max_radius
        #             min_sphere = np.empty(box_dim, dtype=bool)
        #             for x in range(box_dim[0]):
        #                 for y in range(box_dim[1]):
        #                     for z in range(box_dim[2]):
        #                         coord = np.array([x, y, z]) * step + box_origin
        #                         if np.linalg.norm(coord-min_coord) < max_radius and\
        #                            np.linalg.norm(coord-min_coord) > min_radius:
        #                             min_sphere[x, y, z] = True
        #                         else:
        #                             min_sphere[x, y, z] = False
        #             # Variance in the valence of the sphere, centered on the electron
        #             tot_var += np.var(orb_box[min_sphere])
        #         else:
        #             all_orient.append([max_orient])
        #         logger.info(f"Total orbital variance of {tot_var} around atom {atm} in fragment {self.chemid}")
        #         # A lowest variance means more spherical symmetry
        #         if tot_var < max_var:
        #             max_var = tot_var
        #             best_atm = to_integrate[0]

        #     # Check if min and max of the orbital are aligned
        #     # The cutoff angle is defined as arccos(1-x),
        #     # where x is the float at the end of this line.
        #     # 0.2 corresponds to approx. 37 degree
        #     logger.info(f"Model atom selected from variance for fragment {self.chemid}: {best_atm}")
        #     max_orient = all_orient[to_integrate.index(best_atm)][0]
        #     min_orient = all_orient[to_integrate.index(best_atm)][1]
        #     if np.linalg.norm(
        #         np.subtract(min_orient, max_orient)) > 0.2:
        #         orient = np.array([max_orient, min_orient], dtype=float32)
        #     else:
        #         orient = np.array([max_orient], dtype=float32)

        # # Reorientate the chosen orientation depending on the neighbours orientation
        # best_grp: list[bool] = []
        # best_grp_atmids = {}
        # # Detect order of atomids for the rotation
        # for num, bond in enumerate(self.bond[best_atm]):
        #     if bond:
        #         if self.atomid[num] not in best_grp_atmids:
        #             best_grp_atmids[self.atomid[num]] = [[int(np.sum(best_grp)+1)], 0]
        #         else:
        #             best_grp_atmids[self.atomid[num]][0].append(np.sum(best_grp)+1)
        #         best_grp.append(True)
        #     else:
        #         best_grp.append(False)
        # grp1_centered = np.empty((int(np.sum(best_grp))+1, 3), dtype=float32)
        # grp1_centered[0] = cube_geom[best_atm] - cube_geom[best_atm]
        # grp1_centered[1:] = cube_geom[best_grp] - cube_geom[best_atm]
        # for atm in to_integrate:
        #     for atom_id in best_grp_atmids:
        #         best_grp_atmids[atom_id][1] = 0
        #     grp2_centered = np.empty((np.sum(best_grp)+1,3), dtype=float32)
        #     grp2_centered[0] = cube_geom[atm] - cube_geom[atm]
        #     for num, bond in enumerate(self.bond[atm]):
        #         if bond:
        #             same_atomid, order_index = best_grp_atmids[self.atomid[num]]
        #             grp2_centered[same_atomid[order_index]] = cube_geom[num] - cube_geom[atm]
        #             best_grp_atmids[self.atomid[num]][1] += 1
        #     U = rmsd.kabsch(grp1_centered, grp2_centered)
        #     self.orient[str(atm)] = np.dot(orient, U)
        return orient

    def pp_aligned_with_bond(self,
                             index: int,
                             angle: float = 0.0,
                             last_neighbor: int | None = None
                             ) -> NDArray[float32]:
        """Create pivot point aligned with the bond

        Args:
            index (int): Index of atom in fragment.
            angle (float, optional): Angle deviation (degree) in the plane
                                     bond+COM or bond+last_neighbor.
                                     Defaults to 0.0.
            last_neighbor (int, optional): Index for last point to define
                                           the plane in which
                                           to move the pivot point.
                                           Defaults to None.

        Returns:
            list[float]: coordinate of the pivot point.
        """
        ra_pos: NDArray[Any] = np.array(self.geom[index], dtype=float)
        neighbor_pos: list[list[float]] = []
        for neighbor_index, this_bond in enumerate(self.bond[index]):
            if this_bond != 0:
                neighbor_pos.append(self.geom[neighbor_index])
                break

        if angle == 0.0 or not isinstance(angle, (float, int)) \
           or last_neighbor is None or not isinstance(last_neighbor, int):
            try:
                pp_orient: NDArray[Any] = np.subtract(ra_pos, neighbor_pos[0])
            except NameError:
                logger.warning("Could not find any bond for atom {}.\
                                Setting it to COM".format(self.atom[index]))
                pp_orient = np.subtract(ra_pos, self.com) + 0.0000000001
        else:
            neighbor_pos.append(self.geom[last_neighbor])
            v1: NDArray[Any] = np.subtract(neighbor_pos[0], ra_pos)
            v2: NDArray[Any] = np.subtract(neighbor_pos[1], ra_pos)
            axis: NDArray[floating[Any]] = geometry.unit_vector(
                vector=np.cross(v2, v1))
            pp_orient = np.dot(
                geometry.rotation_matrix(axis=axis, theta=math.radians(angle)),
                v1)
        pp_orient = geometry.unit_vector(pp_orient)

        return np.array([pp_orient], dtype=float32)

    def pp_from_homo(self,
                     index: int
                     ) -> NDArray[float32]:
        if os.path.isfile(f'{self.parent}/vrctst/{self.chemid}_vts.cube'):
            with open(f'{self.parent}/vrctst/{self.chemid}_vts.cube', 'r') as f:
                cubefile: list[str] = f.readlines()
        else:
            raise TypeError('HOMO mode selected for pivot point orientation, but cubefile is not available.')
        step: NDArray[float32] = np.zeros(3, dtype=float32)
        dim: NDArray[int16] = np.zeros(3, dtype=int16)
        elements: list[str] = []
        num2str: dict[str, str] = {
            '1': 'H',
            '6': 'C',
            '7': 'N',
            '8': 'O',
            '16': 'S'}
        # Saves info from the header: step size, origin and dimension
        for ln, line in enumerate(cubefile):
            if ln == 2:
                natom, x, y, z = line.split()[:-1]
                origin: NDArray[float32] = np.array([x, y, z],
                                                    dtype=float32)
                natm = abs(int(natom))
                cube_geom: NDArray[float32] = np.empty((natm,3), dtype=float32)
                continue
            elif ln == 3:
                dim[0] = int16(line.split()[0])
                step[0] = float32(line.split()[1])
                continue
            elif ln == 4:
                dim[1] = int16(line.split()[0])
                step[1] = float32(line.split()[2])
                continue
            elif ln == 5:
                dim[2] = int16(line.split()[0])
                step[2] = float32(line.split()[3])
                continue
            # Save the cube geometry
            elif ln > 5 and ln <= 5 + natm:
                num, _, x, y, z = line.split()
                cube_geom[ln-6] = np.array([x, y, z], dtype=float32)
                elements.append(num2str[num])
            elif ln > 5 and len(line.split()) == 2:
                start: int = ln + 1
                break

        cube_geom *= BOHRtoANGSTROM
        origin *= BOHRtoANGSTROM
        step *= BOHRtoANGSTROM

        # Create the box

        # Size of a z-axis block in number of lines
        ysize = int(math.ceil(dim[2]/6))
        xsize = int(ysize*dim[1])

        # Change the geom for the one in the cube to ensure same orientation
        tmp_atm = Atoms(elements, positions=cube_geom)
        tmp_stp = StationaryPoint.from_ase_atoms(tmp_atm)
        tmp_stp.characterize()
        reorder_coord(self, tmp_stp)
        self.geom = tmp_stp.geom
        self.recentre()

        # for eq in self.equiv:
        #     if index in eq:
        #         to_integrate: list[int] = eq
        #         break

        # max_integral = -1
        # max_var: float = np.inf
        # all_orient = []
        

        # define a searchbox for max/min value around atom of interest
        searchbox: NDArray[Any] = np.array(
            [cube_geom[index]-1, cube_geom[index]+1])
        search_origin: NDArray[int16] = np.array(np.trunc(
            np.absolute(origin - searchbox[0]) / step), dtype=int16)
        search_end: NDArray[int16] = np.array(np.trunc(
            np.absolute(origin - searchbox[1]) / step)+1, dtype=int16)

        box_dim = (
                search_end[0] - search_origin[0],
                search_end[1] - search_origin[1],
                search_end[2] - search_origin[2]
                )

        start_line: int = \
            start + \
            search_origin[0] * xsize + \
            search_origin[1] * ysize
        stop_line: int = \
            start + \
            search_end[0] * xsize + \
            search_end[1] * ysize + 1

        box_origin = search_origin * step + origin

        # Defines a sphere in the box
        radius = 0.8
        sphere = np.empty(box_dim, dtype=bool)
        for x in range(box_dim[0]):
            for y in range(box_dim[1]):
                for z in range(box_dim[2]):
                    coord = np.array([x, y, z]) * step + box_origin
                    if np.linalg.norm(coord-cube_geom[index]) < radius:
                        sphere[x, y, z] = True
                    else:
                        sphere[x, y, z] = False

        # Saves the part of interest of the orbital in numpy array
        orb_box: NDArray[float32] = np.empty(
            shape=box_dim,
            dtype=float32
        )

        # Read data for the box only
        save_pos = np.arange(0, box_dim[2])
        for ln, line in enumerate(cubefile[start_line:stop_line]):
            nx = int(np.trunc((ln+start_line-start) / xsize))
            if np.in1d(nx, np.arange(search_origin[0], search_end[0])):
                rest_y = (ln+start_line-start) % xsize
                ny = int(np.trunc(rest_y / ysize))
                if np.in1d(ny, np.arange(search_origin[1], search_end[1])):
                    # rest_z goes from 0 to ysize-1
                    rest_z = rest_y % ysize
                    indexes_in_line: NDArray[Any] = np.arange(rest_z*6, rest_z*6+6)
                    is_in_box: NDArray[bool_] = np.in1d(
                        indexes_in_line,
                        np.arange(search_origin[2], search_end[2]))
                    zidx_in_box = (
                        indexes_in_line[is_in_box] - search_origin[2]
                        ).tolist()
                    # Boolean mask
                    where2save = np.in1d(save_pos, zidx_in_box)
                    which2save = np.where(is_in_box)[0]
                    if not any(where2save):
                        continue
                    orb_box[nx-search_origin[0],
                            ny-search_origin[1],
                            where2save] = np.take(line.split(),
                                                    indices=which2save)

        # find indexes of max value in sphere
        orb_max: float32 = np.max(orb_box[sphere])

        # Search indexes of the value in the box
        tmp_max_idx: NDArray[Any] = np.transpose(np.where(orb_box == orb_max))
        # Take closest value in case several values are minimum
        if len(tmp_max_idx) > 1:
            dist: list[floating[Any]] = [
                np.linalg.norm(i*step+origin-cube_geom[index])
                for i in tmp_max_idx]
            max_idx: NDArray[Any] = tmp_max_idx[np.where(dist == np.min(dist))][0]
        else:
            max_idx: NDArray[Any] = tmp_max_idx[0]

        max_coord = step * max_idx + box_origin

        max_orient: NDArray[float32] = np.array(
            geometry.unit_vector(max_coord - cube_geom[index]),
            dtype=float32)

        # # Define a sphere around highest density in the box
        # max_radius = np.linalg.norm(max_coord - cube_geom[index])
        # min_radius = 0.8 * max_radius
        # max_sphere = np.empty(box_dim, dtype=bool)
        # for x in range(box_dim[0]):
        #     for y in range(box_dim[1]):
        #         for z in range(box_dim[2]):
        #             coord = np.array([x, y, z]) * step + box_origin
        #             if np.linalg.norm(coord-max_coord) < max_radius and\
        #             np.linalg.norm(coord-max_coord) > min_radius:
        #                 max_sphere[x, y, z] = True
        #             else:
        #                 max_sphere[x, y, z] = False

        # # Variance in the valence of the sphere, centered on the electron             
        # tot_var = np.var(orb_box[max_sphere])

        # Check if minimum should be a pivot point
        orb_min: float32 = np.min(orb_box[sphere])
        if abs(orb_min) > 0.8*orb_max:
            # Search indexes of the value in the box
            tmp_min_idx = np.transpose(np.where(orb_box == orb_min))
            # Take closest value in case several values are minimum
            if len(tmp_min_idx) > 1:
                dist = [np.linalg.norm(i*step+origin-cube_geom[index]) for i in tmp_min_idx]
                min_idx = tmp_min_idx[np.where(dist == np.min(dist))][0]
            else:
                min_idx = tmp_min_idx[0]

            min_coord = step * min_idx + box_origin

            min_orient: NDArray[float32] = np.array(
                geometry.unit_vector(min_coord - cube_geom[index]),
                dtype=float32)

            # all_orient.append([max_orient, min_orient])

            # # Define a sphere around lowest density in the box
            # max_radius = np.linalg.norm(min_coord - cube_geom[index])
            # min_radius = 0.8 * max_radius
            # min_sphere = np.empty(box_dim, dtype=bool)
            # for x in range(box_dim[0]):
            #     for y in range(box_dim[1]):
            #         for z in range(box_dim[2]):
            #             coord = np.array([x, y, z]) * step + box_origin
            #             if np.linalg.norm(coord-min_coord) < max_radius and\
            #                 np.linalg.norm(coord-min_coord) > min_radius:
            #                 min_sphere[x, y, z] = True
            #             else:
            #                 min_sphere[x, y, z] = False
            # # Variance in the valence of the sphere, centered on the electron
            # tot_var += np.var(orb_box[min_sphere])
        # else:
        #     all_orient.append([max_orient])
        # logger.info(f"Total orbital variance of {tot_var} around atom {index} in fragment {self.chemid}")
        # A lowest variance means more spherical symmetry
        # if tot_var < max_var:
        #     max_var = tot_var
        #     best_atm = to_integrate[0]

        # Check if min and max of the orbital are aligned
        # The cutoff angle is defined as arccos(1-x),
        # where x is the float at the end of this line.
        # 0.2 corresponds to approx. 37 degree
        # logger.info(f"Model atom selected from variance for fragment {self.chemid}: {best_atm}")
        # max_orient = all_orient[to_integrate.index(best_atm)][0]
        # min_orient = all_orient[to_integrate.index(best_atm)][1]
            if np.linalg.norm(
                np.subtract(min_orient, max_orient)) > 0.2:
                orient = np.array([max_orient, min_orient], dtype=float32)
            else:
                orient = np.array([max_orient], dtype=float32)
        else:
            orient = np.array([max_orient], dtype=float32)
        
        return orient


    def pp_triangle(self,
                    index: int,
                    angle: float | None = None) -> NDArray[float32]:
        """Create pivot point in the middle of the big angle between the 2 neighbors.

        Args:
            index (int): for an angle ABC, index of atom B in fragment.
            angle (float, optional): Deviation angle.
                                     Defaults to None, which leads to
                                     in between the two attached atoms.

        Returns:
            pp_coord (list[float]): 3D coordinate of the pivot point.
        """
        ra_pos: NDArray[Any] = np.array(self.geom[index], dtype=float)
        neighbor_pos: list[NDArray[Any]] = []
        for ni, bond in enumerate(self.bond[index]):
            if bond != 0:
                neighbor_pos.append(np.array(self.geom[ni], dtype=float))
        v1: NDArray[Any] = np.subtract(neighbor_pos[0], ra_pos)
        v2: NDArray[Any] = np.subtract(neighbor_pos[1], ra_pos)
        small_angle: float = geometry.calc_angle(
            a=neighbor_pos[0],
            b=ra_pos,
            c=neighbor_pos[1])
        big_angle: float = 2*pi - small_angle
        if angle is None or not isinstance(angle, (float, int)):
            angle: float = big_angle/2
        axis: NDArray[floating[Any]] = geometry.unit_vector(
            vector=np.cross(a=v2,
                            b=v1))
        pp_orient: NDArray[Any] = np.dot(
            a=geometry.rotation_matrix(axis, angle),
            b=v1)

        return np.array(geometry.unit_vector(vector=pp_orient), dtype=float32)

    def pp_bipyramide_triangle_base(self,
                                    index: int) -> NDArray[float32]:
        """Create one or two pivot points on one or both side of the plane.

        Args:
            index (int): atom having 3 neighbours

        Returns:
            _type_: _description_
        """
        n_pp = 2
        ra_pos: NDArray[Any] = np.array(self.geom[index], dtype=float)
        neighbor_pos = []
        for neighbor_index, this_bond in enumerate(
         np.array(self.bonds)[0, index]):
            if this_bond != 0:
                neighbor_pos.append(np.array(self.geom[neighbor_index],
                                             dtype=float))
        v1: NDArray[Any] = np.subtract(neighbor_pos[0], ra_pos)
        v2: NDArray[Any] = np.subtract(neighbor_pos[1], ra_pos)
        v3: NDArray[Any] = np.subtract(neighbor_pos[2], ra_pos)
        plane: tuple[NDArray[Any], float] = geometry.plane_from_points(
            v0=v1,
            v1=v2,
            v2=v3)
        ra_to_plane: float = geometry.dist_point_to_plane(
            point=ra_pos,
            plane=plane)
        if abs(ra_to_plane) >= .1:
            # If carbon atom out of plane,
            # only place a single pp on other side of the plane
            n_pp = 1
        pp_list = []
        # To know in which direction to place the pivot point
        plane_direction = geometry.unit_vector(np.dot(plane[0], ra_pos))
        for i in range(n_pp):
            pp_orient = np.array(
                geometry.unit_vector(vector=plane[0]) *
                plane_direction*np.power(-1, i),
                dtype=float)
            pp_list.append(geometry.unit_vector(pp_orient))
            
        return np.array(pp_list, dtype=float32)
