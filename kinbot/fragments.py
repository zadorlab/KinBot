from typing import Any
from kinbot import kb_path
from kinbot.stationary_pt import StationaryPoint
from ase.atoms import Atoms
import numpy as np
from numpy.typing import NDArray
from numpy import float64, int64, pi, floating
from kinbot import pp_tables
from kinbot import geometry
import logging
import copy
import math

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
        self.atom = Atoms(symbols=symbols,
                          positions=geom
                          )
        if charge is None:
            self.charge: int = sum(self.atom.get_initial_charges())
        else:
            self.charge: int = charge
        if mult is None and\
           self.atom.calc is None or\
           'mult' not in self.atom.calc.parameters:
            self.mult: int = 1
        elif mult is not None:
            self.mult: int = mult
        else:
            self.mult: int = self.atom.calc.parameters['mult']
        self.formula: str = str(self.atom.get_chemical_symbols())
        self.frag_name: str = self.atom.get_chemical_formula()
        self.geom: NDArray[Any] = np.subtract(geom,
                                              self.atom.get_center_of_mass())
        self.atom = Atoms(symbols=symbols,
                          positions=geom
                          )

        self.frag_name: str
        Fragment.set_fragnames(self)

        super(Fragment, self).__init__(name=self.frag_name,
                                       charge=self.charge,
                                       mult=self.mult,
                                       atom=self.atom.symbols,
                                       geom=self.geom)
        self.characterize()
        self.com: NDArray[float64] = np.zeros(3)

    def __repr__(self) -> str:
        """Return the string representation of the fragment used in rotdPy.

        Returns:
            str: string to rebuild the fragment in rotdPy input.
        """
        # TODO: Modify nonlinear by a variable that detects linearity
        with open(f'{kb_path}/tpl/rotdPy_frag.tpl', 'r') as f:
            tpl: str = f.read()
        rpr: str = tpl.format(
            frag_name=self.frag_name,
            frag_type='Nonlinear',
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

    def get_chemical_formula(self) -> str:
        """Create the string of element in the order of the coordinates

        Returns:
            str: chemical elements
        """
        all_elem: str = ""
        for elem in self.atom:
            all_elem += f"{elem}"
        return all_elem

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
                          ) -> list[list[float]]:
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
        atom_type: str = self.get_atom_type(index)
        coord: list[list[float]] = self.get_pp_coord(
            index=index,
            atom_type=atom_type,
            dist_from_ra=dist_from_ra)

        return np.round(np.asarray(coord),
                        decimals=5).tolist()

    def get_atom_type(self,
                      index: int) -> str:
        """Detect the connectivity of the atom index to
        get its atom type.

        Args:
            index (int): index of atom in fragment

        Returns:
            atom_type (str): atom type name that depends on connectivity.
        """
        element = self.atom[index]
        nconnect = 0
        ndouble = 0
        ntriple = 0
        for this_bond in self.bond[index]:
            if this_bond == 0:
                continue
            nconnect += this_bond
            match this_bond:
                case 2:
                    ndouble += 1
                case 3:
                    ntriple += 1
        atom_type: str = pp_tables.atom_type_table(element=element,
                                                   nconnect=nconnect,
                                                   ndouble=ndouble,
                                                   ntriple=ntriple)
        return atom_type

    def get_pp_coord(self,
                     index: int,
                     atom_type: str,
                     dist_from_ra: float) -> list[list[float]]:
        """Call the appropriate method to create
        the pivot point depending on the atom type.

        Args:
            index (int): index of the atom holding the pivot point
            atom_type (str): string describing the geometry of the
                             environment around the atom
            dist_from_ra (float): length of the pivot point

        Returns:
            list[list[float]]: list of 3D coordinates.
        """
        match atom_type:
            case 'H' | 'C' | 'O' | 'S':
                # Create pivot point on atom
                return [self.get_pp_on_atom(index=index)]
            case 'H_lin':
                pp_coord: list[list[float]] = \
                    [self.create_pp_aligned_with_bond(
                        index=index,
                        length=dist_from_ra)]
                return pp_coord
            case 'C_lin':
                pp_coord: list[list[float]] = \
                    [self.create_pp_aligned_with_bond(
                        index=index,
                        length=dist_from_ra)]
                return pp_coord
            case 'C_tri':
                pp_coord: list[list[float]] = \
                    [self.create_pp_triangle(
                        index=index,
                        length=dist_from_ra)]
                return pp_coord
            case 'C_quad':
                pp_list: list[list[float]] = \
                    self.create_pp_bipyramide_triangle_base(
                        index=index,
                        length=dist_from_ra)
                return pp_list
            case 'N_tri':
                pp_coord: list[list[float]] = \
                    [self.create_pp_aligned_with_bond(
                        index=index,
                        length=dist_from_ra)]
                return pp_coord
            case 'N_pyr':
                pp_coord: list[list[float]] = \
                    [self.create_pp_triangle(
                        index=index,
                        length=dist_from_ra)]
                return pp_coord
            case 'N_quad':
                pp_list: list[list[float]] = \
                    self.create_pp_bipyramide_triangle_base(
                        index=index,
                        length=dist_from_ra)
                return pp_list
            case 'O_tri':
                pp_coord: list[list[float]] = \
                    [self.create_pp_aligned_with_bond(
                        index=index,
                        length=dist_from_ra)]
                return pp_coord
            case 'O_quad':
                pp_coord: list[list[float]] = \
                    [self.create_pp_triangle(
                        index=index,
                        length=dist_from_ra)]
                return pp_coord
            case 'S_tri':
                pp_coord: list[list[float]] = \
                    [self.create_pp_aligned_with_bond(
                        index=index,
                        length=dist_from_ra)]
                return pp_coord
            case 'S_pyr':
                pp_coord: list[list[float]] = \
                    [self.create_pp_triangle(
                        index=index,
                        length=dist_from_ra)]
                return pp_coord
            case 'S_lin':
                pp_coord: list[list[float]] = \
                    [self.create_pp_aligned_with_bond(
                        index=index,
                        length=dist_from_ra)]
                return pp_coord
            case 'S_bip_tri_l':
                pp_list: list[list[float]] = \
                    self.create_pp_bipyramide_triangle_base(
                        index=index,
                        length=dist_from_ra)
                return pp_list
            case 'S_bip_tri':
                return [[]]
            case 'S_quad':
                pp_list: list[list[float]] = \
                    self.create_pp_bipyramide_triangle_base(
                        index=index,
                        length=dist_from_ra)
                return pp_list
            case 'S_bip_quad_t':
                return [[]]
            case _:
                return [[]]

    def create_pp_aligned_with_bond(self,
                                    index: int,
                                    length: float,
                                    angle: float = 0.0,
                                    last_neighbor: int | None = None
                                    ) -> list[float]:
        """Create pivot point aligned with the bond

        Args:
            index (int): Index of atom in fragment.
            length (float): pivot point length.
            angle (float, optional): Angle deviation (degree) in the plane
                                     bond+COM or bond+last_neighborg.
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

        # Multiply unit vector with correct orientation with desired length
        pp_vect = np.asarray(length) * geometry.unit_vector(pp_orient)
        # Place pivots point in molecular frame
        pp_coord: NDArray[Any] = np.add(ra_pos, pp_vect).tolist()
        return pp_coord.tolist()

    def create_pp_triangle(self,
                           index: int,
                           length: float,
                           angle: float | None = None) -> list[float]:
        """Create pivot point in the middle of the big angle between the 2 neighbors.

        Args:
            index (int): for an angle ABC, index of atom B in fragment.
            length (float, optional): Length of pivot point. Defaults to None.
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
        pp_vect: NDArray[Any] = length * geometry.unit_vector(vector=pp_orient)
        pp_coord: NDArray[Any] = np.add(ra_pos, pp_vect)
        return pp_coord.tolist()

    def create_pp_bipyramide_triangle_base(self,
                                           index: int,
                                           length: float) -> list[list[float]]:
        """Create one or two pivot points on one or both side of the plane.

        Args:
            index (int): _description_
            length (float): _description_

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
            pp_vect = length * geometry.unit_vector(pp_orient)
            pp_coord: NDArray[Any] = np.add(ra_pos, pp_vect)
            pp_list.append(pp_coord.tolist())
        return pp_list
