from typing import Any
from kinbot import kb_path
from kinbot.stationary_pt import StationaryPoint
from ase.atoms import Atoms
import numpy as np
from numpy import ndarray
from numpy import pi
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
                 ):
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
        self.ra: ndarray = np.array(ra, dtype=int)
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
        self.geom: ndarray[Any] = np.subtract(geom,
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
        self.com: ndarray[float] = np.zeros(3)

    def __repr__(self) -> str:
        #TODO: Modify nonlinear by a variable that detects linearity
        with open(f'{kb_path}/tpl/rotdPy_frag.tpl', 'r') as f:
            tpl = f.read()
        rpr: str = tpl.format(
            frag_name=self.frag_name,
            frag_type='Nonlinear',
            formula=self.formula,
            positions=np.round(np.array(self.geom),
                               decimals=4).tolist()
        )
        rpr = rpr.replace('], [',
                          '],\n                     [')

        return rpr

    @classmethod
    def set_fragnames(cls, self) -> None:
        cls._fragnames = [inst.frag_name for inst in cls._instances[:-1]]
        if f"{self.frag_name}_0" in cls._fragnames:
            self.frag_name = f"{self.frag_name}_{self.frag_number}"
        elif self.frag_name in cls._fragnames:
            index = cls._fragnames.index(self.frag_name)
            cls._instances[:-1][index].frag_name = f"{self.frag_name}_{index}"
            cls._fragnames = [inst.frag_name for inst in cls._instances[:-1]]
            self.frag_name = f"{self.frag_name}_{self.frag_number}"
        cls._fragnames.append(self.frag_name)

    @classmethod
    def get_fragnames(cls) -> list[str]:
        return cls._fragnames

    def get_chemical_formula(self) -> str:
        all_elem = ""
        for elem in self.atom:
            all_elem += f"{elem}"
        return all_elem

    def get_pp_on_com(self):
            return np.round(copy.copy(self.com),\
                            decimals=5).tolist()

    def get_pp_on_atom(self, index):
            return np.round(copy.copy(self.geom[index]),\
                            decimals=5).tolist()

    def get_pp_next_to_ra(self, index, dist_from_ra=0.0):
        atom_type = self.get_atom_type(index)
        coord = self.get_pp_coord(index,\
                                  atom_type,\
                                  dist_from_ra=dist_from_ra)
        return np.round(coord, decimals=5).tolist()

    def get_atom_type(self, index):
        element = self.atom[index]
        nconnect = 0
        ndouble = 0
        ntriple = 0
        for this_bond in np.array(self.bond)[ index]:
            if this_bond != 0:
                nconnect += 1
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

    def get_pp_coord(self, index, atom_type, dist_from_ra=None):
        if dist_from_ra is None:
            # Return a list of relevant distances to try
            dist_from_ra = pp_tables.pp_length_table(atom_type[0])
        match atom_type:
            case 'H' | 'C' | 'O' | 'S':
                # Create pivot point on atom
                return [self.get_pp_on_atom(index)]
            case 'H_lin':
                pp_coord = self.create_pp_aligned_with_bond(index, length=dist_from_ra)
                return pp_coord
            case 'C_lin':
                pp_coord = self.create_pp_aligned_with_bond(index, length=dist_from_ra)
                return pp_coord
            case 'C_tri':
                pp_coord = self.create_pp_triangle(index, length=dist_from_ra)
                return pp_coord
            case 'C_quad':
                pp_list = self.create_pp_bipyramide_triangle_base(index, length=dist_from_ra)
                return pp_list
            case 'N_tri':
                pp_coord = self.create_pp_aligned_with_bond(index, length=dist_from_ra)
                return pp_coord
            case 'N_pyr':
                pp_coord = self.create_pp_triangle(index, length=dist_from_ra)
                return pp_coord
            case 'N_quad':
                pp_list = self.create_pp_bipyramide_triangle_base(index, length=dist_from_ra)
                return pp_list
            case 'O_tri':
                pp_coord = self.create_pp_aligned_with_bond(index, length=dist_from_ra)
                return pp_coord
            case 'O_quad':
                pp_coord = self.create_pp_triangle(index, length=dist_from_ra)
                return pp_coord
            case 'S_tri':
                pp_coord = self.create_pp_aligned_with_bond(index, length=dist_from_ra)
                return pp_coord
            case 'S_pyr':
                pp_coord = self.create_pp_triangle(index, length=dist_from_ra)
                return pp_coord
            case 'S_lin':
                pp_coord = self.create_pp_aligned_with_bond(index, length=dist_from_ra)
                return pp_coord
            case 'S_bip_tri_l':
                pp_list = self.create_pp_bipyramide_triangle_base(index, length=dist_from_ra)
                return pp_list
            case 'S_bip_tri':
                pass
            case 'S_quad':
                pp_list = self.create_pp_bipyramide_triangle_base(index, length=dist_from_ra)
                return pp_list
            case 'S_bip_quad_t':
                pass

    def create_pp_aligned_with_bond(self,
                                    index,
                                    length=None,
                                    angle=None,
                                    last_neighbours=None):
        #Create pivot point aligned with the bond
        ra_pos = np.array(self.geom[index], dtype=float)
        neighbour_pos = []
        for neighbour_index, this_bond in enumerate(self.bond[index]):
            if this_bond != 0:
                neighbour_pos.append(self.geom[neighbour_index])
                break

        if angle is None or not isinstance(angle, (float, int)) \
        or last_neighbours is None or not isinstance(last_neighbours, int):
            try:
                pp_orient = np.subtract(ra_pos, neighbour_pos[0])
            except NameError:
                logger.warning(f"Could not find any bond for atom {self.atom[index]}. Setting it to COM")
                neighbour_pos = self.com
                pp_orient = np.subtract(ra_pos, neighbour_pos) + 0.0000000001
        else:
            neighbour_pos.append(self.geom[last_neighbours])
            v1 = np.subtract(neighbour_pos[0], ra_pos)
            v2 = np.subtract(neighbour_pos[1], ra_pos)
            axis = geometry.unit_vector(np.cross(v2, v1))
            pp_orient = np.dot(geometry.rotation_matrix(axis, math.radians(angle)), v1)
            
        #Multiply unit vector with correct orientation with desired length
        pp_vect = np.asarray(length) * geometry.unit_vector(pp_orient)
        #Place pivots point in molecular frame
        pp_coord = np.add(ra_pos, pp_vect)
        return pp_coord

    def create_pp_triangle(self, index, length=None, angle=None):
        #Create pivot point in the middle of the big angle between the 2 neighbours
        ra_pos = np.array(self.geom[index], dtype=float)
        neighbour_pos = []
        for neighbour_index, this_bond in enumerate(self.bond[index]):
            if this_bond != 0:
                neighbour_pos.append(np.array(self.geom[neighbour_index], dtype=float))
        v1 = np.subtract(neighbour_pos[0], ra_pos)
        v2 = np.subtract(neighbour_pos[1], ra_pos)
        small_angle = geometry.calc_angle(neighbour_pos[0], ra_pos, neighbour_pos[1])
        big_angle = 2*pi - small_angle
        if angle is None or not isinstance(angle, (float, int)):
            angle = big_angle/2
        axis = geometry.unit_vector(np.cross(v2, v1))
        pp_orient = np.dot(geometry.rotation_matrix(axis, angle), v1)
        if length is None:
            length = pp_tables.pp_length_table(self.atom[index], par=self.par)
        pp_vect = length * geometry.unit_vector(pp_orient)
        pp_coord = np.add(ra_pos, pp_vect)
        return pp_coord

    def create_pp_bipyramide_triangle_base(self, index, length=None):
        #Create one or two pivot points on one or both side of the plane
        n_pp = 2
        ra_pos = np.array(self.geom[index], dtype=float)
        neighbour_pos = []
        for neighbour_index, this_bond in enumerate(np.array(self.bonds)[0,index]):
            if this_bond != 0:
                neighbour_pos.append(np.array(self.geom[neighbour_index], dtype=float))
        v1 = np.subtract(neighbour_pos[0], ra_pos)
        v2 = np.subtract(neighbour_pos[1], ra_pos)
        v3 = np.subtract(neighbour_pos[2], ra_pos)
        plane = geometry.plane_from_points(v1, v2, v3)
        ra_to_plane = geometry.dist_point_to_plane(ra_pos, plane)
        if abs(ra_to_plane) >= .1:
            #If carbon atom out of plane, only place a single pp on other side of the plane
            n_pp = 1
        pp_list = []
        #To know in which direction to place the pivot point
        plane_direction = geometry.unit_vector(np.dot(plane[0], ra_pos))
        if length is None:
            length = pp_tables.pp_length_table(self.atom[index],par=self.par)
        for i in range(n_pp):
            pp_orient = np.array(geometry.unit_vector(plane[0])*plane_direction*np.power(-1,i), dtype=float)
            try:
                pp_vect = length * geometry.unit_vector(pp_orient)
            except NameError:
                logger.warning(f"Length of pivot point not defined yet for atom {self.atom[index]}. Setting it to 0.5A.")
                length = 0.5
                pp_vect = length * geometry.unit_vector(pp_orient)
            pp_coord = np.add(ra_pos, pp_vect)
            pp_list.append(pp_coord)
        return pp_list
