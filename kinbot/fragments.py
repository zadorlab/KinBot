from kinbot.stationary_pt import StationaryPoint
import numpy as np
from numpy import pi
from kinbot.pp_tables import *
from kinbot.kb_trigo import *

class Fragment(StationaryPoint):
    """
    Class that creates fragments in the form of Atom objects from the ASE package.

    """
    def __init__(self, name, charge, mult, smiles='', structure=None, natom=0,\
                    atom=None, geom=None, wellorts=0, fragA=None, fragB=None):

        super(Fragment, self).__init__(self, name, charge, mult, smiles='', structure=None, natom=0,\
                                            atom=None, geom=None, wellorts=0, fragA=None, fragB=None)

    @classmethod
    def from_ase_atoms(cls, atoms, **kwargs):
        super(Fragment, cls).from_ase_atoms(cls, atoms, **kwargs)

    def set_parent_chemid(self, p_chemid):
        self.parent_chemid = p_chemid

    def set_map(self, p_map):
        self.map = np.array(p_map, dtype=int)

    def get_center_of_mass(self):
        super(Fragment, self).get_center_of_mass()

    def get_chemical_formula(self):
        all_elem = ""
        for elem in self.atom:
            all_elem += f"{elem}"
        return all_elem
    
    def set_pivot_points(self, dist, ra_indexes_in_parent):
        self.pivot_points = []
        if dist >= 12:
            self.set_pp_on_com()
        elif dist >= 10 and dist < 12:
            self.set_pp_on_com()
            self.set_ra(ra_indexes_in_parent)
            for index in self.ra:
                self.set_pp_on_ra(index)
        elif dist >= 6 and dist < 10:
            self.set_ra(ra_indexes_in_parent)
            for index in self.ra:
                self.set_pp_on_ra(index)
        elif dist >= 5 and dist < 6:
            self.set_ra(ra_indexes_in_parent)
            for index in self.ra:
                self.set_pp_on_ra(index)
            if self.natom == 1:
                pass
            else:
                self.set_pp_next_to_ra()
        elif dist < 5:
            self.set_ra(ra_indexes_in_parent)
            if self.natom == 1:
                self.set_pp_on_com()
            else:
                self.set_pp_next_to_ra()
                
    def set_ra(self, ra_indexes_in_parent): 
        #Find the atomid of all reactive atoms
        ra_indexes_in_frag = [ i for i, x in enumerate(self.map) if x == ra_indexes_in_parent]
        #Find all equivalent atoms
        if len(self.atom_uniq) != self.natom:
            for atom in ra_indexes_in_frag:
                index = self.atom_uniq.index(atom)
                for equivalent_atom in self.atom_eqv[index]:
                    if equivalent_atom not in ra_indexes_in_frag:
                        ra_indexes_in_frag.append(equivalent_atom)
        self.ra = np.array(ra_indexes_in_frag, dtype=int)

    def set_pp_on_com(self):
        coord = super(Fragment, self).get_center_of_mass()
        self.pivot_points.append(coord)

    def set_pp_on_ra(self, index):
        self.pivot_points.append(self.geom[index])

    def set_pp_next_to_ra(self):
        for index in self.ra:
            atom_type = self.get_atom_type(index)
            coord = self.get_pp_coord(index, atom_type)
            if type(coord) is list:
                for this_pp in coord:
                    self.pivot_points.append(this_pp)
            else:
                self.pivot_points.append(coord)

    def get_atom_type(self, index):
        element = self.atom(index)
        nconnect = 0
        ndouble = 0
        ntriple = 0
        for this_bond in self.bonds[index]:
            nconnect += this_bond
            match this_bond:
                case 2:
                    ndouble += 1
                case 3:
                    ntriple += 1
        atom_type = atom_type_table(element, nconnect, ndouble, ntriple)
        return atom_type
    
    def get_pp_coord(self, index, atom_type):
        match atom_type:
            case 'H' | 'C' | 'O' | 'S':
                #Create pivot point on atom
                self.set_pp_on_ra(index)
            case 'C_lin':
                #Create pivot point aligned with the bond
                ra_pos = np.array(self.geom(index), dtype=float)
                for neighbour_index, this_bond in enumerate(self.bonds[index]):
                    if this_bond != 0:
                        neighbour_pos = np.array(self.geom(neighbour_index), dtype=float)
                        break
                pp_orient = np.subtract(ra_pos, neighbour_pos)
                length = pp_lenght_table(self.atom[index])
                pp_vect = length * unit_vector(pp_orient)
                pp_coord = np.add(ra_pos, pp_vect)
                return pp_coord
            case 'C_tri':
                #Create pivot point in the middle of the big angle between the 2 fragments
                ra_pos = np.array(self.geom(index), dtype=float)
                neighbour_pos = []
                for neighbour_index, this_bond in enumerate(self.bonds[index]):
                    if this_bond != 0:
                        neighbour_pos.append(np.array(self.geom(neighbour_index), dtype=float))
                v1 = np.subtract(neighbour_pos[0], ra_pos)
                v2 = np.subtract(neighbour_pos[1], ra_pos)
                small_angle = angle_between(v1, v2)
                big_angle = 2*pi - small_angle
                axis = unit_vector(np.cross(v2, v1))
                pp_orient = np.dot(rotation_matrix(axis, big_angle/2), v1)
                length = pp_lenght_table(self.atom[index])
                pp_vect = length * unit_vector(pp_orient)
                pp_coord = np.add(ra_pos, pp_vect)
                return pp_coord
            case 'C_quad':
                #Create one or two pivot points on one or both side of the plane
                n_pp = 2
                ra_pos = np.array(self.geom(index), dtype=float)
                neighbour_pos = []
                for neighbour_index, this_bond in enumerate(self.bonds[index]):
                    if this_bond != 0:
                        neighbour_pos.append(np.array(self.geom(neighbour_index), dtype=float))
                v1 = np.subtract(neighbour_pos[0], ra_pos)
                v2 = np.subtract(neighbour_pos[1], ra_pos)
                v3 = np.subtract(neighbour_pos[2], ra_pos)
                plane = plane_from_points(v1, v2, v3)
                ra_to_plane = dist_point_to_plane(ra_pos, plane)
                if abs(ra_to_plane) >= .1:
                    #If carbon atom out of plane, only place a single pp on other side of the plane
                    n_pp = 1
                pp_list = []
                #To know in which direction to place the pivot point
                plane_direction = unit_vector(np.dot(plane[0], ra_pos))
                length = pp_lenght_table(self.atom[index])
                for i in range(n_pp):
                    pp_orient = np.array(unit_vector(plane[0])*plane_direction*np.power(-1,i), dtype=float)
                    pp_vect = length * unit_vector(pp_orient)
                    pp_coord = np.add(ra_pos, pp_vect)
                    pp_list.append(pp_coord)
                return pp_list

            case 'N_tri':
                pass
            case 'N_pyr':
                pass
            case 'O_tri':
                #Create pivot point aligned with the bond
                ra_pos = np.array(self.geom(index), dtype=float)
                for neighbour_index, this_bond in enumerate(self.bonds[index]):
                    if this_bond != 0:
                        neighbour_pos = np.array(self.geom(neighbour_index), dtype=float)
                        break
                pp_orient = np.subtract(ra_pos, neighbour_pos)
                length = pp_lenght_table(self.atom[index])
                pp_vect = length * unit_vector(pp_orient)
                pp_coord = np.add(ra_pos, pp_vect)
                return pp_coord
            case 'S_tri':
                pass
            case 'S_pyr':
                pass
            case 'S_lin':
                pass
            case 'S_bip_tri_l':
                pass
            case 'S_bip_tri':
                pass
            case 'S_quad':
                pass
            case 'S_bip_quad_t':
                pass