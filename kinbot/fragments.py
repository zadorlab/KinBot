from ase.db import connect
from ase import Atoms
from kinbot.stationary_pt import StationaryPoint

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
        super(Fragment, self).from_ase_atoms(cls, atoms, **kwargs)

    def set_parent_chemid(self, p_chemid):
        self.parent_chemid = p_chemid

    def set_map(self, p_map):
        self.map = np.array(p_map, dtype=int)

    def get_center_of_mass(self):
        super(Fragment, self).get_center_of_mass()

    def set_pivot_points(self, distances, ra_indexes_in_parent):
        for d in distances: #Distances must be in Angstrom
            self.pivot_points = []
            if d >= 12:
                self.set_pp_on_com()
            else if d >= 10 and d < 12:
                self.set_pp_on_com()
                self.set_ra(ra_indexes_in_parent)
                self.set_pp_on_ra()
            else if d >= 6 and d < 10:
                self.set_ra(ra_indexes_in_parent)
                self.set_pp_on_ra()
            else if d >= 5 and d < 6:
                self.set_ra(ra_indexes_in_parent)
                self.set_pp_on_ra()
                if self.natom == 1:
                    pass
                else:
                    self.set_pp_on_rc()
            else if d < 5:
                self.set_ra(ra_indexes_in_parent)
                if self.natom == 1:
                    self.set_pp_on_ra()
                else:
                    self.set_pp_on_rc()
                
    def set_ra(self, ra_indexes_in_parent): 
        #Find the atomid of all reactive atoms
        ra_indexes_in_frag = [ i for i, x in enumerate(self.map) if x == ra_indexes_in_parent]
        #Find all equivalent atoms
        if len(self.atom_uniq) != self.natom:
            for atom in initial_ra_indexes_in_frag:
                index = self.atom_uniq.index(atom)
                for equivalent_atom in self.atom_eqv[index]:
                    if equivalent_atom not in ra_indexes_in_frag:
                        ra_indexes_in_frag.append(equivalent_atom)
        self.ra = np.array(ra_indexes_in_frag, dtype=int)

    def set_pp_on_com(self):
        coord = super(Fragment, self).get_center_of_mass()
        self.pivot_point.append(coord)

    def set_pp_on_ra(self):
        for index in self.ra:
            self.pivot_points.append(self.geom[index])

    def set_pp_on_rc(self):
        
