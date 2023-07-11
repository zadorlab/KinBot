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

    def set_parent_chemid(self, p_chemid):
        self.parent_chemid = p_chemid

    def set_map(self, p_map)
        self.map = np.array(p_map, dtype=int)
