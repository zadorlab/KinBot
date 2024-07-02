import numpy as np
from kinbot import constants

class VRC_TST_Surface:
    __id__ = 0
    def __init__(self, pp_coords, dist_matrix, info):
        """Surface object characterized by a distance matrix and pivot points coordinates."""
        self.dist_matrix = (dist_matrix / constants.BOHRtoANGSTROM).tolist()
        self.centers = {"0" : np.ndarray.tolist(np.asarray(pp_coords[0])/ constants.BOHRtoANGSTROM),
                        "1" : np.ndarray.tolist(np.asarray(pp_coords[1])/ constants.BOHRtoANGSTROM)}
        self.overhead = f"""
#Surface id: {VRC_TST_Surface.__id__}
#{info[0]}
#{info[1]}
#{info[2]}"""
        VRC_TST_Surface.__id__ += 1


    def __repr__(self):
        string_rpr = f"{self.overhead}" + "\n" +\
                  f"Surface(pivotpoints={self.centers}" + ",\n" +\
                  f"          distances={self.dist_matrix}),"\
            .replace("[[", "np.array([[")\
            .replace("]]","]])")\
            .replace("]]),","]]),\n                    ")\
            .replace("])),","])),\n\n")
        return string_rpr
