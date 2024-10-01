import numpy as np
from numpy import ndarray
from kinbot import constants
from kinbot import kb_path


class VRC_TST_Surface:
    __id__ = 0

    def __init__(self,
                 pp_coords: list[list[list[float]]],
                 dist_matrix: ndarray,
                 info: list[str]) -> None:
        """Surface object characterized by a
        distance matrix and pivot points coordinates.
        The pivot points' coordinates and the dist matrix
        are received in Angstrom
        and converted to Bohr."""
        self.dist_matrix = dist_matrix
        self.id = VRC_TST_Surface.__id__
        with open(f'{kb_path}/tpl/rotdPy_surf.tpl', 'r') as f:
            tpl = f.read()
        self.repr: str = tpl.format(surf_id=VRC_TST_Surface.__id__,
                                    info0=info[0],
                                    info1=info[1],
                                    info2=info[2],
                                    pp_coords_f1=(np.asarray(pp_coords[0]) /
                                                  constants.BOHRtoANGSTROM)
                                                  .tolist(),
                                    pp_coords_f2=(np.asarray(pp_coords[1]) /
                                                  constants.BOHRtoANGSTROM)
                                                  .tolist(),
                                    dist_matrix=(dist_matrix /
                                                 constants.BOHRtoANGSTROM)
                                                 .tolist()
                                    )
        VRC_TST_Surface.__id__ += 1

    def __repr__(self) -> str:
        return self.repr
