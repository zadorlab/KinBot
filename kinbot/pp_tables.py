import numpy as np
from kinbot import constants


def pp_length_table() -> dict[str, list[float]]:
    """Create a default list of distances to try from a given pivot point.

    Args:
        element (str): chemical symbol

    Returns:
        list[float]: list of distances for the pivot point.
    """
    pp_table: dict[str, list[float]] = {
        'H': (np.array([0.125, 0.25, 0.375, 0.5])*constants.BOHRtoANGSTROM).tolist(),
        'C': (np.array([0.25, 0.5, 0.75, 1.0])*constants.BOHRtoANGSTROM).tolist(),
        'N': (np.array([0.25, 0.5, 0.75, 1.0])*constants.BOHRtoANGSTROM).tolist(),
        'O': (np.array([0.25, 0.5, 0.75, 1.0])*constants.BOHRtoANGSTROM).tolist(),
        'S': (np.array([0.25, 0.5, 0.75, 1.0])*constants.BOHRtoANGSTROM).tolist(),
        'X': (np.array([0.25, 0.5, 0.75, 1.0])*constants.BOHRtoANGSTROM).tolist()}

    return pp_table
