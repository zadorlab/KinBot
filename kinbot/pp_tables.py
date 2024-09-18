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
        'H': [0.125, 0.25, 0.375, 0.5],
        'C': [0.25, 0.5, 0.75, 1.0],
        'N': [0.25, 0.5, 0.75, 1.0],
        'O': [0.25, 0.5, 0.75, 1.0],
        'S': [0.25, 0.5, 0.75, 1.0],
        'X': [0.25, 0.5, 0.75, 1.0]}

    return pp_table
