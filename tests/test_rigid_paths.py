import numpy as np

from kinbot.stationary_pt import StationaryPoint


def test_rigidity_uses_the_first_shortest_path_in_atom_order():
    # The two paths from 0 to 4 have equal length. Only the route through
    # atom 2 contains a substituted double bond. Earlier motif order chose 1.
    species = StationaryPoint('two paths', 0, 1, atom=['C'] * 6,
                              geom=np.zeros((6, 3)))
    species.bond = np.zeros((6, 6), dtype=int)
    for i, j, order in [(0, 1, 1), (1, 4, 1), (0, 2, 1), (2, 4, 1),
                        (2, 3, 2), (3, 5, 1)]:
        species.bond[i, j] = species.bond[j, i] = order
    species.cycle = [0] * 6
    assert not species.rigid_along_path(4, 0)
    species.bond[1, 4] = species.bond[4, 1] = 0
    assert species.rigid_along_path(4, 0)
    species.bond[2, 4] = species.bond[4, 2] = 0
    assert not species.rigid_along_path(4, 0)
