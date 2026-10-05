"""Keep the ring lists used by reaction and symmetry calculations."""
import numpy as np
import pytest

from kinbot.stationary_pt import StationaryPoint


@pytest.mark.parametrize('edges, expected', [
    ([(0, 1), (1, 2)], []),
    ([(0, 1), (1, 2), (2, 0)], [[0, 1, 2]]),
    ([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
      (4, 6), (6, 7), (7, 8), (8, 9), (9, 5)],
     [[0, 1, 2, 3, 4, 5], [4, 5, 9, 8, 7, 6],
      [0, 1, 2, 3, 4, 6, 7, 8, 9, 5]]),
    ([(0, 1), (1, 4), (0, 2), (2, 4), (0, 3), (3, 4)],
     [[0, 1, 4, 2], [0, 1, 4, 3], [0, 2, 4, 3]]),
    ([(i, j) for i in range(4) for j in range(i)],
     [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3], [0, 1, 2, 3]]),
])
def test_closed_cycle_search_keeps_legacy_ring_order(edges, expected):
    n = max(max(edge) for edge in edges) + 1
    point = StationaryPoint('rings', 0, 1, atom=['C'] * n, geom=np.ones((n, 3)))
    point.bond = np.zeros((n, n), dtype=int)
    for a, b in edges:
        point.bond[a, b] = point.bond[b, a] = 1
    point.find_cycle()
    assert point.cycle_chain == expected
    assert point.cycle == [int(any(i in ring for ring in expected)) for i in range(n)]
