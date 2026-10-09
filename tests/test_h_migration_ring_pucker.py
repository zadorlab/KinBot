"""Five-membered H-migration TS rings keep a small pucker through the bond-setting step.

With every bond frozen at its reactant value the flat ring is the constrained
minimum, and a TS search started in a plane cannot leave it (C2H5OO 1,5-H shift
was never found with Sella). The last constrained step therefore keeps the
first ring dihedral at a small angle with its current sign, carried as a
value in the fix list so the job script sets it rather than freezing whatever
the geometry has.
"""

import unittest

import numpy as np

from kinbot.reactions.reac_intra_H_migration import IntraHMigration
from kinbot.stationary_pt import StationaryPoint
from kinbot.utils import get_unique_list_of_lists


SYM = ['O', 'O', 'C', 'C', 'H', 'H', 'H', 'H', 'H']
# CH3CH2OO (O1 bonded to C3, O2 the radical oxygen)
WELL = np.array([[-0.528, 0.4717, -0.003], [-1.7192, -0.0575, 0.0067],
                 [0.5015, -0.5527, -0.0014], [1.8353, 0.1502, 0.0002],
                 [0.3432, -1.1643, 0.8857], [0.347, -1.1629, -0.89],
                 [1.9477, 0.774, -0.8855], [1.9442, 0.7768, 0.8845],
                 [2.6363, -0.5888, 0.0028]])
# end of the dihedral scan for the O2 <- H7 1,5-shift: puckered ring, +52.6 deg
STEP13 = np.array([[0.0134, 0.0458, 0.0217], [-0.4408, 0.7947, -0.9452],
                   [1.4611, 0.0439, -0.0553], [1.9269, 1.4746, -0.0882],
                   [1.766, -0.4938, 0.8387], [1.7453, -0.523, -0.9388],
                   [1.4331, 2.0075, 0.721], [3.0033, 1.5192, 0.0521],
                   [1.6685, 1.9654, -1.0248]])


def reaction(instance):
    species = StationaryPoint('ethylperoxy', 0, 2, atom=SYM, geom=WELL)
    species.characterize()
    return IntraHMigration(species, None, {}, instance, 'test')


class TestRingPucker(unittest.TestCase):
    def test_five_ring_keeps_a_signed_pucker_in_the_bond_setting_step(self):
        rxn = reaction([1, 0, 2, 3, 6])
        step, fix, change, release = rxn.get_constraints(rxn.dihstep + 1, STEP13)
        self.assertEqual(step, rxn.dihstep + 1)
        self.assertIn([2, 1, 3, 4, 15.], fix)
        self.assertNotIn([2, 1, 3, 4], release)
        self.assertIn([1, 3, 4, 7], release)
        self.assertEqual(sorted(change), [[2, 7, 1.2], [4, 7, 1.35]])

    def test_sign_follows_the_current_dihedral(self):
        rxn = reaction([1, 0, 2, 3, 6])
        mirrored = STEP13 * np.array([1., 1., -1.])
        _, fix, _, _ = rxn.get_constraints(rxn.dihstep + 1, mirrored)
        self.assertIn([2, 1, 3, 4, -15.], fix)

    def test_other_ring_sizes_release_every_dihedral(self):
        rxn = reaction([1, 0, 2, 4])  # 1,3-shift, four-membered ring
        _, fix, _, release = rxn.get_constraints(rxn.dihstep + 1, WELL)
        self.assertFalse([f for f in fix if len(f) == 5])
        self.assertEqual(release[-1], [2, 1, 3, 5])

    def test_unique_fix_entries_keep_angle_and_dihedral_order(self):
        unique = get_unique_list_of_lists([[4, 2], [2, 4], [2, 1, 3, 4], [4, 3, 1, 2],
                                           [1, 2, 3, 4], [2, 1, 3, 4, 15.], [2, 1, 3, 4, -15.]])
        self.assertEqual(unique, [[4, 2], [2, 1, 3, 4], [1, 2, 3, 4],
                                  [2, 1, 3, 4, 15.], [2, 1, 3, 4, -15.]])


if __name__ == '__main__':
    unittest.main()
