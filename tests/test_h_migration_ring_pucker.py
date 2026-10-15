"""Five-membered H-migration TS rings keep a small pucker through the bond-setting step.

With every bond frozen at its reactant value the flat ring is the constrained
minimum, and a TS search started in a plane cannot leave it (C2H5OO 1,5-H shift
was never found with Sella). The last constrained step therefore sets the
first ring dihedral to a small angle with its current sign as a geometry
change: the starting geometry is moved there and the four atoms are then
frozen, which every backend holds (a Gaussian ModRedundant line carrying a
value and F does not).
"""

import ast
import os
import unittest
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np

from kinbot import geometry, modify_geom
from kinbot.reac_family import carry_out_reaction
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
        self.assertIn([2, 1, 3, 4, 15.], change)
        self.assertFalse([f for f in fix if len(f) > 4])
        self.assertNotIn([2, 1, 3, 4], release)
        self.assertIn([1, 3, 4, 7], release)
        self.assertEqual(change, [[2, 7, 1.2], [4, 7, 1.35], [2, 1, 3, 4, 15.]])  # pucker last

    def test_sign_follows_the_current_dihedral(self):
        rxn = reaction([1, 0, 2, 3, 6])
        mirrored = STEP13 * np.array([1., 1., -1.])
        _, _, change, _ = rxn.get_constraints(rxn.dihstep + 1, mirrored)
        self.assertIn([2, 1, 3, 4, -15.], change)

    def test_other_ring_sizes_release_every_dihedral(self):
        rxn = reaction([1, 0, 2, 4])  # 1,3-shift, four-membered ring
        _, _, change, release = rxn.get_constraints(rxn.dihstep + 1, WELL)
        self.assertFalse([c for c in change if len(c) == 5])
        self.assertEqual(release[-1], [2, 1, 3, 5])

    def test_the_changes_move_the_geometry_to_the_pucker_and_the_bond_targets(self):
        rxn = reaction([1, 0, 2, 3, 6])
        _, _, change, _ = rxn.get_constraints(rxn.dihstep + 1, STEP13)
        zero_based = [[c - 1 for c in ci[:-1]] + [ci[-1]] for ci in change]
        success, moved = modify_geom.modify_coordinates(
            rxn.species, 'test', STEP13.copy(), zero_based, rxn.species.bond)
        self.assertEqual(success, 1)
        self.assertAlmostEqual(geometry.calc_dihedral(*moved[[1, 0, 2, 3]])[0], 15., delta=0.2)
        self.assertAlmostEqual(np.linalg.norm(moved[1] - moved[6]), 1.2, delta=0.03)
        self.assertAlmostEqual(np.linalg.norm(moved[3] - moved[6]), 1.35, delta=0.03)
        # The same changes with the pucker first lose it to the bond refinement.
        _, flat = modify_geom.modify_coordinates(
            rxn.species, 'test', STEP13.copy(), zero_based[-1:] + zero_based[:-1], rxn.species.bond)
        self.assertLess(abs(geometry.calc_dihedral(*flat[[1, 0, 2, 3]])[0]), 5.)

    def test_gaussian_job_starts_at_the_pucker_and_freezes_it_without_a_value(self):
        rxn = reaction([1, 0, 2, 3, 6])
        rxn.qc = SimpleNamespace(
            qc='gauss', use_sella=False, ppn=1,
            check_qc=lambda name: 'normal',
            is_in_database=lambda name: False,
            get_qc_arguments=lambda *args, **kwargs: {'opt': 'ModRedun,CalcFC'},
            get_qc_geom=lambda name, natom, **kwargs: (0, STEP13.copy()),
            submit_qc=lambda *args, **kwargs: 1)
        rxn.par = {'sella_kwargs': {}, 'sella_fmax': 0.01, 'sella_steps': 100, 'calc_kwargs': {},
                   'fc_model_path': '', 'fc_task_name': '', 'fc_device': ''}
        with TemporaryDirectory() as directory:
            previous = os.getcwd()
            os.chdir(directory)
            try:
                carry_out_reaction(rxn, rxn.dihstep + 1, 'g16')
                script = open(rxn.instance_name + '.py').read()
            finally:
                os.chdir(previous)
        kwargs = geom = None
        for node in ast.parse(script).body:
            if isinstance(node, ast.Assign) and any(getattr(t, 'id', '') == 'kwargs' for t in node.targets):
                kwargs = ast.literal_eval(node.value)
            if (isinstance(node, ast.Assign) and isinstance(node.value, ast.Call)
                    and getattr(node.value.func, 'id', None) == 'Atoms'):
                geom = np.array([ast.literal_eval(k.value) for k in node.value.keywords
                                 if k.arg == 'positions'][0])
        lines = kwargs['addsec'].splitlines()
        self.assertIn('2 1 3 4 F', lines)
        self.assertFalse([line for line in lines if '15.0' in line])
        self.assertIn('1 3 4 7 A', lines)
        self.assertAlmostEqual(geometry.calc_dihedral(*geom[[1, 0, 2, 3]])[0], 15., delta=0.2)
        self.assertAlmostEqual(np.linalg.norm(geom[1] - geom[6]), 1.2, delta=0.03)

    def test_unique_fix_entries_keep_angle_and_dihedral_order(self):
        unique = get_unique_list_of_lists([[4, 2], [2, 4], [2, 1, 3, 4], [4, 3, 1, 2],
                                           [1, 2, 3, 4], [1, 2, 3], [3, 2, 1], [1, 3, 2]])
        self.assertEqual(unique, [[4, 2], [2, 1, 3, 4], [1, 2, 3, 4], [1, 2, 3], [1, 3, 2]])


if __name__ == '__main__':
    unittest.main()
