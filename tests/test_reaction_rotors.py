"""A forming/breaking bond is not interchangeable with a spectator bond."""
import copy
import json
from pathlib import Path
import unittest

import numpy as np
from ase.build import molecule

from kinbot import frequencies, symmetry, zmatrix
from kinbot.molecular_symmetry import graph_equivalence_classes
from kinbot.stationary_pt import StationaryPoint


def departing_ts():
    data = json.loads((Path(__file__).parent / 'reference' / 'departing_atom_ts.json').read_text())
    point = StationaryPoint('departing_atom_ts', 0, 2,
                            atom=data['atom'], geom=np.array(data['geom']), wellorts=1)
    point.characterize()
    for key in ('bond', 'bonds', 'rads', 'reac_bond'):
        setattr(point, key, np.array(data[key]))
    point.find_cycle()
    point.find_conf_dihedral()
    return point


class TestReactionRotors(unittest.TestCase):
    def test_symmetric_exchange_retains_physical_external_rotation(self):
        point = StationaryPoint('H3', 0, 2, atom=['H']*3,
            geom=np.array([[-.9, 0., 0.], [0., 0., 0.], [.9, 0., 0.]]), wellorts=1)
        point.characterize()
        point.reac_bond = np.array([[0, -1, 0], [-1, 0, 1], [0, 1, 0]])
        symmetry.calculate_symmetry(point)
        self.assertEqual(point.sigma_ext, 2)
        self.assertEqual(graph_equivalence_classes(point), [0, 1, 0])

    def test_departing_atom_does_not_supply_methyl_period_but_can_define_angle(self):
        point = departing_ts()
        original_ids = list(point.atomid)
        symmetry.calculate_symmetry(point)
        self.assertEqual(point.sigma_int[0][1], 1)
        self.assertEqual(point.dihed[0], [5, 0, 1, 2])
        self.assertIn(point.dihed[0], point.conf_dihed)
        self.assertEqual(list(point.atomid), original_ids)
        # Either reference defines the angle; the physical top still
        # includes the departing atom: it moves rigidly with its TS fragment.
        old = [6, 0, 1, 2]
        self.assertEqual(frequencies.partition(point, old, point.natom),
                         frequencies.partition(point, point.dihed[0], point.natom))
        self.assertEqual(frequencies.partition(point, point.dihed[0], point.natom)[0],
                         [0, 5, 6, 7])
        _, _, _, order = zmatrix.make_zmat_from_cart(point, 0, point.geom, 0)
        self.assertEqual(order[:4], point.dihed[0])
        self.assertEqual(sorted(order), list(range(point.natom)))
        # The projection depends on the axis and top, not on which outer atom
        # defines the angle. Compare the actual projection with the same Hessian.
        rng = np.random.default_rng(7)
        h = rng.normal(size=(3 * point.natom, 3 * point.natom))
        h = h.T @ h
        old_point = copy.copy(point)
        old_point.dihed = [old] + point.dihed[1:]
        before = frequencies.get_frequencies(old_point, h, point.geom)
        after = frequencies.get_frequencies(point, h, point.geom)
        for left, right in zip(before, after):
            np.testing.assert_allclose(left, right, atol=1.e-6)

    def test_same_rule_applies_to_departing_fluorine(self):
        atoms = molecule('CH3OH')
        point = StationaryPoint('reference', 0, 1,
                                atom=atoms.get_chemical_symbols(), geom=atoms.positions)
        point.characterize()
        bond = point.bond.copy()
        elements = list(point.atom)
        for i in (2, 4, 5):
            elements[i] = 'F'
        point = StationaryPoint('fluorine_departure', 0, 1,
                                atom=elements, geom=atoms.positions, wellorts=1)
        point.characterize(bond_mx=bond)
        symmetry.calculate_symmetry(point)
        self.assertEqual(point.sigma_int[0][1], 3)
        point.reac_bond = np.zeros_like(bond)
        point.reac_bond[0, 2] = point.reac_bond[2, 0] = -1
        point.find_conf_dihedral()
        symmetry.calculate_symmetry(point)
        self.assertEqual(point.sigma_int[0][1], 1)
        self.assertEqual(point.dihed[0], [2, 0, 1, 3])
        classes = graph_equivalence_classes(point)
        self.assertNotEqual(classes[2], classes[4])
        self.assertEqual(classes[4], classes[5])
        # Reaction role does not veto an otherwise valid single-bond axis.
        point.reac_bond[0, 1] = point.reac_bond[1, 0] = 1
        point.find_conf_dihedral()
        self.assertEqual(point.dihed, [[2, 0, 1, 3]])

    def test_forming_and_breaking_roles_cannot_exchange(self):
        point = departing_ts()
        point.reac_bond[0, 6] = point.reac_bond[6, 0] = 1
        classes = graph_equivalence_classes(point)
        self.assertEqual(len({classes[5], classes[6], classes[7]}), 3)


if __name__ == '__main__':
    unittest.main()
