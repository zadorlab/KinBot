"""Keep conformer indices, sources and failed observations associated."""
from dataclasses import FrozenInstanceError
from types import SimpleNamespace
import unittest
import numpy as np
from kinbot.conformer_records import inventory, retain, update_member


class TestConformerRecords(unittest.TestCase):
    def test_inventory_is_immutable_and_keeps_discarded_sources(self):
        species = SimpleNamespace(name='test', atom=['H', 'H'])
        geometries = [np.array([[0., 0., 0.], [0., 0., .74]]), None]
        records = inventory(species, geometries, [-1., -999.], [[4000.], None], [0, 1],
                            {0: {'source_job': 'conf/test_0000',
                                 'electronic_energy_hartree': -1.01, 'zpe_hartree': .01}})
        geometries[0][1, 2] = 99.
        self.assertEqual(records[0].geometry[1][2], .74)
        with self.assertRaises(FrozenInstanceError):
            records[0].status = 'changed'
        retain(species, records, [0])
        self.assertEqual(len(species.conformer_inventory), 2)
        self.assertEqual(species.conformer_records[0].source_job, 'conf/test_0000')
        self.assertIsNone(species.conformer_inventory[1].geometry)

    def test_l2_order_and_failure_do_not_reassign_another_member(self):
        species = SimpleNamespace(name='test', atom=['H'])
        records = inventory(species, [[[0, 0, 0]], [[1, 0, 0]]],
                            [-1., -2.], [[], []], [0, 0])
        retain(species, records, [0, 1])
        species.conformer_index = [1, 0]
        species.conformer_geom = [[[2, 0, 0]], [[3, 0, 0]]]
        species.conformer_energy = [-3., -4.]
        species.conformer_zeroenergy = [-2.9, -3.9]
        species.conformer_freq = [[], []]
        update_member(species, 1, source_job='test_0001_high')
        update_member(species, 0, accepted=False)
        self.assertEqual(species.conformer_records[1].geometry, ((2., 0., 0.),))
        self.assertEqual(species.conformer_records[1].electronic_energy_hartree, -3.)
        self.assertFalse(species.conformer_records[0].retained)
        self.assertEqual(species.conformer_inventory[0].member_id, 'test:conformer:0')

    def test_misaligned_arrays_fail_before_association(self):
        with self.assertRaisesRegex(ValueError, 'different lengths'):
            inventory(SimpleNamespace(atom=['H']), [], [-1.], [], [])

    def test_all_failed_search_keeps_failure_indices_without_properties(self):
        records = inventory(SimpleNamespace(atom=['H']), [], [], [], [1, 1])
        self.assertEqual([record.index for record in records], [0, 1])
        self.assertTrue(all(record.status == 'failed' for record in records))


if __name__ == '__main__':
    unittest.main()
