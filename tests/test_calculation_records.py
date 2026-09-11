"""Keep selected geometry and thermochemistry on one calculation record."""
from types import SimpleNamespace
import unittest

import numpy as np

from kinbot import constants
from kinbot.calculation import load_calculation_record, selected_calculation_job
from kinbot.conformers import Conformers


class TestCalculationRecords(unittest.TestCase):
    def test_cached_lowest_conformer_tracks_its_calculation_source(self):
        conformers = Conformers.__new__(Conformers)
        conformers.species = SimpleNamespace(natom=3)
        conformers.get_name = lambda: 'ts'
        conformers.qc = SimpleNamespace(get_qc_energy=lambda job: (0, -10.),
            get_qc_zpe=lambda job: (0, .02),
            get_qc_geom=lambda job, natom: (0, np.ones((natom, 3))))
        geom, energy, zpe = conformers.lowest_conf_info()
        optimization = SimpleNamespace(species=SimpleNamespace(confs=conformers),
            par={'high_level': 0, 'conformer_search': 1}, shigh=1,
            log_name=lambda high: 'ts')
        self.assertEqual(selected_calculation_job(optimization), 'conf/ts_low')
        self.assertEqual((energy, zpe), (-10., .02))

    def test_selected_record_replaces_all_parent_properties_together(self):
        row = SimpleNamespace(positions=np.arange(9).reshape(3, 3), data={
            'energy': -10. / constants.EVtoHARTREE, 'zpe': .04,
            'frequencies': [-800., 400., 1500.], 'hess': np.eye(9),
            'status': 'normal'})
        requested = []
        qc = SimpleNamespace(db=SimpleNamespace(
            select=lambda name: requested.append(name) or [row]))
        species = SimpleNamespace(natom=3, geom=np.zeros((3, 3)), energy=-9.,
                                  zpe=.03, freq=[-900.], reduced_freqs=[-900.])
        load_calculation_record(species, qc, 'conformer')
        self.assertEqual(requested, ['conformer'])
        np.testing.assert_array_equal(species.geom, row.positions)
        self.assertAlmostEqual(species.energy, -10.)
        self.assertEqual(species.zpe, .04)
        self.assertEqual(species.freq, row.data['frequencies'])
        np.testing.assert_array_equal(species.hess, np.eye(9))
        self.assertEqual(species.source_job, 'conformer')

    def test_calculation_selection_tracks_l1_conformer_and_l2_parent(self):
        optimization = SimpleNamespace(
            species=SimpleNamespace(confs=SimpleNamespace(selected_job='conf/ts_0004')),
            par={'high_level': 0, 'conformer_search': 1}, shigh=1,
            log_name=lambda high: 'ts_high' if high else 'ts')
        self.assertEqual(selected_calculation_job(optimization), 'conf/ts_0004')
        optimization.par['high_level'] = 1
        self.assertEqual(selected_calculation_job(optimization), 'ts_high')
        optimization.par.update(high_level=0, conformer_search=0)
        self.assertEqual(selected_calculation_job(optimization), 'ts')


if __name__ == '__main__':
    unittest.main()
