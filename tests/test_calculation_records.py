"""Keep selected geometry and thermochemistry on one calculation record."""
from types import SimpleNamespace
import unittest

import numpy as np

from kinbot import constants
from kinbot.calculation import load_calculation_record, selected_calculation_job
from kinbot.conformers import Conformers
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import optical_scope
from ase.build import molecule


class TestCalculationRecords(unittest.TestCase):
    def ts_fixture(self):
        point = StationaryPoint.from_ase_atoms(molecule('H2O'))
        point.characterize()
        optical_scope(point)
        point.wellorts = 1
        return point

    def test_cached_lowest_conformer_tracks_its_calculation_source(self):
        conformers = Conformers.__new__(Conformers)
        conformers.species = self.ts_fixture()
        conformers.optical_population = 'specified'
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
        species = self.ts_fixture()
        row = SimpleNamespace(name='conformer', id=7, symbols=['O', 'H', 'H'],
                             positions=species.geom + .2, data={
            'energy': -10. / constants.EVtoHARTREE, 'zpe': .04,
            'frequencies': [-800., 400., 1500.], 'hess': np.eye(9),
            'status': 'normal'})
        requested = []
        qc = SimpleNamespace(db=SimpleNamespace(
            select=lambda name: requested.append(name) or [row]))
        species.energy, species.zpe = -9., .03
        species.freq = species.reduced_freqs = [-900.]
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
