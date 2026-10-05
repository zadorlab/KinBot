"""Keep a specified configuration through selection, refinement and caches."""
from kinbot.species_routing import routing_name
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock
import numpy as np
from ase import Atoms
from ase.db import connect
from kinbot import constants
from kinbot.calculation import load_calculation_record
from kinbot.conformer_counting import writer_members
from kinbot.conformers import Conformers
from kinbot.optimize import Optimize
from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import optical_scope
from kinbot.stereo_routing import StereoRoutingError


class TestOpticalScopeLifecycle(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        Path('conf').mkdir()
        Path('input.json').write_text(json.dumps({'barrier_threshold': 50,
            'multi_conf_tst': 1, 'conformer_search': 1, 'high_level': 0,
            'rotor_scan': 0, 'L3_calc': 0}))
        self.par = Parameters('input.json').par
        p = StationaryPoint('test', 0, 1, smiles='C[C@H](O)CC')
        p.characterize()
        p.name = routing_name(p)
        p.energy, p.zpe = -20., .01
        p.freq = [100.] * (3 * p.natom - 6)
        self.species = p
        self.original = p.geom.copy()
        self.mirror = p.geom * [-1, 1, 1]
        self.db = connect('kinbot.db')
        self.write(f'{p.name}_well', self.original, -20.)
        latest = lambda job: list(self.db.select(name=job))[-1]
        self.qc = SimpleNamespace(qc='fc', db=self.db, publish_result=Mock(),
            get_qc_geom=lambda job, n, **kw: (0, latest(job).positions),
            get_qc_energy=lambda job, **kw: (0, latest(job).data['energy'] * constants.EVtoHARTREE),
            get_qc_zpe=lambda job, **kw: (0, latest(job).data['zpe']),
            get_qc_freq=lambda job, n, **kw: (0, latest(job).data['frequencies']))

    def write(self, job, geom, energy):
        self.db.write(Atoms(self.species.atom, positions=geom), name=job,
            data={'energy': energy / constants.EVtoHARTREE, 'zpe': .01,
                  'frequencies': self.species.freq, 'status': 'normal', 'hess': []})

    def test_lower_mirror_cannot_replace_configured_parent_or_writer_scope(self):
        p = self.species
        opt = Optimize(p, self.par, self.qc)
        reference = optical_scope(p)['identity']['id']
        p.confs = Conformers(p, self.par, self.qc)
        p.confs.conf, p.confs.conf_status = 2, [0, 0]
        self.write(p.confs.get_job_name(0), self.mirror, -20.001)
        self.write(p.confs.get_job_name(1), self.original, -20.)
        opt.scycconf, opt.sconf = 1, 0
        opt.do_optimization()
        np.testing.assert_array_equal(p.geom, self.original)
        self.assertEqual(p.source_job, p.confs.get_job_name(1))
        self.assertEqual(p.conformer_index, [1])
        self.assertEqual(writer_members(p)[0].index, 1)
        self.assertEqual(optical_scope(p)['identity']['id'], reference)
        self.assertEqual(p.conformer_inventory[0].exclusion_reason,
                         'different configured stereoisomer')

    def test_cached_low_selected_source_and_l2_reject_opposite_configuration(self):
        p = self.species
        opt = Optimize(p, self.par, self.qc)
        p.confs = Conformers(p, self.par, self.qc)
        job = f'conf/{p.name}_low'
        self.write(job, self.mirror, -21.)
        with self.assertRaisesRegex(ValueError, 'different configured stereoisomer'):
            load_calculation_record(p, self.qc, job)
        np.testing.assert_array_equal(p.geom, self.original)

        geometry, energy, _ = p.confs.lowest_conf_info()
        np.testing.assert_array_equal(geometry, self.original)
        self.assertEqual(p.confs.selected_job, p.name + '_well')
        self.assertAlmostEqual(energy, -20.)
        p.confs.conf, p.confs.conf_status = 1, [0]
        self.write(p.confs.get_job_name(0), self.original, -20.)
        geometry, energy, zpe = p.confs.lowest_conf_info()
        np.testing.assert_array_equal(geometry, self.original)
        self.assertAlmostEqual(energy, -20.)
        self.assertAlmostEqual(zpe, .01)
        self.write(opt.log_name(1), self.mirror, -21.)
        opt.compare_structures()
        self.assertEqual(opt.shigh, -999)
        np.testing.assert_array_equal(p.geom, self.original)

    def test_no_sampled_recovery_rejects_incomplete_parent_locally(self):
        p = self.species
        p.confs = Conformers(p, self.par, self.qc)
        job = f'conf/{p.name}_low'
        self.write(job, self.mirror, -21.)
        for field in ('energy', 'zpe', 'frequencies'):
            with self.subTest(missing=field):
                data = {'energy': -20./constants.EVtoHARTREE, 'zpe': .01,
                        'frequencies': p.freq, 'status': 'normal'}
                data.pop(field)
                self.db.write(Atoms(p.atom, positions=self.original), name=p.name+'_well', data=data)
                with self.assertRaisesRegex(StereoRoutingError, 'No compatible cached'):
                    p.confs.lowest_conf_info()


if __name__ == '__main__':
    unittest.main()
