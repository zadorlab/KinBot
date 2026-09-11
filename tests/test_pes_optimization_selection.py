"""Accepted calculations reach PES and restarts through conventional names."""
import logging
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from ase.build import molecule
from ase.db import connect

from kinbot import constants
from kinbot.calculation import publish_optimization_result
from kinbot.conformers import Conformers
from kinbot.parameters import Parameters
from kinbot.pes import get_energy
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint


class TestPESOptimizationSelection(unittest.TestCase):
    def setUp(self):
        logger = patch('kinbot.pes.logger', logging.getLogger('KinBot'), create=True)
        logger.start()
        self.addCleanup(logger.stop)
        temp = TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temp.name)
        self.atoms = molecule('CH3OH')
        self.point = StationaryPoint.from_ase_atoms(self.atoms)
        self.point.characterize()
        self.well = str(self.point.chemid)
        Path(self.well).mkdir()
        Path('conf').mkdir()
        self.db = connect(f'{self.well}/kinbot.db')
        Path('input.json').write_text('{"barrier_threshold": 100, "rotor_scan": 0}')
        self.par = Parameters('input.json').par
        self.par.update(high_level=0, conformer_search=0, rotor_scan=1, queuing='local')
        self.qc = QuantumChemistry(self.par)
        self.qc.db = self.db
        self.freq = [100.] * 12

    def assertEnergy(self, actual, expected):
        for value, target in zip(actual, expected):
            self.assertAlmostEqual(value, target, places=12)

    def record(self, job, energy, zpe=.1, status='normal', atoms=None, **extra):
        return self.db.write(self.atoms if atoms is None else atoms, name=job,
            data=dict(energy=energy / constants.EVtoHARTREE, zpe=zpe, status=status,
                      frequencies=self.freq, **extra))

    def selection(self, job, ts=0, high=0, conf=0, just_high=False, publish=True):
        base = job if ts else job + '_well'
        source = base + '_hir_restart_1'
        self.record(source, -100.02, .1002, hess=np.eye(18))
        Path(source + '.log').write_text('accepted calculation B\ndone\n')
        optimization = SimpleNamespace(
            name=job, just_high=just_high,
            species=SimpleNamespace(source_job=source, geom=self.atoms.positions.copy(),
                                    energy=-100.02, zpe=.1002, freq=self.freq),
            qc=self.qc, shigh=1, shir=1, defer_hir=False,
            par={'high_level': high, 'conformer_search': conf, 'rotor_scan': 1},
            log_name=lambda level: base + ('_high' if level else ''))
        if publish:
            publish_optimization_result(optimization)
        return optimization

    def test_selected_result_reaches_conventional_readers_for_wells_and_ts(self):
        for ts in (0, 1):
            for conf in (0, 1):
                for high in (0, 1):
                    with self.subTest(ts=ts, conf=conf, high=high):
                        job = self.well + '_test_ts' if ts else self.well
                        base = job if ts else job + '_well'
                        target = base + '_high' if high else (f'conf/{job}_low' if conf else base)
                        self.record(target, -100.)
                        Path(target + '.log').write_text('old calculation A\ndone\n')
                        opt = self.selection(job, ts, high, conf)
                        self.assertEnergy(get_energy([self.well], job, ts, high, conf=conf), (-100.02, .1002))
                        row = list(self.db.select(name=target))[-1]
                        np.testing.assert_array_equal(row.positions, self.atoms.positions)
                        np.testing.assert_array_equal(row.data['hess'], np.eye(18))
                        self.assertEqual(list(row.data['frequencies']), self.freq)
                        self.assertEqual(Path(target + '.log').read_text(), 'accepted calculation B\ndone\n')
                        self.assertTrue(any(p.read_text() == 'old calculation A\ndone\n'
                                            for p in Path(target).parent.glob(Path(target).name + '.log.restart_*')))
                        count = self.db.count(name=target)
                        publish_optimization_result(opt)
                        self.assertEqual(self.db.count(name=target), count)
        self.assertEqual(list(self.db.select('name=optimization/*')), [])

    def test_lowest_conformer_restart_reads_the_published_result(self):
        self.record(f'conf/{self.well}_low', -100.)
        self.selection(self.well, conf=1)
        confs = Conformers(self.point, self.par, self.qc)
        geom, energy, zpe = confs.lowest_conf_info()
        self.assertEnergy((energy, zpe), (-100.02, .1002))
        np.testing.assert_array_equal(geom, self.atoms.positions)

    def test_reference_and_other_level_calculations_keep_their_energy(self):
        job = self.well + '_test_ts'
        for suffix, energy in (('', -100.), ('_mp2', -101.), ('_bls', -102.), ('_high', -103.)):
            self.record(job + suffix, energy)
        self.selection(job, ts=1)
        self.assertEnergy(get_energy([], job, 1, 0, mp2=1), (-101., .1))
        self.assertEnergy(get_energy([], job, 1, 0, bls=1), (-102., .1))
        self.assertEnergy(get_energy([], job, 1, 1), (-103., .1))
        self.assertEnergy(get_energy([], job, 1, 0), (-100.02, .1002))

    def test_unaccepted_source_changes_leave_last_accepted_result_intact_and_warn(self):
        job = self.well + '_test_ts'
        opt = self.selection(job, ts=1)
        for status in ('normal', 'error'):
            self.record(opt.species.source_job, -110., status=status)
            with self.assertLogs('KinBot', level='WARNING'):
                publish_optimization_result(opt)
            self.assertEnergy(get_energy([], job, 1, 0), (-100.02, .1002))
        self.db.delete([row.id for row in self.db.select(name=opt.species.source_job)])
        with self.assertLogs('KinBot', level='WARNING'):
            publish_optimization_result(opt)
        self.assertEnergy(get_energy([], job, 1, 0), (-100.02, .1002))

    def test_later_conventional_result_and_new_acceptance_follow_last_row_wins(self):
        job = self.well + '_test_ts'
        opt = self.selection(job, ts=1)
        self.record(job, -100.04, .1004)
        self.assertEnergy(get_energy([], job, 1, 0), (-100.04, .1004))
        publish_optimization_result(opt)
        self.assertEnergy(get_energy([], job, 1, 0), (-100.02, .1002))

    def test_pending_deferred_and_failed_optimizations_do_not_publish(self):
        opt = self.selection(self.well, publish=False)
        base = opt.log_name(0)
        for values in ({'shir': 0}, {'shir': 1, 'shigh': -999}, {'shigh': 1, 'defer_hir': True}):
            opt.__dict__.update(values)
            publish_optimization_result(opt)
        self.assertEqual(self.db.count(name=base), 0)

    def test_just_high_does_not_publish_under_conformer_name(self):
        opt = self.selection(self.well, conf=1, just_high=True)
        self.assertEqual(self.db.count(name=f'conf/{self.well}_low'), 0)
        self.assertEqual(list(self.db.select(name=opt.log_name(0)))[-1].data['status'], 'normal')

    def test_native_artifacts_are_replaced_together_and_originals_preserved(self):
        for backend, suffix in (('gauss', '.fchk'), ('qchem', '_freq.out'), ('fc', '_sella.log')):
            with self.subTest(backend=backend):
                job = self.well + '_' + backend + '_ts'
                opt = self.selection(job, ts=1, publish=False)
                source = opt.species.source_job
                Path(source + suffix).write_text('new native data')
                Path(job + suffix).write_text('old native data')
                Path(job + '.chk').write_text('old checkpoint with no replacement')
                publish_optimization_result(opt)
                self.assertEqual(Path(job + suffix).read_text(), 'new native data')
                self.assertEqual(Path(source + suffix).read_text(), 'new native data')
                self.assertFalse(Path(job + '.chk').exists())
                self.assertTrue(list(Path('.').glob(job + '.chk.restart_*')))

    def test_interrupted_publication_does_not_expose_new_properties_with_old_files(self):
        opt = self.selection(self.well, publish=False)
        base = opt.log_name(0)
        self.record(base, -100.)
        Path(base + '.log').write_text('old result\ndone\n')
        with patch('kinbot.qc.copyfile', side_effect=OSError('copy failed')):
            with self.assertRaises(OSError):
                publish_optimization_result(opt)
        self.assertEqual(list(self.db.select(name=base))[-1].data['status'], 0)
        with self.assertRaises(ValueError):
            get_energy([self.well], self.well, 0, 0)
        publish_optimization_result(opt)
        self.assertEnergy(get_energy([self.well], self.well, 0, 0), (-100.02, .1002))

    def test_running_target_is_not_overwritten(self):
        opt = self.selection(self.well, publish=False)
        with patch.object(self.qc, 'check_qc', return_value='running'):
            with self.assertRaises(ValueError):
                publish_optimization_result(opt)
        self.assertEqual(self.db.count(name=opt.log_name(0)), 0)


if __name__ == '__main__':
    unittest.main()
