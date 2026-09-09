"""PES consumes accepted HIR-refinement energies without extra IRCs."""
import os
import logging
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from ase.build import molecule
from ase.db import connect

from kinbot import constants
from kinbot.calculation import (record_optimization_selection,
                               optimization_selection_row)
from kinbot.pes import get_energy
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
        point = StationaryPoint.from_ase_atoms(self.atoms)
        point.characterize()
        self.well = str(point.chemid)
        Path(self.well).mkdir()
        self.db = connect(f'{self.well}/kinbot.db')

    def assertEnergy(self, actual, expected):
        for value, target in zip(actual, expected):
            self.assertAlmostEqual(value, target, places=12)

    def record(self, job, energy, zpe=.1, status='normal', atoms=None):
        return self.db.write(self.atoms if atoms is None else atoms, name=job,
                             data={'energy': energy / constants.EVtoHARTREE,
                                   'zpe': zpe, 'status': status,
                                   'frequencies': [100., 200., 300.]})

    def selection(self, base, high=0, conf=0, source=None):
        source = source or base + '_hir_restart_1'
        self.record(source, -100.02, .1002)
        optimization = SimpleNamespace(
            species=SimpleNamespace(source_job=source, geom=self.atoms.positions.copy(),
                                    energy=-100.02, zpe=.1002, freq=[100., 200., 300.]),
            qc=SimpleNamespace(db=self.db),
            shigh=1, shir=1, defer_hir=False,
            par={'high_level': high, 'conformer_search': conf, 'rotor_scan': 1},
            log_name=lambda high: base)
        record_optimization_selection(optimization)
        return optimization

    def test_selected_l1_refinement_reaches_pes_for_wells_and_ts(self):
        for ts in (0, 1):
            for conf in (0, 1):
                with self.subTest(ts=ts, conf=conf):
                    job = self.well + '_test_ts' if ts else self.well
                    base = job if ts else job + '_well'
                    conventional = f'conf/{job}_low' if conf else base
                    self.record(conventional, -100.)
                    optimization = self.selection(base, conf=conf)
                    self.assertEnergy(get_energy([self.well], job, ts, 0,
                                                conf=conf, rotor_scan=1), (-100.02, .1002))
                    # Ordinary terminal polls do not duplicate publication.
                    record_optimization_selection(optimization)
                    self.assertEqual(self.db.count(name=f'optimization/{base}'), conf + 1)

    def test_reference_calculations_keep_their_own_energy(self):
        job = self.well + '_test_ts'
        self.record(job, -100.)
        self.record(job + '_mp2', -101.)
        self.record(job + '_bls', -102.)
        self.record(job + '_high', -103.)
        self.selection(job)
        self.assertEnergy(get_energy([], job, 1, 0, mp2=1), (-101., .1))
        self.assertEnergy(get_energy([], job, 1, 0, bls=1), (-102., .1))
        self.assertEnergy(get_energy([], job, 1, 1), (-103., .1))
        self.assertEnergy(get_energy([], job, 1, 0, rotor_scan=0), (-100., .1))

    def test_missing_or_superseded_source_does_not_resurrect_an_old_selection(self):
        job = self.well + '_test_ts'
        self.record(job, -100.)
        optimization = self.selection(job)
        source = optimization.species.source_job
        self.record(source, -100.03)
        record_optimization_selection(optimization)
        self.assertEqual(self.db.count(name=f'optimization/{job}'), 1)
        self.assertEnergy(get_energy([], job, 1, 0), (-100., .1))
        self.record(source, -100.03, status='error')
        self.assertEnergy(get_energy([], job, 1, 0), (-100., .1))
        self.db.delete([row.id for row in self.db.select(name=source)])
        self.assertEnergy(get_energy([], job, 1, 0), (-100., .1))

    def test_new_conventional_calculation_supersedes_old_pointer(self):
        job = self.well + '_test_ts'
        self.record(job, -100.)
        self.selection(job)
        self.record(job, -100.04, .1004)
        self.assertEnergy(get_energy([], job, 1, 0), (-100.04, .1004))

    def test_new_acceptance_of_a_cached_source_publishes_after_legacy_updates(self):
        job = self.well + '_test_ts'
        self.record(job, -100.)
        optimization = self.selection(job)
        self.record(job, -100.04, .1004)
        resumed = SimpleNamespace(**{key: value for key, value in vars(optimization).items()
                                     if key != '_recorded_selection'})
        record_optimization_selection(resumed)
        self.assertEnergy(get_energy([], job, 1, 0), (-100.02, .1002))

    def test_selected_well_still_requires_the_requested_identity(self):
        base = self.well + '_well'
        self.record(base, -100.)
        other = base + '_wrong_product'
        self.record(other, -110., atoms=molecule('H2O'))
        optimization = self.selection(base)
        optimization.species.source_job = other
        optimization.species.energy, optimization.species.zpe = -110., .1
        optimization.species.geom = molecule('H2O').positions
        record_optimization_selection(optimization)
        with self.assertRaises(ValueError):
            get_energy([self.well], self.well, 0, 0)

    def test_pending_and_deferred_optimizations_are_not_published(self):
        base = self.well + '_well'
        optimization = self.selection(base)
        for state in ({'shir': 0}, {'shir': 1, 'shigh': -999},
                      {'shigh': 1, 'defer_hir': True}):
            optimization.__dict__.update(state)
            optimization.species.source_job = 'not_accepted'
            record_optimization_selection(optimization)
        self.assertEqual(self.db.count(name=f'optimization/{base}'), 1)
        self.assertIsNone(optimization_selection_row(self.db, base, 0, 1))


if __name__ == '__main__':
    unittest.main()
