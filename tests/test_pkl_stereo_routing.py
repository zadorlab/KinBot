"""Pickle ingestion preserves verified stereochemical job routing."""
import json
import os
import pickle
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from kinbot.qc import QuantumChemistry
from kinbot.species_routing import connect, prepare_qc_routing, resolve_job, routing_name
from kinbot.stationary_pt import StationaryPoint


class TestPklStereoIntegration(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        self.p = StationaryPoint('butanol', 0, 1, smiles='CC[C@H](O)C')
        self.p.characterize()
        self.mirror = StationaryPoint('mirror', 0, 1, atom=self.p.atom.copy(),
                                      geom=self.p.geom * [-1., 1., 1.])
        self.mirror.characterize()
        self.key = routing_name(self.p)
        self.other_key = routing_name(self.mirror)
        self.assertNotEqual(self.key, self.other_key)
        self.legacy = str(self.p.chemid)
        Path(self.legacy + '.json').write_text(json.dumps({
            'charge': 0, 'mult': 1,
            'structure': [v for atom, xyz in zip(self.p.atom, self.p.geom)
                          for v in [str(atom), *map(float, xyz)]]}))
        self.qc = QuantumChemistry.__new__(QuantumChemistry)
        self.qc.db = connect('kinbot.db')
        self.qc.qc, self.qc.queuing, self.qc.job_ids = 'fc', 'local', {}
        self.qc.par = {'error_missing_local': False, 'optical_population': 'specified'}
        prepare_qc_routing(self.qc, self.p)
        prepare_qc_routing(self.qc, self.mirror)

    def write_result(self, job, energy):
        Path(job + '.pkl').write_bytes(pickle.dumps({
            'name': job, 'sym': list(self.p.atom), 'pos': self.p.geom.tolist(),
            'data': {'energy': energy, 'status': 'normal', 'charge': 0, 'multiplicity': 1}}))
        Path(job + '_sella.log').write_text('done\n')

    def test_verified_legacy_pickle_and_log_are_used_under_configured_request(self):
        old, configured = self.legacy + '_well', self.key + '_well'
        self.write_result(old, -1.)
        self.assertEqual(resolve_job(self.qc.db, configured), old)
        self.assertEqual(self.qc.check_qc(configured), 'normal')
        self.assertFalse(Path(old + '.pkl').exists())
        self.assertEqual(self.qc.db.raw.get(name=old).data.energy, -1.)
        self.assertEqual(len(list(self.qc.db.raw.select(name=configured))), 0)
        self.assertEqual(self.qc.check_qc(configured), 'normal')
        self.assertEqual(len(list(self.qc.db.raw.select(name=old))), 1)

    def test_configured_pickle_supersedes_verified_legacy_pickle(self):
        old, configured = self.legacy + '_well', self.key + '_well'
        self.write_result(old, -1.)
        self.write_result(configured, -2.)
        self.assertEqual(self.qc.check_qc(configured), 'normal')
        self.assertEqual(self.qc.db.raw.get(name=configured).data.energy, -2.)
        self.assertEqual(len(list(self.qc.db.raw.select(name=old))), 0)
        self.assertTrue(Path(old + '.pkl').exists())
        self.assertFalse(Path(configured + '.pkl').exists())

    def test_other_enantiomer_does_not_consume_legacy_pickle(self):
        old, other = self.legacy + '_well', self.other_key + '_well'
        self.write_result(old, -1.)
        self.assertEqual(resolve_job(self.qc.db, other), other)
        self.assertEqual(self.qc.check_qc(other), 'error')
        self.assertTrue(Path(old + '.pkl').exists())
        self.assertEqual(len(list(self.qc.db.raw.select(name=old))), 0)
        self.assertEqual(len(list(self.qc.db.raw.select(name=other))), 0)

    def test_aliased_pickle_disappearing_during_database_write_is_tolerated(self):
        old, configured = self.legacy + '_well', self.key + '_well'
        self.write_result(old, -1.)
        raw_write = self.qc.db.raw.write

        def write_and_remove(*args, **kwargs):
            Path(old + '.pkl').unlink()
            return raw_write(*args, **kwargs)

        with patch.object(self.qc.db.raw, 'write', side_effect=write_and_remove):
            self.assertEqual(self.qc.check_qc(configured), 'normal')
        self.assertEqual(self.qc.db.raw.get(name=old).data.energy, -1.)
