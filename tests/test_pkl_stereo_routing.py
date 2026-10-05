"""Pickle ingestion keeps literal R and S calculation names separate."""
import os
import pickle
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from ase.db import connect
from kinbot.qc import QuantumChemistry
from kinbot.run_format import ensure_current_run
from kinbot.species_routing import routing_name
from kinbot.stationary_pt import StationaryPoint


class TestPklStereoIntegration(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        ensure_current_run(create=True)
        self.p = StationaryPoint('butanol', 0, 1, smiles='CC[C@H](O)C')
        self.p.characterize()
        self.mirror = StationaryPoint('mirror', 0, 1, atom=self.p.atom.copy(),
                                      geom=self.p.geom * [-1., 1., 1.])
        self.mirror.characterize()
        self.key = routing_name(self.p)
        self.other_key = routing_name(self.mirror)
        self.assertNotEqual(self.key, self.other_key)
        self.qc = QuantumChemistry.__new__(QuantumChemistry)
        self.qc.db = connect('kinbot.db')
        self.qc.qc, self.qc.queuing, self.qc.job_ids = 'fc', 'local', {}
        self.qc.par = {'error_missing_local': False, 'optical_population': 'specified'}

    def write_result(self, point, energy):
        job = routing_name(point) + '_well'
        Path(job + '.pkl').write_bytes(pickle.dumps({
            'name': job, 'sym': list(point.atom), 'pos': point.geom.tolist(),
            'data': {'energy': energy, 'status': 'normal', 'charge': 0, 'multiplicity': 1}}))
        Path(job + '_sella.log').write_text('done\n')
        return job

    def test_configured_result_is_ingested_once(self):
        job = self.write_result(self.p, -1.)
        self.assertEqual(self.qc.check_qc(job), 'normal')
        self.assertFalse(Path(job + '.pkl').exists())
        self.assertEqual(self.qc.db.get(name=job).data.energy, -1.)
        self.assertEqual(self.qc.check_qc(job), 'normal')
        self.assertEqual(len(list(self.qc.db.select(name=job))), 1)

    def test_configured_request_does_not_consume_the_other_enantiomers_pickle(self):
        first = self.write_result(self.p, -1.)
        second = self.write_result(self.mirror, -2.)
        self.assertEqual(self.qc.check_qc(first), 'normal')
        self.assertTrue(Path(second + '.pkl').exists())
        self.assertEqual(len(list(self.qc.db.select(name=second))), 0)
        self.assertEqual(self.qc.check_qc(second), 'normal')
        self.assertEqual(self.qc.db.get(name=first).data.energy, -1.)
        self.assertEqual(self.qc.db.get(name=second).data.energy, -2.)

    def test_pickle_disappearing_during_database_write_is_tolerated(self):
        job = self.write_result(self.p, -1.)
        write = self.qc.db.write

        def write_and_remove(*args, **kwargs):
            Path(job + '.pkl').unlink()
            return write(*args, **kwargs)

        with patch.object(self.qc.db, 'write', side_effect=write_and_remove):
            self.assertEqual(self.qc.check_qc(job), 'normal')
        self.assertEqual(self.qc.db.get(name=job).data.energy, -1.)
