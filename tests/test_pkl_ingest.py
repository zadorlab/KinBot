"""check_qc must survive a {job}.pkl that vanishes while it is being ingested."""

import os
import pickle
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from ase.db import connect

from kinbot.qc import QuantumChemistry


JOB = '123_Intra_disproportionation_R_3_10'
PAYLOAD = {'sym': ['H', 'H'], 'pos': [[0., 0., 0.], [0.74, 0., 0.]],
           'name': JOB, 'data': {'energy': -1.0}}


class DeletingDB:
    """Wraps a real ase db; deletes the pkl during write, like a compute job
    that restarts and clears its stale handoff while the driver is ingesting."""

    def __init__(self, db):
        self.db = db

    def select(self, **kwargs):
        return self.db.select(**kwargs)

    def write(self, *args, **kwargs):
        os.remove(f'{JOB}.pkl')
        return self.db.write(*args, **kwargs)


class LaggingDB:
    """In-memory db whose first select after a write still shows the old rows,
    as a networked sqlite file can. Each select is a fresh generator."""

    def __init__(self):
        self.rows = []
        self.stale = 0

    def select(self, **kwargs):
        visible = self.rows[:len(self.rows) - self.stale]
        self.stale = 0
        return (row for row in visible)

    def write(self, *args, **kwargs):
        self.rows.append(kwargs)
        self.stale = 1


class TestPklIngest(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)

    def qc(self, db):
        qc = SimpleNamespace(db=db, qc='fc', queuing='local', job_ids={},
                             par={'error_missing_local': False},
                             is_in_database=lambda job: 0)
        qc.ingest_pkl = lambda job: QuantumChemistry.ingest_pkl(qc, job)
        return qc

    def write_pkl(self):
        with open(f'{JOB}.pkl', 'wb') as f:
            pickle.dump(PAYLOAD, f)

    def test_pkl_removed_mid_ingest_does_not_crash_and_row_is_written(self):
        self.write_pkl()
        db = connect('kinbot.db')
        QuantumChemistry.check_qc(self.qc(DeletingDB(db)), JOB)
        self.assertEqual(sum(1 for _ in db.select(name=JOB)), 1)
        self.assertFalse(Path(f'{JOB}.pkl').exists())

    def test_missing_pkl_does_not_reingest_existing_row(self):
        db = connect('kinbot.db')
        from ase import Atoms
        db.write(Atoms(PAYLOAD['sym'], PAYLOAD['pos']), name=JOB, data=PAYLOAD['data'])
        QuantumChemistry.check_qc(self.qc(db), JOB)
        self.assertEqual(sum(1 for _ in db.select(name=JOB)), 1)

    def test_pkl_is_ingested_and_removed_on_the_normal_path(self):
        self.write_pkl()
        db = connect('kinbot.db')
        QuantumChemistry.check_qc(self.qc(db), JOB)
        rows = list(db.select(name=JOB))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].data['energy'], -1.0)
        self.assertFalse(Path(f'{JOB}.pkl').exists())

    def test_row_visibility_wait_requeries_the_db(self):
        # with the old exhausted-generator loop this hangs forever
        self.write_pkl()
        db = LaggingDB()
        with patch('kinbot.qc.time.sleep') as sleep:
            QuantumChemistry.check_qc(self.qc(db), JOB)
        self.assertEqual(len(db.rows), 1)
        self.assertEqual(sleep.call_count, 1)
        self.assertFalse(Path(f'{JOB}.pkl').exists())

    def test_partial_pkl_is_skipped_after_timeout(self):
        Path(f'{JOB}.pkl').write_bytes(b'\x80\x04')  # truncated pickle
        qc = self.qc(connect('kinbot.db'))
        with patch('kinbot.qc.time.sleep'), self.assertLogs('KinBot', 'WARNING'):
            self.assertIsNone(QuantumChemistry.ingest_pkl(qc, JOB, max_wait=3))


if __name__ == '__main__':
    unittest.main()
