"""A done stamp must not report an earlier step's result as the new one.

Each step of a reaction search resubmits the same job name. The job writes
its {job}.pkl and then appends 'done' to the log. Over a network file system
the driver can see the stamp before the pkl. The last database row then still
belongs to the previous step, and reporting it as 'normal' starts the next
step from the wrong geometry.
"""

import os
import pickle
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock, patch

from ase import Atoms
from ase.db import connect

from kinbot.qc import QuantumChemistry


JOB = '611491631250510000002_intra_H_migration_2_7'


def result(step):
    return {'sym': ['H', 'H'], 'pos': [[0., 0., 0.], [0., 0., 0.7 + 0.01 * step]],
            'name': JOB, 'data': {'energy': -1. - step, 'status': 'normal'}}


class TestDoneStampTiming(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        self.qc = QuantumChemistry.__new__(QuantumChemistry)
        self.qc.qc, self.qc.queuing = 'fc', 'slurm'
        self.qc.par = {'error_missing_local': True}
        self.qc.db = connect('kinbot.db')
        self.qc.job_ids = {JOB: '1'}
        self.qc._rows_at_submit = {}
        # the previous step's result is in, its pkl already consumed
        previous = result(13)
        self.qc.db.write(Atoms(previous['sym'], previous['pos']), name=JOB, data=previous['data'])
        # the queue no longer lists the job
        queue = Mock()
        queue.communicate.return_value = (b'JOBID NAME\n', b'')
        patcher = patch('kinbot.qc.subprocess.Popen', return_value=queue)
        patcher.start()
        self.addCleanup(patcher.stop)

    def stamp_done(self):
        Path(f'{JOB}_sella.log').write_text('Sella 11 ... converged\ndone\n')

    def drop_pkl(self, step):
        with open(f'{JOB}.pkl', 'wb') as handle:
            pickle.dump(result(step), handle)

    def test_stamp_before_pkl_is_not_a_result(self):
        self.qc._rows_at_submit[JOB] = 1  # step 14 submitted with one row present
        self.stamp_done()
        self.assertEqual(self.qc._check_qc(JOB), 0)
        self.assertEqual(sum(1 for _ in self.qc.db.select(name=JOB)), 1)

    def test_result_is_reported_once_its_row_is_in(self):
        self.qc._rows_at_submit[JOB] = 1
        self.stamp_done()
        self.assertEqual(self.qc._check_qc(JOB), 0)
        self.drop_pkl(14)
        self.assertEqual(self.qc._check_qc(JOB), 'normal')
        rows = list(self.qc.db.select(name=JOB))
        self.assertEqual(len(rows), 2)
        self.assertAlmostEqual(rows[-1].data['energy'], -15.)
        self.assertFalse(Path(f'{JOB}.pkl').exists())

    def test_wrapper_keeps_callers_waiting_during_the_gap(self):
        self.qc._rows_at_submit[JOB] = 1
        self.stamp_done()
        self.assertEqual(self.qc.check_qc(JOB), 'running')
        self.drop_pkl(14)
        self.assertEqual(self.qc.check_qc(JOB), 'normal')

    def test_restart_without_submission_record_trusts_the_database(self):
        self.stamp_done()
        self.assertEqual(self.qc._check_qc(JOB), 'normal')

    def test_result_written_before_submission_returns_is_accepted(self):
        # Most backends write their row from inside the job. A fast job can
        # do so before the submission command returns, so the count that a
        # result must exceed is taken before launching.
        self.qc.job_ids = {}
        self.qc.queue_job_limit, self.qc.queue_name, self.qc.slurm_feature = 0, 'test', ''
        self.qc.par['queue_template'] = ''

        def launch(command, **kwargs):
            process = Mock()
            if command[0] == 'sbatch':
                finished = result(14)
                self.qc.db.write(Atoms(finished['sym'], finished['pos']), name=JOB, data=finished['data'])
                process.communicate.return_value = (b'Submitted batch job 7\n', b'')
            else:
                process.communicate.return_value = (b'JOBID NAME\n', b'')
            return process

        with patch('kinbot.qc.subprocess.Popen', side_effect=launch):
            self.assertEqual(QuantumChemistry.submit_qc(self.qc, JOB, 1), 1)
            self.assertEqual(self.qc._rows_at_submit[JOB], 1)
            self.stamp_done()
            self.assertEqual(self.qc._check_qc(JOB), 'normal')
            self.assertEqual(self.qc.check_qc(JOB), 'normal')


if __name__ == '__main__':
    unittest.main()
