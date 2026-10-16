"""Gaussian jobs run through Sella stamp 'done' in {job}_sella.log.

check_qc read {job}.log for every Gaussian job, which the Sella templates
leave at Gaussian's own 'Normal termination' line, so a finished Sella
optimization was never reported as done. With use_sella the native
templates still run some jobs, and those keep their stamp in {job}.log.
"""

import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock, patch

from ase import Atoms
from ase.db import connect

from kinbot.qc import QuantumChemistry


JOB = '611491631250510000002_well'


class TestSellaGaussianDoneStamp(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        self.qc = QuantumChemistry.__new__(QuantumChemistry)
        self.qc.qc, self.qc.queuing, self.qc.use_sella = 'gauss', 'slurm', True
        self.qc.par = {'error_missing_local': True}
        self.qc.db = connect('kinbot.db')
        self.qc.job_ids = {JOB: '1'}
        self.qc._rows_at_submit = {JOB: 0}
        self.qc.db.write(Atoms('HH', [[0., 0., 0.], [0., 0., 0.74]]), name=JOB,
                         data={'energy': -1., 'status': 'normal'})
        queue = Mock()
        queue.communicate.return_value = (b'JOBID NAME\n', b'')
        patcher = patch('kinbot.qc.subprocess.Popen', return_value=queue)
        patcher.start()
        self.addCleanup(patcher.stop)
        Path(f'{JOB}.log').write_text(' Normal termination of Gaussian 16.\n')

    def test_sella_job_is_done_when_its_sella_log_is_stamped(self):
        self.assertEqual(self.qc._check_qc(JOB), 0)
        Path(f'{JOB}_sella.log').write_text('Sella 7 ... converged\ndone\n')
        self.assertEqual(self.qc._check_qc(JOB), 'normal')

    def test_native_gaussian_job_keeps_its_own_stamp(self):
        Path(f'{JOB}.log').write_text(' Normal termination of Gaussian 16.\ndone\n')
        self.assertEqual(self.qc._check_qc(JOB), 'normal')

    def test_without_sella_the_sella_log_is_ignored(self):
        self.qc.use_sella = False
        Path(f'{JOB}_sella.log').write_text('stale\ndone\n')
        self.assertEqual(self.qc._check_qc(JOB), 0)


if __name__ == '__main__':
    unittest.main()
