"""Job-writing regressions with all scheduler submission mocked out."""

import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from kinbot.qc import QuantumChemistry


class TestSchedulerJobs(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)

    def test_pbs_uses_requested_processors_and_records_job_id(self):
        for requested, expected in ((4, 4), (0, 1)):
            with self.subTest(requested=requested):
                qc = SimpleNamespace(
                    check_qc=lambda job: 0, queue_job_limit=0,
                    par={'queue_template': ''}, queuing='pbs', qc='gauss',
                    queue_name='test', job_ids={},
                )
                process = Mock()
                process.communicate.return_value = (b'123.server\n', b'')
                with patch('kinbot.qc.subprocess.Popen', return_value=process) as submit:
                    result = QuantumChemistry.submit_qc(qc, 'test_job', requested)
                script = Path('test_job.pbs').read_text()
                self.assertIn(f'nodes=1:ppn={expected}', script)
                self.assertIn(f'OMP_NUM_THREADS={expected}', script)
                self.assertEqual(qc.job_ids, {'test_job': '123'})
                self.assertEqual(result, 1)
                self.assertEqual(submit.call_args.args[0], ['qsub', 'test_job.pbs'])


if __name__ == '__main__':
    unittest.main()
