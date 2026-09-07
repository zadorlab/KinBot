"""Job-writing regressions with all scheduler submission mocked out."""

import ast
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


class TestVRCJobs(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        Path('vrctst').mkdir()
        self.qc = SimpleNamespace(
            qc='gauss', qc_command='g16', ppn=8,
            par={'calc_kwargs': {}, 'vrc_tst_scan_sella': False,
                 'vrc_tst_scan_deviation': 0.1, 'vts_ang_dev': 10,
                 'sella_kwargs': {}},
            get_qc_arguments=Mock(return_value={}), submit_qc=Mock(),
        )
        self.fragment = SimpleNamespace(
            chemid=1, charge=0, mult=2, nel=1, natom=1,
            atom=['H'], geom=[[0., 0., 0.]], bond01=[[0]],
        )

    def test_fragment_job_uses_fragment_electron_count(self):
        job = QuantumChemistry.qc_vts_frag(self.qc, self.fragment)
        self.qc.get_qc_arguments.assert_called_once_with(job, 2, 0, 1, vts=1)
        self.qc.submit_qc.assert_called_once_with(job, 1)
        ast.parse(Path(f'{job}.py').read_text())

    def test_scan_jobs_use_reacting_species_electron_count(self):
        species = SimpleNamespace(
            chemid=2, charge=0, mult=1, nel=2, natom=2,
            atom=['H', 'H'], geom=[[0., 0., 0.], [1., 0., 0.]],
        )
        reaction = SimpleNamespace(
            species=species, instance_name='reaction', scan_coo=[0, 1],
            irc_prod=SimpleNamespace(bondlist=[]), maps=[[0], [1]],
            products=[self.fragment, self.fragment],
        )
        for sella, asymptote in ((False, False), (True, False), (True, True)):
            with self.subTest(sella=sella, asymptote=asymptote):
                self.qc.par['vrc_tst_scan_sella'] = sella
                self.qc.get_qc_arguments.reset_mock()
                self.qc.submit_qc.reset_mock()
                job = QuantumChemistry.qc_vts(
                    self.qc, reaction, species.geom, 0, [], asymptote, species.geom)
                self.qc.get_qc_arguments.assert_called_once_with(job, 1, 0, 2, vts=1)
                self.qc.submit_qc.assert_called_once_with(job, 2)
                ast.parse(Path(f'{job}.py').read_text())


if __name__ == '__main__':
    unittest.main()
