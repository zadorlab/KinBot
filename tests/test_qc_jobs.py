"""Job-writing regressions with all scheduler submission mocked out."""

import ast
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import threading
import time
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch
import numpy as np
from ase import Atoms

from kinbot.qc import (QuantumChemistry, _export_slurm_excluded_nodes,
                       _slurm_job_active, _worker_exit_status)
from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint
from kinbot.utils import clean_files, too_far


class TestSchedulerJobs(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)

    def test_ts_refinement_uses_the_requested_level_in_the_actual_job(self):
        Path('input.json').write_text(json.dumps({'barrier_threshold': 100.,
            'method': 'HF', 'basis': 'sto-3g',
            'high_level_method': 'MP2', 'high_level_basis': '6-31g'}))
        qc = QuantumChemistry(Parameters('input.json', show_warnings=False).par)
        qc.submit_qc = Mock()
        point = SimpleNamespace(name='ts', wellorts=1, mult=1, charge=0, nel=2,
                                atom=['H', 'H'], geom=[[0., 0., 0.], [1., 0., 0.]])
        for level, suffix, method, basis in ((0, '_hir_restart_1', 'hf', 'sto-3g'),
                                             (1, '_high', 'mp2', '6-31g')):
            qc.qc_opt_ts(point, point.geom, high_level=level, ext=suffix)
            script = Path(f'ts{suffix}.py').read_text().lower()
            ast.parse(script)
            self.assertIn(f"'method': '{method}'", script)
            self.assertIn(f"'basis': '{basis}'", script)

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
                self.assertIn(
                    f'{shlex.quote(sys.executable)} test_job.py', script)
                self.assertEqual(qc.job_ids, {'test_job': '123'})
                self.assertEqual(result, 1)
                self.assertEqual(submit.call_args.args[0], ['qsub', 'test_job.pbs'])

    def test_default_slurm_uses_requested_partition_and_submitting_python(self):
        qc = SimpleNamespace(
            check_qc=lambda job: 0, queue_job_limit=0,
            par={'queue_template': ''}, queuing='slurm', qc='fc',
            queue_name='day-long-cpu', slurm_feature='', job_ids={},
        )
        process = Mock()
        process.communicate.return_value = (b'Submitted batch job 123\n', b'')
        discovery = Mock(returncode=0, stdout='04:00:00\n', stderr='')
        with patch('kinbot.qc.subprocess.Popen', return_value=process), \
                patch('kinbot.qc.subprocess.run', return_value=discovery), \
                patch.dict('kinbot.qc._SLURM_WALLTIME_CACHE', {}, clear=True):
            QuantumChemistry.submit_qc(qc, 'test_job', 8)
        script = Path('test_job.sbatch').read_text()
        self.assertIn('#SBATCH --partition=day-long-cpu', script)
        self.assertNotIn('#SBATCH -q ', script)
        self.assertIn('#SBATCH --time=04:00:00', script)
        self.assertIn(f'{shlex.quote(sys.executable)} test_job.py', script)

    def test_ordinary_slurm_worker_sources_automatic_backend_setup(self):
        qc = SimpleNamespace(
            check_qc=lambda job: 0, queue_job_limit=0,
            par={'queue_template': ''}, queuing='slurm', qc='gauss',
            qc_command='g16', queue_name='day-long-cpu', slurm_feature='',
            job_ids={},
        )
        qc._worker_site_setup = QuantumChemistry._worker_site_setup.__get__(qc)
        process = Mock()
        process.communicate.return_value = (b'Submitted batch job 123\n', b'')
        discovery = Mock(returncode=0, stdout='04:00:00\n', stderr='')
        with patch('kinbot.qc.subprocess.Popen', return_value=process), \
                patch('kinbot.qc.subprocess.run', return_value=discovery), \
                patch('kinbot.anl.site.render_site_setup',
                      return_value='# resolved Gaussian setup\n'), \
                patch.dict('kinbot.qc._SLURM_WALLTIME_CACHE', {}, clear=True):
            QuantumChemistry.submit_qc(qc, 'automatic_site_job', 8)
        script = Path('automatic_site_job.sbatch').read_text()
        setup = Path('.kinbot_site_setup_gaussian.sh').resolve()
        assert setup.read_text() == '# resolved Gaussian setup\n'
        assert 'export KINBOT_BACKEND=gaussian' in script
        assert f'source {setup}' in script
        assert script.index(f'source {setup}') < script.index(
            f'{shlex.quote(sys.executable)} automatic_site_job.py')

    def test_site_can_select_bundled_slurm_partition_template(self):
        import kinbot

        template = str(Path(kinbot.__file__).parent / 'tpl' /
                       'slurm_partition.tpl')
        qc = SimpleNamespace(
            check_qc=lambda job: 0, queue_job_limit=0,
            par={'queue_template': template}, queuing='slurm', qc='fc',
            queue_name='day-long-cpu', slurm_feature='', job_ids={},
        )
        process = Mock()
        process.communicate.return_value = (b'Submitted batch job 123\n', b'')
        discovery = Mock(returncode=0, stdout='1-00:00:00\n', stderr='')
        with patch('kinbot.qc.subprocess.Popen', return_value=process), \
                patch('kinbot.qc.subprocess.run', return_value=discovery), \
                patch.dict('kinbot.qc._SLURM_WALLTIME_CACHE', {}, clear=True):
            QuantumChemistry.submit_qc(qc, 'partition_job', 8)
        script = Path('partition_job.sbatch').read_text()
        self.assertIn('#SBATCH --partition=day-long-cpu', script)
        self.assertNotIn('#SBATCH -q ', script)
        self.assertIn('#SBATCH --time=1-00:00:00', script)
        self.assertIn('kinbot_exit_file=partition_job.exitcode', script)
        self.assertIn('${SLURM_JOB_ID:-unknown}', script)
        self.assertIn('exit "$kinbot_exit_code"', script)

    def test_completed_scheduler_job_without_result_becomes_error(self):
        qc = SimpleNamespace(queuing='slurm', job_ids={'lost': '123'})
        with patch.object(QuantumChemistry, '_check_qc', return_value=0), \
                patch('kinbot.qc._slurm_job_record', return_value={}), \
                patch('kinbot.qc.time.monotonic', side_effect=[0., 61.]), \
                self.assertLogs('KinBot', 'ERROR'):
            self.assertEqual(QuantumChemistry.check_qc(qc, 'lost'), 'running')
            self.assertEqual(QuantumChemistry.check_qc(qc, 'lost'), 'error')

    def test_scheduler_query_outage_is_bounded_before_submission(self):
        qc = SimpleNamespace(
            queuing='slurm', username='tester', queue_job_limit=8)
        failed = Mock(returncode=1, stdout='', stderr='controller unavailable')
        with patch('kinbot.qc.subprocess.run', return_value=failed), \
                patch('kinbot.qc.time.monotonic', side_effect=[0., 61.]), \
                patch('kinbot.qc.time.sleep'), \
                self.assertRaisesRegex(RuntimeError, 'continuously'):
            QuantumChemistry.limit_jobs(qc)

    def test_queue_limit_counts_only_the_reserved_slurm_branch(self):
        root = str((Path.cwd() / 'anl_automatic').resolve())
        qc = SimpleNamespace(
            queuing='slurm', username='tester', queue_job_limit=2,
            queue_job_root=root, queue_job_root_recursive=True)
        listing = Mock(returncode=0, stderr='', stdout='\n'.join((
            root + '/reference_l2/a',
            str((Path.cwd() / 'vrctst').resolve()),
        )))
        with patch('kinbot.qc.subprocess.run', return_value=listing) as run:
            self.assertEqual(QuantumChemistry.limit_jobs(qc), 0)
        self.assertEqual(
            run.call_args.args[0],
            ['squeue', '-h', '-u', 'tester', '-o', '%Z'])

    def test_automatic_queue_admission_survives_long_controller_outage(self):
        root = str((Path.cwd() / 'anl_automatic').resolve())
        qc = SimpleNamespace(
            queuing='slurm', username='tester', queue_job_limit=2,
            queue_job_root=root, queue_job_root_recursive=True)
        failed = Mock(returncode=1, stdout='',
                      stderr='Unable to contact slurm controller')
        recovered = Mock(returncode=0, stdout='', stderr='')
        with patch('kinbot.qc.subprocess.run', side_effect=[
                failed, failed, recovered]), \
                patch('kinbot.qc.time.monotonic', side_effect=[0., 120.]), \
                patch('kinbot.qc.time.sleep') as sleep:
            self.assertEqual(QuantumChemistry.limit_jobs(qc), 0)
        self.assertEqual(sleep.call_count, 2)

    def test_automatic_reference_submission_uses_cross_process_lock(self):
        qc = SimpleNamespace()
        with patch.object(QuantumChemistry, '_submit_qc', return_value=1) as submit, \
                patch.dict(os.environ, {
                    'KINBOT_QUEUE_LOCK': str(Path.cwd() / 'shared.lock')}):
            self.assertEqual(
                QuantumChemistry.submit_qc(qc, 'job', 8), 1)
        submit.assert_called_once_with(qc, 'job', 8, 1, None)
        self.assertTrue(Path('shared.lock').is_file())

    def test_concurrent_legacy_threads_have_atomic_queue_admission(self):
        active = 0
        maximum = 0
        guard = threading.Lock()

        def submit(*args):
            nonlocal active, maximum
            with guard:
                active += 1
                maximum = max(maximum, active)
            time.sleep(.01)
            with guard:
                active -= 1
            return 1

        qc = SimpleNamespace()
        with patch.object(QuantumChemistry, '_submit_qc', side_effect=submit):
            workers = [threading.Thread(
                target=QuantumChemistry.submit_qc,
                args=(qc, f'job_{index}', 1)) for index in range(8)]
            for worker in workers:
                worker.start()
            for worker in workers:
                worker.join()
        self.assertEqual(maximum, 1)

    def test_completed_worker_reads_result_after_current_slurm_generation(self):
        Path('input.json').write_text(json.dumps({
            'barrier_threshold': 100., 'queuing': 'slurm'}))
        qc = QuantumChemistry(
            Parameters('input.json', show_warnings=False).par)
        qc.qc = 'fc'
        job = 'cached_ts'
        atoms = Atoms('H2', positions=[[0., 0., 0.], [.74, 0., 0.]])
        qc.db.write(atoms, name=job,
                    data={'energy': -1., 'status': 'normal'})
        Path(job + '_sella.log').write_text('done\n')
        Path(job + '.exitcode').write_text('123 0\n')
        qc.job_ids[job] = '123'
        with patch('kinbot.qc._slurm_job_active', return_value=False) as active:
            self.assertEqual(qc.check_qc(job), 'normal')
        active.assert_not_called()

    def test_legacy_profiled_job_keeps_one_field_exit_marker_compatibility(self):
        Path('input.json').write_text(json.dumps({
            'barrier_threshold': 100., 'queuing': 'slurm'}))
        qc = QuantumChemistry(
            Parameters('input.json', show_warnings=False).par)
        qc.qc = 'fc'
        job = 'legacy_completed'
        atoms = Atoms('H2', positions=[[0., 0., 0.], [.74, 0., 0.]])
        qc.db.write(atoms, name=job,
                    data={'energy': -1., 'status': 'normal'})
        Path(job + '_sella.log').write_text('done\n')
        Path(job + '.exitcode').write_text('0\n')
        qc.job_ids[job] = '123'
        qc.job_generation_tagged[job] = False

        with patch('kinbot.qc._slurm_job_active') as active:
            self.assertEqual(qc.check_qc(job), 'normal')
        active.assert_not_called()

    def test_active_replacement_ignores_previous_worker_error(self):
        Path('input.json').write_text(json.dumps({
            'barrier_threshold': 100., 'queuing': 'slurm'}))
        qc = QuantumChemistry(
            Parameters('input.json', show_warnings=False).par)
        qc.qc = 'fc'
        job = 'reaction_path'
        atoms = Atoms('H2', positions=[[0., 0., 0.], [.74, 0., 0.]])
        qc.db.write(atoms, name=job,
                    data={'energy': -1., 'status': 'error'})
        Path(job + '_sella.log').write_text('done\n')
        Path(job + '.exitcode').write_text('1\n')
        qc.job_ids[job] = '456'

        with patch('kinbot.qc._slurm_job_active', return_value=True):
            self.assertEqual(qc.check_qc(job), 'running')

    def test_released_replacement_ignores_previous_generation_marker(self):
        Path('input.json').write_text(json.dumps({
            'barrier_threshold': 100., 'queuing': 'slurm'}))
        qc = QuantumChemistry(
            Parameters('input.json', show_warnings=False).par)
        qc.qc = 'fc'
        job = 'reaction_path'
        atoms = Atoms('H2', positions=[[0., 0., 0.], [.74, 0., 0.]])
        qc.db.write(atoms, name=job,
                    data={'energy': -1., 'status': 'normal'})
        Path(job + '_sella.log').write_text('done\n')
        Path(job + '.exitcode').write_text('455 0\n')
        qc.job_ids[job] = '456'

        with patch('kinbot.qc._slurm_job_active', return_value=False):
            self.assertEqual(qc._check_qc(job), 0)

    def test_current_generation_marker_accepts_completed_worker(self):
        Path('input.json').write_text(json.dumps({
            'barrier_threshold': 100., 'queuing': 'slurm'}))
        qc = QuantumChemistry(
            Parameters('input.json', show_warnings=False).par)
        qc.qc = 'fc'
        job = 'reaction_path'
        atoms = Atoms('H2', positions=[[0., 0., 0.], [.74, 0., 0.]])
        qc.db.write(atoms, name=job,
                    data={'energy': -1., 'status': 'normal'})
        Path(job + '.exitcode').write_text('456 0\n')
        qc.job_ids[job] = '456'

        with patch('kinbot.qc._slurm_job_active', return_value=False):
            self.assertEqual(qc.check_qc(job), 'normal')

    def test_worker_exit_marker_normalizes_pbs_server_suffix(self):
        Path('worker.exitcode').write_text('123.server 0\n')
        self.assertEqual(
            _worker_exit_status('worker.exitcode', '123'), (0, True))
        self.assertIsNone(_worker_exit_status('worker.exitcode', '124'))

    def test_completed_worker_accepts_database_row_without_legacy_done_marker(self):
        Path('input.json').write_text(json.dumps({
            'barrier_threshold': 100., 'queuing': 'slurm'}))
        qc = QuantumChemistry(
            Parameters('input.json', show_warnings=False).par)
        qc.qc = 'fc'
        job = 'restored_high'
        atoms = Atoms('H2', positions=[[0., 0., 0.], [.74, 0., 0.]])
        qc.db.write(atoms, name=job,
                    data={'energy': -1., 'status': 'normal'})
        Path(job + '.exitcode').write_text('0\n')
        # Restored routes can deliberately have no historical scheduler id.
        # The exit marker still proves that the worker/database handoff ended.
        with patch('kinbot.qc.subprocess.run',
                   side_effect=AssertionError('Slurm must not be queried')):
            self.assertEqual(qc.check_qc(job), 'normal')

    def test_completed_worker_without_database_result_is_terminal_error(self):
        Path('input.json').write_text(json.dumps({
            'barrier_threshold': 100., 'queuing': 'slurm'}))
        qc = QuantumChemistry(
            Parameters('input.json', show_warnings=False).par)
        qc.qc = 'fc'
        Path('empty.exitcode').write_text('0\n')
        with self.assertLogs('KinBot', 'ERROR') as messages:
            self.assertEqual(qc.check_qc('empty'), 'error')
        self.assertIn('did not publish a calculation result',
                      '\n'.join(messages.output))

    def test_vrc_scan_completion_waits_for_paired_frozen_geometry(self):
        Path('input.json').write_text(json.dumps({
            'barrier_threshold': 100., 'queuing': 'slurm'}))
        qc = QuantumChemistry(
            Parameters('input.json', show_warnings=False).par)
        qc.qc = 'gauss'
        Path('vrctst').mkdir()
        job = 'vrctst/reaction_vts_pt_asymptote'
        atoms = Atoms('H2', positions=[[0., 0., 0.], [30., 0., 0.]])
        qc.db.write(atoms, name=job,
                    data={'energy': -1., 'status': 'normal'})
        Path(job + '.log').write_text('done\n')

        self.assertIsNone(qc._database_qc_status(job))

        qc.db.write(atoms, name=job + '_fr',
                    data={'energy': -1., 'status': 'normal'})
        Path(job + '_fr.log').write_text('done\n')
        self.assertEqual(qc._database_qc_status(job), 'normal')

    def test_finished_vrc_worker_without_frozen_geometry_is_error(self):
        Path('input.json').write_text(json.dumps({
            'barrier_threshold': 100., 'queuing': 'slurm'}))
        qc = QuantumChemistry(
            Parameters('input.json', show_warnings=False).par)
        qc.qc = 'gauss'
        Path('vrctst').mkdir()
        job = 'vrctst/reaction_vts_pt_asymptote'
        atoms = Atoms('H2', positions=[[0., 0., 0.], [30., 0., 0.]])
        qc.db.write(atoms, name=job,
                    data={'energy': -1., 'status': 'normal'})
        Path(job + '.exitcode').write_text('0\n')
        with self.assertLogs('KinBot', 'ERROR'):
            self.assertEqual(qc._check_qc(job), 'error')

    def test_slurm_status_query_is_targeted_and_bounded(self):
        result = Mock(returncode=0, stdout='123\n', stderr='')
        with patch('kinbot.qc.subprocess.run', return_value=result) as run:
            self.assertTrue(_slurm_job_active('123'))
        self.assertEqual(
            run.call_args.args[0], ['squeue', '-h', '-j', '123', '-o', '%i'])
        self.assertEqual(run.call_args.kwargs['timeout'], 15)

    def test_slurm_invalid_job_id_is_terminal_not_an_outage(self):
        result = Mock(
            returncode=1, stdout='',
            stderr='slurm_load_jobs error: Invalid job id specified')
        with patch('kinbot.qc.subprocess.run', return_value=result):
            self.assertFalse(_slurm_job_active('52906599'))

    def test_slurm_status_timeout_remains_running(self):
        qc = SimpleNamespace(
            queuing='slurm', job_ids={'pending': '123'},
            ingest_pkl=lambda job: None,
            is_in_database=lambda job: False)
        with patch('kinbot.qc._slurm_job_active', return_value=None), \
                self.assertLogs('KinBot', 'WARNING') as messages:
            self.assertEqual(
                QuantumChemistry._check_qc(qc, 'pending'), 'running')
            self.assertEqual(
                QuantumChemistry._check_qc(qc, 'pending'), 'running')
        self.assertEqual(sum('Slurm status query is unavailable' in message
                             for message in messages.output), 1)

    def test_scheduler_submission_timeout_is_bounded(self):
        qc = SimpleNamespace(
            check_qc=lambda job: 0, queue_job_limit=0,
            par={'queue_template': ''}, queuing='slurm', qc='fc',
            queue_name='day-long-cpu', slurm_feature='', job_ids={},
        )
        process = Mock()
        process.communicate.side_effect = [
            subprocess.TimeoutExpired(['sbatch'], 30), (b'', b'')]
        process.kill = Mock()
        discovery = Mock(returncode=0, stdout='1-00:00:00\n', stderr='')
        with patch('kinbot.qc.subprocess.Popen', return_value=process), \
                patch('kinbot.qc.subprocess.run', return_value=discovery), \
                patch.dict('kinbot.qc._SLURM_WALLTIME_CACHE', {}, clear=True):
            with self.assertRaisesRegex(RuntimeError,
                                        'submission timed out'):
                QuantumChemistry.submit_qc(qc, 'timeout_job', 1)
        process.kill.assert_called_once_with()

    def test_automatic_submission_retries_explicit_controller_failure(self):
        root = str((Path.cwd() / 'anl_automatic').resolve())
        qc = SimpleNamespace(
            check_qc=lambda job: 0, queue_job_limit=0,
            queue_job_root=root,
            par={'queue_template': ''}, queuing='slurm', qc='fc',
            queue_name='day-long-cpu', slurm_feature='', job_ids={})
        unavailable = Mock(returncode=1)
        unavailable.communicate.return_value = (
            b'', b'sbatch: error: Unable to contact slurm controller')
        accepted = Mock(returncode=0)
        accepted.communicate.return_value = (
            b'Submitted batch job 123\n', b'')
        discovery = Mock(returncode=0, stdout='1-00:00:00\n', stderr='')
        with patch('kinbot.qc.subprocess.Popen',
                   side_effect=[unavailable, accepted]) as submit, \
                patch('kinbot.qc.subprocess.run', return_value=discovery), \
                patch('kinbot.qc.time.sleep') as sleep, \
                patch.dict('kinbot.qc._SLURM_WALLTIME_CACHE', {}, clear=True):
            self.assertEqual(
                QuantumChemistry.submit_qc(qc, 'retry_job', 1), 1)
        self.assertEqual(submit.call_count, 2)
        sleep.assert_called_once_with(5)
        self.assertEqual(qc.job_ids['retry_job'], '123')

    def test_automatic_submission_retries_transient_process_launch(self):
        root = str((Path.cwd() / 'anl_automatic').resolve())
        qc = SimpleNamespace(
            check_qc=lambda job: 0, queue_job_limit=0,
            queue_job_root=root,
            par={'queue_template': ''}, queuing='slurm', qc='fc',
            queue_name='day-long-cpu', slurm_feature='', job_ids={})
        accepted = Mock(returncode=0)
        accepted.communicate.return_value = (
            b'Submitted batch job 123\n', b'')
        discovery = Mock(returncode=0, stdout='1-00:00:00\n', stderr='')
        with patch('kinbot.qc.subprocess.Popen', side_effect=[
                OSError(11, 'Resource temporarily unavailable'),
                accepted]) as submit, \
                patch('kinbot.qc.subprocess.run', return_value=discovery), \
                patch('kinbot.qc.time.sleep') as sleep, \
                patch.dict('kinbot.qc._SLURM_WALLTIME_CACHE', {}, clear=True):
            self.assertEqual(
                QuantumChemistry.submit_qc(qc, 'retry_launch', 1), 1)
        self.assertEqual(submit.call_count, 2)
        sleep.assert_called_once_with(5)

    def test_unrecognized_scheduler_success_output_is_hard_error(self):
        qc = SimpleNamespace(
            check_qc=lambda job: 0, queue_job_limit=0,
            par={'queue_template': ''}, queuing='slurm', qc='fc',
            queue_name='day-long-cpu', slurm_feature='', job_ids={})
        process = Mock(returncode=0)
        process.communicate.return_value = (b'unexpected success\n', b'')
        discovery = Mock(returncode=0, stdout='1-00:00:00\n', stderr='')
        with patch('kinbot.qc.subprocess.Popen', return_value=process), \
                patch('kinbot.qc.subprocess.run', return_value=discovery), \
                patch.dict('kinbot.qc._SLURM_WALLTIME_CACHE', {}, clear=True), \
                self.assertRaisesRegex(RuntimeError,
                                       'Cannot parse scheduler job id'):
            QuantumChemistry.submit_qc(qc, 'bad_output', 1)

    def test_slurm_io_launch_failure_is_immediate_and_excludes_node(self):
        qc = SimpleNamespace(queuing='slurm', job_ids={'lost': '123'},
                             _slurm_excluded_nodes=set())
        record = {
            'JobState': 'FAILED',
            'Reason': 'RaisedSignal:53(Real-time_signal_19)',
            'ExitCode': '0:53',
            'NodeList': 'node16',
        }
        with patch.object(QuantumChemistry, '_check_qc', return_value=0), \
                patch('kinbot.qc._slurm_job_record', return_value=record), \
                patch.dict(os.environ, {}, clear=True), \
                self.assertLogs('KinBot', 'ERROR') as messages:
            self.assertEqual(QuantumChemistry.check_qc(qc, 'lost'), 'error')
            self.assertEqual(os.environ['SBATCH_EXCLUDE'], 'node16')
        self.assertIn('never started', '\n'.join(messages.output))
        self.assertEqual(qc._slurm_excluded_nodes, {'node16'})
        saved = json.loads(Path('.kinbot_slurm_nodes.json').read_text())
        self.assertEqual(saved['excluded_nodes'], ['node16'])

    def test_slurm_node_exclusion_preserves_and_reaches_nested_submitters(self):
        with patch.dict(os.environ, {'SBATCH_EXCLUDE': 'node15'}, clear=True):
            _export_slurm_excluded_nodes({'node16', 'node15'})
            self.assertEqual(os.environ['SBATCH_EXCLUDE'], 'node15,node16')

    def test_recorded_slurm_node_is_excluded_from_later_jobs(self):
        qc = SimpleNamespace(
            check_qc=lambda job: 0, queue_job_limit=0,
            par={'queue_template': ''}, queuing='slurm', qc='fc',
            queue_name='day-long-cpu', slurm_feature='', job_ids={},
            _slurm_excluded_nodes={'node16'},
        )
        process = Mock()
        process.communicate.return_value = (b'Submitted batch job 124\n', b'')
        discovery = Mock(returncode=0, stdout='1-00:00:00\n', stderr='')
        with patch('kinbot.qc.subprocess.Popen', return_value=process), \
                patch('kinbot.qc.subprocess.run', return_value=discovery), \
                patch.dict('kinbot.qc._SLURM_WALLTIME_CACHE', {}, clear=True):
            QuantumChemistry.submit_qc(qc, 'replacement', 8)
        script = Path('replacement.sbatch').read_text()
        self.assertIn('#SBATCH --exclude=node16', script)

    def test_nonfinite_geometry_is_never_a_normal_qc_result(self):
        Path('input.json').write_text(json.dumps({
            'barrier_threshold': 100., 'queuing': 'local'}))
        qc = QuantumChemistry(
            Parameters('input.json', show_warnings=False).par)
        atoms = Atoms('H2', positions=[[0., 0., 0.], [np.nan] * 3])
        qc.db.write(atoms, name='bad',
                    data={'energy': np.nan, 'status': 'normal'})
        Path('bad.log').write_text('done\n')
        with self.assertLogs('KinBot', 'ERROR'):
            self.assertEqual(qc._check_qc('bad'), 'error')
        self.assertEqual(too_far(atoms.positions), -1)

    def test_fairchem_cleanup_preserves_unrelated_driver_log(self):
        Path('driver.log').write_text('')
        Path('reaction_sella.log').write_text('')
        Path('pending.pkl').write_bytes(b'pending')
        clean_files(qc='fc')
        self.assertTrue(Path('driver.log').is_file())
        self.assertFalse(Path('reaction_sella.log').exists())
        self.assertTrue(Path('pending.pkl').exists())

        Path('diagnostic.pkl').write_bytes(b'pending')
        clean_files(qc='fc', diagnostic=True)
        self.assertTrue(Path('diagnostic.pkl').is_file())

    def test_replacement_submission_removes_only_its_stale_handoff(self):
        Path('replacement.pkl').write_bytes(b'stale')
        Path('unrelated.pkl').write_bytes(b'pending')
        qc = SimpleNamespace(
            check_qc=lambda job: 0, queue_job_limit=0,
            par={'queue_template': ''}, queuing='slurm', qc='fc',
            queue_name='day-long-cpu', slurm_feature='', job_ids={},
        )
        process = Mock()
        process.communicate.return_value = (b'Submitted batch job 123\n', b'')
        discovery = Mock(returncode=0, stdout='1-00:00:00\n', stderr='')
        with patch('kinbot.qc.subprocess.Popen', return_value=process), \
                patch('kinbot.qc.subprocess.run', return_value=discovery), \
                patch.dict('kinbot.qc._SLURM_WALLTIME_CACHE', {}, clear=True):
            QuantumChemistry.submit_qc(qc, 'replacement', 1)
        self.assertFalse(Path('replacement.pkl').exists())
        self.assertTrue(Path('unrelated.pkl').exists())


class TestVRCJobs(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        Path('vrctst').mkdir()
        self.qc = QuantumChemistry.__new__(QuantumChemistry)
        self.qc.qc, self.qc.qc_command, self.qc.ppn = 'gauss', 'g16', 8
        self.qc.par = {'calc_kwargs': {}, 'vrc_tst_scan_sella': False,
                       'vrc_tst_scan_deviation': 0.1, 'vts_ang_dev': 10,
                       'sella_kwargs': {}}
        self.qc.get_qc_arguments = Mock(return_value={})
        self.qc.submit_qc = Mock()
        self.fragment = StationaryPoint('H', 0, 2, atom=['H'], geom=np.array([[0., 0., 0.]]))
        self.fragment.characterize()

    def test_fragment_job_uses_fragment_electron_count(self):
        job = QuantumChemistry.qc_vts_frag(self.qc, self.fragment)
        self.qc.get_qc_arguments.assert_called_once_with(job, 2, 0, 1, vts=1)
        self.qc.submit_qc.assert_called_once_with(job, 1)
        ast.parse(Path(f'{job}.py').read_text())

    def test_fragment_job_keeps_legacy_calc_kwargs_fallback(self):
        self.qc.par['calc_kwargs'] = {'scf': 'xqc'}
        self.qc.merge_kwargs = Mock(return_value={'scf': 'xqc'})
        QuantumChemistry.qc_vts_frag(self.qc, self.fragment)
        self.qc.merge_kwargs.assert_called_once()
        self.assertEqual(self.qc.merge_kwargs.call_args.args[1], {'scf': 'xqc'})

    def test_fragment_job_prefers_explicit_vrc_kwargs(self):
        self.qc.par['calc_kwargs'] = {'scf': 'legacy'}
        self.qc.par['vrc_tst_calc_kwargs'] = {'scf': 'vrc'}
        self.qc.merge_kwargs = Mock(return_value={'scf': 'vrc'})
        QuantumChemistry.qc_vts_frag(self.qc, self.fragment)
        self.qc.merge_kwargs.assert_called_once()
        self.assertEqual(self.qc.merge_kwargs.call_args.args[1], {'scf': 'vrc'})

    def test_explicit_empty_vrc_kwargs_disable_legacy_fallback(self):
        self.qc.par['calc_kwargs'] = {'scf': 'legacy'}
        self.qc.par['vrc_tst_calc_kwargs'] = {}
        self.qc.merge_kwargs = Mock()
        QuantumChemistry.qc_vts_frag(self.qc, self.fragment)
        self.qc.merge_kwargs.assert_not_called()

    def test_scan_jobs_use_reacting_species_electron_count(self):
        species = StationaryPoint('H2', 0, 1, atom=['H', 'H'],
                                  geom=np.array([[0., 0., 0.], [.74, 0., 0.]]))
        species.characterize()
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
                    self.qc, reaction, species.geom, 0, [], asymptote, species.geom.tolist())
                self.qc.get_qc_arguments.assert_called_once_with(job, 1, 0, 2, vts=1)
                self.qc.submit_qc.assert_called_once_with(job, 2)
                ast.parse(Path(f'{job}.py').read_text())
                script = Path(f'{job}.py').read_text()
                frozen_done = script.index("with open(logfile, 'a') as f:")
                relaxed_done = script.rindex(
                    f"with open('{job}.log', 'a') as f:")
                self.assertGreater(relaxed_done, frozen_done)


if __name__ == '__main__':
    unittest.main()
