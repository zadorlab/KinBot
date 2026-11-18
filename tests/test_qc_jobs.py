"""Job-writing regressions with all scheduler submission mocked out."""

import ast
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch
import numpy as np
import pytest

from kinbot.qc import QuantumChemistry
from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint


@pytest.mark.parametrize('calculator', ['gauss', 'qchem', 'orca', 'fc'])
def test_isolated_proton_needs_no_worker_at_any_geometry_level(tmp_path, monkeypatch, calculator):
    monkeypatch.chdir(tmp_path)
    Path('input.json').write_text(json.dumps({'qc': calculator, 'queuing': 'local',
                                              'barrier_threshold': 100., 'ppn': 8, 'fc_model_path': 'uma-s-1p1'}))
    qc = QuantumChemistry(Parameters('input.json').par)
    qc.submit_qc = Mock(side_effect=AssertionError('submitted a proton calculation'))
    proton = StationaryPoint('proton', 1, 1, atom=['H'], geom=np.array([[1., 2., 3.]]))
    proton.characterize()
    name = str(proton.chemid)
    for level, suffix in ((0, '_product'), (1, '_well_high')):
        job = name + suffix
        for restart in range(2):
            qc.qc_opt(proton, proton.geom, initial_product=not level, high_level=level)
            assert qc.get_qc_energy(job) == (0, 0.)
        assert qc.db.count(name=job) == 1
        assert qc.get_qc_zpe(job) == (0, 0.)
        assert qc.get_qc_freq(job, 1) == (0, [])
        np.testing.assert_equal(qc.get_qc_geom(job, 1)[1], proton.geom)
        np.testing.assert_equal(qc.read_qc_hess(job, 1), np.zeros((3, 3)))
        assert not Path(job + '.py').exists()
    qc.qc_conf(proton, proton.geom, 0, semi_emp=1)
    assert qc.get_qc_energy(f'conf/{name}_semi_emp_0000') == (0, 0.)
    qc.qc_opt(proton, proton.geom, ext='_custom')
    assert qc.check_qc(name + '_custom') == 'normal'
    qc.submit_qc.assert_not_called()


def test_proton_result_survives_older_worker_and_rejects_hydride(tmp_path, monkeypatch):
    import pickle
    from kinbot.stereo_routing import StereoRoutingError

    monkeypatch.chdir(tmp_path)
    Path('input.json').write_text(json.dumps({'qc': 'fc', 'queuing': 'local',
                                              'barrier_threshold': 100., 'fc_model_path': 'uma-s-1p1'}))
    qc = QuantumChemistry(Parameters('input.json').par)
    qc.submit_qc = Mock(side_effect=AssertionError('unexpected calculation'))
    proton = StationaryPoint('proton', 1, 1, atom=['H'], geom=np.zeros((1, 3)))
    proton.characterize()
    hydride = StationaryPoint('hydride', -1, 1, atom=['H'], geom=proton.geom)
    hydride.characterize()
    assert proton.chemid == hydride.chemid
    job = f'conf/{proton.chemid}_0000'
    Path('conf').mkdir()
    Path(job + '.pkl').write_bytes(pickle.dumps(dict(sym=['H'], pos=proton.geom,
        name=job, data={'status': 'error'})))
    qc.qc_conf(proton, proton.geom, 0)
    assert qc.get_qc_energy(job) == (0, 0.)
    assert qc.qc_freq(proton, job) == job
    for call in (lambda: qc.qc_conf(hydride, hydride.geom, 0),
                 lambda: qc.qc_freq(hydride, job)):
        with pytest.raises(StereoRoutingError, match='charge or multiplicity'):
            call()
    assert qc.db.count(name=job) == 1
    qc.submit_qc.assert_not_called()


def test_hydrogen_ionization_submits_only_the_neutral_atom(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('aie').mkdir()
    Path('input.json').write_text(json.dumps({'qc': 'gauss', 'queuing': 'local',
                                              'barrier_threshold': 100.}))
    qc = QuantumChemistry(Parameters('input.json').par)
    qc.submit_qc = Mock()
    atom = StationaryPoint('hydrogen', 0, 2, atom=['H'], geom=np.zeros((1, 3)))
    atom.characterize()
    qc.qc_aie(atom, atom.geom, '0000')
    base = f'aie/{atom.chemid}_AIE'
    qc.submit_qc.assert_called_once_with(base + '0_0000', 1)
    assert qc.get_qc_energy(base + '1_0000') == (0, 0.)
    assert not Path(base + '1_0000.py').exists()


class TestSchedulerJobs(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)

    def test_queue_limit_counts_jobs_instead_of_empty_lines(self):
        for queue in ('slurm', 'pbs'):
            for empty in (b'', b'\n', b'  \n\n'):
                with self.subTest(queue=queue, empty=empty):
                    qc = SimpleNamespace(queuing=queue, username='guest', queue_job_limit=1)
                    with patch('kinbot.qc.subprocess.check_output', return_value=empty), \
                            patch('kinbot.qc.time.sleep') as sleep:
                        self.assertEqual(QuantumChemistry.limit_jobs(qc), 0)
                    sleep.assert_not_called()
            for limit, occupied in ((1, b'123 job\n'), (1, b'123 job'),
                                    (2, b'123 job\n456 job\n')):
                with self.subTest(queue=queue, limit=limit, occupied=occupied):
                    qc = SimpleNamespace(queuing=queue, username='guest', queue_job_limit=limit)
                    with patch('kinbot.qc.subprocess.check_output', side_effect=[occupied, b'']), \
                            patch('kinbot.qc.time.sleep') as sleep:
                        self.assertEqual(QuantumChemistry.limit_jobs(qc), 0)
                    sleep.assert_called_once_with(1)

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
                    db=SimpleNamespace(select=lambda **kwargs: iter([object()])),
                )
                process = Mock()
                process.communicate.return_value = (b'123.server\n', b'')
                with patch('kinbot.qc.subprocess.Popen', return_value=process) as submit:
                    result = QuantumChemistry.submit_qc(qc, 'test_job', requested)
                script = Path('test_job.pbs').read_text()
                self.assertIn(f'nodes=1:ppn={expected}', script)
                self.assertIn(f'OMP_NUM_THREADS={expected}', script)
                self.assertEqual(qc.job_ids, {'test_job': '123'})
                self.assertEqual(qc._rows_at_submit, {'test_job': 1})
                self.assertEqual(result, 1)
                self.assertEqual(submit.call_args.args[0], ['qsub', 'test_job.pbs'])


class TestVRCJobs(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        Path('vrctst').mkdir()
        self.qc = QuantumChemistry.__new__(QuantumChemistry)
        self.qc.qc, self.qc.qc_command, self.qc.ppn = 'gauss', 'g16', 8
        self.qc.theory_level = None
        self.qc.frequency_mode = 'auto'
        self.qc.optimizer = 'native'
        self.qc.profile_explicit = False
        self.qc.use_sella = False
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


if __name__ == '__main__':
    unittest.main()
