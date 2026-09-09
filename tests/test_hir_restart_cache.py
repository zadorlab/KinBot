"""HIR restarts must submit new work instead of accepting a backend's old cache."""
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock, patch

import numpy as np
from ase import Atoms

from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry


class TestHIRRestartCache(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(directory.name)
        Path('hir').mkdir()
        Path('input.json').write_text(json.dumps({'barrier_threshold': 100., 'queuing': 'pbs'}))
        self.qc = QuantumChemistry(Parameters('input.json', show_warnings=False).par)
        self.atoms = Atoms('HH', positions=[[0., 0., 0.], [1., 0., 0.]])

    def test_every_backend_resubmits_invalidated_hir_and_l2_results(self):
        process = Mock()
        process.communicate.return_value = (b'123.server\n', b'')
        # No live scheduler/job is contacted. Actual status and submit methods run.
        with patch('kinbot.qc.subprocess.call', return_value=1), \
             patch('kinbot.qc.subprocess.Popen', return_value=process), \
             patch('kinbot.qc.time.sleep'):
            for backend, suffix in (('gauss', '.log'), ('qchem', '.out'),
                                    ('nwchem', '.out'), ('orca', '_sella.log'),
                                    ('fc', '_sella.log'), ('nn_pes', None)):
                for scan in (True, False):
                    with self.subTest(backend=backend, scan=scan):
                        self.qc.qc = backend
                        job = f'hir/{backend}_hir_0_00' if scan else f'{backend}_well_high'
                        old = self.qc.db.write(self.atoms, name=job, data={
                            'energy': -1., 'zpe': .01, 'frequencies': [100.],
                            'hess': np.eye(6), 'status': 'normal'})
                        if suffix:
                            Path(job + suffix).write_text('old normal output\ndone\n')
                        Path(job + '.fchk').write_text('old Hessian')
                        Path(job + '.chk').write_text('old checkpoint')
                        self.assertEqual(self.qc.check_qc(job), 'normal')
                        self.qc.invalidate_qc(job)
                        self.assertEqual(self.qc.check_qc(job), 0)
                        self.assertFalse(Path(job + '.fchk').exists())
                        self.assertFalse(Path(job + '.chk').exists())
                        self.assertEqual(self.qc.db.get(id=old).data.energy, -1.)
                        archives = list(Path(job).parent.glob(Path(job).name + '.fchk.restart_*'))
                        self.assertEqual(len(archives), 1)
                        self.assertEqual(archives[0].read_text(), 'old Hessian')
                        self.assertEqual(self.qc.submit_qc(job, 2), 1)
                        process.communicate.assert_called()
                        # A new completed calculation supersedes the invalidation marker.
                        self.qc.db.write(self.atoms, name=job, data={
                            'energy': -2., 'zpe': .02, 'frequencies': [200.],
                            'hess': np.eye(6) * 2, 'status': 'normal'})
                        if backend == 'qchem':
                            with open(job + '_freq.out', 'w') as handle:
                                handle.write(' Mass-Weighted Hessian Matrix\n')
                                np.savetxt(handle, np.eye(6) * 2)
                                handle.write(' Translations and Rotations\n')
                        if suffix:
                            Path(job + suffix).write_text('new normal output\ndone\n')
                        self.assertEqual(self.qc.check_qc(job), 'normal')
                        np.testing.assert_allclose(self.qc.read_qc_hess(job, 2), np.eye(6) * 2)
                        self.assertEqual(self.qc.submit_qc(job, 2), 0)

    def test_running_job_is_preserved(self):
        with patch.object(self.qc, 'check_qc', return_value='running'):
            with self.assertRaisesRegex(ValueError, 'running'):
                self.qc.invalidate_qc('still_running')
        self.assertEqual(self.qc.db.count(), 0)


if __name__ == '__main__':
    unittest.main()
