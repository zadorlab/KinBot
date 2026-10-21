"""Generated-input checks for ordinary unimolecular job writers.

References were compared with master 684e513c. The Gaussian IRC
worker also contains an inactive bimolecular branch; ordinary runs retain
master's IRC and endpoint optimization steps.
The temporary working directory is replaced with ``{WORKDIR}``.
"""

import hashlib
import json
import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

from ase.build import molecule

from kinbot.irc import IRC
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint


REFERENCE_HASHES = {
    'gaussian_well':
        'f7f1604a366f9c3389ac9938c33171d0c350fe9cca9c4a9a394b0643c3c2a39a',
    'gaussian_ts':
        'a69025662e6e359e587635315f03ee68a13a2c25000f875069a10a97372bc3d7',
    'gaussian_conformer':
        '3535993f708d4024cee4eab0dea6501133414783599f752c2782dee1d7f12a5b',
    'gaussian_ring':
        'a6ea90d50e762dfa4de55fe4ab42e491301e2a6b96489153f4a5d869c2301b1b',
    # The generated worker selects the ordinary branch with a literal False.
    'gaussian_irc':
        '8a71b090d699ed13d2bfbbb6dbf16ec26c5cf93e417f7fb328bc64368d653c44',
    # Shared selected-Hessian fix requests VIBMAN_PRINT=4 at L1. Upstream
    # now imports KinBot's Q-Chem wrapper rather than the original ASE class;
    # both reference hashes below use that current import.
    'qchem_conformer':
        'a421526c67eef835d50057ea368e6f539519d86195a748cb1fe9fd4fb3b0dc07',
    'qchem_conformer_before_hessian_print':
        'b9a349f25f7421a87869ab825b15f53423d6117a823e3f2303571ec4cd131816',
    # Includes the theory branch's reporting of calculator startup failures.
    'orca_conformer':
        'bed0b68ddd07db32263ba75f64da010402e178d78a6de415a1756fad5fc9b6e2',
}


class TestUnimolecularGeneratedInputs(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.previous = Path.cwd()
        self.workdir = Path(self.temporary.name)
        os.chdir(self.workdir)
        Path('conf').mkdir()

    def tearDown(self):
        os.chdir(self.previous)
        self.temporary.cleanup()

    def parameters(self, **updates):
        path = self.workdir / 'parameters.json'
        path.write_text(json.dumps({'barrier_threshold': 50., **updates}))
        return Parameters(path).par

    @staticmethod
    def species(name='reference', wellorts=0):
        atoms = molecule('CH4')
        point = StationaryPoint(
            name, 0, 1, atom=atoms.get_chemical_symbols(),
            geom=atoms.positions, wellorts=wellorts)
        point.characterize()
        return point

    def digest(self, path):
        value = Path(path).read_text().replace(
            str(self.workdir), '{WORKDIR}')
        return hashlib.sha256(value.encode()).hexdigest()

    @staticmethod
    def disable_submission(qc):
        qc.submit_qc = lambda *args, **kwargs: 0
        for backend in getattr(qc, 'backends', {}).values():
            backend.submit_qc = lambda *args, **kwargs: 0
        return qc

    def test_default_gaussian_jobs_keep_the_reviewed_writer_changes(self):
        par = self.parameters()
        qc = self.disable_submission(QuantumChemistry(par))
        well = self.species()
        ts = self.species('reference_ts', wellorts=1)
        qc.qc_opt(well, well.geom)
        qc.qc_opt_ts(ts, ts.geom)
        qc.qc_conf(well, well.geom, 1)
        qc.qc_ring_conf(
            well, well.geom, [[1, 2, 3, 4]], [60.], 0, 0)
        qc.get_qc_geom = lambda *args, **kwargs: (0, ts.geom)
        reaction = SimpleNamespace(
            qc=qc, species=ts, instance_name='reference_ts', par=par)
        IRC(reaction, par).do_irc_calculations()

        files = {
            'gaussian_well': f'{well.chemid}_well.py',
            'gaussian_ts': 'reference_ts.py',
            'gaussian_conformer': f'conf/{well.chemid}_0001.py',
            'gaussian_ring': f'conf/{well.chemid}_r0000_0000.py',
            'gaussian_irc': 'reference_ts_IRC_F.py',
        }
        for label, filename in files.items():
            with self.subTest(label=label):
                self.assertEqual(self.digest(filename),
                                 REFERENCE_HASHES[label])

    def test_other_backend_conformers_match_master(self):
        for backend, command in (
                ('qchem', 'qchem'),
                ('orca', 'orca')):
            with self.subTest(backend=backend):
                par = self.parameters(qc=backend, qc_command=command)
                qc = self.disable_submission(QuantumChemistry(par))
                point = self.species(backend)
                qc.qc_conf(point, point.geom, 2)
                self.assertEqual(
                    self.digest(f'conf/{point.chemid}_0002.py'),
                    REFERENCE_HASHES[f'{backend}_conformer'])
                if backend == 'qchem':
                    generated = Path(f'conf/{point.chemid}_0002.py').read_text().replace(
                        str(self.workdir), '{WORKDIR}')
                    self.assertIn("'vibman_print': '4'", generated)
                    original = generated.replace("'vibman_print': '4', ", '')
                    self.assertEqual(hashlib.sha256(original.encode()).hexdigest(),
                                     REFERENCE_HASHES['qchem_conformer_before_hessian_print'])

    def test_legacy_json_names_map_to_l0_names(self):
        par = self.parameters(
            semi_emp_conformer_search=1, semi_emp_method='pm7')
        self.assertEqual(par['l0_conformer_search'], 1)
        self.assertEqual(par['l0_method'], 'pm7')

    def test_calculator_kwargs_remain_level_specific(self):
        par = self.parameters(
            calc_kwargs={'mem': '1GB'}, l0_calc_kwargs={'mem': '2GB'})
        qc = self.disable_submission(QuantumChemistry(par))
        point = self.species()
        qc.qc_conf(point, point.geom, 3)
        qc.qc_conf(point, point.geom, 4, semi_emp=1)
        qc.qc_ring_conf(
            point, point.geom, [[1, 2, 3, 4]], [60.], 1, 0)

        l1 = Path(f'conf/{point.chemid}_0003.py').read_text()
        l0 = Path(f'conf/{point.chemid}_semi_emp_0004.py').read_text()
        ring = Path(f'conf/{point.chemid}_r0001_0000.py').read_text()
        self.assertIn("'mem': '700MW,1GB'", l1)
        self.assertNotIn("'mem': '700MW,2GB'", l1)
        for generated in (l0, ring):
            self.assertIn("'mem': '700MW,2GB'", generated)
            self.assertNotIn("'mem': '700MW,1GB'", generated)


if __name__ == '__main__':
    unittest.main()
