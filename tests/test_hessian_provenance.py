"""Check that native QChem and ASE Hessians produce the same frequencies."""

import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from ase.build import molecule
import numpy as np

from kinbot import constants, frequencies
from kinbot.optimize import Optimize
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint


class TestHessianWeighting(unittest.TestCase):
    def test_native_and_sella_qchem_results_project_identically(self):
        with TemporaryDirectory() as directory:
            previous = Path.cwd()
            try:
                os.chdir(directory)
                atoms = molecule('H2O')
                point = StationaryPoint('water', 0, 1,
                                        atom=atoms.get_chemical_symbols(), geom=atoms.positions)
                point.characterize()
                hessian = np.diag(np.arange(1., 10.))
                masses = np.repeat([constants.exact_mass[atom] for atom in point.atom], 3)
                weighted = hessian / np.sqrt(np.outer(masses, masses))
                with open('audit_freq.out', 'w') as handle:
                    handle.write(' Mass-Weighted Hessian Matrix\n')
                    np.savetxt(handle, weighted)
                    handle.write(' Translations and Rotations\n')
                expected, _ = frequencies.get_frequencies(point, hessian, point.geom)

                for use_sella in (False, True):
                    with self.subTest(use_sella=use_sella):
                        qc = QuantumChemistry.__new__(QuantumChemistry)
                        qc.qc = 'qchem'
                        qc.use_sella = use_sella
                        qc.check_qc = lambda job: 'normal'
                        qc.db = SimpleNamespace(select=lambda **kw: [
                            SimpleNamespace(data={'hess': hessian})])
                        optimization = Optimize.__new__(Optimize)
                        optimization.species = point
                        optimization.qc = qc
                        optimization.par = {
                            'conformer_search': 0, 'rotation_restart': -1,
                            'high_level': 0, 'rotor_scan': 1,
                            'multi_conf_tst': 0, 'L3_calc': 0,
                        }
                        optimization.shir = 1
                        optimization.restart = 0
                        optimization.just_high = False
                        optimization.defer_hir = False
                        optimization.wait = 0
                        optimization.log_name = lambda level: 'audit'
                        optimization.do_optimization()
                        np.testing.assert_allclose(point.reduced_freqs, expected, rtol=1.e-8)
            finally:
                os.chdir(previous)


if __name__ == '__main__':
    unittest.main()
