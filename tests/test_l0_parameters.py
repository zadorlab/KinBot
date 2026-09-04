"""L0 input names preserve existing conformer-search settings and jobs."""

import ast
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from ase.build import molecule

from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint


class TestL0Parameters(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        Path('conf').mkdir()

    @staticmethod
    def parameters(**options):
        path = Path('input.json')
        path.write_text(json.dumps({'barrier_threshold': 50., **options}))
        return Parameters(path).par

    def test_defaults_keep_am1_and_disable_l0_search(self):
        par = self.parameters()
        self.assertEqual(par['l0_method'], 'am1')
        self.assertEqual(par['L0_conformer_search'], 0)
        self.assertNotIn('semi_emp_method', par)
        self.assertNotIn('semi_emp_conformer_search', par)

    def test_old_and_new_names_produce_the_same_parameters(self):
        old = self.parameters(semi_emp_method='pm6', semi_emp_conformer_search=1)
        new = self.parameters(l0_method='pm6', L0_conformer_search=1)
        self.assertEqual(old, new)
        both = self.parameters(semi_emp_method='pm6', l0_method='pm6',
                               semi_emp_conformer_search=1, L0_conformer_search=1)
        self.assertEqual(both, new)

    def test_conflicting_aliases_are_rejected(self):
        for options in [{'semi_emp_method': 'pm6', 'l0_method': 'am1'},
                        {'semi_emp_conformer_search': 1, 'L0_conformer_search': 0}]:
            with self.subTest(options=options), self.assertRaisesRegex(
                    IOError, 'specify different values'):
                self.parameters(**options)

    def test_aliases_write_identical_l0_jobs_without_changing_l1(self):
        atoms = molecule('CH4')
        species = StationaryPoint('methane', 0, 1,
                                  atom=atoms.get_chemical_symbols(), geom=atoms.positions)
        species.characterize()
        scripts = []
        for names in [{'semi_emp_method': 'pm6', 'semi_emp_conformer_search': 1},
                      {'l0_method': 'pm6', 'L0_conformer_search': 1}]:
            qc = QuantumChemistry(self.parameters(**names))
            qc.submit_qc = lambda *args, **kwargs: 0
            jobs = []
            for is_l0, prefix, method in [(1, 'semi_emp_', 'pm6'), (0, '', 'b3lyp')]:
                qc.qc_conf(species, species.geom, 1, semi_emp=is_l0)
                script = Path(f'conf/{species.chemid}_{prefix}0001.py').read_text()
                self.assertIn(f"'method': '{method}'", script)
                ast.parse(script)
                jobs.append(script)
            scripts.append(jobs)
        self.assertEqual(scripts[0], scripts[1])


if __name__ == '__main__':
    unittest.main()
