"""Sella job scripts for Gaussian keep the job's path and need a checkpoint for Guess=Read.

Conformer and rotor jobs are named conf/... and hir/...; KinBot reads their
done stamp from {job}_sella.log. The templates used os.path.basename for every
code, which put the Sella log of such a job in the working directory, where
nothing looks for it. ORCA keeps the base name: its script runs in a scratch
directory whose contents are copied back into the job's directory.

The templates also switched Gaussian to Guess=Read after the first energy even
when the chk parameter had been removed, which stops Gaussian in l301.
"""

import os
from textwrap import dedent
from types import SimpleNamespace
import unittest

from kinbot import kb_path


TEMPLATES = ('ase_sella_opt_well', 'ase_sella_hir', 'ase_sella_ts_end')


class FakeCalculator:
    instances = []

    def __init__(self, **kwargs):
        self.kwargs = dict(kwargs)
        FakeCalculator.instances.append(self)

    def get_potential_energy(self):
        return 0.


def render(template, code, label):
    text = open(os.path.join(kb_path, 'tpl', f'{template}.tpl.py')).read()
    text = text.replace('{{', '\x00').replace('}}', '\x01')
    text = text.replace('{Code}', code).replace('{label}', label)
    return text.replace('\x00', '{').replace('\x01', '}')


def run_header(template, code, label, kwargs):
    """Execute the lines from the calculator setup through the basename choice."""
    text = render(template, code, label)
    start = text.index(f'mol.calc = {code}(**kwargs)')
    end = text.index('\n', text.index('basename = ', start))
    mol = SimpleNamespace()
    mol.get_potential_energy = lambda: mol.calc.get_potential_energy()
    namespace = {'os': os, code: FakeCalculator, 'mol': mol, 'kwargs': dict(kwargs)}
    FakeCalculator.instances = []
    exec(text[start:end], namespace)
    if template == 'ase_sella_opt_well':
        # This startup is inside the worker's error handler so a failed
        # Gaussian call produces an error row instead of leaving no result.
        start = text.index("        if '" + code + "' == 'Gaussian' and kwargs.get('chk'):")
        end = text.index('        if len(mol.symbols)', start)
        exec(dedent(text[start:end]), namespace)
    return namespace


class TestSellaTemplatePaths(unittest.TestCase):
    def test_gaussian_jobs_in_subdirectories_keep_their_path(self):
        for template in TEMPLATES:
            with self.subTest(template=template):
                ns = run_header(template, 'Gaussian', 'conf/x_0001', {'method': 'am1'})
                self.assertEqual(ns['basename'], 'conf/x_0001')
                text = render(template, 'Gaussian', 'conf/x_0001')
                self.assertIn("open(f'{basename}_sella.log', 'a')", text)

    def test_orca_jobs_keep_the_base_name_for_the_scratch_copy(self):
        for template in TEMPLATES:
            with self.subTest(template=template):
                ns = run_header(template, 'ORCA', 'hir/x_hir_0_01', {})
                self.assertEqual(ns['basename'], 'x_hir_0_01')

    def test_guess_read_only_with_a_checkpoint(self):
        for template in TEMPLATES:
            with self.subTest(template=template):
                ns = run_header(template, 'Gaussian', 'conf/x_0001', {'method': 'am1'})
                self.assertNotIn('guess', ns['mol'].calc.kwargs)
                self.assertEqual(len(FakeCalculator.instances), 1)
                ns = run_header(template, 'Gaussian', 'x_well', {'method': 'am1', 'chk': 'x_well'})
                self.assertEqual(ns['mol'].calc.kwargs['guess'], 'Read')
                self.assertEqual(len(FakeCalculator.instances), 2)


if __name__ == '__main__':
    unittest.main()
