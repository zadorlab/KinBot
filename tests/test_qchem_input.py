"""Native Q-Chem inputs must carry constraints in their own $opt block (#67)."""

import os
import pathlib
import tempfile
import unittest

from ase import Atoms

from kinbot.ase_modules.calculators.qchem import QChem

TPL = pathlib.Path(__file__).resolve().parent.parent / 'kinbot' / 'tpl'


class TestQChemAddsec(unittest.TestCase):
    def render(self, **extra):
        atoms = Atoms('H2', positions=[[0., 0., 0.], [0., 0., 0.74]])
        kwargs = dict(label='probe', jobtype='opt', method='b3lyp', basis='sto-3g',
                      unrestricted='true', nt=2)
        kwargs.update(extra)
        with tempfile.TemporaryDirectory() as directory:
            cwd = os.getcwd()
            os.chdir(directory)
            try:
                QChem(**kwargs).write_input(atoms)
                return pathlib.Path('probe.inp').read_text()
            finally:
                os.chdir(cwd)

    def test_constraints_are_a_separate_opt_block_not_a_rem_variable(self):
        # the constraint text reac_family.py builds for a Q-Chem TS search
        addsec = '$opt\nCONSTRAINT\nSTRE 1 2 0.74\nENDCONSTRAINT\n$end\n'
        text = self.render(addsec=addsec)
        rem = text.split('$rem')[1].split('$end')[0]
        self.assertNotIn('ADDSEC', rem.upper())
        self.assertNotIn('addsec', text.lower().replace('$opt', ''))
        self.assertIn('$opt\nCONSTRAINT\nSTRE 1 2 0.74\nENDCONSTRAINT\n$end', text)
        # the block follows $molecule, where Q-Chem expects it
        self.assertGreater(text.index('$opt'), text.index('$molecule'))

    def test_template_kwargs_are_accepted(self):
        # what qc.py hands the calculator for a Q-Chem frequency recovery job
        text = self.render(jobtype='freq', vibman_print='4', xc_grid='3')
        jobtype = next(line.split()[-1] for line in text.splitlines() if 'JOBTYPE' in line)
        self.assertEqual(jobtype, 'FREQ')
        self.assertIn('VIBMAN_PRINT', text.upper())
        self.assertIn('SYM_IGNORE', text.upper())

    def test_no_template_uses_the_stock_ase_qchem_calculator(self):
        """The stock class writes every keyword into $rem, including addsec."""
        offenders = [p.name for p in TPL.glob('*.py')
                     if 'from ase.calculators.qchem import QChem' in p.read_text()]
        self.assertEqual(offenders, [])
        users = [p.name for p in TPL.glob('ase_qchem_*.tpl.py')]
        self.assertTrue(users)
        for name in users:
            self.assertIn('from kinbot.ase_modules.calculators.qchem import QChem',
                          (TPL / name).read_text(), name)


if __name__ == '__main__':
    unittest.main()
