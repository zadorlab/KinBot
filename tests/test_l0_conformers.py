"""Keep preliminary well-conformer results separate from L1 artifacts."""

import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from ase.build import molecule
from ase.io import read
import numpy as np

from kinbot import constants
from kinbot.calculation import selected_calculation_job
from kinbot.conformers import Conformers
from kinbot.optimize import Optimize
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint


class TestL0ConformerResults(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        Path('conf').mkdir()
        Path('input.json').write_text(json.dumps({'barrier_threshold': 50.}))
        self.par = Parameters('input.json', show_warnings=False).par
        self.atoms = molecule('CH3CH2OH')
        self.species = StationaryPoint('ethanol', 0, 1,
            atom=self.atoms.get_chemical_symbols(), geom=self.atoms.positions)
        self.species.characterize()
        self.species.name = str(self.species.chemid)
        self.qc = QuantumChemistry(self.par)
        # Read real ASE records, with scheduling disabled.
        self.qc.check_qc = Mock(return_value='normal')
        self.qc.qc_conf = Mock(return_value=0)
        self.search = Conformers(self.species, self.par, self.qc, semi_emp=1)
        self.low_job = f'conf/{self.species.chemid}_low'

    def record(self, job, energy, zpe=0., atoms=None, freq=None):
        if freq is None:
            freq = [100.] * (3 * self.species.natom - 6)
        self.qc.db.write(self.atoms if atoms is None else atoms, name=job,
            data={'energy': energy / constants.EVtoHARTREE, 'zpe': zpe,
                  'frequencies': freq, 'status': 'normal'})









    def test_l1_selection_loads_all_properties_without_an_irc_service(self):
        self.qc.qc = 'fc'
        self.par.update(conformer_search=1, high_level=0, rotor_scan=0,
                        multi_conf_tst=0, L3_calc=0)
        self.species.wellorts = 1
        self.species.name = 'ethanol_ts'
        parent_freq = [-1000.] + [100.] * (3 * self.species.natom - 7)
        selected_freq = [-900.] + [150.] * (3 * self.species.natom - 7)
        self.species.energy, self.species.zpe = -100., .1
        self.species.freq = parent_freq
        self.record(self.species.name, -100., .1, freq=parent_freq)
        optimization = Optimize(self.species, self.par, self.qc)
        search = Conformers(self.species, self.par, self.qc)
        self.species.confs = search
        self.record(search.get_job_name(0), -100., .1, freq=parent_freq)
        selected = self.atoms.copy()
        selected.positions += .01
        self.record(search.get_job_name(1), -100.01, .1001,
                    atoms=selected, freq=selected_freq)
        search.conf, search.conf_status = 2, [0, 0]
        optimization.scycconf = optimization.ssemi_empconf = 1
        optimization.sconf = 0
        optimization.do_optimization()
        self.assertEqual(self.species.source_job, search.get_job_name(1))
        np.testing.assert_allclose(self.species.geom, selected.positions)
        self.assertAlmostEqual(self.species.energy, -100.01)
        self.assertEqual(self.species.zpe, .1001)
        self.assertEqual(self.species.freq, selected_freq)
        self.assertEqual(self.species.reduced_freqs, selected_freq)
        selection = list(self.qc.db.select(name='conf/ethanol_ts_low'))[-1]
        self.assertEqual(selection.data.copied_from_job, search.get_job_name(1))



if __name__ == '__main__':
    unittest.main()
