"""Final ordinary RRHO/HIR output keeps corrected modes separate from raw data."""
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from ase.build import molecule

from kinbot.calculation import load_calculation_record
from kinbot.mess import MESS
from kinbot.optimize import Optimize
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.stationary_pt import StationaryPoint


class TestOrdinaryThermalFrequencies(unittest.TestCase):
    def test_final_writer_corrects_soft_modes_and_preserves_raw_record(self):
        previous = Path.cwd()
        with TemporaryDirectory() as directory:
            try:
                os.chdir(directory)
                Path('input.json').write_text(json.dumps({'barrier_threshold': 100.}))
                par = Parameters('input.json', show_warnings=False).par
                par.update(high_level=0, rotor_scan=0, conformer_search=0,
                           multi_conf_tst=0, L3_calc=0, pes=0)
                qc = QuantumChemistry(par)
                atoms = molecule('CH3OH')
                for saddle in (0, 1):
                    for projected in (False, True):
                        with self.subTest(saddle=saddle, projected=projected):
                            p = StationaryPoint('thermal', 0, 1, wellorts=saddle,
                                atom=atoms.get_chemical_symbols(), geom=atoms.positions)
                            p.characterize()
                            raw = ([-1000., -20.] if saddle else [-20., 100.]) + [200.] * 10
                            job = f'selected_{saddle}_{projected}'
                            qc.db.write(atoms, name=job, data={'energy': -100., 'zpe': .01,
                                'frequencies': raw, 'hess': np.eye(3 * len(atoms)), 'status': 'normal'})
                            load_calculation_record(p, qc, job)
                            opt = Optimize(p, par, qc)
                            opt.selected_job, opt.shigh, opt.shir = job, 1, 1
                            reduced = raw[:-1] if projected else raw
                            par['rotor_scan'] = int(projected)
                            with patch('kinbot.optimize.frequencies.get_frequencies', return_value=(raw, reduced)), \
                                 patch.object(opt, '_ensure_selected_hessian', return_value=True):
                                opt.do_optimization()
                            self.assertEqual(p.freq, raw)
                            self.assertEqual(list(qc.db.get(name=job).data.frequencies), raw)
                            self.assertIn(20., p.reduced_freqs)
                            self.assertNotIn(-20., p.reduced_freqs)
                            writer = MESS(par, p)
                            writer.well_names = {p.chemid: 'w1'}
                            writer.ts_names = {'saddle': 'ts1'}
                            # Projection itself is tested elsewhere; no rotor potential in this fixture.
                            with patch.object(writer, 'make_rotors', return_value=''):
                                if saddle:
                                    p.reac_type = ['test']
                                    reaction = SimpleNamespace(ts=p, products=[p], instance_name='saddle')
                                    output = writer.write_barrier(reaction, 0, 30., 20., 0., 1., 1., 0)[0]
                                    self.assertEqual(p.reduced_freqs[0], -1000.)
                                else:
                                    output = writer.write_well(p, 0., 1., 0)
                            self.assertNotIn('-20.0', output)
                            self.assertIn('20.0', output)
            finally:
                os.chdir(previous)


if __name__ == '__main__':
    unittest.main()
