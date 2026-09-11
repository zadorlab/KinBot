"""Final direct/PES inputs use actual serialized MC endpoint admissibility bounds."""
import json
import logging
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from kinbot import constants
from kinbot.mess import MESS, finalize_mc_mess
from kinbot.parameters import Parameters
from kinbot.pes import create_mess_input
from tests.test_mess_conformers import point, values


class TestFinalMCInput(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(directory.name)
        Path('me').mkdir()
        Path('123').mkdir()
        Path('input.json').write_text(json.dumps({'barrier_threshold': 100.,
            'smiles': 'O', 'multi_conf_tst': 1, 'conformer_search': 1,
            'high_level': 0, 'me': 0, 'uq': 0, 'epsilon': 100., 'sigma': 3.}))
        self.par = Parameters('input.json', show_warnings=False).par

    def render(self, endpoint_shift, ts_energy, corrected, pes, survivors=1):
        well, product, ts = point('water', 123), point('other', 789), point('saddle', 456, True)
        for p, parent in ((well, 0.), (product, 10.), (ts, ts_energy)):
            p.energy += parent / constants.AUtoKCAL
            p.conformer_zeroenergy = [e + parent / constants.AUtoKCAL for e in p.conformer_zeroenergy]
        for p in (well, product):
            p.conformer_index = [-999, -999, 2] if survivors == 1 else [0, -999, 2]
            p.conformer_zeroenergy[2] = p.energy + p.zpe + endpoint_shift / constants.AUtoKCAL
            if survivors == 2:
                p.conformer_zeroenergy[0] = p.conformer_zeroenergy[2] + .01
            if survivors == 0:
                p.conformer_index = [-999] * 3
        well.mass = 18.
        reaction = SimpleNamespace(ts=ts, products=[product], instance_name='saddle',
                                   do_vdW=False, prod_opt=[SimpleNamespace(species=product)], mp2=0)
        well.reac_obj, well.reac_ts_done = [reaction], [-1]
        self.par.update(pes=int(pes), correct_submerged=int(corrected))
        writer = MESS(self.par, well)
        if not pes:
            writer.write_input(None)
        else:
            writer.create_short_names()
            Path('123/123_0000.mess').write_text(writer.write_well(well, 0., 1., 0))
            Path('123/789_0000.mess').write_text(writer.write_well(product, 0., 1., 0))
            Path('123/saddle_0000.mess').write_text(writer.write_barrier(
                reaction, 0, ts_energy, ts_energy - 10., 0., 1., 1., 0)[0])
            with patch('kinbot.pes.logger', logging.getLogger('KinBot'), create=True):
                create_mess_input(self.par, ['123', '789'], [],
                    [['123', 'saddle', ['789'], ts_energy]], [], [],
                    {'123': 0., '789': 10.}, {}, {'123': '123', '789': '123'}, 18., False)
        output = Path('me/mess_0000.inp').read_text()
        # Both routes write one barrier and an optional MC Union.
        return output[output.index('  Barrier'):]

    def test_final_input_corrects_actual_well_union_bounds_in_both_routes(self):
        for pes in (False, True):
            for survivors in (1, 2):
                for corrected in (False, True):
                    with self.subTest(pes=pes, survivors=survivors, corrected=corrected):
                        output = self.render(6.28, 14., corrected, pes, survivors)
                        self.assertEqual(values(output, 'ZeroEnergy'),
                                         [16.28 if corrected else 14., 20.28])
                        # The corrected member at the well floor has no tunnel,
                        # even though its parent-relative depths would be positive.
                        self.assertEqual(output.count('Tunneling   Eckart'), 1)
                        self.assertEqual(values(output, 'WellDepth'), [20.28, 10.28])

    def test_lower_endpoint_member_does_not_inherit_parent_admissibility_floor(self):
        for pes in (False, True):
            output = self.render(-6.28, 4., True, pes)
            self.assertEqual(values(output, 'ZeroEnergy'), [4., 10.28])
            # Retain the selected-parent Eckart approximation away from the bound.
            self.assertEqual(values(output, 'WellDepth'), [10.28, .28])

    def test_parent_fallback_remains_the_endpoint_ground(self):
        for pes in (False, True):
            output = self.render(6.28, 4., True, pes, survivors=0)
            self.assertEqual(values(output, 'ZeroEnergy'), [10., 10.28])
            self.assertEqual(values(output, 'WellDepth'), [10.28, .28])

    def test_final_rounding_and_bimolecular_ground_semantics(self):
        # Fragment RRHO zeros do not replace the explicit product GroundEnergy.
        text = ('Well w\nZeroEnergy[kcal/mol] 10.01\n'
                'Bimolecular p\nZeroEnergy[kcal/mol] 500\nGroundEnergy[kcal/mol] 5\n'
                'Barrier b w p\nZeroEnergy[kcal/mol] 30 ! -60\n'
                'Tunneling Eckart\nImaginaryFrequency[1/cm] 1000\n'
                'CutoffEnergy[kcal/mol] 20\nWellDepth[kcal/mol] 20\n'
                'WellDepth[kcal/mol] 25\nEnd\nEnd ! RRHO\n')
        output = finalize_mc_mess(text, correct_submerged=True)
        self.assertEqual(values(output, 'ZeroEnergy'), [10.01, 500., 10.01])
        self.assertNotIn('Tunneling Eckart', output)
        self.assertEqual(finalize_mc_mess(output, True), output)


if __name__ == '__main__':
    unittest.main()
