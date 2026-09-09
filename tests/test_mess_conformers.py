"""Render conformer thermochemistry without running MESS or QC jobs."""

import json
import logging
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from ase.build import molecule

from kinbot import constants
from kinbot.mess import MESS, apply_conformer_shifts
from kinbot.parameters import Parameters
from kinbot.pes import create_mess_input


def point(name, chemid, saddle=False):
    atoms = molecule('H2O')
    freq = [-1000., 2000., 3000.] if saddle else [1000., 2000., 3000.]
    return SimpleNamespace(
        name=name, chemid=chemid, smiles='O', natom=3, mult=1,
        atom=atoms.get_chemical_symbols(), geom=atoms.positions,
        energy=-76., zpe=.01, sigma_ext=2, nopt=1,
        reduced_freqs=freq, conformer_index=[0, -999, 2],
        conformer_geom=[atoms.positions, None, atoms.positions * 1.1],
        conformer_freq=[freq, None, [freq[0] * 1.1, 2100., 3100.]],
        conformer_zeroenergy=[-75.99, -1000., -75.98],
        reac_type=['test'])


class TestConformerSerialization(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        self.well = point('water', 123)
        self.ts = point('saddle', 456, saddle=True)
        self.writer = MESS({
            'pes': 0, 'multi_conf_tst': 1, 'rotor_scan': 0,
            'freq_uq_ref': 100., 'freq_uq_max_exp': 2.,
        }, self.well)
        self.writer.well_names = {123: 'w1'}
        self.writer.ts_names = {'saddle': 'ts1'}
        self.reaction = SimpleNamespace(
            ts=self.ts, products=[self.well], instance_name='saddle')

    def render(self, barrier=False):
        if barrier:
            return self.writer.write_barrier(
                self.reaction, 0, 30., 20., 0., 1., 1., 0)[0]
        return self.writer.write_well(self.well, 0., 1., 0)

    def test_failed_conformer_is_not_rendered_or_used_as_energy_reference(self):
        for barrier in (False, True):
            with self.subTest(barrier=barrier):
                output = self.render(barrier)
                self.assertEqual(output.count('End ! RRHO'), 2)
                self.assertIn('number of species in union is 2', output)
                energies = values(output, 'ZeroEnergy')
                reference = 30. if barrier else 0.
                self.assertAlmostEqual(energies[0], reference)
                self.assertAlmostEqual(energies[1], reference + 6.28)

    def test_single_surviving_conformer_uses_its_own_properties(self):
        for barrier in (False, True):
            with self.subTest(barrier=barrier):
                species = self.ts if barrier else self.well
                species.conformer_index = [-999, -999, 2]
                output = self.render(barrier)
                self.assertIn('2100.0', output)
                self.assertIn('3100.0', output)

    def test_no_surviving_conformer_uses_valid_parent(self):
        for barrier in (False, True):
            with self.subTest(barrier=barrier):
                species = self.ts if barrier else self.well
                species.conformer_index = [-999] * 3
                output = self.render(barrier)
                self.assertNotIn('Union', output)
                self.assertIn('2000.0', output)

    def test_direct_output_applies_each_ground_energy_and_tunneling_parameters(self):
        output = self.render()
        self.assertEqual(values(output, 'ZeroEnergy'), [0., 6.28])
        output = self.writer.write_barrier(
            self.reaction, 0, 30., 20., 0., 1., 1.1, 0)[0]
        self.assertEqual(values(output, 'ZeroEnergy'), [30., 36.28])
        self.assertEqual(values(output, 'ImaginaryFrequency'), [1100., 1210.])
        self.assertEqual(values(output, 'CutoffEnergy'), [20., 26.28])
        self.assertEqual(values(output, 'WellDepth'), [30., 20., 36.28, 26.28])

    def test_conformer_offsets_use_parent_reference_after_l2_order_reversal(self):
        self.well.conformer_zeroenergy[0] -= .02
        self.assertEqual(values(self.render(), 'ZeroEnergy'), [-12.55, 6.28])

    def test_direct_and_pes_routes_agree_without_double_applying_offsets(self):
        for barrier in (False, True):
            with self.subTest(barrier=barrier):
                self.writer.par['pes'] = 0
                direct = self.render(barrier)
                self.writer.par['pes'] = 1
                intermediate = self.render(barrier)
                resolved = intermediate.format(
                    name='reference', zeroenergy=30. if barrier else 0.,
                    cutoff=20., welldepth1=30., welldepth2=20.)
                final = apply_conformer_shifts(resolved)
                for keyword in ('ZeroEnergy', 'ImaginaryFrequency',
                                'CutoffEnergy', 'WellDepth'):
                    self.assertEqual(values(final, keyword), values(direct, keyword))
                self.assertEqual(apply_conformer_shifts(final), final)

    def test_submerged_conformer_does_not_reuse_parent_tunneling(self):
        self.ts.conformer_zeroenergy[0] -= .1
        output = self.render(barrier=True)
        self.assertEqual(output.count('Tunneling   Eckart'), 1)
        self.assertIn('! barrier is submerged', output)

    def test_complete_pes_output_removes_submerged_member_tunneling_after_shifts(self):
        self.ts.conformer_zeroenergy[0] -= .1
        direct = self.render(barrier=True)
        self.writer.par['pes'] = 1
        barrier = self.render(barrier=True)
        well = self.render()
        Path('123').mkdir()
        Path('123/saddle_0000.mess').write_text(barrier)
        for name in ('123', '789'):
            Path(f'123/{name}_0000.mess').write_text(well)
        Path('input.json').write_text(json.dumps({
            'barrier_threshold': 100., 'smiles': 'O', 'multi_conf_tst': 1,
            'conformer_search': 1,
            'me': 0, 'uq': 0, 'epsilon': 100., 'sigma': 3.,
        }))
        par = Parameters('input.json', show_warnings=False).par
        # The PES command initializes its logger in main().
        with patch('kinbot.pes.logger', logging.getLogger('KinBot'), create=True):
            create_mess_input(par, ['123', '789'], [],
                [['123', 'saddle', ['789'], 30.]], [], [],
                {'123': 0., '789': 10.}, {}, {'123': '123', '789': '123'}, 18., False)
        output = Path('me/mess_0000.inp').read_text().split('# BARRIERS\n', 1)[1]
        self.assertEqual(output.count('End ! RRHO'), 2)
        self.assertEqual(output.count('Tunneling   Eckart'), 1)
        for keyword in ('ZeroEnergy', 'ImaginaryFrequency', 'CutoffEnergy', 'WellDepth'):
            self.assertEqual(values(output, keyword), values(direct, keyword))

    def test_complete_pes_corrects_resolved_member_grounds_only_when_requested(self):
        self.ts.conformer_zeroenergy[0] -= .1
        self.writer.par['pes'] = 1
        Path('123').mkdir()
        Path('123/saddle_0000.mess').write_text(self.render(barrier=True))
        well = self.render()
        for name in ('123', '789'):
            Path(f'123/{name}_0000.mess').write_text(well)
        Path('input.json').write_text(json.dumps({'barrier_threshold': 100.,
            'smiles': 'O', 'multi_conf_tst': 1, 'conformer_search': 1,
            'me': 0, 'uq': 0, 'epsilon': 100., 'sigma': 3.}))
        par = Parameters('input.json', show_warnings=False).par
        for corrected in (0, 1):
            par['correct_submerged'] = corrected
            with patch('kinbot.pes.logger', logging.getLogger('KinBot'), create=True):
                create_mess_input(par, ['123', '789'], [],
                    [['123', 'saddle', ['789'], 30.]], [], [],
                    {'123': 0., '789': 10.}, {}, {'123': '123', '789': '123'}, 18., False)
            output = Path('me/mess_0000.inp').read_text().split('# BARRIERS\n', 1)[1]
            self.assertEqual(values(output, 'ZeroEnergy'),
                             [10. if corrected else -32.75, 36.28])
            self.assertEqual(output.count('Tunneling   Eckart'), 1)
            self.assertEqual(values(output, 'WellDepth'), [36.28, 26.28])

    def test_tunneling_parameters_must_remain_positive_after_serialization(self):
        for mc in (0, 1):
            self.writer.par['multi_conf_tst'] = mc
            for depth in (-.001, 0., .001, .006):
                with self.subTest(mc=mc, depth=depth):
                    # The first conformer has no shift; the second is higher.
                    output = self.writer.write_barrier(
                        self.reaction, 0, 10. + depth, depth, 0., 1., 1., 0)[0]
                    expected = (1 if mc else 0) + int(round(depth, 2) > 0.)
                    self.assertEqual(output.count('Tunneling   Eckart'), expected)
                    self.assertTrue(all(value > 0. for keyword in ('CutoffEnergy', 'WellDepth')
                                        for value in values(output, keyword)))
        # Also exercise an offset that nearly submerges just one MC member.
        self.writer.par['multi_conf_tst'] = 1
        self.ts.conformer_zeroenergy[0] -= 19.999 / constants.AUtoKCAL
        direct = self.render(barrier=True)
        self.assertEqual(direct.count('Tunneling   Eckart'), 1)
        self.writer.par['pes'] = 1
        deferred = apply_conformer_shifts(self.render(barrier=True).format(
            name='barrier', zeroenergy=30., cutoff=20., welldepth1=30., welldepth2=20.))
        # PES removes this zero-cutoff block in its final serialization pass.
        self.assertEqual(values(deferred, 'CutoffEnergy')[0], 0.)

    def test_accepted_small_imaginary_modes_are_corrected_without_mutating_raw_members(self):
        for barrier in (False, True):
            species = self.ts if barrier else self.well
            raw = [-1000., -20., 3000.] if barrier else [-20., 2000., 3000.]
            species.conformer_freq[0] = raw.copy()
            output = self.render(barrier)
            self.assertNotIn('-20.0', output)
            self.assertIn('20.0', output)
            self.assertEqual(species.conformer_freq[0], raw)

    def test_legacy_pes_comments_and_unshifted_blocks_remain_supported(self):
        text = ('ZeroEnergy[kcal/mol] 0\nCutoffEnergy[kcal/mol] 4\nEnd ! RRHO\n'
                'ZeroEnergy[kcal/mol] 10 ! 2\n'
                'ImaginaryFrequency[1/cm] 1000 ! -1200\n'
                'CutoffEnergy[kcal/mol] 5\nEnd ! RRHO\n'
                'ZeroEnergy[kcal/mol] 30\nCutoffEnergy[kcal/mol] 7\n')
        result = apply_conformer_shifts(text)
        self.assertEqual(values(result, 'ZeroEnergy'), [0., 12., 30.])
        self.assertEqual(values(result, 'CutoffEnergy'), [4., 7., 7.])
        self.assertEqual(values(result, 'ImaginaryFrequency'), [1200.])


def values(output, keyword):
    return [float(line.split()[1]) for line in output.splitlines()
            if line.strip().startswith(keyword + '[')]


if __name__ == '__main__':
    unittest.main()
