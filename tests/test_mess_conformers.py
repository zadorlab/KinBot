"""Render conformer thermochemistry without running MESS or QC jobs."""
from kinbot.species_routing import routing_key

import json
import copy
import logging
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from ase.build import molecule

from kinbot import constants
from kinbot.mess import MESS, apply_conformer_shifts, validate_mess_populations
from kinbot.mess_mirrors import _endpoint_references
from kinbot.parameters import Parameters
from kinbot import constants
from kinbot.pes import create_mess_input
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity
from kinbot.stereo_routing import StereoRoutingError
from kinbot.conformer_counting import representative_record
from kinbot.symmetry import calculate_symmetry
from test_molecular_symmetry import methoxy


def point(name, chemid, saddle=False):
    atoms = molecule('H2O')
    freq = [-1000., 2000., 3000.] if saddle else [1000., 2000., 3000.]
    species = StationaryPoint(name, 0, 1, atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    species.characterize()
    species.__dict__.update(
        name=name, chemid=chemid, smiles='O', natom=3, mult=1,
        atom=atoms.get_chemical_symbols(), geom=atoms.positions,
        energy=-76., zpe=.01, sigma_ext=2, nopt=1,
        freq=list(freq), reduced_freqs=freq, conformer_index=[0, -999, 2],
        conformer_geom=[atoms.positions, None, atoms.positions * 1.3],
        conformer_freq=[freq, None, [freq[0] * 1.1, 2100., 3100.]],
        conformer_zeroenergy=[-75.99, -1000., -75.98],
        reac_type=['test'])
    return species


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

    def test_distorted_members_keep_legacy_rotational_symmetry_in_direct_and_pes(self):
        for barrier in (False, True):
            species = self.ts if barrier else self.well
            species.conformer_geom[2] = species.conformer_geom[2].copy()
            species.conformer_geom[2][1] *= 1.6
            direct = self.render(barrier)
            self.assertEqual(values(direct, 'SymmetryFactor'), [2., 2.])
            self.writer.par['pes'] = 1
            deferred = self.render(barrier)
            self.assertEqual(values(deferred, 'SymmetryFactor'), [2., 2.])
            self.writer.par['pes'] = 0

    def test_methoxy_symmetry_is_three_in_final_direct_and_pes_blocks(self):
        p = methoxy()
        calculate_symmetry(p)
        p.conformer_geom = [p.geom]
        p.conformer_freq = [p.freq]
        p.conformer_zeroenergy = [p.energy + p.zpe]
        self.writer.well_names[routing_key(p)] = 'methoxy'
        self.writer.fragment_names = {p.chemid: 'methoxy', self.well.chemid: 'water'}
        key = '_'.join(sorted(map(str, [routing_key(p), routing_key(self.well)])))
        self.writer.bimolec_names = {key: 'products'}
        for mc in (0, 1):
            for pes in (0, 1):
                for indices in ([7], []):
                    with self.subTest(mc=mc, pes=pes, indices=indices):
                        self.writer.par.update(multi_conf_tst=mc, pes=pes)
                        p.conformer_index = indices
                        self.assertEqual(values(self.writer.write_well(p, 0., 1., 0),
                                                'SymmetryFactor'), [3.])
                        fragment = self.writer.write_bimol([p, self.well], 0., 1., 1., 0, 0)
                        self.assertEqual(values(fragment, 'SymmetryFactor')[0], 3.)
        self.assertEqual(p.sigma_ext, 3.)

    def test_fragment_union_preserves_member_offsets_and_total_ground_energy(self):
        other = point('second', 234)
        self.writer.fragment_names = {123: 'f1', 234: 'f2'}
        self.writer.bimolec_names = {'123_234': 'p1'}
        self.well.conformer_zeroenergy[0] -= .02
        direct = self.writer.write_bimol([self.well, other], 0., 1., 1., 0, 0)
        self.assertEqual(direct.count('End ! Union'), 2)
        self.assertEqual(values(direct, 'ZeroEnergy'), [0., 18.83, 0., 6.28])
        self.writer.par['pes'] = 1
        deferred = self.writer.write_bimol([self.well, other], 0., 1., 1., 0, 0)
        reference = ((self.well.energy + self.well.zpe) * constants.AUtoKCAL)
        assembled = apply_conformer_shifts(deferred.format(
            name='p1', fr_name_123='f1', fr_name_234='f2', ground_energy=round(reference, 2)))
        self.assertEqual(values(assembled, 'ZeroEnergy'), values(direct, 'ZeroEnergy'))
        self.assertAlmostEqual(values(assembled, 'GroundEnergy')[0],
                               values(direct, 'GroundEnergy')[0], delta=.01)
        # The retained selected-parent sum must survive independently of the
        # fragment minimum, even if no attached saddle supplies an Eckart depth.
        self.assertEqual(_endpoint_references([direct])['p1'], round(reference, 2))
        self.assertEqual(_endpoint_references([assembled])['p1'], round(reference, 2))

    def test_mc_well_keeps_parent_reference_when_no_tunneling_path_exists(self):
        self.well.conformer_zeroenergy[0] -= .02
        direct = self.render()
        self.assertLess(min(values(direct, 'ZeroEnergy')), 0.)
        self.assertEqual(_endpoint_references([direct])['w1'], 0.)
        self.writer.par['pes'] = 1
        deferred = self.render()
        assembled = apply_conformer_shifts(deferred.format(name='w1', zeroenergy=5.))
        self.assertLess(min(values(assembled, 'ZeroEnergy')), 5.)
        self.assertEqual(_endpoint_references([assembled])['w1'], 5.)

    def test_no_surviving_conformer_uses_valid_parent(self):
        for barrier in (False, True):
            with self.subTest(barrier=barrier):
                species = self.ts if barrier else self.well
                species.conformer_index = [-999] * 3
                output = self.render(barrier)
                self.assertNotIn('Union', output)
                self.assertIn('2000.0', output)

    def test_specified_chiral_mc_fallback_uses_member_counting_in_all_blocks(self):
        p = StationaryPoint('fixed', 0, 1, smiles='C[C@H](O)CC')
        p.characterize()
        p.energy, p.zpe, p.sigma_ext, p.nopt = -76., .01, 1, 2
        p.freq = p.reduced_freqs = [100.] * (3 * p.natom - 6)
        p.conformer_index = []
        self.writer.well_names[p.name] = 'fixed'
        self.writer.well_names[routing_key(p)] = 'fixed'
        self.assertEqual(values(self.writer.write_well(p, 0., 1., 0), 'SymmetryFactor'), [1.])
        ts = copy.copy(p)
        ts.name, ts.wellorts = 'fixed_ts', 1
        ts.freq = ts.reduced_freqs = [-1000.] + p.freq[1:]
        ts.optical_reference = canonical_identity(p)
        self.writer.ts_names[ts.name] = 'ts2'
        reaction = SimpleNamespace(ts=ts, products=[self.well], instance_name=ts.name)
        p.reac_type = self.well.reac_type
        self.writer.species = p  # actual configured reactant of this channel
        output = self.writer.write_barrier(reaction, 0, 30., 20., 0., 1., 1., 0)[0]
        self.assertEqual(values(output, 'SymmetryFactor'), [1.])
        self.writer.fragment_names = {routing_key(p): 'f1', self.well.chemid: 'f2'}
        self.writer.bimolec_names = {'_'.join(sorted(map(str, [routing_key(p), routing_key(self.well)]))): 'p1'}
        self.well.conformer_index = []
        output = self.writer.write_bimol([p, self.well], 0., 1., 1., 0, 0)
        self.assertEqual(values(output, 'SymmetryFactor'), [1., 2.])
        self.assertEqual(p.rrho_representative_counting['remaining_optical_weight'], 1.)
        self.writer.par['multi_conf_tst'] = 0
        # A specified enantiomer excludes its global mirror with MC off too.
        self.assertEqual(values(self.writer.write_well(p, 0., 1., 0), 'SymmetryFactor'), [1.])
        output = self.writer.write_barrier(reaction, 0, 30., 20., 0., 1., 1., 0)[0]
        self.assertEqual(values(output, 'SymmetryFactor'), [1.])
        self.writer.par['optical_population'] = 'racemic'
        self.assertEqual(values(self.writer.write_well(p, 0., 1., 0), 'SymmetryFactor'), [.5])

    def test_mc_variational_barrier_uses_legacy_representative_with_warning(self):
        self.well.reac_type[0] = 'barrierless_saddle'
        self.reaction.prod_opt = [SimpleNamespace(species=self.well), SimpleNamespace(species=self.well)]
        output = self.render(barrier=True)
        self.assertIn('additional conformers are not included', output)
        self.assertIn('Variational', output)
        self.assertNotIn('End ! Union', output)
        self.assertEqual(output.count('End ! RRHO'), 2)  # inner RRHO and outer PST
        self.assertEqual(values(output, 'ZeroEnergy')[-1], 30.)

    def test_nonmc_separately_named_mirrors_cannot_each_claim_the_same_racemate(self):
        self.writer.par.update(multi_conf_tst=0, optical_population='racemic')
        blocks = []
        for smiles in ('C[C@H](O)CC', 'C[C@@H](O)CC'):
            p = StationaryPoint('fixed', 0, 1, smiles=smiles)
            p.characterize()
            p.energy, p.zpe, p.sigma_ext, p.nopt = -76., .01, 1, 2
            p.freq = p.reduced_freqs = [100.] * (3*p.natom-6)
            self.writer.well_names[routing_key(p)] = str(routing_key(p))
            blocks.append(self.writer.write_well(p, 0., 1., 0))
        with self.assertRaisesRegex(ValueError, 'Overlapping racemic populations'):
            validate_mess_populations('\n'.join(blocks))

    def test_ts_mirror_population_cannot_be_routed_to_one_specified_chiral_product(self):
        def prepared(smiles):
            p = StationaryPoint('configured', 0, 1, smiles=smiles)
            p.characterize()
            p.energy, p.zpe = -20., .01
            p.freq = p.reduced_freqs = [100.] * (3 * p.natom - 6)
            p.reac_type, p.conformer_index = ['test'], []
            return p
        reactant, product = prepared('CCCCO'), prepared('C[C@H](O)CC')
        ts = copy.copy(product)
        ts.name, ts.wellorts = 'counting_fixture', 1
        ts.freq = ts.reduced_freqs = [-1000.] + product.freq[1:]
        ts.optical_reference = canonical_identity(reactant)
        self.writer.species = reactant
        reaction = SimpleNamespace(ts=ts, products=[product], instance_name=ts.name)
        for explicit in (False, True):
            for fragmented in (False, True):
                reaction.products = [product, self.well] if fragmented else [product]
                if explicit:
                    ts.conformer_index = [0, 1]
                    ts.conformer_geom = [ts.geom, ts.geom * [-1, 1, 1]]
                    ts.conformer_freq = [ts.freq, ts.freq]
                    ts.conformer_zeroenergy = [-19.99, -19.99]
                with self.assertRaisesRegex(StereoRoutingError, 'TS mirror orbit'):
                    self.writer._validate_mc_endpoints(reaction, [representative_record(ts)])
        # The global mirror of an unordered R+S endpoint is the same pair.
        mirror_product = copy.copy(product)
        mirror_product.geom = product.geom * [-1, 1, 1]
        reaction.products = [product, mirror_product]
        self.writer._validate_mc_endpoints(reaction, [representative_record(ts)])
        # A fixed-chiral reactant excludes its mirrored TS under specified scope.
        ts.optical_reference = canonical_identity(product)
        reaction.products = [product]
        self.writer._validate_mc_endpoints(reaction, [representative_record(ts)])
        # A self-mirror saddle does not introduce an omitted optical partner.
        self.ts.wellorts = 1
        self.ts.optical_reference = canonical_identity(reactant)
        reaction.ts = self.ts
        self.writer._validate_mc_endpoints(reaction, [representative_record(self.ts)])
        self.writer.par['optical_population'] = 'racemic'
        with self.assertRaisesRegex(StereoRoutingError, 'independent racemates'):
            self.writer._mc_product_identities([product, mirror_product])
        self.writer._mc_product_identities([product, self.well])
        self.writer.par['multi_conf_tst'] = 0
        with self.assertRaisesRegex(StereoRoutingError, 'independent racemates'):
            self.writer.write_bimol([product, mirror_product], 0., 1., 1., 0, bless=False)

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

    def test_full_pes_keeps_declared_parent_eckart_depths_below_ensemble_ground_shifts(self):
        self.well.conformer_zeroenergy[0] -= .02
        self.well.source_job = 'selected_well'
        self.ts.source_job = 'selected_ts'
        direct = self.render(barrier=True)
        self.assertEqual(values(direct, 'WellDepth'), [30., 20., 36.28, 26.28])
        reference = self.ts.mess_tunneling_reference
        self.assertEqual(reference['left_parent_observation']['source_job'], 'selected_well')
        self.assertFalse(reference['ensemble_minima_redefine_depths'])
        self.writer.par['pes'] = 1
        barrier, well = self.render(barrier=True), self.render()
        Path('123').mkdir()
        Path('123/saddle_0000.mess').write_text(barrier)
        for name in ('123', '789'):
            Path(f'123/{name}_0000.mess').write_text(well)
        Path('input.json').write_text(json.dumps({'barrier_threshold': 100., 'smiles': 'O',
            'multi_conf_tst': 1, 'conformer_search': 1, 'me': 0, 'uq': 0,
            'epsilon': 100., 'sigma': 3.}))
        par = Parameters('input.json', show_warnings=False).par
        with patch('kinbot.pes.logger', logging.getLogger('KinBot'), create=True):
            create_mess_input(par, ['123', '789'], [], [['123', 'saddle', ['789'], 30.]],
                [], [], {'123': 0., '789': 10.}, {}, {'123': '123', '789': '123'}, 18., False)
        output = Path('me/mess_0000.inp').read_text()
        self.assertIn('selected-parent endpoint approximation', output)
        self.assertIn(-12.55, values(output, 'ZeroEnergy'))
        self.assertEqual(values(output, 'WellDepth'), values(direct, 'WellDepth'))
        self.assertIsNone(self.ts.mess_tunneling_reference['input_depths_kcal_mol'])

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
            if line.strip() and line.split()[0].split('[')[0] == keyword]


if __name__ == '__main__':
    unittest.main()
