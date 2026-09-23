"""Direct/PES complex selection and saddle-free dissociation output."""
import copy
from dataclasses import replace
import itertools
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from ase.build import molecule
from kinbot import constants, pes, symmetry
from kinbot.mess import MESS, apply_conformer_shifts
from kinbot.conformer_counting import representative_record
from test_mess_conformers import values
from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint
from kinbot.species_routing import routing_name
from kinbot.product_complex import reassess_product_complex


def species(formula, mult=1):
    atoms = molecule(formula)
    p = StationaryPoint(formula, 0, mult, atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    p.characterize()
    p.name = routing_name(p)
    symmetry.calculate_symmetry(p)
    p.energy, p.zpe = -40., .01
    p.freq = p.reduced_freqs = [1000.] * max(0, 3*p.natom - (5 if p.natom == 2 else 6))
    return p


class TestComplexMESS(unittest.TestCase):
    def setUp(self):
        tmp = TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(tmp.name)
        Path('me').mkdir()
        Path('input.json').write_text(json.dumps(dict(barrier_threshold=100.,
            high_level=0, rotor_scan=0, multi_conf_tst=0, conformer_search=0,
            me=0, uq=0, epsilon=100., sigma=3., queuing='local')))
        self.par = Parameters('input.json', show_warnings=False).par

    def reactions(self, homolysis=False):
        parent = species('CH4')
        products = [species('CH3', 2), species('H', 2)]
        products[0].energy, products[1].energy = -39., -1.
        reactions = []
        for name, barrier, depth in [('high', 35., -1.), ('low', 25., -1.01)]:
            ts = copy.deepcopy(parent)
            ts.name = name
            ts.energy += barrier/constants.AUtoKCAL
            ts.wellorts = 1
            ts.freq = ts.reduced_freqs = [-1000.] + [1000.] * 8
            complex_species = copy.deepcopy(parent)
            complex_species.name = name + '_IRC_F_prod'
            complex_species.energy += depth/constants.AUtoKCAL
            reactions.append(SimpleNamespace(instance_name=name, ts=ts,
                products=products, prod_opt=[SimpleNamespace(species=p) for p in products],
                do_vdW=not homolysis, mp2=0, irc_prod=complex_species,
                irc_prod_opt=SimpleNamespace(species=complex_species)))
        if homolysis:
            reactions = reactions[:1]
            reactions[0].instance_name = 'test_hom_sci'
            reactions[0].ts = copy.deepcopy(parent)
        parent.reac_obj = reactions
        parent.reac_ts_done = [-1] * len(reactions)
        parent.reac_type = ['hom_sci' if homolysis else 'test'] * len(reactions)
        return parent, reactions

    def test_direct_selects_lowest_ordinary_barrier_without_complex(self):
        parent, reactions = self.reactions()
        for order in (reactions, reactions[::-1]):
            parent.reac_obj = order
            writer = MESS(self.par, parent)
            writer.write_input(None)
            output = Path('me/mess_0000.inp').read_text()
            self.assertNotIn(' high\n', output)
            self.assertNotIn('{wellname}', output)
            headers = [line.split('!')[0].split() for line in output.splitlines()
                       if line.strip().startswith('Barrier ')]
            low = writer.ts_names['low']
            self.assertEqual(len(headers), 1)
            self.assertEqual(headers[0][:3], ['Barrier', low, 'w_1'])
            self.assertTrue(headers[0][3].startswith('b_'))
            self.assertEqual(sum(l.strip().startswith('Well ') for l in output.splitlines()), 1)
            self.assertNotIn('Core PhaseSpaceTheory', output)
            self.assertEqual(values(output, 'WellDepth'), [25., 26.01])
            self.assertTrue(all(r.do_vdW for r in reactions))
            self.assertNotIn('high_IRC_F_prod', output)

    def test_final_unbound_complex_is_omitted_without_losing_inner_saddle(self):
        parent, reactions = self.reactions()
        separated = sum(p.energy + p.zpe for p in reactions[0].products)
        saved = []
        for reaction in reactions:
            p = reaction.irc_prod_opt.species
            p.energy = separated - p.zpe + .02/constants.AUtoKCAL
            reaction.vdW_depth = 3.04  # initial-well comparison, now obsolete
            saved.append(p.geom.copy())
            reassess_product_complex(reaction, self.par)
        writer = MESS(self.par, parent)
        writer.write_input(None)
        text = Path('me/mess_0000.inp').read_text()
        self.assertEqual(sum(l.strip().startswith('Well ') for l in text.splitlines()), 1)
        self.assertEqual(sum(l.strip().startswith('Barrier ') for l in text.splitlines()), 1)
        self.assertNotIn('Core PhaseSpaceTheory', text)
        self.assertIn(' low\n', text)
        for reaction, geom in zip(reactions, saved):
            self.assertFalse(reaction.do_vdW)
            self.assertAlmostEqual(reaction.vdW_depth, -.02)
            np.testing.assert_array_equal(reaction.irc_prod_opt.species.geom, geom)

    def test_final_complex_check_uses_retained_mc_ground_not_selected_parent(self):
        parent, reactions = self.reactions()
        r = reactions[0]
        r.irc_prod_opt.species.energy = (sum(p.energy+p.zpe for p in r.products)
            - r.irc_prod_opt.species.zpe - 1./constants.AUtoKCAL)
        fragment = r.products[0]
        reference = representative_record(fragment)
        lower = replace(reference, index=2,
                        zero_energy_hartree=reference.zero_energy_hartree-4/constants.AUtoKCAL)
        fragment.conformer_index = [2]
        fragment.conformer_geom = [lower.geometry]
        fragment.conformer_zeroenergy = [lower.zero_energy_hartree]
        fragment.conformer_freq = [lower.frequencies_cm1]
        fragment.conformer_inventory = (lower,)
        fragment.conformer_records = {2: lower}
        reassess_product_complex(r, dict(self.par, multi_conf_tst=1))
        self.assertFalse(r.do_vdW)
        expected = (lower.zero_energy_hartree + r.products[1].energy + r.products[1].zpe
                    - r.irc_prod_opt.species.energy - r.irc_prod_opt.species.zpe)*constants.AUtoKCAL
        self.assertAlmostEqual(r.vdW_depth, expected)

    def test_bound_complex_with_large_imaginary_mode_keeps_direct_saddle(self):
        parent, reactions = self.reactions()
        reference = parent.energy + parent.zpe
        products = reactions[0].products
        products[0].energy += (reference + 43.19/constants.AUtoKCAL
                              - sum(p.energy + p.zpe for p in products))
        saved = []
        for reaction, barrier in zip(reactions, (52.01, 42.01)):
            point = reaction.irc_prod_opt.species
            point.energy = reference + 31.93/constants.AUtoKCAL - point.zpe
            point.freq = [-71.0308, -10.5093] + [1000.] * 7
            point.reduced_freqs = [-71.0308, 10.5093] + [1000.] * 7
            reaction.ts.energy = reference + barrier/constants.AUtoKCAL - reaction.ts.zpe
            saved.append((point, point.geom.copy(), list(point.freq), reaction.ts))
        writer = MESS(self.par, parent)
        with self.assertLogs('KinBot', level='WARNING') as logs:
            for reaction in reactions:
                reassess_product_complex(reaction, self.par)
            writer.write_input(None)
        output = Path('me/mess_0000.inp').read_text()
        headers = [line.split('!')[0].split() for line in output.splitlines()
                   if line.strip().startswith('Barrier ')]
        self.assertEqual(len(headers), 1)
        self.assertEqual(headers[0][1:3], [writer.ts_names['low'], 'w_1'])
        self.assertTrue(headers[0][3].startswith('b_'))
        self.assertEqual(sum(l.strip().startswith('Well ') for l in output.splitlines()), 1)
        self.assertNotIn('Core PhaseSpaceTheory', output)
        self.assertNotIn('Tunneling   Eckart', output)
        self.assertIn('barrier is submerged', output)
        self.assertIn('-71.0308', '\n'.join(logs.output))
        for reaction, (point, geom, raw, ts) in zip(reactions, saved):
            self.assertFalse(reaction.do_vdW)
            self.assertFalse(reaction.final_vdW_assessment['retained'])
            self.assertAlmostEqual(reaction.vdW_depth, 11.26)
            self.assertIs(reaction.irc_prod_opt.species, point)
            self.assertIs(reaction.ts, ts)
            self.assertEqual(point.freq, raw)
            np.testing.assert_array_equal(point.geom, geom)
        self.assertEqual(parent.reac_ts_done, [-1, -1])

    def test_complex_frequency_check_uses_existing_configured_tolerance(self):
        for mode, threshold, retained in ((100., 50., True), (-10.5, 50., True),
                                         (-50., 50., True), (-71., 50., False),
                                         (-71., 80., True)):
            with self.subTest(mode=mode, threshold=threshold):
                _, reactions = self.reactions()
                reaction = reactions[0]
                point = reaction.irc_prod_opt.species
                point.freq = [mode] + [1000.] * 8
                raw = list(point.freq)
                reassess_product_complex(reaction, dict(self.par, imagfreq_threshold=threshold))
                self.assertEqual(reaction.do_vdW, retained)
                self.assertEqual(point.freq, raw)

    def test_distinct_stereopaths_reach_products_and_keep_own_tunneling_reference(self):
        parent, reactions = self.reactions()
        original = reactions[0].irc_prod_opt.species
        par = dict(self.par, multi_conf_tst=1)
        with patch('kinbot.mess.reaction_path_id', side_effect=lambda r: r.instance_name):
            writer = MESS(par, parent)
            writer.write_input(None)
        text = Path('me/mess_0000.inp').read_text()
        self.assertIn('Union ! 2 stereochemical pathways', text)
        self.assertEqual(sum(line.strip().startswith('Well ') for line in text.splitlines()), 1)
        self.assertNotIn('Core PhaseSpaceTheory', text)
        # Master's inner-TS reference is retained; it is not a kinetic well.
        high = reactions[0].ts.mess_tunneling_reference
        self.assertAlmostEqual(high['right_parent_observations'][0]['electronic_energy_hartree'], original.energy)
        self.assertAlmostEqual(high['input_depths_kcal_mol'][1], 36.)
        from kinbot.mess_mirrors import population_keys
        keys = population_keys(reactions[0].products)
        self.assertIn(' '.join(keys), text)
        self.assertNotIn('IRC_F_prod', text)

    def test_unmodeled_complex_optics_cannot_block_direct_output(self):
        for mc in (0, 1):
            parent, reactions = self.reactions()
            complexes = [r.irc_prod_opt.species for r in reactions]
            writer = MESS(dict(self.par, multi_conf_tst=mc), parent)
            original = writer._parent_symmetry
            def checked(point, **kwargs):
                self.assertTrue(all(point is not p for p in complexes))
                return original(point, **kwargs)
            from kinbot.mess import writer_members as members
            def checked_members(point, *args, **kwargs):
                self.assertTrue(all(point is not p for p in complexes))
                return members(point, *args, **kwargs)
            with patch.object(writer, '_parent_symmetry', side_effect=checked), \
                    patch('kinbot.mess.writer_members', side_effect=checked_members):
                writer.write_input(None)
            self.assertTrue(all(r.do_vdW for r in reactions))

    def test_submerged_refined_ts_keeps_calculated_energies_without_complex_well(self):
        for mc in (0, 1):
            parent, reactions = self.reactions()
            reference = parent.energy + parent.zpe
            products = reactions[0].products
            products[0].energy += (reference + 60.11/constants.AUtoKCAL
                                  - sum(p.energy+p.zpe for p in products))
            for reaction, barrier in zip(reactions, (61.042, 57.340)):
                reaction.ts.energy = reference + barrier/constants.AUtoKCAL - reaction.ts.zpe
                point = reaction.irc_prod_opt.species
                point.energy = reference + 59.577/constants.AUtoKCAL - point.zpe
            writer = MESS(dict(self.par, multi_conf_tst=mc), parent)
            writer.write_input(None)
            text = Path('me/mess_0000.inp').read_text()
            # MC fragment RRHO blocks also carry a relative zero energy.
            self.assertEqual(values(text, 'ZeroEnergy')[-1], 57.34)
            self.assertTrue(all(value == 0. for value in values(text, 'ZeroEnergy')[:-1]))
            self.assertEqual(values(text, 'GroundEnergy'), [60.11])
            self.assertNotIn('Tunneling   Eckart', text)
            self.assertNotIn('Core PhaseSpaceTheory', text)
            self.assertNotIn('59.58', text)
            self.assertTrue(all(r.do_vdW for r in reactions))

    def test_pes_workers_still_write_complex_and_phase_space_models(self):
        from kinbot.mess_mirrors import population_keys
        from kinbot.species_routing import mess_filename
        for mc in (0, 1):
            parent, reactions = self.reactions()
            writer = MESS(dict(self.par, pes=1, multi_conf_tst=mc), parent)
            with patch('kinbot.mess.reaction_path_id', side_effect=lambda r: r.instance_name):
                writer.write_input(None)
            for reaction in reactions:
                name = reaction.irc_prod_opt.species.name
                complex_text = Path(name + '_0000.mess').read_text()
                self.assertRegex(complex_text, r'Well\s+\{name\}')
                self.assertIn('ZeroEnergy[kcal/mol]', complex_text)
                barrier = Path(reaction.instance_name + '_0000.mess').read_text()
                self.assertIn(' '.join(population_keys([reaction.irc_prod_opt.species])), barrier)
            key = '_'.join(sorted(routing_name(p) for p in reactions[0].products))
            text = Path(mess_filename(key, 0)).read_text()
            self.assertIn('Core PhaseSpaceTheory', text)
            self.assertIn('{wellname}', text)
            self.assertIn('{prodname}', text)

    def test_homolysis_never_counts_or_renders_placeholder_saddle(self):
        parent, reactions = self.reactions(True)
        writer = MESS(self.par, parent)
        original = writer._parent_symmetry
        def checked(p, **kwargs):
            self.assertIsNot(p, reactions[0].ts)
            return original(p, **kwargs)
        with patch.object(writer, '_parent_symmetry', side_effect=checked), \
             patch.object(writer, 'write_barrier', side_effect=AssertionError('placeholder TS')):
            writer.write_input(None)
        output = Path('me/mess_0000.inp').read_text()
        self.assertEqual(output.count('Core PhaseSpaceTheory'), 1)
        self.assertNotIn('Tunneling', output)
        self.assertNotIn('{blessname}', output)

    def test_mc_phase_space_uses_product_combinations_and_parent_normalization(self):
        parent, reactions = self.reactions(True)
        fragment = reactions[0].products[0]
        base = representative_record(fragment)
        # Associated calculated properties, including a tolerated soft mode.
        first = replace(base, index=0, zero_energy_hartree=base.zero_energy_hartree-1./constants.AUtoKCAL,
                        frequencies_cm1=(-20., *base.frequencies_cm1[1:]))
        second = replace(base, index=2, zero_energy_hartree=base.zero_energy_hartree+2./constants.AUtoKCAL,
                         sigma_ext=base.sigma_ext/2, remaining_optical_weight=2.)
        ensemble = {0: first, 2: second}
        from kinbot.mess import writer_members as original
        def members(point, population, **kwargs):
            return ensemble if point is fragment else original(point, population, **kwargs)
        par = dict(self.par, multi_conf_tst=1)
        with patch('kinbot.mess.writer_members', side_effect=members):
            writer = MESS(par, parent)
            writer.write_input(None)
            direct = Path('me/mess_0000.inp').read_text()
            writer.par['pes'] = 1
            deferred = writer.write_bimol(reactions[0].products, 0., 1., 1., 0, 1)
        expected_parent = (sum(p.energy+p.zpe for p in reactions[0].products)
                           - parent.energy-parent.zpe)*constants.AUtoKCAL
        replacements = dict(name='b_1', ground_energy=expected_parent, blessname='bl',
                            wellname='w_1', prodname='b_1')
        replacements.update({'fr_name_'+routing_name(p): 'fr_'+str(i)
                             for i,p in enumerate(reactions[0].products)})
        deferred = apply_conformer_shifts(deferred.format(**replacements))
        for text in (direct, deferred):
            barrier = text[text.index('  Barrier'):]
            self.assertIn('Union ! product conformer combinations', barrier)
            self.assertEqual(values(barrier, 'ZeroEnergy'),
                             [round(expected_parent-1., 2), round(expected_parent+2., 2)])
            self.assertEqual(values(barrier, 'SymmetryFactor'), [1., .25])
            self.assertIn('20.0', barrier)
            self.assertNotIn('-20.0', barrier)
            self.assertEqual(values(text, 'GroundEnergy'), [round(expected_parent-1., 2)])
        self.assertEqual(first.frequencies_cm1[0], -20.)


    def test_two_mc_fragment_ensembles_form_four_phase_space_terms(self):
        parent, _ = self.reactions(True)
        products = [species('H2O'), species('NH3')]
        choices = []
        for p, offset, divisor in zip(products, (2., 3.), (.5, 1./3.)):
            base = representative_record(p)
            choices.append({0: base, 1: replace(base, index=1,
                zero_energy_hartree=base.zero_energy_hartree+offset/constants.AUtoKCAL,
                sigma_ext=base.sigma_ext*divisor)})
        writer = MESS(dict(self.par, multi_conf_tst=1), parent)
        def members(p, population):
            return choices[next(i for i, point in enumerate(products) if point is p)]
        with patch('kinbot.mess.writer_members', side_effect=members):
            direct = writer._phase_space_models(products, 'H5N1O1', 50., 0., 1., 1.)
            writer.par['pes'] = 1
            deferred = writer._phase_space_models(products, 'H5N1O1', '{ground_energy}', 0., 1., 1.)
        deferred = apply_conformer_shifts(deferred.format(ground_energy=50.))
        for output in (direct, deferred):
            self.assertEqual(output.count('Core PhaseSpaceTheory'), 4)
            self.assertEqual(values(output, 'ZeroEnergy'), [50., 53., 52., 55.])
            self.assertEqual(values(output, 'SymmetryFactor'), [1., 1./3., .5, 1./6.])
            self.assertEqual(values(output, 'Frequencies'), [9.]*4)

    def test_pes_complex_choice_does_not_depend_on_discovery_order(self):
        observations = [ ['p', name, prods, 30., None, 'vdW_IRC_F_prod']
                         for name, prods in [('high', ['a', 'b']), ('low', ['b', 'a']),
                                             ('stereo', ['a', 'b-s-other'])] ]
        energies = {'high_IRC_F_prod': -1., 'low_IRC_F_prod': -2., 'stereo_IRC_F_prod': -3.}
        expected = {'high_IRC_F_prod': 'low_IRC_F_prod', 'low_IRC_F_prod': 'low_IRC_F_prod',
                    'stereo_IRC_F_prod': 'stereo_IRC_F_prod'}
        for order in itertools.permutations(observations):
            self.assertEqual(pes.find_min_vdW(list(order), energies), expected)


if __name__ == '__main__':
    unittest.main()
