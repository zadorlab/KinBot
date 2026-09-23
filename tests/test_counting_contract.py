"""Representation counting from real methanol evidence and small state models."""
import copy
import json
import pickle
import unittest
from unittest.mock import patch

import numpy as np

from kinbot import constants, frequencies, symmetry
from kinbot.conformer_records import inventory, retain
from kinbot.counting_contract import site_counting
from kinbot.stereo_identity import optical_scope
from kinbot.stationary_pt import StationaryPoint
from kinbot.thermochemistry import thermochemistry_evidence
from tests.counting_fixtures import methanol_data, saved_point
from tests.test_conformer_counting import peroxide


class TestCountingContract(unittest.TestCase):
    def methyl(self):
        return saved_point(methanol_data()['saddles'][0])

    def test_saved_both_channels_and_fragment_rotors(self):
        for run in ('MeOH_H', 'MeOH_H_rotor'):
            fixture = methanol_data(run)
            for data in fixture['saddles']:
                species = saved_point(data)
                np.testing.assert_array_equal(species.bond, data['bond'])
                record = thermochemistry_evidence(species)
                self.assertEqual(record['sigma_ext'], 1)
                methyl = data['instance'][1] == 3
                expected = 2 if methyl and run == 'MeOH_H' else 1
                self.assertEqual(record['optical_counting']['remaining_multiplier'], expected)
                self.assertEqual(record['source_row_id'], data['source_row_id'])
                self.assertEqual(record['raw_harmonic_frequencies_cm-1'], data['raw_frequencies_cm-1'])
                self.assertEqual(record['zpe_hartree'], data['zpe_hartree'])
                self.assertEqual(record['electronic_energy_hartree'], data['electronic_energy_hartree'])
                count = site_counting(data['site_degeneracy'], data['site_members'],
                                      data['site_relations'], data['stereopath_id'])
                self.assertEqual(count['unlabelled_model_conversion']['additional_site_multiplier'], 1)
                self.assertEqual(count['degeneracy'], 3 if methyl else 1)
                self.assertEqual(count['physical_path_completeness'], 'not_established')
                if methyl and run.endswith('_rotor'):
                    rotor = record['hir']['rotors'][0]
                    self.assertEqual(rotor['dihedral'], [4, 1, 2, 6])
                    self.assertEqual(rotor['sigma_int'], 1)
                    self.assertEqual(rotor['mirror_coverage']['witnesses'][0]['mirror_point'], 7)
                    self.assertEqual(len(record['thermochemical_frequencies_cm-1']), 14)
                json.dumps(record, allow_nan=False)
        fixture = methanol_data()
        methanol = thermochemistry_evidence(saved_point(fixture['fragments']['methanol']))
        self.assertEqual(methanol['hir']['rotors'][0]['sigma_int'], 3)
        self.assertEqual(len(methanol['thermochemical_frequencies_cm-1']), 11)
        ch2oh = thermochemistry_evidence(saved_point(fixture['fragments']['CH2OH']))
        self.assertEqual(ch2oh['hir']['rotors'][0]['sigma_int'], 2)
        self.assertTrue(ch2oh['hir']['rotors'][0]['usable'])
        self.assertEqual(ch2oh['optical_counting']['remaining_multiplier'], 1.)
        self.assertEqual(ch2oh['optical_counting']['status'], 'resolved')
        self.assertIsNone(ch2oh['hir']['rotors'][0]['potential_periodicity_verified'])
        self.assertEqual(len(ch2oh['thermochemical_frequencies_cm-1']), 8)
        self.assertEqual(saved_point(fixture['fragments']['CH3O']).sigma_ext, 3)

    def test_partial_scan_uses_measured_witnesses_but_not_fitted_points(self):
        p = self.methyl()
        h = p.hir
        h.hir_status[0][3] = 1
        h.hir_energies[0][3] = -1.
        self.assertTrue(h.fourier_fit('partial', np.arange(12) * np.pi / 6, 0))
        self.assertTrue(h.is_valid_rotor(0))
        self.assertEqual(thermochemistry_evidence(p)['optical_counting']['remaining_multiplier'], 1)
        h.hir_status[0][7] = 1
        h.hir_energies[0][7] = -1.
        self.assertTrue(h.fourier_fit('partial', np.arange(12) * np.pi / 6, 0))
        record = thermochemistry_evidence(p)
        self.assertTrue(record['hir']['rotors'][0]['usable'])
        self.assertEqual(record['optical_counting']['remaining_multiplier'], 1)
        self.assertNotIn('measured_coverage', record['optical_counting'])
        self.assertIn('coordinate_coverage', record['optical_counting'])
        self.assertIsNotNone(record['hir']['rotors'][0]['points'][7]['fitted_relative_energy_kcal_mol'])
        self.assertIsNotNone(record['hir']['rotors'][0]['points'][7]['observation']['electronic_energy_hartree'])

    def test_stale_or_absent_reference_does_not_certify_coverage(self):
        for change in ('scan', 'projection', 'row', 'geometry', 'hessian'):
            p = self.methyl()
            if change == 'scan': p.hir.scan_reference = None
            elif change == 'projection': p.rotor_projection.pop('reference')
            elif change == 'row': p.source_row_id += 1
            elif change == 'geometry': p.geom[6, 2] += .2
            else: p.hess[0][0] += .1
            self.assertIsNone(thermochemistry_evidence(p)['optical_counting']['remaining_multiplier'])

    def test_scan_energy_comparison_does_not_mix_selected_and_hir_levels(self):
        p = self.methyl()
        p.energy -= 1.  # A different selected-level zero is not a scan contradiction.
        self.assertEqual(thermochemistry_evidence(p)['optical_counting']['remaining_multiplier'], 1)
        p.hir.hir_raw_energies[0][7] += 1. / constants.AUtoKCAL
        record = thermochemistry_evidence(p)
        self.assertEqual(record['optical_counting']['remaining_multiplier'], 1.)
        self.assertTrue(record['optical_counting']['warnings'])

    def test_measured_hir_mirror_accepts_sub_kcal_energy_disagreement(self):
        for delta, resolved in ((.5, True), (.99, True), (1.01, False)):
            with self.subTest(delta_kcal_mol=delta):
                p = self.methyl()
                p.hir.hir_raw_energies[0][7] = (p.hir.hir_raw_energies[0][0]
                                               + delta / constants.AUtoKCAL)
                counting = thermochemistry_evidence(p)['optical_counting']
                self.assertEqual(counting['mirror_energy_tolerance_kcal_mol'], 1.)
                self.assertEqual(counting['status'] == 'resolved', resolved)
                self.assertEqual(counting['remaining_multiplier'], 1.)
                self.assertEqual(bool(counting.get('warnings')), not resolved)

    def test_full_turn_symmetry_quotient_does_not_repeat_optical_weight(self):
        p = self.methyl()
        axis = p.dihed[0][1:3]
        p.sigma_int[min(axis)][max(axis)] = 3
        record = thermochemistry_evidence(p)
        self.assertTrue(record['hir']['rotors'][0]['mirror_coverage']['mirror_observed_in_full_scan'])
        self.assertEqual(record['optical_counting']['remaining_multiplier'], 1.)
        coverage = record['hir']['rotors'][0]['mirror_coverage']
        self.assertEqual(coverage['counting_domain_degrees'], [0., 360.])
        self.assertEqual(coverage['potential_period_degrees'], 120.)
        self.assertEqual(coverage['symmetry_quotient'], 3)
        self.assertFalse(coverage['symmetry_periodicity_verified'])
        self.assertTrue(any(not w['within_serialized_potential_period'] for w in coverage['witnesses']))
        p = self.methyl()
        p.dihed.append(list(p.dihed[0]))
        p.hir.scan_reference['dihedrals'] = copy.deepcopy(p.dihed)
        for attr in ('hir_status', 'hir_energies', 'hir_raw_energies', 'hir_geoms',
                     'hir_fourier', 'hir_fit_diagnostics', 'scan_jobs', 'point_observations'):
            getattr(p.hir, attr).append(copy.deepcopy(getattr(p.hir, attr)[0]))
        p.hir.projection_status.append({'rotor_index': 1, 'projected': True, 'reason': None})
        p.rotor_projection['reference']['dihedrals'] = copy.deepcopy(p.dihed)
        p.rotor_projection['internal_rank'] = 2
        p.reduced_freqs = p.reduced_freqs[:-1]
        record = thermochemistry_evidence(p)
        self.assertIsNone(record['optical_counting']['remaining_multiplier'])
        self.assertIn('multiple', record['optical_counting']['reason'])

    def test_projection_decisions_and_current_frequencies_must_agree(self):
        for change in ('missing', 'duplicate', 'rank', 'axes', 'demoted', 'frequencies'):
            with self.subTest(change=change):
                p = self.methyl()
                if change == 'missing': p.rotor_projection['rotors'] = []
                elif change == 'duplicate': p.rotor_projection['rotors'] *= 2
                elif change == 'rank': p.rotor_projection['internal_rank'] = 0
                elif change == 'axes': p.rotor_projection['reference']['dihedrals'][0][1:3] = [0, 1]
                elif change == 'demoted': p.hir.hir_status[0] = [1] * 12
                else: p.reduced_freqs = list(p.freq)
                before = pickle.dumps(vars(p))
                result = thermochemistry_evidence(p)
                self.assertEqual(result['optical_counting']['status'], 'unresolved')
                self.assertIsNone(result['optical_counting']['remaining_multiplier'])
                self.assertEqual(len(result['hir']['rotors'][0]['points']), 12)
                self.assertEqual(pickle.dumps(vars(p)), before)

    def test_demoted_rotor_is_a_harmonic_fallback(self):
        p = self.methyl()
        p.hir.hir_status[0] = [1] * 12
        _, p.reduced_freqs = frequencies.get_frequencies(p, p.hess, p.geom)
        record = thermochemistry_evidence(p)
        self.assertEqual(record['optical_counting']['representation'], 'RRHO')
        self.assertEqual(record['optical_counting']['remaining_multiplier'], 2)
        self.assertEqual(len(record['thermochemical_frequencies_cm-1']), 15)

    def test_unknown_scope_and_forbidden_measured_mirror_are_distinct(self):
        p = self.methyl()
        scope = optical_scope(copy.copy(p))
        scope['mirror_allowed'] = False
        # The same geometry/scan under a specified chiral input population.
        with patch('kinbot.optical.optical_scope', return_value=scope):
            result = thermochemistry_evidence(p)
        self.assertIn('outside the specified', result['optical_counting']['reason'])
        self.assertIsNone(result['optical_counting']['remaining_multiplier'])
        self.assertTrue(result['hir']['rotors'][0]['usable'])
        p.optical_reference = {'status': 'unavailable'}
        result = thermochemistry_evidence(p)['optical_counting']
        self.assertEqual(result['status'], 'legacy_unverified')
        self.assertEqual(result['remaining_multiplier'], p.nopt)

    def mc(self, geoms):
        p = peroxide()
        p.energy, p.zpe = -1.01, .01
        p.freq = p.reduced_freqs = [100.] * 6
        p.conformer_representation = 'MC-RRHO'
        p.conformer_geom = geoms
        p.conformer_index = list(range(len(geoms)))
        p.conformer_zeroenergy = [-1.] * len(geoms)
        p.conformer_freq = [[100.] * 6 for _ in geoms]
        return p

    def test_explicit_mc_mirrors_and_archives_are_counted_without_export_side_effects(self):
        geom = peroxide().geom
        for explicit in (False, True):
            p = self.mc([geom, geom * [-1, 1, 1]] if explicit else [geom])
            observations = inventory(p, [geom, geom * [-1, 1, 1]], [-1., -1.], [[100.] * 6] * 2, [0, 0])
            retain(p, observations, p.conformer_index)
            before = pickle.dumps(vars(p))
            with patch('kinbot.conformer_counting.preserve_counting_error', side_effect=AssertionError('no writes')):
                first, second = thermochemistry_evidence(p), thermochemistry_evidence(p)
            self.assertEqual(first, second)
            self.assertEqual(pickle.dumps(vars(p)), before)
            counting = first['optical_counting']
            self.assertIsNone(counting['remaining_multiplier'])
            self.assertEqual([m['remaining_optical_weight'] for m in counting['members']], [1., 1.] if explicit else [2.])
            self.assertEqual(len(first['conformer_inventory']), 2)
            if explicit:
                self.assertEqual(counting['members'][0]['mirror_partner_ids'], [counting['members'][1]['member_id']])

    def test_mixed_mc_weights_and_unsupported_populations_remain_explicit(self):
        geom = peroxide().geom
        other = geom.copy()
        other[3] = [1.75, -.65, .62]
        p = self.mc([geom, geom * [-1, 1, 1], other])
        result = thermochemistry_evidence(p)['optical_counting']
        self.assertEqual(len(result['members']), 3)
        self.assertEqual([m['remaining_optical_weight'] for m in result['members']], [1., 1., 2.])
        self.assertIsNone(result['remaining_multiplier'])
        p.optical_reference = {'status': 'unsupported', 'reason': 'unassigned configuration'}
        before = pickle.dumps(vars(p))
        with patch('kinbot.conformer_counting.preserve_counting_error', side_effect=AssertionError('no writes')):
            result = thermochemistry_evidence(p)
        self.assertEqual(pickle.dumps(vars(p)), before)
        self.assertEqual(len(result['conformer_inventory']), 3)
        self.assertIsNone(result['optical_counting']['remaining_multiplier'])
        self.assertIsNone(result['conformer_counting_error'])
        self.assertTrue(all(member['optical_evidence']['status'] == 'legacy_unverified'
                            for member in result['optical_counting']['members']))

    def test_missing_member_arrays_are_preserved_with_an_error(self):
        p = self.mc([peroxide().geom])
        p.conformer_zeroenergy = []
        before = pickle.dumps(vars(p))
        result = thermochemistry_evidence(p)
        self.assertEqual(pickle.dumps(vars(p)), before)
        self.assertTrue(result['conformer_counting_error'])
        self.assertEqual(result['unassociated_conformer_arrays']['conformer_index'], [0])
        self.assertIsNone(result['optical_counting']['remaining_multiplier'])

    def test_duplicate_indices_and_extra_member_properties_are_not_silently_joined(self):
        for change in ('duplicate', 'extra_energy', 'extra_geometry'):
            p = self.mc([peroxide().geom, peroxide().geom * [-1, 1, 1]])
            if change == 'duplicate': p.conformer_index = [0, 0]
            elif change == 'extra_energy': p.conformer_zeroenergy.append(-2.)
            else: p.conformer_geom.append(peroxide().geom)
            before = pickle.dumps(vars(p))
            result = thermochemistry_evidence(p)
            self.assertTrue(result['conformer_counting_error'])
            self.assertEqual(result['optical_counting']['status'], 'unresolved')
            self.assertIsNotNone(result['unassociated_conformer_arrays'])
            self.assertEqual(pickle.dumps(vars(p)), before)

    def test_nonfinite_member_energy_cannot_certify_mirror_compatibility(self):
        from kinbot.conformer_counting import writer_members, CountingError
        for energy in (float('nan'), float('inf'), -float('inf')):
            p = self.mc([peroxide().geom, peroxide().geom * [-1, 1, 1]])
            p.conformer_zeroenergy[1] = energy
            before = pickle.dumps(vars(p))
            result = thermochemistry_evidence(p)
            self.assertEqual(pickle.dumps(vars(p)), before)
            self.assertEqual(result['optical_counting']['status'], 'unresolved')
            self.assertIn('finite E + ZPE', result['conformer_counting_error'])
            self.assertEqual(len(result['conformer_inventory']), 2)
            self.assertFalse(result['optical_counting']['members'])
            with patch('kinbot.conformer_counting.preserve_counting_error') as preserve:
                with self.assertRaises(CountingError):
                    writer_members(p)
                preserve.assert_called_once()

    def test_specified_configuration_and_global_racemic_pair_scope(self):
        p = StationaryPoint('chiral', 0, 1, smiles='F[C@](Cl)(Br)I')
        p.characterize()
        p.energy, p.zpe = -100., .01
        p.freq = p.reduced_freqs = [100.] * 9
        symmetry.calculate_symmetry(p)
        specified = thermochemistry_evidence(p)['optical_counting']
        self.assertEqual(specified['remaining_multiplier'], 1)
        p.optical_population = 'racemic'
        racemic = thermochemistry_evidence(p)['optical_counting']
        self.assertEqual(racemic['remaining_multiplier'], 2)
        self.assertIsNone(racemic['total_configurational_states'])
        self.assertIn('not independent fragment racemates', racemic['scope'])
        pair = StationaryPoint('pair', 0, 1, atom=list(p.atom) * 2,
                               geom=np.vstack((p.geom, p.geom + [10., 0., 0.])))
        pair.characterize()
        pair.energy, pair.zpe = -200., .02
        pair.freq = pair.reduced_freqs = [100.] * 24
        pair.optical_population = 'racemic'
        scope = thermochemistry_evidence(pair)['optical_counting']
        self.assertEqual(scope['allowed_global_mirror_states'], 2)
        self.assertIsNone(scope['total_configurational_states'])  # Not all RR/RS/SR/SS populations.

    def test_same_model_explicit_and_representative_counts_agree(self):
        from scipy.integrate import quad
        harmonic = saved_point(methanol_data('MeOH_H')['saddles'][0])
        harmonic_weight = thermochemistry_evidence(harmonic)['optical_counting']['remaining_multiplier']
        rotor_weight = thermochemistry_evidence(self.methyl())['optical_counting']['remaining_multiplier']
        for temperature in (200., 500., 1500.):
            beta = 1. / (.001987204258 * temperature)
            one_saddle = np.exp(-7. * beta)
            # Two distinct rigid mirror structures, independently of labels.
            self.assertAlmostEqual(harmonic_weight * one_saddle, sum([one_saddle, one_saddle]))
            reactant = quad(lambda a: np.exp(-beta * .5 * (1 - np.cos(3*a))), 0, 2*np.pi)[0] / 3
            density = lambda a: np.exp(-beta * 2 * (1 - np.cos(2*a)))
            full = quad(density, 0, 2*np.pi)[0]
            explicit = quad(density, 0, np.pi)[0] + quad(density, np.pi, 2*np.pi)[0]
            self.assertAlmostEqual(rotor_weight * full / reactant, explicit / reactant)
            # O-H has one rigid TS and no extra optical or site factor.
            oh = saved_point(methanol_data()['saddles'][1])
            self.assertEqual(thermochemistry_evidence(oh)['optical_counting']['remaining_multiplier'], 1)
        point = self.methyl()
        data = methanol_data()['saddles'][0]
        args = [data[k] for k in ('site_degeneracy', 'site_members', 'site_relations', 'stereopath_id')]
        original = site_counting(*args)
        point.sigma_int[1][2] = 7  # A symmetry divisor is not a label-class proof.
        self.assertEqual(site_counting(*args), original)
        methane = site_counting(4, {'hydrogen': [1, 2, 3, 4]}, {'hydrogen': 'homotopic'}, 'methane')
        self.assertEqual((12/3) * methane['unlabelled_model_conversion']['additional_site_multiplier'], 4)

    def test_site_conversion_does_not_claim_completeness_or_accept_conflicts(self):
        base = (3, {'donor': [2, 3, 4]}, {'donor': 'homotopic'}, 'site')
        count = site_counting(*base)
        self.assertEqual(count['unlabelled_model_conversion']['labelled_replica_normalization'], 1/3)
        self.assertNotIn('remaining_rate_weight', count)
        for args, extra in [(base, {'conflict': True}), ((4, *base[1:]), {}),
                            ((3, None, None, None), {})]:
            self.assertIsNone(site_counting(*args, **extra)['unlabelled_model_conversion']['additional_site_multiplier'])

    def test_malformed_site_members_are_preserved_without_certifying_a_conversion(self):
        for atoms in ('123', 3, [True], [-1], [[1], [2, 3]], [1, 1, 2], [1., 2., 3.]):
            result = site_counting(3, {'H': atoms}, {'H': 'homotopic'}, 'site')
            self.assertIsNone(result['unlabelled_model_conversion']['additional_site_multiplier'])
            self.assertEqual(result['site_members']['H'], atoms)
        result = site_counting(np.int64(3), {'H': np.array([1, 2, 3])}, {'H': 'homotopic'}, 'site')
        self.assertEqual(result['unlabelled_model_conversion']['additional_site_multiplier'], 1)


if __name__ == '__main__':
    unittest.main()
