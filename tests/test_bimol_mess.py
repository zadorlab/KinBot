"""Saved methanol calculations and analytic limits of bimolecular MESS models."""
import copy
import math
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
from tempfile import TemporaryDirectory
import unittest
from types import SimpleNamespace

import numpy as np

from kinbot import constants
from kinbot.bimol_mess import BimolecularMESS, RecordRenderer
from kinbot.bimol_rate_models import ensemble, population_route, prepare_models
from kinbot.parameters import Parameters
from kinbot.bimol_network import BimolecularNetwork, _species_data
from tests.test_bimol_counting_contract import saved_network
from tests.counting_fixtures import methanol_data, saved_point


def parameters(directory):
    path = Path(directory) / 'input.json'
    path.write_text(json.dumps(dict(barrier_threshold=200., me=0, high_level=0)))
    par = Parameters(str(path)).par
    par.update(TemperatureList=[300, 500, 1000], PressureList=[760],
               rotor_scan=1, high_level=0, epsilon=100., sigma=3.,
               queuing='local', me=1, run_me=1, uq_n=1, me_skip_vdW=0,
               mess_command=shutil.which('mess') or 'mess',
               mess_abstraction_command=shutil.which('mess-abstraction') or 'mess-abstraction')
    return par


def methanol_records(mc=False, hir=True):
    network = saved_network()
    network.par.update(multi_conf_tst=int(mc), rotor_scan=int(hir and not mc))
    points = [network.species.fragA, network.species.fragB]
    for reaction in network.species.reac_obj:
        # The portable fixture predates the shared pathway annotation. Supply
        # its already saved site class; do not infer a new reaction here.
        reaction.ts.stereopath_id = reaction.stereopath_id
        reaction.ts.stereopath_metadata = {'id': reaction.stereopath_id}
        points.extend([reaction.ts, *(opt.species for opt in reaction.prod_opt)])
    if mc or not hir:
        for point in points:
            point.hir = None
            point.rotor_projection = None
            point.reduced_freqs = list(point.freq)
            point.conformer_representation = 'MC-RRHO' if mc else 'selected conformer'
    return network.as_dict()


def rates(path):
    return np.array([[float(word) for word in line.split()] for line in path.read_text().splitlines()
                     if line.strip() and line.split()[0].isdigit()])


def add_complex(records, identifier, role, endpoint, energy):
    """A positive-mode analytic well for testing MESS, not a computed complex."""
    data = copy.deepcopy(records['rate_models']['pathways'][0]['stationary_point'])
    data['name'], data['energy'] = identifier, energy - data['zpe']
    thermo = data['thermochemistry']
    data['frequencies'] = [abs(f) for f in thermo['raw_harmonic_frequencies_cm-1']]
    thermo.update(electronic_energy_hartree=data['energy'],
                  hir={'rotors': []}, frequency_projection=None,
                  conformer_representation='selected conformer')
    thermo['raw_harmonic_frequencies_cm-1'] = data['frequencies']
    thermo['thermochemical_frequencies_cm-1'] = data['frequencies']
    thermo['optical_counting']['remaining_multiplier'] = 1.
    records['rate_models']['complexes'].append(dict(id=identifier, role=role,
        asymptote=endpoint, stationary_point=data, sampling='analytic well fixture'))


def reverse_atoms(data):
    order = list(reversed(range(len(data['atom']))))
    inverse = {old: new for new, old in enumerate(order)}
    for key in ('atom', 'geometry'):
        data[key] = [data[key][i] for i in order]
    data['bond'] = np.asarray(data['bond'])[np.ix_(order, order)].tolist()
    thermo = data['thermochemistry']
    thermo['atoms'], thermo['geometry_angstrom'] = data['atom'], data['geometry']
    members = list(thermo['conformer_inventory'])
    if thermo.get('rrho_representative_counting'):
        members.append(thermo['rrho_representative_counting'])
    for member in members:
        if member.get('geometry'):
            member['geometry'] = [member['geometry'][i] for i in order]
    for rotor in thermo['hir']['rotors']:
        for key in ('axis', 'dihedral'):
            rotor[key] = [inverse[i] for i in rotor[key]]
        rotor['partition'] = [[inverse[i] for i in part] for part in rotor['partition']]


class TestBimolecularMESS(unittest.TestCase):
    def test_invalid_product_frequencies_do_not_remove_unrelated_channel(self):
        previous = Path.cwd()
        with TemporaryDirectory() as directory:
            try:
                records = methanol_records(hir=False)
                path = records['rate_models']['pathways'][0]
                node = next(n for n in records['nodes'] if n['id'] == path['product'])
                point = next(d for d in node['fragment_data'] if len(d['atom']) > 2)
                point['thermochemistry']['thermochemical_frequencies_cm-1'][0] = -1000.
                os.chdir(directory)
                writer = BimolecularMESS(parameters(directory), records)
                writer.write_input()
                self.assertEqual(len(writer.mess_jobs), 1)
                self.assertTrue(any('invalid harmonic frequencies' in w for w in writer.warnings))
            finally:
                os.chdir(previous)

    def test_submerged_reverse_endpoint_does_not_receive_direct_tst(self):
        previous = Path.cwd()
        with TemporaryDirectory() as directory:
            try:
                par = parameters(directory)
                records = methanol_records(hir=False)
                path = records['rate_models']['pathways'][0]
                node = next(n for n in records['nodes'] if n['id'] == path['product'])
                for data in node['fragment_data']:
                    data['energy'] += 1.  # Unsupported energy ordering, deliberately large.
                os.chdir(directory)
                writer = BimolecularMESS(par, records)
                writer.write_input()
                self.assertEqual(len(writer.mess_jobs), 1)
                self.assertTrue(any('at or below a separated endpoint' in warning for warning in writer.warnings))
            finally:
                os.chdir(previous)

    @unittest.skipUnless(shutil.which('mess-abstraction'), 'Argonne abstraction companion is not installed')
    def test_fragment_and_atom_order_leave_all_three_rate_models_unchanged(self):
        previous = Path.cwd()
        with TemporaryDirectory() as directory:
            try:
                par = parameters(directory)
                for mode in ({'hir': False}, {}, {'mc': True}):
                    original = methanol_records(**mode)
                    permuted = copy.deepcopy(original)
                    for node in permuted['nodes']:
                        node['fragment_data'].reverse()
                        for point in node['fragment_data']:
                            reverse_atoms(point)
                    for path in permuted['rate_models']['pathways']:
                        reverse_atoms(path['stationary_point'])
                    observed = []
                    for i, records in enumerate((original, permuted)):
                        root = Path(directory, str(len(observed)) + str(mode))
                        root.mkdir(exist_ok=True)
                        os.chdir(root)
                        writer = BimolecularMESS(par, records)
                        writer.write_input()
                        writer.run()
                        self.assertFalse(writer.warnings)
                        observed.append({tuple(job['pathways']): rates(Path('me', job['stem'] + '.out'))
                                         for job in writer.mess_jobs})
                    self.assertEqual(set(observed[0]), set(observed[1]))
                    for path in observed[0]:
                        np.testing.assert_allclose(observed[0][path], observed[1][path], rtol=1.e-6)
            finally:
                os.chdir(previous)

    def test_distinct_stereochemical_paths_are_summed_without_enabling_mc(self):
        previous = Path.cwd()
        with TemporaryDirectory() as directory:
            try:
                records = methanol_records()
                first = records['rate_models']['pathways'][0]
                second = copy.deepcopy(first)
                second.update(id='other_stereochemical_path', path_id='other_stereochemical_class')
                records['rate_models']['pathways'].append(second)
                os.chdir(directory)
                writer = BimolecularMESS(parameters(directory), records)
                writer.write_input()
                job = next(j for j in writer.mess_jobs if first['id'] in j['pathways'])
                self.assertEqual(set(job['pathways']), {first['id'], second['id']})
                text = Path('me', job['stem'] + '.inp').read_text()
                self.assertIn('Barrier abstraction\nUnion', text)
                self.assertEqual(text.count('Tunneling   Eckart'), 2)
            finally:
                os.chdir(previous)
    def test_racemic_routes_preserve_relative_configurations(self):
        r = (('R',), ('S',))
        s = (('S',), ('R',))
        self.assertEqual(population_route(r, r), population_route(s, s))
        self.assertEqual(population_route(r, s), population_route(s, r))
        self.assertNotEqual(population_route(r, r), population_route(r, s))

    def test_specified_chiral_entrance_does_not_create_other_enantiomer(self):
        records = methanol_records()
        reactant = next(n for n in records['nodes'] if n['id'] == records['reference_node'])
        reactant['fragment_data'][0]['stereochemical_identity'].update(id='R', mirror_id='S')
        product = next(n for n in records['nodes'] if n['id'] == records['rate_models']['pathways'][0]['product'])
        product['fragment_data'][0]['stereochemical_identity'].update(id='P_R', mirror_id='P_S')
        with TemporaryDirectory() as directory:
            writer = BimolecularMESS(parameters(directory), records)
            self.assertEqual(len(writer.models['pathways']), 2)

    def test_achiral_entrance_completes_specified_mirror_product(self):
        records = methanol_records()
        product = next(n for n in records['nodes'] if n['id'] == records['rate_models']['pathways'][0]['product'])
        product['fragment_data'][0]['stereochemical_identity'].update(id='P_R', mirror_id='P_S')
        with TemporaryDirectory() as directory:
            writer = BimolecularMESS(parameters(directory), records)
            self.assertEqual(len(writer.models['pathways']), 3)
            mirror = next(p for p in writer.models['pathways'] if p['id'].endswith('_mirror'))
            self.assertTrue(mirror['stationary_point']['reflected_model'])
            original = next(p for p in writer.models['pathways'] if p['id'] == mirror['id'][:-7])
            self.assertEqual(writer._depths(original), writer._depths(mirror))

    def test_saved_hir_and_mc_models_use_their_own_frequency_treatment(self):
        with TemporaryDirectory() as directory:
            par = parameters(directory)
            for mc in (False, True):
                data = methanol_records(mc=mc)
                writer = BimolecularMESS(par, data)
                self.assertEqual(len(writer.models['pathways']), 2)
                for path in writer.models['pathways']:
                    text = writer.rrho(path['stationary_point'], writer.reference,
                                       depths=writer._depths(path))
                    if mc:
                        self.assertNotIn('Rotor     ', text)
                    else:
                        expected = sum(r['usable'] for r in path['stationary_point']['thermochemistry']['hir']['rotors'])
                        self.assertEqual(text.count('Rotor     '), expected)
                    self.assertIn('Tunneling   Eckart', text)
                    self.assertEqual(path['site_counting']['unlabelled_model_conversion']
                                     ['additional_site_multiplier'], 1.)

    def test_conformer_union_deduplicates_atom_permutations(self):
        point = saved_point(methanol_data()['fragments']['methanol'])
        point.hir, point.rotor_projection = None, None
        point.reduced_freqs = list(point.freq)
        other = copy.deepcopy(point)
        order = list(reversed(range(point.natom)))
        other.atom, other.geom = np.asarray(other.atom)[order], other.geom[order]
        other.bond = other.bond[np.ix_(order, order)]
        other.bonds = [b[np.ix_(order, order)] for b in other.bonds]
        other.rads = [np.asarray(r)[order] for r in other.rads]
        selected = ensemble([point, other], 'specified')
        self.assertEqual(len(selected.conformer_records), 1)
        self.assertEqual(len(selected.conformer_inventory), 2)

    def test_only_lowest_complex_is_used_with_harmonic_hir_and_mc(self):
        network = saved_network()
        first = copy.copy(network.species.reac_obj[0])
        second = copy.copy(first)
        point = saved_point(methanol_data()['fragments']['methanol'])
        point.hir, point.rotor_projection = None, None
        point.reduced_freqs = list(point.freq)
        other = copy.deepcopy(point)
        point.name, other.name = 'first_complex', 'second_complex'
        point.energy, point.zpe = -100., .01
        # Higher electronic energy, but lower E + ZPE: this must be selected.
        other.energy, other.zpe = -99.999, .005
        first.do_reac_vdW = second.do_reac_vdW = True
        first.do_vdW = second.do_vdW = True
        first.irc_reactant_opt = SimpleNamespace(species=point)
        second.irc_reactant_opt = SimpleNamespace(species=other)
        first.irc_prod_opt = SimpleNamespace(species=copy.deepcopy(other))
        second.irc_prod_opt = SimpleNamespace(species=copy.deepcopy(point))
        network.species.reac_obj = [first, second]
        network.species.reac_ts_done = [-1, -1]
        for transition in network.transitions:
            reaction = transition['reaction']
            reaction.ts.stereopath_id = reaction.stereopath_id
            reaction.ts.stereopath_metadata = {'id': reaction.stereopath_id}
        expected = _species_data(other)
        with TemporaryDirectory() as directory:
            renderer = RecordRenderer(parameters(directory), None)
            for mc, hir in ((False, False), (False, True), (True, False)):
                for reverse in (False, True):
                    with self.subTest(mc=mc, hir=hir, reverse=reverse):
                        network.par.update(multi_conf_tst=int(mc), rotor_scan=int(hir))
                        network.species.reac_obj = [second, first] if reverse else [first, second]
                        models = prepare_models(network)
                        self.assertEqual({c['role'] for c in models['complexes']}, {'rc', 'pc'})
                        self.assertEqual(len(models['complexes']), 2)
                        for complex_model in models['complexes']:
                            data = complex_model['stationary_point']
                            self.assertEqual(data, expected)
                            self.assertEqual(len(renderer.members(data)), 1)
                            self.assertNotIn('Union', renderer.rrho(data, point.energy + point.zpe))
                        self.assertEqual(len(models['pathways']), 2)
                        self.assertTrue(all(p['stationary_point']['thermochemistry']
                            ['conformer_representation'] == ('MC-RRHO' if mc else 'selected conformer')
                            for p in models['pathways']))

    def test_complex_energy_disagreement_preserves_distinct_reaction_paths(self):
        network = saved_network()
        # Use the saved methanol + H endpoint geometry with positive test
        # frequencies; these complex energies are deliberately synthetic.
        first = network.species.reac_obj[0]
        for index, reaction in enumerate(network.species.reac_obj):
            point = copy.deepcopy(first.irc_reactant)
            point.name = f'complex_{index}'
            point.energy = network.reference - (1. + .5 * index) / constants.AUtoKCAL
            point.zpe = 0.
            point.freq = [100.] * (3 * point.natom - 6)
            point.reduced_freqs = list(point.freq)
            reaction.do_reac_vdW = True
            reaction.irc_reactant_opt = SimpleNamespace(species=point)
            reaction.ts.stereopath_id = reaction.stereopath_id
            reaction.ts.stereopath_metadata = {'id': reaction.stereopath_id}
        network = BimolecularNetwork(network.species, network.par)
        data = network.as_dict()
        complexes = [node for node in data['nodes'] if node['kind'] == 'complex']
        self.assertEqual(len(complexes), 2)
        self.assertTrue(all(node['energy_conflicts'] for node in complexes))
        self.assertEqual(len(data['rate_models']['complexes']), 1)
        selected = data['rate_models']['complexes'][0]
        self.assertEqual(selected['stationary_point']['name'], 'complex_1')
        self.assertEqual(len(data['rate_models']['pathways']), 2)
        self.assertFalse(data['rate_models']['warnings'])
        for path in data['rate_models']['pathways']:
            self.assertEqual(path['reactant_complex'], selected['id'])
            self.assertEqual(path['site_counting']['unlabelled_model_conversion']
                             ['additional_site_multiplier'], 1.)

    def test_repeated_ts_with_different_complex_energies_has_one_rate_path(self):
        for role in ('rc', 'pc'):
            for mc, hir in ((False, False), (False, True), (True, False)):
                with self.subTest(role=role, mc=mc, hir=hir):
                    network = saved_network()
                    first = network.species.reac_obj[0]
                    first.ts.stereopath_id = first.stereopath_id
                    first.ts.stereopath_metadata = {'id': first.stereopath_id}
                    second = copy.copy(first)
                    second.instance_name += '_repeat'
                    endpoint = first.irc_reactant if role == 'rc' else first.irc_prod
                    for index, reaction in enumerate((first, second)):
                        point = copy.deepcopy(endpoint)
                        point.name = f'complex_{index}'
                        point.energy = network.reference - (1. + .5 * index) / constants.AUtoKCAL
                        point.zpe = 0.
                        point.freq = [100.] * (3 * point.natom - 6)
                        point.reduced_freqs = list(point.freq)
                        setattr(reaction, 'do_reac_vdW' if role == 'rc' else 'do_vdW', True)
                        setattr(reaction, 'irc_reactant_opt' if role == 'rc' else 'irc_prod_opt',
                                SimpleNamespace(species=point))
                    network.species.reac_obj = [first, second]
                    network.species.reac_ts_done = [-1, -1]
                    network.par.update(multi_conf_tst=int(mc), rotor_scan=int(hir))
                    network = BimolecularNetwork(network.species, network.par)
                    data = network.as_dict()
                    self.assertEqual(len(data['transitions']), 2)
                    self.assertTrue(all(not t['endpoint_conflicts'] and not t['energy_conflicts']
                                        for t in data['transitions']))
                    complexes = [node for node in data['nodes'] if node['kind'] == 'complex']
                    self.assertEqual(len(complexes), 2)
                    self.assertTrue(all(node['energy_conflicts'] for node in complexes))
                    models = data['rate_models']
                    self.assertFalse(models['warnings'])
                    self.assertEqual(len(models['pathways']), 1)
                    self.assertEqual(len(models['pathways'][0]['observations']), 2)
                    self.assertEqual(len(models['complexes']), 1)
                    self.assertEqual(models['complexes'][0]['stationary_point']['name'], 'complex_1')

                    # A disagreement about equivalent reacting sites must
                    # remain an error even when the lower TS record is valid.
                    second.ts = copy.deepcopy(first.ts)
                    second.ts.energy += 1.e-8
                    second.site_degeneracy += 1
                    conflicting = BimolecularNetwork(network.species, network.par)
                    self.assertFalse(prepare_models(conflicting)['pathways'])
                    second.site_degeneracy = first.site_degeneracy

                    # Agreeing sites cannot excuse conflicting TS energies.
                    second.ts.energy = first.ts.energy + .5 / constants.AUtoKCAL
                    conflicting = BimolecularNetwork(network.species, network.par)
                    self.assertTrue(all(t['energy_conflicts'] for t in conflicting.transitions))
                    self.assertFalse(prepare_models(conflicting)['pathways'])

    def test_unprojected_rotor_is_not_written_twice(self):
        data = methanol_records()
        point = next(n['fragment_data'][1] for n in data['nodes'] if n['id'] == data['reference_node'])
        if len(point['atom']) == 1:
            point = next(n['fragment_data'][0] for n in data['nodes'] if n['id'] == data['reference_node'])
        point = copy.deepcopy(point)
        thermo = point['thermochemistry']
        thermo['thermochemical_frequencies_cm-1'] = thermo['raw_harmonic_frequencies_cm-1']
        for rotor in thermo['frequency_projection']['rotors']:
            rotor['projected'] = False
        with TemporaryDirectory() as directory:
            writer = RecordRenderer(parameters(directory), None)
            self.assertNotIn('Rotor     Hindered', writer.rrho(point, point['energy'] + point['zpe']))

    def test_wrong_units_and_missing_projected_rotor_are_not_silently_used(self):
        point = _species_data(saved_point(methanol_data()['fragments']['methanol']))
        bad = copy.deepcopy(point)
        bad['energy_unit'] = 'eV'
        with self.assertRaisesRegex(ValueError, 'hartree'):
            RecordRenderer.members(bad)
        bad = copy.deepcopy(point)
        bad['thermochemistry']['hir']['rotors'][0]['usable'] = False
        with self.assertRaisesRegex(ValueError, 'projected frequencies'):
            RecordRenderer.members(bad)

    @unittest.skipUnless(shutil.which('mess-abstraction'), 'Argonne abstraction companion is not installed')
    def test_actual_mess_matches_independent_saved_hir_calculation(self):
        previous = Path.cwd()
        with TemporaryDirectory() as directory:
            try:
                par = parameters(directory)
                data = methanol_records()
                os.chdir(directory)
                writer = BimolecularMESS(par, data)
                writer.write_input()
                writer.run()
                self.assertFalse(writer.warnings)
                forward = sorted((rates(Path('me') / (job['stem'] + '.out'))[:, 1]
                                  for job in writer.mess_jobs), key=lambda a: a[0])
                np.testing.assert_allclose(forward, [[5.37149e-18, 8.36044e-16, 2.04624e-13],
                                                    [2.03750e-15, 6.33755e-14, 2.89865e-12]], rtol=2.e-5)
                for job in writer.mess_jobs:
                    self.assertTrue(np.all(rates(Path('me') / (job['stem'] + '.out')) > 0))
            finally:
                os.chdir(previous)

    @unittest.skipUnless(shutil.which('mess-abstraction'), 'Argonne abstraction companion is not installed')
    def test_labelled_replicas_equal_one_unlabelled_contribution(self):
        with TemporaryDirectory() as directory:
            par = parameters(directory)
            writer = BimolecularMESS(par, methanol_records())
            path = next(p for p in writer.models['pathways'] if p['site_counting']['degeneracy'] == 3)
            fragments = ''.join(writer.bimolecular(dict(writer.nodes[p], id=name), writer.reference)
                                for p, name in zip((path['reactant'], path['product']), ('R', 'P')))
            block = writer.rrho(path['stationary_point'], writer.reference, depths=writer._depths(path))
            replica = re.sub(r'(SymmetryFactor\s+)([0-9.eE+-]+)',
                             lambda m: m[1] + str(float(m[2]) * 3.), block, count=1)
            for name, ts in [('one', block), ('three', 'Union\n' + replica * 3 + '\nEnd\n')]:
                text = (f'TemperatureList[K] 300 500 1000\nModelEnergyLimit[kcal/mol] 150\n'
                        f'AtomDistanceMin[angstrom] 0.5\nRateOutput {name}.out\nLogOutput {name}.log\n'
                        + fragments + 'Barrier TS\n' + ts)
                Path(directory, name + '.inp').write_text(text)
                run = subprocess.run([par['mess_abstraction_command'], name + '.inp'], cwd=directory,
                                     capture_output=True, text=True)
                self.assertEqual(run.returncode, 0, run.stderr)
            np.testing.assert_allclose(rates(Path(directory, 'one.out')),
                                       rates(Path(directory, 'three.out')), rtol=1.e-5)

    @unittest.skipUnless(shutil.which('mess'), 'MESS is not installed')
    def test_complex_networks_run_and_disconnected_product_wells_are_separate(self):
        previous = Path.cwd()
        with TemporaryDirectory() as directory:
            try:
                par = parameters(directory)
                par.update(TemperatureList=[500], ModelEnergyLimit=100.)
                for mode in ('rc_only', 'rc_pc', 'pc_only'):
                    root = Path(directory, mode)
                    root.mkdir()
                    records = methanol_records()
                    writer = BimolecularMESS(par, records)
                    paths = records['rate_models']['pathways']
                    if mode != 'pc_only':
                        add_complex(records, 'rc', 'rc', records['reference_node'],
                                    writer.reference - 5 / constants.AUtoKCAL)
                        for path in paths:
                            path['reactant_complex'] = 'rc'
                    if mode != 'rc_only':
                        for i, path in enumerate(paths):
                            zero = sum(writer.ground(d) for d in writer.nodes[path['product']]['fragment_data'])
                            add_complex(records, f'pc{i}', 'pc', path['product'], zero - 2 / constants.AUtoKCAL)
                            path['product_complex'] = f'pc{i}'
                    os.chdir(root)
                    writer = BimolecularMESS(par, records)
                    writer.write_input()
                    self.assertFalse(writer.warnings)
                    self.assertEqual(len(writer.mess_jobs), 2 if mode == 'pc_only' else 1)
                    writer.run()
                    for job in writer.mess_jobs:
                        text = Path('me', job['stem'] + '.inp').read_text()
                        self.assertIn(f'MicroRateOutput {job["stem"]}.micro.out', text)
                        self.assertTrue(Path('me', job['stem'] + '.out').stat().st_size)
                    if mode == 'rc_only':
                        self.assertEqual(text.count('Core PhaseSpaceTheory'), 1)
            finally:
                os.chdir(previous)

    @unittest.skipUnless(shutil.which('mess-abstraction'), 'Argonne abstraction companion is not installed')
    def test_absolute_capture_matches_c6_formula_for_hir_and_mc(self):
        from ase import Atoms
        from scipy.constants import physical_constants, k
        raw = methanol_data()['fragments']
        fragments = [_species_data(saved_point(raw[key])) for key in ('methanol', 'H')]
        with TemporaryDirectory() as directory:
            writer = RecordRenderer(parameters(directory), None)
            for mc in (False, True):
                pair = copy.deepcopy(fragments)
                if mc:
                    t = pair[0]['thermochemistry']
                    first = copy.deepcopy(writer.members(pair[0])[0])
                    second = dict(first, zero_energy_hartree=first['zero_energy_hartree'] + .5/constants.AUtoKCAL)
                    t.update(conformer_representation='MC-RRHO', hir={'rotors': []},
                             conformer_inventory=[dict(r, retained=True, status='valid') for r in (first, second)])
                ref = sum(writer.ground(d) for d in pair)
                text = ('TemperatureList[K] 300 500 1000\nAtomDistanceMin[angstrom] 0.5\n'
                        'RateOutput capture.out\nLogOutput capture.log\n')
                text += ''.join(writer.bimolecular(dict(id=name, fragment_data=pair), ref) for name in ('R', 'P'))
                text += 'Barrier capture\n' + writer.pst(pair, ref, 2)
                Path(directory, 'capture.inp').write_text(text)
                run = subprocess.run([shutil.which('mess-abstraction'), 'capture.inp'], cwd=directory,
                                     capture_output=True, text=True)
                self.assertEqual(run.returncode, 0, run.stderr)
                masses = [Atoms(d['atom']).get_masses().sum() for d in pair]
                mu = masses[0] * masses[1] / sum(masses)
                mu *= physical_constants['atomic mass constant'][0] / physical_constants['electron mass'][0]
                conversion = (100 * physical_constants['Bohr radius'][0])**3 / physical_constants['atomic unit of time'][0]
                for temperature, forward, reverse in rates(Path(directory, 'capture.out')):
                    expected = (6 * math.sqrt(math.pi) * 2**(-1/6) * math.gamma(5/3) * .6**(1/3) *
                                (k * temperature / physical_constants['Hartree energy'][0])**(1/6) / math.sqrt(mu) * conversion)
                    self.assertAlmostEqual(forward / expected, 1., delta=4.e-4)
                    self.assertEqual(forward, reverse)


if __name__ == '__main__':
    unittest.main()
