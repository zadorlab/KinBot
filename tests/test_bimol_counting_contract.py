"""Final canonical writers consume complete counting records without job files."""
import copy
import json
import os
import pickle
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

from kinbot.bimol_network import BimolecularNetwork
from kinbot.optimize import Optimize
from kinbot.reactions.reac_abstraction import Abstraction
from kinbot.stationary_pt import StationaryPoint
from kinbot.thermochemistry import thermochemistry_evidence
from tests.counting_fixtures import methanol_data, saved_point


def saved_network():
    fixture = methanol_data()
    fragments = {key: saved_point(raw) for key, raw in fixture['fragments'].items()}
    a, b = fragments['H'], fragments['methanol']
    parent = StationaryPoint('reactants', 0, 2, atom=list(a.atom) + list(b.atom),
                             geom=np.vstack((a.geom + [10., 0., 0.], b.geom)), fragA=a, fragB=b)
    parent.characterize()
    parent.energy, parent.zpe = a.energy + b.energy, a.zpe + b.zpe
    par = {'bimol_energy_tolerance': .25, 'bimol_ts_rmsd': .1, 'bimol_basin_rmsd': .1}
    for data in fixture['saddles']:
        reaction = Abstraction(parent, None, par, np.array(data['instance'], dtype=np.int64), data['source_job'])
        reaction.ts = saved_point(data)
        for key in ('site_members', 'site_relations', 'stereopath_id'):
            setattr(reaction, key, copy.deepcopy(data[key]))
        reaction.site_degeneracy = np.int64(data['site_degeneracy'])
        reaction.ts.rotor_projection = getattr(reaction.ts, 'rotor_projection', None)
        if reaction.ts.rotor_projection:
            reaction.ts.rotor_projection['internal_rank'] = np.int64(1)
        endpoints = []
        for raw in data['endpoints']:
            point = StationaryPoint('endpoint', 0, 2, atom=raw['atoms'], geom=np.array(raw['geometry_angstrom']))
            point.characterize()
            endpoints.append(point)
        reaction.irc_reactant, reaction.irc_prod = endpoints
        for key in ('H2', 'CH2OH' if data['instance'][1] == 3 else 'CH3O'):
            optimization = Optimize.__new__(Optimize)
            optimization.species = fragments[key]
            reaction.prod_opt.append(optimization)
        reaction.actual_bond_changes = {'formed': [reaction.instance[:2]], 'broken': [reaction.instance[1:]]}
        archive = thermochemistry_evidence(saved_point(methanol_data('MeOH_H')['saddles'][0]))
        reaction.final_conformer_archive = [{'observation': 'previous-observation', 'index': np.int64(7),
                                            'thermochemistry': archive, 'independent_rate_contribution': False}]
        parent.reac_obj.append(reaction)
        parent.reac_ts_done.append(-1)
    return BimolecularNetwork(parent, par)


class TestBimolecularCountingContract(unittest.TestCase):
    def test_saved_numeric_types_and_all_canonical_records_round_trip(self):
        network = saved_network()
        before = pickle.dumps(network.species)
        expected = network.as_dict()
        archive = copy.deepcopy(expected['transitions'][0]['conformer_archive'])
        with TemporaryDirectory() as directory:
            previous = Path.cwd()
            try:
                os.chdir(directory)
                # No worker scripts, database or IRC trajectory files exist here.
                with patch('sqlite3.connect', side_effect=AssertionError('No DB reconstruction')):
                    network.write_unique_ts()
                    network.write_json()
                    first = Path('bimol_network.json').read_bytes()
                    network.write_json()
                    network.write_unique_ts()
                    self.assertEqual(Path('bimol_network.json').read_bytes(), first)
                data = json.loads(first)
                index = json.loads(Path('bimol/unique_ts/index.json').read_text())
                self.assertEqual(index['schema_version'], 3)
                for transition in data['transitions']:
                    meta = json.loads((Path('bimol/unique_ts') / transition['id'] / 'metadata.json').read_text())
                    self.assertEqual(meta['schema_version'], 3)
                    self.assertEqual(meta, next(m for m in index['unique_ts'] if m['id'] == transition['id']))
                    self.assertEqual(meta['thermochemistry'], transition['stationary_point']['thermochemistry'])
                    self.assertEqual(meta['site_counting'], transition['site_counting'])
                    self.assertEqual(meta['site_counting']['unlabelled_model_conversion']['additional_site_multiplier'], 1)
                    for label in ('reactant', 'product'):
                        node = next(n for n in data['nodes'] if n['id'] == transition[label])
                        self.assertEqual(meta['endpoints'][label]['fragment_data'], node['fragment_data'])
                    self.assertEqual(meta['thermochemistry']['optical_counting']['remaining_multiplier'], 1)
                    self.assertFalse(meta['thermochemistry']['rate_model_ready'])
                self.assertEqual(data, expected)
                self.assertEqual(data['transitions'][0]['conformer_archive'], archive)
            finally:
                os.chdir(previous)
        self.assertEqual(pickle.dumps(network.species), before)

    def test_conflicts_do_not_turn_the_site_class_into_assigned_flux(self):
        network = saved_network()
        for transition in network.transitions:
            transition['endpoint_conflicts'] = ['incompatible-observation']
        data = network.as_dict()
        for transition in data['transitions']:
            self.assertIsNone(transition['site_counting']['unlabelled_model_conversion']['additional_site_multiplier'])
            self.assertEqual(transition['site_counting']['physical_path_completeness'], 'not_established')
        self.assertFalse(data['rate_model_ready'])

    def test_equivalent_site_numbers_do_not_change_counting(self):
        network = saved_network()
        transition = network.transitions[0]
        reaction = transition['reaction']
        reaction.site_members = {'abstractor': [0], 'transferred_atom': [3]}
        reaction.site_relations = {'abstractor': 'homotopic', 'transferred_atom': 'homotopic'}
        reaction.site_degeneracy = transition['degeneracy'] = 1
        observation = copy.copy(reaction)
        observation.site_members = {'abstractor': [0], 'transferred_atom': [4]}
        transition['_observations'] = [reaction, observation]

        expected = network._site_counting(transition)
        self.assertEqual(expected['unlabelled_model_conversion']['additional_site_multiplier'], 1.)
        self.assertEqual(expected['site_members'], reaction.site_members)
        with TemporaryDirectory() as directory:
            network.write_unique_ts(Path(directory) / 'unique')
            meta = json.loads((Path(directory) / 'unique' / transition['id'] / 'metadata.json').read_text())
        self.assertEqual(meta['site_counting'], expected)
        self.assertEqual(next(t for t in network.as_dict()['transitions']
                              if t['id'] == transition['id'])['site_counting'], expected)

        # Different or malformed classes must still be reported as unresolved.
        for change in (
                {'stereopath_id': 'different-pathway'},
                {'site_relations': {'abstractor': 'homotopic', 'transferred_atom': 'diastereotopic'}},
                {'site_degeneracy': 2},
                {'site_members': {'abstractor': [0], 'transferred_atom': [4, 4]}},
                {'site_members': {'abstractor': [0], 'transferred_atom': [-1]}}):
            with self.subTest(change=change):
                conflicting = copy.copy(observation)
                conflicting.__dict__.update(change)
                transition['_observations'] = [reaction, conflicting]
                self.assertIsNone(network._site_counting(transition)[
                    'unlabelled_model_conversion']['additional_site_multiplier'])

    def test_complex_endpoint_keeps_its_own_record(self):
        network = saved_network()
        selected = network.transitions[0]
        # Exercise the complex-record writer without revising RC classification.
        point = saved_point(methanol_data()['fragments']['methanol'])
        network.nodes['complex-fixture'] = {'id': 'complex-fixture', 'kind': 'complex',
                                            'species': point, 'energy': -1.}
        selected['reactant'] = 'complex-fixture'
        with TemporaryDirectory() as directory:
            network.write_unique_ts(Path(directory) / 'unique')
            meta = json.loads((Path(directory) / 'unique' / selected['id'] / 'metadata.json').read_text())
        node = next(n for n in network.as_dict()['nodes'] if n['id'] == 'complex-fixture')
        self.assertEqual(meta['endpoints']['reactant']['stationary_point'], node['stationary_point'])
        self.assertNotEqual(node['stationary_point']['thermochemistry']['source_job'], meta['selected_ts_job'])


if __name__ == '__main__':
    unittest.main()
