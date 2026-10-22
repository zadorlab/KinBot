"""Observation/group IDs preserve evidence without assigning reaction flux."""
import copy
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace as NS
import unittest
from unittest.mock import Mock, patch
import numpy as np
from ase.build import molecule
from kinbot import bimol
from kinbot.bimol_network import BimolecularNetwork
from kinbot.stationary_pt import StationaryPoint
from kinbot.conformer_records import ConformerRecord
from kinbot.stereo_routing import StereoRoutingError


class TestBimolecularObservations(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        atoms = molecule('CH4')
        self.point = StationaryPoint('fixture', 0, 1,
            atom=atoms.get_chemical_symbols(), geom=atoms.positions)
        self.point.characterize()
        self.point.energy, self.point.zpe = -3., .01
        self.par = {'multi_conf_tst': 0, 'bimol_energy_tolerance': .25,
                    'bimol_ts_rmsd': .1, 'bimol_basin_rmsd': .1}

    def reaction(self, index, point=None, product=0):
        point = copy.deepcopy(self.point) if point is None else point
        return NS(bimol_search=True, family_name='abstraction', stereopath_id='same-site',
            instance_name=f'candidate_{index}', instance=[0, 1, 2], ts=point,
            rc_basin_id='rc', pc_basin_id='pc', site_degeneracy=1,
            reactant_asymptote_id='reactants', product_asymptote_id=f'products-{product}',
            irc_reactant=self.point, irc_prod=self.point,
            l0_validation={'endpoints': [
                {'state': 'reactant', 'configured_identity': {'status': 'assigned', 'id': 'reactant'}},
                {'state': 'product', 'bonds': [[0, 1+product]],
                 'configured_identity': {'status': 'assigned', 'id': f'product-{product}'}}]})

    def network(self):
        network = object.__new__(BimolecularNetwork)
        network.reference, network.par = 0., self.par
        network.nodes, network.transitions, network.connections = {}, [], []
        network.reference_node = 'rc'
        return network

    def test_distinct_product_groups_keep_their_own_final_irc_fallback(self):
        reactions = [self.reaction(i, product=i//2) for i in range(4)]
        reactions[3].ts.geom[1, 0] += .001
        species = NS(reac_obj=reactions, reac_ts_done=[6]*4)
        bimol.release_final_ts_clusters(species, self.par)
        self.assertEqual(species.reac_ts_done, [7, 8, 7, 8])
        self.assertNotEqual(reactions[0].final_ts_cluster_id, reactions[2].final_ts_cluster_id)
        initial = {r.instance_name: r.final_ts_cluster_id for r in reactions}
        species.reac_ts_done[0], species.reac_ts_done[2] = 2, -999
        bimol.release_final_ts_clusters(species, self.par)
        self.assertEqual(species.reac_ts_done, [2, -998, -999, 7])
        self.assertTrue(reactions[3].final_irc_attempted)
        permuted = [self.reaction(i, product=i//2) for i in range(4)]
        permuted[3].ts.geom[1, 0] += .001
        second = NS(reac_obj=permuted[::-1], reac_ts_done=[6]*4)
        bimol.release_final_ts_clusters(second, self.par)
        self.assertEqual({r.instance_name:r.final_ts_cluster_id for r in permuted}, initial)

    def test_unknown_fragment_identity_preserves_evidence_without_creating_a_node(self):
        network = self.network()
        for status in ('unsupported', 'unavailable'):
            with self.subTest(status=status), patch('kinbot.bimol_network.canonical_identity',
                    return_value={'status': status, 'reason': 'identity fixture'}):
                with self.assertRaisesRegex(StereoRoutingError, 'Observations preserved'):
                    network._add_products([self.point])
        self.assertEqual(network.nodes, {})
        reports = [json.loads(p.read_text()) for p in Path('unsupported_observations').glob('*.json')]
        self.assertEqual(len(reports), 2)
        self.assertTrue(all(not r['rate_model_ready'] for r in reports))
        np.testing.assert_array_equal(reports[0]['observations'][0]['geometry_angstrom'], self.point.geom)

    def test_same_site_selection_does_not_call_conflicting_energies_duplicates(self):
        reactions = [self.reaction(i) for i in range(2)]
        reactions[1].ts.energy += .01
        species = NS(reac_obj=reactions, reac_ts_done=[6, 6])
        bimol.release_final_ts_clusters(species, self.par)
        self.assertEqual([r.final_conformer_role for r in reactions], ['selected', 'archived conformer'])
        species.reac_ts_done[0] = 2
        bimol.release_final_ts_clusters(species, self.par)
        self.assertEqual(reactions[1].final_conformer_archive[1]['role'], 'archived conformer')

    def test_restart_reads_completed_duplicate_endpoints_at_requested_level(self):
        for high_level, complete_endpoint in ((0, False), (1, False), (1, True), (1, 'worker_file')):
            with self.subTest(high_level=high_level, complete_endpoint=complete_endpoint):
                reactions = [self.reaction(i) for i in range(2)]
                species = NS(reac_obj=reactions, reac_ts_done=[6, 6])
                par = dict(self.par, high_level=high_level)
                bimol.release_final_ts_clusters(species, par)
                species.reac_ts_done[0] = -1
                duplicate = reactions[1]
                job = bimol.final_irc_job_name(duplicate)
                suffixes = ('', '_prod', '_prod_high') if high_level and not complete_endpoint else ('', '_prod')
                completed = {f'{job}_IRC_{d}{suffix}' for d in ('F', 'R') for suffix in suffixes}
                row = NS(data=dict(status='normal', energy=-1., zpe=.01, frequencies=[100.],
                                   endpoint_converged=True))
                rows = [row] if complete_endpoint is True else []
                def check_qc(name):
                    if complete_endpoint == 'worker_file' and name.endswith('_prod'):
                        rows[:] = [row]
                    return 'normal' if name in completed else 0
                duplicate.qc = NS(check_qc=Mock(side_effect=check_qc),
                    for_job=lambda name: NS(theory_level='l2', qc='fc' if complete_endpoint == 'worker_file' else 'gauss'),
                    db=NS(select=lambda **kw: rows))
                bimol.release_final_ts_clusters(species, par)
                self.assertEqual(species.reac_ts_done, [-1, 7])
                self.assertEqual(duplicate.final_irc_name, job)
                self.assertTrue(duplicate.final_irc_attempted)
                self.assertEqual(duplicate.final_conformer_role, 'duplicate')
                self.assertEqual({call.args[0] for call in duplicate.qc.check_qc.call_args_list}, completed)
                # Reading the existing path is requested once, not on every poll.
                duplicate.qc.check_qc.reset_mock()
                bimol.release_final_ts_clusters(species, par)
                duplicate.qc.check_qc.assert_not_called()

    def test_restart_does_not_start_missing_failed_or_obsolete_duplicate_jobs(self):
        for case in ('missing_endpoint', 'failed_irc', 'running_endpoint', 'missing_high', 'changed_ts'):
            with self.subTest(case=case):
                reactions = [self.reaction(i) for i in range(2)]
                species = NS(reac_obj=reactions, reac_ts_done=[6, 6])
                par = dict(self.par, high_level=int(case == 'missing_high'))
                bimol.release_final_ts_clusters(species, par)
                species.reac_ts_done[0] = -1
                duplicate = reactions[1]
                job = bimol.final_irc_job_name(duplicate)
                statuses = {f'{job}_IRC_{d}{suffix}': 'normal'
                            for d in ('F', 'R') for suffix in ('', '_prod')}
                if case == 'missing_endpoint':
                    statuses.pop(job + '_IRC_R_prod')
                elif case == 'failed_irc':
                    statuses[job + '_IRC_F'] = 'error'
                elif case == 'running_endpoint':
                    statuses[job + '_IRC_R_prod'] = 'running'
                elif case == 'changed_ts':
                    duplicate.ts.geom[1, 0] += .001
                duplicate.qc = NS(check_qc=lambda name: statuses.get(name, 0),
                    for_job=lambda name: NS(theory_level='l2', qc='gauss'), db=NS(select=lambda **kw: []))
                bimol.release_final_ts_clusters(species, par)
                self.assertEqual(species.reac_ts_done, [-1, -998])

    def test_retained_same_geometry_observations_have_unique_ids_and_export_paths(self):
        reactions = [self.reaction(i) for i in range(2)]
        reactions[1].ts.energy += .01
        network = self.network()
        network._cluster_transitions(reactions)
        self.assertEqual(len({t['id'] for t in network.transitions}), 2)
        self.assertTrue(all(t['energy_conflicts'] and 'unresolved' in t['flux_assignment']
                            for t in network.transitions))
        network.nodes = {key: {'kind': 'well', 'species': self.point} for key in ('rc', 'pc')}
        network.write_unique_ts('unique_ts')
        self.assertEqual(len(json.loads(Path('unique_ts/index.json').read_text())['unique_ts']), 2)
        network.nodes = {}
        network._cluster_complexes([{'fragments': ('same',), 'species': r.ts, 'reaction': r}
                                    for r in reactions], 'rc')
        self.assertEqual(len(network.nodes), 2)
        self.assertTrue(all(node['energy_conflicts'] for node in network.nodes.values()))
        self.assertTrue(all(node['observations'] for node in network.nodes.values()))
        self.assertNotEqual(reactions[0].rc_basin_id, reactions[1].rc_basin_id)
        self.assertEqual({m for node in network.nodes.values() for m in node['members']},
                         {'candidate_0', 'candidate_1'})

    def test_endpoint_conflicts_include_compatible_nonrepresentative_observations(self):
        reactions = [self.reaction(i) for i in range(3)]
        for i, reaction in enumerate(reactions):
            reaction.ts.geom[1, 0] += .2*i
            reaction.ts.energy += i * 1.e-6
        reactions[2].pc_basin_id = 'other_pc'
        reactions[2].product_asymptote_id = 'other_products'
        network = self.network()
        network._cluster_transitions(reactions)
        self.assertEqual(len(network.transitions), 2)
        for transition in network.transitions:
            self.assertEqual(len(transition['endpoint_conflicts']), 1)
            pairs = [set(w['observations']) for e in transition['conflict_evidence'] for w in e['witnesses']]
            self.assertIn({'candidate_1', 'candidate_2'}, pairs)
        json.dumps(network.as_dict(), allow_nan=False)

    def test_archives_keep_nonselected_inventories_and_union_merged_observations(self):
        reactions = [self.reaction(i) for i in range(4)]
        for i, reaction in enumerate(reactions):
            reaction.ts.energy += i * 1.e-6
            reaction.ts.conformer_inventory = (ConformerRecord('local:17', 17,
                f'conf/source_{i}', 'failed', attempted_source_job=f'source_{i}_high'),)
        # Two independent selection groups, later numerically merged by the network.
        for pair in (reactions[:2], reactions[2:]):
            bimol.release_final_ts_clusters(NS(reac_obj=pair, reac_ts_done=[6, 6]), self.par)
        network = self.network()
        network._cluster_transitions([reactions[0], reactions[2]])
        data = network.as_dict()['transitions'][0]
        self.assertEqual(len(data['observations']), 2)
        self.assertEqual(len(data['conformer_archive']), 4)
        records = [r for observation in data['conformer_archive'] for r in observation['conformer_inventory']]
        self.assertEqual({r['source_job'] for r in records}, {f'conf/source_{i}' for i in range(4)})
        self.assertEqual(len({r['observation_member_id'] for r in records}), 4)
        self.assertTrue(all(r['status'] == 'failed' and r['attempted_source_job'] for r in records))
        network.nodes = {key: {'kind': 'well', 'species': self.point} for key in ('rc', 'pc')}
        network.write_unique_ts('unique_ts')
        exported = json.loads(Path('unique_ts/index.json').read_text())['unique_ts'][0]
        self.assertEqual(exported['conformer_archive'], data['conformer_archive'])


if __name__ == '__main__':
    unittest.main()
