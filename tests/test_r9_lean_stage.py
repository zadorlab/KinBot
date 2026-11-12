"""Verify the selected R9 drafts against native atom order and discovery."""
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import tempfile
import unittest

from kinbot.parameters import Parameters
from kinbot.reaction_finder import ReactionFinder
from kinbot.stationary_pt import StationaryPoint
from kinbot.vrc_multistate import validate_surface_model


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(path.read_text())


@contextmanager
def temporary_directory():
    previous = Path.cwd()
    with tempfile.TemporaryDirectory() as directory:
        try:
            os.chdir(directory)
            yield
        finally:
            os.chdir(previous)


class NoCalculation:
    def __getattr__(self, key):
        raise AssertionError(f'Review-stage discovery attempted QC access through {key}')


class TestLeanR9Stage(unittest.TestCase):
    def test_exact_scope_and_normal_discovery(self):
        expected_counts = {'ch3s': (2, 2, 4, 2), 'ch3ss': (5, 4, 9, 5)}
        for case, counts in expected_counts.items():
            folder = ROOT / 'examples/anl' / (case + '_anl1_f12') / 'production_r9'
            scope = read(folder / 'scope.json')
            parameters = Parameters(str(folder / 'input.draft.json')).par
            populations = list(scope['populations'].values())
            observed = (
                sum(p['kind'] == 'Well' for p in populations),
                sum(p['kind'] == 'Bimolecular' for p in populations),
                len(scope['connections']),
                sum(c['kind'] == 'vrc' for c in scope['connections']),
            )
            self.assertEqual(observed, counts)
            self.assertEqual(set(parameters['keep_chemids']), set(scope['discovery_workers']))
            self.assertEqual(hashlib.sha256((folder / 'lean_scout.inp').read_bytes()).hexdigest(),
                             scope['source_mess_sha256'])
            self.assertEqual(parameters['queue_job_limit'], 20)
            self.assertEqual(parameters['queue_exclusive'], 1)
            self.assertEqual(parameters['ppn'], 8)
            self.assertEqual(parameters['l1_pes_energy_ceiling'], 150.)
            self.assertTrue(parameters['l1_vdw_depth_check'])
            self.assertEqual(parameters['profiled_l1_reuse_only'], 1)
            self.assertEqual(parameters['composite_method'], 'ANL1-F12')
            self.assertTrue(parameters['electronic_state_models'])
            self.assertEqual(parameters['run_me'], 0)
            for geometry in read(folder / 'selection_geometries.json'):
                worker = geometry['worker']
                with self.subTest(case=case, worker=worker), temporary_directory():
                    species = StationaryPoint(worker, geometry['charge'],
                                              geometry['multiplicity'],
                                              structure=geometry['structure'])
                    species.characterize()
                    species.name = worker
                    finder = ReactionFinder(species, parameters, NoCalculation())
                    finder.find_reactions()
                    expected = {c['reaction'] for c in scope['connections'] if c['owner'] == worker}
                    companions = {name.replace('_hom_sci_', '_dissociation_saddle_')
                                  for name in expected if '_hom_sci_' in name}
                    self.assertEqual(set(species.reac_name), expected | companions)

    def test_model_readiness_is_explicit_for_every_selected_vrc(self):
        status_counts = {}
        for case in ('ch3s', 'ch3ss'):
            folder = ROOT / 'examples/anl' / (case + '_anl1_f12') / 'production_r9'
            data = read(folder / 'input.draft.json')
            scope = read(folder / 'scope.json')
            readiness = read(folder / 'readiness.json')
            self.assertFalse(readiness['submission_ready'])
            selected = {c['reaction'] for c in scope['connections'] if c['kind'] == 'vrc'}
            self.assertEqual(selected, {c['reaction'] for c in readiness['channels']})
            models = [data['vrc_tst_surface_model'], *data.get('vrc_tst_surface_models', {}).values()]
            installed = set()
            for model in models:
                validate_surface_model(model)
                installed.update(c['reaction'] for c in model['channels'])
            for channel in readiness['channels']:
                status = channel['status']
                status_counts[status] = status_counts.get(status, 0) + 1
                self.assertEqual(channel['reaction'] in installed,
                                 status in ('ready', 'installed_selection_supported_seed_binding_pending'))
            requests = read(folder / 'new_h_loss_reviews.json')
            self.assertEqual(len(requests), 1)
            self.assertIsNone(requests[0]['selected_cas'])
            self.assertIn(requests[0]['reaction'], selected)
        self.assertEqual(status_counts, {
            'ready': 2,
            'new_h_loss_review': 2,
            'installed_selection_supported_seed_binding_pending': 1,
            'selection_supported_installation_pending': 2,
        })


if __name__ == '__main__':
    unittest.main()
