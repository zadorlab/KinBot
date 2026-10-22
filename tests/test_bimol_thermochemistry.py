"""Canonical bimolecular artifacts contain the shared HIR evidence records."""
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from kinbot.bimol_network import BimolecularNetwork
from tests.test_thermochemistry_evidence import scanned_point


class TestCanonicalThermochemistry(unittest.TestCase):
    def test_final_ts_and_fragment_evidence_survive_both_canonical_exports(self):
        point, angles = scanned_point()
        point.hir.fourier_fit('test', angles, 0)
        point.optical_isomers = 2
        reaction = SimpleNamespace(
            ts=point, instance_name='candidate', instance=[0, 1, 2],
            irc_reactant=point, irc_prod=point)
        network = BimolecularNetwork.__new__(BimolecularNetwork)
        network.par = {}
        network.reference_node = 'reactants'
        network.nodes = {
            name: {'id': name, 'kind': 'bimolecular', 'energy': 0.,
                   'fragment_species': [point]}
            for name in ('reactants', 'products')}
        network.transitions = [{
            'id': 'ts', 'energy': 10., 'reactant': 'reactants',
            'product': 'products', 'species': point, 'reaction': reaction,
            'members': ['candidate'], 'degeneracy': 3}]
        network.connections = []
        with TemporaryDirectory() as directory:
            previous = Path.cwd()
            try:
                os.chdir(directory)
                network.write_json()
                network.write_unique_ts()
                data = json.loads(Path('bimol_network.json').read_text())
                index = json.loads(Path('bimol/unique_ts/index.json').read_text())
                metadata = json.loads(Path('bimol/unique_ts/ts/metadata.json').read_text())
                self.assertEqual(metadata['irc_level'], 'L1')
                network.par['high_level'] = 1
                network.write_unique_ts('bimol/unique_ts_l2')
                l2_metadata = json.loads(Path('bimol/unique_ts_l2/ts/metadata.json').read_text())
                l2_index = json.loads(Path('bimol/unique_ts_l2/index.json').read_text())
                self.assertEqual(l2_metadata['level'], 'L2')
                self.assertEqual(l2_metadata['irc_level'], 'L2')
                self.assertEqual(l2_index['unique_ts'][0]['irc_level'], 'L2')
            finally:
                os.chdir(previous)
        canonical = data['transitions'][0]['stationary_point']['thermochemistry']
        self.assertEqual(canonical, metadata['thermochemistry'])
        self.assertEqual(metadata, index['unique_ts'][0])
        self.assertEqual(data['nodes'][0]['fragment_data'][0]['thermochemistry'], canonical)
        self.assertEqual(canonical['hir']['rotors'][0]['dihedral'], point.dihed[0])
        self.assertEqual(len(canonical['hir']['rotors'][0]['points']), 12)
        self.assertEqual(data['transitions'][0]['degeneracy'], 3)
        self.assertIsNone(canonical['optical_counting']['remaining_multiplier'])
        self.assertFalse(data['rate_model_ready'])
        self.assertEqual(data['schema_version'], 3)
        self.assertIn('raw harmonic', metadata['frequencies_semantics'])
        self.assertEqual(len(metadata['frequencies_cm-1']), 12)
        self.assertEqual(len(data['transitions'][0]['stationary_point']['frequencies']), 11)


if __name__ == '__main__':
    unittest.main()
