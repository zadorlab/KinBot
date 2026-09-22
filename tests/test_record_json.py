"""Canonical numeric records must work with the standard strict JSON writer."""
import json
import unittest

import numpy as np

from kinbot.thermochemistry import thermochemistry_evidence, json_record
from tests.test_thermochemistry_evidence import scanned_point


class TestRecordJSON(unittest.TestCase):
    def test_numpy_values_and_nonfinite_observations(self):
        data = json_record({'index': np.int64(3), 'axes': np.array([1, 2]),
                            'nested': {'rank': np.int64(1), 'value': np.float64(np.nan)}})
        self.assertIsInstance(data['index'], int)
        self.assertEqual(data['axes'], [1, 2])
        self.assertIsNone(data['nested']['value'])
        self.assertEqual(data['data_issues'][0]['path'], '$.nested.value')
        json.dumps(data, allow_nan=False)
        species, _ = scanned_point()
        species.freq[1] = np.nan
        species.rotor_projection = {'internal_rank': np.int64(1)}
        data = thermochemistry_evidence(species)
        self.assertIsNone(data['raw_harmonic_frequencies_cm-1'][1])
        self.assertTrue(data['data_issues'])
        json.dumps(data, allow_nan=False)


if __name__ == '__main__':
    unittest.main()
