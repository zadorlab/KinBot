"""Values pasted into generated job scripts must be builtin Python types (#84)."""

import pathlib
import re
import unittest

import numpy as np

from kinbot.utils import plain_symbols, plain_geometry

KINBOT = pathlib.Path(__file__).resolve().parent.parent / 'kinbot'


class TestPlainTemplateValues(unittest.TestCase):
    def test_symbols_from_numpy_array_are_builtin_str(self):
        symbols = plain_symbols(np.array(['F', 'C', 'O']))
        self.assertEqual(symbols, ['F', 'C', 'O'])
        self.assertTrue(all(type(s) is str for s in symbols))
        # repr must be usable in a script that never imports numpy
        self.assertNotIn('np.', repr(symbols))
        self.assertEqual(eval(repr(symbols)), symbols)

    def test_geometry_from_numpy_array_is_builtin_floats(self):
        geom = plain_geometry(np.array([[1.0606767, -0.68427962, 0.23015213],
                                        [0.02739081, 0.06808482, -0.12572072]]))
        self.assertTrue(all(type(x) is float for row in geom for x in row))
        self.assertNotIn('np.', repr(geom))
        np.testing.assert_allclose(eval(repr(geom)), geom)

    def test_geometry_accepts_lists_of_numpy_rows(self):
        rows = [np.array([0., 0., 0.]), np.array([1.1, 0., 0.])]
        self.assertEqual(plain_geometry(rows), [[0., 0., 0.], [1.1, 0., 0.]])

    def test_no_template_site_pastes_raw_numpy_values(self):
        """Guard against reintroducing list(species.atom) / list([list(gi) ...])."""
        offenders = []
        for name in ('qc.py', 'irc.py', 'reac_family.py'):
            text = (KINBOT / name).read_text()
            for match in re.finditer(r'atom=list\(|=list\(\[list\(gi\)', text):
                line = text.count('\n', 0, match.start()) + 1
                offenders.append(f'{name}:{line}')
        self.assertEqual(offenders, [])


if __name__ == '__main__':
    unittest.main()
