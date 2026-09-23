"""Linearity decision for the external-mode projection and for the MESS geometry.

A linear molecule left slightly bent by an optimiser must keep 3N-5 vibrations;
a genuinely bent molecule, a clearly bent one, and a prolate molecule with
off-axis atoms must keep 3N-6. Species judged linear are written to MESS as
exactly linear, with a comment; stored geometries are never modified. A
quasi-linear species that KinBot judges non-linear but MESS would call linear
gets explicit rotational constants.
"""

import unittest

from ase import Atoms
from ase.build import molecule
import numpy as np

from kinbot import frequencies
from kinbot.mess import MESS
from kinbot.stationary_pt import StationaryPoint

A2B = 1.0 / 0.529177210903
R0 = 1.16 * A2B


# --- a rotation/translation-invariant model potential for CO2 (Hartree, Bohr) ----
def model_energy(x, ks, kb, d0):
    """Two harmonic C-O stretches plus a harmonic bend on the carbon's distance
    from the O...O line; d0 = 0 gives a linear minimum, d0 > 0 a bent one."""
    o1, c, o2 = x[0:3], x[3:6], x[6:9]
    energy = 0.5 * ks * ((np.linalg.norm(c - o1) - R0) ** 2 + (np.linalg.norm(c - o2) - R0) ** 2)
    u = o2 - o1
    u = u / np.linalg.norm(u)
    v = c - o1
    d = np.linalg.norm(v - np.dot(v, u) * u)
    return energy + 0.5 * kb * (d - d0) ** 2


def model_hessian(geom_angstrom, d0=0., ks=0.573, kb=0.147, h=1e-4):
    """Numerical Cartesian Hessian (Hartree/Bohr^2); kb reproduces the 667 cm-1
    CO2 bend at the linear geometry."""
    x = np.asarray(geom_angstrom, dtype=float).ravel() * A2B
    n = len(x)
    hess = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            pts = []
            for si, sj in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                y = x.copy()
                y[i] += si * h
                y[j] += sj * h
                pts.append(model_energy(y, ks, kb, d0))
            hess[i, j] = (pts[0] - pts[1] - pts[2] + pts[3]) / (4 * h * h)
    return 0.5 * (hess + hess.T)


def co2(angle):
    """CO2 with the given O-C-O angle, carbon at the origin."""
    half = np.radians((180. - angle) / 2.)
    return np.array([[-1.16 * np.cos(half), 1.16 * np.sin(half), 0.],
                     [0., 0., 0.],
                     [1.16 * np.cos(half), 1.16 * np.sin(half), 0.]])


def bent_minimum_offset(angle):
    return 1.16 * np.sin(np.radians((180. - angle) / 2.)) * A2B


def species_from(symbols, geom, name='probe'):
    point = StationaryPoint(name, 0, 1, atom=list(symbols), geom=np.asarray(geom, dtype=float))
    point.characterize()
    return point


def count(point, hessian=None):
    natom = point.natom
    if hessian is None:
        hessian = np.eye(3 * natom)
    freqs, _ = frequencies.get_frequencies(point, hessian, point.geom)
    return len(freqs), freqs


def methyl_cyanodiyne():
    x = [0.0, 1.46, 2.67, 4.05, 5.26, 6.64, 7.80]
    symbols = ['C', 'C', 'C', 'C', 'C', 'C', 'N']
    positions = [[xi, 0., 0.] for xi in x]
    for k in range(3):
        angle = 2 * np.pi * k / 3
        positions.append([-0.36, 1.03 * np.cos(angle), 1.03 * np.sin(angle)])
        symbols.append('H')
    return Atoms(symbols, positions=positions)


def writer():
    mess = MESS.__new__(MESS)
    mess.par = {}
    return mess


class TestProjectionCount(unittest.TestCase):
    def test_exactly_linear_molecules_keep_3n_minus_5(self):
        for name in ('CO2', 'HCN', 'C2H2', 'N2'):
            with self.subTest(molecule=name):
                atoms = molecule(name)
                point = species_from(atoms.get_chemical_symbols(), atoms.positions)
                self.assertEqual(count(point)[0], 3 * len(atoms) - 5)

    def test_linear_minimum_left_bent_by_optimiser_keeps_the_bend(self):
        """Case A: the potential's minimum is linear, the geometry is a residual."""
        for angle in (179.9, 179.5, 179., 178.5):
            with self.subTest(oco_angle=angle):
                geom = co2(angle)
                point = species_from('OCO', geom)
                n, freqs = count(point, model_hessian(geom, d0=0.))
                self.assertEqual(n, 4)
                # both bend components present at ~667 cm-1, nothing near zero
                self.assertGreater(min(freqs), 600.)

    def test_genuinely_bent_minimum_is_not_linear(self):
        """Case B: the same angles, but the potential's minimum is at that angle."""
        for angle in (179.9, 179.5, 179.):
            with self.subTest(oco_angle=angle):
                geom = co2(angle)
                point = species_from('OCO', geom)
                n, freqs = count(point, model_hessian(geom, d0=bent_minimum_offset(angle)))
                self.assertEqual(n, 3)
                self.assertGreater(min(freqs), 600.)   # the ~0 cm-1 rotation is gone

    def test_angle_gate_rejects_clearly_bent_geometries(self):
        for angle in (175., 170., 160.):
            with self.subTest(oco_angle=angle):
                geom = co2(angle)
                point = species_from('OCO', geom)
                # even with a linear-minimum Hessian: 175 deg is not a residual
                self.assertEqual(count(point, model_hessian(geom, d0=0.))[0], 3)

    def test_prolate_molecule_with_off_axis_hydrogens_is_not_linear(self):
        atoms = methyl_cyanodiyne()
        point = species_from(atoms.get_chemical_symbols(), atoms.positions)
        self.assertEqual(count(point)[0], 3 * len(atoms) - 6)

    def test_genuinely_bent_molecules_keep_3n_minus_6(self):
        for name in ('H2O', 'NH3', 'CH3OH', 'C2H6'):
            with self.subTest(molecule=name):
                atoms = molecule(name)
                point = species_from(atoms.get_chemical_symbols(), atoms.positions)
                self.assertEqual(count(point)[0], 3 * len(atoms) - 6)

    def test_decision_is_logged(self):
        geom = co2(179.)
        point = species_from('OCO', geom, name='co2')
        with self.assertLogs('KinBot', level='INFO') as logs:
            count(point, model_hessian(geom, d0=0.))
        self.assertIn('treated as a linear rotor', '\n'.join(logs.output))

    def test_frequencies_are_real_for_degenerate_spectra(self):
        atoms = molecule('CO2')
        point = species_from(atoms.get_chemical_symbols(), atoms.positions)
        rng = np.random.default_rng(7)
        block = rng.normal(size=(9, 9))
        n, freqs = count(point, block @ block.T)
        self.assertEqual(n, 4)
        self.assertTrue(all(isinstance(f, float) for f in freqs), freqs)


class TestLinearize(unittest.TestCase):
    def test_returns_collinear_copy_and_leaves_input_alone(self):
        geom = co2(179.)
        before = geom.copy()
        linear = frequencies.linearize(geom, ['O', 'C', 'O'])
        np.testing.assert_array_equal(geom, before)
        self.assertEqual(frequencies.max_bend_deviation(linear, species_from('OCO', geom).bond), 0.)
        # bond lengths change only in second order
        for i in (0, 2):
            self.assertLess(abs(np.linalg.norm(linear[i] - linear[1]) - np.linalg.norm(geom[i] - geom[1])), 1e-3)


class TestMessGeometry(unittest.TestCase):
    def test_linear_species_is_written_exactly_linear_with_a_note(self):
        point = species_from('OCO', co2(179.), name='co2')
        point.reduced_freqs = [667., 667., 1333., 2349.]           # 3N-5
        block = writer().rotor_geom(point)
        lines = [l for l in block.splitlines() if not l.strip().startswith('!')]
        coords = np.array([[float(v) for v in l.split()[1:]] for l in lines])
        self.assertEqual(frequencies.max_bend_deviation(coords, point.bond), 0.)
        self.assertIn('! geometry linearised for the rigid-rotor model', block)
        np.testing.assert_allclose(point.geom, co2(179.))        # stored geometry untouched

    def test_nonlinear_species_is_written_as_calculated_without_note(self):
        point = species_from('OCO', co2(179.), name='co2')
        point.reduced_freqs = [667., 1333., 2349.]                 # 3N-6
        block = writer().rotor_geom(point)
        self.assertNotIn('!', block)
        self.assertEqual(block, writer().make_geom(point.geom, point.atom))

    def test_inconsistent_count_and_geometry_gets_a_warning_not_a_snap(self):
        point = species_from('OCO', co2(170.), name='co2')
        point.reduced_freqs = [667., 667., 1333., 2349.]
        with self.assertLogs('KinBot', level='WARNING'):
            block = writer().rotor_geom(point)
        self.assertIn('! WARNING: 3N-5 frequencies but the geometry deviates', block)
        self.assertIn(writer().make_geom(point.geom, point.atom), block)

    def test_diatomic_and_ordinary_species_untouched(self):
        atoms = molecule('H2O')
        point = species_from(atoms.get_chemical_symbols(), atoms.positions, name='h2o')
        point.reduced_freqs = [1595., 3657., 3756.]
        self.assertEqual(writer().rotor_geom(point), writer().make_geom(point.geom, point.atom))


class TestReverseMismatch(unittest.TestCase):
    """KinBot non-linear (3N-6) but MESS's I_min/I_mid < 1e-5 test would say linear."""

    def test_quasi_linear_species_gets_explicit_rotational_constants(self):
        geom = co2(179.9)                                   # I_min/I_mid ~ 2e-7
        point = species_from('OCO', geom, name='co2')
        point.reduced_freqs = [667., 1333., 2349.]          # 3N-6: a bent minimum
        self.assertLess(frequencies.moment_ratio(geom, point.atom), frequencies.MESS_LINEAR_MOMENT_RATIO)
        with self.assertLogs('KinBot', level='WARNING'):
            line = writer().rotor_core_line(point)
        self.assertIn('RotationalConstants[1/cm]', line)
        values = [float(v) for v in line.splitlines()[0].split()[1:]]
        self.assertEqual(len(values), 3)
        np.testing.assert_allclose(values, frequencies.rotational_constants(geom, point.atom), rtol=1e-5)
        self.assertGreater(values[0], 1e4)                  # the near-axis rotation is frozen out
        self.assertIn('! quasi-linear species', line)

    def test_ordinary_species_get_no_extra_line(self):
        cases = [(('O', 'C', 'O'), co2(179.), [667., 1333., 2349.]),          # bent, MESS also says non-linear
                 (('O', 'C', 'O'), co2(179.9), [667., 667., 1333., 2349.]),   # linear for both
                 (molecule('H2O').get_chemical_symbols(), molecule('H2O').positions, [1595., 3657., 3756.])]
        for symbols, geom, freqs in cases:
            point = species_from(symbols, geom, name='x')
            point.reduced_freqs = freqs
            self.assertEqual(writer().rotor_core_line(point), '')

    def test_rotational_constants_of_linear_co2(self):
        # exactly linear CO2: one vanishing moment, two equal ones, B = h/(8 pi^2 c) / (2 m_O r^2)
        atoms = molecule('CO2')
        symbols = atoms.get_chemical_symbols()
        r = np.linalg.norm(atoms.positions[1] - atoms.positions[0])
        expected = frequencies.ROTATIONAL_CONSTANT_FACTOR / (2 * 15.994915 * r ** 2)
        b = frequencies.rotational_constants(atoms.positions, symbols)
        self.assertAlmostEqual(b[1], expected, places=4)
        self.assertAlmostEqual(b[2], expected, places=4)
        self.assertAlmostEqual(b[1], 0.379, places=2)   # 0.390 cm-1 experimentally, at r = 1.162 A


class TestMCLinearityIntegration(unittest.TestCase):
    def test_member_renderer_uses_member_geometry_and_frequency_count(self):
        from kinbot.conformer_records import ConformerRecord
        renderer = MESS({'multi_conf_tst': 1, 'freq_uq_ref': 1000., 'freq_uq_max_exp': 1.}, None)
        point = species_from('OCO', co2(170.), name='co2')
        point.mult = 1
        for angle, modes, expected in [
                (179., [667., 667., 1333., 2349.], 'geometry linearised'),
                (179.9, [667., 1333., 2349.], 'RotationalConstants[1/cm]')]:
            record = ConformerRecord('c1', 1, 'c1', 'valid',
                geometry=tuple(map(tuple, co2(angle))),
                frequencies_cm1=tuple(modes), sigma_ext=2., remaining_optical_weight=1.)
            block = renderer._member_rrho(point, record, 1., 0.)
            self.assertIn(expected, block)
            np.testing.assert_allclose(point.geom, co2(170.))


if __name__ == '__main__':
    unittest.main()
