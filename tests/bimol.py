"""Unit tests for bimolecular geometry and force helpers."""

import copy
import ast
import unittest
import os
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.calculators.calculator import Calculator, all_changes

from kinbot import bimol, constants, frequencies, symmetry, zmatrix
from kinbot.bimol_network import BimolecularNetwork, symmetry_rmsd
from kinbot.irc import BimolecularIRCResult, IRC
from kinbot.optimize import Optimize
from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.reaction_finder_bimol import ReactionFinderBimol
from kinbot.reaction_generator import ReactionGenerator
from kinbot.reactions.reac_abstraction import (
    Abstraction, _roles, generate_encounters)
from kinbot.stationary_pt import StationaryPoint
from kinbot.species_routing import routing_name
from kinbot.stereochemistry import (assign_optical_metadata,
                                     conformer_symmetry, optical_isomers,
                                     refine_equivalence_group,
                                     stereotopic_hydrogen_equivalence)


def compatible_l0_validation():
    # Fixed known endpoint evidence for tests isolating selection policy.
    return {'endpoints': [
        {'state': state, 'configured_identity': {'status': 'assigned', 'id': state}}
        for state in ('reactant', 'product')]}


class ZeroCalculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def calculate(self, atoms=None, properties=('energy', 'forces'),
                  system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.results = {'energy': 0., 'forces': np.zeros((len(atoms), 3))}


def stationary(atoms, mult):
    species = StationaryPoint('test', 0, mult,
                              atom=atoms.get_chemical_symbols(),
                              geom=atoms.positions.copy())
    species.characterize()
    return species


def combined(first, second, mult=2):
    atoms = list(first.atom) + list(second.atom)
    geom = np.vstack((first.geom, second.geom + np.array([6., 0., 0.])))
    species = StationaryPoint('bimol', 0, mult, atom=atoms, geom=geom,
                              fragA=first, fragB=second)
    species.characterize()
    return species


def default_parameters():
    with TemporaryDirectory() as directory:
        path = Path(directory, 'input.json')
        path.write_text(json.dumps({'barrier_threshold': 50.}))
        return Parameters(path).par


class TestBimol(unittest.TestCase):
    def setUp(self):
        # A theory change in one test must not restart another test's jobs.
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(directory.name)

    @classmethod
    def setUpClass(cls):
        cls.methane = stationary(molecule('CH4'), 1)
        cls.hydroxyl = stationary(molecule('OH'), 2)

    def test_fragment_roles_are_order_independent(self):
        forward = combined(self.methane, self.hydroxyl)
        reverse = combined(self.hydroxyl, self.methane)
        self.assertEqual(_roles(forward, [5, 1, 0]),
                         ((1, 0), (0, 1), (0, 0)))
        self.assertEqual(_roles(reverse, [0, 3, 2]),
                         ((0, 0), (1, 1), (1, 0)))

    def test_encounter_preserves_canonical_order(self):
        settings = bimol.DEFAULTS.copy()
        forward = combined(self.methane, self.hydroxyl)
        reverse = combined(self.hydroxyl, self.methane)
        first = next(generate_encounters(
            forward, self.methane.geom, self.hydroxyl.geom,
            [5, 1, 0], settings))[0]
        second = next(generate_encounters(
            reverse, self.hydroxyl.geom, self.methane.geom,
            [0, 3, 2], settings))[0]
        self.assertEqual(first.get_chemical_symbols(), forward.atom)
        self.assertEqual(second.get_chemical_symbols(), reverse.atom)
        np.testing.assert_allclose(first.positions[:5], self.methane.geom)
        np.testing.assert_allclose(second.positions[2:], self.methane.geom)

    def test_open_direction_environments(self):
        atomic = stationary(Atoms('Cl', positions=[[0., 0., 0.]]), 2)
        self.assertEqual(len(bimol.open_directions(
            atomic, atomic.geom, 0)), 1)
        oh_direction = bimol.open_directions(
            self.hydroxyl, self.hydroxyl.geom, 0)[0]
        oh_bond = bimol.unit_vector(
            self.hydroxyl.geom[1] - self.hydroxyl.geom[0])
        angle = np.degrees(np.arccos(np.dot(oh_direction, oh_bond)))
        self.assertAlmostEqual(angle, 104.5, places=6)

        methyl = stationary(molecule('CH3'), 2)
        directions = bimol.open_directions(methyl, methyl.geom, 0)
        self.assertEqual(len(directions), 2)
        np.testing.assert_allclose(directions[0], -directions[1])
        for neighbor in np.flatnonzero(methyl.bond01[0]):
            bond = bimol.unit_vector(methyl.geom[neighbor] - methyl.geom[0])
            self.assertAlmostEqual(np.dot(directions[0], bond), 0., places=6)

    def test_approach_schedule(self):
        settings = bimol.DEFAULTS.copy()
        directions = list(bimol.approach_directions(
            np.array([0., 0., 1.]), settings))
        self.assertEqual(len(directions), 13)
        angles = [round(np.degrees(np.arccos(direction[2])))
                  for direction in directions]
        self.assertEqual(angles.count(0), 1)
        self.assertEqual(angles.count(20), 4)
        self.assertEqual(angles.count(40), 4)
        self.assertEqual(angles.count(60), 4)

    def test_linear_bias_force_signs(self):
        atoms = Atoms(['C', 'H', 'O'], positions=[[0., 0., 0.],
                                                  [1., 0., 0.],
                                                  [3., 0., 0.]])
        atoms.calc = ZeroCalculator()
        settings = {'break_scale': 1., 'form_scale': 1.}
        atoms.set_constraint(bimol._bias_constraints(
            1., [[2, 1]], [[1, 0]], settings))
        forces = atoms.get_forces()
        np.testing.assert_allclose(forces, [[-1., 0., 0.],
                                           [2., 0., 0.],
                                           [-1., 0., 0.]])
        self.assertAlmostEqual(atoms.get_potential_energy(), 1.)
        self.assertAlmostEqual(
            atoms.get_potential_energy(apply_constraint=False), 0.)

    def test_force_path_frames_do_not_retain_bias(self):
        atoms = Atoms(['O', 'H', 'C'], positions=[[3., 0., 0.],
                                                  [1., 0., 0.],
                                                  [0., 0., 0.]])
        settings = bimol.DEFAULTS.copy()
        settings['path_steps'] = 1
        with TemporaryDirectory() as directory:
            trial = bimol._trial(
                atoms, .1, [[0, 1]], [[1, 2]], ZeroCalculator(),
                settings, Path(directory))
        self.assertTrue(trial['frames'])
        self.assertTrue(all(not frame.constraints
                            for frame in trial['frames']))

    def test_force_path_starts_at_first_nonzero_force(self):
        species = combined(self.hydroxyl, self.methane)
        encounter = Atoms(species.atom, positions=species.geom)
        settings = bimol.DEFAULTS.copy()
        settings.update({'force_max': .2, 'force_levels': 3})
        forces = []
        directories = []

        def trial(start, force, formed, broken, calculator, settings,
                  directory):
            forces.append(force)
            directories.append(directory.name)
            return {'force': force, 'atoms': start.copy(), 'frames': [],
                    'product': False, 'first_product': None,
                    'record': {'force': force}}

        with TemporaryDirectory() as directory, \
                patch.object(bimol, '_trial', side_effect=trial):
            bimol._force_path(
                encounter, [[0, 3]], [[3, 2]], ZeroCalculator(), settings,
                Path(directory))
        self.assertEqual(forces, [.1, .2])
        self.assertEqual(directories, ['001_coarse', '002_coarse'])

    def test_expected_graphs_are_order_independent(self):
        forward = combined(self.methane, self.hydroxyl)
        reverse = combined(self.hydroxyl, self.methane)
        reactant, product = bimol._expected_graphs(
            forward, [[5, 1]], [[1, 0]])
        self.assertIn((0, 1), reactant)
        self.assertNotIn((0, 1), product)
        self.assertIn((1, 5), product)
        reactant, product = bimol._expected_graphs(
            reverse, [[0, 3]], [[3, 2]])
        self.assertIn((2, 3), reactant)
        self.assertNotIn((2, 3), product)
        self.assertIn((0, 3), product)

    def test_reaction_finder_is_order_independent(self):
        parameters = {
            'bimol_families': ['abstraction'],
            'skip_families': ['none'],
            'specific_reaction': 0,
            'break_bonds': [],
            'form_bonds': [],
            'ringrange': [3, 9],
            'one_reaction_comb': 0,
            'one_reaction_fam': 0,
            'scan_step': 20,
        }
        forward = combined(self.methane, self.hydroxyl)
        reverse = combined(self.hydroxyl, self.methane)
        ReactionFinderBimol(forward, parameters, None).find_reactions()
        ReactionFinderBimol(reverse, parameters, None).find_reactions()
        self.assertEqual(forward.reac_inst, [[5, 1, 0]])
        self.assertEqual(reverse.reac_inst, [[0, 3, 2]])
        self.assertTrue(all(type(index) is int
                            for index in forward.reac_inst[0]
                            + reverse.reac_inst[0]))
        self.assertEqual(forward.reac_obj[0].stereopath_id,
                         reverse.reac_obj[0].stereopath_id)
        self.assertEqual(forward.reac_obj[0].site_degeneracy, 4)
        order = [0, 3, 1, 4, 2]
        permuted_methane = stationary(Atoms(
            [self.methane.atom[index] for index in order],
            positions=self.methane.geom[order]), 1)
        permuted = combined(permuted_methane, self.hydroxyl)
        ReactionFinderBimol(permuted, parameters, None).find_reactions()
        self.assertEqual(forward.reac_obj[0].stereopath_id,
                         permuted.reac_obj[0].stereopath_id)

    def test_isobutane_site_degeneracies_are_one_and_nine(self):
        parameters = {
            'bimol_families': ['abstraction'],
            'skip_families': ['none'],
            'specific_reaction': 0,
            'break_bonds': [],
            'form_bonds': [],
            'ringrange': [3, 9],
            'one_reaction_comb': 0,
            'one_reaction_fam': 0,
            'scan_step': 20,
        }
        cluster_parameters = {'bimol_energy_tolerance': 0.25,
                              'bimol_basin_rmsd': 0.15,
                              'bimol_ts_rmsd': 0.10}
        isobutane = stationary(molecule('isobutane'), 1)
        for first, second in ((self.hydroxyl, isobutane),
                              (isobutane, self.hydroxyl)):
            species = combined(first, second)
            species.energy = species.zpe = 0.
            ReactionFinderBimol(species, parameters, None).find_reactions()
            network = BimolecularNetwork(species, cluster_parameters)
            degeneracies = sorted(network._site_degeneracy([reaction])
                                  for reaction in species.reac_obj)
            self.assertEqual(degeneracies, [1, 9])

    def test_oxygen_abstraction_keeps_spin_and_groups_equivalent_sites(self):
        parameters = default_parameters()
        parameters['bimol_families'] = ['abstraction']
        ethane = stationary(molecule('C2H6'), 1)
        for mult in (1, 3):
            oxygen = stationary(molecule('O2'), mult)
            self.assertEqual(oxygen.bond[0][1], 2)
            self.assertFalse(np.any(oxygen.rads))
            original_bonds = oxygen.bond.copy()
            original_radicals = copy.deepcopy(oxygen.rads)
            path_ids = []
            for first, second in ((ethane, oxygen), (oxygen, ethane)):
                with self.subTest(mult=mult, oxygen_first=first is oxygen):
                    species = combined(first, second, mult=mult)
                    ReactionFinderBimol(species, parameters, None).find_reactions()
                    self.assertEqual(len(species.reac_obj), 1)
                    reaction = species.reac_obj[0]
                    path_ids.append(reaction.stereopath_id)
                    a, h, d = reaction.instance
                    self.assertEqual([species.atom[i] for i in (a, h, d)],
                                     ['O', 'H', 'C'])
                    self.assertEqual(len(reaction.site_members['abstractor']), 2)
                    self.assertEqual(len(reaction.site_members['transferred_atom']), 6)
                    self.assertEqual(reaction.site_degeneracy, 12)
                    self.assertEqual(_roles(species, reaction.instance)[0][0],
                                     0 if first is oxygen else 1)
                    self.assertEqual(species.mult, mult)
                    self.assertEqual(oxygen.mult, mult)
                    np.testing.assert_array_equal(oxygen.bond, original_bonds)
                    np.testing.assert_array_equal(oxygen.rads, original_radicals)
                    encounter, _ = next(generate_encounters(
                        species, first.geom, second.geom, reaction.instance,
                        bimol.DEFAULTS.copy()))
                    self.assertTrue(np.all(np.isfinite(encounter.positions)))
                    oxygen_atoms = reaction.site_members['abstractor']
                    self.assertAlmostEqual(encounter.get_distance(*oxygen_atoms),
                                           np.linalg.norm(oxygen.geom[0] - oxygen.geom[1]))
                    data = reaction.get_bimol_search_data()
                    self.assertEqual(data['formed_bonds'], [[a, h]])
                    self.assertEqual(data['broken_bonds'], [[h, d]])
            self.assertEqual(path_ids[0], path_ids[1])

    def test_oxygen_exception_does_not_enable_other_closed_shell_abstractors(self):
        parameters = default_parameters()
        parameters['bimol_families'] = ['abstraction']
        for formula in ('N2', 'H2O2'):
            with self.subTest(formula=formula):
                abstractor = stationary(molecule(formula), 1)
                species = combined(self.methane, abstractor, mult=1)
                ReactionFinderBimol(species, parameters, None).find_reactions()
                self.assertEqual(species.reac_inst, [])

    def test_search_retains_every_saddle_without_running_irc(self):
        species = combined(self.hydroxyl, self.methane)
        encounter = Atoms(species.atom, positions=species.geom)
        candidate = {'atoms': encounter, 'mode': np.ones(3 * species.natom)}
        poses = [(encounter, {'pose': index, 'angle': 180., 'r_ah': 2.})
                 for index in range(2)]
        result = {'saddle': True, 'converged': True}
        with TemporaryDirectory() as directory, \
                patch.object(bimol, '_collision_free', return_value=True), \
                patch.object(bimol, '_force_path',
                             return_value=([candidate], {'status': 'ok'})), \
                patch.object(bimol, '_refine_candidate',
                             return_value=(result, encounter)), \
                patch.object(bimol, '_irc_endpoints',
                             side_effect=AssertionError):
            status, geometries, metadata = bimol.find_ts_candidates(
                species, [self.hydroxyl.geom], [self.methane.geom],
                [0, 3, 2], [[0, 3]], [[3, 2]],
                lambda *args: iter(poses), ZeroCalculator(), directory,
                {'sequential': False})
        self.assertEqual(status, 0)
        self.assertEqual(len(geometries), 2)
        self.assertEqual(len(metadata['candidates']), 2)

    def test_one_failed_pose_does_not_discard_other_saddles(self):
        species = combined(self.hydroxyl, self.methane)
        encounter = Atoms(species.atom, positions=species.geom)
        candidate = {'atoms': encounter, 'mode': np.ones(3 * species.natom)}
        poses = [(encounter.copy(), {'pose': index, 'angle': 180.})
                 for index in range(2)]
        with TemporaryDirectory() as directory, \
                patch.object(bimol, '_collision_free', return_value=True), \
                patch.object(bimol, '_force_path', side_effect=[
                    RuntimeError('one bad pose'),
                    ([candidate], {'status': 'ok'}),
                ]), \
                patch.object(bimol, '_refine_candidate', return_value=(
                    {'saddle': True, 'converged': True}, encounter)):
            status, geometries, metadata = bimol.find_ts_candidates(
                species, [self.hydroxyl.geom], [self.methane.geom],
                [0, 3, 2], [[0, 3]], [[3, 2]],
                lambda *args: iter(poses), ZeroCalculator(), directory,
                {'sequential': False})
        self.assertEqual(status, 0)
        self.assertEqual(len(geometries), 1)
        self.assertEqual(metadata['attempts'][0]['status'], 'error')
        self.assertEqual(metadata['attempts'][1]['conformers'], [0, 0])

    def test_l0_validation_retains_an_unexpected_product(self):
        species = combined(self.hydroxyl, self.methane)
        endpoints = [
            {'direction': 'forward', 'state': 'reactant',
             'converged': True, 'bonds': [], 'geometry_angstrom': species.geom.tolist()},
            {'direction': 'reverse', 'state': 'unexpected',
             'converged': True, 'bonds': [[0, 2]], 'geometry_angstrom': species.geom.tolist()},
        ]
        with TemporaryDirectory() as directory, patch.object(
                bimol, '_irc_endpoints', return_value=endpoints):
            result = bimol.validate_candidate(
                species, species.geom, [[0, 3]], [[3, 2]],
                ZeroCalculator(), directory)
        self.assertTrue(result['valid'])
        self.assertFalse(result['product_matches_expected'])

    @patch('kinbot.qc.connect')
    def test_l0_command_is_independent_from_l1(self, _):
        parameters = default_parameters()
        parameters.update({'bimol': 1, 'qc_command': '/opt/qchem',
                           'l0_qc_command': '/opt/qchem-l0'})
        qc = QuantumChemistry(parameters)
        self.assertEqual(qc.qc_command, '/opt/qchem')
        self.assertEqual(qc.for_level('l0').qc_command, '/opt/qchem-l0')

        parameters['l0_qc_command'] = ''
        Path('different_command').mkdir()
        os.chdir('different_command')
        qc = QuantumChemistry(parameters)
        self.assertEqual(qc.for_level('l0').qc_command, 'g16')

    def test_l0_backend_default_is_scoped_to_bimolecular_runs(self):
        with TemporaryDirectory() as directory:
            inherited = Path(directory, 'inherited.json')
            inherited.write_text(json.dumps({
                'qc': 'qchem', 'qc_command': 'qchem',
                'barrier_threshold': 50.,
                'semi_emp_conformer_search': 1,
                'semi_emp_method': 'pm7',
            }))
            parameters = Parameters(inherited).par
            self.assertEqual(parameters['l0_qc'], '')
            from kinbot.theory import resolve_profiles
            self.assertEqual(resolve_profiles(parameters), {})
            bimol_profiles = resolve_profiles({**parameters, 'bimol': 1})
            self.assertEqual(bimol_profiles['l0'].calculator, 'gaussian')
            self.assertEqual(parameters['l0_conformer_search'], 1)
            self.assertEqual(parameters['l0_method'], 'pm7')

            explicit = Path(directory, 'explicit.json')
            explicit.write_text(json.dumps({
                'qc': 'qchem', 'qc_command': 'qchem',
                'barrier_threshold': 50.,
                'l0_qc': 'orca', 'l0_qc_command': 'orca',
            }))
            self.assertEqual(Parameters(explicit).par['l0_qc'], 'orca')

    def test_fairchem_l1_rotors_ignore_method_basis_placeholders(self):
        with TemporaryDirectory() as directory:
            path = Path(directory, 'fairchem.json')
            path.write_text(json.dumps({
                'barrier_threshold': 50.,
                'qc': 'fc',
                'qc_command': 'python',
                'fc_model_path': 'uma-m-1p1',
                'rotor_scan': 1,
            }))
            parameters = Parameters(path).par
        self.assertEqual(parameters['qc'], 'fc')
        self.assertNotEqual(parameters['method'].lower(),
                            parameters['high_level_method'].lower())

    def test_l2_without_conformer_search_requires_expert_opt_in(self):
        with TemporaryDirectory() as directory:
            path = Path(directory, 'l2.json')
            path.write_text(json.dumps({
                'barrier_threshold': 50.,
                'high_level': 1,
                'allow_l2_without_conf': 1,
            }))
            parameters = Parameters(path).par
        self.assertEqual(parameters['allow_l2_without_conf'], 1)
        self.assertEqual(parameters['conformer_search'], 0)

    @patch('kinbot.qc.connect')
    def test_l0_sella_defaults_do_not_mutate_l1_settings(self, _):
        parameters = default_parameters()
        parameters.update({'qc': 'gauss', 'use_sella': 0,
                           'l0_qc': 'orca', 'sella_kwargs': {}})
        qc = QuantumChemistry(parameters)
        level = qc.for_level('l0')
        self.assertEqual(qc.par['sella_kwargs'], {})
        self.assertTrue(level.par['sella_kwargs']['internal'])
        self.assertFalse(qc.use_sella)
        self.assertTrue(level.use_sella)

    @patch('kinbot.qc.connect')
    def test_l0_conformer_submission_uses_the_l0_backend(self, _):
        parameters = default_parameters()
        parameters.update({'qc': 'qchem', 'qc_command': 'qchem',
                           'l0_qc': 'orca', 'l0_qc_command': 'orca'})
        qc = QuantumChemistry(parameters)
        species = copy.deepcopy(self.methane)
        with TemporaryDirectory() as directory:
            old = Path.cwd()
            try:
                os.chdir(directory)
                Path('conf').mkdir()
                with patch.object(
                        QuantumChemistry, 'submit_qc', autospec=True,
                        return_value=0) as submit:
                    qc = QuantumChemistry(parameters)
                    qc.qc_conf(species, species.geom, 0, semi_emp=1)
                self.assertEqual(submit.call_args.args[0].qc, 'orca')
            finally:
                os.chdir(old)

    @patch('kinbot.qc.connect')
    def test_l0_conformer_state_is_initialized(self, _):
        parameters = default_parameters()
        qc = QuantumChemistry(parameters)
        optimization = Optimize(
            copy.deepcopy(self.methane), parameters, qc)
        self.assertEqual(optimization.ssemi_empconf, -1)

    def test_bimolecular_ts_can_defer_hir_until_after_irc(self):
        parameters = default_parameters()
        parameters['rotor_scan'] = 1
        reactant = combined(self.hydroxyl, self.methane)
        product = copy.deepcopy(reactant)
        product.bond[2, 3] = product.bond[3, 2] = 0
        product.bond[0, 3] = product.bond[3, 0] = 1
        product.bond01 = (product.bond > 0).astype(int)
        product.bonds = [product.bond.copy()]
        product.characterize(bond_mx=product.bond)
        transition_state = copy.deepcopy(reactant)
        transition_state.wellorts = 1
        reaction = SimpleNamespace(ts=transition_state,
                                   irc_reactant=reactant, irc_prod=product)
        bimol.update_ts_topology_from_irc(reaction)
        transition_state.freq = [-300., 100.]
        row = SimpleNamespace(
            id=1, name=transition_state.name, positions=transition_state.geom.copy(),
            data={'energy': -100. / constants.EVtoHARTREE, 'zpe': .1,
                  'frequencies': [-400., 200.], 'status': 'normal'})
        qc = SimpleNamespace(db=SimpleNamespace(select=lambda **kwargs: [row]))
        optimization = Optimize(
            transition_state, parameters, qc, defer_hir=True)
        with patch('kinbot.optimize.HIR', side_effect=AssertionError):
            optimization.do_optimization()
        self.assertEqual(optimization.shir, 1)
        self.assertIsNone(transition_state.hir)
        self.assertEqual(transition_state.source_job, row.name)
        self.assertAlmostEqual(transition_state.energy, -100.)
        self.assertEqual(transition_state.freq, [-400., 200.])
        self.assertEqual(transition_state.reduced_freqs, [-400., 200.])

    @patch('kinbot.qc.connect')
    def test_final_gaussian_irc_uses_selected_geometry_and_fresh_hessian(self, _):
        parameters = default_parameters()
        parameters.update(bimol=1, qc='gauss', use_sella=0)
        qc = QuantumChemistry(parameters)
        species = combined(self.hydroxyl, self.methane)
        reaction = Abstraction(species, qc, parameters, [0, 3, 2], 'candidate')
        reaction.ts = copy.deepcopy(species)
        reaction.ts.geom = species.geom + .123
        reaction.final_irc_name = 'candidate_final_0123456789'
        with TemporaryDirectory() as directory:
            old = Path.cwd()
            try:
                os.chdir(directory)
                Path('candidate.chk').write_text('discovery checkpoint')
                with patch.object(qc.for_level('l1'), 'submit_qc', return_value=0):
                    IRC(reaction, parameters).do_irc_calculations()
                for direction in ('F', 'R'):
                    job = f'{reaction.final_irc_name}_IRC_{direction}'
                    tree = ast.parse(Path(job + '.py').read_text())
                    assignments = {node.targets[0].id: node.value for node in tree.body
                                   if isinstance(node, ast.Assign)
                                   and isinstance(node.targets[0], ast.Name)}
                    kwargs = ast.literal_eval(assignments['kwargs'])
                    options = kwargs['irc'].split(',')
                    self.assertIn('CalcFC', options)
                    self.assertNotIn('RCFC', options)
                    self.assertNotIn('geom', kwargs)
                    self.assertNotIn('guess', kwargs)
                    positions = next(key.value for key in assignments['mol'].keywords
                                     if key.arg == 'positions')
                    np.testing.assert_allclose(ast.literal_eval(positions), reaction.ts.geom)
                    self.assertFalse(Path(job + '.chk').exists())
            finally:
                os.chdir(old)

    @patch('kinbot.qc.connect')
    def test_qc_opt_ts_honors_requested_theory_level(self, _):
        parameters = default_parameters()
        qc = QuantumChemistry(parameters)
        species = combined(self.hydroxyl, self.methane)
        species.name = 'selected_ts'
        with TemporaryDirectory() as directory:
            old = Path.cwd()
            try:
                os.chdir(directory)
                with patch.object(qc, 'get_qc_arguments', return_value={}) \
                        as arguments, patch.object(
                            qc, 'submit_qc', return_value=0):
                    qc.qc_opt_ts(species, species.geom, high_level=0)
                    self.assertFalse(arguments.call_args.kwargs['high_level'])
                    qc.qc_opt_ts(species, species.geom, high_level=1)
                    self.assertTrue(arguments.call_args.kwargs['high_level'])
            finally:
                os.chdir(old)

    @patch('kinbot.qc.connect')
    def test_final_sella_irc_script_uses_selected_geometry(self, _):
        parameters = default_parameters()
        parameters['bimol'] = 1
        parameters['use_sella'] = 1
        parameters['sella_kwargs']['internal'] = False
        qc = QuantumChemistry(parameters)
        qc.use_sella = True
        species = combined(self.hydroxyl, self.methane)
        reaction = Abstraction(
            species, qc, parameters, [0, 3, 2], 'candidate')
        reaction.ts = copy.deepcopy(species)
        reaction.ts.geom = species.geom + .123
        reaction.final_irc_name = 'candidate_final_0123456789'
        with TemporaryDirectory() as directory:
            old = Path.cwd()
            try:
                os.chdir(directory)
                with patch.object(qc.for_level('l1'), 'submit_qc', return_value=0):
                    IRC(reaction, parameters).do_irc_calculations()
                for direction in ('F', 'R'):
                    script = Path(
                        f'{reaction.final_irc_name}_IRC_{direction}.py')
                    compile(script.read_text(), str(script), 'exec')
                    self.assertIn('calc_vibrations', script.read_text())
                self.assertIn('.123', Path(
                    f'{reaction.final_irc_name}_IRC_F.py').read_text())
            finally:
                os.chdir(old)

    @patch('kinbot.qc.connect')
    def test_l0_search_job_files_are_nested_under_bimol(self, _):
        parameters = default_parameters()
        qc = QuantumChemistry(parameters)
        species = combined(self.hydroxyl, self.methane)
        reaction = Abstraction(species, qc, parameters, [0, 3, 2], 'search')
        with TemporaryDirectory() as directory:
            old = Path.cwd()
            try:
                os.chdir(directory)
                with patch.object(qc, 'submit_qc', return_value=1):
                    qc.qc_bimol_search(reaction)
                script = Path(bimol.search_job('search', sequential=True) + '.py')
                self.assertTrue(script.is_file())
                self.assertTrue(Path('perm/bimol/search').is_dir())
                compile(script.read_text(), str(script), 'exec')
                self.assertIn("_bimol_guess_sequential.done'", script.read_text())
                self.assertIn('os.replace(pkl_temporary, pkl_file)',
                              script.read_text())
                self.assertIn('os.replace(done_temporary, done_file)',
                              script.read_text())
                candidate = {
                    'validation_directory':
                        'bimol/unique_candidates/candidate_000',
                    'candidate': {'geometry': species.geom.tolist()},
                }
                with patch.object(qc, 'submit_qc', return_value=1):
                    qc.qc_bimol_irc(reaction, candidate)
                irc_script = Path(candidate['validation_directory'],
                                  'candidate_000_bimol_irc.py')
                compile(irc_script.read_text(), str(irc_script), 'exec')
                self.assertIn("_bimol_irc.done'", irc_script.read_text())
                self.assertIn('os.replace(pkl_temporary, pkl_file)',
                              irc_script.read_text())
                self.assertIn('os.replace(done_temporary, done_file)',
                              irc_script.read_text())
            finally:
                os.chdir(old)

    @patch('kinbot.qc.connect')
    def test_queue_limit_one_allows_the_first_job(self, _):
        parameters = default_parameters()
        parameters.update({'queuing': 'slurm', 'queue_job_limit': 1,
                           'username': 'guest'})
        qc = QuantumChemistry(parameters)
        with patch('kinbot.qc.subprocess.check_output', return_value=b''):
            self.assertEqual(qc.limit_jobs(), 0)

    def test_failed_l1_refinement_keeps_valid_l0_ts(self):
        with TemporaryDirectory() as directory:
            old = Path.cwd()
            path = Path(directory, 'bimol/search/c0000_0001/p002/candidate_03')
            try:
                os.chdir(directory)
                path.mkdir(parents=True)
                Path(path, 'result.json').write_text(
                    json.dumps({'saddle': True}))
                reaction = SimpleNamespace(
                    search_parent='search',
                    search_candidate={'directory': str(path)})
                bimol.update_bimol_status(reaction, l1_valid=False,
                                          final_status='failed')
            finally:
                os.chdir(old)
            result = json.loads(Path(path, 'result.json').read_text())
            self.assertFalse(result['l1_valid'])
            old = Path.cwd()
            try:
                os.chdir(directory)
                bimol.update_bimol_status(reaction, l1_valid=True,
                                          final_status='final_irc_valid')
            finally:
                os.chdir(old)
            result = json.loads(Path(path, 'result.json').read_text())
            self.assertTrue(result['l1_valid'])

    def test_irc_returns_both_bimolecular_endpoints(self):
        species = combined(self.hydroxyl, self.methane)
        reactant = species.geom.copy()
        product = reactant.copy()
        product[3] = product[0] + np.array([0., 0., 0.96])
        product[2:] += np.array([2., 0., 0.])
        product[3] = product[0] + np.array([0., 0., 0.96])

        class QC:
            def for_level(self, level):
                return self

            def get_qc_geom(self, name, natom, allow_error=1):
                return 0, reactant if '_F_' in name else product

        reaction = Abstraction(
            species, QC(), {'scan_step': 20}, [0, 3, 2], 'test')
        result = IRC(reaction, {'bimol': 1}).irc2stationary_pt()
        self.assertIsInstance(result, BimolecularIRCResult)
        self.assertEqual(result.reactant_direction, 'F')
        self.assertEqual(result.product_direction, 'R')

    def test_irc_retains_an_unexpected_bimolecular_product(self):
        species = combined(self.hydroxyl, self.methane)
        reactant = species.geom.copy()
        unexpected = reactant.copy()
        unexpected[3] = unexpected[0] + np.array([0., .96, 0.])
        unexpected[2] = unexpected[0] + np.array([0., 0., 1.35])

        class QC:
            def for_level(self, level):
                return self

            def get_qc_geom(self, name, natom, allow_error=1):
                return 0, reactant if '_F_' in name else unexpected

        reaction = Abstraction(
            species, QC(), {'scan_step': 20}, [0, 3, 2], 'test')
        result = IRC(reaction, {'bimol': 1}).irc2stationary_pt()
        self.assertIsInstance(result, BimolecularIRCResult)
        self.assertFalse(result.product_matches_expected)

    def test_abstraction_rotor_bypasses_transferred_hydrogen(self):
        atom = ['C', 'O', 'O', 'H', 'H', 'H',
                'C', 'O', 'H', 'H', 'H', 'H']
        geom = np.array([
            [3.4893204, .6036982, -.4948179],
            [3.1696910, -.6841238, .0100114],
            [2.2738556, -.5269153, 1.0593594],
            [4.2394292, .4348445, -1.2661997],
            [3.9002077, 1.2307904, .2975766],
            [2.6070159, 1.0776537, -.9300171],
            [.0254758, -.5559212, -.0919514],
            [-.1666782, .6737621, -.6835423],
            [1.1902319, -.5706644, .5402898],
            [.1591444, -1.3220532, -.8503892],
            [-.6882728, -.8151854, .6908030],
            [-.4462984, 1.3115785, -.0191018],
        ])
        reactant_bonds = [
            (0, 1), (1, 2), (0, 3), (0, 4), (0, 5),
            (6, 7), (6, 8), (6, 9), (6, 10), (7, 11),
        ]

        def point(name, bonds, wellorts=0):
            species = StationaryPoint(
                name, 0, 2, atom=atom, geom=geom, wellorts=wellorts)
            matrix = np.zeros((len(atom), len(atom)), dtype=int)
            for first, second in bonds:
                matrix[first, second] = matrix[second, first] = 1
            species.bond = matrix
            species.bond01 = matrix.copy()
            species.bonds = [matrix.copy()]
            species.characterize(bond_mx=matrix)
            return species

        reactant = point('reactant', reactant_bonds)
        product = point(
            'product',
            [bond for bond in reactant_bonds if bond != (7, 11)]
            + [(2, 11)])
        reaction = Abstraction(
            reactant, None, {'scan_step': 20}, [2, 11, 7], 'test')
        reaction.irc_reactant = reactant
        reaction.irc_prod = product
        reaction.ts = point('ts', reactant_bonds + [(2, 11)], wellorts=1)

        rotor = reaction.prepare_hir_rotors()

        self.assertEqual(rotor, [6, 7, 2, 1])
        self.assertTrue(all(11 not in item for item in reaction.ts.dihed))
        _, _, _, order = zmatrix.make_zmat_from_cart(
            reaction.ts, reaction.ts.dihed.index(rotor), geom, 0)
        self.assertEqual(sorted(order), list(range(len(atom))))

    def test_methane_abstraction_rotor_uses_another_hydrogen(self):
        reactant = combined(self.hydroxyl, self.methane)
        product = copy.deepcopy(reactant)
        product.bond[2, 3] = product.bond[3, 2] = 0
        product.bond[0, 3] = product.bond[3, 0] = 1
        product.bond01 = (product.bond > 0).astype(int)
        product.bonds = [product.bond.copy()]
        product.characterize(bond_mx=product.bond)
        reaction = Abstraction(
            reactant, None, {'scan_step': 20}, [0, 3, 2], 'test')
        reaction.irc_reactant = reactant
        reaction.irc_prod = product
        reaction.ts = copy.deepcopy(reactant)
        bimol.update_ts_topology_from_irc(reaction)

        rotor = reaction.prepare_hir_rotors()

        self.assertEqual(rotor[1:3], [2, 0])
        self.assertEqual(reaction.ts.atom[rotor[0]], 'H')
        self.assertEqual(reaction.ts.atom[rotor[3]], 'H')
        self.assertNotEqual(rotor[0], 3)
        self.assertNotIn(3, rotor)

    def test_terminal_transferred_hydrogen_preserves_rotor_axis(self):
        atom = ['H', 'C', 'O', 'H', 'H', 'H', 'H']
        geom = np.array([
            [1.1835208457, -1.9348098265, 1.2141312249],
            [-.1311901453, -.6527442162, -.1266704555],
            [.0478598973, .6695002198, .2355193476],
            [.5991135594, -1.4151544987, .6269793071],
            [.1867250011, -.9039052140, -1.1396078245],
            [-1.1487191912, -.9610381242, .0923166505],
            [.9109040165, .9707515896, -.0637243247],
        ])

        def point(name, bonds):
            species = StationaryPoint(
                name, 0, 2, atom=atom, geom=geom, wellorts=1)
            matrix = np.zeros((len(atom), len(atom)), dtype=int)
            for first, second in bonds:
                matrix[first, second] = matrix[second, first] = 1
            species.bond = matrix
            species.bond01 = matrix.copy()
            species.bonds = [matrix.copy()]
            species.characterize(bond_mx=matrix)
            return species

        methanol = [(1, 2), (1, 3), (1, 4), (1, 5), (2, 6)]
        reactant = point('reactant', methanol)
        product = point('product', [bond for bond in methanol
                                    if bond != (1, 3)] + [(0, 3)])
        reaction = Abstraction(
            reactant, None, {'scan_step': 20}, [0, 3, 1], 'test')
        reaction.irc_reactant = reactant
        reaction.irc_prod = product
        reaction.ts = point('ts', methanol + [(0, 3)])

        rotor = reaction.prepare_hir_rotors()

        self.assertEqual(rotor[1:3], [1, 2])
        self.assertIn(rotor[0], [4, 5])
        self.assertEqual(rotor[3], 6)
        self.assertNotIn(3, rotor)

    def test_stereotopic_hydrogen_classes_are_refined_selectively(self):
        ethanol = StationaryPoint('ethanol', 0, 1, smiles='CCO')
        ethanol.characterize()
        ethanol_ch2 = next(group for group in ethanol.atom_eqv
                           if len(group) == 2
                           and ethanol.atom[group[0]] == 'H')
        classes = refine_equivalence_group(ethanol, ethanol_ch2)
        self.assertEqual(len(classes), 1)
        self.assertEqual(classes[0]['relation'], 'enantiotopic')

        chiral = StationaryPoint('chiral', 0, 1, smiles='CC(O)CC')
        chiral.characterize()
        chiral_ch2 = next(group for group in chiral.atom_eqv
                          if len(group) == 2
                          and chiral.atom[group[0]] == 'H')
        classes = refine_equivalence_group(chiral, chiral_ch2)
        self.assertEqual([len(item['members']) for item in classes], [1, 1])
        self.assertTrue(all(item['relation'] == 'diastereotopic'
                            for item in classes))
        refined = stereotopic_hydrogen_equivalence(chiral)
        self.assertTrue(all(not set(chiral_ch2).issubset(group)
                            for group in refined))

    def test_abstraction_site_degeneracy_uses_stereotopic_class(self):
        hydrogen = stationary(Atoms('H', positions=[[0., 0., 0.]]), 2)
        species = combined(hydrogen, stationary(molecule('CH3OH'), 1))
        par = default_parameters()
        par.update({'bimol_families': ['abstraction'],
                    'skip_families': []})
        ReactionFinderBimol(species, par, None).find_reactions()
        methyl = next(reaction for reaction in species.reac_obj
                      if species.atom[reaction.instance[2]] == 'C')
        hydroxyl = next(reaction for reaction in species.reac_obj
                        if species.atom[reaction.instance[2]] == 'O')
        self.assertEqual(methyl.site_degeneracy, 3)
        self.assertEqual(hydroxyl.site_degeneracy, 1)
        self.assertNotEqual(methyl.stereopath_id, hydroxyl.stereopath_id)

    def test_configured_reactant_pairs_have_distinct_reaction_names(self):
        hydrogen = stationary(Atoms('H', positions=[[0., 0., 0.]]), 2)
        pairs = []
        for smiles in ('C[C@H](O)[C@H](F)C', 'C[C@H](O)[C@@H](F)C'):
            fragment = StationaryPoint('configured', 0, 1, smiles=smiles)
            fragment.characterize()
            species = combined(hydrogen, fragment)
            par = default_parameters()
            par.update({'bimol_families': ['abstraction'], 'skip_families': []})
            ReactionFinderBimol(species, par, None).find_reactions()
            self.assertTrue(species.reac_name)
            self.assertTrue(all(name.startswith(routing_name(species) + '_abstraction_')
                                for name in species.reac_name))
            pairs.append(species)
        self.assertEqual(pairs[0].chemid, pairs[1].chemid)
        self.assertTrue(set(pairs[0].reac_name).isdisjoint(pairs[1].reac_name))

    def test_abstraction_enumerates_diastereotopic_hydrogens_separately(self):
        hydrogen = stationary(Atoms('H', positions=[[0., 0., 0.]]), 2)
        donor = StationaryPoint('chiral', 0, 1, smiles='CC(O)CC')
        donor.characterize()
        ch2_hydrogens = next(group for group in donor.atom_eqv
                             if len(group) == 2
                             and donor.atom[group[0]] == 'H')
        species = combined(hydrogen, donor)
        par = default_parameters()
        par.update({'bimol_families': ['abstraction'],
                    'skip_families': []})
        ReactionFinderBimol(species, par, None).find_reactions()
        reactions = [reaction for reaction in species.reac_obj
                     if reaction.instance[1] - 1 in ch2_hydrogens]
        self.assertEqual(len(reactions), 2)
        self.assertEqual({reaction.site_degeneracy for reaction in reactions},
                         {1})
        self.assertEqual(len({reaction.stereopath_id
                              for reaction in reactions}), 2)
        self.assertTrue(all(
            reaction.site_relations['transferred_atom'] == 'diastereotopic'
            for reaction in reactions))

    def test_chiral_partner_splits_an_enantiotopic_donor_pair(self):
        abstractor = StationaryPoint(
            'chiral_radical', 0, 2, smiles='[CH2]C(O)C')
        abstractor.characterize()
        donor = StationaryPoint('ethanol', 0, 1, smiles='CCO')
        donor.characterize()
        ch2_hydrogens = next(group for group in donor.atom_eqv
                             if len(group) == 2
                             and donor.atom[group[0]] == 'H')
        species = combined(abstractor, donor)
        par = default_parameters()
        par.update({'bimol_families': ['abstraction'],
                    'skip_families': []})
        ReactionFinderBimol(species, par, None).find_reactions()
        reactions = [reaction for reaction in species.reac_obj
                     if reaction.instance[1] - abstractor.natom
                     in ch2_hydrogens]
        self.assertEqual(len(reactions), 2)
        self.assertTrue(all(not reaction.allow_reflection
                            for reaction in reactions))
        self.assertTrue(all(reaction.site_degeneracy == 1
                            for reaction in reactions))

    def test_stereopath_labels_protect_candidate_clustering(self):
        species = combined(self.hydroxyl, self.methane)
        first = SimpleNamespace(
            family_name='abstraction', instance=[0, 3, 2],
            instance_name='first', stereopath_id='path-a')
        second = SimpleNamespace(
            family_name='abstraction', instance=[0, 3, 2],
            instance_name='second', stereopath_id='path-b')
        candidate = {'geometry': species.geom.tolist(),
                     'result': {'energy': -1.}}
        groups = bimol.cluster_candidates(species, [
            {'reaction': first, 'candidate': candidate},
            {'reaction': second, 'candidate': copy.deepcopy(candidate)},
        ], .1, .25)
        self.assertEqual(len(groups), 2)

    def test_meoh_h_optical_pair_is_detected_from_one_ts(self):
        atom = ['H', 'C', 'O', 'H', 'H', 'H', 'H']
        geometries = [
            np.array([
                [1.1835208457, -1.9348098265, 1.2141312249],
                [-.1311901453, -.6527442162, -.1266704555],
                [.0478598973, .6695002198, .2355193476],
                [.5991135594, -1.4151544987, .6269793071],
                [.1867250011, -.9039052140, -1.1396078245],
                [-1.1487191912, -.9610381242, .0923166505],
                [.9109040165, .9707515896, -.0637243247],
            ]),
            np.array([
                [1.4808507459, 1.5145045743, -1.0093746777],
                [-.0208370909, -.6393750782, -.1409279628],
                [-.1506027260, .7441295391, .0534346826],
                [.9039927188, -1.0380903984, .2830089931],
                [-.1093756215, -.9246895406, -1.1919996296],
                [-.8610744344, -1.0810138823, .4022870813],
                [.7693132164, 1.2952935227, -.5305860609],
            ]),
        ]
        expected = [2, 1]
        transfers = [3, 6]
        for geom, count, transferred in zip(
                geometries, expected, transfers):
            species = StationaryPoint(
                'ts', 0, 2, atom=atom, geom=geom, wellorts=1)
            bonds = [(1, 2), (1, 3), (1, 4), (1, 5), (2, 6),
                     (0, transferred)]
            matrix = np.zeros((len(atom), len(atom)), dtype=int)
            for first, second in bonds:
                matrix[first, second] = matrix[second, first] = 1
            species.bond = matrix
            species.bond01 = matrix.copy()
            species.bonds = [matrix.copy()]
            species.characterize(bond_mx=matrix)
            optical, mirror_rmsd = optical_isomers(species)
            self.assertEqual(optical, count)
            self.assertGreater(mirror_rmsd, .4 if count == 2 else 0.)
            if count == 1:
                self.assertLess(mirror_rmsd, .01)
            symmetry.calculate_symmetry(species)
            metadata = assign_optical_metadata(species)
            self.assertEqual(metadata['nopt'], 1)
            self.assertEqual(metadata['optical_isomers'], count)
            if count == 2:
                species.conformer_geom = [geom, geom * [-1., 1., 1.]]
                species.conformer_index = [2, 7]
                conformers = conformer_symmetry(species)
                self.assertEqual([item['index'] for item in conformers],
                                 [2, 7])
                self.assertTrue(all(item['optical_isomers'] == 2
                                    for item in conformers))

    def test_selected_properties_come_from_one_database_row(self):
        species = combined(self.hydroxyl, self.methane)
        selected_geom = species.geom + .25
        row = SimpleNamespace(
            name='selected', id=1,
            positions=selected_geom,
            data={'energy': -27.211386245988,
                  'zpe': .0123,
                  'frequencies': [-321., 45.],
                  'hess': np.eye(3 * species.natom)})
        qc = SimpleNamespace(
            db=SimpleNamespace(select=lambda name: iter([row])))
        bimol.load_calculation_record(species, qc, 'selected')
        np.testing.assert_allclose(species.geom, selected_geom)
        self.assertAlmostEqual(species.energy, -1.)
        self.assertEqual(species.zpe, .0123)
        self.assertEqual(species.freq, [-321., 45.])
        self.assertEqual(species.reduced_freqs, species.freq)
        self.assertTrue(species.hess)
        self.assertEqual(species.source_job, 'selected')

    def test_ts_seed_topology_uses_l0_endpoints_not_encounter_contacts(self):
        species = combined(self.hydroxyl, self.methane)
        geometry = species.geom.copy()
        geometry[0] = geometry[2] + np.array([0., 0., 1.])
        _, product_graph = bimol._expected_graphs(
            species, [[0, 3]], [[3, 2]])
        row = SimpleNamespace(
            name='test', id=1,
            positions=geometry,
            data={'energy': -27.211386245988, 'zpe': .01,
                  'frequencies': [-300., 50.],
                  'hess': np.eye(3 * species.natom)})
        qc = SimpleNamespace(
            db=SimpleNamespace(select=lambda name: iter([row])),
            get_qc_geom=lambda name, natom: (0, geometry))
        reaction = Abstraction(
            species, qc, {'scan_step': 20}, [0, 3, 2], 'test')
        from kinbot.stereo_identity import canonical_identity
        product_geometry = species.geom.copy()
        product_geometry[3] = product_geometry[0] + np.array([1., 0., 0.])
        product_endpoint = StationaryPoint(
            'product', 0, 2, atom=species.atom, geom=product_geometry)
        product_endpoint.characterize()
        reaction.l0_validation = {'endpoints': [
            {'state': 'reactant',
             'configured_identity': canonical_identity(species),
             'geometry_angstrom': species.geom.tolist()},
            {'state': 'product',
             'configured_identity': canonical_identity(product_endpoint),
             'geometry_angstrom': product_geometry.tolist(),
             'bonds': [list(bond) for bond in sorted(product_graph)]},
        ]}
        ts = bimol.build_bimolecular_ts(reaction, qc)
        self.assertEqual(ts.bond01[0][2], 0)
        self.assertEqual(ts.bond01[0][3], 1)

    def test_symmetry_rmsd_handles_equivalent_atom_permutations(self):
        first = stationary(molecule('CH4'), 1)
        order = [0, 2, 1, 3, 4]
        second = stationary(Atoms(
            [first.atom[index] for index in order],
            positions=first.geom[order]), 1)
        self.assertLess(symmetry_rmsd(first, second), 1.e-5)

    def test_bimolecular_projection_respects_hessian_weighting(self):
        from kinbot.calculation import geometry_reference
        from kinbot.hindered_rotors import HIR
        species = stationary(molecule('H2O'), 1)
        hessian = np.diag(np.arange(1., 10.))
        masses = np.repeat([constants.exact_mass[atom] for atom in species.atom], 3)
        expected, _ = frequencies.get_frequencies(species, hessian, species.geom)
        par = default_parameters()
        par['rotor_scan'] = 1
        for use_sella in (False, True):
            for stage in ('selected_properties', 'final_hir'):
                with self.subTest(use_sella=use_sella, stage=stage):
                    species = stationary(molecule('H2O'), 1)
                    qc = QuantumChemistry.__new__(QuantumChemistry)
                    qc.par = par
                    qc.qc = 'qchem'
                    qc.use_sella = use_sella
                    native_hessian = (hessian if use_sella else
                                      hessian / np.sqrt(np.outer(masses, masses)))
                    row = SimpleNamespace(name='selected', id=1,
                        symbols=list(species.atom), positions=species.geom, data={
                            'energy': -10., 'zpe': .02, 'frequencies': expected})
                    if use_sella:
                        row.data['hess'] = hessian
                    qc.db = SimpleNamespace(select=lambda **kwargs: [row])
                    qc.read_qc_hess = lambda *args: native_hessian
                    optimization = Optimize.__new__(Optimize)
                    optimization.species = species
                    optimization.par = par
                    optimization.qc = qc
                    optimization.name = species.name
                    optimization.selected_job = 'selected'
                    optimization.just_high = optimization.defer_hir = False
                    optimization.shigh = optimization.shir = 1
                    optimization.restart = optimization.wait = 0
                    qc.publish_result = lambda *args: None
                    species.hir = HIR(species, qc, par)
                    species.hir.scan_reference = geometry_reference(species)
                    species.hir.scan_reference['dihedrals'] = []
                    if stage == 'selected_properties':
                        bimol.load_selected_species_properties(optimization, qc)
                    else:
                        from kinbot.calculation import load_calculation_record
                        load_calculation_record(species, qc, 'selected')
                        species.hir.check_hir = lambda: 1
                        generator = ReactionGenerator(species, par, qc, None)
                        reaction = SimpleNamespace(
                            ts=species, ts_opt=optimization)
                        self.assertTrue(generator._advance_bimol_hir(reaction))
                    np.testing.assert_allclose(species.reduced_freqs, expected, rtol=1.e-8)

    def test_symmetry_rmsd_handles_reflections(self):
        species = combined(self.hydroxyl, self.methane)
        geom = species.geom + np.random.default_rng(7).normal(
            0., .13, species.geom.shape)
        reflected = geom * [-1., 1., 1.]
        self.assertLess(symmetry_rmsd(
            species, species, geom, reflected,
            [0, 3, 2], [0, 3, 2]), 1.e-5)
        self.assertGreater(symmetry_rmsd(
            species, species, geom, reflected,
            [0, 3, 2], [0, 3, 2], allow_reflection=False), .05)

    def test_symmetry_rmsd_preserves_bond_order_labels(self):
        first = stationary(molecule('CH4'), 1)
        second = copy.deepcopy(first)
        second.bond[0][1] = second.bond[1][0] = 2
        for matrix in second.bonds:
            matrix[0][1] = matrix[1][0] = 2
        self.assertTrue(np.isinf(symmetry_rmsd(first, second)))

    def test_l0_clustering_collapses_symmetry_equivalent_sites(self):
        species = combined(self.hydroxyl, self.methane)
        first = SimpleNamespace(instance_name='first', instance=[0, 3, 2],
                                family_name='abstraction')
        second = SimpleNamespace(instance_name='second', instance=[0, 4, 2],
                                 family_name='abstraction')
        permuted = species.geom.copy()
        permuted[[3, 4]] = permuted[[4, 3]]
        entries = [
            {'reaction': first, 'candidate': {
                'geometry': species.geom.tolist(),
                'result': {'energy': -10.}}},
            {'reaction': second, 'candidate': {
                'geometry': permuted.tolist(),
                'result': {'energy': -10.}}},
        ]
        groups = bimol.cluster_candidates(
            species, entries, .1, .25)
        self.assertEqual(len(groups), 1)
        self.assertEqual(len(groups[0]), 2)

    def test_l0_clustering_collapses_mirrored_candidates(self):
        species = combined(self.hydroxyl, self.methane)
        reaction = SimpleNamespace(instance_name='first',
                                   instance=[0, 3, 2],
                                   family_name='abstraction')
        geom = species.geom + np.random.default_rng(7).normal(
            0., .13, species.geom.shape)
        entries = [
            {'reaction': reaction, 'candidate': {
                'geometry': geom.tolist(), 'result': {'energy': -10.}}},
            {'reaction': SimpleNamespace(
                instance_name='second', instance=[0, 3, 2],
                family_name='abstraction'), 'candidate': {
                    'geometry': (geom * [-1., 1., 1.]).tolist(),
                    'result': {'energy': -10.}}},
        ]
        groups = bimol.cluster_candidates(species, entries, .1, .25)
        self.assertEqual(len(groups), 1)

    def test_network_ids_ignore_fragment_and_atom_order(self):
        forward = combined(self.hydroxyl, self.methane)
        order = list(range(self.hydroxyl.natom, forward.natom)) \
            + list(range(self.hydroxyl.natom))
        reverse = StationaryPoint(
            'reverse', 0, 2,
            atom=[forward.atom[index] for index in order],
            geom=forward.geom[order], fragA=self.methane,
            fragB=self.hydroxyl)
        reverse.characterize()
        for species in (forward, reverse):
            species.energy = species.zpe = 0.
            species.reac_obj = []
            species.reac_ts_done = []
        parameters = {'bimol_energy_tolerance': 0.25,
                      'bimol_basin_rmsd': 0.15,
                      'bimol_ts_rmsd': 0.10}
        first = BimolecularNetwork(forward, parameters)
        second = BimolecularNetwork(reverse, parameters)
        self.assertEqual(set(first.nodes), set(second.nodes))

        for network, species, name in ((first, forward, 'first'),
                                       (second, reverse, 'second')):
            network._cluster_complexes([{
                'fragments': network._fragment_key(
                    [species.fragA, species.fragB]),
                'species': species,
                'reaction': SimpleNamespace(instance_name=name),
            }], 'rc')
        first_ids = {node_id for node_id in first.nodes
                     if node_id.startswith('rc_')}
        second_ids = {node_id for node_id in second.nodes
                      if node_id.startswith('rc_')}
        self.assertEqual(first_ids, second_ids)

    def test_fragment_identity_includes_charge_and_multiplicity(self):
        neutral = copy.deepcopy(self.hydroxyl)
        ion = copy.deepcopy(self.hydroxyl)
        ion.charge = 1
        self.assertNotEqual(
            BimolecularNetwork._fragment_key([neutral]),
            BimolecularNetwork._fragment_key([ion]))

    def test_fragment_identity_preserves_fixed_enantiomers(self):
        point = StationaryPoint('chiral', 0, 1, smiles='C[C@H](O)CC')
        point.characterize()
        mirror = copy.copy(point)
        mirror.geom = point.geom * [-1., 1., 1.]
        self.assertNotEqual(BimolecularNetwork._fragment_key([point]),
                            BimolecularNetwork._fragment_key([mirror]))

    def test_network_preserves_conflicting_endpoints_without_assigning_extra_flux(self):
        network = object.__new__(BimolecularNetwork)
        network.reference = 0.
        network.par = {'bimol_energy_tolerance': .25, 'bimol_ts_rmsd': .1}
        network.transitions = []
        ts = copy.deepcopy(self.methane)
        ts.energy, ts.zpe = -3., .01
        reactions = [SimpleNamespace(ts=ts, instance_name=f'candidate_{i}',
            family_name='abstraction', stereopath_id='same-site', instance=[0, 1, 2],
            rc_basin_id=f'rc_{i}', pc_basin_id='products', site_degeneracy=3,
            reactant_asymptote_id=f'reactants_{i}', product_asymptote_id='products')
            for i in range(2)]
        network._cluster_transitions(reactions)
        self.assertEqual(len(network.transitions), 2)
        for transition in network.transitions:
            self.assertEqual(len(transition['endpoint_conflicts']), 1)
            self.assertIn('unresolved', transition['flux_assignment'])
            self.assertEqual(transition['degeneracy'], 3)

    def test_final_ts_clustering_precedes_final_basin_identity(self):
        species = combined(self.hydroxyl, self.methane)
        first_ts = copy.deepcopy(species)
        second_ts = copy.deepcopy(species)
        second_ts.geom *= [-1., 1., 1.]
        for point in (first_ts, second_ts):
            point.energy = -3.
            point.zpe = .01
        reactions = [SimpleNamespace(
            bimol_search=True, family_name='abstraction',
            instance_name=f'candidate_{index}', instance=[0, 3, 2],
            l0_validation=compatible_l0_validation(),
            ts=point) for index, point in enumerate((first_ts, second_ts))]
        species.reac_obj = reactions
        species.reac_ts_done = [6, 6]
        parameters = {'bimol_energy_tolerance': .25,
                      'bimol_ts_rmsd': .1}
        self.assertTrue(bimol.release_final_ts_clusters(
            species, parameters))
        self.assertEqual(sorted(species.reac_ts_done), [7, 8])
        representative = reactions[species.reac_ts_done.index(7)]
        self.assertTrue(representative.final_irc_name.startswith(
            f'{representative.instance_name}_final_'))
        self.assertFalse(hasattr(representative, 'rc_basin_id'))

    def test_single_conf_tst_keeps_only_lowest_channel_conformer(self):
        species = combined(self.hydroxyl, self.methane)
        points = [copy.deepcopy(species), copy.deepcopy(species)]
        for index, point in enumerate(points):
            point.energy = -3. - index * 1.e-4
            point.zpe = .01
        validation = compatible_l0_validation()
        reactions = [SimpleNamespace(
            bimol_search=True, family_name='abstraction',
            instance_name=f'candidate_{index}', instance=[0, 3, 2],
            l0_validation=validation, ts=point)
            for index, point in enumerate(points)]
        species.reac_obj = reactions
        species.reac_ts_done = [6, 6]
        parameters = {'bimol_energy_tolerance': .25,
                      'bimol_ts_rmsd': .1,
                      'multi_conf_tst': 0}
        with patch('kinbot.bimol_network.symmetry_rmsd', return_value=1.):
            bimol.release_final_ts_clusters(species, parameters)
        self.assertEqual(species.reac_ts_done, [8, 7])
        self.assertEqual(reactions[1].final_ts_rank, 0)
        archive = reactions[1].final_conformer_archive
        self.assertEqual(len(archive), 2)
        self.assertEqual(reactions[0].final_conformer_role, 'archived conformer')
        self.assertEqual(archive[1]['geometry_angstrom'], points[0].geom.tolist())
        self.assertFalse(archive[1]['independent_rate_contribution'])

    def test_l0_endpoint_records_its_geometry_and_configured_identity(self):
        atoms = molecule('CH4')
        atoms.info.update(charge=0, spin=1)
        graph = bimol._bond_graph(atoms)
        with TemporaryDirectory() as directory, patch('kinbot.bimol.Sella') as optimizer:
            optimizer.return_value.run.return_value = False
            record, final = bimol._endpoint(atoms, ZeroCalculator(),
                {'endpoint_fmax': .05, 'endpoint_steps': 2},
                (graph, set(), [], []), Path(directory))
        self.assertFalse(record['converged'])  # no stricter convergence gate
        self.assertEqual(record['configured_identity']['status'], 'assigned')
        self.assertEqual(record['configured_identity']['charge'], 0)
        self.assertEqual(record['configured_identity']['multiplicity'], 1)
        np.testing.assert_array_equal(record['geometry_angstrom'], final.positions)
        atoms.info.clear()
        self.assertEqual(bimol._endpoint_configuration(atoms)['status'], 'unavailable')

    def test_unknown_or_opposite_configured_l0_products_survive_to_final_irc(self):
        species = combined(self.hydroxyl, self.methane)
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.AddHs(Chem.MolFromSmiles('CC[C@H](C)CCC'))
        AllChem.EmbedMolecule(mol, randomSeed=2)
        atoms = Atoms([a.GetSymbol() for a in mol.GetAtoms()],
                      positions=mol.GetConformer().GetPositions())
        atoms.info.update(charge=0, spin=1)
        identity = bimol._endpoint_configuration(atoms)
        reflected = atoms.copy(); reflected.positions *= [-1., 1., 1.]
        mirror = bimol._endpoint_configuration(reflected)
        self.assertEqual(identity['status'], 'assigned')
        self.assertNotEqual(identity['id'], mirror['id'])
        for observations in ((identity, mirror), ({}, {})):
            with self.subTest(observations=observations):
                reactions = []
                for index, product in enumerate(observations):
                    ts = copy.deepcopy(species); ts.energy, ts.zpe = -3., .01
                    validation = compatible_l0_validation()
                    validation['endpoints'][1]['configured_identity'] = product
                    reactions.append(SimpleNamespace(bimol_search=True, family_name='abstraction',
                        instance_name=f'configured_{index}', instance=[0, 3, 2],
                        ts=ts, l0_validation=validation))
                species.reac_obj, species.reac_ts_done = reactions, [6, 6]
                # Even indistinguishable TS coordinates cannot overrule endpoints.
                with patch('kinbot.bimol_network.symmetry_rmsd', return_value=0.):
                    bimol.release_final_ts_clusters(species, {'multi_conf_tst': 0,
                        'bimol_energy_tolerance': .25, 'bimol_ts_rmsd': .1})
                self.assertEqual(species.reac_ts_done, [7, 7])
                self.assertFalse(any(hasattr(r, 'rc_basin_id') for r in reactions))

    def test_multi_conf_tst_keeps_distinct_channel_conformers(self):
        species = combined(self.hydroxyl, self.methane)
        points = [copy.deepcopy(species), copy.deepcopy(species)]
        for index, point in enumerate(points):
            point.energy = -3. - index * 1.e-4
            point.zpe = .01
        reactions = [SimpleNamespace(
            bimol_search=True, family_name='abstraction',
            instance_name=f'candidate_{index}', instance=[0, 3, 2],
            ts=point) for index, point in enumerate(points)]
        species.reac_obj = reactions
        species.reac_ts_done = [6, 6]
        parameters = {'bimol_energy_tolerance': .25,
                      'bimol_ts_rmsd': .1,
                      'multi_conf_tst': 1}
        with patch('kinbot.bimol_network.symmetry_rmsd', return_value=1.):
            bimol.release_final_ts_clusters(species, parameters)
        self.assertEqual(species.reac_ts_done, [7, 7])

    def test_distinct_ring_geometry_survives_single_conformer_selection(self):
        point = StationaryPoint('ring', 0, 1, smiles='C1CCCCC1')
        point.characterize()
        other = copy.copy(point)
        other.geom = point.geom.copy()
        other.geom[0, 2] += 1.
        for item in (point, other):
            item.energy, item.zpe = -3., .01
        species = SimpleNamespace(reac_obj=[SimpleNamespace(
            bimol_search=True, family_name='abstraction', instance_name=f'ring_{i}',
            instance=[0, 1, 2], l0_validation=compatible_l0_validation(),
            ts=item) for i, item in enumerate((point, other))],
            reac_ts_done=[6, 6])
        par = {'multi_conf_tst': 0, 'bimol_energy_tolerance': .25, 'bimol_ts_rmsd': .1}
        bimol.release_final_ts_clusters(species, par)
        archive = species.reac_obj[0].final_conformer_archive
        self.assertTrue(all(item['contains_graph_cycle'] for item in archive))
        self.assertEqual(archive[1]['role'], 'archived conformer')
        self.assertEqual(archive[1]['geometry_angstrom'], other.geom.tolist())
        self.assertIsNone(archive[1]['equilibration_established'])
        species.reac_ts_done[0] = 2
        bimol.release_final_ts_clusters(species, par)
        self.assertEqual(species.reac_ts_done, [2, -998])
        self.assertEqual(len(species.reac_obj[0].final_conformer_archive), 2)

    def test_single_conf_tst_preserves_distinct_observed_products(self):
        species = combined(self.hydroxyl, self.methane)
        points = [copy.deepcopy(species), copy.deepcopy(species)]
        for point in points:
            point.energy = -3.
            point.zpe = .01
        reactions = []
        for index, point in enumerate(points):
            reactions.append(SimpleNamespace(
                bimol_search=True, family_name='abstraction',
                instance_name=f'candidate_{index}', instance=[0, 3, 2],
                l0_validation={'endpoints': [
                    {'state': 'reactant', 'bonds': [[0, 1]]},
                    {'state': 'product',
                     'bonds': [[0, 3], [1, 2 + index]]},
                ]}, ts=point))
        species.reac_obj = reactions
        species.reac_ts_done = [6, 6]
        parameters = {'bimol_energy_tolerance': .25,
                      'bimol_ts_rmsd': .1,
                      'multi_conf_tst': 0}
        with patch('kinbot.bimol_network.symmetry_rmsd', return_value=0.):
            bimol.release_final_ts_clusters(species, parameters)
        self.assertEqual(species.reac_ts_done, [7, 7])

    def test_failed_final_irc_promotes_next_unique_cluster_member(self):
        species = combined(self.hydroxyl, self.methane)
        points = [copy.deepcopy(species), copy.deepcopy(species)]
        points[1].geom[0, 0] += .002
        for index, point in enumerate(points):
            point.energy = -3. + index * 1.e-6
            point.zpe = .01
        reactions = [SimpleNamespace(
            bimol_search=True, family_name='abstraction',
            instance_name=f'candidate_{index}', instance=[0, 3, 2],
            l0_validation=compatible_l0_validation(), ts=point) for index, point in enumerate(points)]
        species.reac_obj = reactions
        species.reac_ts_done = [6, 6]
        parameters = {'bimol_energy_tolerance': .25,
                      'bimol_ts_rmsd': .1}
        bimol.release_final_ts_clusters(species, parameters)
        failed = species.reac_ts_done.index(7)
        waiting = species.reac_ts_done.index(8)
        species.reac_ts_done[failed] = -999
        reactions[failed].bimol_final_status = 'final_irc_failed'
        bimol.release_final_ts_clusters(species, parameters)
        self.assertEqual(species.reac_ts_done[waiting], 7)
        self.assertTrue(reactions[waiting].final_irc_attempted)

    def test_network_uses_bimolecular_reactant_and_deduplicates_ts(self):
        reactants = combined(self.hydroxyl, self.methane)
        reactants.fragA.energy = -1.0
        reactants.fragA.zpe = 0.01
        reactants.fragB.energy = -2.0
        reactants.fragB.zpe = 0.02
        reactants.energy = -3.0
        reactants.zpe = 0.03

        water = stationary(molecule('H2O'), 1)
        methyl = stationary(molecule('CH3'), 2)
        water.energy, water.zpe = -1.2, 0.01
        methyl.energy, methyl.zpe = -2.1, 0.02
        rc = combined(self.hydroxyl, self.methane)
        pc = combined(water, methyl)
        ts = combined(self.hydroxyl, self.methane)
        rng = np.random.default_rng(7)
        for point, energy in ((rc, -3.01), (pc, -3.31), (ts, -3.005)):
            point.geom += rng.normal(0., .13, point.geom.shape)
            point.energy = energy
            point.zpe = 0.03
            point.reduced_freqs = [-500., 100.] if point is ts else [50.]
        reflected = [copy.deepcopy(point) for point in (rc, pc, ts)]
        for point in reflected:
            point.geom *= [-1., 1., 1.]

        reactions = []
        for index in range(2):
            reaction_rc, reaction_pc, reaction_ts = (
                (rc, pc, ts) if index == 0 else reflected)
            reaction = SimpleNamespace(
                instance_name=f'abstraction_{index}', family_name='abstraction',
                instance=[0, 3, 2], do_reac_vdW=True, do_vdW=True,
                irc_reactant=reaction_rc, irc_prod=reaction_pc,
                irc_reactant_direction='F', irc_prod_direction='R',
                ts=reaction_ts,
                irc_reactant_opt=SimpleNamespace(species=reaction_rc),
                irc_prod_opt=SimpleNamespace(species=reaction_pc),
                prod_opt=[SimpleNamespace(species=water),
                          SimpleNamespace(species=methyl)])
            reactions.append(reaction)
        reactants.reac_obj = reactions
        reactants.reac_ts_done = [-1, -1]
        parameters = {'bimol_energy_tolerance': 0.25,
                      'bimol_basin_rmsd': 0.15,
                      'bimol_ts_rmsd': 0.10}
        network = BimolecularNetwork(reactants, parameters)

        self.assertEqual(sum(node['kind'] == 'bimolecular'
                             for node in network.nodes.values()), 2)
        self.assertEqual(sum(node['kind'] == 'complex'
                             for node in network.nodes.values()), 2)
        self.assertEqual(len(network.transitions), 1)
        self.assertEqual(network.transitions[0]['degeneracy'], 4)

        self.assertEqual(len(network.transitions[0]['members']), 2)
        reactant = network.nodes[reactions[0].reactant_asymptote_id]
        self.assertEqual(reactant['kind'], 'bimolecular')
        self.assertAlmostEqual(reactant['energy'], 0.)

        with TemporaryDirectory() as directory:
            old = Path.cwd()
            try:
                os.chdir(directory)
                network.write_json()
                network.write_pesviewer()
                Path('bimol/unique_ts/stale').mkdir(parents=True)
                network.write_unique_ts()
            finally:
                os.chdir(old)
            text = Path(directory, 'pesviewer.inp').read_text()
            self.assertIn(reactions[0].reactant_asymptote_id, text)
            self.assertNotIn(f'\n{reactants.chemid} 0.0', text)
            unique = json.loads(Path(
                directory, 'bimol/unique_ts/index.json').read_text())
            self.assertEqual(len(unique['unique_ts']), 1)
            self.assertEqual(unique['unique_ts'][0]['level'], 'L1')
            self.assertEqual(
                unique['unique_ts'][0]['endpoints']['reactant'][
                    'irc_direction'], 'F')
            self.assertEqual(
                unique['unique_ts'][0]['endpoints']['product'][
                    'irc_direction'], 'R')
            self.assertEqual(len(list(Path(
                directory, 'bimol/unique_ts').glob('*/ts.xyz'))), 1)
            self.assertFalse(Path(
                directory, 'bimol/unique_ts/stale').exists())
            original_index = Path(
                directory, 'bimol/unique_ts/index.json').read_text()
            old = Path.cwd()
            try:
                os.chdir(directory)
                with patch('kinbot.bimol_network.write',
                           side_effect=RuntimeError('staging failure')):
                    with self.assertRaises(RuntimeError):
                        network.write_unique_ts()
            finally:
                os.chdir(old)
            self.assertEqual(Path(
                directory, 'bimol/unique_ts/index.json').read_text(),
                original_index)

        reactions[1].instance = [0, 4, 2]
        with patch('kinbot.bimol_network.symmetry_rmsd', return_value=0.):
            network = BimolecularNetwork(reactants, parameters)
        self.assertEqual(network.transitions[0]['degeneracy'], 4)

        reactions[0].search_instances = [[0, 3, 2], [0, 4, 2]]
        reactants.reac_obj = [reactions[0]]
        reactants.reac_ts_done = [-1]
        network = BimolecularNetwork(reactants, parameters)
        self.assertEqual(network.transitions[0]['degeneracy'], 4)
        self.assertEqual(len(network.transitions[0]['site_instances']), 2)

    def test_unexpected_single_product_is_exported_as_a_well(self):
        reactants = combined(self.hydroxyl, self.methane)
        reactants.fragA.energy = -1.
        reactants.fragA.zpe = .01
        reactants.fragB.energy = -2.
        reactants.fragB.zpe = .02
        reactants.energy = -3.
        reactants.zpe = .03
        product = copy.deepcopy(reactants)
        product.energy = -3.1
        product.zpe = .03
        ts = copy.deepcopy(reactants)
        ts.energy = -2.95
        ts.zpe = .03
        reaction = SimpleNamespace(
            instance_name='unexpected_association',
            family_name='abstraction', instance=[0, 3, 2],
            do_reac_vdW=False, do_vdW=False,
            irc_reactant=reactants, irc_prod=product, ts=ts,
            prod_opt=[SimpleNamespace(species=product)])
        reactants.reac_obj = [reaction]
        reactants.reac_ts_done = [-1]
        network = BimolecularNetwork(reactants, {
            'bimol_energy_tolerance': .25,
            'bimol_basin_rmsd': .15,
            'bimol_ts_rmsd': .1,
        })
        node = network.nodes[reaction.product_asymptote_id]
        self.assertEqual(node['kind'], 'well')
        serialized = next(item for item in network.as_dict()['nodes']
                          if item['id'] == node['id'])
        self.assertIn('stationary_point', serialized)
        self.assertNotIn('fragment_data', serialized)

    def test_validation_runs_once_after_unique_candidate_selection(self):
        species = combined(self.hydroxyl, self.methane)
        parameters = {'scan_step': 20, 'bimol_ts_rmsd': .1,
                      'bimol_energy_tolerance': .25,
                      'bimol_search_kwargs': {'sequential': False}}

        class Row:
            def __init__(self, data):
                self.data = data

        class DB:
            validation_queries = []

            def select(self, name):
                if name.endswith('_bimol_irc'):
                    self.validation_queries.append(name)
                    return iter([Row({'validation': {
                        'valid': True,
                        'endpoints': [
                            {'direction': 'forward', 'state': 'reactant'},
                            {'direction': 'reverse', 'state': 'product'},
                        ],
                    }})])
                candidates = [
                    {'geometry': species.geom.tolist(),
                     'directory': 'source/first',
                     'result': {'energy': -10.}},
                    {'geometry': (species.geom + .1).tolist(),
                     'directory': 'source/second',
                     'result': {'energy': -10.}},
                ]
                return iter([Row({'search': {'candidates': candidates}})])

        class QC:
            db = DB()

            def check_qc(self, name):
                return 'normal'

        qc = QC()
        name = 'bimol_abstraction'
        reaction = Abstraction(species, qc, parameters, [0, 3, 2], name)
        species.reac_type = ['abstraction']
        species.reac_inst = [[0, 3, 2]]
        species.reac_step = [0]
        species.reac_scan_energy = [[]]
        species.reac_ts_done = [0]
        species.reac_ts_geom = [0]
        species.reac_ts_freq = [0]
        species.reac_name = [name]
        species.reac_obj = [reaction]
        with TemporaryDirectory() as directory:
            old = Path.cwd()
            try:
                os.chdir(directory)
                bimol.prepare_bimolecular_reactions(
                    species, parameters, qc)
                unique = json.loads(Path(
                    'bimol/validated_l0/index.json').read_text())
                self.assertFalse(Path('bimol/unique_ts').exists())
            finally:
                os.chdir(old)
        self.assertEqual(len(species.reac_obj), 1)
        self.assertTrue(species.reac_name[0].startswith(f'{name}_l0_ts_'))
        np.testing.assert_allclose(species.reac_obj[0].guess_geom,
                                   species.geom)
        self.assertEqual(len(unique['validated_l0']), 1)
        self.assertEqual(len(unique['validated_l0'][0]['members']), 2)
        self.assertEqual(len(set(qc.db.validation_queries)), 1)


if __name__ == '__main__':
    unittest.main()
