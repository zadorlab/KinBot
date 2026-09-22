"""Stereochemical H-transfer classes and their statistical boundaries.

Geometries below are deterministic structural fixtures, not optimized TSs.
The supplied endpoint bond graphs isolate the fleeting-stereo comparison.
"""
import copy
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from kinbot import symmetry
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity, configured_geometry_allowed, optical_scope
from kinbot.reaction_path import (prepare_stereopath, path_geometry_allowed,
    reaction_path_id, same_path_class, summary_path_line, read_summary_paths,
    endpoint_snapshot, StereoAssignmentUnavailable)
from kinbot.conformer_counting import evaluate_members
from kinbot.conformer_records import ConformerRecord
from kinbot.hindered_rotors import HIR
from kinbot.optimize import Optimize
from kinbot.parameters import Parameters
from kinbot.reaction_finder import ReactionFinder
from kinbot.species_routing import input_species


def peroxy(smiles='CC[C@H](C)O[O]'):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=312) == 0
    p = StationaryPoint('peroxy', 0, 2,
        atom=[a.GetSymbol() for a in mol.GetAtoms()], geom=mol.GetConformer().GetPositions())
    p.characterize()
    p.energy, p.zpe = -100., .1
    p.freq = p.reduced_freqs = [500.] * (3*p.natom-6)
    symmetry.calculate_symmetry(p)
    return p


def transfer(p, hydrogen, donor=1, acceptor=5):
    q = copy.deepcopy(p)
    q.bond[donor, hydrogen] = q.bond[hydrogen, donor] = 0
    q.bond[acceptor, hydrogen] = q.bond[hydrogen, acceptor] = 1
    q.bonds = [q.bond.copy()]
    q.rads = [np.zeros(p.natom, int)]
    q.rads[0][donor] = 1
    q.calc_chemid()
    q.find_atom_eqv()
    ts = copy.deepcopy(p)
    ts.name = f'peroxy_transfer_{hydrogen}'
    ts.wellorts = 1
    ts.optical_reference = canonical_identity(p)
    ts.bond = np.maximum(p.bond, q.bond)
    ts.bonds = [ts.bond.copy()]
    ts.calc_chemid()
    ts.find_atom_eqv()
    ts.find_cycle()
    ts.find_conf_dihedral()
    symmetry.calculate_symmetry(ts)
    ts.freq = ts.reduced_freqs = [-1000.] + [500.] * (3*p.natom-7)
    prepare_stereopath(ts, p, q)
    return ts, q


class TestReactionPaths(unittest.TestCase):
    def setUp(self):
        self.p = peroxy()
        self.hydrogens = [i for i in range(self.p.natom)
                          if self.p.atom[i] == 'H' and self.p.bond[1, i]]

    def test_two_diastereotopic_transfers_have_different_direction_independent_classes(self):
        keys = []
        for hydrogen in self.hydrogens:
            ts, q = transfer(self.p, hydrogen)
            keys.append(ts.stereopath_id)
            reverse = copy.deepcopy(ts)
            prepare_stereopath(reverse, q, self.p)
            self.assertEqual(reverse.stereopath_id, ts.stereopath_id)
        self.assertEqual(len(set(keys)), 2)

    def test_mirror_reaching_another_specified_product_is_not_in_this_channel(self):
        ts, _ = transfer(self.p, self.hydrogens[0])
        # Isolate population logic: achiral entrance, configured chiral exit.
        ts.optical_reference = canonical_identity(peroxy('CCCO[O]'))
        self.assertFalse(optical_scope(ts)['mirror_allowed'])
        self.assertTrue(optical_scope(ts, 'racemic')['mirror_allowed'])

    def test_atom_order_and_global_reflection_do_not_make_new_path_classes(self):
        ts, q = transfer(self.p, self.hydrogens[0])
        original = ts.stereopath_id
        order = np.arange(self.p.natom)[::-1]
        inverse = {int(old): new for new, old in enumerate(order)}
        p = copy.deepcopy(self.p)
        for item in (p, q, ts):
            item.atom = np.asarray(item.atom)[order]
            item.geom = item.geom[order] * [-1., 1., 1.]
            item.bond = item.bond[np.ix_(order, order)]
            item.bonds = [matrix[np.ix_(order, order)] for matrix in item.bonds]
            item.rads = [np.asarray(rad)[order] for rad in item.rads]
            item.atom_eqv = [[inverse[i] for i in group] for group in item.atom_eqv]
        prepare_stereopath(ts, p, q)
        self.assertEqual(ts.stereopath_id, original)

    def test_conformer_cannot_move_to_other_route_or_unrequested_mirror(self):
        ts, _ = transfer(self.p, self.hydrogens[0])
        swapped = ts.geom.copy()
        swapped[self.hydrogens] = swapped[self.hydrogens[::-1]]
        self.assertFalse(configured_geometry_allowed(ts, swapped))
        mirror = ts.geom * [-1., 1., 1.]
        self.assertFalse(path_geometry_allowed(ts, mirror))
        self.assertTrue(path_geometry_allowed(ts, mirror, 'racemic'))
        records = [ConformerRecord(str(i), i, source_job=None, geometry=geom.tolist(),
                     zero_energy_hartree=-99.9, frequencies_cm1=tuple(ts.freq), status='valid')
                   for i, geom in enumerate((ts.geom, swapped))]
        observed, groups = evaluate_members(ts, records)
        self.assertEqual(groups, [[0]])
        self.assertEqual(observed[1].exclusion_reason, 'different stereochemical pathway')
        self.assertEqual(observed[0].remaining_optical_weight, 1.)

    def test_l2_geometry_and_mode_similarity_do_not_admit_another_pathway(self):
        ts, _ = transfer(self.p, self.hydrogens[0])
        swapped = ts.geom.copy()
        swapped[self.hydrogens] = swapped[self.hydrogens[::-1]]
        qc = SimpleNamespace(qc='fc', use_sella=True,
            get_qc_geom=lambda *a, **kw: (0, swapped),
            get_qc_freq=lambda *a: (0, ts.freq))
        with TemporaryDirectory() as directory:
            filename = Path(directory) / 'input.json'
            filename.write_text(json.dumps({'barrier_threshold': 100., 'conformer_search': 0}))
            par = Parameters(str(filename)).par
        opt = Optimize(ts, par, qc)
        mode = np.arange(ts.natom*3).reshape(-1, 3)
        with patch('kinbot.optimize.geometry.equal_geom', return_value=True), \
             patch('kinbot.optimize.reader_sella.read_imag_mode', return_value=mode):
            opt.compare_structures()
        self.assertEqual(opt.shigh, -999)
        np.testing.assert_array_equal(ts.geom, self.p.geom)

    def test_hir_must_not_interpolate_over_another_explicit_route(self):
        ts, _ = transfer(self.p, self.hydrogens[0])
        swapped = ts.geom.copy()
        swapped[self.hydrogens] = swapped[self.hydrogens[::-1]]
        ts.dihed = [[6, 0, 1, 2]]
        qc = SimpleNamespace(get_qc_geom=lambda job, *a: (0, swapped if job.endswith('_03') else ts.geom),
                             get_qc_energy=lambda job: (0, -100.))
        hir = HIR(ts, qc, {'nrotation': 12, 'plot_hir_profiles': 0, 'rotor_0_test': 0})
        hir.hir_status = [[-1]*12]
        hir.hir_energies = [[0.]*12]
        hir.hir_geoms = [[None]*12]
        with patch('kinbot.hindered_rotors.geometry.equal_geom', return_value=True):
            hir.test_hir()
        self.assertEqual(hir.hir_status[0], [0, 0, 0, 1]+[0]*8)
        np.testing.assert_array_equal(hir.point_observations[0][3]['geometry_angstrom'], swapped)
        self.assertFalse(hir.is_valid_rotor(0))
        self.assertIn('different stereochemical pathway', hir.invalid_rotor_reason(0))

    def test_summary_round_trip_and_no_unknown_route_double_counting(self):
        ts, q = transfer(self.p, self.hydrogens[0])
        reaction = SimpleNamespace(ts=ts, instance_name=ts.name)
        line = summary_path_line(reaction)
        paths = read_summary_paths(['SUCCESS 30.0 test product', line])
        self.assertEqual(paths[ts.name], reaction_path_id(reaction))
        self.assertTrue(same_path_class(None, None))
        self.assertTrue(same_path_class(ts.stereopath_id, ts.stereopath_id))
        self.assertFalse(same_path_class('a', 'b'))
        self.assertFalse(same_path_class(ts.stereopath_id, None))
        ts.geom[self.hydrogens] = ts.geom[self.hydrogens[::-1]]
        with self.assertRaisesRegex(ValueError, 'changed stereochemical pathway'):
            summary_path_line(reaction)

    def test_ordinary_enantiotopic_transfer_keeps_original_selection_policy(self):
        p = peroxy('CCCO[O]')
        hydrogens = [i for i in range(p.natom) if p.atom[i] == 'H' and p.bond[1, i]]
        for hydrogen in hydrogens:
            ts, _ = transfer(p, hydrogen, acceptor=4)
            self.assertEqual(ts.stereopath_id, 'ordinary')
            line = summary_path_line(SimpleNamespace(ts=ts, instance_name=ts.name))
            self.assertEqual(read_summary_paths([line])[ts.name], 'ordinary')

    def test_ordinary_transfer_still_works_without_optional_rdkit(self):
        p = peroxy('CCCO[O]')
        hydrogen = next(i for i in range(p.natom) if p.atom[i] == 'H' and p.bond[1, i])
        ts, q = transfer(p, hydrogen, acceptor=4)
        with patch('kinbot.stereo_identity._strings', side_effect=ImportError('no RDKit')):
            prepare_stereopath(ts, p, endpoint_snapshot(q))
            line = summary_path_line(SimpleNamespace(ts=ts, instance_name=ts.name))
            self.assertEqual(read_summary_paths([line])[ts.name], 'ordinary')

    def test_demonstrated_split_requires_canonical_assignment(self):
        ts, q = transfer(self.p, self.hydrogens[0])
        with patch('kinbot.stereo_identity._strings', side_effect=ImportError('no RDKit')):
            with self.assertRaisesRegex(StereoAssignmentUnavailable, 'canonical stereo assignment'):
                prepare_stereopath(ts, self.p, endpoint_snapshot(q))

    def test_runnable_inputs_enumerate_two_diastereotopic_routes_and_one_control(self):
        directory = Path(__file__).resolve().parents[1] / 'examples/stereochemical_h_transfer'
        for parent, path, expected in [('secbutylperoxy', [5, 4, 2, 1], {9, 10}),
                                       ('propylperoxy_control', [4, 3, 2, 1], {8})]:
            for treatment, mc in [('hir', 0), ('mc_rrho', 1)]:
                filename = directory / f'{parent}_{treatment}.json'
                par = Parameters(str(filename), show_warnings=False).par
                p = input_species(filename)
                ReactionFinder(p, par, None).find_reactions()
                self.assertEqual({int(r.instance[-1]) for r in p.reac_obj
                                  if list(r.instance[:-1]) == path}, expected)
                self.assertEqual(par['multi_conf_tst'], mc)
                self.assertEqual(par['rotor_scan'], 1-mc)
                self.assertEqual(par['me'], 2)


if __name__ == '__main__':
    unittest.main()
