"""Current saved inputs and results must retain their requested stereoisomer."""
import copy
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch, Mock
import numpy as np
from ase import Atoms
from ase.db import connect
from kinbot.reaction_generator import ReactionGenerator
from kinbot.stereo_identity import canonical_identity, UnsupportedStereochemistry
from kinbot.stereo_routing import (guard_well_job, guard_pes_input,
    require_same_configuration, StereoRoutingError, _row_species, preserve_observations)
from kinbot.stationary_pt import StationaryPoint
from kinbot.species_routing import routing_name
from kinbot.qc import QuantumChemistry
from kinbot.conformer_counting import representative_record, CountingError


class TestStereoRouting(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.addCleanup(os.chdir, Path.cwd())
        os.chdir(temporary.name)
        self.first = StationaryPoint('R', 0, 1, smiles='C[C@H](O)CC')
        self.first.characterize()
        self.second = copy.copy(self.first)
        self.second.geom = self.first.geom * [-1, 1, 1]
        self.second.name = 'S'

    def test_supported_configurations_remain_independent_objects(self):
        first, second = self.first, self.second
        self.assertEqual(first.chemid, second.chemid)
        fragments = [first, second]
        ReactionGenerator.equate_identical(None, fragments)
        self.assertIs(fragments[1], second)
        unique = [first]
        ReactionGenerator.equate_unique(None, [second], unique)
        self.assertEqual(len(unique), 2)
        self.assertIs(unique[1], second)

    def test_opposite_stereo_cache_is_rejected(self):
        db = connect('kinbot.db')
        job = routing_name(self.first) + '_well'
        db.write(Atoms(self.second.atom, positions=self.second.geom), name=job,
                 data={'status': 'normal', 'energy': -100., 'zpe': .01,
                       'charge': self.second.charge, 'multiplicity': self.second.mult})
        qc = SimpleNamespace(db=db)
        with self.assertRaisesRegex(StereoRoutingError, 'different stereoisomer'):
            guard_well_job(qc, self.first, self.first.geom, job)
        self.assertEqual(len(list(db.select(name=job))), 1)
        self.assertEqual(len(list(db.select(name=f'stereochemistry/{job}'))), 0)

    def test_input_scope_survives_a_restart_before_the_first_result(self):
        db = connect('kinbot.db')
        job = f'{self.first.chemid}_well'
        guard_well_job(SimpleNamespace(db=db), self.first, self.first.geom, job)
        with self.assertRaisesRegex(StereoRoutingError, 'cached input'):
            guard_well_job(SimpleNamespace(db=connect('kinbot.db')), self.second, self.second.geom, job)

    def test_raw_labels_graph_and_optical_scope_survive_reference_and_refusal(self):
        point = self.first
        point.isotopes = [13] + [0] * (point.natom-1)
        point.formal_charges = [0] * point.natom
        point.optical_reference = canonical_identity(point)
        point.optical_population = 'specified'
        db = connect('kinbot.db')
        guard_well_job(SimpleNamespace(db=db), point, point.geom, 'labelled_well')
        row = db.get(name='stereochemistry/labelled_well')
        restored = _row_species(point, row, 'labelled_well')
        report = json.loads(Path(preserve_observations('fixture', [restored])).read_text())
        context = report['observations'][0]['chemical_context']
        self.assertEqual(context['isotopes'], point.isotopes)
        self.assertEqual(context['formal_charges'], point.formal_charges)
        self.assertEqual(context['optical_reference'], json.loads(json.dumps(point.optical_reference)))
        self.assertIsNone(context['optical_counting_scope'])
        for name in ('bond', 'bond01', 'bonds', 'rads'):
            np.testing.assert_array_equal(context[name], getattr(point, name))

    def test_labelled_completed_result_uses_only_its_trusted_input_labels(self):
        point = self.first
        point.isotopes = [13] + [0] * (point.natom-1)
        db = connect('kinbot.db')
        qc = SimpleNamespace(db=db)
        guard_well_job(qc, point, point.geom, 'labelled')
        db.write(Atoms(point.atom, positions=point.geom), name='labelled', data={'status': 'normal'})
        guard_well_job(qc, point, point.geom, 'labelled')
        # An old result without requested-input provenance cannot acquire the
        # isotope tags from the new caller just to make its identity match.
        db.write(Atoms(point.atom, positions=point.geom), name='old', data={'status': 'normal'})
        with self.assertRaises(StereoRoutingError):
            guard_well_job(qc, point, point.geom, 'old')

    def test_explicit_result_state_cannot_be_replaced_with_the_callers_charge_or_spin(self):
        point = StationaryPoint('oxide', -2, 1, atom=['O'], geom=np.zeros((1, 3)))
        point.characterize()
        db = connect('kinbot.db')
        qc = SimpleNamespace(db=db, par={'multi_conf_tst': 1})
        for index, data in enumerate(({'charge': 0, 'mult': 1}, {'charge': -2, 'multiplicity': 3})):
            job = f'state_{index}'
            # Test both old unreferenced caches and contradicting results after
            # a trusted request. Neither may replace the reported state.
            if index:
                guard_well_job(qc, point, point.geom, job)
            db.write(Atoms(point.atom, positions=point.geom), name=job, data={'status': 'normal', **data})
            with self.assertRaisesRegex(StereoRoutingError, 'charge or multiplicity'):
                guard_well_job(qc, point, point.geom, job)
        row = SimpleNamespace(data={'status': 'normal'}, symbols=['O'], positions=point.geom,
            calculator_parameters={'charge': 0, 'multiplicity': 3}, id=123)
        observed = _row_species(point, row, 'calculator state')
        self.assertEqual((observed.charge, observed.mult), (0, 3))
        self.assertEqual(observed.calculation_state_evidence['charge']['source'], 'calculator_parameters.charge')
        row.calculator_parameters = {}
        unknown = _row_species(point, row, 'unreported state')
        self.assertIn('unreported', unknown.calculation_state_evidence['charge']['source'])
        self.assertEqual(unknown.calculation_state_evidence['charge']['observations'], [])

    def test_compatible_cache_and_existing_connectivity_change_path_remain_available(self):
        db = connect('kinbot.db')
        job = f'{self.first.chemid}_well'
        db.write(Atoms(self.first.atom, positions=self.first.geom), name=job, data={'status': 'normal'})
        guard_well_job(SimpleNamespace(db=db), self.first, self.first.geom, job)
        guard_well_job(SimpleNamespace(db=db), self.first, self.first.geom, job)
        self.assertEqual(len(list(db.select(name=f'stereochemistry/{job}'))), 1)
        # Existing fragment dissociation handling is not a stereo-cache collision.
        separated = self.first.geom.copy()
        separated[-1] += 20.
        db.write(Atoms(self.first.atom, positions=separated), name=job, data={'status': 'normal'})
        guard_well_job(SimpleNamespace(db=db), self.first, self.first.geom, job)

    def test_pes_input_is_not_overwritten_with_another_configuration(self):
        path = Path('well.json')
        structure = [value for atom, xyz in zip(self.first.atom, self.first.geom)
                     for value in [str(atom), *map(float, xyz)]]
        text = json.dumps({'charge': 0, 'mult': 1, 'structure': structure})
        path.write_text(text)
        with self.assertRaises(StereoRoutingError):
            guard_pes_input(self.second, path)
        self.assertEqual(path.read_text(), text)
        guard_pes_input(self.first, path)

    def test_unsupported_well_is_rejected_before_writing_a_reference(self):
        db = connect('kinbot.db')
        qc = SimpleNamespace(db=db)
        for smiles in ('FC=C=CF', 'Cc1ccc2cc3ccccc3cc2c1'):
            point = StationaryPoint('unsupported', 0, 1, smiles=smiles)
            point.characterize()
            self.assertEqual(canonical_identity(point)['status'], 'unsupported')
            job = f'{point.chemid}_well'
            with self.assertRaises(UnsupportedStereochemistry):
                guard_well_job(qc, point, point.geom, job)
            self.assertEqual(len(list(db.select(name=f'stereochemistry/{job}'))), 0)
            point.optical_reference = canonical_identity(point)
            with self.assertRaises(UnsupportedStereochemistry):
                routing_name(point)
