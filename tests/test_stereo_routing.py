"""Stereo conflicts cannot be hidden by legacy chemid aliases or disk caches."""
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
from kinbot.stereo_identity import canonical_identity
from kinbot.stereo_routing import (guard_well_job, guard_pes_input,
    require_same_configuration, StereoRoutingError, _row_species, preserve_observations)
from kinbot.stationary_pt import StationaryPoint
from kinbot.species_routing import routing_name, prepare_qc_routing, resolve_job
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

    def test_preexisting_opposite_stereo_cache_is_checked_without_memory_alias(self):
        db = connect('kinbot.db')
        job = f'{self.first.chemid}_well'
        db.write(Atoms(self.second.atom, positions=self.second.geom), name=job,
                 data={'status': 'normal', 'energy': -100., 'zpe': .01,
                       'charge': self.second.charge, 'multiplicity': self.second.mult})
        qc = QuantumChemistry.__new__(QuantumChemistry)
        qc.db, qc.submit_qc = db, Mock()
        prepare_qc_routing(qc, self.first)
        configured = routing_name(self.first) + '_well'
        self.assertNotEqual(configured, job)
        self.assertEqual(resolve_job(db, configured), configured)
        qc.submit_qc.assert_not_called()
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

    def test_no_rdkit_preserves_explicitly_unverified_legacy_alias_behavior(self):
        with patch('kinbot.stereo_routing.canonical_identity', return_value={'status': 'unavailable'}):
            require_same_configuration(self.first, self.second, 'legacy')
        self.assertIn('unverified', self.first.stereo_routing_status)

    def test_fresh_unsupported_ordinary_jobs_and_identical_polls_remain_available(self):
        db = connect('kinbot.db')
        qc = SimpleNamespace(db=db, par={'multi_conf_tst': 0})
        # Substituted PAHs remain outside the narrow anthracene/phenanthrene class.
        for smiles in ('FC=C=CF', 'Cc1ccc2cc3ccccc3cc2c1'):
            point = StationaryPoint('unsupported', 0, 1, smiles=smiles)
            point.characterize()
            self.assertEqual(canonical_identity(point)['status'], 'unsupported')
            job = f'{point.chemid}_well'
            guard_well_job(qc, point, point.geom, job)
            # Normal late optical bookkeeping does not change a QC request.
            point.optical_reference = canonical_identity(point)
            point.optical_population = 'specified'
            db.write(Atoms(point.atom, positions=point.geom), name=job, data={'status': 'normal'})
            guard_well_job(qc, point, point.geom, job)
            self.assertIn('unverified', point.stereo_routing_status)
            self.assertEqual(len(list(db.select(name=f'stereochemistry/{job}'))), 1)
            guard_well_job(SimpleNamespace(db=db, par={'multi_conf_tst': 1}), point, point.geom, job)
            for change in ('coordinates', 'isotopes', 'charge'):
                other = copy.copy(point)
                other.geom = point.geom.copy()
                if change == 'coordinates':
                    other.geom[0, 0] += .01
                elif change == 'isotopes':
                    other.isotopes = [13] + [0] * (other.natom-1)
                else:
                    other.charge = 1
                with self.subTest(smiles=smiles, change=change):
                    if change == 'coordinates':
                        guard_well_job(qc, other, other.geom, job)
                    else:
                        with self.assertRaises(StereoRoutingError):
                            guard_well_job(qc, other, other.geom, job)
            db.write(Atoms(point.atom, positions=point.geom), name=job+'_old', data={'status': 'normal'})
            guard_well_job(qc, point, point.geom, job+'_old')
            self.assertEqual(len(list(db.select(name=f'stereochemistry/{job}_old'))), 1)
            Path(job+'_pending.py').write_text('# existing job without an input reference\n')
            guard_well_job(qc, point, point.geom, job+'_pending')


if __name__ == '__main__':
    unittest.main()
