"""Ordinary PAHs and biaryls use the shared stereo and counting workflow."""
import copy
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import numpy as np
from ase import Atoms
from ase.db import connect
from rdkit import Chem
from rdkit.Chem import AllChem
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity
from kinbot.species_routing import routing_key, same_species
from kinbot.stereo_routing import guard_well_job
from kinbot.conformer_counting import writer_members
from kinbot.calculation import load_calculation_record
from kinbot import constants
from tests.conformer_fixtures import record_conformers


ANTHRACENE = 'c1ccc2cc3ccccc3cc2c1'
PHENANTHRENE = 'c1ccc2c(c1)ccc1ccccc12'
# PubChem CID 98863: https://pubchem.ncbi.nlm.nih.gov/compound/Hexahelicene
HELICENE = 'C1=CC=C2C(=C1)C=CC3=C2C4=C(C=C3)C=CC5=C4C6=CC=CC=C6C=C5'
PYRENE = 'c1cc2ccc3cccc4ccc(c1)c2c34'


def molecule(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 17
    parameters.useRandomCoords = True
    parameters.useBasicKnowledge = False  # allow the helicene ring framework to twist
    if AllChem.EmbedMolecule(mol, parameters) != 0:
        raise ValueError('Failed to embed fixture')
    if AllChem.UFFOptimizeMolecule(mol, maxIters=1000) != 0:
        raise ValueError('Failed to optimize fixture')
    Chem.Kekulize(mol, clearAromaticFlags=True)
    return mol


def observation(mol):
    return SimpleNamespace(atom=[a.GetSymbol() for a in mol.GetAtoms()],
        geom=mol.GetConformer().GetPositions(), bond=Chem.GetAdjacencyMatrix(mol, useBO=True),
        charge=0, mult=1)


class TestPAHStereoScope(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.points = []
        for smiles in (ANTHRACENE, PHENANTHRENE, PYRENE, 'C'+ANTHRACENE,
                       'c1ccc2c(c1)Cc1ccccc1-2', '[c]1ccc2cc3ccccc3cc2c1',
                       'c1ccccc1-c1ccccc1'):
            obs = observation(molecule(smiles))
            point = StationaryPoint('PAH', 0, 2 if smiles.startswith('[c]') else 1,
                                    atom=obs.atom, geom=obs.geom)
            point.characterize()
            cls.points.append(point)

    def test_supported_skeletons_keep_ordinary_names_and_member_weights(self):
        for original in self.points:
            p = copy.copy(original)
            self.assertEqual(canonical_identity(p)['status'], 'assigned')
            self.assertEqual(routing_key(p), p.chemid)
            other = StationaryPoint('independent', 0, p.mult, atom=p.atom, geom=p.geom.copy())
            other.characterize()
            self.assertTrue(same_species(p, other))
            p.conformer_index = [7]
            p.conformer_geom = [p.geom.copy()]
            p.conformer_zeroenergy = [-1.]
            p.conformer_freq = [[100.] * (3*p.natom-6)]
            record_conformers(p)
            members = writer_members(p)
            self.assertEqual(members[0].index, 7)
            self.assertIsNotNone(members[0].stereo_identity)
            self.assertIn(members[0].remaining_optical_weight, (1., 2.))

    def test_bending_rigid_motion_and_permutation_do_not_change_support_or_key(self):
        for original in self.points:
            expected = canonical_identity(original)['id']
            p = copy.copy(original)
            centered = p.geom - p.geom.mean(axis=0)
            _, _, axes = np.linalg.svd(centered, full_matrices=False)
            # Deliberately exceed the rejected 0.05 A planarity proposal.
            p.geom = p.geom + .25*np.sin(centered @ axes[0])[:, None]*axes[-1]
            p.geom = p.geom @ np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]) + 4.
            self.assertEqual(canonical_identity(p)['id'], expected)
            order = np.arange(p.natom)[::-1]
            p.atom, p.geom = np.asarray(p.atom)[order], p.geom[order]
            p.bond = p.bond[np.ix_(order, order)]
            p.bonds = [matrix[np.ix_(order, order)] for matrix in p.bonds]
            p.rads = [np.asarray(rad)[order] for rad in p.rads]
            self.assertEqual(canonical_identity(p)['id'], expected)
            self.assertEqual(routing_key(p), original.chemid)

    def test_current_well_results_remain_reusable(self):
        with TemporaryDirectory() as directory:
            qc = SimpleNamespace(db=connect(str(Path(directory)/'kinbot.db')), par={}, job_ids={})
            for original in self.points:
                p = copy.copy(original)
                job = f'{p.chemid}_well'
                freq = [100.] * (3*p.natom-6)
                qc.db.write(Atoms(p.atom, positions=p.geom), name=job,
                    data={'status': 'normal', 'energy': -1./constants.EVtoHARTREE,
                          'zpe': .01, 'frequencies': freq, 'charge': 0, 'multiplicity': p.mult})
                guard_well_job(qc, p, p.geom, job)
                load_calculation_record(p, qc, job)
                self.assertEqual(p.source_job, job)
                self.assertEqual(p.energy, -1.)
                self.assertEqual(len(list(qc.db.select(name=job))), 1)

    def test_assignment_preserves_virtual_and_real_isotope_labels(self):
        p = copy.copy(self.points[0])
        hydrogen = next(i for i, atom in enumerate(p.atom) if atom == 'H')
        self.assertEqual(canonical_identity(p, tagged_atom=hydrogen)['status'], 'assigned')
        p.isotopes = [0] * p.natom
        p.isotopes[hydrogen] = 2
        assigned = canonical_identity(p)
        self.assertEqual(assigned['status'], 'assigned')
        self.assertNotEqual(assigned['id'], canonical_identity(self.points[0])['id'])

    def test_virtual_labels_do_not_create_a_physical_biaryl_axis(self):
        from kinbot.stereochemistry import motif_identity, virtually_labelled
        p = self.points[-1]
        original = canonical_identity(p)
        labelled = virtually_labelled(p, {i: i for i in range(p.natom)})
        self.assertEqual(canonical_identity(labelled)['status'], 'assigned')
        self.assertEqual(motif_identity(labelled, list(range(p.natom)))['status'], 'assigned')
        self.assertEqual(canonical_identity(p, tagged_atom=0)['status'], 'assigned')
        self.assertEqual(canonical_identity(p), original)
        self.assertFalse(hasattr(p, 'isotopes'))

    def test_helicity_is_not_encoded_as_a_supported_configuration(self):
        # Removing the PAH ban does not add a helical stereoisomer identifier.
        from kinbot.optical import compare_rigid
        obs = observation(molecule(HELICENE))
        identity = canonical_identity(obs)
        self.assertEqual(identity['status'], 'assigned')
        self.assertEqual(identity['id'], identity['mirror_id'])
        self.assertEqual(compare_rigid(obs, obs.geom, obs.geom)['status'], 'distinct')

    def test_existing_biaryl_and_coordination_graph_guards_remain_in_force(self):
        for smiles, reason in [('Cc1cccc(F)c1-c1c(C)cccc1F', 'biaryl'),
                               ('FS(F)(F)(F)(F)F', 'coordination')]:
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            Chem.Kekulize(mol, clearAromaticFlags=True)
            # These graph-only exclusions precede any 3-D stereo assignment.
            mol.AddConformer(Chem.Conformer(mol.GetNumAtoms()))
            result = canonical_identity(observation(mol))
            self.assertEqual(result['status'], 'unsupported')
            self.assertIn(reason, result['reason'])


if __name__ == '__main__':
    unittest.main()
