"""A narrow ordinary PAH class must not trigger the unsupported-helicity guard."""
import copy
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import numpy as np
from ase import Atoms
from rdkit import Chem
from rdkit.Chem import AllChem
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity, _three_ring_benzenoid
from kinbot.species_routing import routing_key, same_species, connect
from kinbot.stereo_routing import guard_well_job
from kinbot.conformer_counting import writer_members
from kinbot.calculation import load_calculation_record
from kinbot import constants


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
        for smiles in (ANTHRACENE, PHENANTHRENE):
            obs = observation(molecule(smiles))
            point = StationaryPoint('PAH', 0, 1, atom=obs.atom, geom=obs.geom)
            point.characterize()
            cls.points.append(point)

    def test_supported_skeletons_keep_ordinary_names_and_member_weights(self):
        for original in self.points:
            p = copy.copy(original)
            self.assertEqual(canonical_identity(p)['status'], 'assigned')
            self.assertEqual(routing_key(p), p.chemid)
            other = StationaryPoint('independent', 0, 1, atom=p.atom, geom=p.geom.copy())
            other.characterize()
            self.assertTrue(same_species(p, other))
            p.conformer_index = [7]
            p.conformer_geom = [p.geom.copy()]
            p.conformer_zeroenergy = [-1.]
            p.conformer_freq = [[100.] * (3*p.natom-6)]
            members = writer_members(p)
            self.assertEqual(members[0].index, 7)
            self.assertIsNotNone(members[0].stereo_identity)
            self.assertEqual(members[0].remaining_optical_weight, 1.)

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

    def test_legacy_well_results_remain_reusable(self):
        with TemporaryDirectory() as directory:
            qc = SimpleNamespace(db=connect(str(Path(directory)/'kinbot.db')), par={}, job_ids={})
            for original in self.points:
                p = copy.copy(original)
                job = f'{p.chemid}_well'
                freq = [100.] * (3*p.natom-6)
                qc.db.write(Atoms(p.atom, positions=p.geom), name=job,
                    data={'status': 'normal', 'energy': -1./constants.EVtoHARTREE,
                          'zpe': .01, 'frequencies': freq, 'charge': 0, 'multiplicity': 1})
                guard_well_job(qc, p, p.geom, job)
                load_calculation_record(p, qc, job)
                self.assertEqual(p.source_job, job)
                self.assertEqual(p.energy, -1.)
                self.assertEqual(len(list(qc.db.select(name=job))), 1)

    def test_exemption_preserves_virtual_and_real_isotope_labels(self):
        p = copy.copy(self.points[0])
        hydrogen = next(i for i, atom in enumerate(p.atom) if atom == 'H')
        self.assertEqual(canonical_identity(p, tagged_atom=hydrogen)['status'], 'assigned')
        p.isotopes = [0] * p.natom
        p.isotopes[hydrogen] = 2
        assigned = canonical_identity(p)
        self.assertEqual(assigned['status'], 'assigned')
        self.assertNotEqual(assigned['id'], canonical_identity(self.points[0])['id'])

    def test_larger_substituted_saturated_and_disconnected_graphs_stay_outside_class(self):
        cases = [(PYRENE, 'unsupported'), ('C'+ANTHRACENE, 'unsupported'),
                 # Fewer pi rings already bypass the guard; they do not enter this exemption.
                 ('C1CCc2cc3ccccc3cc2C1', 'assigned'), (ANTHRACENE+'.C', 'unsupported')]
        for smiles, status in cases:
            with self.subTest(smiles=smiles):
                mol = molecule(smiles)
                self.assertFalse(_three_ring_benzenoid(mol))
                self.assertEqual(canonical_identity(observation(mol))['status'], status)
        # An extra bridge must not pass merely because all carbon atoms remain in the core.
        mol = Chem.RWMol(molecule(ANTHRACENE))
        mol.AddBond(0, 4, Chem.BondType.SINGLE)
        mol = mol.GetMol()
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)
        self.assertFalse(_three_ring_benzenoid(mol))

    def test_actual_helicene_and_mirror_remain_unsupported(self):
        obs = observation(molecule(HELICENE))
        carbon = obs.geom[np.asarray(obs.atom) == 'C']
        centered = carbon-carbon.mean(axis=0)
        _, _, axes = np.linalg.svd(centered, full_matrices=False)
        self.assertGreater(np.max(abs(centered @ axes[-1])), .2)
        for geom in (obs.geom, obs.geom * [-1., 1., 1.]):
            result = canonical_identity(obs, geom)
            self.assertEqual(result['status'], 'unsupported')
            self.assertIn('helical/planar', result['reason'])

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
