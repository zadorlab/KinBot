import copy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
import numpy as np
from types import SimpleNamespace
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity, optical_scope
from kinbot.stereochemistry import refine_equivalence_group, stereotopic_hydrogen_equivalence
from kinbot.parameters import Parameters
from kinbot.reaction_finder import ReactionFinder
from kinbot import stereochemistry


def point(smiles):
    species = StationaryPoint('test', 0, 1, smiles=smiles)
    species.characterize()
    return species


class TestStereoIdentity(unittest.TestCase):
    def test_axial_chirality_is_not_declared_achiral_by_missing_rdkit_tags(self):
        species = SimpleNamespace(atom=['C', 'C', 'C', 'H', 'F', 'H', 'Cl'],
            geom=np.array([[-1.3, 0, 0], [0, 0, 0], [1.3, 0, 0],
                [-1.8, .9, 0], [-1.8, -.9, 0], [1.8, 0, .9], [1.8, 0, -.9]]),
            bond=np.zeros((7, 7), int), charge=0, mult=1)
        for i, j, order in [(0, 1, 2), (1, 2, 2), (0, 3, 1), (0, 4, 1), (2, 5, 1), (2, 6, 1)]:
            species.bond[i, j] = species.bond[j, i] = order
        self.assertEqual(canonical_identity(species)['status'], 'unsupported')
        self.assertIn('axial', canonical_identity(species)['reason'])
        self.assertFalse(optical_scope(species)['mirror_allowed'])

    def test_canonical_atom_permutation_and_mirror(self):
        species = point('C[C@H](O)CC')
        original = canonical_identity(species)
        mirror = canonical_identity(species, species.geom * [-1, 1, 1])
        self.assertEqual(original['status'], 'assigned')
        self.assertNotEqual(original['id'], mirror['id'])
        self.assertEqual(original['mirror_family_id'], mirror['mirror_family_id'])
        order = np.arange(species.natom)[::-1]
        other = copy.copy(species)
        other.atom = np.asarray(species.atom)[order]
        other.geom = species.geom[order]
        other.bond = species.bond[np.ix_(order, order)]
        other.bonds = [matrix[np.ix_(order, order)] for matrix in species.bonds]
        other.rads = [np.asarray(rad)[order] for rad in species.rads]
        self.assertEqual(canonical_identity(other)['id'], original['id'])

    def test_specified_stereoisomer_and_explicit_racemate(self):
        species = point('C[C@H](O)CC')
        self.assertFalse(optical_scope(species)['mirror_allowed'])
        self.assertTrue(optical_scope(species, 'racemic')['mirror_allowed'])
        self.assertTrue(optical_scope(point('CCO'))['mirror_allowed'])

    def test_targeted_hydrogen_refinement_preserves_global_equivalence(self):
        for smiles, expected in [('CCO', 1), ('C[C@H](O)CC', 2)]:
            species = point(smiles)
            original = copy.deepcopy(species.atom_eqv)
            group = next(g for g in original if len(g) == 2 and species.atom[g[0]] == 'H')
            self.assertEqual(len(refine_equivalence_group(species, group)), expected)
            stereotopic_hydrogen_equivalence(species)
            self.assertEqual(species.atom_eqv, original)

    def test_stereotopic_refinement_substitutes_each_hydrogen_only_once(self):
        species = point('CC')
        group = next(group for group in species.atom_eqv if len(group) == 6)
        with patch.object(stereochemistry, 'canonical_identity', wraps=canonical_identity) as assigned:
            classes = refine_equivalence_group(species, group)
        self.assertEqual(assigned.call_count, len(group))
        self.assertEqual(classes, [{'members': group, 'representative': group[0],
                                    'relation': 'homotopic'}])

    def test_single_site_needs_no_stereochemical_assignment(self):
        with patch.object(stereochemistry, 'canonical_identity') as assigned:
            self.assertEqual(refine_equivalence_group(None, []), [])
            self.assertEqual(refine_equivalence_group(None, [3]),
                             [{'members': [3], 'representative': 3, 'relation': 'homotopic'}])
        assigned.assert_not_called()

    def test_stereotopic_metadata_and_representatives_survive_shuffling_and_fallback(self):
        for smiles, size, relation in [('CC', 6, 'homotopic'),
                                      ('CCO', 2, 'enantiotopic'),
                                      ('C[C@H](O)CC', 2, 'diastereotopic')]:
            species = point(smiles)
            group = next(g for g in species.atom_eqv
                         if len(g) == size and species.atom[g[0]] == 'H')
            for order in (group, group[::-1]):
                members = [[atom] for atom in order] if relation == 'diastereotopic' else [order]
                expected = [{'members': atoms, 'representative': atoms[0], 'relation': relation}
                            for atoms in members]
                for status in ('assigned', 'unavailable', 'unsupported'):
                    with self.subTest(smiles=smiles, order=order, status=status):
                        assignment = ({'wraps': canonical_identity} if status == 'assigned'
                                      else {'return_value': {'status': status}})
                        with patch.object(stereochemistry, 'canonical_identity', **assignment):
                            self.assertEqual(refine_equivalence_group(species, order), expected)

    def test_diastereotopic_hydrogens_survive_actual_transfer_family_objects(self):
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / 'input.json'
            for family in ('intra_H_migration', 'ketoenol'):
                species = point('CC(=O)C[C@H](O)C')
                original = copy.deepcopy(species.atom_eqv)
                hydrogens = next(g for g in original if len(g) == 2 and species.atom[g[0]] == 'H')
                donor = next(i for i in range(species.natom) if species.bond[i][hydrogens[0]])
                carbonyl = next(i for i in range(species.natom) if species.atom[i] == 'C'
                               and any(species.atom[j] == 'O' and species.bond[i][j] == 2
                                       for j in range(species.natom)))
                oxygen = next(j for j in range(species.natom) if species.atom[j] == 'O'
                              and species.bond[carbonyl][j] == 2)
                path.write_text(json.dumps({'barrier_threshold': 50, 'families': [family]}))
                finder = ReactionFinder(species, Parameters(path).par, None)
                finder.find_reactions()
                observed = [reaction for reaction in species.reac_obj
                            if list(reaction.instance[:-1]) == [oxygen, carbonyl, donor]
                            and species.atom[reaction.instance[-1]] == 'H']
                self.assertEqual({r.instance[-1] for r in observed}, set(hydrogens))
                self.assertEqual(len({r.instance_name for r in observed}), 2)
                self.assertEqual(species.atom_eqv, original)

    def test_charge_multiplicity_and_double_bond_identity(self):
        species = point('F/C=C/F')
        trans = canonical_identity(species)['id']
        self.assertNotEqual(trans, canonical_identity(point('F/C=C\\F'))['id'])
        species.charge = 1
        self.assertNotEqual(trans, canonical_identity(species)['id'])
        species.charge = 0
        species.mult = 3
        self.assertNotEqual(trans, canonical_identity(species)['id'])


if __name__ == '__main__':
    unittest.main()
