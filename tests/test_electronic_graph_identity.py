"""Structural stereo identity must not interpret valence deficits as electrons.

The peroxide coordinates are saved calculations; the thermochemical numbers in
the writer test are synthetic and exercise serialization only.
"""
import copy
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import numpy as np
from ase.db import connect
from rdkit import Chem
from rdkit.Chem import AllChem

from kinbot.conformer_counting import representative_record
from kinbot.irc import IRC
from kinbot.mess import MESS
from kinbot.reaction_generator import ReactionGenerator
from kinbot.species_routing import routing_key, same_species
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity, optical_scope
from kinbot.stereo_routing import guard_well_job, require_same_configuration
from kinbot.symmetry import calculate_symmetry


FIXTURES = json.loads((Path(__file__).parent / 'reference' /
                       'peroxy_irc_endpoints.json').read_text())


def endpoint(record, side):
    p = StationaryPoint(side, record['charge'], record['multiplicity'],
                       atom=record['atoms'], geom=np.array(record[side + '_geometry']))
    p.characterize()
    return p


def from_smiles(smiles, charge, multiplicity):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=312) == 0
    p = StationaryPoint(smiles, charge, multiplicity,
        atom=[a.GetSymbol() for a in mol.GetAtoms()],
        geom=mol.GetConformer().GetPositions())
    p.characterize()
    return p


class TestElectronicGraphIdentity(unittest.TestCase):
    def test_saved_peroxide_endpoints_keep_valid_identity_and_irc_product(self):
        for record in FIXTURES:
            with self.subTest(run=record['run']):
                reactant, product = (endpoint(record, side)
                                     for side in ('reactant', 'product'))
                self.assertTrue(any(np.any(rad < 0) for rad in product.rads))
                arrays = {key: copy.deepcopy(getattr(product, key))
                          for key in ('geom', 'atom', 'bond', 'bonds', 'rads')}
                identity = canonical_identity(product)
                self.assertEqual(identity['status'], 'assigned')
                self.assertEqual(identity['formal_charge_localization'], 'not supplied')
                self.assertIn('not inferred', identity['electronic_localization'])
                self.assertEqual(canonical_identity(reactant)['id'], record['reactant_identity'])
                self.assertNotEqual(identity['id'], record['reactant_identity'])
                # Virtual substitution must also work on the observed product.
                hydrogen = int(record['channel'].split('_')[1]) - 1
                self.assertEqual(canonical_identity(product, tagged_atom=hydrogen)['status'], 'assigned')
                qc = SimpleNamespace(get_qc_geom=lambda job, *a, **kw:
                    (0, product.geom if '_F_' in job else reactant.geom))
                reaction = SimpleNamespace(species=reactant, qc=qc, instance_name='saved')
                observed = IRC(reaction, {'bimol': 0, 'optical_population': 'specified'}).irc2stationary_pt()
                self.assertIsInstance(observed, StationaryPoint)
                self.assertEqual(canonical_identity(observed)['id'], identity['id'])
                for key, before in arrays.items():
                    np.testing.assert_array_equal(getattr(product, key), before)

    def test_product_reuse_qc_guard_and_mc_mess_need_no_electronic_guess(self):
        with TemporaryDirectory() as directory:
            original = Path.cwd()
            self.addCleanup(os.chdir, original)
            os.chdir(directory)
            for record in FIXTURES:
                with self.subTest(run=record['run']):
                    product = endpoint(record, 'product')
                    products, _ = product.start_multi_molecular()
                    self.assertEqual(len(products), 1)
                    other = endpoint(record, 'product')
                    self.assertTrue(same_species(product, other))
                    require_same_configuration(product, other, 'product reuse')
                    objects = [product, other]
                    ReactionGenerator.equate_identical(None, objects)
                    self.assertIs(objects[0], objects[1])
                    qc = SimpleNamespace(db=connect(record['run'] + '.db'), par={'multi_conf_tst': 1})
                    guard_well_job(qc, product, product.geom, str(routing_key(product)) + '_well')
                    product.energy, product.zpe = -100., .1
                    product.freq = product.reduced_freqs = [500.] * (3*product.natom - 6)
                    product.conformer_index = [0]
                    product.conformer_geom = [product.geom.copy()]
                    product.conformer_freq = [product.freq[:]]
                    product.conformer_zeroenergy = [product.energy + product.zpe]
                    calculate_symmetry(product)
                    counted = representative_record(product)
                    # Sec-butyl retains a specified stereocentre. Propyl has
                    # only conformational chirality, so its mirror is allowed.
                    self.assertEqual(counted.remaining_optical_weight,
                                     1. if record['run'].startswith('secbutyl') else 2.)
                    writer = MESS({'pes': 0, 'multi_conf_tst': 1, 'rotor_scan': 0,
                                   'optical_population': 'specified',
                                   'freq_uq_ref': 100., 'freq_uq_max_exp': 2.}, product)
                    writer.well_names = {routing_key(product): 'product'}
                    text = writer.write_well(product, 0., 1., 0)
                    self.assertIn('End ! RRHO', text)
                    self.assertNotIn('nan', text.lower())
                    divisor = next(float(line.split()[1]) for line in text.splitlines()
                                   if line.strip().startswith('SymmetryFactor'))
                    self.assertEqual(divisor, counted.sigma_ext / counted.remaining_optical_weight)

    def test_ordinary_radicals_ions_and_resonance_graphs_retain_structural_keys(self):
        cases = [('CO', 0, 1), ('C[O]', 0, 2), ('[CH2]O', 0, 2),
                 ('[CH2]', 0, 3), ('[O][O]', 0, 3), ('C[O-]', -1, 1),
                 ('C[O+](C)C', 1, 1), ('[NH4+]', 1, 1),
                 ('C[N+](=O)[O-]', 0, 1), ('CS(=O)(=O)C', 0, 1),
                 ('C/C=C/C', 0, 1), ('C[C@H](O)CC', 0, 1),
                 ('[CH2]C=C', 0, 2)]
        for smiles, charge, mult in cases:
            with self.subTest(smiles=smiles):
                p = from_smiles(smiles, charge, mult)
                original = canonical_identity(p)
                self.assertEqual(original['status'], 'assigned')
                # The rads field is not an independently assigned electron state.
                p.rads = []
                self.assertEqual(canonical_identity(p)['id'], original['id'])
                p.mult += 2
                self.assertNotEqual(canonical_identity(p)['id'], original['id'])
                p.mult = mult
                p.charge += 1
                self.assertNotEqual(canonical_identity(p)['id'], original['id'])

    def test_configured_ions_keep_supplied_charges_stereo_and_atom_order_invariance(self):
        p = from_smiles('C[O+](C)[C@H](F)Cl', 1, 1)
        self.assertTrue(np.any(p.rad < 0))
        unknown = canonical_identity(p)
        p.formal_charges = [1 if atom == 'O' else 0 for atom in p.atom]
        before = canonical_identity(p)
        self.assertNotEqual(before['id'], unknown['id'])
        self.assertEqual(before['formal_charge_localization'], 'specified')
        self.assertTrue(before['is_chiral_configuration'])
        self.assertFalse(optical_scope(p)['mirror_allowed'])
        self.assertTrue(optical_scope(p, 'racemic')['mirror_allowed'])
        order = np.arange(p.natom)[::-1]
        q = copy.copy(p)
        for name in ('atom', 'geom', 'formal_charges'):
            setattr(q, name, np.asarray(getattr(p, name))[order])
        q.bond = p.bond[np.ix_(order, order)]
        q.bonds = [matrix[np.ix_(order, order)] for matrix in p.bonds]
        q.rads = [np.asarray(rad)[order] for rad in p.rads]
        self.assertEqual(canonical_identity(q)['id'], before['id'])
        mirror = canonical_identity(q, q.geom * [-1., 1., 1.])
        self.assertEqual(mirror['id'], before['mirror_id'])
        self.assertEqual(mirror['mirror_family_id'], before['mirror_family_id'])


if __name__ == '__main__':
    unittest.main()
