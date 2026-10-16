"""SMILES starts choose the same stereoisomer without changing explicit input."""
import json
import os
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

import numpy as np
import pytest
from rdkit import Chem

from kinbot import cheminfo
from kinbot.species_routing import routing_name
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import canonical_identity


def point(smiles):
    species = StationaryPoint('input', 0, 1, smiles=smiles)
    species.characterize()
    return species


def test_smiles_names_and_geometries_repeat_in_separate_processes():
    # Includes unspecified centres, partial stereo, two centres, E/Z and axes.
    smiles = ['CCC(O)C', 'CC(Cl)C(Cl)C', 'C[C@H](F)C(Cl)C',
              'CC=CC', 'CC=C=CC', 'CC=C=CC=C=CC', 'c1ccccc1',
              'Fc1ccccc1-c1ccccc1Cl', 'CCO']
    script = '''
import json
from kinbot.stationary_pt import StationaryPoint
from kinbot.species_routing import routing_name
results = []
for smiles in json.loads(input()):
    species = StationaryPoint('input', 0, 1, smiles=smiles)
    species.characterize()
    results.append([routing_name(species), species.geom.tolist()])
print(json.dumps(results))
'''
    env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]),
               PYTHONDONTWRITEBYTECODE='1')
    runs = [json.loads(subprocess.check_output(
        [sys.executable, '-c', script], input=json.dumps(smiles), text=True,
        env=env).splitlines()[-1]) for _ in range(2)]
    assert runs[0] == runs[1]
    for name, geometry in runs[0]:
        assert name and np.isfinite(geometry).all()


@pytest.mark.parametrize('first,second,mirrors', [
    ('C[C@H](O)CC', 'C[C@@H](O)CC', True),
    ('C/C=C/C', 'C/C=C\\C', False),
])
def test_explicit_smiles_and_supplied_coordinates_keep_their_stereochemistry(first, second, mirrors):
    a, b = point(first), point(second)
    for smiles, species in ((first, a), (second, b)):
        observed = Chem.Mol(species.obmol)
        Chem.AssignStereochemistryFrom3D(observed, replaceExistingTags=True)
        assert Chem.MolToSmiles(Chem.RemoveHs(observed)) == Chem.MolToSmiles(
            Chem.MolFromSmiles(smiles))
    one, two = canonical_identity(a), canonical_identity(b)
    assert routing_name(a) != routing_name(b)
    assert (one['mirror_id'] == two['id']) == mirrors
    # A supplied geometry overrides SMILES and must never be embedded again.
    structure = [value for atom, xyz in zip(b.atom, b.geom)
                 for value in (atom, *xyz)]
    with patch.object(cheminfo, 'generate_3d_structure', side_effect=AssertionError):
        supplied = StationaryPoint('input', 0, 1, smiles=first, structure=structure)
        supplied.characterize()
    np.testing.assert_array_equal(supplied.geom, b.geom)
    assert routing_name(supplied) == routing_name(b)


def test_failed_embedding_has_no_random_fallback_and_aromatic_bonds_are_retained():
    with patch.object(cheminfo.AllChem, 'EmbedMolecule', return_value=-1), \
            patch.object(cheminfo.AllChem, 'MMFFOptimizeMolecule', side_effect=AssertionError):
        with pytest.raises(ValueError, match='supply coordinates'):
            point('CC(O)CC')
    _, _, bonds = cheminfo.generate_3d_structure('c1ccccc1', obabel=0)
    assert np.count_nonzero(bonds == 2) == 6


@pytest.mark.parametrize('smiles', ['C.O', 'CC.O', '[CH3].[CH3]', 'C.O.N'])
def test_disconnected_smiles_keep_their_fragments_separate(smiles):
    molecule, structure, bonds = cheminfo.generate_3d_structure(smiles, obabel=0)
    species = StationaryPoint('fragments', 0, 1, structure=structure)
    species.characterize()
    np.testing.assert_array_equal(species.bond01, bonds > 0)
    fragments, _ = species.start_multi_molecular()
    assert len(fragments) == len(Chem.GetMolFrags(molecule))
    _, repeated, _ = cheminfo.generate_3d_structure(smiles, obabel=0)
    assert repeated == structure
