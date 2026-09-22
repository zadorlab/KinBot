import copy
import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from kinbot.stationary_pt import StationaryPoint
from kinbot import symmetry
from kinbot.conformer_counting import evaluate_members
from kinbot.conformer_records import ConformerRecord
from test_molecular_symmetry import methoxy


def optimized(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
    assert AllChem.MMFFOptimizeMolecule(mol, maxIters=1000) == 0
    p = StationaryPoint(smiles, 0, 1, atom=[a.GetSymbol() for a in mol.GetAtoms()],
                        geom=mol.GetConformer().GetPositions())
    p.characterize()
    return p


@pytest.mark.parametrize('meso,same', [
    ('C[C@H](Cl)[C@H](Cl)C', 'C[C@H](Cl)[C@@H](Cl)C'),
    ('C[C@H](Cl)C[C@H](Cl)C', 'C[C@H](Cl)C[C@@H](Cl)C')])
def test_atom_and_bond_centred_rules_preserve_fixed_configuration(meso, same):
    for smiles, expected in [(meso, 1), (same, 2)]:
        p = optimized(smiles)
        atomid = list(p.atomid)
        eqv = copy.deepcopy(p.atom_eqv)
        symmetry._calculate_symmetry(p)
        internal = copy.deepcopy(p.sigma_int)
        symmetry.calculate_symmetry(p)
        assert p.sigma_ext == expected
        assert p.sigma_int == internal
        assert list(p.atomid) == atomid and p.atom_eqv == eqv


def test_meso_mirror_pair_is_counted_once_per_geometry():
    p = optimized('C[C@H](Cl)[C@H](Cl)C')
    records = [ConformerRecord(str(i), i, str(i), 'valid', geometry=g.tolist(),
                               zero_energy_hartree=-100., frequencies_cm1=(100.,))
               for i, g in enumerate((p.geom, p.geom * [-1., 1., 1.]))]
    one, _ = evaluate_members(p, records[:1])
    both, groups = evaluate_members(p, records)
    assert one[0].remaining_optical_weight / one[0].sigma_ext == 2.
    assert sum(both[i].remaining_optical_weight / both[i].sigma_ext for g in groups for i in g) == 2.


def test_ordinary_graph_conventions_stay_unchanged():
    p = methoxy()
    symmetry.calculate_symmetry(p)
    assert p.sigma_ext == 3
    for smi in ['CO', '[CH2]O', 'CC', 'CCCC', 'C1CCCCC1', 'C1=CC=CC=C1']:
        p = StationaryPoint(smi, 0, 2 if smi.startswith('[CH2]') else 1, smiles=smi)
        p.characterize()
        symmetry._calculate_symmetry(p)
        old = p.sigma_ext, copy.deepcopy(p.sigma_int)
        symmetry.calculate_symmetry(p)
        assert (p.sigma_ext, p.sigma_int) == old
