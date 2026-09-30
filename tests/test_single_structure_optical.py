import copy
import numpy as np
from kinbot import symmetry
from kinbot.mess import MESS
from kinbot.stereo_identity import canonical_identity
from test_conformer_counting import peroxide
from test_configured_rotational_symmetry import optimized
from tests.counting_fixtures import saved_point, methanol_data


def test_harmonic_conformational_mirrors_receive_the_missing_weight():
    writer = MESS({'multi_conf_tst': 0, 'optical_population': 'specified'}, None)
    for p in (peroxide(), optimized('CCO')):
        symmetry.calculate_symmetry(p)
        p.freq = p.reduced_freqs = [500.] * (3*p.natom-6)
        assert not canonical_identity(p)['is_chiral_configuration']
        # Force ethanol into a gauche OH conformation if the embedded minimum is anti.
        if len(p.atom) > 4:
            from rdkit.Chem import rdMolTransforms
            from kinbot.molecular_symmetry import geometric_mirror_states
            if geometric_mirror_states(p) == 1:
                # A defined structural fixture avoids depending on MMFF conformer selection.
                from rdkit import Chem
                from rdkit.Chem import AllChem
                mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
                AllChem.EmbedMolecule(mol, randomSeed=1)
                h = next(a.GetIdx() for a in mol.GetAtomWithIdx(2).GetNeighbors() if a.GetSymbol() == 'H')
                rdMolTransforms.SetDihedralDeg(mol.GetConformer(), 0, 1, 2, h, 60.)
                p.geom = mol.GetConformer().GetPositions()
        assert writer._parent_symmetry(p) == p.sigma_ext / 2.


def test_meso_achiral_geometry_does_not_keep_legacy_nopt_two():
    p = optimized('C[C@H](Cl)[C@H](Cl)C')
    symmetry.calculate_symmetry(p)
    p.freq = p.reduced_freqs = [500.]*(3*p.natom-6)
    writer = MESS({'multi_conf_tst': 0}, p)
    assert p.sigma_ext == 1
    assert writer._parent_symmetry(p) == .5  # two gauche mirror wells / sigma1


def test_completed_methanol_ts_hir_does_not_receive_a_second_optical_two():
    data = methanol_data()
    # Same saved reference fixture used by the existing counting contract.
    entry = data['saddles'][0]
    p = saved_point(entry)
    writer = MESS({'multi_conf_tst': 0}, p)
    assert writer._parent_symmetry(p) == p.sigma_ext
