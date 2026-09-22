import copy
import json
from pathlib import Path
from kinbot.parameters import Parameters
from kinbot.reaction_finder import ReactionFinder
from kinbot.stereochemistry import reaction_atom_equivalence
from test_reaction_paths import peroxy


def search(point, tmp_path, families, ringrange=(5, 7)):
    path = tmp_path / 'input.json'
    path.write_text(json.dumps({'families': families, 'family_restrict': 1,
                               'ringrange': list(ringrange), 'barrier_threshold': 100.}))
    original = copy.deepcopy(point.atom_eqv)
    finder = ReactionFinder(point, Parameters(path, show_warnings=False).par, None)
    finder.find_reactions()
    assert point.atom_eqv == original
    return {tuple(int(i)+1 for i in reaction.instance) for reaction in point.reac_obj}


def test_ho2_elimination_searches_both_diastereotopic_hydrogens(tmp_path):
    paths = search(peroxy(), tmp_path, ['HO2_Elimination_from_PeroxyRadical'])
    assert (10, 2, 3, 5, 6) in paths
    assert (11, 2, 3, 5, 6) in paths


def test_heavy_arm_is_not_removed_before_reaching_the_transfer_site(tmp_path):
    p = peroxy('C[C@H](O)[C@H](O[O])[C@@H](O)C')
    paths = search(p, tmp_path, ['intra_H_migration'])
    assert (5, 4, 2, 1, 10) in paths
    assert (5, 4, 2, 3, 14) in paths


def test_nonhydrogen_migration_uses_the_same_full_motif_refinement(tmp_path):
    # Two Cl atoms are diastereotopic in the presence of the remote fixed centre.
    p = peroxy('[CH2]C[C@H](F)C(Cl)Cl')
    chlorine = [i for i, atom in enumerate(p.atom) if atom == 'Cl']
    groups = reaction_atom_equivalence(p)
    assert not any(set(chlorine) <= set(group) for group in groups)
    paths = search(p, tmp_path, ['intra_R_migration'])
    assert {path[-1]-1 for path in paths} >= set(chlorine)


def test_achiral_control_keeps_existing_atom_equivalence():
    p = peroxy('CCCO[O]')
    assert {frozenset(g) for g in reaction_atom_equivalence(p)} == {frozenset(g) for g in p.atom_eqv}


def test_joint_h2_elimination_sites_survive_without_duplicate_reverse_searches(tmp_path):
    point = peroxy('CCCCC')
    point.mult = 1
    original = copy.deepcopy(point.atom_eqv)
    paths = search(point, tmp_path, ['h2_elim'], ringrange=(4, 5))
    # Terminal elimination and two relative configurations for internal H2
    # elimination. Each single hydrogen belongs to an unsplit class initially.
    assert paths == {(6, 1, 2, 9), (9, 2, 3, 11), (9, 2, 3, 12)}
    assert point.atom_eqv == original


def test_joint_selection_refines_meso_molecule_without_splitting_global_atoms(tmp_path):
    point = peroxy('C[C@H](O)C[C@@H](O)C')
    point.mult = 1
    paths = search(point, tmp_path, ['h2_elim'], ringrange=(4, 5))
    assert (11, 2, 4, 13) in paths
    assert (11, 2, 4, 14) in paths
    assert len(paths) == 4  # two other site classes, with reverse duplicates removed


def test_final_positional_filter_preserves_demonstrated_joint_stereo(tmp_path):
    # Two ring routes share the same end atoms, but the fixed sidechain
    # stereocentre makes the intervening ring arms diastereotopic.
    point = peroxy('C[C@H](O)C1CCC1')
    filename = tmp_path / 'filter.json'
    filename.write_text(json.dumps({'barrier_threshold': 100.}))
    finder = ReactionFinder(point, Parameters(filename, show_warnings=False).par, None)
    first, second = [12, 3, 4, 5, 15], [12, 3, 6, 5, 15]
    assert finder._motif_graph(first) == finder._motif_graph(second)
    assert finder._motif_key(first) != finder._motif_key(second)
    finder.reactions['h2_elim'] = []
    finder.new_reaction([first, second, first[::-1]], 'h2_elim', a=0, b=-1, cross=True)
    assert finder.reactions['h2_elim'] == [first, second]


def test_ring_discovery_names_keep_motifs_results_and_exclusions_separate(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from unittest.mock import Mock, patch
    from ase import Atoms
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from kinbot.stationary_pt import StationaryPoint
    from kinbot.reaction_generator import ReactionGenerator
    from kinbot.species_routing import connect, matches_name, resolve_job

    monkeypatch.chdir(tmp_path)
    filename = tmp_path / 'input.json'
    filename.write_text(json.dumps({'barrier_threshold': 100., 'delete_intermediate_files': 0}))
    par = Parameters(filename, show_warnings=False).par
    db = connect('kinbot.db')
    # Complete discovery must preserve both H choices and both internal ring arms.
    cases = [
        ('CC1CCCCC1', 'Intra_RH_Add_Endocyclic_R', [(4, 5, 6, 19), (4, 5, 6, 20)]),
        ('OC1CCCC1', 'Intra_RH_Add_Endocyclic_R', [(3, 4, 5, 14), (3, 4, 5, 15)]),
        ('CC1CCO1', 'Intra_RH_Add_Endocyclic_R', [(4, 3, 2, 9), (4, 3, 2, 10)]),
        ('C[C@H](O)C1CCC1', 'intra_H_migration', [(3, 4, 5, 15), (3, 6, 5, 15)]),
    ]
    for smiles, family, instances in cases:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        assert AllChem.EmbedMolecule(mol, randomSeed=312) == 0
        point = StationaryPoint('ring', 0, 1, atom=[a.GetSymbol() for a in mol.GetAtoms()],
                                geom=mol.GetConformer().GetPositions())
        point.characterize()
        finder = ReactionFinder(point, par, SimpleNamespace(db=db))
        finder.find_reactions()
        names = {(kind, tuple(r.instance)): r.instance_name
                 for kind, r in zip(point.reac_type, point.reac_obj)}
        assert len(point.reac_name) == len(set(point.reac_name))
        assert point.reac_name == [r.instance_name for r in point.reac_obj]
        first, second = [names[family, instance] for instance in instances]
        base = first.rsplit('_m', 1)[0]
        assert first != second and matches_name(second, [base])
        assert matches_name(first, [first]) and not matches_name(second, [first])

        # An old unsplit result cannot become either new candidate's result.
        atoms = Atoms(point.atom, positions=point.geom)
        db.write(atoms, name=base, data={'energy': -3.})
        for job, energy in [(first, -1.), (second, -2.)]:
            assert resolve_job(db, job) == job
            assert not list(db.select(name=job))
            db.write(atoms, name=job, data={'energy': energy})
            assert db.get(name=job).data.energy == energy
        assert db.get(name=base).data.energy == -3.

        # Reordering retained searches must not exchange their calculation names.
        other = StationaryPoint('ring', 0, 1, atom=point.atom, geom=point.geom)
        other.characterize()
        reordered = ReactionFinder(other, par, SimpleNamespace(db=db))
        for kind, motifs in reversed(list(finder.reactions.items())):
            reordered.reaction_matrix(list(reversed(motifs)), kind)
        reordered._name_distinct_motifs()
        assert names == {(kind, tuple(r.instance)): r.instance_name
                         for kind, r in zip(other.reac_type, other.reac_obj)}

        # A legacy exclusion removes every candidate before any QC access.
        skipped = dict(par, skip_reactions=[base] + [n for n in point.reac_name
                                                    if n not in (first, second)])
        qc = Mock()
        with patch('kinbot.reaction_generator.time.sleep'):
            ReactionGenerator(point, skipped, qc, str(filename)).generate()
        assert qc.mock_calls == []
        assert all(status == -999 for status in point.reac_ts_done)
