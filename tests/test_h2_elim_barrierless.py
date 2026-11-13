"""Declared cation H2 losses reuse fragment optimization without a TS job."""
import json

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from kinbot.parameters import Parameters
from kinbot.qc import QuantumChemistry
from kinbot.reaction_finder import ReactionFinder
from kinbot.reaction_generator import ReactionGenerator
from kinbot.reactions.reac_h2_elim_barrierless import H2ElimBarrierless
from kinbot.reactions.reac_homolytic_scission import HS
from kinbot.stationary_pt import StationaryPoint


def benzenium(charge=1):
    molecule = Chem.AddHs(Chem.MolFromSmiles('C1C=CC=C[CH+]1'))
    assert AllChem.EmbedMolecule(molecule, randomSeed=312) == 0
    point = StationaryPoint('benzenium', charge, 1 if charge else 2,
        atom=[atom.GetSymbol() for atom in molecule.GetAtoms()],
        geom=molecule.GetConformer().GetPositions())
    point.characterize()
    return point


def find(point, tmp_path, **settings):
    filename = tmp_path / 'input.json'
    options = dict(barrier_threshold=100.,
        families=['h2_elim_barrierless'], ringrange=[4, 9])
    options.update(settings)
    filename.write_text(json.dumps(options))
    par = Parameters(filename, show_warnings=False).par
    ReactionFinder(point, par, None).find_reactions()
    return par


@pytest.mark.parametrize('charge, count', [(1, 1), (0, 0), (-1, 0)])
def test_same_carbon_pair_only_for_cations(tmp_path, charge, count):
    point = benzenium(charge)
    find(point, tmp_path)
    assert len(point.reac_obj) == count
    if count:
        reaction = point.reac_obj[0]
        first, carbon, second = reaction.instance
        assert point.atom[carbon] == 'C'
        assert point.atom[first] == point.atom[second] == 'H'
        assert point.bond[first, carbon] == point.bond[second, carbon] == 1
        assert reaction.scan == reaction.max_step == 0


def test_declared_family_can_be_disabled(tmp_path):
    point = benzenium()
    find(point, tmp_path, families=['all'], skip_families=['h2_elim_barrierless'])
    assert 'h2_elim_barrierless' not in point.reac_type



def test_declared_exit_does_not_remove_the_saddle_search(tmp_path):
    point = benzenium()
    find(point, tmp_path, families=['h2_elim', 'h2_elim_barrierless'], ringrange=[3, 4])
    assert set(point.reac_type) == {'h2_elim', 'h2_elim_barrierless'}
    instances = {kind: set(map(tuple, (reaction.instance for reaction in point.reac_obj
                                     if reaction.family_name == family)))
                 for kind, family in [('h2_elim', 'h2elim'),
                                      ('h2_elim_barrierless', 'h2_elim_barrierless')]}
    assert instances['h2_elim'] == instances['h2_elim_barrierless']

def test_h2_geometry_charge_choices_and_parent_are_preserved(tmp_path):
    point = benzenium()
    find(point, tmp_path)
    original_geom = point.geom.copy()
    original_bond = point.bond.copy()
    reaction = point.reac_obj[0]
    product = reaction.get_barrierless_product()
    first, carbon, second = reaction.instance
    assert product.name == reaction.instance_name + '_prod'
    assert product.bond[first, carbon] == product.bond[second, carbon] == 0
    assert product.bond[first, second] == 1
    assert np.linalg.norm(product.geom[first] - product.geom[second]) == pytest.approx(.74)
    fragments, _ = product.start_multi_molecular(vary_charge=True)
    assert sorted((frag.natom, frag.charge, frag.mult) for frag in fragments) == [
        (2, 0, 1), (2, 1, 2), (11, 0, 2), (11, 1, 1)]
    np.testing.assert_array_equal(point.geom, original_geom)
    np.testing.assert_array_equal(point.bond, original_bond)


def test_homolytic_scission_still_only_deletes_its_bond(tmp_path):
    point = benzenium()
    par = find(point, tmp_path)
    reaction = HS(point, None, par, [0, 1], 'hom_sci_1_2')
    product = reaction.get_barrierless_product()
    expected = point.bond.copy()
    expected[0, 1] = expected[1, 0] = 0
    np.testing.assert_array_equal(product.bond, expected)
    np.testing.assert_array_equal(product.geom, point.geom)


def test_controller_starts_product_optimization_without_a_saddle(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    point = benzenium()
    par = find(point, tmp_path)
    qc = QuantumChemistry(par)
    qc.check_qc = lambda name: 0

    class ProductSubmitted(Exception):
        pass

    def submit(fragment, geometry, **kwargs):
        assert fragment.natom in (2, 11)
        assert kwargs['initial_product']
        raise ProductSubmitted

    qc.qc_opt = submit
    with pytest.raises(ProductSubmitted):
        ReactionGenerator(point, par, qc, 'input.json').generate()
    assert point.reac_ts_done == [2]
    assert point.reac_obj[0].ts is None
    assert isinstance(point.reac_obj[0], H2ElimBarrierless)
