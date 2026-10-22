"""Classify observed H transfers without changing their original search targets."""
import copy
import json
from types import SimpleNamespace as NS

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from kinbot import bimol
from kinbot.bimol_rate_models import prepare_models
from kinbot.irc import IRC, BimolecularIRCResult
from kinbot.parameters import Parameters
from kinbot.reaction_finder_bimol import ReactionFinderBimol
from kinbot.reaction_path import prepare_stereopath
from kinbot.stationary_pt import StationaryPoint
from tests.test_bimol_counting_contract import saved_network


def _reactants(first_smiles, second_smiles, first_mult, second_mult):
    fragments = []
    for i, (smiles, mult) in enumerate(((first_smiles, first_mult),
                                       (second_smiles, second_mult))):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        assert AllChem.EmbedMolecule(mol, randomSeed=312) == 0
        point = StationaryPoint(f'fragment_{i}', 0, mult,
            atom=[atom.GetSymbol() for atom in mol.GetAtoms()],
            geom=mol.GetConformer().GetPositions() + [10. * i, 0., 0.])
        point.characterize()
        point.energy, point.zpe = -1., .01
        fragments.append(point)
    a, b = fragments
    parent = StationaryPoint('reactants', 0, 2,
        atom=list(a.atom) + list(b.atom), geom=np.vstack([a.geom, b.geom]),
        fragA=a, fragB=b)
    parent.characterize()
    parent.energy, parent.zpe = -2., .02
    return parent


def _enumerate(species, tmp_path):
    source = tmp_path / 'input.json'
    source.write_text(json.dumps({'barrier_threshold': 50.}))
    parameters = Parameters(str(source)).par
    parameters['bimol'] = 1
    finder = ReactionFinderBimol(species, parameters, None)
    finder.find_reactions()
    return species.reac_obj


def _transferred_product(species, abstractor, hydrogen):
    geom = species.geom.copy()
    neighbors = np.flatnonzero(species.bond01[abstractor])
    if len(neighbors):
        existing = geom[neighbors[0]] - geom[abstractor]
        axis = np.array([1., 0., 0.])
        if abs(np.dot(existing / np.linalg.norm(existing), axis)) > .8:
            axis = np.array([0., 1., 0.])
        direction = np.cross(existing, axis)
        direction /= np.linalg.norm(direction)
    else:
        direction = np.array([0., 1., 0.])
    distance = .74 if species.atom[abstractor] == 'H' else .97
    geom[hydrogen] = geom[abstractor] + distance * direction
    product = StationaryPoint('observed_product', species.charge, species.mult,
                              atom=species.atom, geom=geom)
    product.characterize()
    return product


def test_equivalent_methane_h_is_accepted_by_l0_and_final_irc(tmp_path, monkeypatch):
    species = _reactants('C', '[OH]', 1, 2)
    reaction, = _enumerate(species, tmp_path)
    target = list(reaction.instance)
    data = reaction.get_bimol_search_data()
    alternative = next(choice for choice in data['equivalent_bond_changes']
                       if choice['instance'] != target)
    abstractor, hydrogen, donor = alternative['instance']
    product = _transferred_product(species, abstractor, hydrogen)
    assert donor == target[2] and hydrogen != target[1]
    assert reaction.site_degeneracy == 4
    assert not bimol.matches_expected_product(species, product.geom,
        data['formed_bonds'], data['broken_bonds'])
    assert bimol.matches_expected_product(species, product.geom,
        data['formed_bonds'], data['broken_bonds'], data['equivalent_bond_changes'])

    endpoints = [dict(state=state, converged=True, geometry_angstrom=point.geom.tolist())
                 for state, point in (('reactant', species), ('unexpected', product))]
    monkeypatch.setattr(bimol, '_irc_endpoints', lambda *args: copy.deepcopy(endpoints))
    checked = bimol.validate_candidate(species, species.geom,
        data['formed_bonds'], data['broken_bonds'], None, tmp_path / 'l0',
        search_kwargs={'equivalent_bond_changes': data['equivalent_bond_changes']})
    assert checked['valid'] and checked['product_matches_expected']

    reaction.qc = NS(for_level=lambda level: reaction.qc,
        get_qc_geom=lambda name, *args, **kwargs: (
            0, species.geom if '_IRC_F_' in name else product.geom))
    result = IRC(reaction, {'bimol': 1}).irc2stationary_pt()
    assert isinstance(result, BimolecularIRCResult)
    assert result.product_matches_expected
    reaction.irc_reactant, reaction.irc_prod = result.reactant, result.product
    reaction.ts = copy.deepcopy(species)
    bimol.update_ts_topology_from_irc(reaction)
    assert reaction.instance == target
    assert reaction.observed_site['instance'] == [abstractor, hydrogen, donor]
    assert reaction.observed_site['site_degeneracy'] == 4


def test_other_enumerated_transfer_uses_its_own_site_counting(tmp_path):
    species = _reactants('CCO', '[H]', 1, 2)
    reactions = _enumerate(species, tmp_path)
    reaction = next(r for r in reactions if r.instance[2] == 0)
    other = next(r for r in reactions if species.atom[r.instance[2]] == 'O')
    target = list(reaction.instance)
    original_class = copy.deepcopy((reaction.stereopath_id, reaction.site_degeneracy,
                                   reaction.site_members, reaction.site_relations))
    assert reaction.site_degeneracy == 3 and other.site_degeneracy == 1
    abstractor, hydrogen, donor = other.instance
    reaction.irc_reactant = species
    reaction.irc_prod = _transferred_product(species, abstractor, hydrogen)
    reaction.ts = copy.deepcopy(species)
    bimol.update_ts_topology_from_irc(reaction)
    assert reaction.instance == target
    assert (reaction.stereopath_id, reaction.site_degeneracy,
            reaction.site_members, reaction.site_relations) == original_class
    assert reaction.observed_site['instance'] == [abstractor, hydrogen, donor]
    assert reaction.observed_site['site_degeneracy'] == other.site_degeneracy
    assert reaction.observed_site['stereopath_id'] == other.stereopath_id
    assert reaction.observed_site['site_members'] == other.site_members


def test_unenumerated_product_survives_network_hir_and_rate_preparation():
    network = saved_network()
    reaction = network.species.reac_obj[0]
    target = list(reaction.instance)
    original_degeneracy = reaction.site_degeneracy
    completed = network.species.reac_obj
    for candidate in completed:
        prepare_stereopath(candidate.ts, candidate.irc_reactant, candidate.irc_prod)
    # A valid route with no enumerated reaction must not inherit another
    # search's number of equivalent hydrogens.
    network.species.reac_obj = [NS(bimol_search=False)]
    bimol.assign_observed_site(reaction)
    network.species.reac_obj = completed
    assert reaction.observed_site['family_name'] == 'bimolecular'
    assert reaction.observed_site['site_degeneracy'] == 1
    assert list(reaction.instance) == target
    assert reaction.site_degeneracy == original_degeneracy
    reaction.product_matches_expected = False
    network.transitions.clear()
    network._cluster_transitions(completed)
    observed = next(t for t in network.transitions if not t['product_matches_expected'])
    assert observed['family'] == 'bimolecular' and observed['degeneracy'] == 1
    expected_atoms = sorted({atom for pairs in reaction.actual_bond_changes.values()
                             for pair in pairs for atom in pair})
    assert observed['reaction'].instance == expected_atoms
    assert list(reaction.instance) == target
    network.par.update(rotor_scan=1, multi_conf_tst=0)
    network._select_hir_models()
    assert observed['hir_selection']['selected']
    models = prepare_models(network)
    assert len(models['pathways']) == 2 and not models['warnings']
    assert observed['id'] in {path['id'] for path in models['pathways']}
