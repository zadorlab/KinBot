"""HIR chooses one TS per pathway without deleting calculated observations."""
import copy
import json

import pytest

from kinbot import bimol, constants
from kinbot.bimol_network import BimolecularNetwork
from tests.test_bimol_counting_contract import saved_network


def network_with_alternative():
    network = saved_network()
    network.par.update(rotor_scan=1, multi_conf_tst=0)
    for transition in network.transitions:
        bimol.update_ts_topology_from_irc(transition['reaction'])
    original = network.transitions[0]
    alternative = copy.deepcopy(original)
    alternative['id'] = 'alternative_ts'
    reaction = alternative['reaction']
    reaction.instance_name += '_alternative'
    alternative['members'] = [reaction.instance_name]
    reaction.ts.energy += 1. / constants.AUtoKCAL
    alternative['energy'] = network._relative_energy(reaction.ts)
    network.transitions.append(alternative)
    return network, original, alternative


def selected_ids(network):
    return set(network.hir_selection['selected_transition_ids'])


def test_completed_network_build_selects_both_methanol_channels():
    prepared, _, _ = network_with_alternative()
    network = BimolecularNetwork(prepared.species, prepared.par)
    assert len(selected_ids(network)) == 2
    assert all(t['hir_selection']['selected'] for t in network.transitions)


def test_hir_selection_uses_chemical_endpoints_not_complex_geometry(monkeypatch):
    network, original, alternative = network_with_alternative()
    reaction = alternative['reaction']
    for side, species in [('reactant', reaction.irc_reactant), ('product', reaction.irc_prod)]:
        node_id = 'different_' + side + '_complex'
        network.nodes[node_id] = dict(id=node_id, kind='complex', species=species,
                                     energy=network._relative_energy(species))
        alternative[side] = node_id
    # Selection does not repeat the geometric duplicate test or require the
    # TSs to fall within the numerical-duplicate energy tolerance.
    monkeypatch.setattr('kinbot.bimol_network.symmetry_rmsd',
                        lambda *a, **kw: pytest.fail('Selection must use pathway classes'))
    before = copy.deepcopy(network.as_dict())
    network._select_hir_models()
    assert len(network.transitions) == 3
    assert len(selected_ids(network)) == 2  # Methanol C-H and O-H abstraction.
    assert original['id'] in selected_ids(network)
    assert alternative['hir_selection']['representative'] == original['id']
    assert not alternative['hir_selection']['selected']
    after = network.as_dict()
    for transition in after['transitions']:
        transition.pop('hir_selection')
    after['hir_selection'] = None
    assert after == before  # Keep all energies, HIR, endpoints and archives.


def test_selection_uses_energy_plus_zpe_and_is_order_independent():
    network, original, alternative = network_with_alternative()
    alternative['species'].energy = original['species'].energy - .01
    alternative['species'].zpe = original['species'].zpe + .02
    alternative['energy'] = network._relative_energy(alternative['species'])
    network._select_hir_models()
    first = copy.deepcopy(network.hir_selection)
    assert original['id'] in selected_ids(network)
    network.transitions.reverse()
    network._select_hir_models()
    assert network.hir_selection == first
    alternative['species'].zpe = original['species'].zpe
    alternative['energy'] = network._relative_energy(alternative['species'])
    network._select_hir_models()
    assert alternative['id'] in selected_ids(network)
    assert original['hir_selection']['representative'] == alternative['id']


def test_distinct_stereochemical_paths_survive_with_the_same_products():
    network, original, alternative = network_with_alternative()
    for transition, path in [(original, 'htransfer:first'), (alternative, 'htransfer:second')]:
        transition['species'].stereopath_id = path
        transition['species'].stereopath_metadata['id'] = path
    network._select_hir_models()
    assert len(selected_ids(network)) == 3
    assert original['hir_selection']['selected']
    assert alternative['hir_selection']['selected']


def test_configured_product_identities_remain_distinct():
    network, original, alternative = network_with_alternative()
    # Nodes already distinguish configured fragment identities. Two different
    # such nodes must not collapse even for the same ordinary TS path label.
    product_id = alternative['reaction'].product_asymptote_id
    other_id = 'other_configured_product'
    network.nodes[other_id] = dict(network.nodes[product_id], id=other_id)
    alternative['reaction'].product_asymptote_id = other_id
    network._select_hir_models()
    assert len(selected_ids(network)) == 3


@pytest.mark.parametrize('invalid', ['missing_path', 'nonfinite_energy'])
def test_invalid_observation_does_not_replace_a_valid_hir_model(invalid):
    network, original, alternative = network_with_alternative()
    alternative['energy'] = original['energy'] - 10.
    if invalid == 'missing_path':
        del alternative['species'].stereopath_metadata
    else:
        alternative['energy'] = float('nan')
    network._select_hir_models()
    assert original['id'] in selected_ids(network)
    assert not alternative['hir_selection']['selected']
    assert alternative['hir_selection']['representative'] is None
    assert alternative['hir_selection']['reason']
    assert alternative in network.transitions


def test_observed_other_product_can_supply_the_selected_hir_model():
    network, original, alternative = network_with_alternative()
    alternative['energy'] = original['energy'] - .1
    alternative['product_matches_expected'] = False
    network._select_hir_models()
    assert alternative['id'] in selected_ids(network)
    assert original['id'] not in selected_ids(network)


def test_mc_and_harmonic_runs_keep_their_existing_selection():
    network, _, _ = network_with_alternative()
    network._select_hir_models()
    for rotor_scan, mc in [(0, 0), (0, 1), (1, 1)]:
        network.par.update(rotor_scan=rotor_scan, multi_conf_tst=mc)
        network._select_hir_models()
        assert network.hir_selection is None
        assert all('hir_selection' not in t for t in network.transitions)
        assert len(network.transitions) == 3


def test_selection_does_not_clear_observation_conflicts():
    network, original, alternative = network_with_alternative()
    original['endpoint_conflicts'] = [alternative['id']]
    original['flux_assignment'] = 'unresolved observation conflicts; no independent flux assigned'
    network._select_hir_models()
    assert original['hir_selection']['selected']
    assert original['endpoint_conflicts'] == [alternative['id']]
    assert original['flux_assignment'].startswith('unresolved')
    assert not network.as_dict()['rate_model_ready']


def test_both_exports_keep_all_observations_and_the_same_hir_selection(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    network, original, alternative = network_with_alternative()
    network._select_hir_models()
    network.write_json()
    network.write_unique_ts()
    data = json.loads((tmp_path / 'bimol_network.json').read_text())
    directory = tmp_path / 'bimol/unique_ts'
    index = json.loads((directory / 'index.json').read_text())
    assert data['hir_selection'] == index['hir_selection'] == network.hir_selection
    assert not data['rate_model_ready']
    assert len(data['transitions']) == len(index['unique_ts']) == 3
    for transition in data['transitions']:
        metadata = json.loads((directory / transition['id'] / 'metadata.json').read_text())
        assert metadata['hir_selection'] == transition['hir_selection']
        assert metadata['thermochemistry'] == transition['stationary_point']['thermochemistry']
        assert (directory / transition['id'] / 'ts.xyz').is_file()
        assert (directory / transition['id'] / 'irc_reactant_endpoint.xyz').is_file()
        assert (directory / transition['id'] / 'irc_product_endpoint.xyz').is_file()
    assert (directory / alternative['id'] / 'ts.xyz').is_file()
