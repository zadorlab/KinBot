"""Product-side VRC discovery must survive both dispatch passes and retries."""
from copy import deepcopy
from types import SimpleNamespace

import numpy as np
import pytest

from kinbot.stationary_pt import StationaryPoint
from kinbot.vrc_tst_scan import VTS


def setup_exit():
    endpoint = StationaryPoint('endpoint', 0, 2, atom=['H', 'H', 'H'],
                              geom=np.asarray([[0., 0., 0.], [0., 0., .74], [0., 0., 3.]]))
    endpoint.characterize()
    products, _ = endpoint.start_multi_molecular()
    reaction = SimpleNamespace(instance_name='parent_exit', irc_prod=endpoint,
                               products=products, do_vdW=True)
    well = SimpleNamespace(name='parent', chemid='parent', reac_obj=[reaction])
    par = dict(vrc_tst_scan={'parent': ['parent_exit']}, vrc_tst_noscan={},
               pes=1, pp_orient='geometric', rotdpy_run=0)
    scanner = VTS(well, par, None)
    scanner.scan_reac[reaction.instance_name] = reaction
    return scanner, reaction


def test_coordinate_discovery_never_joins_fragments_or_changes_endpoint():
    scanner, reaction = setup_exit()
    before = deepcopy(reaction.irc_prod)
    for _ in range(3):
        scanner.find_scan_coos(['parent_exit'])
        assert reaction.scan_coo == [1, 2]
        assert len(reaction.irc_prod.start_multi_molecular()[0]) == 2
        assert np.array_equal(reaction.irc_prod.geom, before.geom)
        assert np.array_equal(reaction.irc_prod.bond, before.bond)
        assert np.array_equal(reaction.irc_prod.bonds, before.bonds)


def test_anchor_uses_endpoint_indices_after_product_atom_reordering(monkeypatch):
    scanner, reaction = setup_exit()
    original = scanner.match_order
    def reordered(name):
        original(name)
        reaction.maps[0] = reaction.maps[0][::-1]
    monkeypatch.setattr(scanner, 'match_order', reordered)
    scanner.find_scan_coos(['parent_exit'])
    assert reaction.scan_coo == [1, 2]


def test_full_selection_dispatch_and_sampling_reentry_preserves_endpoint(tmp_path, monkeypatch):
    from kinbot import vrc_orbitals, vrc_multistate, vrc_tst_scan, pes
    scanner, reaction = setup_exit()
    monkeypatch.chdir(tmp_path)
    retained = tmp_path/'existing_native.out'
    retained.write_bytes(b'previous accepted calculation')
    before = deepcopy(reaction.irc_prod)
    monkeypatch.setattr(vrc_tst_scan, 'configured_selection',
                        lambda mapping, well, default: mapping.get('parent', default))
    monkeypatch.setattr(vrc_tst_scan, 'routing_name', lambda species: species.name)
    monkeypatch.setattr(vrc_orbitals, 'asymptotic_fragments', lambda p: False)
    monkeypatch.setattr(vrc_orbitals, 'pivot_backend', lambda p: 'none')
    monkeypatch.setattr(scanner, 'configured_reactions', lambda p: ['parent_exit'] if p else [])
    monkeypatch.setattr(scanner, 'save_products', lambda names: None)
    monkeypatch.setattr(scanner, 'find_equiv', lambda names: None)
    monkeypatch.setattr(vrc_multistate, 'channel_for_reaction',
                        lambda p, r: ({'reaction': r.instance_name}, {}))
    dispatch, sampling = [], []
    def calculate(vts, models):
        assert list(models) == ['parent_exit']
        assert reaction.scan_coo == [1, 2]
        assert np.array_equal(reaction.irc_prod.geom, before.geom)
        assert len(reaction.irc_prod.start_multi_molecular()[0]) == 2
        dispatch.append(1)
    monkeypatch.setattr(vrc_multistate, 'calculate_channels', calculate)
    monkeypatch.setattr(pes, 'create_rotdpy_inputs',
                        lambda p, rows, *a, **k: sampling.append(rows))
    for _ in range(2):
        scanner.calculate_correction_potentials()
    assert len(dispatch) == len(sampling) == 2
    assert all(rows[0][1] == 'parent_exit' for rows in sampling)
    assert np.array_equal(reaction.irc_prod.bond, before.bond)
    assert retained.read_bytes() == b'previous accepted calculation'


def test_connected_endpoint_remains_a_fragment_identity_hold():
    from kinbot.vrc_multistate import FragmentIdentityError
    scanner, reaction = setup_exit()
    reaction.irc_prod.make_extra_bond(*reaction.irc_prod.start_multi_molecular())
    with pytest.raises(FragmentIdentityError, match='does not contain two fragments'):
        scanner.find_scan_coos(['parent_exit'])


@pytest.mark.parametrize('name', ['hcs2_h2_exit', 'ch2s_sh_exit'])
def test_previously_published_product_exit_anchors_survive_reentry(name):
    import json
    from test_product_exit_recon import GEOMS, PLAN
    match = json.loads((PLAN/(name+'.json')).read_text())['request_match']
    endpoint = StationaryPoint('endpoint', 0, 2, atom=np.array(match['symbols']),
                              geom=np.array(GEOMS[name]))
    endpoint.characterize()
    products, _ = endpoint.start_multi_molecular()
    reaction = SimpleNamespace(instance_name=match['reaction'], irc_prod=endpoint,
                               products=products, do_vdW=True)
    scanner = VTS(None, {}, None)
    scanner.scan_reac[reaction.instance_name] = reaction
    before = deepcopy(endpoint)
    for _ in range(3):
        scanner.find_scan_coos([reaction.instance_name])
        assert sorted(i+1 for i in reaction.scan_coo) == sorted(match['bond'])
        assert np.array_equal(endpoint.bond, before.bond)
        assert np.array_equal(endpoint.geom, before.geom)
