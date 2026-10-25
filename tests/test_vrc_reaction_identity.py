"""Same atom pair can label distinct accepted reactions, never the same model."""
from copy import deepcopy
import json
from types import SimpleNamespace

import pytest

from kinbot.execution import atomic_json
from kinbot.vrc_multistate import (validate_surface_model, select_surface_channel,
    extend_surface_model, channel_for_reaction, FragmentIdentityError, SurfaceModelRequired)
from test_vrc_multistate import model


def scoped():
    data = model()
    channel = data['channels'][0]
    channel['reaction'] = 'well_hom_sci_1_2'
    other = deepcopy(channel)
    other['reaction'] = 'well_rearrangement_exit'
    other['bottlenecks'][0]['states'][0]['energy_offset_cm_inverse'] = 42.
    data['channels'].append(other)
    return data


def test_same_bond_scoped_dispatch_is_exact_and_additive(tmp_path):
    from test_rotdpy_profiles import sulfur_correction
    data = scoped()
    validated = validate_surface_model(data)
    first, second = validated['channels']
    assert select_surface_channel(validated, [2,1], first['reaction']) == first
    assert select_surface_channel(validated, [1,2], second['reaction']) == second
    for name in (None, 'another_reaction'):
        with pytest.raises(SurfaceModelRequired):
            select_surface_channel(validated, [1,2], name)
    previous = deepcopy(data); previous['channels'] = previous['channels'][:1]
    assert extend_surface_model(previous, data) == validated
    changed = deepcopy(data); changed['channels'][0]['reaction'] = 'renamed'
    with pytest.raises(ValueError, match='migration'):
        extend_surface_model(previous, changed)
    par,_,_ = sulfur_correction(tmp_path)
    par.update(vrc_tst_surface_model=data, vrc_tst_surface_model_file='')
    products = [SimpleNamespace(mult=2, charge=0, smiles='[CH3]'),
                SimpleNamespace(mult=3, charge=0, smiles='[S]')]
    reaction = SimpleNamespace(scan_coo=[0,1], instance_name=second['reaction'], products=products)
    assert channel_for_reaction(par, reaction)[0] == second


@pytest.mark.parametrize('change', [
    lambda c: c[1].pop('reaction'),
    lambda c: c[1].update(reaction=c[0]['reaction']),
    lambda c: c[1].update(reaction=''),
    lambda c: c[1].update(reaction='../escape'),
])
def test_ambiguous_or_unsafe_scoping_refused(change):
    data = scoped(); change(data['channels'])
    with pytest.raises(ValueError):
        validate_surface_model(data)


def test_only_matching_request_acknowledged_and_identity_hold_stays(tmp_path):
    from test_vrc_review_status import setup_review
    from kinbot.vrc_review_status import acknowledge_installed
    par, species, directory, request = setup_review(tmp_path)
    par['vrc_tst_surface_model']['channels'][0]['reaction'] = request.parent.name
    other = request.parent.parent/'123_product_exit'/'request.json'
    record = json.loads(request.read_text()); record['reaction'] = other.parent.name
    atomic_json(other, record)
    results = acknowledge_installed(par, species, directory)
    assert {r['reaction']:r['status'] for r in results} == {
        request.parent.name:'model_accepted', other.parent.name:'needs_review'}
    record = json.loads(request.read_text())
    record.update(status='needs_review', review_kind='fragment_identity')
    atomic_json(request, record)
    assert all(r['status']=='needs_review' for r in acknowledge_installed(par, species, directory))


def test_reference_probe_resolves_nonhomolysis_exact_model(tmp_path):
    from test_rotdpy_profiles import sulfur_correction
    from kinbot.vrc_reference_probes import leaf_parameters
    par,_,_ = sulfur_correction(tmp_path)
    par.update(vrc_tst_surface_model=scoped(), vrc_tst_surface_model_file='')
    well = tmp_path/'well'; well.mkdir()
    atomic_json(well/'well.json', par)
    name = 'well_rearrangement_exit'
    atomic_json(well/'vrctst/orbital_review'/name/'request.json',
                dict(reaction=name, bond=[1,2], products=[dict(multiplicity=2),dict(multiplicity=3)]))
    leaf, fragments, mult, root = leaf_parameters(well, name+'__inner_m2_r1')
    assert (fragments,mult,root)==([2,3],2,1)
    assert leaf['_vrc_dividing_coordinate']=='reactive'
    with pytest.raises(ValueError):
        leaf_parameters(well, 'well_unknown__inner_m2_r1')


def test_legacy_model_cannot_acknowledge_two_reactions_on_same_bond(tmp_path):
    from test_vrc_review_status import setup_review
    from kinbot.vrc_review_status import acknowledge_installed
    par, species, directory, request = setup_review(tmp_path)
    record = json.loads(request.read_text())
    other = request.parent.parent/'123_product_exit'/'request.json'
    record['reaction'] = other.parent.name
    atomic_json(other, record)
    before = {p: p.read_bytes() for p in (request, other)}
    reports = acknowledge_installed(par, species, directory)
    assert len(reports) == 2
    assert all(r['status'] == 'needs_review' and 'reaction-specific' in r['error'] for r in reports)
    assert all(p.read_bytes() == content for p, content in before.items())


@pytest.mark.parametrize('mismatch', ['identity','charge','spin','count'])
def test_both_fragment_identities_checked_before_reordering(mismatch, monkeypatch):
    from kinbot import vrc_tst_scan as scan
    monkeypatch.setattr(scan, 'same_species', lambda a,b: a.identity==b.identity)
    monkeypatch.setattr(scan, 'reorder_coord', lambda **k: pytest.fail('Must fail before reordering'))
    parts = [SimpleNamespace(identity=i, charge=0, mult=2) for i in (1,2)]
    products = deepcopy(parts)
    if mismatch=='count': products.pop()
    else: setattr(products[1], {'identity':'identity','charge':'charge','spin':'mult'}[mismatch], 9)
    reaction = SimpleNamespace(parts=parts, products=products, maps=[[0],[1]])
    vts = scan.VTS.__new__(scan.VTS); vts.scan_reac={'reaction':reaction}
    before = deepcopy([p.__dict__ for p in products])
    with pytest.raises(FragmentIdentityError): vts.match_order('reaction')
    assert [p.__dict__ for p in products]==before


def test_reversed_valid_products_are_matched_both_ways(monkeypatch):
    from kinbot import vrc_tst_scan as scan
    monkeypatch.setattr(scan, 'same_species', lambda a,b: a.identity==b.identity)
    monkeypatch.setattr(scan, 'reorder_coord', lambda **k: k['map_B'])
    products=[SimpleNamespace(identity=i,charge=0,mult=2,reset_order=lambda:None) for i in (1,2)]
    rx=SimpleNamespace(parts=products,products=products[::-1],maps=[[0],[1]])
    vts=scan.VTS.__new__(scan.VTS);vts.scan_reac={'reaction':rx}
    vts.match_order('reaction')
    assert rx.products==products


def test_scoped_orbitals_install_into_only_one_channel(tmp_path):
    from kinbot.orbital_review import install_selection
    from test_orbital_review import _selected
    source=tmp_path/'input.json'; atomic_json(source, dict(mult=2,vrc_tst_surface_model=scoped()))
    selected=_selected(tmp_path)
    with pytest.raises(ValueError, match='match'):
        install_selection(source, selected, tmp_path/'bad.json', bond=[1,2], stage='inner')
    path=install_selection(source, selected, tmp_path/'good.json', bond=[1,2], stage='inner',
                           reaction='well_rearrangement_exit')
    channels=json.loads(path.read_text())['vrc_tst_surface_model']['channels']
    assert not channels[0]['bottlenecks'][0].get('parameters')
    assert channels[1]['bottlenecks'][0]['parameters']['vrc_tst_active_space']['active_orbitals']==4
