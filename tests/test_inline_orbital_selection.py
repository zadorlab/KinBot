"""User-selected, molecule-independent orbital indices reach native inputs."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from kinbot.execution import file_sha256
from kinbot.molpro import (active_orbital_swaps, continued_vrc_state_model,
                          _vrc_multireference_method, validate_vrc_state_model)
from kinbot.molpro_reconnaissance import select_record
from kinbot.orbital_review import install_selection
from kinbot.parameters import Parameters
from kinbot.vrc_correction import resolved_state_model, expected_levels


def selection(tmp_path, indices):
    checkpoint = tmp_path/'seed.wfu'
    checkpoint.write_bytes(b'reviewed water orbitals fixture')
    return {'orbital_selection': {'mode': 'molpro_record', 'path': str(checkpoint),
        'record': '2140.2', 'sha256': file_sha256(checkpoint), 'orbital_type': 'natural',
        'align_active_orbitals': False, 'active_indices': indices},
        'state_average': [{'multiplicity': 1, 'roots': 1}], 'selected_roots': {1: 1}}


def test_arbitrary_source_indices_are_packed_once_and_alignment_is_preserved(tmp_path):
    model = selection(tmp_path, [3, 6, 8])
    deck, _ = _vrc_multireference_method('caspt2-f12(4,3)', 10, 0, state_model=model)
    identities = list(range(1, 10))
    for a, b in active_orbital_swaps([3, 6, 8], 3):
        identities[a-1], identities[b-1] = identities[b-1], identities[a-1]
        assert f'rotate,{a}.1,{b}.1,0' in deck
    assert identities[3:6] == [3, 6, 8]
    assert 'occ,6\nclosed,3' in deck
    checkpoint = tmp_path/'new.wfu'; checkpoint.write_bytes(b'optimized target CAS')
    new = continued_vrc_state_model(model, checkpoint, nelectron=10, active_electrons=4, active_orbitals=3)
    assert new['orbital_selection']['active_indices'] == [4, 5, 6]
    assert new['orbital_selection']['align_active_orbitals'] is False
    following, _ = _vrc_multireference_method('caspt2-f12(4,3)', 10, 0, state_model=new)
    assert 'rotate,' not in following and 'diab,' not in following
    assert model['orbital_selection']['active_indices'] == [3, 6, 8]


@pytest.mark.parametrize('indices', [[1, 1, 2], [0, 2, 3], [True, 2, 3], [1, 2], [1.0, 2, 3]])
def test_invalid_orbital_indices_fail_before_native_execution(tmp_path, indices):
    with pytest.raises(ValueError, match='active_indices'):
        validate_vrc_state_model(selection(tmp_path, indices), active_electrons=4, active_orbitals=3)


def test_repacking_cannot_align_against_wrong_original_block(tmp_path):
    model = selection(tmp_path, [3, 6, 8]); model['orbital_selection']['align_active_orbitals'] = True
    with pytest.raises(ValueError, match='DIAB'):
        _vrc_multireference_method('caspt2-f12(4,3)', 10, 0, state_model=model)


def test_recon_selection_installs_one_self_contained_input_for_different_molecule(tmp_path):
    seed = tmp_path/'seed.wfu'; seed.write_bytes(b'water selection')
    selected = select_record(seed, tmp_path/'selection.json', record='2140.2', basis='avdz-f12',
        active_electrons=4, active_orbitals=3, active_indices=[3, 6, 8],
        state_average=[{'multiplicity': 1, 'roots': 1}], selected_roots={1: 1},
        review_note='user reviewed orbitals', orbital_type='natural')
    source = tmp_path/'input.json'
    source.write_text(json.dumps({'smiles': 'O', 'mult': 1, 'barrier_threshold': 1000, 'vrc_tst_sample_method': 'caspt2-f12(2,2)',
        'vrc_tst_sample_basis': 'avdz-f12', 'vrc_tst_scan_method': 'caspt2-f12(2,2)',
        'vrc_tst_scan_basis': 'avdz-f12', 'vrc_tst_high_method': 'caspt2-f12(2,2)',
        'vrc_tst_high_basis': 'avdz-f12'}))
    output = install_selection(source, selected, tmp_path/'complete.json')
    selected.unlink()  # The installed run must not read the selection JSON.
    par = Parameters(str(output)).par
    assert par['vrc_tst_sample_method'] == 'caspt2-f12(4,3)'
    assert not par['vrc_tst_state_model_file']
    model = resolved_state_model(par, ['O', 'H', 'H'], 0, None, basis='avdz-f12')
    assert model['orbital_selection']['active_indices'] == [3, 6, 8]
    assert expected_levels(par, symbols=['O', 'H', 'H'], charge=0)['state_models_by_basis']['avdz-f12'] == model
    child = tmp_path/'nested'/'worker.json'; child.parent.mkdir(); child.write_text(json.dumps(par))
    worker = Parameters(str(child)).par
    assert resolved_state_model(worker, ['O', 'H', 'H'], 0, None, basis='avdz-f12') == model
    with pytest.raises(ValueError, match='no orbital record'):
        resolved_state_model(worker, ['O', 'H', 'H'], 0, None, basis='avtz-f12')


def test_unreviewed_inline_chemistry_can_be_edited_but_not_launched(tmp_path):
    from test_vrc_multistate import model
    from kinbot.vrc_multistate import load_surface_model
    surface = model(); surface['reviewed'] = False
    path = tmp_path/'input.json'
    path.write_text(json.dumps({'smiles': 'O', 'mult': 1, 'barrier_threshold': 1000,
                                'vrc_tst_surface_model': surface}))
    par = Parameters(str(path)).par
    with pytest.raises(ValueError, match='reviewed=true'):
        load_surface_model(par)


def test_tight_preset_defaults_and_explicit_force_override(tmp_path):
    path = tmp_path/'input.json'
    base = {'smiles': 'O', 'mult': 1, 'barrier_threshold': 1000,
            'theory_preset': 'b3lyp-b2plyp-anl'}
    path.write_text(json.dumps(base))
    assert Parameters(str(path)).par['sella_fmax'] == 5e-4
    path.write_text(json.dumps({**base, 'sella_fmax': 5e-5, 'anl_geometry_fmax': 2e-4}))
    par = Parameters(str(path)).par
    assert par['sella_fmax'] == 5e-5 and par['anl_geometry_fmax'] == 2e-4
