"""Spin-splitting VRC correction arithmetic and input contracts."""

import os
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from kinbot.anl.extrapolation import two_point_cbs
from kinbot.molpro import (Molpro, _vrc_multireference_method,
                           automatic_vrc_state_model,
                           compact_molpro_titles,
                           full_valence_active_space,
                           stage_vrc_state_model,
                           validate_vrc_state_model,
                           vrc_multireference_energy_contract)
from kinbot.vrc_correction import (
    active_space_feasibility, active_space_settings,
    apply_vrc_accuracy_preset,
    assemble_sampling_reference, assemble_spin_splitting_reference,
    correction_components, coupled_multiplicities, expected_levels,
    resolved_state_model, sampling_settings, spin_splitting_settings,
    spin_splitting_weights, validate_vrc_levels,
)
from kinbot.vrc_tst_scan import VTS


def _parameters(**updates):
    values = {
        'vrc_tst_correction_model': 'spin_splitting',
        'vrc_tst_spin_splitting': None,
        'vrc_tst_sample_method': 'caspt2(2,2)',
        'vrc_tst_sample_basis': 'vdz',
        'vrc_tst_high_method': 'mrci+q(2,2)',
        'vrc_tst_high_basis': 'avtz',
    }
    values.update(updates)
    return values


def test_default_spin_splitting_uses_tz_qz_cbs_for_all_three_branches():
    settings = spin_splitting_settings(_parameters())
    assert settings['triplet_reference']['basis'] == [
        'avtz-f12', 'avqz-f12']
    assert [level['basis'] for level in settings['gap_levels']] == [
        ['avtz-f12', 'avqz-f12'], ['avtz-f12', 'avqz-f12']]
    assert all(level['extrapolation_power'] == pytest.approx(3.7)
               for level in [settings['triplet_reference'],
                             *settings['gap_levels']])
    assert [level['method'] for level in settings['gap_levels']] == [
        'caspt2-f12(2,2)', 'mrci-f12+q(2,2)']
    components = correction_components(_parameters())
    assert len(components) == 11
    assert {component['multiplicity'] for component in components[1:]} == {
        1, 3}
    assert expected_levels(_parameters())['trusted_correction']['model'] == \
        'spin_splitting'


def test_general_triplet_pair_uses_highest_available_reference_spin():
    settings = spin_splitting_settings(
        _parameters(), parent_multiplicity=3, fragment_multiplicities=[3, 3])
    assert settings['target_multiplicity'] == 3
    assert settings['reference_multiplicity'] == 5


def test_direct_trusted_cbs_has_distinct_components_and_correct_arithmetic():
    from kinbot.vrc_correction import assemble_trusted_reference
    parameters = _parameters(vrc_tst_correction_model='direct',
        vrc_tst_high_basis=['avtz-f12', 'avqz-f12'],
        vrc_tst_high_extrapolation_power=3.7)
    components = correction_components(parameters)
    assert [(item['role'], item['basis']) for item in components[1:]] == [
        ('trusted:0', 'avtz-f12'), ('trusted:1', 'avqz-f12')]
    assert len({item['suffix'] for item in components}) == len(components)
    result, metadata = assemble_trusted_reference([3., 30.],
        {'trusted:0': [-10., -9.], 'trusted:1': [-10.1, -9.1]}, parameters)
    assert result == pytest.approx([
        two_point_cbs(-10., -10.1, upper_cardinal=4, power=3.7),
        two_point_cbs(-9., -9.1, upper_cardinal=4, power=3.7)])
    assert metadata['basis'] == ['avtz-f12', 'avqz-f12']


def test_declared_fragment_context_also_selects_the_high_spin_anchor():
    settings = spin_splitting_settings(_parameters(vrc_tst_spin_splitting={
        'target_multiplicity': 3, 'fragment_multiplicities': [3, 3]}))
    assert settings['reference_multiplicity'] == 5


def test_already_high_spin_target_does_not_silently_use_low_spin_cc_anchor():
    with pytest.raises(ValueError, match='already the highest'):
        spin_splitting_settings(
            _parameters(), parent_multiplicity=4, fragment_multiplicities=[2, 3])
    # An explicit expert reference remains possible, with the parent fixed.
    settings = spin_splitting_settings(_parameters(vrc_tst_spin_splitting={
        'reference_multiplicity': 2}), parent_multiplicity=4,
        fragment_multiplicities=[2, 3])
    assert settings['reference_multiplicity'] == 2


@pytest.mark.parametrize('override,message', [
    ({'target_multiplicity': 4}, 'accepted parent'),
    ({'fragment_multiplicities': [2, 5]}, 'accepted products'),
])
def test_explicit_spin_settings_cannot_relabel_the_accepted_channel(override, message):
    with pytest.raises(ValueError, match=message):
        spin_splitting_settings(_parameters(vrc_tst_spin_splitting=override),
                                parent_multiplicity=2, fragment_multiplicities=[2, 3])


def test_distance_switch_prefers_mrci_short_and_caspt2_long_range():
    settings = spin_splitting_settings(_parameters(vrc_tst_spin_splitting={
        'weighting': {'model': 'distance_switch',
                     'distance_knots_angstrom': [2.5, 3.125, 3.75, 4.375, 5.0],
                     'weights': {'caspt2_f12': [.2, .4, .5, .6, .8],
                                 'mrci_f12_q': [.8, .6, .5, .4, .2]}}}))
    weights = spin_splitting_weights([2., 2.5, 3.125, 3.75, 4.375, 5., 8.],
                                     settings)
    assert weights['caspt2_f12'] == pytest.approx(
        [.2, .2, .4, .5, .6, .8, .8])
    assert weights['mrci_f12_q'] == pytest.approx(
        [.8, .8, .6, .5, .4, .2, .2])
    assert np.asarray(weights['caspt2_f12']) + np.asarray(
        weights['mrci_f12_q']) == pytest.approx(np.ones(7))


def test_default_weights_have_no_molecule_specific_distance_preference():
    settings = spin_splitting_settings(_parameters())
    assert settings['weighting'] == {'model': 'arithmetic_mean'}
    weights = spin_splitting_weights([1., 5., 50.], settings)
    assert weights['caspt2_f12'] == pytest.approx([.5, .5, .5])
    with pytest.raises(ValueError):
        spin_splitting_settings(_parameters(vrc_tst_spin_splitting={
            'weighting': {'model': 'distance_switch'}}))


def test_sulfur_radical_spin_couplings_and_high_accuracy_active_spaces():
    assert coupled_multiplicities(2, 3) == (2, 4)
    assert coupled_multiplicities(2, 2) == (1, 3)
    settings = spin_splitting_settings(_parameters(
        vrc_tst_sample_method='caspt2-f12b(13,8)',
        vrc_tst_sample_basis=['avtz-f12', 'avqz-f12'],
        vrc_tst_spin_splitting={
            'target_multiplicity': 2,
            'reference_multiplicity': 4,
            'fragment_multiplicities': [2, 3],
            'reference': {
                'method': 'uccsd(t)-f12b',
                'basis': ['vqz-f12', 'v5z-f12'],
            },
            'gap_levels': [
                {'label': 'caspt2_f12', 'method': 'caspt2-f12b(13,8)',
                 'basis': ['vqz-f12', 'v5z-f12'], 'shift': .3,
                 'ipea': .25},
                {'label': 'mrci_f12_q', 'method': 'mrci-f12b+q(13,8)',
                 'basis': ['vqz-f12', 'v5z-f12']},
            ],
        }))
    assert settings['target_multiplicity'] == 2
    assert settings['reference_multiplicity'] == 4
    assert [item['method'] for item in settings['gap_levels']] == [
        'caspt2-f12b(13,8)', 'mrci-f12b+q(13,8)']


def test_whole_molecule_full_valence_spaces_and_feasibility_are_derived():
    ch3s = full_valence_active_space(['C', 'H', 'H', 'H', 'S'])
    ch3ss = full_valence_active_space(['C', 'H', 'H', 'H', 'S', 'S'])
    assert (ch3s['active_electrons'], ch3s['active_orbitals']) == (13, 8)
    assert (ch3ss['active_electrons'], ch3ss['active_orbitals']) == (19, 12)
    assert [center['atom'] for center in ch3s['orbital_selection']['centers']] \
        == [1, 5]
    assert [center['atom'] for center in ch3ss['orbital_selection']['centers']] \
        == [1, 5, 6]
    assert active_space_feasibility(ch3s, [2, 4]) == {
        '2': {'alpha_electrons': 7, 'beta_electrons': 6,
              'determinants': 224},
        '4': {'alpha_electrons': 8, 'beta_electrons': 5,
              'determinants': 56},
    }
    assert active_space_feasibility(ch3ss, [2, 4]) == {
        '2': {'alpha_electrons': 10, 'beta_electrons': 9,
              'determinants': 14520},
        '4': {'alpha_electrons': 11, 'beta_electrons': 8,
              'determinants': 5940},
    }

    ch3s_rydberg = full_valence_active_space(
        ['C', 'H', 'H', 'H', 'S'],
        additional_centers=[{'atom': 1, 'orbitals': ['3s']}])
    assert (ch3s_rydberg['active_electrons'],
            ch3s_rydberg['active_orbitals']) == (13, 9)
    assert ch3s_rydberg['orbital_selection']['centers'][0] == {
        'atom': 1, 'orbitals': ['2s', '2p', '3s']}


def test_radical_pair_is_a_guarded_two_doublet_cas22_shortcut():
    parameters = _parameters(
        vrc_tst_active_space='radical_pair',
        vrc_tst_sample_method='caspt2-f12b(9,7)',
        vrc_tst_sample_basis='avtz-f12',
        vrc_tst_scan_method='caspt2-f12b(9,7)',
        vrc_tst_scan_basis='avtz-f12',
        vrc_tst_spin_splitting={'fragment_multiplicities': [2, 2]},
    )
    space = active_space_settings(parameters, ['C', 'C', *['H'] * 6])
    assert space == {
        'active_electrons': 2,
        'active_orbitals': 2,
        'orbital_selection': None,
        'selection': 'two_doublet_radical_pair',
        'validation': 'accepted_products_must_be_two_doublets',
    }
    validate_vrc_levels(parameters)
    model = resolved_state_model(
        parameters, ['C', 'C', *['H'] * 6], 0,
        spin_splitting_settings(parameters))
    assert 'orbital_selection' not in model
    assert [entry['multiplicity'] for entry in model['state_average']] == [1, 3]
    assert expected_levels(
        parameters, symbols=['C', 'C', *['H'] * 6]
    )['sampling']['method'] == 'caspt2-f12b(2,2)'

    with pytest.raises(ValueError, match='two doublet'):
        validate_vrc_levels({
            **parameters,
            'vrc_tst_spin_splitting': {'fragment_multiplicities': [2, 3]},
        })
    with pytest.raises(ValueError, match='spin_splitting'):
        validate_vrc_levels({
            **parameters,
            'vrc_tst_correction_model': 'direct',
        })

    holder = VTS.__new__(VTS)
    holder.par = parameters
    accepted = SimpleNamespace(
        products=[SimpleNamespace(mult=2, atom=['C', 'H', 'H', 'H']),
                  SimpleNamespace(mult=2, atom=['C', 'H', 'H', 'H'])],
        species=SimpleNamespace(atom=['C', 'C', *['H'] * 6],
                                charge=0, mult=1))
    settings, _ = VTS._spin_context(holder, accepted)
    assert (settings['target_multiplicity'],
            settings['reference_multiplicity']) == (1, 3)
    rejected = SimpleNamespace(
        products=[SimpleNamespace(mult=2, atom=['C', 'H', 'H', 'H']),
                  SimpleNamespace(mult=3, atom=['S'])],
        species=SimpleNamespace(atom=['C', 'H', 'H', 'H', 'S'],
                                charge=0, mult=2))
    with pytest.raises(ValueError, match='two doublet fragments'):
        VTS._spin_context(holder, rejected)


@pytest.mark.parametrize('stem', ['point', 'Long_Sulfur_Channel_' * 5])
def test_hash_verified_molpro_orbitals_are_staged_reproducibly(tmp_path, stem):
    source = tmp_path / 'expert.wfu'
    source.write_bytes(b'molpro orbital record fixture')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    model = validate_vrc_state_model({
        'orbital_selection': {
            'mode': 'molpro_record',
            'path': str(source),
            'record': '2140.2',
            'sha256': digest,
        },
        'state_average': [
            {'multiplicity': 1, 'roots': 1},
            {'multiplicity': 3, 'roots': 1},
        ],
    })
    staged = stage_vrc_state_model(model, tmp_path / 'calculation', stem)
    staged_path = Path(staged['orbital_selection']['path'])
    assert staged_path.is_absolute()
    assert len(staged_path.name) <= 31 and staged_path.name.islower()
    assert stage_vrc_state_model(model, tmp_path / 'calculation', stem) == staged
    other = stage_vrc_state_model(model, tmp_path / 'calculation', stem + '_other')
    assert other['orbital_selection']['path'] != str(staged_path)
    assert staged_path.read_bytes() == source.read_bytes()
    program, _ = _vrc_multireference_method(
        'caspt2-f12b(2,2)', 18, spin=0, state_model=staged)
    assert f'file,2,{staged_path.name}\n' in program
    assert str(staged_path.parent) not in program
    assert 'start,2140.2' in program
    assert 'occ,10\nclosed,8\n' in program
    assert 'diab,2140.2' in program
    automatic = automatic_vrc_state_model(
        ['C', 'C', *['H'] * 6], 0, [1, 3],
        active_space={
            'active_electrons': 2,
            'active_orbitals': 2,
            'orbital_selection': None,
        })
    assert 'orbital_selection' not in automatic

    source.write_bytes(b'changed')
    with pytest.raises(ValueError, match='SHA-256'):
        validate_vrc_state_model(model)


def test_selected_space_can_load_a_relative_reconnaissance_record(tmp_path):
    wfu = tmp_path / 'reviewed.wfu'
    wfu.write_bytes(b'reviewed orbitals')
    digest = hashlib.sha256(wfu.read_bytes()).hexdigest()
    payload = {
        'schema': 1,
        'active_space': {
            'mode': 'selected', 'active_electrons': 5,
            'active_orbitals': 4, 'source': 'reviewed C-S manifold',
        },
        'state_model': {
            'orbital_selection': {
                'mode': 'molpro_record', 'path': str(wfu),
                'record': '2140.2', 'sha256': digest,
            },
            'state_average': [
                {'multiplicity': 2, 'roots': 3},
                {'multiplicity': 4, 'roots': 3},
            ],
        },
    }
    (tmp_path / 'state.json').write_text(json.dumps(payload))
    (tmp_path / 'input.json').write_text('{}')
    parameters = _parameters(
        input_file=str(tmp_path / 'input.json'),
        vrc_tst_active_space={
            'mode': 'selected', 'active_electrons': 5,
            'active_orbitals': 4, 'source': 'reviewed C-S manifold',
        },
        vrc_tst_state_model_file='state.json',
        vrc_tst_sample_method='caspt2-f12b(2,2)',
        vrc_tst_sample_basis='avtz-f12',
        vrc_tst_scan_method='caspt2-f12b(2,2)',
        vrc_tst_scan_basis='avtz-f12',
        vrc_tst_spin_splitting={
            'target_multiplicity': 2,
            'reference_multiplicity': 4,
            'fragment_multiplicities': [2, 3],
        },
    )
    validate_vrc_levels(parameters)
    settings = spin_splitting_settings(parameters)
    model = resolved_state_model(
        parameters, ['C', 'S', 'H', 'H', 'H'], 0, settings,
        roots={2: 1, 4: 1})
    assert model['orbital_selection']['sha256'] == digest
    assert model['state_average'] == [
        {'multiplicity': 2, 'roots': 1, 'weights': [1.0]},
        {'multiplicity': 4, 'roots': 1, 'weights': [1.0]},
    ]
    assert expected_levels(
        parameters, symbols=['C', 'S', 'H', 'H', 'H'], settings=settings
    )['sampling']['method'] == 'caspt2-f12b(5,4)'


def test_external_state_model_selects_a_basis_specific_wfu(tmp_path):
    models = {}
    for basis in ('avtz-f12', 'avqz-f12'):
        wfu = tmp_path / f'{basis}.wfu'
        wfu.write_bytes(basis.encode())
        models[basis] = {
            'orbital_selection': {
                'mode': 'molpro_record', 'path': str(wfu),
                'record': '2140.2',
                'sha256': hashlib.sha256(wfu.read_bytes()).hexdigest(),
            },
            'state_average': [{'multiplicity': 2, 'roots': 1}],
        }
    (tmp_path / 'state.json').write_text(json.dumps({
        'schema': 2,
        'active_space': {
            'mode': 'selected', 'active_electrons': 5,
            'active_orbitals': 4, 'source': 'reviewed manifold',
        },
        'basis_models': models,
    }))
    parameters = _parameters(
        input_file=str(tmp_path / 'input.json'),
        vrc_tst_active_space={
            'mode': 'selected', 'active_electrons': 5,
            'active_orbitals': 4, 'source': 'reviewed manifold',
        },
        vrc_tst_state_model_file='state.json',
        vrc_tst_sample_method='caspt2-f12b(5,4)',
        vrc_tst_sample_basis=['avtz-f12', 'avqz-f12'],
        vrc_tst_scan_method='caspt2-f12b(5,4)',
        vrc_tst_scan_basis='avtz-f12',
        vrc_tst_spin_splitting={
            'target_multiplicity': 2,
            'reference_multiplicity': 4,
            'fragment_multiplicities': [2, 3],
            'reference': {'method': 'uccsd(t)-f12b',
                          'basis': ['avtz-f12', 'avqz-f12']},
            'gap_levels': [
                {'label': 'caspt2_f12', 'method': 'caspt2-f12b(5,4)',
                 'basis': ['avtz-f12', 'avqz-f12']},
            ],
            'weighting': {'model': 'arithmetic_mean'},
        },
    )
    components = correction_components(
        parameters, symbols=['C', 'S', 'H', 'H', 'H'], settings=None,
        state_roots={2: 1, 4: 1})
    selected = {
        item['basis']: Path(item['state_model']['orbital_selection']['path']).name
        for item in components
    }
    assert selected == {
        'avtz-f12': 'avtz-f12.wfu',
        'avqz-f12': 'avqz-f12.wfu',
    }
    with pytest.raises(ValueError, match='no orbital record'):
        resolved_state_model(
            parameters, ['C', 'S', 'H', 'H', 'H'], 0,
            spin_splitting_settings(parameters), basis='av5z-f12')


def test_anl1_f12_vrc_preset_is_augmented_and_preserves_overrides():
    parameters = {'vrc_tst_accuracy_preset': 'anl1-f12',
                  'vrc_tst_scan_basis': 'aug-cc-pVQZ'}
    apply_vrc_accuracy_preset(parameters, {'vrc_tst_accuracy_preset',
                                           'vrc_tst_scan_basis'})
    assert parameters['vrc_tst_active_space'] == 'full_valence'
    assert parameters['vrc_tst_scan_basis'] == 'aug-cc-pVQZ'
    assert parameters['vrc_tst_scan_gradient'] == 'numerical'
    assert parameters['vrc_tst_sample_basis'] == 'avtz-f12'
    assert parameters['vrc_tst_spin_splitting']['reference']['basis'] == 'avtz-f12'
    assert parameters['vrc_tst_high_basis'] == 'avtz-f12'
    assert all(level['basis'] == 'avtz-f12'
               for level in parameters['vrc_tst_spin_splitting']['gap_levels'])


def test_sampling_surface_accepts_an_adjacent_cbs_pair():
    parameters = _parameters(
        vrc_tst_sample_method='caspt2-f12b(9,7)',
        vrc_tst_sample_basis=['avtz-f12', 'avqz-f12'],
        vrc_tst_sample_extrapolation_power=3.7)
    settings = sampling_settings(parameters)
    assert settings == {
        'method': 'caspt2-f12b(9,7)',
        'basis': ['avtz-f12', 'avqz-f12'],
        'extrapolation_power': 3.7,
    }
    components = correction_components(parameters)
    assert [item['role'] for item in components[:2]] == [
        'sampling:0', 'sampling:1']
    curve, metadata = assemble_sampling_reference(
        [2.5, 30.],
        {'sampling:0': [-10., -9.], 'sampling:1': [-10.1, -9.1]},
        settings)
    assert curve == pytest.approx([
        two_point_cbs(-10., -10.1, upper_cardinal=4),
        two_point_cbs(-9., -9.1, upper_cardinal=4)])
    assert metadata['power'] == pytest.approx(3.7)


def test_spin_splitting_extrapolates_triplet_and_both_gap_states():
    settings = spin_splitting_settings(_parameters())
    distances = [2.5, 5.0]
    energies = {
        'triplet_reference:0': [-10.0, -9.0],
        'triplet_reference:1': [-10.1, -9.1],
        'gap:caspt2_f12:singlet:0': [-5.30, -4.20],
        'gap:caspt2_f12:singlet:1': [-5.35, -4.25],
        'gap:caspt2_f12:triplet:0': [-5.00, -4.00],
        'gap:caspt2_f12:triplet:1': [-5.02, -4.02],
        'gap:mrci_f12_q:singlet:0': [-6.40, -5.30],
        'gap:mrci_f12_q:singlet:1': [-6.45, -5.35],
        'gap:mrci_f12_q:triplet:0': [-6.00, -5.00],
        'gap:mrci_f12_q:triplet:1': [-6.03, -5.03],
    }
    result, metadata = assemble_spin_splitting_reference(
        distances, energies, settings)
    triplet = np.asarray([
        two_point_cbs(-10.0, -10.1, upper_cardinal=4),
        two_point_cbs(-9.0, -9.1, upper_cardinal=4),
    ])
    cas_gap = np.asarray([
        two_point_cbs(-5.30, -5.35, upper_cardinal=4)
        - two_point_cbs(-5.00, -5.02, upper_cardinal=4),
        two_point_cbs(-4.20, -4.25, upper_cardinal=4)
        - two_point_cbs(-4.00, -4.02, upper_cardinal=4),
    ])
    mrci_gap = np.asarray([
        two_point_cbs(-6.40, -6.45, upper_cardinal=4)
        - two_point_cbs(-6.00, -6.03, upper_cardinal=4),
        two_point_cbs(-5.30, -5.35, upper_cardinal=4)
        - two_point_cbs(-5.00, -5.03, upper_cardinal=4),
    ])
    expected = triplet + .5 * (cas_gap + mrci_gap)
    assert result == pytest.approx(expected)
    assert metadata['triplet_cbs']['power'] == pytest.approx(3.7)
    assert set(metadata['gap_cbs']) == {'caspt2_f12', 'mrci_f12_q'}


def test_arithmetic_mean_reproduces_the_published_two_gap_average():
    settings = spin_splitting_settings(_parameters(
        vrc_tst_spin_splitting={
            'weighting': {'model': 'arithmetic_mean'},
            'triplet_reference': {'basis': 'avtz-f12'},
            'gap_levels': [
                {'label': 'caspt2_f12', 'method': 'caspt2-f12(2,2)',
                 'basis': 'avtz-f12', 'shift': .2, 'ipea': .25},
                {'label': 'mrci_f12_q', 'method': 'mrci-f12+q(2,2)',
                 'basis': 'avtz-f12'},
            ],
        }))
    energies = {
        'triplet_reference': [-10.],
        'gap:caspt2_f12:singlet': [-5.3],
        'gap:caspt2_f12:triplet': [-5.0],
        'gap:mrci_f12_q:singlet': [-6.4],
        'gap:mrci_f12_q:triplet': [-6.0],
    }
    result, metadata = assemble_spin_splitting_reference([3.], energies, settings)
    assert result == pytest.approx([-10.35])
    assert metadata['method_sensitivity']['absolute_spread'] == pytest.approx([.1])
    assert metadata['method_sensitivity']['asymptote_referenced_spread'] == pytest.approx([0.])
    assert 'not a statistical error bar' in metadata['method_sensitivity']['interpretation']


def test_asymptote_is_input_driven_and_outside_the_scan_grid():
    from kinbot.vrc_correction import asymptote_distance
    assert asymptote_distance({'vrc_tst_asymptote_distance': 75.,
                               'vrc_tst_scan_points': [3., 45.]}) == 75.
    with pytest.raises(ValueError, match='exceed every'):
        asymptote_distance({'vrc_tst_scan_points': [30.]})
    with pytest.raises(ValueError, match='positive and finite'):
        asymptote_distance({'vrc_tst_asymptote_distance': True})


def test_molpro_spin_splitting_program_syntax_and_energy_variables(tmp_path):
    cas, cas_energy = _vrc_multireference_method(
        'caspt2-f12(2,2)', 18, spin=0, shift=.2, ipea=.25)
    assert '{rs2-f12,shift=0.2,ipea=0.25,thrden=1.d-8,thrvar=1.d-10}' in cas
    assert 'KBVRCSEL=energy(1)' in cas
    assert cas_energy == 'KBVRCSEL'
    mrci, mrci_energy = _vrc_multireference_method(
        'mrci-f12+q(2,2)', 18, spin=2)
    assert 'wf,18,1,2' in mrci
    assert '{mrci-f12,scalf12=1,cabs_singles=1,thrden=1.d-8,thrvar=1.d-10;option,cluster=0;}' in mrci
    assert 'KBVRCVAR=energy(1)' in mrci
    assert 'KBVRCSEL=energd(1)' in mrci
    assert mrci_energy == 'KBVRCSEL'

    for method in ('caspt2-f12a(2,2)', 'mrci-f12a+q(2,2)'):
        with pytest.raises(ValueError, match='no F12A/B variants'):
            _vrc_multireference_method(method, 18)
    legacy, _ = _vrc_multireference_method('caspt2-f12b(2,2)', 18, shift=.2, ipea=.25)
    assert legacy == cas

    contract = vrc_multireference_energy_contract('mrci-f12+q(2,2)')
    assert contract == {
        'method': 'MRCI-F12+Q',
        'active_electrons': 2,
        'active_orbitals': 2,
        'f12_variant': None,
        'explicitly_correlated': True,
        'legacy_alias': None,
        'energy_variable': 'ENERGD(1)',
        'cabs_singles_included': True,
        'davidson_correction': 'relaxed_reference',
        'molpro_cluster_option': 0,
    }

    species = SimpleNamespace(
        atom=['C', 'H', 'H', 'H', 'H'],
        geom=np.zeros((5, 3)), natom=5, name='methane', charge=0, mult=1)
    parameters = {
        'vrc_tst_scan_molpro_tpl': '',
        'vrc_tst_scan_molpro_key': 'MYENERGY',
        'vrc_tst_sample_method': 'caspt2(2,2)',
        'vrc_tst_sample_basis': 'vdz',
        'vrc_tst_high_method': 'mrci+q(2,2)',
        'vrc_tst_high_basis': 'avtz',
    }
    previous = Path.cwd()
    os.chdir(tmp_path)
    try:
        Path('vrctst/molpro').mkdir(parents=True)
        Molpro(species, parameters).create_molpro_input(
            name='triplet', VTS=True,
            vrc_level={'method': 'uccsd(t)-f12b', 'basis': 'avqz-f12',
                       'multiplicity': 3})
        generated = Path('vrctst/molpro/triplet.inp').read_text()
    finally:
        os.chdir(previous)
    assert 'basis = avqz-f12' in generated
    assert '{rhf;wf,10,1,2}' in generated
    assert 'uccsd(t)-f12b' in generated
    assert 'MYENERGY = energy' in generated


def test_molpro_reads_audited_davidson_f12b_assignment(tmp_path):
    species = SimpleNamespace(
        atom=['C', 'H', 'H', 'H', 'H'], geom=np.zeros((5, 3)), natom=5,
        name='methane', charge=0, mult=1)
    parameters = {
        'vrc_tst_scan_molpro_tpl': '',
        'vrc_tst_scan_molpro_key': 'MYENERGY',
        'vrc_tst_sample_method': 'caspt2(2,2)',
        'vrc_tst_sample_basis': 'vdz',
        'vrc_tst_high_method': 'mrci+q(2,2)',
        'vrc_tst_high_basis': 'avtz',
    }
    previous = Path.cwd()
    os.chdir(tmp_path)
    try:
        directory = Path('vrctst/molpro')
        directory.mkdir(parents=True)
        molpro = Molpro(species, parameters)
        molpro.create_molpro_input(
            name='mrci', VTS=True,
            vrc_level={'method': 'mrci-f12b+q(2,2)',
                       'basis': 'avtz-f12', 'multiplicity': 1})
        generated = (directory / 'mrci.inp').read_text()
        assert 'KBVRCSEL=energd(1)' in generated
        assert 'MYENERGY = KBVRCSEL' in generated
        assert 'put,xml,mrci.orbitals.xml' in generated
        assert 'put,molden' not in generated
        assert generated.index('maxit,50;}') < generated.index('put,xml,') < generated.index('{mrci-f12')

        # Values follow the structure of Molpro's published MRCI-F12 output:
        # ENERGY is variational; ENERGD is the CABS-inclusive +Q result.
        (directory / 'mrci.out').write_text(
            ' !CI-F12 STATE 1.1 Energy            -114.538939848235\n'
            ' Cluster corrected energies          -114.55264645 '
            '(Davidson, relaxed reference)\n'
            ' SETTING KBVRCVAR = -114.545957200000 AU\n'
            ' SETTING KBVRCSEL = -114.559663800000 AU\n'
            ' SETTING MYENERGY = -114.559663800000 AU\n'
            ' Molpro calculation terminated\n')
        assert molpro.get_molpro_energy(
            'MYENERGY', name='mrci', VTS=True) == (
                1, pytest.approx(-114.5596638))

        good = (directory / 'mrci.out').read_text()
        for broken in (good.replace('Molpro calculation terminated', ''),
                       'PROGRAM * MULTI\nNO CONVERGENCE! MAXIMAL NUMBER OF ITERATIONS REACHED!\n' + good):
            (directory / 'mrci.out').write_text(broken)
            assert molpro.get_molpro_energy('MYENERGY', name='mrci', VTS=True) == (0, -1)

        (directory / 'mrci.out').write_text(
            ' SETTING KBVRCSEL = -114.559663800000 AU\n'
            ' SETTING MYENERGY = -114.545957200000 AU\n'
            ' Molpro calculation terminated\n')
        assert molpro.get_molpro_energy(
            'MYENERGY', name='mrci', VTS=True) == (0, -1)

        molpro.create_molpro_input(
            name='caspt2', VTS=True,
            vrc_level={'method': 'caspt2-f12b(2,2)',
                       'basis': 'avtz-f12', 'multiplicity': 1,
                       'shift': .2, 'ipea': .25})
        assert 'KBVRCSEL=energy(1)' in (
            directory / 'caspt2.inp').read_text()
        (directory / 'caspt2.out').write_text(
            ' SETTING F12A     = -39.100000000000 AU\n'
            ' SETTING F12B     = -39.200000000000 AU\n'
            ' SETTING KBVRCSEL = -39.200000000000 AU\n'
            ' SETTING MYENERGY = -39.200000000000 AU\n'
            ' Molpro calculation terminated\n')
        assert molpro.get_molpro_energy(
            'MYENERGY', name='caspt2', VTS=True) == (
                1, pytest.approx(-39.2))
    finally:
        os.chdir(previous)


def test_single_reference_anchor_cannot_silently_ignore_an_excited_cas_root():
    from kinbot.molpro import validate_vrc_single_reference_root
    model = {'selected_roots': {2: 1, 4: 3}}
    with pytest.raises(ValueError, match='cannot select excited CAS root 3'):
        validate_vrc_single_reference_root('uccsd(t)-f12b', model, 4)
    validate_vrc_single_reference_root('mrci-f12b+q(13,8)', model, 4)
    validate_vrc_single_reference_root('uccsd(t)-f12b', {'selected_roots': {4: 1}}, 4)


@pytest.mark.parametrize('name', ['ion', 'general_channel_name_' * 5])
def test_vrc_avastart_and_native_export_follow_the_requested_charged_state(tmp_path, monkeypatch, name):
    monkeypatch.chdir(tmp_path)
    Path('vrctst/molpro').mkdir(parents=True)
    species = SimpleNamespace(atom=['Li', 'H'], geom=np.asarray([[0., 0., 0.], [1.6, 0., 0.]]),
        natom=2, name='ion', charge=1, mult=2)
    model = automatic_vrc_state_model(['Li', 'H'], 1, [2],
        active_space={'active_electrons': 1, 'active_orbitals': 4,
            'orbital_selection': {'mode': 'avas', 'threshold': .1, 'minao_basis': 'MINAO',
                'centers': [{'atom': 1, 'orbitals': ['2s', '2p']}]}})
    parameters = {'vrc_tst_scan_molpro_tpl': '', 'vrc_tst_scan_molpro_key': 'MYENERGY'}
    Molpro(species, parameters).create_molpro_input(name=name, VTS=True,
        vrc_level={'method': 'caspt2-f12b(1,4)', 'basis': 'av5z-f12', 'state_model': model})
    deck = Path(f'vrctst/molpro/{name}.inp').read_text()
    assert len(deck.splitlines()[0].split(',', 1)[1]) <= 32
    assert 'set,charge=1\nset,spin=1' in deck
    assert deck.index('{avas,') < deck.index('maxit,50;}') < deck.index('put,xml,') < deck.index('{rs2-f12')
    assert 'put,molden' not in deck


@pytest.mark.parametrize('title', [
    'ch3ss_orbital_reconnaissance_pilot orbital reconnaissance: avd__cs_r_1p636013716',
    'arbitrary user input ' * 100,
    '\u03c3 orbital diagnostic',
    ' ' * 100 + 'short label',
])
def test_native_xml_titles_are_short_deterministic_and_leave_chemistry_intact(title):
    body = '\r\nbasis=avtz-f12\r\n{rhf;wf,23,1,1}\r\nput,xml,orbitals.xml\r\n'
    original = '***,' + title + body
    shortened = compact_molpro_titles(original)
    header = shortened.split('\r\n', 1)[0]
    assert header.isascii() and len(header.split(',', 1)[1]) <= 32
    assert shortened == compact_molpro_titles(original)
    assert shortened == compact_molpro_titles(shortened)
    assert shortened[shortened.index('\r\n'):] == body
    assert shortened != compact_molpro_titles('***,' + title + 'different' + body)


def test_short_legacy_molpro_titles_are_unchanged():
    deck = '***,KinBot\n{rhf}\n***, second job\n{rhf}\n'
    assert compact_molpro_titles(deck) == deck


def test_spin_splitting_rejects_nonadjacent_cbs_pairs_and_bad_options():
    with pytest.raises(ValueError, match='adjacent'):
        spin_splitting_settings(_parameters(vrc_tst_spin_splitting={
            'triplet_reference': {'basis': ['avdz-f12', 'avqz-f12']}}))
    with pytest.raises(ValueError, match='shift'):
        spin_splitting_settings(_parameters(vrc_tst_spin_splitting={
            'gap_levels': [
                {'label': 'a', 'method': 'caspt2-f12(2,2)',
                 'basis': 'avtz-f12', 'shift': 'bad'},
                {'label': 'b', 'method': 'mrci-f12+q(2,2)',
                 'basis': 'avtz-f12'},
            ]}))
