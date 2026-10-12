"""Spin-splitting VRC correction arithmetic and input contracts."""

import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from kinbot.anl.extrapolation import two_point_cbs
from kinbot.molpro import Molpro, _vrc_multireference_method
from kinbot.vrc_correction import (
    assemble_spin_splitting_reference, correction_components,
    expected_levels, spin_splitting_settings, spin_splitting_weights,
)


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
    components = correction_components(_parameters())
    assert len(components) == 11
    assert {component['multiplicity'] for component in components[1:]} == {
        1, 3}
    assert expected_levels(_parameters())['trusted_correction']['model'] == \
        'spin_splitting'


def test_distance_switch_prefers_mrci_short_and_caspt2_long_range():
    settings = spin_splitting_settings(_parameters())
    weights = spin_splitting_weights([2., 2.5, 3.125, 3.75, 4.375, 5., 8.],
                                     settings)
    assert weights['caspt2_f12'] == pytest.approx(
        [.2, .2, .4, .5, .6, .8, .8])
    assert weights['mrci_f12_q'] == pytest.approx(
        [.8, .8, .6, .5, .4, .2, .2])
    assert np.asarray(weights['caspt2_f12']) + np.asarray(
        weights['mrci_f12_q']) == pytest.approx(np.ones(7))


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
    expected = triplet + np.asarray([.2, .8]) * cas_gap \
        + np.asarray([.8, .2]) * mrci_gap
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
    result, _ = assemble_spin_splitting_reference([3.], energies, settings)
    assert result == pytest.approx([-10.35])


def test_molpro_spin_splitting_program_syntax_and_energy_variables(tmp_path):
    cas, cas_energy = _vrc_multireference_method(
        'caspt2-f12(2,2)', 18, spin=0, shift=.2, ipea=.25)
    assert '{rs2-f12,shift=0.2,ipea=0.25}' in cas
    assert cas_energy == 'energy(1)'
    mrci, mrci_energy = _vrc_multireference_method(
        'mrci-f12+q(2,2)', 18, spin=2)
    assert 'wf,18,1,2' in mrci
    assert '{mrci-f12}' in mrci
    assert mrci_energy == 'energd(1)'

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
