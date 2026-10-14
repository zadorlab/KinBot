"""High-accuracy sulfur examples exercise open-shell VRC input contracts."""

from pathlib import Path

import pytest

from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint
from kinbot.vrc_correction import (correction_components,
                                   coupled_multiplicities,
                                   spin_splitting_settings)


@pytest.mark.parametrize(
    'name,smiles,active', [
        ('ch3s', '[CH3][S]', '(13,8)'),
        ('ch3ss', 'C[S][S]', '(9,7)'),
    ])
def test_sulfur_examples_are_hands_off_anl1_f12_open_shell_inputs(
        name, smiles, active):
    path = (Path(__file__).parents[1] / 'examples' / 'anl' /
            f'{name}_anl1_f12' / f'{name}.json')
    parsed = Parameters(str(path), show_warnings=False)
    par = parsed.par
    assert par['smiles'] == smiles
    assert par['mult'] == 2
    assert par['composite_method'] == 'ANL1-F12'
    assert par['anl_f12_augmented'] == 1
    assert par['cbh_rung'] == 'auto'
    assert par['families'] == ['all']
    assert par['selected_reactions'] == {}
    assert par['vrc_tst_scan_gradient'] == 'numerical'
    assert active in par['vrc_tst_scan_method']
    assert par['vrc_tst_sample_basis'] == ['avtz-f12', 'avqz-f12']
    settings = spin_splitting_settings(par)
    assert settings['fragment_multiplicities'] == [2, 3]
    assert coupled_multiplicities(*settings['fragment_multiplicities']) == (
        2, 4)
    assert settings['target_multiplicity'] == 2
    assert settings['reference_multiplicity'] == 4
    assert settings['reference']['basis'] == ['vqz-f12', 'v5z-f12']
    assert all(level['basis'] == ['vqz-f12', 'v5z-f12']
               for level in settings['gap_levels'])
    assert [item['role'] for item in correction_components(par)[:2]] == [
        'sampling:0', 'sampling:1']

    species = StationaryPoint(name, 0, 2, smiles=smiles)
    species.characterize()
    assert species.mult == 2
    assert 'S' in species.atom
