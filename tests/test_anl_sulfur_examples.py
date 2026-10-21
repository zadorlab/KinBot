"""High-accuracy sulfur examples exercise open-shell VRC input contracts."""

from pathlib import Path
import importlib.util
import json

import pytest

from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint
from kinbot.vrc_correction import (active_space_settings,
                                   coupled_multiplicities,
                                   spin_splitting_settings)


@pytest.mark.parametrize(
    'name,smiles,active', [
        ('ch3s', '[CH3][S]', (5, 4)),
        ('ch3ss', 'C[S][S]', (7, 5)),
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
    assert par['vrc_tst_active_space']['mode'] == 'selected'
    assert par['vrc_tst_state_model_file'] == par['vrc_tst_surface_model_file'] == ''
    assert par['vrc_tst_surface_model']['channels']
    assert par['homolytic_bonds'] == {'initial': ([[0, 1]] if name == 'ch3s' else [[0, 1], [1, 2]])}
    assert par['vrc_tst_sample_basis'] == 'avdz-f12'
    assert par['vrc_tst_scan_basis'] == 'avdz-f12'
    assert par['vrc_tst_fragment_geometry'] == 'asymptotic_scan'
    assert par['vrc_tst_pivot_level'] == 'L2'
    assert par['vrc_tst_scan_shift'] == .1
    assert par['vrc_tst_scan_ipea'] == .25
    species = StationaryPoint(name, 0, 2, smiles=smiles)
    species.characterize()
    space = active_space_settings(par, species.atom, species.charge)
    assert (space['active_electrons'], space['active_orbitals']) == active
    assert space['selection'] == 'selected_orbital_subspace'
    settings = spin_splitting_settings(
        par, parent_multiplicity=2, fragment_multiplicities=[2, 3])
    assert settings['fragment_multiplicities'] == [2, 3]
    assert coupled_multiplicities(*settings['fragment_multiplicities']) == (
        2, 4)
    assert settings['target_multiplicity'] == 2
    assert settings['reference_multiplicity'] == 4
    assert settings['reference']['basis'] == ['avdz-f12', 'avtz-f12']
    assert settings['reference']['extrapolation_power'] == 3.7
    gap, = settings['gap_levels']
    assert gap['basis'] == ['avdz-f12', 'avtz-f12']
    assert gap['extrapolation_power'] == 3.7
    assert gap['shift'] == .1
    assert gap['ipea'] == .25
    assert species.mult == 2
    assert 'S' in species.atom

    reconnaissance = path.parent/'orbital_reviews/orbital_reconnaissance.json'
    assert reconnaissance.is_file()


def test_validation_recipe_preserves_comparison_level_and_requested_scan_grid():
    root = Path(__file__).parents[1] / 'examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history'
    spec = importlib.util.spec_from_file_location('example_validation', root / 'prepare_validation.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    recipe = json.loads((root / 'validation_candidates.json').read_text())
    assert recipe['comparison_level'] == {'method': 'caspt2-f12(2,2)', 'basis': 'avtz-f12'}
    assert recipe['comparison_settings'] == {'shift': .3, 'ipea': .25}
    definition = module.plan(recipe, root)
    assert len(definition['models']) == 39
    for case in recipe['cases']:
        rows = [m for m in definition['models'] if m['name'].startswith(case['label'] + '_')]
        assert len(rows) == 13
        parameters = Parameters(str((root / case['kinbot_input']).resolve())).par
        assert all(m['target_distance_angstrom'] in parameters['vrc_tst_scan_points']
                   for m in rows if 'target_distance_angstrom' in m)
        for row in rows:
            check = row['correlation_checks'][0]
            assert check['method'] == f'caspt2-f12({row["active_electrons"]},{row["active_orbitals"]})'
            assert check['shift'] == recipe['comparison_settings']['shift']
            assert check['ipea'] == recipe['comparison_settings']['ipea']
            assert row['region'] == 'diagnostic'  # not an inferred inner/outer selection


def test_cas75_addition_keeps_pi_pair_and_matches_the_existing_potential_stencil(monkeypatch):
    root = Path(__file__).parents[1] / 'examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history'
    monkeypatch.syspath_prepend(str(root))
    from prepare_cas75 import definition
    from prepare_validation import plan
    result = definition(root / 'cas75_candidates.json')['models']
    original = plan(json.loads((root / 'validation_candidates.json').read_text()), root)['models']
    reference = {m['name'].rsplit('_', 1)[-1]: m for m in original
                 if m['name'].startswith('ch3ss_cs_cas9_7_')}
    assert len(result) == 5
    for item in result:
        point = item['name'].rsplit('_', 1)[-1]
        assert item.get('target_distance_angstrom') == reference[point].get('target_distance_angstrom')
        assert (item['active_electrons'], item['active_orbitals']) == (7, 5)
        assert len(item['source_wfu_sha256']) == 64
        assert item['close_lowest'] == 6 and item['discard_highest'] == 1
        assert any(check['method'] == 'mrci-f12+q(7,5)' for check in item['correlation_checks'])
        assert {check.get('ipea') for check in item['correlation_checks']} == {None, 0., .25}
        if point == 'tilted':
            assert 'orbital_swaps' not in item
        else:
            assert item['source_directory'].endswith('/avt_f12__cs_r_30')
            # Source indices 12..23: close sigma #19, discard sigma* #23;
            # retain the pi/pi* quartet and the methyl radical (#22).
            order = list(range(12, 24))
            for a, b in item['orbital_swaps']:
                order[a-1], order[b-1] = order[b-1], order[a-1]
            assert set(order[6:-1]) == {17, 18, 20, 21, 22}
