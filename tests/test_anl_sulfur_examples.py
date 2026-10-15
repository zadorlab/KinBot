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
        ('ch3s', '[CH3][S]', (13, 8)),
        ('ch3ss', 'C[S][S]', (9, 7)),
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
    assert par['vrc_tst_state_model_file'] == f'{name}_state_model.json'
    assert par['vrc_tst_surface_model_file'] == f'{name}_surface_model.json'
    assert par['homolytic_bonds'] == {'initial': ([[0, 1]] if name == 'ch3s' else [[0, 1], [1, 2]])}
    assert par['vrc_tst_sample_basis'] == 'avtz-f12'
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
    assert settings['reference']['basis'] == 'avtz-f12'
    assert all(level['basis'] == 'avtz-f12'
               for level in settings['gap_levels'])
    assert species.mult == 2
    assert 'S' in species.atom

    reconnaissance = path.with_name('orbital_reconnaissance.json')
    assert reconnaissance.is_file()


def test_validation_recipe_uses_requested_sampling_and_scan_grid():
    root = Path(__file__).parents[1] / 'examples/anl/orbital_reviews'
    spec = importlib.util.spec_from_file_location('example_validation', root / 'prepare_validation.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    recipe = json.loads((root / 'validation_candidates.json').read_text())
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
            assert check['method'] == f'caspt2-f12b({row["active_electrons"]},{row["active_orbitals"]})'
            assert check['shift'] == parameters['vrc_tst_scan_shift']
            assert check['ipea'] == parameters['vrc_tst_scan_ipea']
            assert row['region'] == 'diagnostic'  # not an inferred inner/outer selection
