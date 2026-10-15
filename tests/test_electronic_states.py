"""Electronic-state contracts shared by graph routing and MESS."""

from types import SimpleNamespace

import pytest

from kinbot.electronic_states import (
    apply_electronic_state_model, electronic_levels_for_species,
    electronic_state_models, mess_electronic_levels,
)
from kinbot.mess import MESS


def _parameters():
    return {
        'electronic_state_models': [{
            'smiles': '[S]',
            'charge': 0,
            'ground_multiplicity': 3,
            'source': 'test fine-structure model',
            'levels': [
                {'energy_cm_inverse': 0., 'multiplicity': 3,
                 'orbital_degeneracy': 1},
                {'energy_cm_inverse': 396.06, 'multiplicity': 3,
                 'orbital_degeneracy': 1},
            ],
        }],
    }


def test_species_model_controls_stable_mess_levels_but_not_reactive_surface():
    species = SimpleNamespace(smiles='[S]', charge=0, mult=3)
    model = apply_electronic_state_model(species, _parameters())
    assert model['source'] == 'test fine-structure model'
    assert electronic_levels_for_species(species, _parameters()) == [
        {'energy_cm_inverse': 0., 'multiplicity': 3,
         'orbital_degeneracy': 1, 'degeneracy': 3},
        {'energy_cm_inverse': 396.06, 'multiplicity': 3,
         'orbital_degeneracy': 1, 'degeneracy': 3},
    ]
    assert electronic_levels_for_species(
        species, _parameters(), ground_only=True) == [
            {'energy_cm_inverse': 0., 'degeneracy': 3}]

    count, lines = mess_electronic_levels(species, _parameters())
    assert count == 2
    assert lines.splitlines() == ['0 3', '396.06 3']


def test_mess_atomic_and_rrho_blocks_render_explicit_counts():
    species = SimpleNamespace(smiles='[S]', charge=0, mult=3)
    writer = MESS.__new__(MESS)
    writer.par = _parameters()
    fields = writer._electronic_fields(species, '      ')
    assert fields == {
        'nelec': 2,
        'electronic_levels': '      0 3\n      396.06 3',
    }
    levels = electronic_levels_for_species(species, writer.par)
    text = writer._atomic_partition_function_input(
        'sulfur', 'S', levels, [300., 1000.])
    assert 'ElectronicLevels[1/cm] 2' in text
    assert '0 3' in text
    assert '396.06 3' in text


def test_no_configured_model_preserves_legacy_internal_smiles_behavior():
    species = SimpleNamespace(
        smiles='[H][C]([H])([O]([H])[O])[H]', charge=0, mult=2)
    assert electronic_levels_for_species(species, {}) == [
        {'energy_cm_inverse': 0., 'degeneracy': 2}]


def test_electronic_degeneracy_must_equal_spin_times_orbital_factor():
    parameters = _parameters()
    parameters['electronic_state_models'][0]['levels'][0]['degeneracy'] = 2
    with pytest.raises(ValueError, match='multiplicity times'):
        electronic_state_models(parameters)


def test_resolved_j_levels_and_centroid_correction_are_general():
    from ase.units import Hartree, invcm
    from kinbot.electronic_states import state_average_roots
    parameters = _parameters()
    record = parameters['electronic_state_models'][0]
    record.update(term_orbital_degeneracy=3, spin_free_reference='term_centroid')
    record['levels'] = [dict(energy_cm_inverse=energy, multiplicity=3, degeneracy=deg)
                        for energy, deg in [(0., 5), (396.055, 3), (573.640, 1)]]
    species = SimpleNamespace(smiles='[S]', charge=0, mult=1)
    model = apply_electronic_state_model(species, parameters)
    assert species.mult == 3
    assert [level['degeneracy'] for level in model['levels']] == [5, 3, 1]
    assert species.spin_orbit_correction['value_hartree'] == pytest.approx(
        -(3 * 396.055 + 573.640) / 9 * invcm / Hartree)
    assert state_average_roots([species], [2, 4], parameters) == {2: 3, 4: 3}
    record.pop('spin_free_reference')
    apply_electronic_state_model(species, parameters)
    assert not hasattr(species, 'spin_orbit_correction')


def test_alias_updates_lewis_identity_without_aliasing_different_connectivity():
    parameters = {'electronic_state_models': [{
        'smiles': '[S][S]', 'smiles_aliases': ['S=S'], 'charge': 0,
        'ground_multiplicity': 3,
        'levels': [{'energy_cm_inverse': 0, 'multiplicity': 3}],
    }]}
    species = SimpleNamespace(smiles='S=S', charge=0, mult=1)
    apply_electronic_state_model(species, parameters)
    assert species.mult == 3 and species.smiles == '[S][S]'
    parameters['electronic_state_models'][0]['smiles_aliases'] = ['S.S']
    with pytest.raises(ValueError, match='same atoms'):
        electronic_state_models(parameters)


def test_spin_parity_and_incomplete_centroid_fail_early():
    parameters = _parameters()
    parameters['electronic_state_models'][0]['levels'][1]['multiplicity'] = 2
    with pytest.raises(ValueError, match='electron count'):
        electronic_state_models(parameters)
    parameters = _parameters()
    parameters['electronic_state_models'][0].update(
        term_orbital_degeneracy=3, spin_free_reference='term_centroid')
    with pytest.raises(ValueError, match='complete ground-term'):
        electronic_state_models(parameters)


def test_spatial_term_declaration_does_not_invent_unresolved_mess_levels():
    from kinbot.electronic_states import state_average_roots
    parameters = _parameters()
    parameters['electronic_state_models'].append({
        'smiles': '[CH3][S]', 'charge': 0, 'ground_multiplicity': 2,
        'term_orbital_degeneracy': 2,
    })
    fragments = [SimpleNamespace(smiles=smiles, charge=0, mult=mult)
                 for smiles, mult in [('[CH3][S]', 2), ('[S]', 3)]]
    parameters['electronic_state_models'][0]['term_orbital_degeneracy'] = 3
    assert state_average_roots(fragments, [2, 4], parameters) == {2: 6, 4: 6}
    assert electronic_levels_for_species(fragments[0], parameters) == [
        {'energy_cm_inverse': 0., 'degeneracy': 2}]
