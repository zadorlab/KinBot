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


@pytest.mark.parametrize('case', ['ch3s', 'ch3ss'])
def test_reviewed_molecular_corrections_preserve_cas_terms(case):
    import json
    from pathlib import Path
    from ase.units import invcm, Hartree
    from kinbot.electronic_states import state_average_roots
    path = Path(__file__).parents[1] / f'examples/anl/{case}_anl1_f12/{case}.json'
    parameters = json.loads(path.read_text())
    models = electronic_state_models(parameters)
    for smiles, gap in [('C[S]', 255.), ('[SH]', 376.835)]:
        model = models[(smiles, 0)]
        assert [(x['energy_cm_inverse'], x['degeneracy']) for x in model['levels']] == [(0., 2), (gap, 2)]
        assert model['spin_orbit_correction']['value_hartree'] == pytest.approx(-gap / 2 * invcm / Hartree)
    assert {smiles for (smiles, charge), model in models.items()
            if 'spin_orbit_correction' in model} == {'[S]', '[SH]', 'C[S]'}
    methyl = SimpleNamespace(smiles='[CH3]', charge=0, mult=2)
    sulfur = SimpleNamespace(smiles='[S]', charge=0, mult=3)
    methylthiyl = SimpleNamespace(smiles='C[S]', charge=0, mult=2)
    disulfur = SimpleNamespace(smiles='[S][S]', charge=0, mult=3)
    assert state_average_roots([methyl, sulfur], [2, 4], parameters) == {2: 3, 4: 3}
    if case == 'ch3ss':
        assert state_average_roots([methyl, disulfur], [2, 4], parameters) == {2: 1, 4: 1}
        assert state_average_roots([methylthiyl, sulfur], [2, 4], parameters) == {2: 6, 4: 6}


def test_molecular_centroid_is_general_and_can_be_explicitly_overridden():
    from ase.units import invcm, Hartree
    # An illustrative OH input exercises the same generic Pi-term interface;
    # this synthetic gap is not added to the sulfur examples as spectroscopic data.
    item = dict(smiles='[OH]', charge=0, ground_multiplicity=2,
                term_orbital_degeneracy=2, spin_free_reference='term_centroid',
                levels=[dict(energy_cm_inverse=x, multiplicity=2, degeneracy=2)
                        for x in (0., 140.)])
    par = {'electronic_state_models': [item]}
    assert electronic_state_models(par)[('[OH]', 0)]['spin_orbit_correction']['value_hartree'] == pytest.approx(-70 * invcm / Hartree)
    item['spin_orbit_correction_hartree'] = -.001
    with pytest.raises(ValueError, match='not both'):
        electronic_state_models(par)
    item.pop('spin_free_reference')
    assert electronic_state_models(par)[('[OH]', 0)]['spin_orbit_correction']['value_hartree'] == -.001


@pytest.mark.parametrize('configured', [False, True])
def test_cached_reference_receives_current_mess_levels_without_changing_qc(tmp_path, configured):
    import json
    from copy import deepcopy
    from pathlib import Path
    from kinbot import kb_path
    writer = MESS.__new__(MESS)
    writer.par = {'TemperatureList': [300., 2000.],
                  'freq_uq_ref': 100., 'freq_uq_max_exp': 4.}
    if configured:
        path = Path(__file__).parents[1] / 'examples/anl/ch3s_anl1_f12/ch3s.json'
        writer.par.update(electronic_state_models=json.loads(path.read_text())['electronic_state_models'])
    writer.pftpl = Path(kb_path, 'tpl/mess_pf.tpl').read_text()
    model = {'smiles': '[SH]', 'name': 'SH', 'charge': 0, 'multiplicity': 2,
             'symbols': ['S', 'H'], 'positions': [[0., 0., 0.], [1.34, 0., 0.]],
             'symmetry_factor': 1., 'frequencies_cm_inverse': [2600.],
             'electronic_levels': [{'energy_cm_inverse': 0., 'degeneracy': 2}],
             'source': {'record_sha256': 'original QC provenance', 'energy_hartree': -398.}}
    original = deepcopy(model)
    text = writer.write_partition_function_model(model, tmp_path/'SH.inp').read_text()
    assert model == original
    assert f'ElectronicLevels[1/cm] {2 if configured else 1}' in ' '.join(text.split())
    assert ('376.835 2' in text) == configured
    assert '2600' in text
    if configured:
        model['multiplicity'] = 4
        with pytest.raises(ValueError, match='multiplicity disagrees'):
            writer.write_partition_function_model(model, tmp_path/'bad.inp')
