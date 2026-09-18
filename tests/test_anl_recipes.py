"""Complete-expression arithmetic and provenance gates use synthetic energies."""

from dataclasses import replace

import pytest

from kinbot.anl.model import ComponentResult, IncompleteRecipeError
from kinbot.anl.recipes import recipe


L2 = '1' * 64
L3 = '2' * 64
SOURCE = 'a' * 64


def components_for(equation, *, multiplicity=1):
    values = {
        'reference_cbs': -100.0,
        'harmonic_zpe': 0.1,
        'vpt2_correction': -0.01,
        'hoe_high': -99.8,
        'hoe_low': -99.7,
        'hoe_tz_high': -99.8,
        'hoe_tz_low': -99.7,
        'hoe_dz_high': -99.9,
        'hoe_dz_low': -99.85,
        'core_valence_cbs': -0.02,
        'scalar_relativistic': -0.005,
        'dboc': 0.002,
        'spin_orbit': 0.0,
    }
    result = {}
    for requirement in equation.requirements:
        backend = requirement.backends[0]
        if requirement.key in ('hoe_high', 'hoe_tz_high') and multiplicity > 1:
            backend = 'mrcc'
        if requirement.key == 'hoe_dz_low' and multiplicity > 1:
            backend = 'mrcc'
        result[requirement.key] = ComponentResult(
            key=requirement.key, value_hartree=values[requirement.key],
            quantity=requirement.quantity, method=requirement.method,
            basis=requirement.basis, backend=backend, state_id='state-A',
            charge=0, multiplicity=multiplicity,
            geometry_sha256=(None if requirement.geometry_role == 'state'
                             else {'l2': L2, 'l3': L3}[requirement.geometry_role]),
            source_sha256=SOURCE, source='synthetic calculation record',
            settings=dict(requirement.settings),
        )
    return result


def evaluate(equation, components, *, multiplicity=1):
    return equation.evaluate(components, state_id='state-A', charge=0,
                             multiplicity=multiplicity,
                             geometry_hashes={'l2': L2, 'l3': L3})


@pytest.mark.parametrize('name', ['ANL0', 'ANL0-F12'])
def test_complete_anl0_expressions_include_zpe_once(name):
    equation = recipe(name)
    components = components_for(equation)
    result = evaluate(equation, components)
    assert result.electronic_hartree == pytest.approx(-100.123)
    assert result.zero_point_hartree == pytest.approx(0.09)
    assert result.zero_k_hartree == pytest.approx(-100.033)
    assert equation.name == (
        'profiled:ANL0-F12:scaled-triples' if name == 'ANL0-F12'
        else name)
    if name == 'ANL0-F12':
        assert components['reference_cbs'].method == 'CCSD(T)-F12b'
        assert components['reference_cbs'].settings['scale_trip'] == 1
    else:
        assert components['reference_cbs'].method == 'CCSD(T)'


def test_anl1_cross_program_higher_order_pairs_and_profiled_label():
    equation = recipe('ANL1', vpt2_method='B2PLYP-D3BJ')
    components = components_for(equation)
    result = evaluate(equation, components)
    assert result.electronic_hartree == pytest.approx(-100.173)
    assert result.zero_k_hartree == pytest.approx(-100.083)
    assert result.recipe == 'profiled:ANL1:B2PLYP-D3BJ'
    assert components['hoe_tz_high'].backend == 'cfour'
    assert components['hoe_tz_low'].backend == 'molpro'
    assert components['hoe_dz_high'].backend == 'mrcc'
    assert components['vpt2_correction'].settings['dispersion'] == 'GD3BJ'


def test_optional_vpt2_cbs_remains_one_zero_point_correction():
    equation = recipe('ANL1', vpt2_method='B2PLYP-D3BJ', vpt2_cbs=True)
    components = components_for(equation)
    result = evaluate(equation, components)
    assert result.recipe == 'profiled:ANL1:B2PLYP-D3BJ:VPT2-CBS'
    assert result.zero_point_hartree == pytest.approx(0.09)
    assert components['vpt2_correction'].basis == 'CBS(cc-pVTZ,cc-pVQZ)'
    assert components['vpt2_correction'].backend == 'composite'


def test_incomplete_or_unreviewed_components_cannot_be_called_anl():
    equation = recipe('ANL0-F12')
    components = components_for(equation)
    del components['hoe_high']
    with pytest.raises(IncompleteRecipeError, match='hoe_high'):
        evaluate(equation, components)
    components = components_for(equation)
    components['vpt2_correction'] = replace(components['vpt2_correction'],
                                            review_required=True)
    with pytest.raises(IncompleteRecipeError, match='quality review'):
        evaluate(equation, components)


def test_geometry_method_and_source_identity_are_enforced():
    equation = recipe('ANL0-F12')
    components = components_for(equation)
    components['vpt2_correction'] = replace(components['vpt2_correction'],
                                            geometry_sha256=L3)
    with pytest.raises(ValueError, match='geometry differs'):
        evaluate(equation, components)
    components = components_for(equation)
    components['reference_cbs'] = replace(components['reference_cbs'],
                                          method='CCSD(T)-F12a')
    with pytest.raises(ValueError, match='method'):
        evaluate(equation, components)
    components = components_for(equation)
    components['reference_cbs'] = replace(components['reference_cbs'],
                                          settings={'scale_trip': 0})
    with pytest.raises(ValueError, match='settings'):
        evaluate(equation, components)
    components = components_for(equation)
    components['dboc'] = replace(components['dboc'], source_sha256='')
    with pytest.raises(ValueError, match='provenance'):
        evaluate(equation, components)


def test_open_shell_higher_order_uses_mrcc_and_spin_orbit_is_explicit():
    equation = recipe('ANL0')
    components = components_for(equation, multiplicity=2)
    assert evaluate(equation, components, multiplicity=2).zero_k_hartree == \
        pytest.approx(-100.033)
    components['hoe_high'] = replace(components['hoe_high'], backend='cfour')
    with pytest.raises(ValueError, match='backend'):
        evaluate(equation, components, multiplicity=2)
    components = components_for(equation, multiplicity=2)
    del components['spin_orbit']
    with pytest.raises(IncompleteRecipeError, match='spin_orbit'):
        evaluate(equation, components, multiplicity=2)

    extended = recipe('ANL1')
    extended_components = components_for(extended, multiplicity=2)
    assert evaluate(extended, extended_components, multiplicity=2).recipe == 'ANL1'
    extended_components['hoe_dz_low'] = replace(
        extended_components['hoe_dz_low'], backend='cfour')
    with pytest.raises(ValueError, match='backend'):
        evaluate(extended, extended_components, multiplicity=2)
