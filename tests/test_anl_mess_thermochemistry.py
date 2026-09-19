"""Scientific contracts for MESSPF and the 0 K/298 K thermochemistry split."""

from io import StringIO

import pytest

from kinbot.anl.thermochemistry import (
    ELEMENT_REFERENCE_HEAT_CONTENT_298_KJ_MOL_ATOM,
    PartitionFunctionPoint,
    formation_enthalpy_298_kj_mol,
    fit_nasa7,
    nasa7_enthalpy_kj_mol,
    nasa7_entropy_j_mol_k,
    parse_messpf_output,
    point_at_temperature,
)


def test_ram_c2f6_mess_energy_is_cbh_formation_enthalpy_difference():
    # Ram et al. SI: C2F6 well relative to the 2 CF3 channel.
    well = -318.333
    fragments = 2 * -111.164
    assert well - fragments == pytest.approx(-96.005)


def test_messpf_derivative_gives_ideal_gas_heat_content():
    output = StringIO('''
Natural log of the partition function, its derivatives, entropy, and thermal capacity:
T, K       methane      methane      methane      methane      methane
                   Z_0          Z_1          Z_2 S, cal/mol/K C, cal/mol/K
200.0       10.0         0.009        0.0          40.0         8.0
298.15      11.0         0.010        0.0          45.0         9.0
''')
    point = point_at_temperature(parse_messpf_output(output))
    expected = 0.00831446261815324 * (298.15 ** 2 * .010 + 298.15)
    assert point.heat_content_0k_kj_mol == pytest.approx(expected)


def test_hf298_uses_element_standard_state_heat_contents():
    elemental = ELEMENT_REFERENCE_HEAT_CONTENT_298_KJ_MOL_ATOM
    value = formation_enthalpy_298_kj_mol(
        -66.5, 10., ['C', 'H', 'H', 'H', 'H'])
    expected = -66.5 + 10. - elemental['C'] - 4 * elemental['H']
    assert value == pytest.approx(expected)


def test_messpf_requires_explicit_298_15_row():
    output = StringIO('''
Z_0 Z_1 Z_2 S, cal/mol/K C, cal/mol/K
298.2 1 0.01 0 10 5
''')
    with pytest.raises(ValueError, match='298.15'):
        point_at_temperature(parse_messpf_output(output))


def test_nasa7_fit_uses_hf298_and_messpf_entropy_anchors():
    rows = []
    for temperature in (200., 250., 298.15, 400., 600., 800., 1000.,
                        1200., 1600., 2200., 3000.):
        cp_cal = (3.5 + 1.e-4 * temperature) * 8.31446261815324 / 4.184
        rows.append(PartitionFunctionPoint(
            temperature, 1., .01, 0., 50. + .001 * temperature, cp_cal))
    fit = fit_nasa7(rows, -123.456)
    assert nasa7_enthalpy_kj_mol(
        fit['low_coefficients'], 298.15) == pytest.approx(-123.456)
    assert nasa7_entropy_j_mol_k(
        fit['low_coefficients'], 298.15) == pytest.approx(
            rows[2].entropy_cal_mol_k * 4.184)
    assert nasa7_enthalpy_kj_mol(
        fit['low_coefficients'], 1000.) == pytest.approx(
            nasa7_enthalpy_kj_mol(fit['high_coefficients'], 1000.))
    assert nasa7_entropy_j_mol_k(
        fit['low_coefficients'], 1000.) == pytest.approx(
            nasa7_entropy_j_mol_k(fit['high_coefficients'], 1000.))
    assert fit['low_cp_rmse_j_mol_k'] < 1.e-8
    assert fit['high_cp_rmse_j_mol_k'] < 1.e-8
