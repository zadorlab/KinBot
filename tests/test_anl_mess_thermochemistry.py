"""Scientific contracts for MESSPF and the 0 K/298 K thermochemistry split."""

from io import StringIO

import pytest

from kinbot.anl.atct import ATcTRecord
from kinbot.anl.cbh import FormationEnthalpy
from kinbot.anl.thermochemistry import (
    PartitionFunctionPoint,
    cbh_formation_enthalpy_at_temperature,
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


def _atct(smiles, h0, h298, ident):
    return ATcTRecord(
        ident, smiles, f'{smiles} (g)', 'g', smiles, h0, h298, .01,
        False, '1.220', 'a' * 64)


def test_hf298_uses_the_same_cbh0_reaction_and_atct_references():
    methane = _atct('C', -66.5, -74.5, 'methane*0')
    hydrogen = _atct('[H][H]', 0., 0., 'hydrogen*0')
    formation = FormationEnthalpy(
        target_smiles='CC', rung=0, method='ANL0-F12',
        reaction_energy_0k_kj_mol=-64., formation_0k_kj_mol=-69.,
        references={'C': methane, '[H][H]': hydrogen},
        atct_version='1.220', atct_source_sha256='a' * 64,
        energy_sources={'CC': 'ethane', 'C': 'methane', '[H][H]': 'hydrogen'},
        stoichiometry={'C': 2, 'CC': -1, '[H][H]': -1})
    heat = {'CC': 12., 'C': 10., '[H][H]': 8.}
    value, trace = cbh_formation_enthalpy_at_temperature(formation, heat)
    thermal_reaction = 2 * 10. - 12. - 8.
    expected = 2 * -74.5 - (-64. + thermal_reaction)
    assert value == pytest.approx(expected)
    assert trace['thermal_reaction_correction_kj_mol'] == pytest.approx(
        thermal_reaction)
    assert trace['references']['C']['kind'] == 'atct'


def test_higher_cbh_rung_recurses_through_computed_lower_rung_values():
    methane = _atct('C', -66.5, -74.5, 'methane*0')
    hydrogen = _atct('[H][H]', 0., 0., 'hydrogen*0')
    ethane = FormationEnthalpy(
        target_smiles='CC', rung=0, method='ANL0-F12',
        reaction_energy_0k_kj_mol=-64., formation_0k_kj_mol=-69.,
        references={'C': methane, '[H][H]': hydrogen},
        atct_version='1.220', atct_source_sha256='a' * 64,
        energy_sources={'CC': 'ethane', 'C': 'methane', '[H][H]': 'hydrogen'},
        stoichiometry={'C': 2, 'CC': -1, '[H][H]': -1})
    propane = FormationEnthalpy(
        target_smiles='CCC', rung=1, method='F12',
        reaction_energy_0k_kj_mol=2., formation_0k_kj_mol=-140.,
        references={'CC': ethane, 'C': FormationEnthalpy(
            target_smiles='C', rung=0, method='ANL0-F12',
            reaction_energy_0k_kj_mol=0., formation_0k_kj_mol=-66.5,
            references={'C': methane}, atct_version='1.220',
            atct_source_sha256='a' * 64, energy_sources={},
            stoichiometry={})},
        atct_version='1.220', atct_source_sha256='a' * 64,
        energy_sources={'CCC': 'propane', 'CC': 'ethane', 'C': 'methane'},
        stoichiometry={'CCC': -1, 'C': -1, 'CC': 2})
    heat = {'CCC': 15., 'CC': 12., 'C': 10., '[H][H]': 8.}
    value, trace = cbh_formation_enthalpy_at_temperature(propane, heat)
    ethane_298 = 2 * -74.5 - (-64. + 2 * 10. - 12. - 8.)
    expected = 2 * ethane_298 - -74.5 - (2. + 2 * 12. - 10. - 15.)
    assert value == pytest.approx(expected)
    assert trace['references']['CC']['kind'] == 'cbh_reaction'
    assert trace['references']['C']['kind'] == 'atct_terminal'

    invalid = FormationEnthalpy(
        target_smiles='CCC', rung=1, method='F12',
        reaction_energy_0k_kj_mol=2., formation_0k_kj_mol=-140.,
        references={'CC': _atct('CC', -69., -84., 'ethane*0'),
                    'C': methane}, atct_version='1.220',
        atct_source_sha256='a' * 64, energy_sources={},
        stoichiometry={'CCC': -1, 'C': -1, 'CC': 2})
    with pytest.raises(ValueError, match='Only CBH-0'):
        cbh_formation_enthalpy_at_temperature(invalid, heat)


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
