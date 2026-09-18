"""CBH generation, ATcT ingestion, and 0 K thermochemistry workflow tests."""

from io import BytesIO
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot.anl.atct import load_atct, parse_atct_html
from kinbot.anl.cbh import (HARTREE_TO_KJ_MOL, ZeroKEnergy,
                            generate_cbh_reaction, generate_for_stationary_point,
                            select_cbh_ladder, solve_formation_enthalpy)


def _row(number, atct_id, name, formula, smiles, h0, uncertainty='0.1'):
    return (f'<tr id="s1n1c0 i{number} CASx">'
            f'<td><span class="Name"><a href="species/?species_number={number}">'
            f'{name}</a></span></td>'
            f'<td class="bkgFormula> <span class="Formula"><button>'
            f'{formula}</button></span></td>'
            f'<td><img alt="{smiles}" /></td>'
            f'<td><span class="DHf0">{h0}</span></td>'
            f'<td><span class="DHf298">0</span></td>'
            f'<td><span class="Uncert">{uncertainty}</span></td>'
            f'<td><span class="Units">'
            f'{"" if uncertainty == "exact" else "kJ/mol"}</span></td>'
            f'<td><span class="ATcTID">{atct_id}</span></td></tr>')


def _table():
    raw = ('ATcT enthalpies of formation based on version 1.222 of the '
           'Thermochemical Network' +
           _row(1, '1333-74-0*0', 'Dihydrogen', 'H2 (g)', '[H][H]',
                '0', 'exact') +
           _row(2, '74-82-8*0', 'Methane', 'CH4 (g)', 'C', '-66.544') +
           _row(3, '74-84-0*0', 'Ethane', 'C2H6 (g)', 'CC', '-68.400') +
           _row(4, '7732-18-5*0', 'Water', 'H2O (g)', 'O', '-238.920') +
           _row(5, '7732-18-5*1', 'Water', 'H2O (g, ortho)', 'O',
                '-238.920')).encode()
    return raw, parse_atct_html(raw, '1.222')


def test_atct_pinned_reader_and_state_matching(tmp_path, monkeypatch):
    raw, table = _table()
    assert table.by_id('74-82-8*0').formation_0k_kj_mol == -66.544
    assert table.gas_by_smiles('O', formula={'H': 2, 'O': 1}).atct_id == '7732-18-5*0'
    assert table.by_id('7732-18-5*1').preferred_formula.endswith('(g, ortho)')
    with pytest.raises(LookupError):
        table.gas_by_smiles('N', formula={'N': 1, 'H': 3})
    with pytest.raises(ValueError, match='header'):
        parse_atct_html(raw, '1.220')
    calls = []

    def fetch(url, timeout):
        calls.append(url)
        return BytesIO(raw)

    monkeypatch.setattr('kinbot.anl.atct.urlopen', fetch)
    first = load_atct('1.222', tmp_path)
    second = load_atct('1.222', tmp_path)
    assert first.source_sha256 == second.source_sha256
    assert len(calls) == 1
    assert (Path(tmp_path) / 'atct_1.222.html').read_bytes() == raw


def test_original_cbh_rungs_from_smiles_connectivity():
    assert generate_cbh_reaction('C', 0) is None
    assert generate_cbh_reaction('CC', 0).stoichiometry == {
        'CC': -1, '[H][H]': -1, 'C': 2}
    assert generate_cbh_reaction('CCC', 1).stoichiometry == {
        'CCC': -1, 'C': -1, 'CC': 2}
    assert generate_cbh_reaction('CCCC', 2).stoichiometry == {
        'CCCC': -1, 'CC': -1, 'CCC': 2}
    assert generate_cbh_reaction('CCCCC', 3).stoichiometry == {
        'CCCCC': -1, 'CCC': -1, 'CCCC': 2}
    assert generate_cbh_reaction('CCO', 1).stoichiometry == {
        'CCO': -1, 'C': -1, 'CC': 1, 'CO': 1}
    assert generate_cbh_reaction('C1CCCCC1', 3) is None


def test_cbh_rejects_states_without_validated_fragment_rules():
    with pytest.raises(NotImplementedError, match='closed-shell'):
        generate_cbh_reaction('[CH3]', 0, multiplicity=2)
    with pytest.raises(NotImplementedError, match='nonaromatic'):
        generate_cbh_reaction('c1ccccc1', 0)
    with pytest.raises(NotImplementedError, match='unusual-valence'):
        generate_cbh_reaction('[CH3]', 0)
    species = SimpleNamespace(smiles='CCC', charge=0, mult=1, wellorts=0,
                              atom=['C'] * 3 + ['H'] * 8)
    assert generate_for_stationary_point(species)[1].rung == 1
    species.atom.pop()
    with pytest.raises(ValueError, match='atoms disagree'):
        generate_for_stationary_point(species)
    species.wellorts = 1
    with pytest.raises(ValueError, match='transition state'):
        generate_for_stationary_point(species)


def test_zero_k_formation_uses_accepted_energies_and_atct_references():
    _, table = _table()
    reaction = generate_cbh_reaction('CC', 0)
    energies = {
        'CC': ZeroKEnergy('CC', -79.0, 'ANL0-F12', 'ethane-output'),
        'C': ZeroKEnergy('C', -40.0, 'ANL0-F12', 'methane-output'),
        '[H][H]': ZeroKEnergy('[H][H]', -1.0, 'ANL0-F12', 'hydrogen-output'),
    }
    result = solve_formation_enthalpy(reaction, energies, table)
    assert result.reaction_energy_0k_kj_mol == pytest.approx(0)
    assert result.formation_0k_kj_mol == pytest.approx(2 * -66.544)
    assert result.atct_version == '1.222'
    energies['CC'] = ZeroKEnergy('CC', -79.001, 'ANL0-F12', 'new-output')
    result = solve_formation_enthalpy(reaction, energies, table)
    assert result.reaction_energy_0k_kj_mol == pytest.approx(
        0.001 * HARTREE_TO_KJ_MOL)
    assert result.formation_0k_kj_mol == pytest.approx(
        2 * -66.544 - 0.001 * HARTREE_TO_KJ_MOL)
    energies['C'] = ZeroKEnergy('C', -40, 'ANL1', 'different-method')
    with pytest.raises(ValueError, match='same'):
        solve_formation_enthalpy(reaction, energies, table)


def test_ladder_uses_highest_complete_common_tier_for_available_rung():
    _, table = _table()
    species = SimpleNamespace(smiles='CCC', charge=0, mult=1, wellorts=0,
                              atom=['C'] * 3 + ['H'] * 8)
    energies = {
        smiles: {tier: ZeroKEnergy(smiles, value, tier, smiles + tier)
                 for tier in ('ANL0', 'L2')}
        for smiles, value in [('CCC', -118.), ('CC', -79.), ('C', -40.)]
    }
    selected = select_cbh_ladder(species, energies, table)
    assert selected.reaction.rung == 1
    assert selected.formation.method == 'ANL0'
    del energies['CC']['ANL0']
    selected = select_cbh_ladder(species, energies, table)
    assert selected.formation.method == 'L2'
    assert any('ANL0' in reason for reason in selected.skipped)
