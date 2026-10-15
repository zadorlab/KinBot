"""CBH generation, ATcT ingestion, and 0 K thermochemistry workflow tests."""

from dataclasses import asdict
from io import BytesIO
from pathlib import Path
import json
import os
import sys
from types import ModuleType, SimpleNamespace

import pytest

from kinbot.anl.atct import (_fetch_reference_records, load_atct,
                             load_atct_references,
                             parse_atct_api, parse_atct_html)
from kinbot.anl.cbh import (HARTREE_TO_KJ_MOL, ZeroKEnergy,
                            formation_from_dict,
                            generate_cbh_reaction, generate_for_stationary_point,
                            select_cbh_ladder, solve_composite_records,
                            solve_formation_enthalpy,
                            zero_k_from_composite_record)


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


def _api_payload():
    return json.dumps([
        {'ATcT_TN_Version': '1.220', 'ATcT_ID': '1333-74-0*0',
         'Name': 'Dihydrogen', 'Formula': 'H2  (g)',
         'Delta_Hf_0K': '0', 'Delta_Hf_298K': '0',
         'Delta_Hf298K_uncertainty': 'exact', 'unit': None,
         'SMILES': '[H][H]'},
        {'ATcT_TN_Version': '1.220', 'ATcT_ID': '74-82-8*0',
         'Name': 'Methane', 'Formula': 'CH4  (g)',
         '∆fH_0K': '-66.543', '∆fH_298K': '-74.513',
         '∆fH_298K_uncertainty': '0.043', 'units': 'kJ/mol',
         'SMILES': 'C'},
    ], ensure_ascii=False).encode()


def test_atct_pinned_reader_and_state_matching(tmp_path, monkeypatch):
    raw, table = _table()
    assert table.by_id('74-82-8*0').formation_0k_kj_mol == -66.544
    assert table.gas_by_smiles('O', formula={'H': 2, 'O': 1}).atct_id == '7732-18-5*0'
    assert table.by_id('7732-18-5*1').preferred_formula.endswith('(g, ortho)')
    with pytest.raises(LookupError):
        table.gas_by_smiles('N', formula={'N': 1, 'H': 3})
    with pytest.raises(ValueError, match='header'):
        parse_atct_html(raw, '1.220')
    api_raw = _api_payload()
    api_table = parse_atct_api(api_raw, '1.220')
    assert api_table.by_id('74-82-8*0').formation_0k_kj_mol == -66.543
    assert api_table.by_id('1333-74-0*0').exact
    with pytest.raises(ValueError, match='expected 1.222'):
        parse_atct_api(api_raw, '1.222')

    calls = []

    def fetch(url, timeout):
        calls.append(url)
        return BytesIO(api_raw)

    monkeypatch.setattr('kinbot.anl.atct.urlopen', fetch)
    first = load_atct('1.220', tmp_path)
    second = load_atct('1.220', tmp_path)
    assert first.source_sha256 == second.source_sha256
    assert len(calls) == 1
    assert (Path(tmp_path) / 'atct_1.220.json').read_bytes() == api_raw


def test_atct_client_reference_query_is_cached(tmp_path, monkeypatch):
    calls = []

    def fetch(smiles, reference_ids):
        calls.append((tuple(smiles), dict(reference_ids)))
        return json.loads(_api_payload())

    monkeypatch.setattr('kinbot.anl.atct._fetch_reference_records', fetch)
    first = load_atct_references(
        ['C', '[H][H]'], '1.220', tmp_path,
        reference_ids={'[H][H]': '1333-74-0*0'})
    second = load_atct_references(
        ['[H][H]', 'C'], '1.220', tmp_path,
        reference_ids={'[H][H]': '1333-74-0*0'})
    assert first.source_sha256 == second.source_sha256
    assert first.gas_by_smiles('C', formula={'C': 1, 'H': 4}).atct_id == '74-82-8*0'
    assert len(calls) == 1


def test_atct_client_ignores_invalid_generic_proxy(monkeypatch):
    observed = {}

    class Species:
        @staticmethod
        def to_dict():
            return {'ATcT_ID': '74-82-8*0',
                    '∆fH_298K_uncertainty': '0.043'}

    def by_smiles(*args, **kwargs):
        observed['ALL_PROXY'] = os.environ.get('ALL_PROXY')
        observed['all_proxy'] = os.environ.get('all_proxy')
        observed['HTTPS_PROXY'] = os.environ.get('HTTPS_PROXY')
        return SimpleNamespace(total=1, items=[Species()])

    api = ModuleType('atct.api')
    api.get_species_by_atctid = lambda *args, **kwargs: Species()
    api.get_species_by_smiles = by_smiles
    package = ModuleType('atct')
    package.api = api
    monkeypatch.setitem(sys.modules, 'atct', package)
    monkeypatch.setitem(sys.modules, 'atct.api', api)
    monkeypatch.setenv('ALL_PROXY', 'socks://proxy.example:80')
    monkeypatch.setenv('all_proxy', 'socks://proxy.example:80')
    monkeypatch.setenv('HTTPS_PROXY', 'http://proxy.example:80')

    records = _fetch_reference_records(['C'], {})

    assert records[0]['ATcT_ID'] == '74-82-8*0'
    assert observed == {'ALL_PROXY': None, 'all_proxy': None,
                        'HTTPS_PROXY': 'http://proxy.example:80'}
    assert os.environ['ALL_PROXY'] == 'socks://proxy.example:80'
    assert os.environ['all_proxy'] == 'socks://proxy.example:80'


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
    assert generate_cbh_reaction('C1CCCCC1', 3).stoichiometry == {
        'C1CCCCC1': -1, 'CCC': -6, 'CCCC': 6}


def test_cbh_preserves_radicals_aromatic_bonds_ions_and_unusual_valence():
    # A one-fragment radical uses the literature CBH-0 parent-hydride/H
    # terminal. A larger radical retains its center in each capped fragment.
    methyl = generate_cbh_reaction('[CH3]', 0, multiplicity=2)
    assert methyl.stoichiometry == {
        '[CH3]': -1, '[H][H]': -1, 'C': 1, '[H]': 1}
    assert methyl.states == {'[CH3]': (0, 2), '[H]': (0, 2)}
    ethyl = generate_cbh_reaction('C[CH2]', 0, multiplicity=2)
    assert ethyl.stoichiometry == {
        'C': 1, '[CH2]C': -1, '[CH3]': 1, '[H][H]': -1}
    assert ethyl.states == {
        '[CH2]C': (0, 2), '[CH3]': (0, 2)}

    benzene = generate_cbh_reaction('c1ccccc1', 1)
    assert benzene.stoichiometry == {
        'C': -6, 'C=C': 3, 'CC': 3, 'c1ccccc1': -1}
    phenyl = generate_cbh_reaction('[c]1ccccc1', 1, multiplicity=2)
    assert phenyl.states['[c]1ccccc1'] == (0, 2)
    assert all(phenyl.states[fragment] == (0, 2)
               for fragment in ('[CH]=C', '[CH2]C', '[CH3]'))

    formate = generate_cbh_reaction('[O-]C=O', 1, charge=-1)
    assert formate.stoichiometry == {
        'C': -1, 'C=O': 1, 'C[O-]': 1, 'O=C[O-]': -1}
    assert formate.states == {
        'O=C[O-]': (-1, 1), 'C[O-]': (-1, 1)}
    zwitterion = generate_cbh_reaction('[NH3+]CC(=O)[O-]', 1)
    assert zwitterion.states['C[NH3+]'] == (1, 1)
    assert zwitterion.states['C[O-]'] == (-1, 1)

    sulfone = generate_cbh_reaction('CS(=O)(=O)C', 1)
    assert sulfone.stoichiometry == {
        'CS': 2, 'CS(C)(=O)=O': -1, 'O=S': 2, 'S': -3,
        '[H][H]': 2}
    nitro = generate_cbh_reaction('C[N+](=O)[O-]', 1)
    assert nitro.states['[NH4+]'] == (1, 1)

    with pytest.raises(ValueError, match='formal charge'):
        generate_cbh_reaction('[NH4+]', 0, charge=0)
    with pytest.raises(ValueError, match='explicitly identify radical'):
        generate_cbh_reaction('O=O', 0, multiplicity=3)
    species = SimpleNamespace(smiles='CCC', charge=0, mult=1, wellorts=0,
                              atom=['C'] * 3 + ['H'] * 8)
    assert generate_for_stationary_point(species)[1].rung == 1
    species.atom.pop()
    with pytest.raises(ValueError, match='atoms disagree'):
        generate_for_stationary_point(species)
    species.wellorts = 1
    with pytest.raises(ValueError, match='transition state'):
        generate_for_stationary_point(species)


@pytest.mark.parametrize('left,right,charge,multiplicity', [
    ('c1ccc2ccccc2c1', 'c1cccc2ccccc12', 0, 1),
    ('n1ccccc1', 'c1ccncc1', 0, 1),
    ('C[N+](=O)[O-]', '[O-][N+](=O)C', 0, 1),
    ('[NH3+]CC(=O)[O-]', '[O-]C(=O)C[NH3+]', 0, 1),
    ('[CH2]c1ccccc1', 'c1ccccc1[CH2]', 0, 2),
])
def test_cbh_is_invariant_to_equivalent_smiles_traversal(
        left, right, charge, multiplicity):
    for rung in range(4):
        assert generate_cbh_reaction(
            left, rung, charge=charge, multiplicity=multiplicity) == \
            generate_cbh_reaction(
                right, rung, charge=charge, multiplicity=multiplicity)


@pytest.mark.parametrize('smiles,charge,multiplicity,rung', [
    ('c1occc1', 0, 1, 2),
    ('c1sccc1', 0, 1, 2),
    ('O=S(=O)([O-])[O-]', -2, 1, 1),
    ('O=P([O-])([O-])[O-]', -3, 1, 1),
    ('FS(F)(F)(F)(F)F', 0, 1, 1),
    ('ClP(Cl)(Cl)(Cl)Cl', 0, 1, 1),
    ('C=C=C', 0, 1, 1),
    ('[O][O]', 0, 3, 0),
    ('[N]=O', 0, 2, 0),
])
def test_cbh_balances_heteroaromatic_ionic_hypervalent_and_radical_corpus(
        smiles, charge, multiplicity, rung):
    reaction = generate_cbh_reaction(
        smiles, rung, charge=charge, multiplicity=multiplicity)
    assert reaction is not None
    assert reaction.stoichiometry[reaction.target_smiles] == -1


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


def test_zero_k_formation_propagates_cbs_half_width_and_atct_uncertainty():
    _, table = _table()
    reaction = generate_cbh_reaction('CC', 0)
    series = 'reference_cbs:f12'
    energies = {
        'CC': ZeroKEnergy(
            'CC', -79., 'ANL0-F12', 'ethane-output',
            basis_uncertainty_series={series: {
                'lower_hartree': -79., 'upper_hartree': -79.001}}),
        'C': ZeroKEnergy(
            'C', -40., 'ANL0-F12', 'methane-output',
            basis_uncertainty_series={series: {
                'lower_hartree': -40., 'upper_hartree': -40.0002}}),
        '[H][H]': ZeroKEnergy(
            '[H][H]', -1., 'ANL0-F12', 'hydrogen-output',
            basis_uncertainty_series={series: {
                'lower_hartree': -1., 'upper_hartree': -1.}}),
    }

    result = solve_formation_enthalpy(reaction, energies, table)

    cbs = .5 * .0006 * HARTREE_TO_KJ_MOL
    atct = 2 * .1
    assert result.uncertainty_components_kj_mol == pytest.approx({
        f'cbs:{series}': cbs,
        'atct:74-82-8*0': atct,
        'atct:1333-74-0*0': 0.,
    })
    assert result.uncertainty_kj_mol == pytest.approx(
        (cbs * cbs + atct * atct) ** .5)
    assert formation_from_dict(asdict(result)) == result


def test_zero_k_formation_rejects_incomplete_cbs_uncertainty_series():
    _, table = _table()
    reaction = generate_cbh_reaction('CC', 0)
    energies = {
        smiles: ZeroKEnergy(
            smiles, value, 'ANL0-F12', smiles,
            basis_uncertainty_series=(
                {'reference_cbs': {
                    'lower_hartree': value, 'upper_hartree': value}}
                if smiles != '[H][H]' else {}))
        for smiles, value in {'CC': -79., 'C': -40., '[H][H]': -1.}.items()
    }
    with pytest.raises(ValueError, match='missing CBS uncertainty series'):
        solve_formation_enthalpy(reaction, energies, table)


def _composite_record(path, electronic, zpe=0., *,
                      recipe='profiled:ANL0-F12:scaled-triples:B2PLYP-D3BJ'):
    components = {
        'electronic': {
            'key': 'electronic', 'value_hartree': electronic,
            'review_required': False, 'source_sha256': 'a' * 64,
        },
        'zpe': {
            'key': 'zpe', 'value_hartree': zpe,
            'review_required': False, 'source_sha256': 'b' * 64,
        },
    }
    path.write_text(json.dumps({
        'schema': 1, 'status': 'complete',
        'recipe': recipe,
        'recipe_version': 1, 'charge': 0, 'multiplicity': 1,
        'components': components,
        'electronic_terms': [{'component': 'electronic', 'coefficient': 1}],
        'zero_point_terms': [{'component': 'zpe', 'coefficient': 1}],
        'electronic_hartree': electronic, 'zero_point_hartree': zpe,
        'zero_k_hartree': electronic + zpe,
    }))


def test_composite_record_handoff_solves_ethane_cbh0(tmp_path, monkeypatch):
    files = {}
    for smiles, name, value in (
            ('CC', 'ethane.json', -79.), ('C', 'methane.json', -40.),
            ('[H][H]', 'hydrogen.json', -1.)):
        path = tmp_path / name
        _composite_record(path, value)
        files[smiles] = path
    energy = zero_k_from_composite_record('CC', files['CC'])
    assert energy.hartree == -79.
    assert energy.method.startswith('profiled:ANL0-F12:')
    assert '#sha256=' in energy.source
    _, table = _table()
    monkeypatch.setattr('kinbot.anl.cbh.load_atct_references',
                        lambda *args, **kwargs: table)
    reaction, formation = solve_composite_records(
        'CC', 0, files, atct_version='1.222', atct_cache=tmp_path)
    assert reaction.stoichiometry == {'C': 2, 'CC': -1, '[H][H]': -1}
    assert formation.reaction_energy_0k_kj_mol == pytest.approx(0.)
    assert formation.formation_0k_kj_mol == pytest.approx(2 * -66.544)
    record = json.loads(files['CC'].read_text())
    record['zero_k_hartree'] += .1
    files['CC'].write_text(json.dumps(record))
    with pytest.raises(ValueError, match='disagrees with components'):
        zero_k_from_composite_record('CC', files['CC'])


def test_ladder_uses_highest_complete_common_tier_for_available_rung():
    _, table = _table()
    species = SimpleNamespace(smiles='CCC', charge=0, mult=1, wellorts=0,
                              atom=['C'] * 3 + ['H'] * 8)
    energies = {
        smiles: {tier: ZeroKEnergy(smiles, value, tier, smiles + tier)
                 for tier in ('ANL0', 'L2')}
        for smiles, value in [('CCC', -118.), ('CC', -79.), ('C', -40.),
                              ('[H][H]', -1.)]
    }
    selected = select_cbh_ladder(species, energies, table)
    assert selected.reaction.rung == 1
    assert selected.formation.method == 'ANL0'
    del energies['CC']['ANL0']
    selected = select_cbh_ladder(species, energies, table)
    assert selected.formation.method == 'L2'
    assert any('ANL0' in reason for reason in selected.skipped)


def test_ladder_prefers_exact_reduced_recipe_before_lower_tier_or_rung():
    _, table = _table()
    species = SimpleNamespace(smiles='CCC', charge=0, mult=1, wellorts=0,
                              atom=['C'] * 3 + ['H'] * 8)
    full = 'profiled:ANL0-F12:scaled-triples:B2PLYP-D3BJ'
    reduced = full + ':reduced[without=higher_order]'
    profiles = {
        smiles: {
            reduced: ZeroKEnergy(smiles, value, reduced, smiles + '-reduced'),
            'L3': ZeroKEnergy(smiles, value + .2, 'L3', smiles + '-l3'),
        }
        for smiles, value in [('CCC', -118.), ('CC', -79.), ('C', -40.),
                              ('[H][H]', -1.)]
    }
    # A complete target profile is reusable, but the CBH-1 reaction cannot
    # use it unless every participant has that exact expression.
    profiles['CCC'][full] = ZeroKEnergy('CCC', -118.1, full, 'propane-full')
    selected = select_cbh_ladder(
        species, profiles, table, max_rung=1,
        method_ladder=(full, reduced, 'L3'))
    assert selected.reaction.rung == 1
    assert selected.formation.method == reduced
    assert selected.formation.references['CC'].method == reduced
    assert any(full in reason for reason in selected.skipped)


def test_ladder_descends_rung_when_no_common_profile_exists_at_higher_rung():
    _, table = _table()
    species = SimpleNamespace(smiles='CCC', charge=0, mult=1, wellorts=0,
                              atom=['C'] * 3 + ['H'] * 8)
    method = 'profiled:ANL0-F12:scaled-triples:B2PLYP-D3BJ'
    profiles = {
        smiles: {method: ZeroKEnergy(smiles, value, method, smiles + '-full')}
        for smiles, value in [('CCC', -118.), ('C', -40.), ('[H][H]', -1.)]
    }
    # Ethane is needed only by CBH-1. Its incompatible profile makes that
    # rung infeasible while every completed CBH-0 value remains reusable.
    profiles['CC'] = {
        'L2': ZeroKEnergy('CC', -79., 'L2', 'ethane-l2')}
    selected = select_cbh_ladder(
        species, profiles, table, max_rung=1,
        method_ladder=(method, 'L2'))
    assert selected.reaction.rung == 0
    assert selected.formation.method == method
    assert any('CBH-1' in reason for reason in selected.skipped)


def test_higher_rung_uses_computed_lower_rung_formation_from_its_own_tier():
    _, table = _table()
    species = SimpleNamespace(smiles='CCC', charge=0, mult=1, wellorts=0,
                              atom=['C'] * 3 + ['H'] * 8)
    energies = {
        'CCC': {'F12': ZeroKEnergy('CCC', -118., 'F12', 'propane-f12')},
        'CC': {
            'F12': ZeroKEnergy('CC', -79., 'F12', 'ethane-f12'),
            'ANL0': ZeroKEnergy('CC', -79.1, 'ANL0', 'ethane-anl0'),
        },
        'C': {
            'F12': ZeroKEnergy('C', -40., 'F12', 'methane-f12'),
            'ANL0': ZeroKEnergy('C', -40.1, 'ANL0', 'methane-anl0'),
        },
        '[H][H]': {
            'ANL0': ZeroKEnergy('[H][H]', -1., 'ANL0', 'hydrogen-anl0')},
    }
    selected = select_cbh_ladder(
        species, energies, table, max_rung=1,
        method_ladder=('F12', 'ANL0'))
    assert selected.formation.method == 'F12'
    assert selected.formation.references['CC'].method == 'ANL0'

    saved = asdict(selected.formation)
    assert formation_from_dict(saved) == selected.formation
    saved['formulas']['CCC']['H'] = 7
    with pytest.raises(ValueError, match='not element balanced'):
        formation_from_dict(saved)


def test_composite_file_profiles_expose_independent_rung_laddering(
        tmp_path, monkeypatch):
    files = {}
    for smiles, name, value, method in (
            ('CCC', 'propane-f12.json', -118., 'F12'),
            ('CC', 'ethane-f12.json', -79., 'F12'),
            ('C', 'methane-f12.json', -40., 'F12'),
            ('CC', 'ethane-anl0.json', -79.1, 'ANL0'),
            ('C', 'methane-anl0.json', -40.1, 'ANL0'),
            ('[H][H]', 'hydrogen-anl0.json', -1., 'ANL0')):
        path = tmp_path / name
        _composite_record(path, value, recipe=method)
        files.setdefault(smiles, {})[method] = path
    _, table = _table()
    monkeypatch.setattr('kinbot.anl.cbh.load_atct_references',
                        lambda *args, **kwargs: table)
    reaction, formation = solve_composite_records(
        'CCC', 1, files, atct_version='1.222', atct_cache=tmp_path)
    assert reaction.rung == 1
    assert formation.method == 'F12'
    assert formation.references['CC'].method == 'ANL0'
