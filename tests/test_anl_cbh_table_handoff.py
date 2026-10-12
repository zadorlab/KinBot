"""Local workflow checks using rounded CBH rows supplied with the ANL plan.

These are deliberately tiny fixtures. They check reaction signs, tier selection,
and KinBot-to-MESS serialization without claiming to validate a QC method.
"""

import json
from pathlib import Path
import shlex
import sys
from types import SimpleNamespace
from unittest.mock import patch

from ase.build import molecule
import pytest

from kinbot import constants
from kinbot.anl.atct import parse_atct_html
from kinbot.anl.cbh import (CBHReaction, DEFAULT_METHOD_LADDER, FormationEnthalpy,
                            HARTREE_TO_KJ_MOL,
                            generate_cbh_reaction, select_cbh_ladder,
                            solve_formation_enthalpy)
from kinbot.anl.thermochemistry import (_count_hindered_rotor_blocks,
                                        audit_messpf_run)
from kinbot.energy import ZeroKEnergy, attach_formation_enthalpy
from kinbot.mess import MESS
from kinbot.parameters import Parameters
from tests.test_anl_cbh_atct import _row
from tests.test_mess_conformers import values
from tests.test_mess_conformers import point


KCAL_TO_KJ = 4.184


def _reference_model_library(path):
    ch4 = molecule('CH4')
    h2 = molecule('H2')
    models = {
        'C': {
            'name': 'methane', 'smiles': 'C',
            'symbols': ch4.get_chemical_symbols(),
            'positions': ch4.get_positions().tolist(),
            'charge': 0, 'multiplicity': 1, 'symmetry_factor': 12.,
            'frequencies_cm_inverse': [1300.] * 9,
            'source': {'kind': 'test-fixture'},
        },
        '[H][H]': {
            'name': 'hydrogen', 'smiles': '[H][H]',
            'symbols': h2.get_chemical_symbols(),
            'positions': h2.get_positions().tolist(),
            'charge': 0, 'multiplicity': 1, 'symmetry_factor': 2.,
            'frequencies_cm_inverse': [4400.],
            'source': {'kind': 'test-fixture'},
        },
    }
    path.write_text(json.dumps({'schema': 1, 'species': models}))
    return path


def _figure_references():
    # The paper table uses ATcT 1.202 and prints only rounded kcal/mol values.
    raw = ('ATcT enthalpies of formation based on version 1.202 of the '
           'Thermochemical Network' +
           _row(1, '1333-74-0*0', 'Dihydrogen', 'H2 (g)', '[H][H]', '0', 'exact') +
           _row(2, '74-82-8*0', 'Methane', 'CH4 (g)', 'C',
                str(-15.903 * KCAL_TO_KJ)) +
           _row(6, '7664-39-3*0', 'Hydrogen fluoride', 'HF (g)', 'F',
                str(-65.172 * KCAL_TO_KJ)) +
           _row(3, '74-84-0*0', 'Ethane', 'C2H6 (g)', 'CC',
                str(-16.461 * KCAL_TO_KJ)) +
           _row(4, '593-53-3*0', 'Fluoromethane', 'CH3F (g)', 'CF',
                str(-54.346 * KCAL_TO_KJ)) +
           _row(5, '12385-13-6*0', 'Atomic hydrogen', 'H (g)', '[H]',
                str(51.633 * KCAL_TO_KJ))).encode()
    return parse_atct_html(raw, '1.202')


def _energies_for_delta(reaction, delta_kcal, method):
    """Choose one target E0 that makes the signed reaction delta exact."""
    carbon, hydrogen, fluorine = -40., -1., -100.
    ethane_delta = -15.345 * KCAL_TO_KJ / HARTREE_TO_KJ_MOL
    fluoromethane_delta = -26.729 * KCAL_TO_KJ / HARTREE_TO_KJ_MOL
    reference_energies = {
        'C': carbon, '[H][H]': hydrogen, 'F': fluorine,
        'CC': 2 * carbon - hydrogen - ethane_delta,
        'CF': carbon + fluorine - hydrogen - fluoromethane_delta,
    }
    reference_sum = sum(coefficient * reference_energies[smiles]
                        for smiles, coefficient in reaction.stoichiometry.items()
                        if smiles != reaction.target_smiles)
    target = reference_sum - delta_kcal * KCAL_TO_KJ / HARTREE_TO_KJ_MOL
    return {smiles: ZeroKEnergy(smiles,
                                target if smiles == reaction.target_smiles
                                else reference_energies[smiles],
                                method, f'{smiles}:{method}:native-output')
            for smiles in reaction.stoichiometry}


def test_ethane_cbh0_table_components_and_heat_of_formation():
    reaction = generate_cbh_reaction('CC', 0)
    assert reaction.stoichiometry == {'C': 2, 'CC': -1, '[H][H]': -1}
    # The table prints the reverse sign of the displayed forward CBH arrow.
    printed_components = [18.080, .073, .006, .098, .006, .022,
                          0., -3.033, .093]
    assert sum(printed_components) == pytest.approx(15.345)
    energies = _energies_for_delta(reaction, -15.345, 'ANL0-F12')
    formation = solve_formation_enthalpy(reaction, energies, _figure_references())
    assert formation.reaction_energy_0k_kj_mol / KCAL_TO_KJ == pytest.approx(-15.345)
    assert formation.formation_0k_kj_mol / KCAL_TO_KJ == pytest.approx(-16.461)
    assert formation.atct_version == '1.202'


def test_fluorinated_cbh1_table_reaction_and_full_fallback_ladder():
    reaction = generate_cbh_reaction('CC(F)(F)F', 1)
    assert reaction.stoichiometry == {
        'CC(F)(F)F': -1, 'C': -3, 'CC': 1, 'CF': 3}
    table = _figure_references()
    species = SimpleNamespace(smiles='CC(F)(F)F', charge=0, mult=1,
                              wellorts=0, atom=['C'] * 2 + ['H'] * 3 + ['F'] * 3)
    labels = DEFAULT_METHOD_LADDER[:-2] + ('L3:small-points', 'L2')
    pools = {}
    for label in labels:
        energies = _energies_for_delta(reaction, 44.938, label)
        for child in (generate_cbh_reaction('CC', 0),
                      generate_cbh_reaction('CF', 0)):
            child_energies = _energies_for_delta(
                child, -15.345 if child.target_smiles == 'CC' else -26.729,
                label)
            energies.update(child_energies)
        for smiles, energy in energies.items():
            pools.setdefault(smiles, {})[label] = energy
    # A single missing accepted reference energy steps the entire reaction down.
    for label in labels:
        selected = select_cbh_ladder(species, pools, table, max_rung=1)
        assert selected.reaction.rung == 1
        assert selected.formation.method == label
        assert selected.formation.formation_0k_kj_mol / KCAL_TO_KJ == pytest.approx(
            -176.728)
        del pools['CF'][label]
    selected = select_cbh_ladder(species, pools, table, max_rung=1)
    assert selected.reaction.rung == 0
    assert any('CBH-1' in reason for reason in selected.skipped)


def test_perfluoroalkane_cbh2_and_cbh3_table_reactions():
    c2f6 = generate_cbh_reaction('FC(F)(F)C(F)(F)F', 2)
    assert c2f6.stoichiometry == {
        c2f6.target_smiles: -1, 'CC': -1, 'CC(F)(F)F': 2}
    c3f8 = generate_cbh_reaction('FC(F)(F)C(F)(F)C(F)(F)F', 3)
    assert c3f8.stoichiometry == {
        c3f8.target_smiles: -1, 'CC(F)(F)C': -1,
        'FC(C(F)(F)C)(F)F': 2}


def test_methyl_radical_table_reaction_requires_explicit_spin_state():
    reaction = generate_cbh_reaction('[CH3]', 0, multiplicity=2)
    assert reaction.stoichiometry == {
        '[CH3]': -1, '[H][H]': -1, 'C': 1, '[H]': 1}
    assert reaction.states == {'[CH3]': (0, 2), '[H]': (0, 2)}
    energy_values = {'[H][H]': -1., 'C': -40., '[H]': -.5}
    reference_sum = sum(coefficient * energy_values[smiles]
                        for smiles, coefficient in reaction.stoichiometry.items()
                        if smiles != '[CH3]')
    target = reference_sum + .096 * KCAL_TO_KJ / HARTREE_TO_KJ_MOL
    energies = {smiles: ZeroKEnergy(
        smiles, target if smiles == '[CH3]' else energy_values[smiles],
        'ANL0-F12', f'{smiles}:output', 0,
        2 if smiles in ('[CH3]', '[H]') else 1)
        for smiles in reaction.stoichiometry}
    with pytest.raises(ValueError, match='explicit ID'):
        solve_formation_enthalpy(reaction, energies, _figure_references())
    result = solve_formation_enthalpy(
        reaction, energies, _figure_references(),
        reference_ids={'[H]': '12385-13-6*0'})
    assert result.reaction_energy_0k_kj_mol / KCAL_TO_KJ == pytest.approx(-.096)
    assert result.formation_0k_kj_mol / KCAL_TO_KJ == pytest.approx(35.826)


def test_rounded_ethane_formation_is_handed_to_fake_mess(tmp_path, monkeypatch):
    reaction = generate_cbh_reaction('CC', 0)
    energies = _energies_for_delta(reaction, -15.345, 'ANL0-F12')
    formation = solve_formation_enthalpy(reaction, energies, _figure_references())
    geometry = molecule('C2H6')
    well = SimpleNamespace(
        name='ethane', chemid=123, smiles='CC', charge=0, mult=1,
        natom=len(geometry), atom=geometry.get_chemical_symbols(),
        geom=geometry.positions, energy=-79., zpe=.01, mass=30.07,
        sigma_ext=1, nopt=1, reduced_freqs=[1000.] * 18,
        conformer_index=[], reac_obj=[], reac_ts_done=[], reac_type=[],
        final_zero_k_energy=energies['CC'])
    attach_formation_enthalpy(well, formation)
    monkeypatch.chdir(tmp_path)
    Path('me').mkdir()
    Path('input.json').write_text(json.dumps({
        'smiles': 'CC', 'me': 0, 'uq': 0, 'rotor_scan': 1,
        'multi_conf_tst': 0, 'epsilon': 100., 'sigma': 3.,
        'barrier_threshold': 100., 'high_level': 1,
        'allow_l2_without_conf': 1,
        'cbh_thermochemistry_models': str(
            _reference_model_library(tmp_path / 'cbh_models.json'))}))
    par = Parameters('input.json', show_warnings=False).par
    writer = MESS(par, well)
    with patch.object(MESS, 'make_rotors', return_value='! L2 HIR fixture') as rotors:
        writer.write_input(None)
    # The writer still calls its L2 HIR path while using the accepted ANL E0.
    rotors.assert_called()
    output = Path('me/mess_0000.inp').read_text()
    assert '! L2 HIR fixture' in output
    assert 'CBH/ANL formation enthalpies and barriers' in output
    assert values(output, 'ZeroEnergy') == [0.]
    partition_input = Path('me/partition_functions/123.inp').read_text()
    assert '298.15' in partition_input
    assert 'ZeroEnergy[kcal/mol]               0.0' in partition_input
    assert '! L2 HIR fixture' in partition_input
    assert '    +++\n' in partition_input
    assert partition_input.rstrip().endswith('  End')
    assert not any(line == 'End' for line in partition_input.splitlines())
    atomic_input = writer.write_partition_function_input(
        SimpleNamespace(natom=1, name='hydrogen', atom=['H'], mult=2),
        Path('me/partition_functions/H.inp')).read_text()
    assert atomic_input.rstrip().endswith('  End')
    assert not any(line == 'End' for line in atomic_input.splitlines())
    runner = Path('me/partition_functions/run_messpf.sh')
    assert runner.stat().st_mode & 0o111
    assert 'KINBOT_MESSPF_COMMAND' in runner.read_text()
    sidecar = json.loads(Path('me/formation_0k.json').read_text())
    record = sidecar['species']['123']
    assert record['formation_0k_kj_mol'] / KCAL_TO_KJ == pytest.approx(-16.461)
    assert record['zero_k_energy_hartree'] == energies['CC'].hartree
    assert record['method'] == 'ANL0-F12'
    assert record['cbh_rung'] == 0
    assert record['reference_ids']['C'] == '74-82-8*0'
    pf_output = Path('me/partition_functions/123.dat')
    rows = '\n'.join(
        f'{temperature:.2f} 11.0 0.010 0.0 54.0 10.0'
        for temperature in (200., 250., 298.15, 400., 600., 800., 1000.,
                            1400., 1800., 2400., 3000.))
    pf_output.write_text('''
Natural log of the partition function, its derivatives, entropy, and thermal capacity:
T, K       ethane       ethane       ethane       ethane       ethane
                   Z_0          Z_1          Z_2 S, cal/mol/K C, cal/mol/K
''' + rows + '\n')
    manifest = json.loads(
        Path('me/partition_functions/manifest.json').read_text())
    for ident, entry in manifest['species'].items():
        if ident != '123':
            Path(entry['output']).write_text(
                'Z_0 Z_1 Z_2 S, cal/mol/K C, cal/mol/K\n' + rows + '\n')
    thermo = writer.read_partition_function_outputs()
    assert thermo['species']['123']['formation_0k_kj_mol'] == pytest.approx(
        formation.formation_0k_kj_mol)
    assert Path('me/thermochemistry_298.json').is_file()
    assert thermo['species']['123']['nasa7']['format'] == 'NASA7'


def test_messpf_executor_replaces_stale_output_and_records_hashes(
        tmp_path, monkeypatch):
    reaction = generate_cbh_reaction('CC', 0)
    energies = _energies_for_delta(reaction, -15.345, 'ANL0-F12')
    formation = solve_formation_enthalpy(
        reaction, energies, _figure_references())
    geometry = molecule('C2H6')
    well = SimpleNamespace(
        name='ethane', chemid=123, smiles='CC', charge=0, mult=1,
        natom=len(geometry), atom=geometry.get_chemical_symbols(),
        geom=geometry.positions, energy=-79., zpe=.01, mass=30.07,
        sigma_ext=1, nopt=1, reduced_freqs=[1000.] * 18,
        conformer_index=[], reac_obj=[], reac_ts_done=[], reac_type=[],
        final_zero_k_energy=energies['CC'])
    attach_formation_enthalpy(well, formation)
    monkeypatch.chdir(tmp_path)
    Path('me').mkdir()
    rows = '\n'.join(
        f'{temperature:.2f} 11.0 0.010 0.0 54.0 10.0'
        for temperature in (200., 250., 298.15, 400., 600., 800., 1000.,
                            1400., 1800., 2400., 3000.))
    binary = tmp_path / 'fake messpf.py'
    binary.write_text(
        f'#!{sys.executable}\n'
        'from pathlib import Path\n'
        'import sys\n'
        'source = Path(sys.argv[-1])\n'
        'source.with_suffix(".dat").write_text(' + repr(
            'Z_0 Z_1 Z_2 entropy heat_capacity other\n' + rows + '\n') + ')\n'
        'print("fake MESSPF completed")\n')
    binary.chmod(0o755)
    Path('input.json').write_text(json.dumps({
        'smiles': 'CC', 'me': 0, 'uq': 0, 'rotor_scan': 0,
        'multi_conf_tst': 0, 'epsilon': 100., 'sigma': 3.,
        'barrier_threshold': 100., 'high_level': 1,
        'allow_l2_without_conf': 1,
        'messpf_command': shlex.join([str(binary), '--fixture']),
        'cbh_thermochemistry_models': str(
            _reference_model_library(tmp_path / 'cbh_models.json'))}))
    par = Parameters('input.json', show_warnings=False).par
    writer = MESS(par, well)
    with patch.object(MESS, 'make_rotors', return_value=''):
        writer.write_input(None)
    output = Path('me/partition_functions/123.dat')
    output.write_text('stale output must not be accepted\n')
    thermo = writer.run_partition_functions()
    assert thermo['species']['123']['formation_0k_kj_mol'] == pytest.approx(
        formation.formation_0k_kj_mol)
    execution = json.loads(
        Path('me/partition_functions/execution.json').read_text())
    record = execution['species']['123']
    assert record['status'] == 'complete'
    assert record['command'][-2:] == ['--fixture', '123.inp']
    assert len(record['input_sha256']) == len(record['output_sha256']) == 64
    assert Path('me/partition_functions/123.stdout').read_text() == (
        'fake MESSPF completed\n')
    audit = audit_messpf_run(tmp_path)
    assert audit['status'] == 'passed'
    assert audit['species']['123']['output_sha256'] == record['output_sha256']

    binary.write_text(
        f'#!{sys.executable}\n'
        'from pathlib import Path\n'
        'import sys\n'
        'source = Path(sys.argv[-1])\n'
        'source.with_suffix(".log").write_text("native parser detail\\n")\n'
        'print("native stderr detail", file=sys.stderr)\n'
        'raise SystemExit(6)\n')
    with pytest.raises(RuntimeError) as caught:
        writer.run_partition_functions()
    message = str(caught.value)
    assert 'native parser detail' in message
    assert 'native stderr detail' in message
    assert '123.log (last 40 lines)' in message


def test_messpf_audit_recognizes_native_hindered_rotor_grammar():
    text = '''
      Rotor     Hindered
        Group 3 4 5
      End ! Rotor
      Rotor Free
      End ! Rotor
    '''
    assert _count_hindered_rotor_blocks(text) == 1


def test_direct_mess_network_accepts_per_species_ladder_fallback(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('me').mkdir()
    Path('input.json').write_text(json.dumps({
        'smiles': 'O', 'me': 0, 'uq': 0, 'rotor_scan': 0,
        'multi_conf_tst': 0, 'epsilon': 100., 'sigma': 3.,
        'barrier_threshold': 100., 'high_level': 1,
        'allow_l2_without_conf': 1}))
    par = Parameters('input.json', show_warnings=False).par
    root, product, ts = (point('water', 123), point('other', 789),
                         point('saddle', 456, True))
    for item, energy, method in (
            (root, -75.99, 'ANL1-F12'),
            (product, -75.98, 'ANL0-F12'),
            (ts, -75.95, 'ANL1-F12')):
        item.charge = 0
        item.final_zero_k_energy = ZeroKEnergy('O', energy, method,
                                              item.name + ':native-output')
    for item, formation_kj in ((root, -100.), (product, -50.)):
        final = item.final_zero_k_energy
        attach_formation_enthalpy(item, FormationEnthalpy(
            final.smiles, 0, final.method, 0., formation_kj, {}, 'test',
            'fixture', {final.smiles: final.source}))
    root.mass = 18.
    reaction = SimpleNamespace(ts=ts, products=[product],
                               instance_name='saddle', do_vdW=False,
                               prod_opt=[SimpleNamespace(species=product)], mp2=0)
    root.reac_obj, root.reac_ts_done = [reaction], [-1]
    with patch.object(MESS, 'make_rotors', return_value=''):
        MESS(par, root).write_input(None)
    output = Path('me/mess_0000.inp').read_text()
    assert sorted(values(output, 'ZeroEnergy')) == pytest.approx(
        sorted([0., round(.04 * constants.AUtoKCAL, 2),
                round(50. / KCAL_TO_KJ, 2)]))
    # The product well is 50 kJ/mol above the root from Hf(0), independent of
    # the different accepted electronic-energy ladder rung.
    assert round(50. / KCAL_TO_KJ, 2) in values(output, 'ZeroEnergy')
    product.formation_enthalpy_0k = None
    with pytest.raises(ValueError, match='formation enthalpy'):
        MESS(par, root).write_input(None)


def test_cbh_anl_homolytic_channel_uses_fragment_threshold_without_saddle(
        tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path('me').mkdir()
    Path('input.json').write_text(json.dumps({
        'smiles': 'O', 'me': 0, 'uq': 0, 'rotor_scan': 0,
        'multi_conf_tst': 0, 'epsilon': 100., 'sigma': 3.,
        'barrier_threshold': 100., 'high_level': 1,
        'allow_l2_without_conf': 1}))
    par = Parameters('input.json', show_warnings=False).par
    root, product_a, product_b = (point('water', 123), point('a', 789),
                                  point('b', 790))
    for item, energy, method, formation_kj in (
            (root, -75.99, 'ANL1-F12', -100.),
            (product_a, -37.1, 'ANL0-F12', -40.),
            (product_b, -38.2, 'L3', -40.)):
        item.charge = 0
        item.final_zero_k_energy = ZeroKEnergy(
            'O', energy, method, item.name + ':native-output')
        attach_formation_enthalpy(item, FormationEnthalpy(
            'O', 0, method, 0., formation_kj, {}, 'test', 'fixture',
            {'O': item.final_zero_k_energy.source}))
    root.mass = 18.
    ts = point('unused_saddle', 456, True)
    reaction = SimpleNamespace(
        ts=ts, products=[product_a, product_b], instance_name='hom_sci_fixture',
        do_vdW=False,
        prod_opt=[SimpleNamespace(species=product_a),
                  SimpleNamespace(species=product_b)], mp2=0)
    root.reac_obj, root.reac_ts_done, root.reac_type = [reaction], [-1], ['hom_sci']
    with patch.object(MESS, 'make_rotors', return_value=''):
        MESS(par, root).write_input(None)
    output = Path('me/mess_0000.inp').read_text()
    assert '{blessname}' not in output
    assert 'Barrier       bl_ts_1 w_1 b_1' in output
    assert 'Barrier       ts_1 ' not in output
    assert values(output, 'GroundEnergy') == [round(20. / KCAL_TO_KJ, 2)]
    sidecar = json.loads(Path('me/formation_0k.json').read_text())
    threshold = sidecar['barriers']['hom_sci_fixture']
    assert threshold['kind'] == 'fragment_channel_threshold'
    assert threshold['relative_0k_kj_mol'] == pytest.approx(20.)
