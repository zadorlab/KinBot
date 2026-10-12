"""Final CBH/ANL, ROTD_py, MESSPF, and MESS orchestration checks."""

import json
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot.anl.end_to_end import (_audit_handoff, _literature_checks,
                                   _write_mess_input, component_breakdown)


def _composite(path, offset):
    values = {
        'reference_cbs': -40.0 + offset,
        'hoe_high': -39.8 + offset,
        'hoe_low': -39.7 + offset,
        'core_valence_cbs': -0.02 + offset,
        'scalar_relativistic': -0.005 + offset,
        'dboc': 0.002 + offset,
        'spin_orbit': offset,
        'harmonic_zpe': 0.1 + offset,
        'vpt2_correction': -0.01 + offset,
    }
    path.write_text(json.dumps({
        'schema': 1, 'status': 'complete',
        'components': {
            key: {'value_hartree': value} for key, value in values.items()},
    }))
    return values


def test_component_breakdown_uses_each_anl_term_once(tmp_path):
    left_path, right_path = tmp_path / 'left.json', tmp_path / 'right.json'
    left = _composite(left_path, 0.)
    right = _composite(right_path, .001)
    reaction = SimpleNamespace(stoichiometry={'left': -1, 'right': 1})
    result = component_breakdown(
        reaction, {'left': left_path, 'right': right_path})
    assert result['reference_cbs']['hartree'] == pytest.approx(.001)
    # The identical offset cancels inside each species' high-low correction.
    assert result['higher_order']['hartree'] == pytest.approx(0.)
    expected = ((right['reference_cbs'] - left['reference_cbs'])
                + ((right['hoe_high'] - right['hoe_low'])
                   - (left['hoe_high'] - left['hoe_low']))
                + sum(right[key] - left[key] for key in (
                    'core_valence_cbs', 'scalar_relativistic', 'dboc',
                    'spin_orbit', 'harmonic_zpe', 'vpt2_correction')))
    assert result['total_zero_k_reaction']['hartree'] == pytest.approx(expected)


def test_component_breakdown_reads_direct_rank_exact_zero(tmp_path):
    left_path, right_path = tmp_path / 'left.json', tmp_path / 'right.json'
    left = _composite(left_path, 0.)
    right = _composite(right_path, .001)
    for path in (left_path, right_path):
        record = json.loads(path.read_text())
        record['components'].pop('hoe_high')
        record['components'].pop('hoe_low')
        record['components']['higher_order_correction'] = {
            'value_hartree': 0.}
        path.write_text(json.dumps(record))
    reaction = SimpleNamespace(stoichiometry={'left': -1, 'right': 1})
    result = component_breakdown(
        reaction, {'left': left_path, 'right': right_path})
    assert result['higher_order']['hartree'] == 0.
    expected = ((right['reference_cbs'] - left['reference_cbs'])
                + sum(right[key] - left[key] for key in (
                    'core_valence_cbs', 'scalar_relativistic', 'dboc',
                    'spin_orbit', 'harmonic_zpe', 'vpt2_correction')))
    assert result['total_zero_k_reaction']['hartree'] == pytest.approx(expected)


def test_mess_continuation_disables_new_vrc_work_and_requires_production(
        tmp_path):
    run = tmp_path / 'run'
    run.mkdir()
    source = run / 'production.json'
    source.write_text(json.dumps({
        'title': 'fixture', 'me': 0, 'run_me': 0,
        'vrc_tst_scan': {'parent': ['reaction']}, 'rotdpy_run': 1,
    }))
    rotdpy = run / 'rotdPy' / 'reaction.py'
    rotdpy.parent.mkdir()
    rotdpy.write_text('# fixture\n')
    handoff = tmp_path / 'handoff.json'
    handoff.write_text('{}\n')
    manifest = {
        'run_dir': str(run), 'kinbot_input': str(source),
        'mess': {'mess_command': 'mess --flag',
                 'messpf_command': 'messpf', 'collider': 'He',
                 'epsilon': 252.3, 'epsilon_unit': 'K', 'sigma': 4.302,
                 'queue_job_limit': 8},
    }
    actual_run, output = _write_mess_input(
        tmp_path, manifest, handoff.resolve(), rotdpy.resolve())
    generated = json.loads(output.read_text())
    assert actual_run == run.resolve()
    assert generated['me'] == generated['run_me'] == 1
    assert generated['run_messpf'] == 1
    assert generated['anl_handoff_manifest'] == str(handoff.resolve())
    assert generated['rotdpy_directory'] == str(rotdpy.parent.resolve())
    assert generated['rotdpy_mess_mode'] == 'production'
    assert generated['rotdpy_run'] == 0
    assert generated['vrc_tst_scan'] == 0
    assert generated['vrc_tst_noscan'] == 0
    assert generated['high_level'] == 1
    assert generated['allow_l2_without_conf'] == 1
    assert generated['rotor_scan'] == 1
    assert generated['collider'] == 'He'
    assert generated['epsilon'] == 252.3
    assert generated['epsilon_unit'] == 'K'
    assert generated['sigma'] == 4.302
    assert generated['queue_job_limit'] == 8


def test_production_restart_script_does_not_reenter_completed_kinbot():
    root = Path(__file__).resolve().parents[1]
    script = (root / 'examples/anl/ethane_profiled_hpc/'
              'run_vrc_production.sh').read_text()
    assert 'if ! test -f "$rotdpy_result"; then\n' in script
    assert '"$python_bin" -m kinbot.kb ethane_vrc_production.json\nfi' in script


def test_final_manifest_defaults_exclude_superseded_ethane_graphs():
    root = Path(__file__).resolve().parents[1]
    script = (root / 'examples/anl/ethane_profiled_hpc/'
              'full_verification_manifest.py').read_text()
    assert "base / 'anl_current_base_ethane_449cfbd'" in script
    assert "base / 'anl_post_geometry_ethane_449cfbd'" in script
    assert 'anl_current_base_ethane_fb3eab2' not in script
    assert 'base=methane_interface' not in script
    assert 'base=hydrogen_interface' not in script


def test_final_manifest_assigns_each_completed_graph_to_a_valid_role(tmp_path):
    root = Path(__file__).resolve().parents[1]
    source = (root / 'examples/anl/ethane_profiled_hpc/'
              'full_verification_manifest.py')
    spec = importlib.util.spec_from_file_location(
        'kinbot_full_verification_manifest_test', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    base = tmp_path / 'completed'
    refs = base / 'cbh0_reference_graphs'
    run = base / 'vrc_production_run'
    (run / 'rotdPy').mkdir(parents=True)
    (run / 'ethane_vrc_production.json').write_text('{}\n')
    (run / 'rotdPy' / f'{module.REACTION}.py').write_text('# completed\n')
    workflows = {
        base / 'anl_interface':
            'anl1-f12-non-mrcc-interface-validation',
        base / 'anl_current_base_ethane_449cfbd':
            'anl-current-base-validation',
        base / 'anl_post_geometry_ethane_449cfbd':
            'anl-post-geometry-validation',
        refs / 'methane_interface':
            'anl1-f12-non-mrcc-interface-validation',
        refs / 'methane_post': 'anl-post-geometry-validation',
        refs / 'hydrogen_interface':
            'anl1-f12-non-mrcc-interface-validation',
        refs / 'hydrogen_post': 'anl-post-geometry-validation',
        refs / 'methyl_composite': 'anl-composite-validation',
        refs / 'atomic_h_composite': 'anl-composite-validation',
        refs / 'methane_source_qz':
            'anl1-higher-order-interface-validation',
        refs / 'methyl_source_qz':
            'anl1-higher-order-interface-validation',
    }
    for directory, name in workflows.items():
        directory.mkdir(parents=True, exist_ok=True)
        (directory / 'workflow.json').write_text(json.dumps({
            'schema': 1, 'name': name,
        }))
    output = base / 'manifest.json'
    module.main([str(base), '--output', str(output)])
    generated = json.loads(output.read_text())
    assert generated['composites']['ethane']['assembly'][
        'base_run'].endswith('anl_current_base_ethane_449cfbd')
    assert 'base_run' not in generated['composites']['methane']['assembly']
    assert 'base_run' not in generated['composites']['hydrogen']['assembly']


def test_applied_handoff_audit_verifies_every_source_hash(tmp_path):
    run = tmp_path / 'run'
    (run / 'me').mkdir(parents=True)
    composite = tmp_path / 'composite.json'
    formation = tmp_path / 'formation.json'
    composite.write_text('{"composite": true}\n')
    formation.write_text('{"formation": true}\n')
    manifest = tmp_path / 'handoff.json'
    manifest.write_text(json.dumps({
        'schema': 1, 'mode': 'cbh-anl', 'strict': True,
        'species': {'species-A': {}}, 'transition_states': {},
    }))
    from kinbot.anl.end_to_end import _sha256
    (run / 'me' / 'anl_handoff.json').write_text(json.dumps({
        'schema': 1, 'status': 'complete', 'strict': True,
        'manifest_sha256': _sha256(manifest),
        'species': {'species-A': {
            'method': 'profiled:ANL0-F12', 'zero_k_hartree': -1.,
            'formation_0k_kj_mol': 2.,
            'composite': str(composite),
            'composite_sha256': _sha256(composite),
            'formation': str(formation),
            'formation_sha256': _sha256(formation),
        }}, 'transition_states': {},
    }))
    data = {'handoff': {'output': str(manifest)}}
    result = _audit_handoff(tmp_path, data, run)
    assert result['status'] == 'passed'
    formation.write_text('changed\n')
    with pytest.raises(ValueError, match='formation file failed hash'):
        _audit_handoff(tmp_path, data, run)


def test_literature_gate_requires_every_declared_component(
        tmp_path, monkeypatch):
    monkeypatch.setattr(
        'kinbot.anl.end_to_end.compare_run',
        lambda benchmark, source: {
            'status': 'passed', 'checks': {'ccsdt_dz': {'passed': True}}})
    with pytest.raises(ValueError, match='lacks required comparisons'):
        _literature_checks(tmp_path, [{
            'kind': 'run', 'benchmark': 'fixture', 'source': str(tmp_path),
            'required_checks': ['ccsdt_dz', 'ccsdtq_dz'],
        }])
