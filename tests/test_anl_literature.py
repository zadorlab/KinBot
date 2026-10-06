"""Published numbers remain comparison targets rather than recipe inputs."""

import json
from pathlib import Path

import pytest

from kinbot.anl import literature
from kinbot.anl.dispatch import validate_spec
from kinbot.anl.literature import BENCHMARKS, compare_values


def test_ethane_tz_component_and_difference_regression_pass():
    benchmark = BENCHMARKS['ethane-tz-2017']
    observed = {
        task_id: target['expected']
        for task_id, target in benchmark['tasks'].items()
    }
    result = compare_values('ethane-tz-2017', observed)
    assert result['status'] == 'passed'
    assert result['checks']['delta_q_dz']['passed'] is True
    assert result['checks']['delta_q_dz']['derived_from'] == [
        'ccsdtq_dz', 'ccsdt_dz']
    assert 'inserted' in result['note']


def test_methane_targeted_mrcc_regression_detects_wrong_energy():
    benchmark = BENCHMARKS['methane-qz-2017']
    observed = {
        'ccsdtq_dz': benchmark['tasks']['ccsdtq_dz']['expected'],
        'ccsdtqp_dz': benchmark['tasks']['ccsdtqp_dz']['expected'] + 1e-4,
    }
    result = compare_values('methane-qz-2017', observed)
    assert result['status'] == 'failed'
    assert result['checks']['ccsdtq_dz']['passed'] is True
    assert result['checks']['ccsdtqp_dz']['passed'] is False
    assert result['checks']['delta_p_dz']['passed'] is False


def test_literature_comparison_requires_a_shared_completed_component():
    with pytest.raises(ValueError, match='no completed components'):
        compare_values('methyl-qz-2017', {'unrelated': -1.0})


def test_published_small_species_geometries_are_state_specific():
    methane = BENCHMARKS['methane-qz-2017']['molecule']
    methyl = BENCHMARKS['methyl-qz-2017']['molecule']
    assert (len(methane['symbols']), methane['multiplicity']) == (5, 1)
    assert (len(methyl['symbols']), methyl['multiplicity']) == (4, 2)
    assert methane['positions'][1][0] == pytest.approx(0.8882878215)
    assert methyl['positions'][1][0] == pytest.approx(1.0777376714)


def test_prepare_cli_builds_only_the_default_small_species_pair_plus_low(
        monkeypatch, tmp_path):
    captured = {}

    def fake_prepare(spec, run_dir):
        captured['spec'] = spec
        return Path(run_dir)

    monkeypatch.setattr('kinbot.anl.dispatch.prepare', fake_prepare)
    run_dir = tmp_path / 'methane'
    assert literature.main([
        'prepare-higher-order', 'methane-qz-2017', str(run_dir)]) == 0
    spec = captured['spec']
    assert [task['id'] for task in spec['tasks']] == [
        'ccsdt_dz', 'ccsdtq_dz', 'ccsdtqp_dz']
    assert spec['molecule']['multiplicity'] == 1
    qp = next(task for task in spec['tasks'] if task['id'] == 'ccsdtqp_dz')
    assert qp['resources']['walltime'] == '24:00:00'
    assert qp['resources']['max_cores'] == 4
    assert spec['intent']['literature_benchmark']['name'] == \
        'methane-qz-2017'
    json.dumps(spec)


def test_prepare_methyl_source_benchmark_uses_its_published_uhf_reference(
        monkeypatch, tmp_path):
    captured = {}

    def fake_prepare(spec, run_dir):
        captured['spec'] = spec
        return Path(run_dir)

    monkeypatch.setattr('kinbot.anl.dispatch.prepare', fake_prepare)
    assert literature.main([
        'prepare-higher-order', 'methyl-qz-2017', str(tmp_path / 'methyl'),
        '--task', 'ccsdt_dz', '--task', 'ccsdtq_dz']) == 0
    tasks = {task['id']: task for task in captured['spec']['tasks']}
    assert tasks['ccsdt_dz']['result_parser']['reference'] == 'ROHF'
    assert tasks['ccsdtq_dz']['result_parser']['reference'] == 'UHF'
    assert 'scftype=UHF' in tasks['ccsdtq_dz']['input_template']
    assert 'rohftype=semicanonical' not in \
        tasks['ccsdtq_dz']['input_template']
    for task in captured['spec']['tasks']:
        task['resources'].update(
            cores=4, memory_mb=64000, partition='test')
    validate_spec(captured['spec'])


def test_run_comparison_rejects_a_different_reference(monkeypatch):
    task = {
        'id': 'ccsdtqp_dz', 'backend': 'mrcc',
        'result_parser': {
            'method': 'CCSDTQ(P)', 'basis': 'cc-pVDZ',
            'reference': 'UHF', 'correlation': 'unrestricted'}}
    monkeypatch.setattr(
        literature, '_load',
        lambda run_dir: (Path(run_dir), {
            'molecule':
                literature.BENCHMARKS['methane-qz-2017']['molecule'],
            'tasks': [task]}, {
            'tasks': {'ccsdtqp_dz': {'status': 'complete'}}}))
    with pytest.raises(ValueError, match='does not match'):
        literature.compare_run('methane-qz-2017', '/synthetic/run')


def test_formation_comparison_uses_computed_cbh_result(tmp_path):
    expected = BENCHMARKS['methyl-qz-2017'][
        'formation_enthalpy_0k_kcal_mol']['ANL0-F12']
    result = tmp_path / 'methyl-cbh.json'
    result.write_text(json.dumps({
        'schema': 1, 'status': 'complete',
        'formation': {
            'target_smiles': '[CH3]', 'method': 'ANL0-F12',
            'formation_0k_kj_mol': (expected + 0.02) * 4.184,
        },
    }))
    comparison = literature.compare_formation(
        'methyl-qz-2017', result, absolute_tolerance_kcal_mol=0.1)
    assert comparison['status'] == 'passed'
    assert comparison['error_kcal_mol'] == pytest.approx(0.02)
    assert 'inserted' in comparison['note']

    profile = 'profiled:ANL0-F12:scaled-triples:B2PLYP-D3BJ'
    payload = json.loads(result.read_text())
    payload['formation']['method'] = profile
    result.write_text(json.dumps(payload))
    comparison = literature.compare_formation(
        'methyl-qz-2017', result, method=profile,
        absolute_tolerance_kcal_mol=0.1)
    assert comparison['status'] == 'passed'
    assert comparison['benchmark_method_family'] == 'ANL0-F12'
    assert comparison['profiled_variant'] is True


def test_formation_comparison_rejects_wrong_species_or_method(tmp_path):
    result = tmp_path / 'wrong-cbh.json'
    result.write_text(json.dumps({
        'schema': 1, 'status': 'complete',
        'formation': {
            'target_smiles': 'C', 'method': 'ANL0',
            'formation_0k_kj_mol': 0.,
        },
    }))
    with pytest.raises(ValueError, match='target'):
        literature.compare_formation('methyl-qz-2017', result)
