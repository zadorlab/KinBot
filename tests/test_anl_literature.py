"""Published numbers remain comparison targets rather than recipe inputs."""

import json
from pathlib import Path

import pytest

from kinbot.anl import literature
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
