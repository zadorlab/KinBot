"""Slurm inventory and atomic live capacity changes, with retained allocations."""
from types import SimpleNamespace

import pytest

from kinbot import run_capacity, run_policy
from tests.test_run_policy import pool, submit


NODES = '''node1|drained*|96|126000|up
node2|idle|96|126000|up
node3|allocated|96|126000|up
node4|mixed|96|126000|up
node5|idle|4|126000|up
node6|idle|96|64000|up
node7|idle|96|126000|down
node8|idle*|96|126000|up
node9|idle|96|126000|up
'''


def test_unique_compatible_nodes_preserve_exclusions_and_flags():
    result = run_capacity.node_inventory(NODES+NODES, ['node9'], 8, 98000)
    assert result['unique_nodes'] == 9
    assert result['eligible'] == ['node2', 'node3', 'node4']
    assert result['idle'] == ['node2']
    assert result['rejected'] == dict(node1='drained*', node5='resources', node6='resources',
        node7='partition_unavailable', node8='idle*', node9='excluded')


@pytest.mark.parametrize('text', ['', 'node[1-3]|idle|96|126000|up', 'node2|idle|x|126000|up',
    NODES+'node2|allocated|96|126000|up\n'])
def test_bad_inventory_never_becomes_capacity(text):
    with pytest.raises(ValueError): run_capacity.node_inventory(text, [], 8, 98000)


def test_fit_keeps_live_work_and_other_policy_fields_and_is_repeatable(pool, monkeypatch):
    run_policy.update(pool.root, exclude=['node9'], timeout_retries=1, optimization_step_retries=2)
    before = run_policy.read(pool.root)
    pool.jobs.update({'old-a':1, 'old-b':1})
    assert submit(pool, 'blocked') is None
    real_run = run_capacity.subprocess.run
    def query(command, **kwargs):
        if command[0] == 'sinfo': return SimpleNamespace(stdout=NODES)
        return real_run(command, **kwargs)
    monkeypatch.setattr(run_capacity.subprocess, 'run', query)
    plan = run_capacity.fit(pool.root, partition='month-long-cpu', cores=8, memory_mb=98000)
    assert plan['max_jobs'] == 3 and run_policy.read(pool.root) == before
    result = run_capacity.fit(pool.root, partition='month-long-cpu', cores=8, memory_mb=98000, apply=True)
    assert result['applied'] and result['effective_max_jobs'] == 3
    after = run_policy.read(pool.root)
    for key, value in before.items():
        if key not in ('max_jobs', 'revision', 'updated_at'): assert after[key] == value
    assert submit(pool, 'blocked') is not None
    assert submit(pool, 'over-cap') is None
    assert set(pool.jobs) == {'old-a', 'old-b', '101'}
    again = run_capacity.fit(pool.root, partition='month-long-cpu', cores=8, memory_mb=98000, apply=True)
    assert not again['applied'] and run_policy.read(pool.root) == after


def test_policy_changed_during_sinfo_is_not_overwritten(tmp_path, monkeypatch):
    run_policy.update(tmp_path, max_jobs=2)
    def query(*a, **k):
        run_policy.update(tmp_path, exclude=['node2'])
        return SimpleNamespace(stdout=NODES)
    monkeypatch.setattr(run_capacity.subprocess, 'run', query)
    with pytest.raises(RuntimeError, match='changed'):
        run_capacity.fit(tmp_path, partition='cpu', cores=8, memory_mb=98000, apply=True)
    assert run_policy.read(tmp_path)['max_jobs'] == 2
    assert run_policy.read(tmp_path)['excluded_nodes'] == ['node2']


@pytest.mark.parametrize('error', ['outage', 'empty', 'no_eligible'])
def test_inventory_failure_preserves_policy(tmp_path, monkeypatch, error):
    run_policy.update(tmp_path, max_jobs=20)
    before = (tmp_path/run_policy.NAME).read_bytes()
    def query(*a, **k):
        if error == 'outage': raise run_capacity.subprocess.TimeoutExpired('sinfo', 30)
        return SimpleNamespace(stdout='' if error == 'empty' else 'node1|down|96|126000|up')
    monkeypatch.setattr(run_capacity.subprocess, 'run', query)
    with pytest.raises((ValueError, run_capacity.subprocess.TimeoutExpired)):
        run_capacity.fit(tmp_path, partition='cpu', cores=8, memory_mb=98000, apply=True)
    assert (tmp_path/run_policy.NAME).read_bytes() == before
