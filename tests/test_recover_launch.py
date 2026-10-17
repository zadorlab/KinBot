"""Cluster-interrupted handoffs must not duplicate drivers or lose old evidence."""
import importlib.util
import json
from pathlib import Path

import pytest


@pytest.fixture
def recovery(monkeypatch):
    package = Path(__file__).parents[1]/'examples/anl/tools'
    monkeypatch.syspath_prepend(str(package))
    spec = importlib.util.spec_from_file_location('recover_launch_test', package/'recover_launch.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module.handoff, 'preflight', lambda root: {'passed': True})
    monkeypatch.setattr(module.handoff, 'archive_legacy_example_layout', lambda root: [])
    monkeypatch.setattr(module, 'old_handoffs', lambda root: [])
    return module


def test_preparation_manifest_excludes_new_runtime_inputs_and_wfu(recovery, tmp_path):
    root = recovery.handoff.prepare(tmp_path/'r7')
    before = (root/'prepared.json').read_bytes()
    runtime = root/'ch3s/123/anl_automatic/tasks'
    runtime.mkdir(parents=True)
    (runtime/'input.json').write_text('{"runtime":true}')
    (runtime/'step.wfu').write_bytes(b'calculated orbitals')
    recovery.handoff.prepare(root)
    assert (root/'prepared.json').read_bytes() == before


def test_partial_driver_launch_retries_only_missing_case(recovery, monkeypatch, tmp_path):
    root = tmp_path/'r7'
    old = recovery.handoff.prepare(tmp_path/'r6')
    (old/'ch3s/native.out').write_text('evidence')
    live = {}
    monkeypatch.setattr(recovery.run_cleanup, 'inspect', lambda p, **kw:
                        {'jobs': [], 'processes': live.get(str(p), [])})
    stopped = []
    monkeypatch.setattr(recovery.run_cleanup, 'stop', lambda p, **kw:
                        stopped.append(p) or {'status': 'stopped'})
    monkeypatch.setattr(recovery, 'argv_for', lambda p: p['argv'])
    spawned = []

    def spawn(inputs):
        result = []
        for filename in inputs:
            directory = Path(filename).parent
            if len(spawned) == 1 and directory.name == 'ch3ss':
                spawned.append('simulated interruption')
                raise OSError('connection lost between launches')
            process = {'pid': 100 + len(spawned), 'argv': ['python', '-m', 'kinbot.pes', filename]}
            live[str(directory)] = [process]
            (directory/'launch.json').write_text(json.dumps({'pid': process['pid']}))
            spawned.append(directory.name)
            result.append({'pid': process['pid']})
        return result

    monkeypatch.setattr(recovery.handoff, 'start_inputs', spawn)
    with pytest.raises(OSError, match='between launches'):
        recovery.recover(root, [old])
    assert stopped == [old/'ch3s', old/'ch3ss']
    assert not old.exists()
    archives = list(tmp_path.glob('r6.stopped.*'))
    assert len(archives) == 1 and (archives[0]/'ch3s/native.out').read_text() == 'evidence'
    second = recovery.recover(root, [old])
    assert len(second['launched']) == 1
    assert second['retained'][0]['state'] == 'running'
    third = recovery.recover(root, [old])
    assert not third['launched']
    assert [item['state'] for item in third['retained']] == ['running', 'running']
    assert spawned == ['ch3s', 'simulated interruption', 'ch3ss']
    assert len(list(tmp_path.glob('r6.stopped.*'))) == 1


def test_scheduler_outage_prevents_retirement_and_launch(recovery, monkeypatch, tmp_path):
    old = recovery.handoff.prepare(tmp_path/'r6')
    def unavailable(*args, **kwargs):
        raise RuntimeError('Slurm unavailable')
    monkeypatch.setattr(recovery.run_cleanup, 'inspect', unavailable)
    monkeypatch.setattr(recovery.run_cleanup, 'stop', lambda *a, **k: pytest.fail('must not stop'))
    monkeypatch.setattr(recovery.handoff, 'start_inputs', lambda *a: pytest.fail('must not launch'))
    with pytest.raises(RuntimeError, match='Slurm unavailable'):
        recovery.recover(tmp_path/'r7', [old])
    assert old.exists()
    status = json.loads((tmp_path/'r7/restart_status.json').read_text())
    assert 'Slurm unavailable' in status['error']


def test_partial_retirement_can_be_retried(recovery, monkeypatch, tmp_path):
    old = recovery.handoff.prepare(tmp_path/'r6')
    root = tmp_path/'r7'
    monkeypatch.setattr(recovery.run_cleanup, 'inspect', lambda *a, **k: {'jobs': [], 'processes': []})
    stopped = []
    def stop(path, **kw):
        stopped.append(path.name)
        if len(stopped) == 2:
            raise RuntimeError('scheduler temporarily down')
        return {'status': 'stopped'}
    monkeypatch.setattr(recovery.run_cleanup, 'stop', stop)
    monkeypatch.setattr(recovery.handoff, 'start_inputs', lambda paths: list(paths))
    with pytest.raises(RuntimeError, match='temporarily down'):
        recovery.recover(root, [old])
    assert old.exists()
    result = recovery.recover(root, [old])
    assert len(result['launched']) == 2
    assert not old.exists()


def test_stale_new_run_is_not_blindly_restarted(recovery, monkeypatch, tmp_path):
    root = recovery.handoff.prepare(tmp_path/'r7')
    (root/'ch3s/launch.json').write_text('{"pid":123}')
    monkeypatch.setattr(recovery.run_cleanup, 'inspect', lambda *a, **k: {'jobs': [], 'processes': []})
    monkeypatch.setattr(recovery.handoff, 'start_inputs', lambda *a: pytest.fail('must not launch'))
    with pytest.raises(RuntimeError, match='previously started but no live'):
        recovery.recover(root, [])


def test_completed_new_cases_are_retained(recovery, monkeypatch, tmp_path):
    root = recovery.handoff.prepare(tmp_path/'r7')
    for case in recovery.handoff.CASES:
        (root/case/'pes_status.json').write_text('{"status":"complete"}')
    monkeypatch.setattr(recovery.run_cleanup, 'inspect', lambda *a, **k: {'jobs': [], 'processes': []})
    monkeypatch.setattr(recovery.handoff, 'start_inputs', lambda paths: list(paths))
    result = recovery.recover(root, [])
    assert not result['launched']
    assert all(item['state'] == 'complete' for item in result['retained'])


def test_active_handoff_lock_prevents_second_preparation(recovery, monkeypatch, tmp_path):
    import fcntl
    root = tmp_path/'r7'
    root.mkdir()
    monkeypatch.setattr(recovery.handoff, 'prepare', lambda *a: pytest.fail('must not prepare'))
    with (root/'.handoff.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert recovery.recover(root, []) == {'status': 'handoff_already_active'}


def test_launcher_selection_uses_exact_script_root_and_action(monkeypatch, tmp_path):
    # Load without the fixture's old_handoffs stub.
    package = Path(__file__).parents[1]/'examples/anl/tools'
    monkeypatch.syspath_prepend(str(package))
    import recover_launch as module
    old = tmp_path/'r6'
    script = str(package/'run_examples.py')
    commands = [
        ['python', '-u', script, 'launch', str(old)],
        ['python', script, 'launch', str(tmp_path/'other')],
        ['python', script, 'monitor', str(old)],
        ['python', '/unrelated/run_examples.py', 'launch', str(old)],
    ]
    processes = [{'pid': i, 'starttime': 'original', 'argv': args} for i, args in enumerate(commands)]
    monkeypatch.setattr(module.run_cleanup, '_processes', lambda *a: processes)
    monkeypatch.setattr(module, 'argv_for', lambda p: p['argv'])
    assert module.old_handoffs(old) == [processes[0]]


def test_retirement_includes_verified_launcher_identity(recovery, monkeypatch, tmp_path):
    old = recovery.handoff.prepare(tmp_path/'r6')
    root = recovery.handoff.prepare(tmp_path/'r7')
    identity = {'pid': 123, 'starttime': 'original'}
    monkeypatch.setattr(recovery.run_cleanup, 'inspect', lambda *a, **k: {'jobs': [], 'processes': []})
    calls = []
    monkeypatch.setattr(recovery.run_cleanup, 'stop', lambda p, **kw: calls.append(kw) or {})
    recovery.handoff.retire_production(root, old, known=[identity])
    assert calls[0]['known'] == [identity]
    assert 'known' not in calls[1]


def test_retire_partial_prepare_preserves_unfinished_artifacts(recovery, monkeypatch, tmp_path):
    from kinbot.run_format import ensure_current_run
    old, root = tmp_path/'r6', tmp_path/'r7'
    ensure_current_run(old/'ch3s', create=True)
    root.mkdir()
    (old/'ch3s/orbital.wfu.partial').write_bytes(b'partial input')
    monkeypatch.setattr(recovery.run_cleanup, 'inspect', lambda *a, **k: {'jobs': [], 'processes': []})
    monkeypatch.setattr(recovery.run_cleanup, 'stop', lambda *a, **k: {})
    result = recovery.handoff.retire_production(root, old, allow_partial=True)
    assert (Path(result['archive'])/'ch3s/orbital.wfu.partial').read_bytes() == b'partial input'


def test_remote_controller_host_is_not_assumed_inactive(recovery, monkeypatch, tmp_path):
    old = recovery.handoff.prepare(tmp_path/'r6')
    root = recovery.handoff.prepare(tmp_path/'r7')
    (old/'ch3s/pes_status.json').write_text('{"host":"another-host","status":"running"}')
    monkeypatch.setattr(recovery.run_cleanup, 'inspect', lambda *a, **k: pytest.fail('wrong host'))
    with pytest.raises(RuntimeError, match='retire from controller host'):
        recovery.handoff.retire_production(root, old)
    assert old.exists()
