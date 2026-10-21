"""Fresh handoff interruption boundaries must never launch duplicate drivers."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot.execution import atomic_json
from kinbot.run_format import ensure_current_run

PACKAGE = Path(__file__).parents[1]/'examples/anl/tools'


@pytest.fixture
def fresh(monkeypatch, tmp_path):
    monkeypatch.syspath_prepend(str(PACKAGE))
    spec = importlib.util.spec_from_file_location('fresh_examples', PACKAGE/'fresh_run.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    config = json.loads((PACKAGE/'fresh_run.json').read_text())
    old = tmp_path/'r6'
    for case in ('ch3s', 'ch3ss'):
        ensure_current_run(old/case, create=True)
        (old/case/'native.out').write_text('expensive evidence')
    events = []
    monkeypatch.setattr(module, 'node_plan', lambda *a: dict(partition='month-long-cpu',
                        states={}, reserved_healthy_nodes=['n1'], max_jobs=2, excluded_nodes=['n1']))
    monkeypatch.setattr(module.run_cleanup, 'inspect', lambda d: {'jobs': [], 'processes': []})
    monkeypatch.setattr(module.run_cleanup, 'stop', lambda d, **kw: events.append(('stop', d.name)))
    monkeypatch.setattr(module.handoff, 'preflight', lambda r: events.append(('preflight', r.name)))
    def clean(archive, root, options, **kw):
        assert not old.exists() and Path(archive).is_dir()
        assert options['empty_home_wfu'] is True
        events.append(('cleanup', root.name))
        return {'status': 'complete'}
    monkeypatch.setattr(module.run_cache_cleanup, 'cleanup', clean)
    def launch(inputs, *, code):
        case = Path(inputs[0]).parent
        events.append(('launch', case.name))
        record = dict(pid=42, input=inputs[0], code=code)
        atomic_json(case/'launch.json', record)
        return [record]
    monkeypatch.setattr(module.handoff, 'start_inputs', launch)
    return SimpleNamespace(module=module, root=tmp_path/'r7', old=old, config=config, events=events)


def test_review_then_apply_retains_evidence_and_repeat_never_relaunches(fresh):
    f = fresh
    result = f.module.run(f.root, f.old, f.config)
    assert result['phase'] == 'reviewed' and f.old.is_dir()
    assert not any(x[0] in ('stop', 'launch') for x in f.events)
    result = f.module.run(f.root, f.old, f.config, apply=True)
    assert result['phase'] == 'launched' and not f.old.exists()
    assert (Path(result['archive'])/'ch3s/native.out').read_text() == 'expensive evidence'
    inputs = json.loads((f.root/'ch3s/input.json').read_text())
    assert inputs['vdW_detection'] == 0 and inputs['anl_augmented'] is True
    assert inputs['queue_live_policy']['max_jobs'] == 2
    assert inputs['queue_live_policy']['orbital_review_reserve'] == 3
    assert inputs['queue_live_control'] == {'hold_failures': False, 'hold_task_failures': True}
    assert not (f.root/'ch3s/.kinbot_resource_policy.json').exists()
    from kinbot import live_control
    assert live_control.find(f.root/'ch3s') == live_control.find(f.root/'ch3ss')
    assert live_control.find(f.root/'ch3s').parent == f.root
    calls = list(f.events)
    f.module.run(f.root, f.old, f.config, apply=True)
    assert calls == f.events


def test_failed_preflight_never_stops_old_run(fresh, monkeypatch):
    def fail(root):
        raise RuntimeError('missing executable')
    monkeypatch.setattr(fresh.module.handoff, 'preflight', fail)
    with pytest.raises(RuntimeError, match='missing executable'):
        fresh.module.run(fresh.root, fresh.old, fresh.config, apply=True)
    assert fresh.old.exists() and not fresh.events


def test_r7_scan_and_sampling_controls_reach_all_reviewed_surfaces(fresh):
    from kinbot.parameters import Parameters
    from kinbot.vrc_correction import asymptote_distance
    from kinbot.vrc_multistate import load_surface_model, state_parameters

    f = fresh
    f.module.handoff.prepare(f.root, overrides=f.config['input_overrides'])
    expected = ([i / 10 for i in range(20, 81, 2)]
                + [i / 2 for i in range(17, 25)] + list(range(13, 21)))
    assert len(expected) == len(set(expected)) == 47
    count = 0
    for case in f.module.handoff.CASES:
        source = PACKAGE.parent/f'{case}_anl1_f12'/f'{case}.json'
        original = Parameters(str(source), show_warnings=False).par
        par = Parameters(str(f.root/case/'input.json'), show_warnings=False).par
        assert par['rotdpy_flux_parameters'] == dict(
            original['rotdpy_flux_parameters'], flux_rel_err=5)
        assert par['vrc_tst_scan_points'] == expected
        assert par['rotdpy_dist'] == original['rotdpy_dist']
        surface, _ = load_surface_model(par)
        for ci, channel in enumerate(surface['channels']):
            for si, stage in enumerate(channel['bottlenecks']):
                previous = original['vrc_tst_surface_model']['channels'][ci]['bottlenecks'][si]
                assert stage['distances_angstrom'] == previous['distances_angstrom']
                assert stage['states'] == previous['states']
                for state in stage['states']:
                    leaf = state_parameters(par, stage, state, [2, 3])
                    assert leaf['vrc_tst_scan_points'] == expected
                    assert leaf['rotdpy_flux_parameters']['flux_rel_err'] == 5
                    assert leaf['rotdpy_dist'] == previous['distances_angstrom']
                    assert asymptote_distance(leaf) == asymptote_distance(original) == 30.
                    count += 1
    assert count == 4


def test_interrupted_after_archive_move_recovers_without_losing_evidence(fresh, monkeypatch):
    original = fresh.module.retire
    def interrupt(*args):
        original(*args)
        raise RuntimeError('connection lost')
    monkeypatch.setattr(fresh.module, 'retire', interrupt)
    with pytest.raises(RuntimeError, match='connection lost'):
        fresh.module.run(fresh.root, fresh.old, fresh.config, apply=True)
    monkeypatch.setattr(fresh.module, 'retire', original)
    result = fresh.module.run(fresh.root, fresh.old, fresh.config, apply=True)
    assert result['phase'] == 'launched'
    assert fresh.events.count(('stop', 'ch3s')) == 1


def test_uncertain_launch_blocks_duplicate_and_keeps_other_case_unstarted(fresh, monkeypatch):
    def uncertain(*args, **kwargs):
        raise RuntimeError('process launch outcome unknown')
    monkeypatch.setattr(fresh.module.handoff, 'start_inputs', uncertain)
    with pytest.raises(RuntimeError, match='launch outcome'):
        fresh.module.run(fresh.root, fresh.old, fresh.config, apply=True)
    with pytest.raises(RuntimeError, match='launch outcome is uncertain'):
        fresh.module.run(fresh.root, fresh.old, fresh.config, apply=True)
    journal = json.loads((fresh.root/'fresh_run.json').read_text())
    assert set(journal['launches']) == {'ch3s'}


def test_partial_launch_receipt_is_adopted_and_only_second_case_starts(fresh, monkeypatch):
    original = fresh.module.handoff.start_inputs
    def interrupt(*args, **kwargs):
        original(*args, **kwargs)
        raise RuntimeError('receipt saved then interrupted')
    monkeypatch.setattr(fresh.module.handoff, 'start_inputs', interrupt)
    with pytest.raises(RuntimeError):
        fresh.module.run(fresh.root, fresh.old, fresh.config, apply=True)
    monkeypatch.setattr(fresh.module.handoff, 'start_inputs', original)
    fresh.module.run(fresh.root, fresh.old, fresh.config, apply=True)
    assert fresh.events.count(('launch', 'ch3s')) == 1
    assert fresh.events.count(('launch', 'ch3ss')) == 1


def test_cache_guard_blocks_launch_and_retry_does_not_repeat_retirement(fresh, monkeypatch):
    f = fresh
    clean = f.module.run_cache_cleanup.cleanup
    def busy(*args, **kwargs):
        raise RuntimeError('ALL account allocations')
    monkeypatch.setattr(f.module.run_cache_cleanup, 'cleanup', busy)
    with pytest.raises(RuntimeError, match='ALL account allocations'):
        f.module.run(f.root, f.old, f.config, apply=True)
    assert not any(x[0] == 'launch' for x in f.events)
    assert json.loads((f.root/'fresh_run.json').read_text())['phase'] == 'retired'
    monkeypatch.setattr(f.module.run_cache_cleanup, 'cleanup', clean)
    f.module.run(f.root, f.old, f.config, apply=True)
    assert f.events.count(('stop', 'ch3s')) == 1
    assert f.events.index(('cleanup', 'r7')) < f.events.index(('launch', 'ch3s'))


def test_node_inventory_counts_nodes_once_and_never_admits_unhealthy(monkeypatch):
    monkeypatch.syspath_prepend(str(PACKAGE))
    import fresh_run
    text = 'n1|idle\nn1|idle\nn2|allocated\nn3|mixed\nn4|idle\nn5|idle\nn6|down\nn7|drained\nn8|unknown\nn9|idle*\n'
    monkeypatch.setattr(fresh_run.subprocess, 'run', lambda *a, **k: SimpleNamespace(stdout=text))
    plan = fresh_run.node_plan('long', 2, preferred=['n2'])
    assert plan['max_jobs'] == 3
    assert plan['reserved_healthy_nodes'] == ['n2', 'n1']
    assert set(plan['excluded_nodes']) == {'n1', 'n2', 'n6', 'n7', 'n8', 'n9'}
