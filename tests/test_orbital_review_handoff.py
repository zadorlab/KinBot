"""Refresh/reconciliation/panel ordering; native/site operations are simulated."""
import fcntl
import json

import pytest

from kinbot import orbital_review_handoff as handoff
from kinbot.execution import atomic_json
from test_remaining_character_panel import tree, plan, PLAN


def test_full_orchestration_preflights_before_refresh_and_repeats_without_restaging(tree, monkeypatch):
    root, output = tree
    source = root/'ch3ss/input.json'
    native = root/'existing-native.out'; native.write_bytes(b'accepted QC')
    atomic_json(root/'existing-state.json', {'job_id': '901', 'status': 'submitted'})
    before = {p:p.read_bytes() for p in root.rglob('*') if p.is_file()}
    events = []
    batch = handoff.review.batch
    def prepare(*args, **kwargs):
        events.append('preflight' if kwargs.get('preflight') else 'prepare')
        return batch(*args, **kwargs)
    monkeypatch.setattr(handoff.review, 'batch', prepare)
    def refresh(production, **options):
        assert events[-1] == 'preflight'
        assert options == dict(apply=True, wait_idle=600.)
        with production.with_suffix('.salvage.lock').open('a') as probe:
            with pytest.raises(BlockingIOError):
                fcntl.flock(probe, fcntl.LOCK_EX | fcntl.LOCK_NB)
        events.append('refresh')
        return dict(phase='launched', code='pinned', jobs=[dict(job_id='901')])
    monkeypatch.setattr(handoff.run_refresh, 'refresh_when_idle', refresh)
    readiness = iter([['owner starting'], [], []])
    monkeypatch.setattr(handoff, '_ready', lambda inputs: next(readiness))
    monkeypatch.setattr(handoff.time, 'sleep', lambda seconds: events.append('wait'))
    def submit(path, **options):
        assert path == source and options['reconcile_models']
        events.append('apply-models' if options.get('apply') else 'preview-models')
        return dict(status='queued', drivers_started=0)
    monkeypatch.setattr(handoff.input_refresh, 'submit', submit)
    def run(panel, **options):
        assert panel == output and options == dict(report_failures=True)
        assert len(handoff.review._saved(panel)['systems']) == 4
        events.append('run')
        return dict(status='complete_review_required')
    monkeypatch.setattr(handoff.review, 'run', run)
    result = handoff.perform(root, PLAN/'screens.json', output, [source],
                             apply=True, expected_max_jobs=23)
    assert events == ['preflight','refresh','wait','preview-models','apply-models','prepare','run']
    assert not result['production_models_accepted']
    staged = {p:p.read_bytes() for p in output.rglob('*') if p.is_file()}
    events.clear()
    handoff.perform(root, PLAN/'screens.json', output, [source], apply=True, expected_max_jobs=23)
    assert events == ['preflight','refresh','preview-models','apply-models','prepare','run']
    assert staged == {p:p.read_bytes() for p in output.rglob('*') if p.is_file()}
    assert all(p.read_bytes() == content for p, content in before.items())


def test_bad_later_recipe_refuses_before_refresh(tree, tmp_path, monkeypatch):
    root, output = tree
    invalid = tmp_path/'invalid'; invalid.mkdir()
    for path in PLAN.glob('*.json'): (invalid/path.name).write_bytes(path.read_bytes())
    recipe = json.loads((invalid/'hcs_h2s_exit.json').read_text())
    recipe['geometry_source']['record_name'] = 'nonexistent_endpoint'
    atomic_json(invalid/'hcs_h2s_exit.json', recipe)
    monkeypatch.setattr(handoff.run_refresh, 'refresh_when_idle',
                        lambda *a,**k: pytest.fail('No refresh after failed panel preflight'))
    with pytest.raises(ValueError):
        handoff.perform(root, invalid/'screens.json', output, [root/'ch3ss/input.json'], apply=True)
    assert not output.exists()


def test_owner_timeout_or_concurrent_input_refusal_does_not_launch_review(tree, monkeypatch):
    root, output = tree
    monkeypatch.setattr(handoff.run_refresh, 'refresh_when_idle', lambda *a,**k: dict(phase='launched'))
    monkeypatch.setattr(handoff, '_ready', lambda inputs: ['PES owner not ready'])
    monkeypatch.setattr(handoff.review, 'run', lambda *a,**k: pytest.fail('No review launch'))
    with pytest.raises(RuntimeError, match='readiness wait expired'):
        handoff.perform(root, PLAN/'screens.json', output, [root/'ch3ss/input.json'],
                        apply=True, wait_owner=0)
    assert not output.exists()
    monkeypatch.setattr(handoff, '_ready', lambda inputs: [])
    def changed(*a, **k): raise ValueError('unaccepted scientific input edit')
    monkeypatch.setattr(handoff.input_refresh, 'submit', changed)
    with pytest.raises(ValueError, match='unaccepted scientific'):
        handoff.perform(root, PLAN/'screens.json', output, [root/'ch3ss/input.json'], apply=True)
    assert not output.exists()


def test_preview_and_duplicate_owner_do_not_replace_controllers(tree, monkeypatch):
    root, output = tree
    monkeypatch.setattr(handoff.run_refresh, 'refresh_when_idle', lambda *a,**k: pytest.fail('No refresh'))
    before = {p:p.read_bytes() for p in root.rglob('*') if p.is_file()}
    assert handoff.perform(root, PLAN/'screens.json', output, [root/'ch3ss/input.json'])['status']=='reviewed'
    assert before == {p:p.read_bytes() for p in root.rglob('*') if p.is_file()}
    with (root/'.kinbot_orbital_handoff.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(RuntimeError, match='Another orbital handoff'):
            handoff.perform(root, PLAN/'screens.json', output, [root/'ch3ss/input.json'], apply=True)
    with pytest.raises(ValueError, match='job limit'):
        handoff.perform(root, PLAN/'screens.json', output, [root/'ch3ss/input.json'],
                        expected_max_jobs=24, apply=True)


def test_readiness_distinguishes_startup_pending_blocked_and_prior_generation(tmp_path, monkeypatch):
    source = tmp_path/'input.json'; source.write_text('{}')
    monkeypatch.setattr(handoff.input_refresh, 'owner_directory', lambda p: None)
    assert handoff._ready([source])
    monkeypatch.setattr(handoff.input_refresh, 'owner_directory', lambda p: tmp_path)
    atomic_json(tmp_path/handoff.pes_phase.NAME,
                dict(token='new', phase='exploring', workers={'worker':{'status':'vrc'}}))
    monkeypatch.setattr(handoff.pes_phase, 'active_worker', lambda p: False)
    assert 'waiting for workers' in handoff._ready([source])[0]
    monkeypatch.setattr(handoff.pes_phase, 'active_worker', lambda p: True)
    assert handoff._ready([source]) == []
    atomic_json(tmp_path/handoff.input_refresh.JOURNAL, dict(token='new'))
    monkeypatch.setattr(handoff.input_refresh, 'receipt', lambda p: dict(status='queued'))
    assert 'revision queued' in handoff._ready([source])[0]
    monkeypatch.setattr(handoff.input_refresh, 'receipt', lambda p: dict(status='blocked'))
    with pytest.raises(RuntimeError, match='blocked'): handoff._ready([source])
    atomic_json(tmp_path/handoff.input_refresh.JOURNAL, dict(token='old'))
    assert handoff._ready([source]) == []
