"""The missing-passwd regression is repaired once, without duplicate QC work."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from kinbot.anl import account_recovery, dispatch
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task
from tests.test_anl_dispatch import _mock_execution


def failed(graph, ident, *, error="'getpwuid(): uid not found: 155382'", trace=True):
    state = dispatch._load(graph)[2]
    atomic_json(graph/'tasks'/ident/'execution.json', {
        'schema': 1, 'task_id': ident, 'geometry_sha256': state['tasks'][ident]['geometry_sha256'],
        'status': 'failed', 'failure_phase': 'validation', 'error': error,
        'traceback': 'File "/kinbot/qc_scratch.py", line 44\n pwd.getpwuid(os.getuid())\nKeyError' if trace else ''})


@pytest.fixture
def graph(tmp_path, monkeypatch):
    spec = dispatch_spec()
    spec['tasks'] = [molpro_task(n, 'rhf\n', geometry_from='initial')
                     for n in ('bad', 'live', 'done', 'child')]
    spec['tasks'][-1]['depends_on'] = ['bad']
    root = dispatch.prepare(spec, tmp_path/'graph')
    failed(root, 'bad')
    _mock_execution(root, 'done')
    state = dispatch._load(root)[2]
    state['tasks']['bad'].update(status='failed', job_id='42')
    state['tasks']['live'].update(status='submitted', job_id='43')
    state['tasks']['child'] = {'status': 'blocked', 'blocked_by': ['bad']}
    atomic_json(root/'state.json', state)
    monkeypatch.setattr(dispatch, '_job_active', lambda job: str(job) == '43')
    return root


def test_legacy_failure_archived_and_retried_without_touching_live_or_success(graph):
    spec_before = (graph/'workflow.json').read_bytes()
    untouched = {p: file_sha256(p) for name in ('live', 'done')
                 for p in (graph/'tasks'/name).iterdir() if p.is_file()}
    state = dispatch.advance(graph)
    assert state['tasks']['bad']['status'] == 'staged'
    assert state['tasks']['bad']['attempt'] == 2
    assert state['tasks']['live']['job_id'] == '43'
    assert state['tasks']['live']['status'] == 'submitted'
    assert state['tasks']['done']['status'] == 'complete'
    assert state['tasks'].get('child', {}).get('status') != 'blocked'
    archive = Path(state['tasks']['bad']['previous_attempt'])
    assert json.loads((archive/'execution.json').read_text())['error'].startswith("'getpwuid")
    assert (graph/'workflow.json').read_bytes() == spec_before
    assert untouched == {p: file_sha256(p) for p in untouched}
    assert not (graph/'tasks/bad/execution.json').exists()
    assert dispatch.advance(graph) == state


def test_same_bug_cannot_trigger_unbounded_retries(graph):
    dispatch.advance(graph)
    failed(graph, 'bad')
    state = dispatch.advance(graph)
    assert state['tasks']['bad']['status'] == 'failed'
    assert state['tasks']['bad']['attempt'] == 2


@pytest.mark.parametrize('change', ['other_error', 'missing_trace', 'wrong_provenance'])
def test_other_failures_not_retried(graph, change):
    if change == 'other_error':
        failed(graph, 'bad', error='SCF did not converge')
    elif change == 'missing_trace':
        failed(graph, 'bad', trace=False)
    else:
        p = graph/'tasks/bad/execution.json'
        v = json.loads(p.read_text()); v['geometry_sha256'] = 'wrong'; atomic_json(p, v)
        with pytest.raises(RuntimeError, match='execution record does not match'):
            dispatch.advance(graph)
        return
    assert dispatch.advance(graph)['tasks']['bad']['status'] == 'failed'
    assert not (graph/'account_lookup_recovery.json').exists()


def test_live_failure_record_is_not_reset_and_is_recovered_after_queue_exit(graph, monkeypatch):
    state = dispatch._load(graph)[2]
    state['tasks']['bad'].update(status='submitted')
    atomic_json(graph/'state.json', state)
    monkeypatch.setattr(dispatch, '_job_active', lambda job: True)
    assert dispatch.advance(graph)['tasks']['bad']['status'] == 'submitted'
    assert not (graph/'account_lookup_recovery.json').exists()
    monkeypatch.setattr(dispatch, '_job_active', lambda job: str(job) == '43')
    assert dispatch.advance(graph)['tasks']['bad']['status'] == 'staged'


def test_interrupted_archive_resumes_before_loading_changed_ledger(graph, monkeypatch):
    from kinbot import run_recover
    original = run_recover.atomic_json
    stopped = False
    def interrupted(path, data):
        nonlocal stopped
        if Path(path) == graph/'state.json' and not stopped:
            stopped = True
            raise OSError('interrupted after archiving task')
        return original(path, data)
    monkeypatch.setattr(run_recover, 'atomic_json', interrupted)
    with pytest.raises(OSError, match='interrupted'):
        dispatch.advance(graph)
    assert not (graph/'tasks/bad').exists()
    monkeypatch.setattr(run_recover, 'atomic_json', original)
    state = dispatch.advance(graph)
    assert state['tasks']['bad']['status'] == 'staged' and state['tasks']['bad']['attempt'] == 2
    assert state['tasks']['live']['job_id'] == '43'
    assert len(list((graph/'attempts').glob('account_lookup_*'))) == 1
