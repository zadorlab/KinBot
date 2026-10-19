"""Outage diagnostics are conservative and never advance the calculation."""
import json
from pathlib import Path

from kinbot import run_recap
from kinbot.anl import dispatch
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task
from tests.test_anl_dispatch import _mock_execution


def test_inventory_verifies_success_and_detects_corruption_without_mutation(tmp_path):
    spec = dispatch_spec()
    task = molpro_task('energy', 'rhf\n', geometry_from='initial')
    spec['tasks'] = [task]
    graph = dispatch.prepare(spec, tmp_path/'graph')
    _mock_execution(graph, 'energy')
    # Real native completion markers are required, in addition to hashes.
    work = graph/'tasks/energy'
    native = work/task['required_outputs'][0]
    native.write_text('Molpro calculation terminated\n')
    execution = json.loads((work/'execution.json').read_text())
    execution['artifacts'][native.name] = file_sha256(native)
    atomic_json(work/'execution.json', execution)
    before = {str(p): file_sha256(p) for p in graph.rglob('*') if p.is_file()}
    rows = run_recap.graph_recap(graph, {}, None, True)
    assert rows[0]['disposition'] == 'reusable'
    assert before == {str(p): file_sha256(p) for p in graph.rglob('*') if p.is_file()}
    native.write_text('truncated')
    assert run_recap.graph_recap(graph, {}, None, True)[0]['disposition'] == 'review_required'


def test_live_missing_and_unknown_scheduler_are_distinct(tmp_path):
    graph = dispatch.prepare(dispatch_spec(), tmp_path/'graph')
    state = dispatch._load(graph)[2]
    state['tasks']['l2_geometry'].update(status='submitted', job_id='123')
    atomic_json(graph/'state.json', state)
    live = {'123': {'job_id': '123', 'state': 'RUNNING', 'directory': str(graph/'tasks/l2_geometry')}}
    assert run_recap.graph_recap(graph, live, None, False)[0]['disposition'] == 'active'
    row = run_recap.graph_recap(graph, {}, None, False)[0]
    assert row['disposition'] == 'unknown' and 'disappearance' in row['reason']
    assert 'Scheduler unavailable' in run_recap.graph_recap(graph, {}, 'timeout', False)[0]['reason']
    live['123']['directory'] = str(tmp_path/'unrelated')
    assert 'different working directory' in run_recap.graph_recap(graph, live, None, False)[0]['reason']


def test_job_found_by_directory_even_if_driver_died_before_recording_id(tmp_path):
    graph = dispatch.prepare(dispatch_spec(), tmp_path/'graph')
    live = {'99': {'job_id': '99', 'state': 'RUNNING', 'directory': str(graph/'tasks/l2_geometry')}}
    row = run_recap.graph_recap(graph, live, None, False)[0]
    assert row['disposition'] == 'active' and row['job_id'] == '99'


def test_two_field_worker_markers_and_compact_diagnostics(tmp_path, monkeypatch):
    root = tmp_path/'run'
    case = root/'example'
    worker = case/'123'
    worker.mkdir(parents=True)
    atomic_json(root/'prepared.json', {'inputs': [str(case/'input.json')]})
    (worker/'a.exitcode').write_text('91 1\n')
    (worker/'b.exitcode').write_text('92 0\n')
    (worker/'a.log').write_text('Error termination: known native failure\n'*100)
    (worker/'b.log').write_text('Error: unrelated old attempt\n')
    monkeypatch.setattr('kinbot.run_cleanup._processes', lambda root: [])
    report = run_recap.recap(root, scheduler=([], None))
    case_report = report['cases'][0]
    assert len(case_report['worker_failures']) == len(case_report['worker_successes']) == 1
    assert case_report['worker_failures'][0]['native_diagnostics'][0]['lines'] == ['Error termination: known native failure']
    assert 'unrelated' not in run_recap.render(report)
    assert len(run_recap.render(report).splitlines()) < 15
    assert not report['restart_ready']
