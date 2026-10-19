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


def test_failed_geometry_is_primary_failure_and_dependents_are_not_ready(tmp_path):
    graph = dispatch.prepare(dispatch_spec(), tmp_path/'graph')
    state = dispatch._load(graph)[2]
    state['tasks']['l2_geometry'].update(status='unavailable', unavailable_reason='native failure')
    state['tasks']['l3_geometry'] = {'status': 'complete'}
    state['tasks']['f12_tz'] = {'status': 'staged'}
    atomic_json(graph/'state.json', state)
    rows = {r['task']: r for r in run_recap.graph_recap(graph, {}, None, False)}
    assert rows['l2_geometry']['disposition'] == 'failed'
    assert rows['l2_geometry']['geometry_invalid']
    assert rows['l3_geometry']['disposition'] == 'invalid_geometry'
    assert rows['f12_tz']['disposition'] == 'waiting'
    assert 'upstream geometry' in rows['f12_tz']['reason']


def test_uncertain_submission_with_output_still_requires_reconciliation(tmp_path):
    graph = dispatch.prepare(dispatch_spec(), tmp_path/'graph')
    _mock_execution(graph, 'l2_geometry', geometry='final.xyz')
    state = dispatch._load(graph)[2]
    state['tasks']['l2_geometry']['status'] = 'submitting'
    atomic_json(graph/'state.json', state)
    assert run_recap.graph_recap(graph, {}, None, True)[0]['disposition'] == 'unknown'
    live = {'7': {'job_id': '7', 'state': 'RUNNING', 'directory': str(graph/'tasks/l2_geometry')}}
    row = run_recap.graph_recap(graph, live, None, False)[0]
    assert row['disposition'] == 'active' and row['job_id'] == '7'
    live['8'] = {**live['7'], 'job_id': '8'}
    row = run_recap.graph_recap(graph, live, None, False)[0]
    assert row['disposition'] == 'unknown'
    assert row['live_job_ids'] == ['7', '8']


def test_evidence_retains_rank_failure_and_inputs_without_copying_wfu(tmp_path):
    graph = dispatch.prepare(dispatch_spec(), tmp_path/'graph')
    state = dispatch._load(graph)[2]
    state['tasks']['l2_geometry']['status'] = 'unavailable'
    atomic_json(graph/'state.json', state)
    work = graph/'tasks/l2_geometry'
    failure = {'returncode': 255, 'files': {'scratch/rank4.log':
               {'tail': '? Error\n? Error reading integral file\nInput/output error'}}}
    atomic_json(work/'l2_geometry_step_0001.failure.json', failure)
    atomic_json(work/'l2_geometry_step_0001.runtime.json', {'host': 'node6', 'scratch_directory': '/scratch/user/work'})
    (work/'l2_geometry_step_0001.inp').write_text('rhf;uccsd(t)\n')
    (work/'l2_geometry_step_0001.out').write_text('header\n' + 'normal progress\n'*20000 + '? Error reading integral file\n')
    (work/'integrals.wfu').write_text('do not copy this')
    (work/'fort.16').write_text('do not copy this either')
    before = {str(p): file_sha256(p) for p in graph.rglob('*') if p.is_file()}
    row = run_recap.graph_recap(graph, {}, None, False, evidence=True)[0]
    files = row['evidence']['files']
    assert files['l2_geometry_step_0001.failure.json']['json'] == failure
    assert files['l2_geometry_step_0001.inp']['head'] == 'rhf;uccsd(t)\n'
    assert files['l2_geometry_step_0001.out']['middle_omitted']
    assert len(files['l2_geometry_step_0001.out']['tail']) <= 65536
    assert not {'integrals.wfu', 'fort.16'} & files.keys()
    assert row['local_recovery_locations'][0]['host'] == 'node6'
    assert before == {str(p): file_sha256(p) for p in graph.rglob('*') if p.is_file()}


def test_render_prioritizes_native_failure_over_geometry_cascade_and_stale_status():
    tasks = [dict(graph='/run/example/123/states/a/composite', task=f'dependent{i}',
                  kind='external', disposition='invalid_geometry', reason='retry geometry') for i in range(40)]
    tasks += [dict(graph='/run/example/123/states/z/composite', task='l3_geometry',
                   kind='ase_optimize', disposition='failed', reason='force failed',
                   native_diagnostics=[
                       {'file': 'slurm.stderr', 'lines': ['Traceback (most recent call last):', 'RuntimeError: native failed']},
                       {'file': 'step.out', 'lines': ['? Error', '? Input/output error']}])]
    case = {'case': 'example', 'counts': {'failed': 1, 'invalid_geometry': 40},
            'live_processes': [], 'worker_failures': [], 'worker_successes': [], 'l2_models': [],
            'branches': [{'branches': {'cbh_anl': {'status': 'running'}}}], 'tasks': tasks, 'errors': []}
    report = {'time': 'now', 'root': '/run', 'scheduler_error': None, 'cases': [case],
              'live_jobs': [{'job_id': '31', 'state': 'RUNNING', 'nodes': 'node6',
                             'elapsed': '1:00', 'directory': '/run/example/123/tasks/gaussian_vpt2'}]}
    text = run_recap.render(report, limit=1)
    assert 'LIVE 31 RUNNING node6' in text
    assert 'no local controller found' in text
    assert 'requiring review/recompute: 40' in text
    assert 'z/l3_geometry: failed' in text and 'Input/output error' in text
    assert 'Traceback' not in text and 'additional issues' not in text
    assert run_recap.short_reason(dict(disposition='unknown', reason='Missing execution record',
                                      native_diagnostics=[{'file': 'step.out', 'lines': ['iterations'], 'tail_only': True}])) == 'Missing execution record'


def test_evidence_cli_writes_only_external_report(tmp_path, monkeypatch, capsys):
    root = tmp_path/'run'
    case = root/'example'
    graph = dispatch.prepare(dispatch_spec(), case/'123/anl_automatic/states/state_a/composite')
    atomic_json(root/'prepared.json', {'inputs': [str(case/'input.json')]})
    state = dispatch._load(graph)[2]
    state['tasks']['l2_geometry']['status'] = 'submitting'
    atomic_json(graph/'state.json', state)
    monkeypatch.setattr('kinbot.run_cleanup._processes', lambda root: [])
    monkeypatch.setattr('kinbot.run_monitor.queue', lambda: ([], None))
    before = {str(p): file_sha256(p) for p in root.rglob('*') if p.is_file()}
    report = tmp_path/'report.json'
    assert run_recap.main([str(root), '--json', str(report), '--evidence', '--limit', '1']) == 0
    saved = json.loads(report.read_text())
    row = saved['cases'][0]['tasks'][0]
    assert saved['evidence_included'] and row['disposition'] == 'unknown'
    assert row['evidence']['entry']['status'] == 'submitting'
    assert 'job.slurm' in row['evidence']['files']
    assert 'reconcile before any resubmission' in capsys.readouterr().out
    assert before == {str(p): file_sha256(p) for p in root.rglob('*') if p.is_file()}
