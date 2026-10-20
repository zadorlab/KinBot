import hashlib
import json
from pathlib import Path

import pytest

from kinbot import run_progress


def put(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) if isinstance(value, dict) else value)


@pytest.fixture
def campaign(tmp_path):
    root = tmp_path/'campaign'
    worker = root/'case'/'123'
    model = {'channels': [{'bond': [1, 2], 'bottlenecks': [dict(
        name='single', states=[{'multiplicity': 2, 'root': 1},
                               {'multiplicity': 2, 'root': 2}],
        parameters={}, distances_angstrom=[3., 4., 5.])]}]}
    put(root/'case/input.json', {'vrc_tst_surface_model': model})
    put(worker/'123.json', {})
    put(worker/'kinbot.log', '    vrc_tst_scan_points: [2.5, 2.7, 2.9]\n'
                             '    vrc_tst_scan_fmax: 0.0005\n')
    graph = worker/'anl_automatic/states/molecule/composite'
    put(graph/'workflow.json', {'molecule': {'smiles': 'C'}, 'tasks': [
        dict(id='l3_geometry', kind='ase_optimize', optimizer={'fmax': .0005}),
        dict(id='anl1_6z', kind='external', depends_on=['l3_geometry'])]})
    put(graph/'state.json', {'tasks': {'l3_geometry': {'status': 'submitted'}}})
    put(graph/'tasks/l3_geometry/optimization.log',
        'Sella 0 10:00:00 -10.0 0.1\nSella 1 10:20:00 -10.1 0.02\n')
    put(root/'case/pes_summary.txt', 'Total 2\nRunning 1\nFinished 0\n')
    return root, worker, graph


def snapshot(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob('*') if p.is_file()}


def test_report_reads_only_and_preserves_unknown_durations(campaign):
    root, worker, graph = campaign
    before = snapshot(root)
    report = run_progress.collect(root, scheduler=([dict(job_id='12', state='RUNNING',
        elapsed='1:20:00', directory=str(graph/'tasks/l3_geometry'))], None))
    assert snapshot(root) == before
    w = report['cases'][0]['workers'][0]
    assert [p['relaxed_points'] for p in w['surface_plan']] == [3, 3]
    task, downstream = w['anl'][0]['tasks']
    assert task['last_fmax'] == .02
    assert task['live']['elapsed'] == '1:20:00'
    assert downstream['depends_on'] == ['l3_geometry']
    assert downstream['status'] == 'not_materialized'
    assert 'eta' not in task
    assert 'planned relaxed points=6' in run_progress.render(report)


def test_no_default_grid_invented_and_stage_override(campaign):
    root, worker, graph = campaign
    (worker/'kinbot.log').write_text('')
    report = run_progress.collect(root, scheduler=([], None))
    assert report['cases'][0]['workers'][0]['surface_plan'][0]['relaxed_points'] is None
    params = {'vrc_tst_scan_points': [2., 3.], 'vrc_tst_surface_model': {
        'channels': [{'bottlenecks': [{'parameters': {'vrc_tst_scan_points': [5.]},
            'states': [{'multiplicity': 2, 'root': 1}]}]}]}}
    assert run_progress.surface_plan(params)[0]['relaxed_points'] == 1


def test_current_graph_pointer_only_and_no_retry_archives(campaign):
    root, worker, graph = campaign
    for name in ('dispatch', 'dispatch_new'):
        put(worker/f'vrctst/molpro/{name}/workflow.json', {'tasks': [{'id': name}]})
        put(worker/f'vrctst/molpro/{name}/state.json', {'tasks': {name: {'status': 'complete'}}})
    put(worker/'vrctst/molpro/current_dispatch.json',
        {'run_dir': str(worker/'vrctst/molpro/dispatch_new')})
    put(worker/'anl_automatic/states/old/archive/composite/workflow.json', {'tasks': [{'id': 'archived'}]})
    report = run_progress.collect(root, scheduler=([], None))
    w = report['cases'][0]['workers'][0]
    assert len(w['anl']) == 1
    assert w['correction_graph']['tasks'][0]['task'] == 'dispatch_new'


def test_bounded_native_reads_and_scheduler_failure(campaign):
    root, worker, graph = campaign
    native = graph/'tasks/l3_geometry/step.live.log'
    put(native, 'not the tail\n' + 'x' * (run_progress.LIMIT * 2) + '\nREAL TIME * 20 SEC\n')
    report = run_progress.collect(root, scheduler=([], 'scheduler unavailable'))
    assert report['scheduler_error'] == 'scheduler unavailable'
    record = report['cases'][0]['workers'][0]['anl'][0]['tasks'][0]
    assert 'not the tail' not in '\n'.join(record['native']['tail'])
    assert record['native']['timing_lines'] == ['REAL TIME * 20 SEC']


def test_output_cannot_modify_run_or_existing_file(campaign, monkeypatch, tmp_path):
    root, _, _ = campaign
    monkeypatch.setattr(run_progress, 'collect', lambda *a, **kw: pytest.fail('Must validate first'))
    with pytest.raises(SystemExit):
        run_progress.main([str(root), '--output', str(root/'report.json')])
    output = tmp_path/'already_exists.json'
    output.write_text('keep')
    with pytest.raises(SystemExit):
        run_progress.main([str(root), '--output', str(output)])
    assert output.read_text() == 'keep'


def test_read_errors_are_visible(campaign):
    root, worker, graph = campaign
    (graph/'state.json').write_text('{')
    report = run_progress.collect(root, scheduler=([], None))
    assert report['read_errors']
    assert 'read errors=1' in run_progress.render(report)
