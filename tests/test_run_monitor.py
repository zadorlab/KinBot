import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import run_monitor as monitor
ORIGINAL_QUEUE = monitor.queue


def put(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) if isinstance(value, dict) else value)


@pytest.fixture
def production(tmp_path, monkeypatch):
    root = tmp_path/'runs'
    case = root/'example'
    put(root/'prepared.json', {'inputs': [str(case/'input.json')]})
    put(case/'input.json', {'queue_job_limit': 10, 'slurm_exclude_nodes': ['node16']})
    monkeypatch.setattr(monitor, 'queue', lambda: ([], None))
    from kinbot import run_cleanup
    monkeypatch.setattr(run_cleanup, '_processes', lambda *a: [])
    return root, case


def test_initial_failure_visible_even_without_accounting_or_driver_output(production):
    root, case = production
    put(case/'driver.log', 'Copyright ...\n')
    put(case/'pes_summary.txt', 'Total 1\nRunning 0\nFinished 1\n\nFinished:\n  well\n')
    put(case/'well/kinbot.log', '16-Oct-26 18:35:07-INFO: Starting high level optimization\n'
        '16-Oct-26 18:35:17-ERROR: Error with high level optimization of initial structure.\n')
    put(case/'well/well_high.log', 'QPErr --- A syntax error was detected\nError termination via Lnk1e\n')
    put(case/'well/well_high.exitcode', '42 1\n')
    put(case/'well/perm/well_high.err', 'Error: segmentation violation\nTraceback (most recent call last):\n'
        'RuntimeError: Gaussian exited with error code 1\n')
    before = {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}
    report = monitor.snapshot(root)
    text = monitor.render(report)
    assert 'QPErr' in text and 'FAILED QC EXIT' in text
    assert 'Starting high level optimization' in text
    assert 'Finished means exited, not successful' in text
    assert 'INACTIVE:' in text and 'completion verdict' in text
    assert 'Copyright' not in text
    assert before == {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}


def test_unknown_scheduler_not_reported_as_empty(production, monkeypatch):
    root, _ = production
    monkeypatch.setattr(monitor, 'queue', lambda: ([], 'Slurm unavailable'))
    text = monitor.render(monitor.snapshot(root))
    assert 'Slurm: UNKNOWN' in text
    assert 'INACTIVE:' not in text


def test_active_jobs_scoped_and_resources_visible(production, monkeypatch):
    root, case = production
    raw = ('42|RUNNING|node20|00:30|31-00:00:00|8|10G|node16|well_high|node20|'
           f'{case}/well\n43|RUNNING|node9|00:30|1:00:00|8|10G||other|node9|/unrelated\n')
    # Test the queue parser separately, then inject its result into snapshot.
    monkeypatch.setattr(monitor, 'queue', ORIGINAL_QUEUE)
    monkeypatch.setattr(monitor.subprocess, 'run', lambda *a, **k: SimpleNamespace(returncode=0, stdout=raw))
    jobs, error = monitor.queue()
    assert error is None and jobs[0]['cpus'] == '8'
    monkeypatch.setattr(monitor, 'queue', lambda: (jobs, None))
    from kinbot import run_cleanup
    monkeypatch.setattr(run_cleanup, '_processes', lambda *a: [])
    text = monitor.render(monitor.snapshot(root))
    assert 'well_high' in text and 'Slurm exclusions: node16' in text
    assert '/unrelated' not in text


def test_reaction_and_dispatch_blockers_visible(production):
    root, case = production
    put(case/'well/kinbot_monitor.out', '4\t6\thom_sci\n-1\t5\taccepted_rxn\n-999\t2\tfailed_rxn\n')
    put(case/'well/vrctst/dispatch/workflow.json', {'tasks': [{'id': 'first'}, {'id': 'next', 'depends_on': ['first']}]})
    put(case/'well/vrctst/dispatch/state.json', {'tasks': {'first': {'status': 'failed', 'error': 'QPErr'}}})
    text = monitor.render(monitor.snapshot(root), include_all=True)
    assert 'product/TS optimization' in text
    assert "waiting for={'first': 'failed'}" in text
    assert 'accepted_rxn' in text


def test_malformed_json_is_reported_and_log_view_survives(production):
    root, case = production
    put(case/'well/state.json', '{"tasks":')
    put(case/'well/kinbot.log', '16-Oct-26 18:35:07-INFO: Starting KinBot\n')
    assert 'READ WARNING:' in monitor.render(monitor.snapshot(root))


def test_follower_append_partial_line_rotation_truncation_and_new_worker(tmp_path):
    path = tmp_path/'pes.log'
    path.write_text('first\npartial')
    follower = monitor.LogFollower()
    assert follower.poll(tmp_path)[0] == [('pes.log', ['first'])]
    assert follower.poll(tmp_path)[0] == []
    with path.open('a') as f:
        f.write(' end\nnext\n')
    assert follower.poll(tmp_path)[0] == [('pes.log', ['partial end', 'next'])]
    path.write_text('reset\n')
    assert follower.poll(tmp_path)[0] == [('pes.log', ['reset'])]
    path.rename(tmp_path/'old.txt')
    path.write_text('rotated\n')
    assert follower.poll(tmp_path)[0] == [('pes.log', ['rotated'])]
    put(tmp_path/'new/kinbot.log', 'new worker\n')
    assert follower.poll(tmp_path)[0] == [('new/kinbot.log', ['new worker'])]


def test_archives_and_external_symlinks_not_followed(tmp_path):
    put(tmp_path/'.restart_test/kinbot.log', 'old failure\n')
    put(tmp_path/'archive/pes.log', 'old failure\n')
    (tmp_path/'external').symlink_to(tmp_path.parent, target_is_directory=True)
    assert monitor.inventory(tmp_path, []) == []


def test_rotd_surfaces_and_execution_failures_visible(production):
    root, case = production
    put(case/'well/rotdPy/kb_reaction/Surface_0/large_payload.log', 'not traversed')
    put(case/'well/rotdPy/kb_reaction/output/surface_0.dat', 'flux')
    put(case/'well/rotdPy/reaction.rotdpy.json', {'status': 'failed', 'error': 'sample limit'})
    text = monitor.render(monitor.snapshot(root))
    assert '1/1 surface output files' in text
    assert 'sample limit' in text
    assert 'large_payload' not in text


def test_follower_handles_native_out_and_overlong_lines(tmp_path):
    path = tmp_path/'molpro.out'
    path.write_text('already shown\n')
    follower = monitor.LogFollower()
    assert follower.poll(tmp_path, 0)[0] == []
    with path.open('a') as stream:
        stream.write('x'*monitor.TAIL_BYTES+'\nnew progress\n')
    assert follower.poll(tmp_path)[1]
    assert 'new progress' in follower.poll(tmp_path)[0][0][1]
