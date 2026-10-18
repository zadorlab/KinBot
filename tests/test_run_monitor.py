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


def test_archived_failure_is_cached_outcome_but_new_attempt_is_visible(production):
    from kinbot.execution import file_sha256
    root, case = production
    marker = case / 'well/rejected.exitcode'
    put(marker, '42 1\n')
    put(case / 'resume_provenance.json', {
        'archived_failure_exits': {'well/rejected.exitcode': file_sha256(marker)}})
    text = monitor.render(monitor.snapshot(root))
    assert 'Cached rejected QC outcomes: 1' in text
    assert 'FAILED QC EXIT' not in text
    put(marker, '43 1\n')
    text = monitor.render(monitor.snapshot(root))
    assert 'FAILED QC EXIT' in text


def test_anl_curvature_recovery_report_is_visible(production):
    root, case = production
    put(case/'anl/tasks/l3_geometry/minimum_recovery.json', {
        'status': 'minimum_curvature_accepted', 'controls': {'attempts': 2},
        'attempts': [{'lowest_frequency_cm_inverse': -587.,
                     'selected_energy_eV': -10.1},
                    {'lowest_frequency_cm_inverse': 125.}]})
    text = monitor.render(monitor.snapshot(root), include_all=True)
    assert 'minimum_curvature_accepted' in text
    assert 'curvature_checks' in text and '-587.0' in text


def test_scf_recovery_visible_without_archived_native_error(production):
    root, case = production
    put(case/'well/atom.scf_recovery.json', {'status': 'complete',
        'reason': 'Gaussian QC-SCF failure in Link 508',
        'original_scf': 'xqc,conver=10',
        'retry_scf': 'conver=10,CDIIS,NoVarAcc,NoIncFock,MaxCycle=256'})
    put(case/'well/restart_gaussian_scf/atom.123/atom.log',
        'No lower point found -- run aborted.\nError termination via l508.exe')
    text = monitor.render(monitor.snapshot(root), include_all=True)
    assert 'CDIIS,NoVarAcc,NoIncFock' in text and 'complete' in text
    assert 'No lower point found' not in text


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


def test_sella_progress_and_frequency_rejection_are_visible(production):
    root, case = production
    put(case/'well/well_high_sella.log',
        'Sella   6 19:11:56 -11916.528873 0.0055 0.0000 0.1000 1.2245 6\n'
        'Worker failure: Frequency validation failed for minimum: imaginary '
        'modes [cm^-1] = -587.3340; 1 below -50 cm^-1.\ndone\n')
    put(case/'well/well_high.log',
        ' Frequencies --  -587.3340 727.3799 979.8327\n'
        ' Normal termination of Gaussian 16\n')
    text = monitor.render(monitor.snapshot(root))
    assert 'Sella   6' in text
    assert 'Frequency validation failed for minimum' in text
    assert '-587.3340' in text


def test_older_failed_product_spectrum_not_hidden_by_newer_jobs(production):
    import os
    root, case = production
    put(case/'well/product_well.exitcode', '52908315 1\n')
    failure = case/'well/product_well_sella.log'
    put(failure, 'Converged but frequency check failed (order=0); frequencies (cm-1):\n'
        '[np.float64(-53.062370746758766), np.float64(303.6845755216985)]\ndone\n')
    os.utime(failure, (1., 1.))
    for i in range(45):
        put(case/f'well/recent_{i}.log', 'Normal termination of Gaussian 16\n')
    text = monitor.render(monitor.snapshot(root))
    assert 'FAILED QC EXIT' in text
    assert 'Converged but frequency check failed' in text
    assert '-53.062370746758766' in text


def test_both_pes_remain_visible_during_targeted_minimum_check(production, monkeypatch):
    root, first = production
    second = root/'second_pes'
    put(root/'prepared.json', {'inputs': [str(first/'input.json'), str(second/'input.json')]})
    put(second/'input.json', {'queue_job_limit': 10})
    put(second/'well/kinbot.log', '16-Oct-26 20:00:00-INFO: Starting rotor scan\n')
    check = first/'minimum_recovery_check'
    put(check/'well_high_sella.log', 'Minimum recovery 1/2: selected plus displacement\n'
        'Sella   7 20:00:00 -10.5 0.0001 0.0\n')
    archive = check/'well_high_minimum_recovery/attempt_1_example'
    put(archive/'well_high.log', 'Error termination HISTORICAL ONLY\n')
    put(archive/'recovery.json', {'status': 'descent_selected', 'attempt': 1,
        'selected': 'plus', 'energy_change_eV': -.0002})
    job = {'job_id': '43', 'state': 'RUNNING', 'nodes': 'node20', 'elapsed': '1:00',
           'limit': '24:00:00', 'cpus': '8', 'memory': '10G', 'excluded_nodes': 'node16',
           'name': 'well_high', 'reason': 'node20', 'directory': str(check)}
    monkeypatch.setattr(monitor, 'queue', lambda: ([job, dict(job, job_id='44', directory=str(second))], None))
    report = monitor.snapshot(root)
    assert len(report['cases']) == 2
    assert [c['jobs'][0]['job_id'] for c in report['cases']] == ['43', '44']
    text = monitor.render(report)
    assert 'Starting rotor scan' in text and 'Minimum recovery 1/2' in text
    assert 'descent_selected' in text and 'HISTORICAL ONLY' not in text


def test_dashboard_separates_both_kinbot_tails_from_qc_without_parameter_dump(production):
    root, first = production
    second = root/'second_pes'
    put(root/'prepared.json', {'inputs': [str(first/'input.json'), str(second/'input.json')]})
    for case in (first, second):
        put(case/'input.json', {'queue_job_limit': 10, 'slurm_exclude_nodes': ['node16']})
        put(case/'well/kinbot_status.json', {'status': 'running'})
        put(case/'well/kinbot.log',
            '    vrc_tst_surface_model: '+str({'convergence': 1e-8, 'payload': 'x'*12000})+'\n'
            '    anl_cc_convergence: 1e-8\n'
            '16-Oct-26 21:04:01-INFO: Starting optimization of initial well...\n'
            '16-Oct-26 21:04:02-WARNING: Runtime warning stays visible\n')
        put(case/'well/energy.log', ' SCF Done: E(UB3LYP) = -437.123 A.U.\n'
            ' 3 16 0 1.054795 3.058235 -0.304045\n')
        put(case/'well/energy_sella.log', ' Step Time Energy fmax cmax rtrust rho trjid\n'
            ' Sella 4 21:05:06 -11920.123 0.0521 0.0 0.1 1.2 4\n')
    before = {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}
    report = monitor.snapshot(root, include_all=True)
    text = monitor.render(report, include_all=True, width=100)
    assert text.count('KINBOT LOGS') == text.count('QC OUTPUT') == 2
    assert text.count('Starting optimization of initial well') == 2
    assert text.count('Runtime warning stays visible') == 2
    assert text.count('SCF Done: E(UB3LYP)') == text.count('Sella 4') == 2
    assert text.count('Configured exclusions: node16') == 1
    assert 'kinbot_status.json' not in text  # status already on WORKER row
    assert text.count('WORKER well: running') == 2
    assert 'surface_model' not in text and 'anl_cc_convergence' not in text
    assert '1.054795' not in text  # coordinates are not solver iterations
    assert 'Step Time Energy fmax' in text
    assert len(text) < 5000 and max(map(len, text.splitlines())) <= 100
    assert text.index('KINBOT LOGS') < text.index('QC OUTPUT')
    assert before == {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}


def test_numeric_native_iterations_and_error_reasons_remain_visible(tmp_path):
    path = tmp_path/'molpro.out'
    put(path, ' PROGRAM * MRCI-F12\n'
        ' ITER. STATE ROOT SQ.NORM CORR.ENERGY TOTAL ENERGY ENERGY CHANGE\n'
        ' 8 1 2 1.11667027 -0.33129541 -437.41636132 -0.00000002 0.32D-08\n'
        ' ? Error\n ? dsygv failed\n ? The problem occurs in cif12solve\n'
        ' GLOBAL ERROR fehler on processor 0\n')
    record = monitor.log_record(path, tmp_path, path.stat().st_mtime)
    content = monitor.log_content(record, lines=8)
    assert any('ITER. STATE ROOT' in line for line in content)
    assert any('-437.41636132' in line for line in content)
    assert any('dsygv failed' in line for line in content)
    assert any('cif12solve' in line for line in content)


def test_details_and_raw_tails_retain_full_metadata(production, monkeypatch):
    root, case = production
    command = f'/opt/python -u -m kinbot.pes {case}/input.json'
    from kinbot import run_cleanup
    monkeypatch.setattr(run_cleanup, '_processes', lambda *a: [
        {'pid': 123, 'state': 'S', 'cwd': str(case), 'command': command}])
    put(case/'well/kinbot.log', '16-Oct-26 21:04:01-INFO: Starting KinBot\n')
    put(case/'well/native.log', 'unclassified backend detail\n1 6 0 1.2 3.4 5.6\n')
    put(case/'well/kinbot_status.json', {'status': 'running'})
    report = monitor.snapshot(root, log_lines=20)
    compact = monitor.render(report, width=1000)
    expanded = monitor.render(report, details=True, raw_logs=True, lines=20, width=1000)
    assert 'kinbot.pes PID 123 [S]' in compact and command not in compact
    assert command in expanded
    assert 'kinbot_status.json' in expanded
    assert 'unclassified backend detail' in expanded
    assert '1 6 0 1.2 3.4 5.6' in expanded
    assert report['cases'][0]['processes'][0]['command'] == command
    assert 'unclassified backend detail' in json.dumps(report)


def test_lines_option_applies_to_collection_and_live_cli(production, capsys):
    import argparse
    root, case = production
    put(case/'well/kinbot.log', ''.join(
        f'16-Oct-26 21:04:01-INFO: Event {i}\n' for i in range(25)))
    parser = argparse.ArgumentParser()
    monitor.add_arguments(parser)
    args = parser.parse_args(['--all', '--details', '--raw-logs', '--lines', '20'])
    monitor.run(root, args)
    text = capsys.readouterr().out
    assert 'Event 5\n' in text and 'Event 24\n' in text
    assert 'Event 4\n' not in text


def test_narrow_display_wraps_long_errors_without_losing_tail(production):
    root, case = production
    message = 'RuntimeError: '+('long diagnostic '*30)+'END OF ERROR'
    put(case/'well/perm/failed.err', message+'\n')
    text = monitor.render(monitor.snapshot(root), width=70)
    assert max(map(len, text.splitlines())) <= 70
    assert 'END OF ERROR' in text


def test_many_composite_records_do_not_bury_second_pes(tmp_path, monkeypatch):
    root = tmp_path / 'both'
    cases = [root / name for name in ('first', 'second')]
    put(root / 'prepared.json', {'inputs': [str(c / 'input.json') for c in cases]})
    monkeypatch.setattr(monitor, 'queue', lambda: ([], None))
    from kinbot import run_cleanup
    monkeypatch.setattr(run_cleanup, '_processes', lambda *a: [])
    for case in cases:
        put(case / 'input.json', {'queue_job_limit': 10})
        put(case / 'well/kinbot.log', '17-Oct-26 00:09:24-ERROR: branch failure\n')
        put(case / 'well/.kinbot_post_reaction_branches.json', {
            'status': 'failed', 'branches': {'vrc_tst': {'status': 'failed', 'error': 'geometry mismatch'}}})
        for i in range(40):
            put(case / f'well/task{i}.result.json', {'status': 'normal'})
            put(case / f'well/composite{i}/state.json', {'tasks': {
                f'stage{j}': {'status': 'staged'} for j in range(25)}})
            put(case / f'well/composite{i}/vpt2.log', 'Normal termination of Gaussian\n')
    report = monitor.snapshot(root, include_all=True)
    compact = monitor.render(report, include_all=True, width=120)
    expanded = monitor.render(report, include_all=True, details=True, width=120)
    assert compact.index('OVERVIEW second') < compact.index('=== first ===')
    assert compact.count('KINBOT LOGS') == 2
    assert compact.count('QC OUTPUT') == 2
    assert len(compact.splitlines()) < 150
    assert len(expanded.splitlines()) > 1500
    assert len(report['cases'][0]['artifacts']) == 81  # no data deleted to simplify display


def test_rehomed_run_monitors_both_cases_and_old_slurm_path(tmp_path, monkeypatch):
    root, old = tmp_path/'new', tmp_path/'old'
    first, second = root/'first', root/'second'
    put(root/'prepared.json', {'inputs': [str(first/'input.json'), str(second/'input.json')]})
    for case in (first, second):
        put(case/'input.json', {'queue_job_limit': 10})
    old.symlink_to(root, target_is_directory=True)
    raw = f'42|RUNNING|node20|00:30|31-00:00:00|8|10G||l3|node20|{old}/first/123\n'
    monkeypatch.setattr(monitor.subprocess, 'run', lambda *a, **kw: SimpleNamespace(returncode=0, stdout=raw, stderr=''))
    monkeypatch.setattr(monitor, 'queue', ORIGINAL_QUEUE)
    monkeypatch.setattr('kinbot.run_cleanup._processes', lambda *a: [])
    report = monitor.snapshot(root)
    assert len(report['cases']) == 2
    text = monitor.render(report)
    assert '42' in text and 'first' in text and 'second' in text
