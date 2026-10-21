import json
import os
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


def test_retired_vrc_graph_hidden_only_when_unchanged_and_not_live(production, monkeypatch):
    from kinbot.execution import file_sha256
    root, case = production
    graph = case/'123/vrctst/molpro_scan/old_tz'
    put(graph/'workflow.json', {'tasks': [{'id': 'relaxed_scan'}]})
    put(graph/'state.json', {'tasks': {'relaxed_scan': {'status': 'failed', 'error': 'RETIRED FAILURE'}}})
    put(graph/'tasks/relaxed_scan/old.out', 'Error termination RETIRED FAILURE\n')
    put(root/'vrc_restart.json', {'retired_scans': [{'graph': str(graph),
         'workflow_sha256': file_sha256(graph/'workflow.json')}]})
    before = {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}
    report = monitor.snapshot(root)
    assert 'RETIRED FAILURE' not in json.dumps(report)
    assert 'retained for comparison: 1' in monitor.render(report)
    assert before == {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}
    monkeypatch.setattr(monitor, 'queue', lambda: ([{'job_id': '100', 'state': 'RUNNING',
        'directory': str(graph/'tasks/relaxed_scan')}], None))
    assert 'RETIRED FAILURE' in json.dumps(monitor.snapshot(root))
    monkeypatch.setattr(monitor, 'queue', lambda: ([], None))
    put(graph/'workflow.json', {'tasks': [{'id': 'relaxed_scan', 'changed': True}]})
    assert 'RETIRED FAILURE' in json.dumps(monitor.snapshot(root))


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


def test_live_pool_is_shown_instead_of_stale_per_case_limit(production):
    from kinbot import run_policy
    root, case = production
    run_policy.update(root, max_jobs=23, exclude=['reserve'])
    before = {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}
    text = monitor.render(monitor.snapshot(root))
    assert 'Shared pool: cap 23' in text and '/ shared pool;' in text
    assert '/ cap 10' not in text
    assert 'New-submission exclusions: node16,reserve' in text
    assert before == {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}


def test_selection_display_distinguishes_old_controller_and_current_choice(production):
    root, case = production
    phase = dict(token='new', epoch=1, phase='exploring', workers={})
    put(case/'pes_exploration.json', phase)
    assert 'New ANL work waits' not in monitor.render(monitor.snapshot(root))
    phase['anl_selection_policy']='lowest_l2'
    put(case/'pes_exploration.json', phase)
    put(case/'pes_l2_selection.json',dict(token='old',epoch=1,status='selected',
        workers={'a':dict(selected=['ts'],omitted={})}))
    text=monitor.render(monitor.snapshot(root))
    assert 'New ANL work waits' in text and '1 retained' not in text
    put(case/'pes_l2_selection.json',dict(token='new',epoch=1,status='selected',
        workers={'a':dict(selected=['ts'],omitted={'other':{}})}))
    assert '1 retained; 1 omitted' in monitor.render(monitor.snapshot(root))


def test_exploration_slot_wait_is_distinct_from_qc_pool(production):
    root, case = production
    put(case/'pes_exploration.json', dict(token='new', epoch=0, phase='exploring',
        max_searches=1, waiting_workers=['c'], exploring_workers=['b'],
        workers={'a':dict(status='vrc'), 'b':dict(status='exploring')}))
    text = monitor.render(monitor.snapshot(root))
    assert 'Exploration slots: 1/1; unstarted wells=1' in text
    assert 'Next well waits for an exploring well to finish: b' in text


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


def test_legacy_live_preflight_import_is_identified(production, monkeypatch):
    from kinbot import run_cleanup
    root, case = production
    graph = case / 'well/anl_automatic/states/state_x/composite'
    monkeypatch.setattr(run_cleanup, '_processes', lambda *a: [{
        'pid': 42, 'state': 'R', 'cwd': str(graph),
        'command': '/venv/bin/python -c import ase, sella; import kinbot.anl.dispatch'}])
    text = monitor.render(monitor.snapshot(root))
    assert 'Composite preflight: checking Python environment' in text
    assert 'states/state_x/composite' in text


@pytest.mark.parametrize('alive', [False, True])
def test_preflight_progress_is_visible_and_stale_record_is_labelled(production, monkeypatch, alive):
    from kinbot import run_cleanup
    root, case = production
    put(case / 'well/anl_automatic/composite_preflight.json', {
        'status': 'running', 'controller_pid': 42, 'phase': 'checking_slurm_script',
        'graphs_complete': 3, 'graphs_total': 13, 'graph': 'stable:C|0|1',
        'task': 'l3_geometry', 'updated_at': monitor.time.time()})
    monkeypatch.setattr(run_cleanup, '_processes', lambda *a: [{
        'pid': 42, 'state': 'S', 'cwd': str(case),
        'command': '/venv/bin/python -m kinbot.kb input.json'}] if alive else [])
    text = monitor.render(monitor.snapshot(root))
    assert 'graphs checked 3/13' in text and 'checking_slurm_script' in text
    assert 'task=l3_geometry' in text
    assert ('recorded state only' in text) == (not alive)


@pytest.mark.parametrize('worker_status', ['running', 'failed'])
def test_previous_branch_report_is_not_a_new_failure(production, worker_status):
    root, case = production
    worker = case/'well'
    put(worker/'kinbot_status.json', {
        'status': worker_status, 'pid': 42,
        'started_utc': '2026-10-21T02:31:48+00:00',
        **({'error': 'NEW CONTROLLER FAILURE'} if worker_status == 'failed' else {})})
    branch = worker/'.kinbot_post_reaction_branches.json'
    put(branch, {'status': 'running', 'mess_gate': 'blocked', 'branches': {
        'cbh_anl': {'status': 'running'},
        'vrc_tst': {'status': 'failed', 'error': 'OLD BOUND GRAPH FAILURE'}}})
    preflight = worker/'anl_automatic/composite_preflight.json'
    put(preflight, {'status': 'complete', 'controller_pid': 41,
                    'phase': 'ready_for_dispatch', 'graphs_complete': 19, 'graphs_total': 19})
    previous = monitor.datetime.fromisoformat('2026-10-21T02:30:00+00:00').timestamp()
    for p in (branch, preflight):
        os.utime(p, (previous, previous))
    before = {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}
    report = monitor.snapshot(root)
    text = monitor.render(report)
    assert 'vrc_tst=not_reported_yet' in text
    assert 'cbh_anl=not_reported_yet' in text
    assert 'OLD BOUND GRAPH FAILURE' not in text
    assert 'MESS gate: blocked' not in text
    assert 'Composite preflight: previous controller report' in text
    assert ('NEW CONTROLLER FAILURE' in text) == (worker_status == 'failed')
    detailed = monitor.render(report, details=True)
    assert '.kinbot_post_reaction_branches.json [previous controller]' in detailed
    assert 'OLD BOUND GRAPH FAILURE' in detailed
    assert 'OLD BOUND GRAPH FAILURE' in json.dumps(report)
    assert before == {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}


@pytest.mark.parametrize('started', ['2026-10-21T02:31:48+00:00', None, 'invalid'])
def test_current_branch_failure_remains_visible_after_controller_exits(production, started):
    root, case = production
    worker = case/'well'
    branch = worker/'.kinbot_post_reaction_branches.json'
    put(branch, {'status': 'failed', 'mess_gate': 'blocked', 'branches': {
        'vrc_tst': {'status': 'failed', 'error': 'CURRENT FAILURE'}}})
    written = monitor.datetime.fromisoformat('2026-10-21T02:31:59+00:00').timestamp()
    os.utime(branch, (written, written))
    # This file is rewritten AFTER the branch failure at controller exit.
    put(worker/'kinbot_status.json', {'status': 'failed', 'pid': 42,
                                      'started_utc': started})
    report = monitor.snapshot(root)
    text = monitor.render(report)
    assert 'vrc_tst=failed' in text and 'CURRENT FAILURE' in text
    assert 'MESS gate: blocked' in text
    assert 'previous controller' not in text


def test_new_branch_report_replaces_startup_history(production):
    root, case = production
    worker = case/'well'
    put(worker/'kinbot_status.json', {'status': 'running', 'pid': 42,
        'started_utc': '2026-10-21T02:31:48+00:00'})
    branch = worker/'.kinbot_post_reaction_branches.json'
    put(branch, {'status': 'running', 'mess_gate': 'waiting', 'branches': {
        'cbh_anl': {'status': 'running'}, 'vrc_tst': {'status': 'running'}}})
    written = monitor.datetime.fromisoformat('2026-10-21T02:31:59+00:00').timestamp()
    os.utime(branch, (written, written))
    text = monitor.render(monitor.snapshot(root))
    assert 'cbh_anl=running; vrc_tst=running' in text
    assert 'MESS gate: waiting' in text
    assert 'previous controller' not in text


def test_unstarted_child_branch_uses_replacement_pes_start(production):
    root, case = production
    put(case/'pes_status.json', {'status':'running', 'started_utc':'2026-10-21T17:32:53+00:00'})
    branch = case/'well/.kinbot_post_reaction_branches.json'
    put(branch, {'branches':{'cbh_anl':{'status':'failed','error':'OLD CHILD FAILURE'}}})
    previous = monitor.datetime.fromisoformat('2026-10-21T15:00:00+00:00').timestamp()
    os.utime(branch,(previous,previous))
    text = monitor.render(monitor.snapshot(root))
    assert 'cbh_anl=not_reported_yet' in text
    assert 'OLD CHILD FAILURE' not in text
    assert 'OLD CHILD FAILURE' in monitor.render(monitor.snapshot(root),details=True)
    # An actual branch failure after the parent launch must remain current.
    current = monitor.datetime.fromisoformat('2026-10-21T17:40:00+00:00').timestamp()
    os.utime(branch,(current,current))
    assert 'cbh_anl=failed' in monitor.render(monitor.snapshot(root))


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


def test_archived_dispatch_failure_excluded_but_current_retry_failure_visible(production):
    root, case = production
    graph = case/'well/anl_automatic/states/state_x/composite'
    put(graph/'workflow.json', {'tasks': [{'id': 'gaussian_vpt2',
        'recovery': {'assigned_native_framework': 'C1'}}]})
    put(graph/'state.json', {'tasks': {'gaussian_vpt2': {'status': 'staged'}}})
    archive = graph/'attempts/gaussian_vpt2/1'
    put(archive/'execution.json', {'status': 'failed', 'error': 'OLD segmentation violation'})
    put(archive/'vpt2.log', 'Error termination OLD native failure\n')
    put(archive/'task.exitcode', '12 1\n')
    before = {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}
    report = monitor.snapshot(root)
    for details in (False, True):
        text = monitor.render(report, details=details)
        assert 'OLD' not in text and 'FAILED QC EXIT' not in text
    assert 'RETRY ANL/state_x/gaussian_vpt2: staged' in monitor.render(report)
    assert 'OLD' not in json.dumps(report)
    assert len(report['cases'][0]['archived_attempt_directories']) == 1
    assert not any('attempts/' in path for path, _ in monitor.LogFollower().poll(root)[0])
    assert before == {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}

    put(graph/'state.json', {'tasks': {'gaussian_vpt2': {
        'status': 'unavailable', 'unavailable_reason': 'NEW native failure'}}})
    put(graph/'tasks/gaussian_vpt2/execution.json', {'status': 'failed', 'error': 'NEW native failure'})
    put(graph/'tasks/gaussian_vpt2/vpt2.log', 'QPErr NEW syntax error\n')
    text = monitor.render(monitor.snapshot(root))
    assert 'unavailable=1' in text and 'NEW native failure' in text and 'QPErr' in text
    assert 'OLD' not in text


def test_failure_cascade_is_bounded_without_losing_logs_or_json(production):
    root, case = production
    error = "'getpwuid(): uid not found: 155382'"
    for index in range(20):
        graph = case/f'123/anl_automatic/states/state_{index}/composite'
        tasks = [{'id': 'geometry'}] + [
            {'id': f'sp_{n}', 'depends_on': ['geometry']} for n in range(15)]
        put(graph/'workflow.json', {'tasks': tasks})
        put(graph/'state.json', {'tasks': {'geometry': {'status': 'failed'}, **{
            f'sp_{n}': {'status': 'blocked', 'blocked_by': ['geometry']} for n in range(15)}}})
        put(graph/'tasks/geometry/execution.json', {'status': 'failed', 'error': error})
        put(graph/'tasks/geometry/minimum_recovery.json', {'status': 'failed', 'error': 'KeyError: '+error})
    put(case/'123/kinbot.log', '19-Oct-26 13:52:00-INFO: Continuing ANL tasks\n')
    put(case/'123/current.out', 'PROGRAM * CCSD\n')
    report = monitor.snapshot(root)
    text = monitor.render(report, width=150)
    assert 'failed=20' in text and 'blocked=300' in text
    assert '20 affected' in text and 'Dependency blocks: 300' in text
    assert text.count('getpwuid') == 1
    assert 'Continuing ANL tasks' in text and 'PROGRAM * CCSD' in text
    assert len(text.splitlines()) < 55
    assert len([a for a in report['cases'][0]['artifacts'] if 'tasks' in a]) == 20
    assert 'sp_14' in monitor.render(report, details=True)


def test_intermediate_execution_record_does_not_duplicate_live_retry_failure(production, monkeypatch):
    root, case = production
    graph = case/'123/anl_automatic/states/state_x/composite'
    put(graph/'workflow.json', {'tasks': [{'id': 'geometry'}]})
    put(graph/'state.json', {'tasks': {'geometry': {'status': 'submitted', 'job_id': '42'}}})
    put(graph/'tasks/geometry/execution.json', {'status': 'failed', 'error': 'old attempt'})
    job = {'job_id': '42', 'state': 'RUNNING', 'nodes': 'node20', 'elapsed': '1:00',
           'limit': '24:00:00', 'cpus': '8', 'memory': '10G', 'excluded_nodes': 'node16',
           'name': 'geometry', 'reason': 'node20', 'directory': str(graph/'tasks/geometry')}
    monkeypatch.setattr(monitor, 'queue', lambda: ([job], None))
    report = monitor.snapshot(root)
    text = monitor.render(report)
    assert 'Live tasks with earlier/intermediate failure records: 1' in text
    assert 'Failures:' not in text
    assert 'old attempt' in monitor.render(report, details=True)
    monkeypatch.setattr(monitor, 'queue', lambda: ([], None))
    assert 'old attempt' in monitor.render(monitor.snapshot(root))


def test_only_hash_matched_cached_worker_error_is_collapsed(production):
    from kinbot.execution import file_sha256
    root, case = production
    marker = case/'well/rejected.exitcode'
    message = '18-Oct-26 14:46:41-ERROR: rejected exited with worker status 1.'
    put(marker, '42 1\n')
    put(case/'resume_provenance.json', {
        'archived_failure_exits': {'well/rejected.exitcode': file_sha256(marker)}})
    put(case/'well/kinbot.log', message+'\n18-Oct-26 14:47:00-INFO: VRC started\n'
        '18-Oct-26 14:47:01-ERROR: DIFFERENT failure remains visible\n')
    # An identically named job in another worker is not the archived failure.
    put(case/'other/kinbot.log', message+'\n')
    report = monitor.snapshot(root)
    text = monitor.render(report)
    assert text.count('rejected exited with worker status 1') == 1
    assert 'DIFFERENT failure remains visible' in text and 'VRC started' in text
    assert 'Cached rejected QC outcomes: 1' in text
    assert monitor.render(report, details=True).count('rejected exited with worker status 1') == 2
    assert message in json.dumps(report)
    put(marker, '43 1\n')
    text = monitor.render(monitor.snapshot(root))
    assert text.count('rejected exited with worker status 1') == 2
    assert 'FAILED QC EXIT' in text


def test_completed_reference_and_rejections_have_details_without_crowding_live_logs(production):
    root, case = production
    put(case/'well/kinbot.log', '18-Oct-26 14:47:00-INFO: Live main worker\n')
    put(case/'well/kinbot_monitor.out', '-999 23 rejected_rxn\n-1 5 accepted_rxn\n4 6 live_rxn\n')
    ref = case/'well/anl_automatic/reference_l2/state_x'
    put(ref/'kinbot_status.json', {'status': 'complete'})
    put(ref/'kinbot.log', '18-Oct-26 12:36:05-WARNING: Stale reference advice\n')
    report = monitor.snapshot(root)
    text = monitor.render(report)
    assert "'failed/rejected': 1" in text and 'Reference/subtask workers: 1 complete' in text
    assert 'Live main worker' in text and 'live_rxn' in text
    assert 'Stale reference advice' not in text and 'rejected_rxn' not in text
    details = monitor.render(report, details=True)
    assert 'Stale reference advice' in details and 'rejected_rxn' in details
    put(case/'well/kinbot_status.json', {'status': 'complete'})
    assert 'Live main worker' in monitor.render(monitor.snapshot(root))


def test_each_task_gets_log_budget_before_repeated_force_steps(production, monkeypatch):
    import os
    root, case = production
    busy = case/'well/anl_automatic/states/state_busy/composite/tasks/l3_geometry'
    slow = case/'well/anl_automatic/states/state_slow/composite/tasks/gaussian_vpt2'
    for i in range(60):
        path = busy/f'step_{i:04d}.out'
        put(path, f' !RHF-UCCSD(T) energy {-400-i}.123\n')
        os.utime(path, (100+i, 100+i))
    put(slow/'vpt2.log', ' SCF Done: E(UB2PLYP) = -398.123 A.U.\n')
    os.utime(slow/'vpt2.log', (1, 1))
    job = {'job_id': '42', 'state': 'RUNNING', 'nodes': 'node20', 'elapsed': '1:00',
           'limit': '24:00:00', 'cpus': '8', 'memory': '10G', 'excluded_nodes': 'node16',
           'name': 'vpt2', 'reason': 'node20', 'directory': str(slow)}
    monkeypatch.setattr(monitor, 'queue', lambda: ([job], None))
    report = monitor.snapshot(root)
    text = monitor.render(report)
    assert '-398.123' in text and '-459.123' in text
    assert '-458.123' not in text  # prior step of the same task is details-only
    assert 'READ WARNING' not in text
    assert report['cases'][0]['log_coverage'] == {'read': 40, 'available': 61}


def test_molpro_dataset_numbers_are_not_solver_iterations(tmp_path):
    path = tmp_path/'molpro.out'
    put(path, ' ITER. SQ.NORM CORR.ENERGY TOTAL ENERGY ENERGY CHANGE\n'
        ' 13 1.18076383 -0.41430535 -437.48412474 0.00000079 -0.00000000 0.95D-10\n'
        ' 1100 1400 1410 1200 1210 1080 1600 1650 1700 1380\n'
        ' 3 16 0 1.054795 3.058235 -0.304045\n')
    record = monitor.log_record(path, tmp_path, path.stat().st_mtime)
    text = '\n'.join(monitor.log_content(record))
    assert '-437.48412474' in text and 'ITER.' in text
    assert '1100 1400' not in text and '1.054795' not in text
    assert '1100 1400' in '\n'.join(monitor.log_content(record, raw=True))
