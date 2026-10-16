"""Live monitoring must surface real evidence without touching the run."""
import json
from pathlib import Path
from types import SimpleNamespace
import pytest

from kinbot import orbital_monitor as monitor


@pytest.mark.parametrize('cause', [
    '?TRYING TO OVERWRITE DUMP RECORD 2140.2 WITH DIFFERENT BASIS DIMENSION',
    '? RHF mode requires a single determinant wavefunction!',
    '? not sure if it works with nstatr/=nstati',
])
def test_monitor_preserves_native_cause_before_mpi_teardown(tmp_path, monkeypatch, cause):
    jobs = _campaign(tmp_path)
    (Path(jobs[0]['directory']) / 'probe.out').write_text(
        ' PROGRAM * MULTI\n ? Error\n ?\n ' + cause + '\n'
        ' ? The problem occurs in native_solver\n GLOBAL ERROR fehler on processor 0\n')
    monkeypatch.setattr(monitor, '_queue', lambda: ({}, None))
    result = monitor.snapshot(tmp_path)
    assert cause in monitor.render(result, details=True)


def _campaign(tmp_path):
    config = tmp_path / 'input.json'
    config.write_text('{"probe_method":"casscf"}')
    jobs = []
    for index in range(4):
        directory = tmp_path / str(index)
        directory.mkdir()
        jobs.append({'name': 'probe', 'directory': str(directory),
                     'status': 'complete' if index == 2 else 'incomplete'})
        if index < 3:
            jobs[-1]['job_id'] = str(index + 10)
    path = tmp_path / 'reconnaissance.json'
    path.write_text(json.dumps({'jobs': jobs, 'slurm_excluded_nodes': ['bad-node']}))
    empty = tmp_path / 'full.json'
    empty.write_text('{"jobs":[]}')
    (tmp_path / 'campaign.json').write_text(json.dumps({
        'status': 'running', 'max_jobs': 20, 'selected_probe_method': 'casscf',
        'systems': [{'title': 'example', 'configuration': str(config),
                     'pilot_manifest': str(path), 'manifest': str(empty)}]}))
    return jobs


def test_snapshot_reads_live_native_errors_without_mutating_or_auditing(tmp_path, monkeypatch):
    jobs = _campaign(tmp_path)
    directory = Path(jobs[0]['directory'])
    (directory / 'execution.json').write_text('{"status":"running","host":"old-node"}')
    (directory / 'probe.out').write_text(' PROGRAM * RS2-F12\n'
        ' NO CONVERGENCE IN REFERENCE CI AFTER200 ITERATIONS\n'
        + '\n'.join(' GLOBAL ERROR fehler on processor ' + str(i) for i in range(15)))
    calls = []
    def scheduler(command, **kwargs):
        calls.append(command)
        return SimpleNamespace(returncode=0, stdout=(
            f'10|RUNNING|bad-node|2:10|None|104000M||{directory}\n'
            f'11|PENDING||0:00|Resources|104000M|bad-node|{jobs[1]["directory"]}\n'))
    monkeypatch.setattr(monitor.subprocess, 'run', scheduler)
    before = {p: p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    result = monitor.snapshot(tmp_path)
    assert result['counts']['pilot_manifest'] == {
        'total': 4, 'RUNNING': 1, 'PENDING': 1, 'VALIDATED': 1, 'WAITING': 1}
    first, pending = result['jobs']
    assert first['state'] == 'RUNNING'  # Log evidence does not relabel a live allocation.
    assert first['stage'] == 'RS2-F12'
    assert 'NO CONVERGENCE' in first['log_matches'][0]['line']
    assert any('ASSIGNED TO EXCLUDED NODE' in a for a in first['alerts'])
    assert any('no node exclusion' in a for a in first['alerts'])
    assert not pending['alerts'] and not pending['log_matches']
    assert pending['output_age_seconds'] is None
    rendered = monitor.render(result, details=True)
    assert '104000M' in rendered and 'NO CONVERGENCE' in rendered and 'Resources' in rendered
    assert 'quiet output alone does not prove a stall' in rendered
    assert len(calls) == 1 and calls[0][0] == 'squeue'
    assert {p: p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()} == before


def test_monitor_reads_explicit_single_phase_registry_without_duplicate_jobs(tmp_path, monkeypatch):
    _campaign(tmp_path)
    original = json.loads((tmp_path / 'campaign.json').read_text())
    for system in original['systems']:
        system['manifest'] = system.pop('pilot_manifest')
    registry = tmp_path / 'custom_registry.json'
    registry.write_text(json.dumps(original))
    monkeypatch.setattr(monitor, '_queue', lambda: ({}, None))
    before = registry.read_bytes()
    result = monitor.snapshot(tmp_path, registry=registry)
    assert set(result['counts']) == {'manifest'}
    assert result['counts']['manifest']['total'] == 4
    assert len(result['jobs']) == 2
    assert registry.read_bytes() == before


def test_monitor_uses_latest_requeue_and_ignores_old_guard_failures(tmp_path, monkeypatch):
    jobs = _campaign(tmp_path)
    directory = Path(jobs[0]['directory'])
    journal = directory / '.kinbot_slurm_io_retries/token/10'
    journal.mkdir(parents=True)
    (journal / 'execution.json').write_text('{"status":"failed","attempts":[{}]}')
    (journal / 'attempt-1').mkdir()
    (journal / 'attempt-1/stderr.txt').write_text('Input/output error OLD\n')
    latest = journal / 'restart-2'
    latest.mkdir()
    (latest / 'execution.json').write_text(
        '{"status":"running","slurm_restart_count":2,"attempts":[{}]}')
    (latest / 'attempt-2').mkdir()
    (latest / 'attempt-2/stderr.txt').write_text('Input/output error CURRENT\n')
    monkeypatch.setattr(monitor, '_queue', lambda: ({}, 'scheduler temporarily unavailable'))
    result = monitor.snapshot(tmp_path)
    row = next(row for row in result['jobs'] if row['job_id'] == '10')
    assert row['state'] == 'QUEUE-UNKNOWN'
    assert row['slurm_restarts'] == 2 and row['guard_retries'] == 1
    assert row['guard_status'] == 'running'
    assert [m['line'] for m in row['log_matches']] == ['Input/output error CURRENT']
    assert result['warnings'][0].startswith('SCHEDULER UNAVAILABLE')


def test_monitor_does_not_attach_a_reused_slurm_id_to_wrong_directory(tmp_path, monkeypatch):
    _campaign(tmp_path)
    monkeypatch.setattr(monitor, '_queue', lambda: ({'10': {
        'state': 'RUNNING', 'directory': '/another/run', 'nodes': 'other-node'}}, None))
    row = next(r for r in monitor.snapshot(tmp_path)['jobs'] if r['job_id'] == '10')
    assert not row['scheduler']
    assert any('identity mismatch' in a for a in row['alerts'])


def test_monitor_surfaces_native_input_limit_before_mpi_abort(tmp_path, monkeypatch):
    jobs = _campaign(tmp_path)
    directory = Path(jobs[0]['directory'])
    (directory / 'probe.out').write_text(' PROGRAM * MULTI\n'
        ' MAXITER IS LIMITED TO999\n ? Error\n ? Illegal input value\n'
        ' ? The problem occurs in muinp1\n GLOBAL ERROR fehler on processor 0\n')
    monkeypatch.setattr(monitor, '_queue', lambda: ({}, None))
    row = next(r for r in monitor.snapshot(tmp_path)['jobs'] if r['job_id'] == '10')
    assert row['log_matches'][0]['line'] == 'MAXITER IS LIMITED TO999'
    assert 'Illegal input value' in monitor.render({'jobs':[row], **{
        k:v for k,v in monitor.snapshot(tmp_path).items() if k != 'jobs'}}, details=True)


def test_full_root_progress_ignores_restart_values_and_survives_long_ci_output(tmp_path, monkeypatch):
    jobs = _campaign(tmp_path)
    manifest_path = tmp_path/'reconnaissance.json'
    manifest = json.loads(manifest_path.read_text())
    manifest['jobs'][0].update(correlation_checks=[{'method':'mrci-f12+q(5,4)'}],
        root_states=[{'multiplicity':m,'root':r} for m in (2,4) for r in (1,2)])
    manifest_path.write_text(json.dumps(manifest))
    directory = Path(jobs[0]['directory'])
    (directory/'probe.out').write_text(
        'reconc1=energy(1)\n PROGRAM * RESTART\n RECONC4 = -99.0 AU\n'
        ' PROGRAM * MULTI\n PROGRAM * MRCI-F12\n RECONC1 = -1.2345 AU\n'
        ' PROGRAM * MRCI-F12\n RECONC2 = -1.2344D+00 AU\n'
        ' PROGRAM * MRCI-F12\n MAXDAV=20 MAXVI=60\n'+(' CI working\n'*15000)+
        ' 8 1 1 1.15 -0.48 -835.14 0.00000000\n'
        ' CPU time for F12 integral evaluation (CABS singles) 1.62 seconds\n')
    monkeypatch.setattr(monitor, '_queue', lambda: ({}, None))
    before = {p:p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    data = monitor.snapshot(tmp_path, roots=True)
    row = next(r for r in data['jobs'] if r['job_id']=='10')
    assert row['tail_truncated']
    assert row['stage']=='MRCI-F12'
    assert row['solver_parameters']['MAXDAV']=='20'
    progress = row['root_progress']
    assert progress['expected']==4 and progress['printed']==2
    assert progress['calls_started']==3
    assert progress['current']['multiplicity']==4 and progress['current']['root']==1
    assert progress['last_selected']['energy_hartree']==pytest.approx(-1.2344)
    rendered = monitor.render(data, details=True)
    assert 'multiplicity=4, root=1' in rendered
    assert 'Correlated energy markers: 2/4' in rendered
    assert 'not validated results' in rendered
    assert {p:p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}==before


def test_root_progress_does_not_accept_precabs_energy_or_reuse_previous_multi(tmp_path):
    job={'name':'probe','correlation_checks':[{'method':'mrci-f12+q(7,5)'}],
         'root_states':[{'multiplicity':2,'root':1}]}
    p=tmp_path/'probe.out'
    p.write_text(' PROGRAM * MULTI\n PROGRAM * MRCI-F12\n RECONC1 = -1.0 AU\n'
                 ' PROGRAM * MULTI\n PROGRAM * MRCI-F12\n'
                 ' !CI-F12 STATE 1.1 Energy -835.14\n'
                 ' Cluster corrected energies -835.21 (Davidson, relaxed reference)\n'
                 ' ? strange..\n')
    result=monitor._native_progress(tmp_path,job)
    assert result['printed']==0 and result['last_selected'] is None
    assert result['current']['root']==1 and result['calls_started']==1
    p.write_text(p.read_text()+' RECONC1 = -835.22 AU\n')
    result=monitor._native_progress(tmp_path,job)
    assert result['printed']==1 and result['current'] is None
    # Markers are display evidence, even if an earlier error exists. Only the
    # controller validates completion; the monitor must not relabel a job.
    p.write_text(p.read_text()+' PROGRAM * RESTART\n RECONC1 = -835.22 AU\n')
    result=monitor._native_progress(tmp_path,job)
    assert result['printed']==0 and result['calls_started']==0


def test_all_view_lists_completed_roots_and_blocked_predecessor_without_reading_completed_output(tmp_path, monkeypatch):
    _campaign(tmp_path)
    p=tmp_path/'reconnaissance.json'
    m=json.loads(p.read_text())
    for i,j in enumerate(m['jobs']): j['name']='point'+str(i)
    m['jobs'][0]['terminal_failure']='failed native solve'
    m['jobs'][2]['correlated_energies_hartree']=[{'energy_hartree':-1.0}]*6
    m['jobs'][3]['orbital_predecessor']='point0'
    p.write_text(json.dumps(m))
    monkeypatch.setattr(monitor,'_queue',lambda: ({},None))
    def scan(directory, job):
        assert job['status']!='complete'
        return None
    monkeypatch.setattr(monitor,'_native_progress',scan)
    data=monitor.snapshot(tmp_path,include_all=True,roots=True)
    assert len(data['jobs'])==4
    assert data['chains'][0]['counts']['VALIDATED']==1
    rendered=monitor.render(data,details=True)
    assert '1/4 validated' in rendered
    assert '6 correlated roots recorded' in rendered
    assert 'predecessor=point0 (FAILED)' in rendered
