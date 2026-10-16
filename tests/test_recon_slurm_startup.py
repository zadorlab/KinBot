"""Recover Slurm's pre-shell I/O failure without hiding chemistry failures."""
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import molpro_reconnaissance as recon
from kinbot.execution import atomic_json, file_sha256, slurm_job_record


def _manifest(root, identifier='71'):
    directory = root / 'probe'
    directory.mkdir(parents=True)
    inp = directory / 'probe.inp'
    inp.write_text('immutable chemistry\n')
    (directory / 'probe.sbatch').write_text('#!/bin/bash\n#SBATCH -n 8\ntrue\n')
    job = {'name': 'probe', 'directory': str(directory), 'status': 'incomplete',
           'input_sha256': file_sha256(inp), 'submission_status': 'submitted',
           'missing_since': 0., 'job_id': identifier}
    manifest = {'jobs': [job], 'status': 'incomplete'}
    path = atomic_json(root / 'reconnaissance.json', manifest)
    record = {'JobId': identifier, 'WorkDir': str(directory), 'JobState': 'FAILED',
              'ExitCode': '0:53', 'Reason': 'RaisedSignal:53(Real-time_signal_19)',
              'BatchHost': 'faulty-host', 'RunTime': '00:00:01'}
    return path, manifest, job, record


@pytest.mark.parametrize('reason', ['unknown', 'cancelled', 'signal9', 'wrong_job',
    'wrong_directory', 'no_node', 'native', 'output', 'execution', 'stderr', 'guard',
    'guard_started', 'live'])
def test_no_automatic_retry_without_proof_of_prestart_failure(tmp_path, monkeypatch, reason):
    path, manifest, job, record = _manifest(tmp_path)
    if reason == 'unknown':
        record = {}
    elif reason == 'cancelled':
        record['JobState'] = 'CANCELLED'
    elif reason == 'signal9':
        record['ExitCode'] = '0:9'
    elif reason == 'wrong_job':
        record['JobId'] = '999'
    elif reason == 'wrong_directory':
        record['WorkDir'] = '/another/campaign'
    elif reason == 'no_node':
        record['BatchHost'] = '(null)'
    elif reason == 'guard_started':
        (Path(job['directory']) / '.kinbot_slurm_io_retries/token/71').mkdir(parents=True)
    elif reason in ('native', 'output', 'execution', 'stderr', 'guard'):
        filename = {'native': 'probe.native.inp', 'output': 'probe.out',
                    'execution': 'execution.json', 'stderr': 'probe.stderr',
                    'guard': '.kinbot_slurm_io_retries/token/71/execution.json'}[reason]
        artifact = Path(job['directory']) / filename
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text('native or shell has started')
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: record)
    live = {'71'} if reason == 'live' else set()
    recon._recover_slurm_startup_failures(path, manifest, live)
    assert job['job_id'] == '71'
    assert not job.get('startup_retry_history')
    assert not manifest.get('slurm_excluded_nodes')


def test_startup_retry_journals_evidence_and_survives_controller_restart(tmp_path, monkeypatch):
    path, manifest, job, record = _manifest(tmp_path)
    job.pop('missing_since')
    now = [100.]
    monkeypatch.setattr(recon.time, 'time', lambda: now[0])
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: record.copy())
    recon._recover_slurm_startup_failures(path, manifest, set())
    assert job['job_id'] == '71'  # Publication grace; save scheduler evidence now.
    assert job['scheduler_record'] == record
    # A new controller can recover even after scontrol purges the old record.
    manifest = json.loads(path.read_text())
    job = manifest['jobs'][0]
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: {})
    now[0] = 161.
    recon._recover_slurm_startup_failures(path, manifest, set())
    assert 'job_id' not in job
    assert job['submission_status'] == 'retry_ready'
    assert manifest['slurm_excluded_nodes'] == ['faulty-host']
    assert job['startup_retry_history'][0]['scheduler_record'] == record
    assert json.loads(path.read_text()) == manifest
    assert (Path(job['directory']) / 'probe.inp').read_text() == 'immutable chemistry\n'


def test_startup_retries_are_bounded_and_do_not_lose_final_failure(tmp_path, monkeypatch):
    path, manifest, job, record = _manifest(tmp_path)
    monkeypatch.setattr(recon, 'slurm_job_record', lambda identifier: dict(record, JobId=identifier))
    for attempt in range(3):
        job.update(job_id=str(71 + attempt), missing_since=0.)
        recon._recover_slurm_startup_failures(path, manifest, set())
    assert len(job['startup_retry_history']) == 2
    assert job['job_id'] == '73'
    with pytest.raises(recon.ProbeFailure, match=r'retry limit \(2\) exhausted'):
        recon._finished(path, manifest, set())


def test_requeued_job_is_reserved_even_before_squeue_sees_it(tmp_path, monkeypatch):
    path, manifest, job, record = _manifest(tmp_path)
    record.update(JobState='PENDING', Reason='BeginTime', ExitCode='0:0')
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: record)
    live = set()
    recon._recover_slurm_startup_failures(path, manifest, live)
    assert live == {'71'} and 'missing_since' not in job
    assert not recon._finished(path, manifest, live)
    assert job['job_id'] == '71'


def test_startup_retry_refuses_changed_input(tmp_path, monkeypatch):
    path, manifest, job, record = _manifest(tmp_path)
    (Path(job['directory']) / 'probe.inp').write_text('changed chemistry')
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: record)
    with pytest.raises(ValueError, match='input changed'):
        recon._recover_slurm_startup_failures(path, manifest, set())
    assert job['job_id'] == '71'


def test_cached_active_record_cannot_keep_a_missing_job_live_forever(tmp_path, monkeypatch):
    path, manifest, job, record = _manifest(tmp_path)
    job['scheduler_record'] = dict(record, JobState='RUNNING', ExitCode='0:0')
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: {})
    live = set()
    recon._recover_slurm_startup_failures(path, manifest, live)
    assert not live and job['job_id'] == '71'
    with pytest.raises(recon.ProbeFailure):
        recon._finished(path, manifest, live)


@pytest.mark.parametrize('state', ['OUT_OF_MEMORY', 'CANCELLED', 'PENDING'])
def test_started_job_retains_scheduler_evidence_without_unproven_retry(tmp_path, monkeypatch, state):
    path, manifest, job, record = _manifest(tmp_path)
    execution = Path(job['directory']) / 'execution.json'
    execution.write_text('{"status":"running"}')
    record.update(JobState=state, ExitCode='0:9')
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: record)
    live = set()
    recon._recover_slurm_startup_failures(path, manifest, live)
    assert job['scheduler_record'] == record
    assert job['job_id'] == '71' and not job.get('startup_retry_history')
    assert execution.read_text() == '{"status":"running"}'
    if state == 'PENDING':
        assert live == {'71'}
    else:
        with pytest.raises(recon.ProbeFailure, match=state):
            recon._finished(path, manifest, live)


@pytest.mark.parametrize('kind,expected', [
    ('interrupted', True), ('io', True), ('chemistry', False),
    ('cancelled', False), ('oom', False), ('live', False), ('grace', False)])
def test_allocation_recovery_is_bounded_and_preserves_native_evidence(tmp_path, monkeypatch, kind, expected):
    path, manifest, job, _ = _manifest(tmp_path)
    directory = Path(job['directory'])
    output = 'interrupted during integral setup\n'
    status = 'running'
    if kind == 'io':
        status, output = 'failed', 'ERROR WRITING 32768 WORDS AT OFFSET 1. TO FILE21\n'
    elif kind == 'chemistry':
        output = 'NO CONVERGENCE! MAXIMAL NUMBER OF ITERATIONS REACHED!\n'
    elif kind in ('cancelled', 'oom'):
        job['scheduler_record'] = {'JobState': 'CANCELLED' if kind == 'cancelled' else 'OUT_OF_MEMORY'}
    elif kind == 'grace':
        job['missing_since'] = 100.
    (directory / 'execution.json').write_text(json.dumps({'status': status}))
    (directory / 'probe.out').write_text(output)
    atomic_json(path, manifest)
    monkeypatch.setattr(recon.time, 'time', lambda: 120.)
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: {})
    result = recon._retry_allocation_failures(path, manifest, {'71'} if kind == 'live' else set(), 1)
    job = result['jobs'][0]
    if expected:
        assert 'job_id' not in job
        assert job['submission_status'] == 'retry_ready'
        assert (directory / 'attempts/1/probe.out').read_text() == output
        assert (directory / 'probe.inp').read_text() == 'immutable chemistry\n'
        assert len(job['allocation_retry_history']) == 1
        # Resuming a controller does not replenish the retry allowance.
        job.update(job_id='72', missing_since=0.)
        (directory / 'execution.json').write_text('{"status":"running"}')
        atomic_json(path, result)
        after = recon._retry_allocation_failures(path, result, set(), 1)
        assert after['jobs'][0]['job_id'] == '72'
        assert (directory / 'execution.json').exists()
    else:
        assert job['job_id'] == '71' and not job.get('allocation_retry_history')


@pytest.mark.parametrize('case', ['confirmed', 'same_host', 'unknown_host', 'no_requeue', 'guard_started', 'cancelled'])
def test_requeued_startup_io_preserves_old_native_files_and_retries_only_confirmed_case(tmp_path, monkeypatch, case):
    path, manifest, job, scheduler = _manifest(tmp_path)
    directory = Path(job['directory'])
    scheduler.update(Restarts='1', RunTime='00:00:00')
    execution = {'status': 'running', 'host': 'previous-host'}
    if case == 'same_host':
        execution['host'] = scheduler['BatchHost']
    elif case == 'unknown_host':
        execution.pop('host')
    elif case == 'no_requeue':
        scheduler['Restarts'] = '0'
    elif case == 'cancelled':
        scheduler['JobState'] = 'CANCELLED'
    elif case == 'guard_started':
        (directory / '.kinbot_slurm_io_retries/token/71/restart-1').mkdir(parents=True)
    job['scheduler_record'] = scheduler
    (directory / 'execution.json').write_text(json.dumps(execution))
    output = 'Program * Restricted Hartree-Fock\nNO CONVERGENCE!\nPROGRAM * MULTI\nCONVERGENCE REACHED!\nPROGRAM * RS2-F12\n'
    (directory / 'probe.out').write_text(output)
    atomic_json(path, manifest)
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: {})
    result = recon._retry_allocation_failures(path, manifest, set(), 2)
    if case == 'confirmed':
        assert result['jobs'][0]['allocation_retry_history'][0]['reason'] == 'requeued_startup_io'
        assert (directory / 'attempts/1/probe.out').read_text() == output
        assert result['slurm_excluded_nodes'] == ['faulty-host']
        assert 'job_id' not in result['jobs'][0]
    else:
        assert result['jobs'][0]['job_id'] == '71'
        assert not (directory / 'attempts').exists()
        assert not (directory / 'attempts').exists()


def test_retry_selected_interruption_preserves_a_running_peer(tmp_path, monkeypatch):
    from copy import deepcopy
    path, manifest, job, _ = _manifest(tmp_path)
    (Path(job['directory']) / 'execution.json').write_text('{"status":"running"}')
    peer = deepcopy(job)
    peer.update(name='peer', directory=str(tmp_path / 'peer'), job_id='72')
    Path(peer['directory']).mkdir()
    artifact = Path(peer['directory']) / 'execution.json'
    artifact.write_text('{"status":"running","owner":"peer"}')
    manifest['jobs'].append(peer)
    atomic_json(path, manifest)
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: {})
    result = recon._retry_failed(path, include_interrupted=True, job_ids={'71'}, live_job_ids={'72'})
    assert result['retry_ready'] == ['probe']
    assert artifact.read_text() == '{"status":"running","owner":"peer"}'
    assert json.loads(path.read_text())['jobs'][1]['job_id'] == '72'


@pytest.mark.parametrize('kind', ['io', 'interrupted'])
def test_controller_recovers_allocation_without_aborting_independent_work(tmp_path, monkeypatch, kind):
    paths = []
    for index in range(3):
        path, manifest, job, _ = _manifest(tmp_path / str(index), str(71 + index))
        paths.append(path)
        if index == 0:
            job['missing_since'] = 100. if kind == 'io' else 0.
            directory = Path(job['directory'])
            (directory / 'execution.json').write_text(json.dumps({
                'status': 'failed' if kind == 'io' else 'running'}))
            (directory / 'probe.out').write_text('Input/output error' if kind == 'io' else 'integral setup')
        elif index == 2:
            job.pop('job_id')
        atomic_json(path, manifest)
    live, submitted = {'72'}, []
    monkeypatch.setattr(recon.time, 'time', lambda: 120.)
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: {})
    monkeypatch.setattr(recon, '_audit', lambda path, **kw: json.loads(Path(path).read_text()))
    def scheduler(command, **kwargs):
        if command[0] == 'squeue':
            return SimpleNamespace(stdout='\n'.join(live), returncode=0)
        assert command[0] == 'sbatch'
        identifier = str(100 + len(submitted))
        submitted.append(kwargs['cwd'])
        live.add(identifier)
        assert len(live) <= 2
        return SimpleNamespace(stdout=identifier, stderr='', returncode=0)
    def complete_tick(_):
        for path in paths:
            manifest = json.loads(path.read_text())
            job = manifest['jobs'][0]
            if job.get('job_id') in live:
                job['status'] = 'complete'
                manifest['status'] = 'complete_review_required'
                atomic_json(path, manifest)
        live.clear()
    monkeypatch.setattr(recon.subprocess, 'run', scheduler)
    monkeypatch.setattr(recon.time, 'sleep', complete_tick)
    reports = recon.run_many(paths, max_jobs=2, startup_checked=True, allocation_retries=2)
    assert len(submitted) == 2 and len(set(submitted)) == 2
    assert reports[1]['jobs'][0]['job_id'] == '72'
    assert reports[0]['jobs'][0]['allocation_retry_history'][0]['reason'] == (
        'native_io' if kind == 'io' else 'interrupted_allocation')
    assert (paths[0].parent / 'probe/attempts/1/execution.json').exists()


def test_exclusions_preserve_site_hostlists_without_changing_parent_environment(monkeypatch):
    monkeypatch.setenv('SBATCH_EXCLUDE', 'site[01-03]')
    options = recon._slurm_submission_options(
        '#!/bin/bash\n#SBATCH -x original\n#SBATCH --exclude=other[04-06]\ntrue\n',
        ['faulty-host'])
    assert options == ['--exclude=site[01-03],faulty-host,original,other[04-06]']
    assert recon.os.environ['SBATCH_EXCLUDE'] == 'site[01-03]'


@pytest.mark.parametrize('limit', [2, 20])
def test_shared_controller_retries_prestart_failure_adopts_live_jobs_and_honors_cap(
        tmp_path, monkeypatch, limit):
    paths, records = [], {}
    for index in range(limit + 2):
        path, manifest, job, record = _manifest(tmp_path / str(index), str(71 + index))
        paths.append(path)
        if index == 0:
            records[job['job_id']] = record
        elif index > 1:
            job.pop('job_id')
            atomic_json(path, manifest)
    live, submitted, tick = {'72'}, [], [0]
    protected = {p: p.read_bytes() for p in tmp_path.rglob('*.inp')}
    monkeypatch.setattr(recon, 'slurm_job_record', lambda identifier: records.get(identifier, {}))
    monkeypatch.setattr(recon, '_audit', lambda path, **kw: json.loads(Path(path).read_text()))

    def scheduler(command, **kwargs):
        if command[0] == 'squeue':
            return SimpleNamespace(stdout='\n'.join(live), returncode=0)
        assert command[0] == 'sbatch'
        assert '--exclude=faulty-host' in command
        identifier = str(100 + len(submitted))
        submitted.append(kwargs['cwd'])
        live.add(identifier)
        assert len(live) <= limit
        return SimpleNamespace(stdout=identifier, stderr='', returncode=0)

    def complete_tick(_):
        tick[0] += 1
        assert tick[0] < 5
        for path in paths:
            manifest = json.loads(path.read_text())
            job = manifest['jobs'][0]
            if job.get('job_id') in live:
                job['status'] = 'complete'
                manifest['status'] = 'complete_review_required'
                atomic_json(path, manifest)
        live.clear()

    monkeypatch.setattr(recon.subprocess, 'run', scheduler)
    monkeypatch.setattr(recon.time, 'sleep', complete_tick)
    reports = recon.run_many(paths, max_jobs=limit, startup_checked=True)
    assert reports[1]['jobs'][0]['job_id'] == '72'  # Healthy job was adopted.
    assert len(submitted) == len(paths) - 1 and len(set(submitted)) == len(submitted)
    assert len(reports[0]['jobs'][0]['startup_retry_history']) == 1
    assert all(m['slurm_excluded_nodes'] == ['faulty-host'] for m in reports)
    assert all(p.read_bytes() == original for p, original in protected.items())


def test_scontrol_evidence_needs_no_accounting(monkeypatch):
    calls = []
    def scheduler(command, **kwargs):
        calls.append(command)
        assert kwargs['timeout'] == 15
        return SimpleNamespace(returncode=0, stdout=(
            'JobId=71 JobState=FAILED Reason=RaisedSignal:53(Real-time_signal_19) '
            'ExitCode=0:53 WorkDir=/campaign/job BatchHost=faulty-host'))
    monkeypatch.setattr(recon.subprocess, 'run', scheduler)
    assert slurm_job_record('71')['ExitCode'] == '0:53'
    assert calls == [['scontrol', 'show', 'job', '-o', '71']]


def test_documented_incident_repair_preserves_live_and_started_jobs(tmp_path, monkeypatch):
    import sys
    from kinbot import execution
    records, systems, paths = {}, [], []
    for index, identifier in enumerate(range(52907719, 52907726)):
        path, manifest, job, record = _manifest(tmp_path / str(index), str(identifier))
        paths.append(path)
        systems.append({'pilot_manifest': str(path), 'manifest': str(path)})
        if index == 0:
            record = {}  # Retained scheduler record has expired: explicit incident retry.
        elif index == 3:
            (Path(job['directory']) / 'probe.native.inp').write_text('worker started')
        elif index == 4:
            record['JobState'] = 'CANCELLED'
        elif index == 5:
            record['WorkDir'] = '/another/campaign'
        elif index == 6:
            job['status'] = 'complete'
            atomic_json(path, manifest)
        records[str(identifier)] = record
    atomic_json(tmp_path / 'campaign.json', {'systems': systems})
    monkeypatch.setattr(sys, 'argv', ['-', str(tmp_path)])
    monkeypatch.setattr(execution, 'slurm_job_record', lambda jid: records[jid])
    monkeypatch.setattr(recon.subprocess, 'run', lambda *a, **kw:
                        SimpleNamespace(stdout='52907721\n', returncode=0))
    retried = []
    def retry(path, *, include_interrupted):
        assert include_interrupted
        manifest = json.loads(path.read_text())
        retried.append(manifest['jobs'][0].pop('job_id'))
        atomic_json(path, manifest)
        return {}
    monkeypatch.setattr(recon, 'retry_failed', retry)
    document = (Path(__file__).parents[1] /
                'docs/SULFUR_REMAINING_CHECKS_2026-10-15.md').read_text()
    section = document.split('## Resume Batch A after the node21 startup failures')[1]
    snippet = section.split('<<\'PY\'\n', 1)[1].split('\nPY\n', 1)[0]
    exec(compile(snippet, '<documented incident repair>', 'exec'), {})
    exec(compile(snippet, '<repeated incident repair>', 'exec'), {})
    assert retried == ['52907719', '52907720']
    assert all(json.loads(p.read_text())['slurm_excluded_nodes'] == ['node21'] for p in paths)
    assert json.loads(paths[2].read_text())['jobs'][0]['job_id'] == '52907721'


def test_explicit_numerical_retry_preserves_live_peer_and_binds_actual_native_input(tmp_path, monkeypatch):
    from copy import deepcopy
    path, manifest, job, _ = _manifest(tmp_path)
    directory = Path(job['directory'])
    inp = directory / 'probe.inp'
    original = ('***,test\nmemory,1000,m\nfile,2,probe.wfu\n'
                'gthresh,energy=1.d-10,coeff=1.d-10\n'
                '{rs2-f12,shift=0.3,ipea=0.25;state,1,3;option,maxiti=200;}\n')
    inp.write_text(original)
    job['input_sha256'] = file_sha256(inp)
    job['ppn'] = 8
    (directory / 'execution.json').write_text('{"status":"failed","error":"reference CI"}')
    peer = deepcopy(job)
    peer.update(name='peer', directory=str(tmp_path / 'peer'), job_id='72')
    Path(peer['directory']).mkdir()
    artifact = Path(peer['directory']) / 'execution.json'
    artifact.write_text('{"status":"running","owner":"peer"}')
    manifest['jobs'].append(peer)
    atomic_json(path, manifest)
    monkeypatch.setattr(recon, 'slurm_job_record', lambda _: {})
    result = recon._retry_failed(path, job_ids={'71'}, live_job_ids={'72'}, max_reference_vectors=60,
                                 max_reference_iterations=2000)
    assert result['retry_ready'] == ['probe']
    assert artifact.read_text() == '{"status":"running","owner":"peer"}'
    job = json.loads(path.read_text())['jobs'][0]
    assert job['retry_max_reference_vectors'] == 60
    assert inp.read_text() == original
    assert json.loads((directory / 'attempts/1/execution.json').read_text())['status'] == 'failed'
    monkeypatch.setattr(recon, 'qc_runtime_environment', lambda *a, **kw: ({}, {'scratch_directory': str(tmp_path)}))
    monkeypatch.setattr(recon, 'cleanup_qc_runtime', lambda *a: None)
    def native(command, **kwargs):
        deck = (directory / command[-1]).read_text()
        assert 'gthresh,energy=1.d-10,coeff=1.d-10' in deck
        assert '{rs2-f12,shift=0.3,ipea=0.25;state,1,3;option,maxiti=2000;option,maxvi=60;}' in deck
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(recon.subprocess, 'run', native)
    recon._execute(job)
    execution = json.loads((directory / 'execution.json').read_text())
    assert execution['retry_max_reference_vectors'] == 60
    assert execution['retry_max_reference_iterations'] == 2000
    assert recon._checkpoint_evidence(directory, job, execution)[1] is None
    job['retry_max_reference_vectors'] = 61
    assert 'override differs' in recon._checkpoint_evidence(directory, job, execution)[1]
    job['retry_max_reference_vectors'] = 60
    # Even a rehashed native file must actually contain the declared override.
    native_path = Path(execution['native_input'])
    native_path.write_text(native_path.read_text().replace('maxvi=60', 'maxvi=9'))
    execution['native_input_sha256'] = file_sha256(native_path)
    assert 'does not contain' in recon._checkpoint_evidence(directory, job, execution)[1]
    assert inp.read_text() == original


@pytest.mark.parametrize('value', [0, -1, 1.5, True])
def test_numerical_retry_rejects_invalid_capacity_before_touching_files(tmp_path, value):
    with pytest.raises(ValueError):
        recon._retry_failed(tmp_path / 'missing.json', job_ids={'1'}, max_reference_vectors=value)
    with pytest.raises(ValueError):
        recon._retry_failed(tmp_path / 'missing.json', job_ids={'1'}, max_reference_iterations=value)


def test_reference_capacity_override_cannot_change_an_unrelated_program():
    text = '{multi;state,3;maxit,150;}\n{mrci-f12;state,1,3;option,maxvi=9,cluster=0;}\n'
    expected = text.replace('maxvi=9', 'maxvi=60')
    assert recon._reference_solver_override(text, 60) == expected
    assert recon._reference_solver_override(expected, 60) == expected
    with pytest.raises(ValueError, match='No supported'):
        recon._reference_solver_override('{uccsd(t)-f12;}', 60)
    with pytest.raises(ValueError, match='explicit MAXDAV'):
        recon._reference_solver_override('{rs2;maxdav,4,9;}', 60)
    with pytest.raises(ValueError, match='explicit --job-id'):
        recon._retry_failed(Path('/missing'), max_reference_vectors=60)
    with pytest.raises(ValueError, match='explicit --job-id'):
        recon._retry_failed(Path('/missing'), max_reference_iterations=2000)
    with pytest.raises(ValueError, match='explicit MAXITER'):
        recon._reference_solver_override('{rs2;maxiter,30,200;}', iterations=2000)
    assert recon._reference_solver_override(text, iterations=2000).startswith('{multi;state,3;maxit,150;}')


def test_manifest_solver_settings_reach_submission_without_mutating_live_work(tmp_path, monkeypatch):
    from copy import deepcopy
    path, manifest, job, _ = _manifest(tmp_path)
    job.pop('job_id')
    job.pop('missing_since', None)
    job['correlation_checks'] = [{'method': 'caspt2-f12(5,4)'}]
    deck = Path(job['directory']) / 'probe.inp'
    deck.write_text('{multi;maxit,150;}\n{rs2-f12;state,1,3;option,maxiti=200;}\n')
    job['input_sha256'] = file_sha256(deck)
    peer = deepcopy(job)
    peer.update(name='live', job_id='72', directory=str(tmp_path / 'live'))
    manifest['jobs'].append(peer)
    manifest['reference_solver_options'] = {'max_reference_vectors': 60, 'max_reference_iterations': 2000}
    atomic_json(path, manifest)
    monkeypatch.setattr(recon, '_audit', lambda p, **kw: json.loads(Path(p).read_text()))
    commands = []
    def submit(command, **kwargs):
        assert command[0] == 'sbatch'
        commands.append(command)
        return SimpleNamespace(returncode=0, stdout='99', stderr='')
    monkeypatch.setattr(recon.subprocess, 'run', submit)
    recon._submit(path, max_jobs=20, live_job_ids={'72'}, startup_checked=True)
    saved = json.loads(path.read_text())['jobs']
    assert len(commands) == 1 and saved[0]['job_id'] == '99'
    assert saved[0]['retry_max_reference_iterations'] == 2000
    assert saved[0]['retry_max_reference_vectors'] == 60
    assert saved[1] == peer
    assert file_sha256(deck) == job['input_sha256']


def test_campaign_explicit_recovery_covers_all_exhausted_startups_without_touching_qc_or_live(tmp_path, monkeypatch):
    from kinbot import orbital_campaign as campaign
    systems, paths, records = [], [], {}
    for index in range(25):
        path, manifest, job, record = _manifest(tmp_path / str(index), str(100 + index))
        job['startup_retry_history'] = [{'job_id': 'old-1'}, {'job_id': 'old-2'}]
        directory = Path(job['directory'])
        if index in (21, 22):
            (directory / 'execution.json').write_text('{"status":"failed","error":"reference CI"}')
            (directory / 'probe.out').write_text('NO CONVERGENCE IN REFERENCE CI AFTER200 ITERATIONS')
        elif index == 23:
            record['JobState'] = 'RUNNING'
            (directory / 'execution.json').write_text('{"status":"running"}')
        elif index == 24:
            record['JobState'] = 'CANCELLED'
        records[str(100 + index)] = record
        atomic_json(path, manifest)
        paths.append(path)
        systems.append({'pilot_manifest': str(path), 'manifest': str(path)})
    atomic_json(tmp_path / 'campaign.json', {'systems': systems})
    monkeypatch.setattr(campaign, '_validate_inputs', lambda record: None)
    monkeypatch.setattr(recon, 'slurm_job_record', lambda identifier: records.get(identifier, {}))
    monkeypatch.setattr(campaign.subprocess, 'run', lambda *a, **kw:
                        SimpleNamespace(stdout='123\n', returncode=0))
    protected = {p: p.read_bytes() for i in (21, 22, 23) for p in paths[i].parent.rglob('*')
                 if p.is_file() and p != paths[i]}
    result = campaign.retry_startup_failures(tmp_path)
    assert len(result['manual_retries']) == 21
    for path in paths[:21]:
        job = json.loads(path.read_text())['jobs'][0]
        assert 'job_id' not in job and job['submission_status'] == 'retry_ready'
        assert len(job['startup_retry_history']) == 2 and job['startup_retry_limit'] == 4
        assert len(job['manual_startup_retry_history']) == 1
        assert (path.parent / 'probe/attempts/1/probe.inp').read_text() == 'immutable chemistry\n'
    for i in (21, 22, 23, 24):
        assert json.loads(paths[i].read_text())['jobs'][0]['job_id'] == str(100+i)
    assert all(p.read_bytes() == original for p, original in protected.items())
    # An unchanged second invocation cannot reset budgets or retry these jobs twice.
    assert campaign.retry_startup_failures(tmp_path)['manual_retries'] == []
    path = paths[0]
    manifest = json.loads(path.read_text())
    job = manifest['jobs'][0]
    assert job['startup_retry_limit'] == 4 and len(job['manual_startup_retry_history']) == 1
    for attempt in range(3):
        job.update(job_id=str(200+attempt), missing_since=0.)
        records[str(200+attempt)] = dict(records['100'], JobId=str(200+attempt))
        recon._recover_slurm_startup_failures(path, manifest, set())
    assert len(job['startup_retry_history']) == 4 and job['job_id'] == '202'
    assert 'retry limit (4) exhausted' in job['startup_failure']
