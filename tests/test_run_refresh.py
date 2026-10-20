"""Controller handoff preserves owned native work and rejects unsafe pauses."""
from copy import deepcopy
import ctypes
import errno
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
from types import SimpleNamespace

import pytest

from kinbot import run_cleanup, run_refresh, run_salvage
from kinbot.anl import dispatch
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task


@pytest.fixture
def live_tree(tmp_path, monkeypatch):
    root = tmp_path / 'run'; root.mkdir()
    pin = root / '.salvage_code.new'; pin.mkdir()
    inputs, processes, jobs, graphs = [], [], [], []
    for index, case in enumerate(('ch3s', 'ch3ss')):
        directory = root / case; directory.mkdir()
        path = directory / 'input.json'; path.write_text('{}')
        inputs.append(str(path))
        worker = directory / '123'; worker.mkdir()
        (worker / '.kinbot_run.json').write_text('{}')
        (worker / 'kinbot.db').write_bytes(b'accepted L1/L2 results')
        for folder, module, name in ((directory, 'pes', 'pes_status.json'),
                                      (worker, 'kb', 'kinbot_status.json')):
            pid = 100 + len(processes)
            processes.append(dict(pid=pid, ppid=1, starttime=str(pid * 10),
                                  cwd=str(folder), state='S',
                                  command=f'/python -m kinbot.{module} input.json'))
            atomic_json(folder / name, dict(pid=pid, host=socket.gethostname(), status='running'))
        atomic_json(directory / 'launch.json', dict(pid=processes[-2]['pid']))
        (directory / 'driver.log').write_text('original driver log')
        for branch in ('anl_automatic/states/test/composite', 'vrctst/molpro_scan/scan'):
            graph = worker / branch
            graph.parent.mkdir(parents=True)
            spec = dispatch_spec()
            spec['tasks'] = [molpro_task('energy', 'rhf\n', geometry_from='initial')]
            dispatch.prepare(spec, graph)
            state = dispatch._load(graph)[2]
            job = str(500 + len(jobs))
            state['tasks']['energy'].update(status='submitted', job_id=job)
            atomic_json(graph / 'state.json', state)
            jobs.append(dict(job_id=job, directory=str(graph / 'tasks/energy'), state='RUNNING'))
            graphs.append(graph)
    atomic_json(root / 'prepared.json', dict(inputs=inputs, sha256={
        str(Path(p).relative_to(root)): file_sha256(p) for p in inputs}))
    atomic_json(root / 'salvage.json', dict(status='launched', inputs=inputs,
                                         code=str(root / '.salvage_code.old')))
    monkeypatch.setattr(run_salvage, '_pin_code', lambda *a, **kw: pin)
    monkeypatch.setattr(run_refresh, '_preflight', lambda root: None)
    monkeypatch.setattr(run_salvage, '_preflight_worker_routes', lambda root: None)
    monkeypatch.setattr('kinbot.anl.automatic.validate_automatic_environment', lambda par: None)
    monkeypatch.setattr('kinbot.parameters.Parameters', lambda path: SimpleNamespace(par={}))
    monkeypatch.setattr(run_cleanup, '_jobs', lambda root: deepcopy(jobs))
    monkeypatch.setattr(run_cleanup, '_processes', lambda folder: [deepcopy(p) for p in processes
                        if Path(p['cwd']).is_relative_to(folder)])
    signals = []
    def send(fd, sig):
        if fd == os.getpid() and sig == 0:
            return
        record = next((p for p in processes if p['pid'] == fd), None)
        if record is None:
            raise ProcessLookupError
        if sig == 0:
            return
        signals.append((fd, sig))
        if sig == signal.SIGKILL:
            processes.remove(record)
        else:
            record['state'] = 'T' if sig == signal.SIGSTOP else 'S'
    monkeypatch.setattr(run_refresh.os, 'pidfd_open', lambda pid: pid, raising=False)
    monkeypatch.setattr(run_refresh.signal, 'pidfd_send_signal', send, raising=False)
    monkeypatch.setattr(run_refresh.os, 'close', lambda fd: None)
    monkeypatch.setattr(run_refresh, '_starttime', lambda pid: str(pid * 10))
    monkeypatch.setattr(run_refresh, '_stopped', lambda pid: all(p['state'] == 'T' for p in processes if p['pid'] == pid))
    launches = []
    def launch(command, **kwargs):
        launches.append((command, kwargs))
        return SimpleNamespace(pid=900 + len(launches))
    monkeypatch.setattr(run_salvage.subprocess, 'Popen', launch)
    return SimpleNamespace(root=root, processes=processes, jobs=jobs, graphs=graphs,
                           signals=signals, launches=launches, send=send)


def scientific_files(root):
    return {str(p.relative_to(root)): file_sha256(p) for p in root.rglob('*') if p.is_file()
            and ('composite' in p.parts or 'molpro_scan' in p.parts or p.name == 'kinbot.db')}


def test_refresh_keeps_live_anl_and_vrc_and_restarts_once(live_tree):
    t = live_tree
    before = scientific_files(t.root)
    original_jobs = deepcopy(t.jobs)
    review = run_refresh.refresh(t.root)
    assert review['phase'] == 'reviewed' and t.signals == []
    result = run_refresh.refresh(t.root, apply=True)
    assert result['phase'] == 'launched'
    assert t.jobs == original_jobs and scientific_files(t.root) == before
    assert len(t.launches) == 2
    assert {pid for pid, sig in t.signals if sig == signal.SIGKILL} == {100, 101, 102, 103}
    assert all(args['env']['PYTHONPATH'].startswith(result['code']) for _, args in t.launches)
    archive = Path(result['archive'])
    assert (archive / 'ch3ss/driver.log').read_text() == 'original driver log'
    assert run_refresh.refresh(t.root, apply=True) == result
    assert len(t.launches) == 2


def test_refresh_preserves_shared_policy_and_all_running_allocations(live_tree):
    from kinbot import run_policy
    t = live_tree
    run_policy.update(t.root, max_jobs=23, exclude=['reserved'])
    policy_bytes = (t.root/run_policy.NAME).read_bytes()
    before = scientific_files(t.root)
    jobs = deepcopy(t.jobs)
    result = run_refresh.refresh(t.root, apply=True)
    assert result['phase'] == 'launched' and t.jobs == jobs
    assert (t.root/run_policy.NAME).read_bytes() == policy_bytes
    assert scientific_files(t.root) == before


@pytest.fixture
def vrc_migration_tree(live_tree, monkeypatch):
    from kinbot import vrc_restart
    import shutil
    t = live_tree
    monkeypatch.setattr(vrc_restart, 'validate_input', lambda *a: None)
    t.commands = []
    t.terminal = {}
    for index, graph in enumerate(t.graphs):
        if 'molpro_scan' not in graph.parts:
            continue
        shutil.rmtree(graph)
        spec = dispatch_spec()
        spec['tasks'] = [molpro_task('relaxed_scan', 'rhf\n', geometry_from='initial')]
        dispatch.prepare(spec, graph)
        state = dispatch._load(graph)[2]
        state['tasks']['relaxed_scan'].update(status='submitted', job_id=t.jobs[index]['job_id'])
        atomic_json(graph/'state.json', state)
        t.jobs[index]['directory'] = str(graph/'tasks/relaxed_scan')
        (graph/'tasks/relaxed_scan/optimization.traj').write_bytes(b'saved optimizer progress')
    def record(ident):
        job = next((j for j in t.jobs if j['job_id'] == ident), t.terminal.get(ident))
        return {} if job is None else dict(JobId=ident, JobState=job['state'], WorkDir=job['directory'])
    def cancel(args):
        t.commands.append(args)
        assert len(args) == 2 and args[0] == 'scancel'
        job = next(j for j in t.jobs if j['job_id'] == args[1])
        t.jobs.remove(job)
        t.terminal[job['job_id']] = dict(job, state='CANCELLED')
        return ''
    t.record = record
    monkeypatch.setattr('kinbot.execution.slurm_job_record', record)
    monkeypatch.setattr(run_refresh, '_pending_request', cancel)
    return t


def test_targeted_vrc_restart_keeps_all_native_files_and_unrelated_jobs(vrc_migration_tree):
    t = vrc_migration_tree
    before = scientific_files(t.root)
    updates = {'vrc_tst_scan_basis': 'avdz-f12',
               'vrc_tst_fragment_geometry': 'asymptotic_scan', 'vrc_tst_pivot_level': 'L2',
               'vrc_tst_spin_splitting': {'reference': {
                   'basis': ['avdz-f12', 'avtz-f12'], 'extrapolation_power': 3.7}}}
    options = dict(vrc_updates=updates, cancel_superseded_vrc_scans=True)
    review = run_refresh.refresh(t.root, **options)
    assert review['phase'] == 'reviewed' and not t.commands and not t.signals
    result = run_refresh.refresh(t.root, apply=True, **options)
    assert result['phase'] == 'launched'
    assert t.commands == [['scancel', '501'], ['scancel', '503']]
    assert {j['job_id'] for j in t.jobs} == {'500', '502'}
    assert before == scientific_files(t.root)
    assert len(t.launches) == 2
    _, inputs = run_salvage._cases(t.root)
    for path in inputs:
        assert json.loads(path.read_text()) == updates
    saved = json.loads((t.root/'vrc_restart.json').read_text())
    assert len(saved['retired_scans']) == 2
    assert run_refresh.refresh(t.root, apply=True, **options) == result
    assert len(t.commands) == 2 and len(t.launches) == 2


def test_vrc_cancel_unknown_outcome_never_launches_duplicates(vrc_migration_tree, monkeypatch):
    t = vrc_migration_tree
    options = dict(vrc_updates={'vrc_tst_scan_basis': 'avdz-f12'}, cancel_superseded_vrc_scans=True)
    monkeypatch.setattr('kinbot.execution.slurm_job_record',
                        lambda ident: {} if ident in t.terminal else t.record(ident))
    with pytest.raises(RuntimeError, match='not yet confirmed'):
        run_refresh.refresh(t.root, apply=True, **options)
    assert not t.launches
    assert json.loads((t.root/'ch3s/input.json').read_text()) == {}
    monkeypatch.setattr('kinbot.execution.slurm_job_record', t.record)
    assert run_refresh.refresh(t.root, apply=True, **options)['phase'] == 'launched'
    assert t.commands == [['scancel', '501'], ['scancel', '503']]


def test_vrc_input_migration_interruption_resumes(vrc_migration_tree, monkeypatch):
    from kinbot import vrc_restart
    t = vrc_migration_tree
    options = dict(vrc_updates={'vrc_tst_scan_basis': 'avdz-f12'}, cancel_superseded_vrc_scans=True)
    original = vrc_restart.atomic_json
    once = True
    def interrupt(path, payload):
        nonlocal once
        if Path(path) == t.root/'prepared.json' and once:
            once = False
            raise OSError('simulated manifest interruption')
        return original(path, payload)
    monkeypatch.setattr(vrc_restart, 'atomic_json', interrupt)
    with pytest.raises(OSError, match='simulated manifest'):
        run_refresh.refresh(t.root, apply=True, **options)
    assert not t.launches
    assert run_refresh.refresh(t.root, apply=True, **options)['phase'] == 'launched'
    assert len(t.commands) == 2 and len(t.launches) == 2


def test_vrc_restart_does_not_touch_other_pending_allocations(vrc_migration_tree):
    t = vrc_migration_tree
    t.jobs[0]['state'] = 'PENDING'
    run_refresh.refresh(t.root, apply=True, vrc_updates={'vrc_tst_scan_basis': 'avdz-f12'},
                        cancel_superseded_vrc_scans=True)
    assert t.jobs[0]['job_id'] == '500' and t.jobs[0]['state'] == 'PENDING'


def test_reference_only_update_retains_running_geometry_jobs(vrc_migration_tree):
    from kinbot import vrc_restart
    t = vrc_migration_tree
    jobs = deepcopy(t.jobs)
    run_refresh.refresh(t.root, apply=True, vrc_updates=vrc_restart.settings(
        reference_bases=['avdz-f12', 'avtz-f12'], reference_power=3.7))
    assert t.jobs == jobs and not t.commands


@pytest.fixture
def pending_tree(legacy_pending_tree):
    return legacy_pending_tree


def test_refresh_replaces_only_confirmed_pending_and_is_idempotent(pending_tree):
    t = pending_tree
    live_before = {p: file_sha256(p) for graph in t.graphs[1:] for p in graph.rglob('*') if p.is_file()}
    result = run_refresh.refresh(t.root, apply=True, replace_pending=True)
    assert result['phase'] == 'launched' and len(t.launches) == 2
    assert result['pending_cancellation_evidence'].keys() == {'500'}
    assert {j['job_id'] for j in result['jobs']} == {'501', '502', '503'}
    state = dispatch._load(t.graphs[0])[2]
    assert 'energy' not in state['tasks'] and state['pending_retries']['energy']['attempt'] == 2
    assert Path(state['pending_retries']['energy']['previous_attempt']).is_dir()
    assert live_before == {p: file_sha256(p) for p in live_before}
    assert run_refresh.refresh(t.root, apply=True, replace_pending=True) == result
    assert len(t.launches) == 2
    assert t.commands.count(['scancel', '500']) == 1


def test_pending_disappearance_requires_cancellation_evidence(pending_tree, monkeypatch):
    t = pending_tree
    monkeypatch.setattr('kinbot.execution.slurm_job_record',
                        lambda job: {} if job in t.cancelled else t.record(job))
    with pytest.raises(RuntimeError, match='scheduler identity/state unavailable'):
        run_refresh.refresh(t.root, apply=True, replace_pending=True)
    assert not t.launches
    assert (t.graphs[0]/'tasks/energy').exists()
    monkeypatch.setattr('kinbot.execution.slurm_job_record', t.record)
    assert run_refresh.refresh(t.root, apply=True, replace_pending=True)['phase'] == 'launched'
    assert t.commands.count(['scancel', '500']) == 1


def test_pending_archive_interruption_resumes_without_reading_partial_ledger(pending_tree, monkeypatch):
    from kinbot import run_recover
    t = pending_tree
    original = run_recover.atomic_json
    stopped = False
    def interrupted(path, data):
        nonlocal stopped
        if Path(path) == t.graphs[0]/'state.json' and not stopped:
            stopped = True
            raise OSError('interrupted pending archive')
        return original(path, data)
    monkeypatch.setattr(run_recover, 'atomic_json', interrupted)
    with pytest.raises(OSError, match='interrupted pending archive'):
        run_refresh.refresh(t.root, apply=True, replace_pending=True)
    assert not t.launches
    monkeypatch.setattr(run_recover, 'atomic_json', original)
    assert run_refresh.refresh(t.root, apply=True, replace_pending=True)['phase'] == 'launched'
    assert len(t.launches) == 2
    assert t.commands.count(['scancel', '500']) == 1


def test_no_pending_cancellation_mechanism_does_not_stop_controllers(live_tree, monkeypatch):
    t = live_tree
    monkeypatch.setattr(run_refresh.subprocess, 'run', lambda *a, **kw:
                        SimpleNamespace(returncode=0, stdout='old client', stderr=''))
    with pytest.raises(RuntimeError, match='scontrol hold/release'):
        run_refresh.refresh(t.root, apply=True, replace_pending=True)
    assert not t.signals and not t.launches


@pytest.fixture
def legacy_pending_tree(live_tree, monkeypatch):
    """Older client with no --ctld; hold changes running requeue priority too."""
    t = live_tree
    t.jobs[0]['state'] = 'PENDING'
    priorities = {j['job_id']: '123' for j in t.jobs}
    t.commands, t.cancelled, t.events = [], {}, {}
    t.priorities = priorities

    def record(ident):
        job = next((j for j in t.jobs if j['job_id'] == ident), None)
        if job:
            return dict(JobId=ident, JobState=job['state'], WorkDir=job['directory'],
                        Priority=priorities[ident])
        return deepcopy(t.cancelled.get(ident, {}))

    def command(args, **kwargs):
        t.commands.append(args)
        if args == ['scancel', '--help']:
            return SimpleNamespace(returncode=0, stdout='-t, --state=STATE', stderr='')
        if args == ['scontrol', '--help']:
            return SimpleNamespace(returncode=0, stdout='', stderr='hold job_list\nrelease job_list')
        assert '--ctld' not in args
        if args[:2] in (['scontrol', 'hold'], ['scontrol', 'release']):
            action, ident = args[1:]
            if action in t.events:
                t.events[action](ident)
            priorities[ident] = '0' if action == 'hold' else '123'
        else:
            assert args[0] == 'scancel' and len(args) == 2
            ident = args[1]
            job = next(j for j in t.jobs if j['job_id'] == ident)
            # Deliberately NO client-side state protection in this simulator:
            # an unsafe bare scancel would kill RUNNING as well as PENDING.
            assert job['state'] == 'PENDING' and priorities[ident] == '0'
            t.cancelled[ident] = {**record(ident), 'JobState': 'CANCELLED'}
            t.jobs.remove(job)
            action = 'cancel'
        if action + '_after' in t.events:
            t.events[action + '_after'](ident)
        return SimpleNamespace(returncode=0, stdout='', stderr='')

    monkeypatch.setattr(run_refresh.subprocess, 'run', command)
    monkeypatch.setattr('kinbot.execution.slurm_job_record', record)
    t.command, t.record = command, record
    return t


def test_older_slurm_refresh_replaces_pending_preserves_running_and_repeats(legacy_pending_tree):
    t = legacy_pending_tree
    before = {p: file_sha256(p) for graph in t.graphs[1:] for p in graph.rglob('*') if p.is_file()}
    result = run_refresh.refresh(t.root, apply=True, replace_pending=True)
    assert result['phase'] == 'launched' and len(t.launches) == 2
    assert result['pending_cancel_method'] == 'hold_then_cancel'
    assert result['pending_cancellation_evidence'].keys() == {'500'}
    assert {j['job_id'] for j in result['jobs']} == {'501', '502', '503'}
    state = dispatch._load(t.graphs[0])[2]
    assert 'energy' not in state['tasks'] and state['pending_retries']['energy']['attempt'] == 2
    assert Path(state['pending_retries']['energy']['previous_attempt']).is_dir()
    assert before == {p: file_sha256(p) for p in before}
    assert run_refresh.refresh(t.root, apply=True, replace_pending=True) == result
    assert t.commands.count(['scontrol', 'hold', '500']) == 1
    assert t.commands.count(['scancel', '500']) == 1
    assert len(t.launches) == 2


def test_older_slurm_start_during_hold_is_released_and_retained(legacy_pending_tree):
    t = legacy_pending_tree
    before = scientific_files(t.root)
    t.events['hold'] = lambda ident: t.jobs[0].update(state='RUNNING')
    result = run_refresh.refresh(t.root, apply=True, replace_pending=True)
    assert result['phase'] == 'launched' and len(result['jobs']) == 4
    assert not t.cancelled and not result['pending_cancellation_evidence']
    assert t.priorities['500'] == '123'
    assert ['scontrol', 'release', '500'] in t.commands
    assert scientific_files(t.root) == before


@pytest.mark.parametrize('operation', ['hold', 'cancel'])
def test_older_slurm_timed_out_ack_is_reconciled_from_scheduler(legacy_pending_tree, operation):
    t = legacy_pending_tree
    def lost_ack(ident):
        raise subprocess.TimeoutExpired(operation, 30)
    t.events[operation + '_after'] = lost_ack
    result = run_refresh.refresh(t.root, apply=True, replace_pending=True)
    assert result['phase'] == 'launched' and len(t.launches) == 2
    assert result['pending_cancellation_evidence'].keys() == {'500'}
    assert t.commands.count(['scancel', '500']) == 1


def test_older_slurm_failed_hold_never_cancels_and_can_resume(legacy_pending_tree):
    t = legacy_pending_tree
    def rejected(ident):
        raise OSError('controller unavailable')
    t.events['hold'] = rejected
    with pytest.raises(RuntimeError, match='hold not confirmed'):
        run_refresh.refresh(t.root, apply=True, replace_pending=True)
    assert not t.cancelled and not t.launches
    assert (t.graphs[0]/'tasks/energy').exists()
    t.events.clear()
    assert run_refresh.refresh(t.root, apply=True, replace_pending=True)['phase'] == 'launched'
    assert t.commands.count(['scancel', '500']) == 1


def test_older_slurm_missing_record_after_hold_never_cancels(legacy_pending_tree, monkeypatch):
    t = legacy_pending_tree
    available = [True]
    monkeypatch.setattr('kinbot.execution.slurm_job_record',
                        lambda ident: t.record(ident) if available[0] else {})
    t.events['hold_after'] = lambda ident: available.__setitem__(0, False)
    with pytest.raises(RuntimeError, match='scheduler identity/state unavailable'):
        run_refresh.refresh(t.root, apply=True, replace_pending=True)
    assert not t.cancelled and not t.launches and t.priorities['500'] == '0'
    available[0] = True
    assert run_refresh.refresh(t.root, apply=True, replace_pending=True)['phase'] == 'launched'
    assert t.commands.count(['scontrol', 'hold', '500']) == 1
    assert t.commands.count(['scancel', '500']) == 1


def test_older_slurm_changed_job_identity_is_not_held(legacy_pending_tree, monkeypatch):
    t = legacy_pending_tree
    monkeypatch.setattr('kinbot.execution.slurm_job_record',
                        lambda ident: {**t.record(ident), 'WorkDir': '/another/run'})
    with pytest.raises(RuntimeError, match='scheduler identity/state unavailable'):
        run_refresh.refresh(t.root, apply=True, replace_pending=True)
    assert not any(a[:2] == ['scontrol', 'hold'] for a in t.commands)
    assert not t.cancelled and not t.launches


def test_older_slurm_archive_interruption_resumes(legacy_pending_tree, monkeypatch):
    from kinbot import run_recover
    t = legacy_pending_tree
    original = run_recover.atomic_json
    def interrupted(path, data):
        if Path(path) == t.graphs[0]/'state.json':
            raise OSError('interrupted pending archive')
        return original(path, data)
    monkeypatch.setattr(run_recover, 'atomic_json', interrupted)
    with pytest.raises(OSError, match='interrupted pending archive'):
        run_refresh.refresh(t.root, apply=True, replace_pending=True)
    monkeypatch.setattr(run_recover, 'atomic_json', original)
    assert run_refresh.refresh(t.root, apply=True, replace_pending=True)['phase'] == 'launched'
    assert len(t.launches) == 2
    assert t.commands.count(['scancel', '500']) == 1


def test_older_slurm_release_failure_is_resumable_without_cancelling_running(legacy_pending_tree):
    t = legacy_pending_tree
    t.events['hold'] = lambda ident: t.jobs[0].update(state='RUNNING')
    t.events['release'] = lambda ident: (_ for _ in ()).throw(OSError('controller unavailable'))
    with pytest.raises(RuntimeError, match='releasing refresh hold not confirmed'):
        run_refresh.refresh(t.root, apply=True, replace_pending=True)
    assert not t.cancelled and not t.launches
    t.events.clear()
    result = run_refresh.refresh(t.root, apply=True, replace_pending=True)
    assert result['phase'] == 'launched' and len(result['jobs']) == 4
    assert not t.cancelled and t.priorities['500'] == '123'


def test_older_slurm_partial_batch_keeps_cancel_evidence_and_running_work(legacy_pending_tree):
    t = legacy_pending_tree
    t.jobs[1]['state'] = 'PENDING'
    def second_hold_fails(ident):
        if ident == '501':
            raise OSError('controller unavailable')
    t.events['hold'] = second_hold_fails
    with pytest.raises(RuntimeError, match='hold not confirmed'):
        run_refresh.refresh(t.root, apply=True, replace_pending=True)
    assert t.cancelled.keys() == {'500'} and not t.launches
    t.cancelled.clear()  # Slurm purged its terminal record before retry.
    t.events.clear()
    result = run_refresh.refresh(t.root, apply=True, replace_pending=True)
    assert result['pending_cancellation_evidence'].keys() == {'500', '501'}
    assert {j['job_id'] for j in result['jobs']} == {'502', '503'}
    assert len(t.launches) == 2 and t.commands.count(['scancel', '500']) == 1


@pytest.mark.parametrize('advertises_ctld', [True, False])
def test_pending_preflight_does_not_assume_ctld_implies_server_filter(monkeypatch, advertises_ctld):
    def command(args, **kw):
        text = ('--ctld --state=STATE' if advertises_ctld else '--state=STATE')
        if args[0] == 'scontrol':
            text = 'hold job_list\nrelease job_list'
        return SimpleNamespace(returncode=0, stdout='', stderr=text)
    monkeypatch.setattr(run_refresh.subprocess, 'run', command)
    assert run_refresh._pending_preflight() == 'hold_then_cancel'


@pytest.mark.parametrize('problem', ['child', 'unowned_job', 'submitting', 'changed_input',
                                     'scheduler', 'staging', 'partial_graph'])
def test_refresh_refuses_uncertain_ownership_without_stopping(live_tree, monkeypatch, problem):
    t = live_tree
    if problem == 'child':
        t.processes.append(dict(pid=105, cwd=str(t.root / 'ch3s/123'), command='/bin/sbatch job.slurm'))
    elif problem == 'unowned_job':
        t.jobs.append(dict(job_id='999', directory=str(t.root / 'ch3s/unowned')))
    elif problem == 'scheduler':
        monkeypatch.setattr(run_cleanup, '_jobs', lambda root: (_ for _ in ()).throw(RuntimeError('outage')))
    elif problem == 'submitting':
        path = t.graphs[0] / 'state.json'
        state = json.loads(path.read_text()); state['tasks']['energy']['status'] = 'submitting'
        atomic_json(path, state)
    elif problem == 'staging':
        (t.graphs[0] / 'tasks/uncommitted').mkdir()
    elif problem == 'partial_graph':
        (t.root / 'ch3s/123/vrctst/molpro_scan/new').mkdir()
    else:
        (t.graphs[0] / 'tasks/energy/energy.inp').write_text('wrong')
    with pytest.raises((RuntimeError, ValueError)):
        run_refresh.refresh(t.root, apply=True)
    assert t.signals == [] and not t.launches


def test_submission_race_after_stop_resumes_original_controllers(live_tree, monkeypatch):
    t = live_tree
    original = run_refresh._owned_jobs
    count = 0
    def inventory(root):
        nonlocal count
        count += 1
        if count == 2:
            raise RuntimeError('submission in flight')
        return original(root)
    monkeypatch.setattr(run_refresh, '_owned_jobs', inventory)
    with pytest.raises(RuntimeError, match='submission in flight'):
        run_refresh.refresh(t.root, apply=True)
    assert all(p['state'] == 'S' for p in t.processes)
    assert not any(sig == signal.SIGKILL for _, sig in t.signals)
    assert not t.launches


def test_pid_reuse_does_not_receive_any_signal(live_tree, monkeypatch):
    monkeypatch.setattr(run_refresh, '_starttime', lambda pid: 'reused')
    with pytest.raises(RuntimeError, match='identity changed'):
        run_refresh.refresh(live_tree.root, apply=True)
    assert not live_tree.signals


def test_interrupted_archive_resumes_without_touching_qc(live_tree, monkeypatch):
    t = live_tree
    before = scientific_files(t.root)
    original = run_salvage._save
    count = 0
    def archive(*args, **kwargs):
        nonlocal count
        count += 1
        if count == 3:
            raise RuntimeError('interrupted archiving')
        return original(*args, **kwargs)
    monkeypatch.setattr(run_salvage, '_save', archive)
    with pytest.raises(RuntimeError, match='interrupted archiving'):
        run_refresh.refresh(t.root, apply=True)
    assert not t.processes and not t.launches
    monkeypatch.setattr(run_salvage, '_save', original)
    assert run_refresh.refresh(t.root, apply=True)['phase'] == 'launched'
    assert scientific_files(t.root) == before
    assert len(t.launches) == 2


@pytest.mark.parametrize('missing', [('open',), ('send',), ('open', 'send')])
def test_missing_python_pidfd_apis_use_linux_syscalls_for_full_handoff(live_tree, monkeypatch, missing):
    t = live_tree
    before = scientific_files(t.root)
    if 'open' in missing:
        monkeypatch.delattr(run_refresh.os, 'pidfd_open')
    if 'send' in missing:
        monkeypatch.delattr(run_refresh.signal, 'pidfd_send_signal')
    monkeypatch.setattr(run_refresh.sys, 'platform', 'linux')
    monkeypatch.setattr(run_refresh.platform, 'machine', lambda: 'x86_64')
    calls = []
    def syscall(number, *args):
        values = [arg.value for arg in args]
        calls.append((number.value, values))
        if number.value == 434:
            assert values[1] == 0
            return values[0]
        assert number.value == 424 and values[2:] == [None, 0]
        try:
            t.send(values[0], values[1])
        except ProcessLookupError:
            ctypes.set_errno(errno.ESRCH)
            return -1
        return 0
    monkeypatch.setattr(run_refresh.ctypes, 'CDLL', lambda *a, **kw: SimpleNamespace(syscall=syscall))
    assert run_refresh.refresh(t.root, apply=True)['phase'] == 'launched'
    assert calls and len(t.launches) == 2
    assert scientific_files(t.root) == before


def test_kernel_probe_failure_happens_before_journal_or_stop(live_tree, monkeypatch):
    t = live_tree
    def unavailable(pid):
        raise OSError(errno.ENOSYS, 'kernel lacks pidfd_open')
    monkeypatch.setattr(run_refresh, '_pidfd_open', unavailable)
    with pytest.raises(RuntimeError, match='No controllers were paused'):
        run_refresh.refresh(t.root, apply=True)
    assert not t.signals and not t.launches
    assert not (t.root / 'controller_refresh.json').exists()


def test_all_handles_checked_before_first_stop(live_tree, monkeypatch):
    t = live_tree
    def unavailable(pid):
        if pid == 102:
            raise PermissionError('second case is inaccessible')
        return pid
    monkeypatch.setattr(run_refresh, '_pidfd_open', unavailable)
    with pytest.raises(PermissionError, match='inaccessible'):
        run_refresh.refresh(t.root, apply=True)
    assert not t.signals and not t.launches


def test_failed_old_preparation_accepts_fixed_code_without_manual_cleanup(live_tree):
    t = live_tree
    before = scientific_files(t.root)
    archive = t.root.parent / 'failed_refresh'; archive.mkdir()
    atomic_json(t.root / 'controller_refresh.json', {
        'phase': 'preparing', 'code': str(t.root / '.salvage_code.bad'),
        'previous_code': str(t.root / '.salvage_code.old'), 'archive': str(archive),
    })
    result = run_refresh.refresh(t.root, apply=True)
    assert result['phase'] == 'launched' and len(t.launches) == 2
    assert result['superseded_preparation_codes'] == [str(t.root / '.salvage_code.bad')]
    assert scientific_files(t.root) == before


@pytest.mark.parametrize('phase', ['stopping', 'stopped', 'prepared'])
def test_code_replacement_rejects_an_already_committed_handoff(live_tree, phase):
    t = live_tree
    archive = t.root.parent / 'failed_refresh'; archive.mkdir()
    atomic_json(t.root / 'controller_refresh.json', {
        'phase': phase, 'code': str(t.root / '.salvage_code.bad'),
        'previous_code': str(t.root / '.salvage_code.old'), 'archive': str(archive),
    })
    with pytest.raises(RuntimeError, match='already changed'):
        run_refresh.refresh(t.root, apply=True)
    assert not t.signals and not t.launches


@pytest.mark.parametrize('host,arch', [('darwin', 'arm64'), ('linux', 'unknown')])
def test_syscall_fallback_refuses_other_platforms_before_loading_libc(monkeypatch, host, arch):
    monkeypatch.setattr(run_refresh.sys, 'platform', host)
    monkeypatch.setattr(run_refresh.platform, 'machine', lambda: arch)
    def forbidden(*args, **kwargs):
        pytest.fail('must not issue a Linux syscall on an unrecognized ABI')
    monkeypatch.setattr(run_refresh.ctypes, 'CDLL', forbidden)
    with pytest.raises(RuntimeError, match='compatibility requires'):
        run_refresh._pidfd_syscall(434, ctypes.c_int(1), ctypes.c_uint(0))


@pytest.mark.skipif(sys.platform != 'linux', reason='requires the real Linux pidfd ABI')
def test_real_linux_process_handles_without_python_wrappers(monkeypatch):
    """Linux CI exercises the native ABI; no controller or scientific job is used."""
    monkeypatch.delattr(os, 'pidfd_open', raising=False)
    monkeypatch.delattr(signal, 'pidfd_send_signal', raising=False)
    run_refresh._probe_process_handles()
    child = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(30)'])
    fd = None
    try:
        fd = run_refresh._pidfd_open(child.pid)
        run_refresh._pidfd_send_signal(fd, signal.SIGSTOP)
        _, status = os.waitpid(child.pid, os.WUNTRACED)
        assert os.WIFSTOPPED(status)
        run_refresh._pidfd_send_signal(fd, signal.SIGCONT)
        run_refresh._pidfd_send_signal(fd, signal.SIGKILL)
        assert child.wait(timeout=5) == -signal.SIGKILL
        with pytest.raises(ProcessLookupError):
            run_refresh._pidfd_send_signal(fd, 0)
    finally:
        if child.poll() is None:
            child.kill(); child.wait(timeout=5)
        if fd is not None:
            os.close(fd)
