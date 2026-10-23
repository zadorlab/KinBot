"""Controller handoff preserves owned native work and rejects unsafe pauses."""
from copy import deepcopy
import ctypes
import errno
import fcntl
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

PIN_CODE = run_salvage._pin_code


@pytest.fixture
def sampling_tree(live_tree, monkeypatch):
    from kinbot import rotdpy_campaign, run_policy
    t = live_tree
    directory = t.root/'ch3s'/'sampling'
    directory.mkdir()
    inp = directory/'channel.py'
    inp.write_text('validation_metadata = {}\n')
    identity = dict(host=socket.gethostname(), pid=700, start='7000', boot='test')
    monkeypatch.setattr('kinbot.pes_phase.process_identity', lambda pid: identity if pid == 700 else None)
    atomic_json(rotdpy_campaign.receipt(inp), dict(input=str(inp), owner_identity=identity,
                status='launched', **rotdpy_campaign.contract(inp)))
    t.processes.extend([
        dict(pid=700, ppid=101, starttime='7000', cwd=str(directory), state='S',
             command=f'/python -u -m kinbot.rotdpy run {inp}'),
        dict(pid=701, ppid=700, starttime='7010', cwd=str(directory), state='S',
             command='/python -m kinbot.rotdpy_policy channel.py'),
        dict(pid=702, ppid=701, starttime='7020', cwd=str(directory), state='S',
             command='/python -m kinbot.run_policy submit job.sh')])
    jobs = directory/'kb_channel/Surface_0/jobs'; jobs.mkdir(parents=True)
    (jobs/'job.sh').write_text('#!/bin/sh\n')
    atomic_json(t.root/run_policy.NAME, {})
    atomic_json(t.root/'.kinbot_policy_submissions/800.json', dict(job_id='800',
                directory=str(jobs), command=['sbatch', 'job.sh'], status='submitted'))
    t.jobs.append(dict(job_id='800', directory=str(jobs), state='RUNNING'))
    t.campaign = inp
    return t


def test_refresh_preserves_campaign_tree_and_allocations(sampling_tree):
    t = sampling_tree
    before = t.campaign.with_suffix('.campaign.json').read_bytes()
    result = run_refresh.refresh(t.root, apply=True)
    assert result['phase'] == 'launched'
    assert {p['pid'] for p in t.processes} == {700, 701, 702}
    assert not any(pid in {700, 701, 702} for pid, sig in t.signals)
    assert t.campaign.with_suffix('.campaign.json').read_bytes() == before
    assert len(t.launches) == 2
    assert any(j['job_id'] == '800' for j in t.jobs)


def test_refresh_imports_terminal_search_outcomes_before_replacement(live_tree):
    from kinbot import search_state
    t = live_tree
    worker = t.root/'ch3s/123'
    inp = worker/'123.json'
    atomic_json(inp, {'smiles': 'CC'})
    monitor = worker/'kinbot_monitor.out'
    monitor.write_text('-999 -1 rejected\n-1 13 accepted product\n3 13 optimizing\n')
    before = (worker/'kinbot.db').read_bytes()
    result = run_refresh.refresh(t.root, apply=True)
    saved = json.loads((worker/search_state.NAME).read_text())
    assert saved['reactions']['rejected']['status'] == -999
    assert saved['reactions']['accepted']['status'] == -1
    assert saved['reactions']['optimizing']['status'] == 3
    archive = Path(result['archive'])/'ch3s/123'
    assert (archive/search_state.NAME).read_bytes() == (worker/search_state.NAME).read_bytes()
    assert (archive/'kinbot_monitor.out').read_bytes() == monitor.read_bytes()
    assert (worker/'kinbot.db').read_bytes() == before
    assert not any(j.get('state') == 'CANCELLED' for j in t.jobs)
    assert len(t.launches) == 2


@pytest.mark.parametrize('source', ['snapshot', 'archived_log'])
def test_refresh_restores_hash_proven_initial_input_and_keeps_jobs(live_tree, source):
    from kinbot import search_state
    t = live_tree
    worker = t.root/'ch3s/123'; inp = worker/'123.json'
    original = dict(smiles='', structure=['H', 0., 0., 0.])
    atomic_json(inp, original)
    (worker/'kinbot_monitor.out').write_text('-999 -1 rejected\n')
    search_state.capture_legacy(inp)
    if source == 'archived_log':
        record = json.loads((worker/search_state.NAME).read_text())
        record.pop('input_snapshot')
        atomic_json(worker/search_state.NAME, record)
        old_archive = t.root.with_name(t.root.name+'.controllers.previous')
        log = old_archive/'ch3s/123/kinbot.log'
        log.parent.mkdir(parents=True)
        log.write_text('INFO: Input parameters\n\tstructure: '+repr(original['structure'])+
                       '\nINFO: Starting KinBot\n')
        atomic_json(t.root/'controller_refresh.json', dict(
            phase='launched', code=str(t.root/'.salvage_code.old'), archive=str(old_archive)))
    checkpoint = (worker/search_state.NAME).read_bytes()
    atomic_json(inp, dict(original, structure=['H', 1., 2., 3.]))
    overwritten = inp.read_bytes()
    chemistry = (worker/'kinbot.db').read_bytes()
    result = run_refresh.refresh(t.root, apply=True)
    assert json.loads(inp.read_text()) == original
    assert (Path(result['archive'])/'ch3s/123/123.json').read_bytes() == overwritten
    assert (worker/search_state.NAME).read_bytes() == checkpoint
    assert (worker/'kinbot.db').read_bytes() == chemistry
    assert all(j['state'] == 'RUNNING' for j in t.jobs)
    assert len(t.launches) == 2


def test_unproven_search_input_change_fails_preflight_before_signalling(live_tree, monkeypatch):
    from kinbot import search_state
    t = live_tree
    worker = t.root/'ch3s/123'; inp = worker/'123.json'
    atomic_json(inp, {'structure': ['H', 0., 0., 0.], 'barrier_threshold': 100.})
    (worker/'kinbot_monitor.out').write_text('-999 -1 rejected\n')
    search_state.capture_legacy(inp)
    before = (worker/search_state.NAME).read_bytes()
    atomic_json(inp, {'structure': ['H', 1., 2., 3.], 'barrier_threshold': 50.})
    monkeypatch.setattr(run_refresh, '_preflight', run_refresh._search_input_plans)
    with pytest.raises(RuntimeError, match='inputs/theory changed'):
        run_refresh.refresh(t.root, apply=True)
    assert not t.signals and not t.launches
    assert (worker/search_state.NAME).read_bytes() == before


@pytest.mark.parametrize('phase_bytes', [b'', b'{"phase":', b'{"phase":"exploring"}'])
def test_refresh_archives_old_phase_before_relaunch_after_owner_loss(live_tree, phase_bytes):
    from kinbot.input_refresh import owner_directory
    t = live_tree
    t.processes.clear()
    before = {p: p.read_bytes() for g in t.graphs for p in (g/'state.json', g/'workflow.json')}
    for case in ('ch3s', 'ch3ss'):
        (t.root/case/'pes_exploration.json').write_bytes(phase_bytes)
        (t.root/case/'chemids').write_text('123\n')
        (t.root/case/'pes_l2_selection.json').write_text('{"status":"selected"}')
    result = run_refresh.refresh(t.root, apply=True)
    archive = Path(result['archive'])
    for case in ('ch3s', 'ch3ss'):
        assert (archive/case/'pes_exploration.json').read_bytes() == phase_bytes
        assert not (t.root/case/'pes_exploration.json').exists()
        assert owner_directory(t.root/case/'input.json', directory=t.root/case) is None
        assert (t.root/case/'chemids').read_text() == '123\n'
        assert (t.root/case/'pes_l2_selection.json').read_text() == '{"status":"selected"}'
    assert all(p.read_bytes() == data for p, data in before.items())
    assert not t.signals
    assert all(j['state'] == 'RUNNING' for j in t.jobs)
    assert len(t.launches) == 2


@pytest.mark.parametrize('fault', ['owner', 'input', 'claim', 'orphan'])
def test_refresh_refuses_unproven_campaign_before_signalling(sampling_tree, fault):
    t = sampling_tree
    if fault == 'owner':
        p = t.campaign.with_suffix('.campaign.json')
        record = json.loads(p.read_text()); record['owner_identity']['start'] = 'wrong'
        atomic_json(p, record)
    elif fault == 'input':
        t.campaign.write_text('validation_metadata = {"changed": True}\n')
    elif fault == 'claim':
        (t.root/'.kinbot_policy_submissions/800.json').unlink()
    else:
        t.processes[:] = [p for p in t.processes if p['pid'] != 700]
    with pytest.raises(RuntimeError):
        run_refresh.refresh(t.root, apply=True)
    assert not t.signals and not t.launches


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


@pytest.fixture
def fresh_tree(live_tree):
    t = live_tree
    pin = t.root / '.fresh_code.original'
    pin.mkdir()
    launches = {}
    for case in ('ch3s', 'ch3ss'):
        directory = t.root / case
        receipt = json.loads((directory/'launch.json').read_text())
        receipt.update(code=str(pin), input=str(directory/'input.json'))
        atomic_json(directory/'launch.json', receipt)
        launches[case] = receipt
    atomic_json(t.root/'fresh_run.json', dict(schema=1, phase='launched', code=str(pin),
                                            launches=launches))
    (t.root/'salvage.json').unlink()
    # One controller tree has exited; the other still owns expensive jobs.
    t.processes[:] = [p for p in t.processes if '/ch3ss' not in p['cwd']]
    return t


def test_fresh_refresh_preserves_qc_and_original_launch_and_restarts_once(fresh_tree):
    t = fresh_tree
    before = scientific_files(t.root)
    original_jobs = deepcopy(t.jobs)
    original_receipt = (t.root/'fresh_run.json').read_bytes()
    prepared = (t.root/'prepared.json').read_bytes()
    assert run_refresh.refresh(t.root)['phase'] == 'reviewed'
    assert not (t.root/'salvage.json').exists()
    assert not t.signals and not t.launches
    result = run_refresh.refresh(t.root, apply=True)
    assert result['phase'] == 'launched'
    assert t.jobs == original_jobs and scientific_files(t.root) == before
    assert (t.root/'fresh_run.json').read_bytes() == original_receipt
    assert (t.root/'prepared.json').read_bytes() == prepared
    assert {pid for pid, sig in t.signals if sig == signal.SIGKILL} == {100, 101}
    assert len(t.launches) == 2
    for _, kwargs in t.launches:
        for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
            assert kwargs['env'][key] == '1'
    assert run_refresh.refresh(t.root, apply=True) == result
    assert len(t.launches) == 2


@pytest.mark.parametrize('problem', ['incomplete', 'receipt', 'code', 'input'])
def test_fresh_refresh_refuses_inconsistent_launch_before_signalling(fresh_tree, problem):
    t = fresh_tree
    path = t.root/'fresh_run.json'
    m = json.loads(path.read_text())
    if problem == 'incomplete':
        m['phase'] = 'cleaned'
    elif problem == 'receipt':
        m['launches']['ch3s']['pid'] += 1
    elif problem == 'code':
        m['code'] = str(t.root.parent/'.fresh_code.other')
    else:
        (t.root/'ch3s/input.json').write_text('{"changed":true}')
    atomic_json(path, m)
    with pytest.raises((RuntimeError, ValueError)):
        run_refresh.refresh(t.root, apply=True)
    assert not t.signals and not t.launches
    assert not (t.root/'salvage.json').exists()


def test_fresh_refresh_recovers_after_launch_receipts_archived(fresh_tree, monkeypatch):
    t = fresh_tree
    original = run_refresh._archive_controllers
    def interrupted(root, archive):
        original(root, archive)
        raise RuntimeError('interrupted after archiving controller receipts')
    monkeypatch.setattr(run_refresh, '_archive_controllers', interrupted)
    before = scientific_files(t.root)
    with pytest.raises(RuntimeError, match='interrupted after'):
        run_refresh.refresh(t.root, apply=True)
    assert not (t.root/'ch3s/launch.json').exists()
    monkeypatch.setattr(run_refresh, '_archive_controllers', original)
    assert run_refresh.refresh(t.root, apply=True)['phase'] == 'launched'
    assert scientific_files(t.root) == before
    assert len(t.launches) == 2


def test_fresh_refresh_pins_current_code_without_rewriting_original(fresh_tree, monkeypatch):
    t = fresh_tree
    monkeypatch.setattr(run_salvage, '_pin_code', PIN_CODE)
    original = (t.root/'fresh_run.json').read_bytes()
    review = run_refresh.refresh(t.root)
    pin = Path(review['code'])
    assert pin.name.startswith('.salvage_code.')
    manifest = json.loads((pin/'manifest.json').read_text())
    assert manifest['kinbot/frequencies.py'] == file_sha256(Path(run_salvage.__file__).with_name('frequencies.py'))
    result = run_refresh.refresh(t.root, apply=True)
    assert result['code'] == str(pin)
    assert (t.root/'fresh_run.json').read_bytes() == original
    assert run_refresh.refresh(t.root, apply=True) == result
    assert len(t.launches) == 2


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


def test_refresh_preserves_profiled_exploration_allocations_and_routes(live_tree):
    t = live_tree
    worker = t.root/'ch3ss/123'
    atomic_json(worker/'123.json', {})
    routes = {}
    for job, phase, name in [('901', 'RUNNING', '123_irc'), ('902', 'PENDING', '123_conf_0001')]:
        routes[name] = dict(level='l1', job_id=job, generation_tagged=True)
        (worker/(name+'.py')).write_text('# expensive worker\n')
        (worker/(name+'.sbatch')).write_text('#!/bin/bash\n')
        t.jobs.append(dict(job_id=job, directory=str(worker), state=phase))
    path = worker/'.kinbot_theory_jobs.json'
    atomic_json(path, dict(schema=1, jobs=routes))
    expected = path.read_bytes()
    old_jobs = deepcopy(t.jobs)
    result = run_refresh.refresh(t.root, apply=True)
    assert result['phase'] == 'launched'
    assert t.jobs == old_jobs
    assert path.read_bytes() == expected
    assert {j['job_id'] for j in result['jobs']} >= {'901', '902'}


def test_refresh_rejects_wrong_directory_in_ordinary_route(live_tree):
    t = live_tree
    worker = t.root/'ch3ss/123'
    atomic_json(worker/'123.json', {})
    atomic_json(worker/'.kinbot_theory_jobs.json', dict(schema=1, jobs={
        '123_irc': dict(level='l1', job_id='901', generation_tagged=True)}))
    t.jobs.append(dict(job_id='901', directory=str(t.root/'ch3s/123'), state='RUNNING'))
    with pytest.raises(RuntimeError, match='ordinary QC'):
        run_refresh.refresh(t.root, apply=True)
    assert not t.signals and not t.launches


def test_network_revisit_handoff_preserves_live_and_pending_graphs(live_tree):
    from kinbot.network_revisit import settings
    t=live_tree
    for case in ('ch3s','ch3ss'):
        worker=t.root/case/'123'
        atomic_json(worker/'123.json',dict(vdW_detection=.5))
        (worker/'kinbot_monitor.out').write_text('-1 0 accepted\n-999 0 retry_this\n')
    t.jobs[-1]['state']='PENDING'
    before=scientific_files(t.root)
    jobs=deepcopy(t.jobs)
    options=settings(0.,True,True)
    preview=run_refresh.refresh(t.root,network_options=options)
    assert preview['phase']=='reviewed' and not t.signals
    result=run_refresh.refresh(t.root,apply=True,network_options=options)
    assert result['phase']=='launched' and t.jobs==jobs
    assert scientific_files(t.root)==before
    for case in ('ch3s','ch3ss'):
        assert json.loads((t.root/case/'input.json').read_text())['vdW_detection']==0
        assert json.loads((t.root/case/'123/.kinbot_revisit.json').read_text())['rejected']==['retry_this']
    assert run_refresh.refresh(t.root,apply=True,network_options=options)==result
    assert len(t.launches)==2


def test_deduplicate_refresh_cancels_only_verified_duplicate_and_retains_chain(live_tree, monkeypatch):
    from tests.test_vrc_graph import graph, spec
    t = live_tree
    scans = t.root/'ch3s/123/vrctst/molpro_scan'
    original = graph(scans/'asymptote_old', spec(), job='1000')
    duplicate = graph(scans/'asymptote_new', spec(walltime='7-00:00:00'), 'submitted', '1001')
    t.jobs.append(dict(job_id='1001', directory=str(duplicate/'tasks/point'), state='RUNNING'))
    preserved = {p: file_sha256(p) for g in (*t.graphs, original) for p in g.rglob('*') if p.is_file()}
    cancelled = []
    terminal = {}
    def record(ident):
        job = next((j for j in t.jobs if j['job_id'] == ident), terminal.get(ident))
        return {} if job is None else dict(JobId=ident, WorkDir=job['directory'], JobState=job['state'])
    def cancel(command):
        assert command == ['scancel', '1001']
        job = next(j for j in t.jobs if j['job_id'] == command[1])
        t.jobs.remove(job)
        terminal['1001'] = dict(job, state='CANCELLED')
        cancelled.append(command[1])
        return ''
    monkeypatch.setattr('kinbot.execution.slurm_job_record', record)
    monkeypatch.setattr(run_refresh, '_pending_request', cancel)
    review = run_refresh.refresh(t.root, deduplicate_vrc=True)
    assert review['phase'] == 'reviewed' and not t.signals and not cancelled
    result = run_refresh.refresh(t.root, apply=True, deduplicate_vrc=True)
    assert result['phase'] == 'launched' and cancelled == ['1001']
    assert len(t.jobs) == 4 and len(t.launches) == 2
    assert preserved == {p: file_sha256(p) for p in preserved}
    assert not duplicate.exists() and original.exists()
    archived = Path(result['archive'])/'vrc_duplicates'/duplicate.relative_to(t.root)
    assert (archived/'workflow.json').is_file()
    assert run_refresh.refresh(t.root, apply=True, deduplicate_vrc=True) == result
    assert cancelled == ['1001'] and len(t.launches) == 2


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


def test_refresh_waits_for_real_worker_lock_release(live_tree, monkeypatch, capsys):
    t = live_tree
    path = t.root/'ch3s/123/.kinbot_run.json'
    before = (path.stat().st_ino, path.read_bytes(), deepcopy(t.jobs))
    with path.open('r+') as holder:
        fcntl.flock(holder, fcntl.LOCK_EX | fcntl.LOCK_NB)
        waited = []
        def release(seconds):
            waited.append(seconds)
            fcntl.flock(holder, fcntl.LOCK_UN)
        monkeypatch.setattr(run_refresh.time, 'sleep', release)
        assert run_refresh.refresh(t.root, apply=True)['phase'] == 'launched'
    assert waited and len(t.launches) == 2
    assert 'Waiting for worker lock release: ch3s/123/.kinbot_run.json' in capsys.readouterr().out
    assert before == (path.stat().st_ino, path.read_bytes(), t.jobs)


@pytest.fixture
def vrc_migration_tree(live_tree, monkeypatch):
    from kinbot import vrc_restart
    import shutil
    t = live_tree
    monkeypatch.setattr(vrc_restart, 'validate_input', lambda *a, **kw: None)
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


@pytest.mark.parametrize('unsafe_change', [None, 'input', 'cancellation', 'archive'])
def test_stopped_lock_failure_can_resume_with_fixed_code_only_before_edits(
        vrc_migration_tree, monkeypatch, unsafe_change):
    from kinbot import vrc_restart
    t = vrc_migration_tree
    options = dict(vrc_updates=vrc_restart.settings('avdz-f12', 'asymptotic_scan', 'L2',
                    ['avdz-f12', 'avtz-f12'], 3.7), cancel_superseded_vrc_scans=True)
    worker_locks = run_refresh._worker_locks
    monkeypatch.setattr(run_refresh, '_worker_locks', lambda root: worker_locks(root, timeout=0))
    before = scientific_files(t.root)
    with (t.root/'ch3s/123/.kinbot_run.json').open('r+') as holder:
        fcntl.flock(holder, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(RuntimeError, match='Worker lock still held:'):
            run_refresh.refresh(t.root, apply=True, **options)
    assert not t.processes and not t.launches and not t.commands
    marker = t.root/'controller_refresh.json'
    journal = json.loads(marker.read_text())
    assert journal['phase'] == 'stopped' and not journal['vrc_migration'].get('cancellations')
    assert scientific_files(t.root) == before
    pin = t.root/'.salvage_code.fixed'; pin.mkdir()
    monkeypatch.setattr(run_salvage, '_pin_code', lambda *a, **kw: pin)
    monkeypatch.setattr(run_refresh, '_worker_locks', worker_locks)
    if unsafe_change == 'input':
        (t.root/'ch3s/input.json').write_text('{"changed":true}')
    elif unsafe_change == 'cancellation':
        journal['vrc_migration']['cancellations'] = {'501': {'phase': 'requested'}}
        atomic_json(marker, journal)
    elif unsafe_change == 'archive':
        (Path(journal['archive'])/'partial_archive').write_text('changed')
    if unsafe_change:
        with pytest.raises(RuntimeError, match='already changed'):
            run_refresh.refresh(t.root, apply=True, **options)
        assert not t.launches and not t.commands
        return
    result = run_refresh.refresh(t.root, apply=True, **options)
    assert result['phase'] == 'launched' and result['code'] == str(pin)
    assert len(t.launches) == 2
    assert t.commands == [['scancel', '501'], ['scancel', '503']]
    assert scientific_files(t.root) == before
    assert run_refresh.refresh(t.root, apply=True, **options) == result
    assert len(t.launches) == 2


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
    archive = t.root.parent / (t.root.name+'.controllers.failed'); archive.mkdir()
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


def test_refresh_adopts_nested_reference_allocations_and_replaces_only_its_controller(live_tree):
    t=live_tree;worker=t.root/'ch3s/123'
    atomic_json(worker/'123.json',{})
    ref=worker/'anl_automatic/reference_l2/state_ref';ref.mkdir(parents=True)
    atomic_json(ref/'input.json',{})
    atomic_json(ref/'.kinbot_theory_jobs.json',dict(schema=1,jobs={'reference_high':dict(level='l2',job_id='908',generation_tagged=True)}))
    (ref/'reference_high.py').write_text('worker')
    (ref/'reference_high.sbatch').write_text('#!/bin/bash\n')
    (ref/'.kinbot_run.json').write_text('{}')
    (ref/'driver.pid').write_text('108')
    (ref/'driver.stderr').write_text('old reference controller log')
    atomic_json(ref/'kinbot_status.json',dict(pid=108,host=socket.gethostname(),status='running'))
    t.processes.append(dict(pid=108,ppid=101,starttime='1080',cwd=str(ref),state='S',
                           command='/python -m kinbot.kb '+str(ref/'input.json')))
    t.jobs.append(dict(job_id='908',directory=str(ref),state='RUNNING'))
    result=run_refresh.refresh(t.root,apply=True)
    assert result['phase']=='launched' and t.jobs[-1]['job_id']=='908'
    assert (108,signal.SIGKILL) in t.signals
    assert (ref/'reference_high.sbatch').read_text()=='#!/bin/bash\n'
    assert (Path(result['archive'])/ref.relative_to(t.root)/'driver.stderr').exists()


def test_direct_cbs_and_reviewed_channel_share_one_handoff(vrc_migration_tree, tmp_path):
    from kinbot import vrc_restart
    from tests.test_input_refresh import surface
    t = vrc_migration_tree
    first = surface()
    first['channels'][0]['bottlenecks'][0]['parameters'].update(
        vrc_tst_scan_method='caspt2-f12(5,4)', vrc_tst_correction_model='spin_splitting')
    existing = surface()
    existing['channels'][0]['bottlenecks'][0]['parameters'].update(
        vrc_tst_scan_method='caspt2-f12(3,4)', vrc_tst_correction_model='spin_splitting',
        vrc_tst_active_space=dict(mode='selected', active_electrons=3, active_orbitals=4))
    payload = dict(vrc_tst_scan_method='caspt2-f12(5,4)', vrc_tst_scan_basis='avdz-f12',
        vrc_tst_fragment_geometry='asymptotic_scan', vrc_tst_correction_model='spin_splitting',
        vrc_tst_spin_splitting={'reference_multiplicity':4},
        vrc_tst_surface_model=first, vrc_tst_surface_models={'456':existing})
    _, paths = run_salvage._cases(t.root)
    for path in paths:
        atomic_json(path, payload)
        atomic_json(path.parent/'resolved_parameters.json', payload)
        atomic_json(path.parent/'123/123.json', payload)
        atomic_json(path.parent/'123/.kinbot_search_options.json',dict(
            vdW_detection=.5,vrc_tst_surface_models=payload['vrc_tst_surface_models']))
    prepared=json.loads((t.root/'prepared.json').read_text())
    prepared['sha256']={str(p.relative_to(t.root)):file_sha256(p) for p in paths}
    atomic_json(t.root/'prepared.json',prepared)
    new = surface()
    new['channels'][0]['bottlenecks'][0]['parameters'].update(
        vrc_tst_scan_method='caspt2-f12(5,4)', vrc_tst_correction_model='direct')
    evidence = tmp_path/'native.out';evidence.write_text('retained native screen')
    proof = tmp_path/'proposal.json'
    atomic_json(proof, dict(schema=1,kind='orbital_model_proposal',case_input=str(paths[0]),
        input_sha256=file_sha256(paths[0]), models={**payload['vrc_tst_surface_models'],'789':new},
        evidence={str(evidence):file_sha256(evidence)}))
    updates = vrc_restart.settings(direct_caspt2_cbs=True, reviewed_models=proof)
    before = scientific_files(t.root)
    reviewed = run_refresh.refresh(t.root, vrc_updates=updates)
    assert not reviewed['vrc_migration']['replace_geometry']
    assert reviewed['vrc_migration']['selected_jobs'] == []
    assert not t.signals
    result = run_refresh.refresh(t.root, apply=True, vrc_updates=updates)
    assert result['phase'] == 'launched' and len(t.jobs)==4 and not t.commands
    assert before == scientific_files(t.root)
    for path in paths:
        for p in (path,path.parent/'resolved_parameters.json',path.parent/'123/123.json'):
            data=json.loads(p.read_text())
            assert data['vrc_tst_correction_model']=='direct'
            assert 'vrc_tst_spin_splitting' not in data
            assert data['vrc_tst_high_basis']==['avdz-f12','avtz-f12']
            assert data['vrc_tst_high_extrapolation_power']==3.7
            par=data['vrc_tst_surface_models']['456']['channels'][0]['bottlenecks'][0]['parameters']
            assert par['vrc_tst_high_method']=='caspt2-f12(3,4)'
            assert par['vrc_tst_correction_model']=='direct'
            assert ('789' in data['vrc_tst_surface_models']) == (path==paths[0])
    for path in paths:
        overrides=json.loads((path.parent/'123/.kinbot_search_options.json').read_text())
        assert set(overrides)=={'vdW_detection','vrc_tst_surface_models'}
        assert overrides['vdW_detection']==.5
        par=overrides['vrc_tst_surface_models']['456']['channels'][0]['bottlenecks'][0]['parameters']
        assert par['vrc_tst_correction_model']=='direct'
        assert ('789' in overrides['vrc_tst_surface_models']) == (path==paths[0])
    assert vrc_restart.audit_direct_inputs(t.root)['input_files']==6
    # Exact journal replay does not reject its now-old proposal or re-edit files.
    assert run_refresh.refresh(t.root,apply=True,vrc_updates=updates)==result
    assert len(t.launches)==2


def test_direct_migration_refuses_changed_proposal_evidence_before_stop(vrc_migration_tree, tmp_path):
    from kinbot import vrc_restart
    t=vrc_migration_tree
    _,paths=run_salvage._cases(t.root)
    proof=tmp_path/'proposal.json';native=tmp_path/'native.out';native.write_text('old')
    atomic_json(proof,dict(schema=1,kind='orbital_model_proposal',case_input=str(paths[0]),
        input_sha256=file_sha256(paths[0]),models={},evidence={str(native):file_sha256(native)}))
    updates=vrc_restart.settings(direct_caspt2_cbs=True,reviewed_models=proof)
    native.write_text('changed')
    with pytest.raises(ValueError,match='evidence changed'):
        run_refresh.refresh(t.root,apply=True,vrc_updates=updates)
    assert not t.signals and len(t.jobs)==4
