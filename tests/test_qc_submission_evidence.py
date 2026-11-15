"""Real shared admission with ordinary QC receipts and a simulated scheduler."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

from kinbot import run_policy, slurm_retry
from kinbot.anl import dispatch
from kinbot.execution import atomic_json
from tests.anl_fixture import dispatch_spec, molpro_task


@pytest.fixture
def pool(tmp_path, monkeypatch):
    run_policy.update(tmp_path, max_jobs=20, fair_submission=False)
    jobs, calls = {}, []
    def scheduler(command, **kwargs):
        if command[0] == 'squeue':
            return SimpleNamespace(returncode=0, stdout='\n'.join(
                f'{job}|1|{directory}' for job, directory in jobs.items()), stderr='')
        if command[:3] == ['scontrol', 'show', 'job']:
            return SimpleNamespace(returncode=1, stdout='',
                                   stderr='slurm_load_jobs error: Invalid job id specified')
        assert command[0] == 'sbatch'
        assert len(jobs) < 20
        job = str(1000+len(calls))
        calls.append((command, kwargs)); jobs[job] = str(kwargs['cwd'])
        return SimpleNamespace(returncode=0, stdout=job+'\n', stderr='')
    monkeypatch.setattr(run_policy.subprocess, 'run', scheduler)
    return SimpleNamespace(root=tmp_path, jobs=jobs, calls=calls, scheduler=scheduler)


def ordinary(pool, job, identity='hir/minimum_hir_0_00', *, code=0, journal=True,
             retry_script=False):
    directory = pool.root/('worker_'+job); directory.mkdir()
    script = directory/('.kinbot_timeout_retries/token/retry.sbatch' if retry_script
                        else identity+'.sbatch')
    script.parent.mkdir(parents=True, exist_ok=True)
    script.write_text(slurm_retry.wrap_script('#!/bin/bash\npython worker.py\n',
        identity=identity, qc_job=identity, exit_marker=identity+'.exitcode'))
    command = ['sbatch', str(script) if retry_script else str(script.relative_to(directory))]
    claim = dict(directory=str(directory), command=command, job_id=job,
                 status='submitted', created_at='2026-10-29T14:19:05+00:00')
    atomic_json(pool.root/'.kinbot_policy_submissions'/(job+'.json'), claim)
    marker=directory/(identity+'.exitcode'); marker.parent.mkdir(exist_ok=True)
    marker.write_text(job+' '+str(code)+'\n')
    token=hashlib.sha256(identity.encode()).hexdigest()[:20]
    folder=directory/'.kinbot_slurm_io_retries'/token/job
    if journal:
        folder.mkdir(parents=True); (folder/'lock').touch()
        atomic_json(folder/'execution.json', dict(schema=1, identity=identity,
            slurm_job_id=job, slurm_restart_count=0, status='failed' if code else 'complete',
            attempts=[dict(returncode=code, failed_result=bool(code))]))
    return claim, marker, folder


@pytest.mark.parametrize('journal', [True, False])
@pytest.mark.parametrize('code', [0, 1])
@pytest.mark.parametrize('retry_script', [True, False])
def test_ordinary_attempt_releases_capacity_without_acceptance_or_retry(pool, journal, code, retry_script):
    claim, marker, folder=ordinary(pool, '123', journal=journal, code=code, retry_script=retry_script)
    atomic_json(pool.root/'.kinbot_policy_admissions.json', {'claim':claim})
    protected={p:p.read_bytes() for p in Path(claim['directory']).rglob('*') if p.is_file()}
    assert run_policy.reconcile_allocations(pool.root)==0
    assert json.loads((pool.root/'.kinbot_policy_admissions.json').read_text())=={}
    proof=json.loads((pool.root/'.kinbot_policy_reconciled/claim.json').read_text())
    assert proof['claim']==claim
    assert proof['execution_status']==('failed' if code else 'complete')
    assert not pool.calls
    assert all(p.read_bytes()==content for p,content in protected.items())
    assert run_policy.reconcile_allocations(pool.root)==0
    assert len(list((pool.root/'.kinbot_policy_reconciled').iterdir()))==1


@pytest.mark.parametrize('fault', ['live','unknown_id','receipt_mismatch','bare_marker',
    'wrong_marker','missing_marker','running_journal','bad_journal','requeue_running',
    'wrong_identity','missing_script','scheduler_down','empty_scheduler_reply'])
def test_uncertain_or_live_attempt_keeps_reservation(pool, monkeypatch, fault):
    claim, marker, folder=ordinary(pool, '123', journal=fault in
        ('running_journal','bad_journal','requeue_running'))
    if fault=='live':pool.jobs['123']=claim['directory']
    elif fault=='unknown_id':claim.pop('job_id')
    elif fault=='receipt_mismatch':claim['created_at']='other'
    elif fault=='bare_marker':marker.write_text('0\n')
    elif fault=='wrong_marker':marker.write_text('122 0\n')
    elif fault=='missing_marker':marker.unlink()
    elif fault in ('running_journal','bad_journal','requeue_running'):
        record=json.loads((folder/'execution.json').read_text())
        if fault=='bad_journal':(folder/'execution.json').write_text('{')
        else:
            record['status']='running'
            if fault=='requeue_running':
                folder=folder/'restart-1';folder.mkdir();record['slurm_restart_count']=1
            atomic_json(folder/'execution.json',record)
    elif fault in ('wrong_identity','missing_script'):
        script=Path(claim['directory'])/claim['command'][-1]
        if fault=='missing_script':script.unlink()
        else:script.write_text(script.read_text().replace('--qc-job hir/', '--qc-job wrong/'))
    elif fault in ('scheduler_down','empty_scheduler_reply'):
        def reply(command,**kwargs):
            if command[0]=='scontrol':
                return SimpleNamespace(returncode=1 if fault=='scheduler_down' else 0,
                    stdout='',stderr='Unable to contact controller' if fault=='scheduler_down' else '')
            return pool.scheduler(command,**kwargs)
        monkeypatch.setattr(subprocess,'run',reply)
    atomic_json(pool.root/'.kinbot_policy_admissions.json',{'claim':claim})
    assert run_policy.reconcile_allocations(pool.root)==1
    claims=json.loads((pool.root/'.kinbot_policy_admissions.json').read_text())
    assert claims==({} if fault=='live' else {'claim':claim})
    assert not pool.calls
    assert not (pool.root/'.kinbot_policy_reconciled').exists()


def test_old_marker_overwritten_by_new_attempt_uses_allocation_journal(pool):
    claim,marker,_=ordinary(pool,'123')
    marker.write_text('999 0\n')
    pool.jobs['999']=claim['directory']
    atomic_json(pool.root/'.kinbot_policy_admissions.json', {'claim':claim})
    assert run_policy.reconcile_allocations(pool.root)==1
    assert marker.read_text()=='999 0\n'
    assert pool.jobs=={'999':claim['directory']}


def test_full_dispatch_fills_freed_capacity_preserves_live_work_and_repeats(pool, monkeypatch):
    """The reported 8 live + 12 old claims must admit exactly 12 ready tasks."""
    claims={}
    for n in range(12):
        claims[str(n)]=ordinary(pool,str(100+n),code=n%2)[0]
    atomic_json(pool.root/'.kinbot_policy_admissions.json',claims)
    pool.jobs.update({str(n):str(pool.root/'live') for n in range(8)})
    protected={p:p.read_bytes() for p in pool.root.glob('worker_*/*/*') if p.is_file()}
    spec=dispatch_spec(); task=molpro_task('first','rhf\n',geometry_from='initial')
    spec['limits']['max_nodes']=20
    spec['tasks']=[dict(deepcopy(task),id='task_'+str(n)) for n in range(15)]
    graph=dispatch.prepare(spec,pool.root/'ready')
    # Existing allocations stay active; dispatch must neither reset nor relaunch.
    monkeypatch.setattr(dispatch,'_job_active',lambda job:job in pool.jobs)
    state=dispatch.advance(graph,submit=True)
    assert len(pool.calls)==12 and len(pool.jobs)==20
    submitted={n:e['job_id'] for n,e in state['tasks'].items() if e['status']=='submitted'}
    assert len(submitted)==12
    assert len(list((pool.root/'.kinbot_policy_reconciled').glob('*.json')))==12
    state=dispatch.advance(graph,submit=True)
    assert len(pool.calls)==12
    assert all(state['tasks'][n]['job_id']==job for n,job in submitted.items())
    assert all(p.read_bytes()==content for p,content in protected.items())


def test_native_exit_marker_before_queue_exit_cannot_release_twice(pool):
    claim,_,_=ordinary(pool,'123')
    pool.jobs['123']=claim['directory']
    atomic_json(pool.root/'.kinbot_policy_admissions.json',{'claim':claim})
    assert run_policy.reconcile_allocations(pool.root)==1
    # The live job is counted from Slurm after its provisional claim is retired.
    assert run_policy.reconcile_allocations(pool.root)==1
    assert not pool.calls


def test_completion_receipt_api_is_unchanged(pool):
    claim,_,_=ordinary(pool,'123',code=1)
    assert slurm_retry.completion_receipt('hir/minimum_hir_0_00','123',claim['directory'])==(1,True)
