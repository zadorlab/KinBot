"""Fresh owners, scheduler outages, requeues, and transaction interruption."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from ase.io import read
from ase.calculators.singlepoint import SinglePointCalculator
import numpy as np

from kinbot import run_policy, run_recover
from kinbot.anl import dispatch, allocation_recovery as recovery
from kinbot.anl.optimization_progress import Progress
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task
from tests.test_anl_dispatch import _mock_execution
from tests.test_optimization_progress import optimizer


@pytest.fixture
def site(tmp_path, monkeypatch):
    run_policy.update(tmp_path, max_jobs=3, allocation_retries=2,
                      allocation_retry_delay=10, allocation_retry_max_delay=30)
    now=[100.];live=set();records={}
    monkeypatch.setattr(recovery.time,'time',lambda:now[0])
    monkeypatch.setattr(recovery,'slurm_job_record',lambda job,**kw:records.get(job,{}))
    monkeypatch.setattr(dispatch,'_job_active',lambda job:job in live)
    return tmp_path,now,live,records


def graph_at(site, geometry=False, downstream=False):
    root,now,live,records=site
    spec=dispatch_spec()
    bad=deepcopy(spec['tasks'][0]) if geometry else molpro_task('bad','rhf\n',geometry_from='initial')
    bad['id']='bad';bad['failure_policy']='skip_component';bad['component_families']=['test']
    spec['tasks']=[bad,molpro_task('live','rhf\n',geometry_from='initial'),
                   molpro_task('done','rhf\n',geometry_from='initial')]
    if downstream:
        child=molpro_task('child','rhf\n',geometry_from='initial');child['depends_on']=['bad']
        spec['tasks'].append(child)
    graph=dispatch.prepare(spec,root/'graph')
    _mock_execution(graph,'done')
    state=dispatch._load(graph)[2]
    state['tasks']['bad'].update(status='submitted',job_id='42',missing_since=0)
    state['tasks']['live'].update(status='submitted',job_id='43');live.add('43')
    atomic_json(graph/'state.json',state)
    records['42']=dict(JobId='42',JobState='NODE_FAIL',WorkDir=str(graph/'tasks/bad'))
    return graph


def test_wait_persists_across_fresh_owner_outage_and_repeated_polls(site,monkeypatch):
    graph=graph_at(site);root,now,live,records=site
    original={p:file_sha256(p) for name in ('live','done') for p in (graph/'tasks'/name).iterdir() if p.is_file()}
    native=graph/'tasks/bad/partial.chk';native.write_bytes(b'old progress')
    first=dispatch.advance(graph)
    assert first['tasks']['bad']['status']=='retry_wait'
    assert first['tasks']['bad']['allocation_retry_wait']['retry_after']==110
    assert first['allocation_retries']['bad']==[]
    assert dispatch.advance(graph)['tasks']['bad']==first['tasks']['bad']
    def outage(job):raise dispatch.SchedulerUnavailableError('Slurm is down for hours')
    monkeypatch.setattr(dispatch,'_job_active',outage);now[0]=10000
    before=(graph/'state.json').read_bytes()
    with pytest.raises(dispatch.SchedulerUnavailableError):dispatch.advance(graph)
    assert (graph/'state.json').read_bytes()==before and native.read_bytes()==b'old progress'
    monkeypatch.setattr(dispatch,'_job_active',lambda job:job in live)
    records.clear()  # Original terminal scontrol record has now expired.
    restarted=dispatch.advance(graph);entry=restarted['tasks']['bad']
    assert entry['status']=='staged' and entry['attempt']==2
    assert len(restarted['allocation_retries']['bad'])==1
    assert (Path(entry['previous_attempt'])/'partial.chk').read_bytes()==b'old progress'
    assert original=={p:file_sha256(p) for p in original}
    assert dispatch.advance(graph)['allocation_retries']==restarted['allocation_retries']


@pytest.mark.parametrize('kind',['CANCELLED','FAILED','TIMEOUT','COMPLETED','unknown','wrong_directory','active'])
def test_no_retry_from_ambiguous_or_scientific_or_live_failure(site,kind):
    graph=graph_at(site);root,now,live,records=site
    if kind=='unknown':records.clear()
    elif kind=='wrong_directory':records['42']['WorkDir']=str(root/'different')
    elif kind=='active':live.add('42')
    else:records['42']['JobState']=kind
    state=dispatch.advance(graph)
    assert not state.get('allocation_retries')
    assert not (graph/'allocation_recovery.json').exists()


def test_scheduler_requeue_is_reattached_and_completed_work_is_adopted(site):
    graph=graph_at(site);root,now,live,records=site
    dispatch.advance(graph);live.add('42');now[0]=10000
    state=dispatch.advance(graph)
    assert state['tasks']['bad']['status']=='submitted'
    assert 'allocation_retry_wait' not in state['tasks']['bad']
    live.remove('42');_mock_execution(graph,'bad')
    state=dispatch.advance(graph)
    assert state['tasks']['bad']['status']=='complete'
    assert not state['allocation_retries']['bad']


def test_bounded_budget_and_policy_disable_never_create_optional_fallback(site):
    graph=graph_at(site);root,now,live,records=site
    for job,delay in [('42',10),('44',20)]:
        records[job]=dict(JobId=job,JobState='PREEMPTED',WorkDir=str(graph/'tasks/bad'))
        state=dispatch._load(graph)[2]
        state['tasks']['bad'].update(status='submitted',job_id=job,missing_since=0)
        atomic_json(graph/'state.json',state)
        wait=dispatch.advance(graph)
        assert wait['tasks']['bad']['allocation_retry_wait']['retry_after']==now[0]+delay
        now[0]+=delay
        assert dispatch.advance(graph)['tasks']['bad']['status']=='staged'
    state=dispatch._load(graph)[2]
    state['tasks']['bad'].update(status='submitted',job_id='45',missing_since=0)
    atomic_json(graph/'state.json',state)
    records['45']=dict(JobId='45',JobState='BOOT_FAIL',WorkDir=str(graph/'tasks/bad'))
    last=dispatch.advance(graph)
    assert last['tasks']['bad']['status']=='failed'
    assert last['tasks']['bad']['failure_phase']=='allocation'
    assert 'limit reached' in last['tasks']['bad']['allocation_retry_blocked']
    assert len(last['allocation_retries']['bad'])==2


def test_disable_pending_wait(site):
    graph=graph_at(site);root,now,live,records=site
    dispatch.advance(graph);run_policy.update(root,allocation_retries=0)
    assert dispatch.advance(graph)['tasks']['bad']['status']=='failed'


def test_live_descendant_prevents_reset(site):
    graph=graph_at(site,downstream=True);root,now,live,records=site
    state=dispatch._load(graph)[2]
    child=next(t for t in dispatch._load(graph)[1]['tasks'] if t['id']=='child')
    dispatch._stage_task(graph,dispatch._load(graph)[1],state,child)
    state['tasks']['child'].update(status='submitted',job_id='55');live.add('55')
    atomic_json(graph/'state.json',state)
    after=dispatch.advance(graph)
    assert after['tasks']['bad']['status']=='failed'
    assert 'downstream' in after['tasks']['bad']['allocation_retry_blocked']
    assert after['tasks']['child']['job_id']=='55'


def test_scontrol_outage_is_not_missing_evidence_or_optional_unavailability(site,monkeypatch):
    graph=graph_at(site);before=(graph/'state.json').read_bytes()
    def unavailable(*a,**kw):raise OSError('controller unavailable')
    monkeypatch.setattr(recovery,'slurm_job_record',unavailable)
    with pytest.raises(dispatch.SchedulerUnavailableError):dispatch.advance(graph)
    assert (graph/'state.json').read_bytes()==before


def test_transaction_interruption_retains_accepted_geometry_checkpoint(site,monkeypatch):
    graph=graph_at(site,geometry=True);root,now,live,records=site
    directory=graph/'tasks/bad';record=json.loads((directory/'task.json').read_text())
    atoms=read(directory/'geometry.xyz');original=atoms.positions.copy()
    task=record['task'];steps=task.get('optimizer',{}).get('steps',100)
    p=Progress(directory,record,atoms,used_chunks=0,chunk_steps=steps,max_chunks=3)
    atoms.positions[1,0]+=.01
    atoms.calc=SinglePointCalculator(atoms,energy=-1.,forces=np.zeros_like(atoms.positions))
    opt=optimizer(nsteps=2);p.attach(opt,atoms);opt.save()
    dispatch.advance(graph);now[0]+=10
    real=run_recover.atomic_json
    def fail(path,data):
        if Path(path)==graph/'state.json':raise OSError('crash after archival')
        return real(path,data)
    monkeypatch.setattr(run_recover,'atomic_json',fail)
    with pytest.raises(OSError,match='crash after'):dispatch.advance(graph)
    monkeypatch.setattr(run_recover,'atomic_json',real)
    state=dispatch.advance(graph);entry=state['tasks']['bad']
    assert entry['status']=='staged' and entry['attempt']==2
    loaded=read(directory/'geometry.xyz')
    assert loaded.positions==pytest.approx(original)
    resumed=Progress(directory,json.loads((directory/'task.json').read_text()),loaded,
                     used_chunks=0,chunk_steps=steps,max_chunks=3)
    assert resumed.steps==2 and loaded.positions==pytest.approx(atoms.positions)
    assert len(state['allocation_retries']['bad'])==1
    assert (Path(entry['previous_attempt'])/'.kinbot_checkpoints').is_dir()


def test_driver_keeps_waiting_until_due_and_handles_submission_outage(site,monkeypatch):
    graph=graph_at(site);root,now,live,records=site
    # Exercise drive_run rather than only the helper. After recovery, emulate
    # an admission outage; it must preserve the single prepared replacement.
    real=dispatch.advance;calls=[]
    def advance(path,**kw):
        state=real(path,submit=False)
        calls.append(state['tasks']['bad']['status'])
        if state['tasks']['bad']['status']=='staged':
            if calls.count('staged')==1:
                raise dispatch.SchedulerUnavailableError('admission service unavailable')
            _mock_execution(graph,'bad');_mock_execution(graph,'live');live.clear()
            return real(path,submit=False)
        return state
    monkeypatch.setattr(dispatch,'advance',advance)
    monkeypatch.setattr(dispatch.time,'sleep',lambda seconds:now.__setitem__(0,now[0]+seconds))
    code,summary=dispatch.drive_run(graph,interval=10)
    assert code==0 and set(summary.values())=={'complete'}
    assert calls==['retry_wait','staged','staged']
    assert len(dispatch._load(graph)[2]['allocation_retries']['bad'])==1


@pytest.mark.parametrize('options',[dict(allocation_retries=True),dict(allocation_retries=11),
    dict(allocation_retry_delay=0),dict(allocation_retry_delay=50,allocation_retry_max_delay=10)])
def test_invalid_allocation_policy_rejected(tmp_path,options):
    with pytest.raises(ValueError):run_policy.update(tmp_path,max_jobs=1,**options)


@pytest.mark.parametrize('reply,raises',[
    (SimpleNamespace(returncode=1,stdout='',stderr='slurm_load_jobs error: Invalid job id specified'),False),
    (SimpleNamespace(returncode=1,stdout='',stderr='Unable to contact slurm controller'),True),
    (SimpleNamespace(returncode=0,stdout='',stderr=''),True)])
def test_strict_terminal_query_distinguishes_expired_record_from_outage(monkeypatch,reply,raises):
    from kinbot import execution
    monkeypatch.setattr(execution.subprocess,'run',lambda *a,**k:reply)
    if raises:
        with pytest.raises(OSError):execution.slurm_job_record('42',strict=True)
    else:assert execution.slurm_job_record('42',strict=True)=={}
