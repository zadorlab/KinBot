"""Queue-all handoff tests: orchestration and scheduler doubles, no native QC."""
from contextlib import contextmanager
import fcntl
import importlib.util
import os
from pathlib import Path
import signal
import sys
from types import SimpleNamespace

import pytest

from kinbot import molpro_reconnaissance as recon, run_policy, run_refresh
from kinbot.execution import atomic_json
from tests.test_ch2_sh_s_followup import tree, followup, screens, relaxation, legacy_metadata_change
from tests.test_requested_orbital_recon import completed_screen

sys.modules.setdefault('continue_ch2_sh_s_review',followup)
spec=importlib.util.spec_from_file_location('queue_review',Path(__file__).resolve().parents[1]/
    'examples/anl/tools/queue_progress_review.py')
queue=importlib.util.module_from_spec(spec);spec.loader.exec_module(queue)


def test_inherited_node_override_preserves_new_failures_and_native_deck(monkeypatch):
    monkeypatch.setenv('SBATCH_EXCLUDE','node16,environment_bad')
    calls=[]
    def expand(args,**kw):
        calls.append(args)
        assert args[:3]==['scontrol','show','hostnames']
        return SimpleNamespace(stdout={'node16,environment_bad':'node16\nenvironment_bad\n',
            'node[21,24,34],script_bad':'node21\nnode24\nnode34\nscript_bad\n'}[args[3]])
    monkeypatch.setattr(recon.subprocess,'run',expand)
    deck='#!/bin/bash\n#SBATCH --exclude=node[21,24,34],script_bad\nrun-native\n'
    result=recon._slurm_submission_options(deck,['new_bad','node16'],
        allow_inherited=['node16','node21','node24','node34'])
    assert result==['--exclude=environment_bad,new_bad,node16,script_bad']
    assert len(calls)==2


def test_all_fixed_probes_admit_together_and_keep_new_node_failures(tree,monkeypatch):
    production,review,destination,*_=tree
    followup.stage(production,review,destination)
    target=destination/'candidates'
    policy_before=(production/run_policy.NAME).read_bytes()
    for s in followup.read(target/'review_launch.json')['systems']:
        m=followup.read(s['manifest']);m['slurm_excluded_nodes']=['badnode','node16']
        atomic_json(s['manifest'],m)
    from kinbot import run_walltime
    monkeypatch.setattr(screens.subprocess,'run',lambda *a,**kw:SimpleNamespace(returncode=0))
    monkeypatch.setattr(run_walltime,'partition_walltime',lambda p:'31-00:00:00')
    calls=[]
    def run_many(paths,**kw):
        calls.append(kw)
        return [followup.read(p) for p in paths]
    monkeypatch.setattr(recon,'run_many',run_many)
    screens.run(target,queue_all=True,allow_nodes=['node16'],report_failures=True)
    assert calls[-1]['max_jobs']==24 and calls[-1]['startup_checked'] is True
    for s in followup.read(target/'review_launch.json')['systems']:
        m=followup.read(s['manifest']);assert m['slurm_excluded_nodes']==['badnode']
        m['slurm_excluded_nodes'].append('node16') # New allocation failure wins on resume.
        atomic_json(s['manifest'],m)
    screens.run(target,queue_all=True,allow_nodes=['node16'])
    assert all('node16' in followup.read(s['manifest'])['slurm_excluded_nodes']
               for s in followup.read(target/'review_launch.json')['systems'])
    assert (production/run_policy.NAME).read_bytes()==policy_before


def test_all_relaxed_graphs_adopt_live_and_repeat_without_restoring_old_exclusions(tree,monkeypatch):
    production,review,destination,*_=tree
    prepared=followup.stage(production,review,destination)
    output=destination/'relaxations'
    run_policy.update(output,exclude=['node16'])
    before=(production/run_policy.NAME).read_bytes()
    states={g:{'tasks':{'relaxed_scan':{'status':'running','job_id':str(40+i)}}}
            for i,g in enumerate(prepared['graphs'])}
    calls=[]
    monkeypatch.setattr(relaxation.dispatch,'preflight_until_scheduler_ready',lambda g:None)
    def advance(g,submit):
        policy=run_policy.read(output)
        assert policy['max_jobs']==10 and submit
        calls.append(states[g]['tasks']['relaxed_scan']['job_id'])
        return {'tasks':{'relaxed_scan':{'status':'complete'}}}
    monkeypatch.setattr(relaxation.dispatch,'advance',advance)
    monkeypatch.setattr(relaxation,'report',lambda p:{'panels':[]})
    relaxation.run(output,queue_all=True,allow_nodes=['node16'])
    assert len(calls)==10 and 'node16' not in run_policy.read(output)['excluded_nodes']
    run_policy.update(output,exclude=['node16'])
    followup.stage(production,review,destination)
    relaxation.run(output,queue_all=True,allow_nodes=['node16'])
    assert calls[:10]==calls[10:] and 'node16' in run_policy.read(output)['excluded_nodes']
    assert (production/run_policy.NAME).read_bytes()==before


def process(pid,script,args,cwd):
    return dict(pid=pid,ppid=1,starttime=str(pid*10),cwd=str(cwd),state='S',
                command=f'/python -B -u {queue.TOOLS/script} '+ ' '.join(map(str,args)))


def test_owner_scope_excludes_production_and_waits_for_scheduler_child(tmp_path,monkeypatch):
    production=tmp_path/'production';review=tmp_path/'review'
    original=process(21,'refresh_progress_20261024.sh',[production],review)
    waiting=process(22,'continue_ch2_sh_s_review.py',[production],tmp_path/'checkout')
    active=process(23,'vrc_relaxation_recon.py',['run',review/'relaxations'],review)
    unrelated=process(24,'vrc_relaxation_recon.py',['run',tmp_path/'other'],tmp_path/'other')
    production_process=dict(pid=25,ppid=1,starttime='50',cwd=str(production),command='python -m kinbot.pes input.json')
    rows=[original,waiting,active,unrelated,production_process]
    monkeypatch.setattr(queue.run_cleanup,'_processes',lambda p:rows)
    assert queue.inventory(production,review)==[original,waiting,active]
    rows.append(dict(pid=26,ppid=23,starttime='60',cwd=str(review),command='sbatch job',state='S'))
    with pytest.raises(run_refresh.RefreshBusyError):queue.inventory(production,review)


def test_freeze_callback_rechecks_review_and_resumes_on_child_race(tmp_path,monkeypatch):
    record=process(21,'vrc_relaxation_recon.py',['run',tmp_path],tmp_path)
    signals=[]
    monkeypatch.setattr(run_refresh,'_pidfd_open',lambda pid:91)
    monkeypatch.setattr(run_refresh,'_starttime',lambda pid:record['starttime'])
    monkeypatch.setattr(run_refresh,'_pidfd_send_signal',lambda fd,sig:signals.append((fd,sig)))
    monkeypatch.setattr(run_refresh,'_stopped',lambda pid:True)
    monkeypatch.setattr(run_refresh.os,'close',lambda fd:None)
    monkeypatch.setattr(run_refresh,'_controllers',lambda _:pytest.fail('Production inventory must not be used'))
    with pytest.raises(run_refresh.RefreshTransientError):
        with run_refresh._freeze(tmp_path,[record],inventory=lambda _:[record,dict(pid=22,starttime='220')]):
            pytest.fail('Race must abort handoff')
    assert signals==[(91,0),(91,signal.SIGSTOP),(91,signal.SIGCONT)]


def test_stop_retries_submission_boundary_and_preserves_jobs(tmp_path,monkeypatch):
    production=tmp_path/'production';review=tmp_path/'review';review.mkdir()
    record=process(21,'vrc_relaxation_recon.py',['run',review/'relaxations'],review)
    monkeypatch.setattr(queue,'inventory',lambda *a:[record])
    monkeypatch.setattr(run_refresh,'_probe_process_handles',lambda:None)
    events=[]
    @contextmanager
    def freeze(root,records,inventory):
        assert inventory(root)==records
        events.append('stop')
        try:yield [(record,91)]
        finally:events.append('resume')
    monkeypatch.setattr(run_refresh,'_freeze',freeze)
    attempts=[]
    def check(review):
        attempts.append(1)
        if len(attempts)==1:raise run_refresh.RefreshTransientError('submitting')
    monkeypatch.setattr(queue,'submission_check',check)
    monkeypatch.setattr(queue.time,'sleep',lambda _:None)
    jobs=[dict(job_id='4242',directory=str(review/'graph'),state='RUNNING')]
    monkeypatch.setattr(queue.run_cleanup,'_jobs',lambda _:jobs)
    monkeypatch.setattr(run_refresh,'_pidfd_send_signal',lambda fd,sig:events.append((fd,sig)))
    assert queue.stop_owners(production,review,review/'receipt.json')==jobs
    assert events==['stop','resume','stop',(91,signal.SIGKILL),'resume']
    assert queue.read(review/'receipt.json')['retained_jobs']==jobs


@pytest.mark.parametrize('replace_owner', [False, True])
def test_bad_preflight_and_duplicate_owner_never_stop_anything(tree,monkeypatch,replace_owner):
    production,review,destination,*_=tree
    monkeypatch.setattr(queue,'stop_owners',lambda *a:pytest.fail('No stopping before valid preflight'))
    def fail(*a):raise ValueError('seed changed')
    monkeypatch.setattr(queue,'preflight',fail)
    with pytest.raises(ValueError,match='seed changed'):
        queue.perform(production,review,Path('/pin'),['node16'],replace_owner=replace_owner)
    with (review/'.queue_all.lock').open('a') as f:
        fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with pytest.raises(RuntimeError,match='lock still held'):
            queue.perform(production,review,Path('/pin'),['node16'])


def test_supervisor_launches_every_phase_before_waiting_and_preserves_completed_evidence(tmp_path,monkeypatch):
    review=tmp_path/'review';review.mkdir();pin=tmp_path/'pin'
    atomic_json(review/'queue_all_refresh.json',{'phase':'owners_starting'})
    launches=[];polls=[]
    class Child:
        def __init__(self,command,**kw):
            assert not polls
            assert '--queue-all' in command and command[-2:]==['--allow-node','node16']
            assert kw['env']['PYTHONPATH']==str(pin)
            launches.append(command);self.pid=len(launches);self.returncode=None
        def poll(self):
            assert len(launches)==4
            polls.append(self.pid);self.returncode=0
            return 0
    for _,_,folder in queue.campaigns(review):folder.mkdir(parents=True)
    monkeypatch.setattr(queue.subprocess,'Popen',Child)
    monkeypatch.setattr(screens,'evidence',lambda _: 'saved native evidence')
    result=queue.supervise(review,pin,['node16'])
    assert result==dict(screens=0,relaxations=0,compact=0,compact_relaxations=0)
    assert queue.read(review/'queue_all_finished.json')['production_model_accepted'] is False


def test_pending_override_retains_other_exclusions_and_running_jobs(tmp_path,monkeypatch):
    review=tmp_path/'review';review.mkdir()
    jobs=[dict(job_id='1',directory=str(review/'a'),state='PENDING'),
          dict(job_id='2',directory=str(review/'b'),state='RUNNING')]
    monkeypatch.setattr(queue.run_cleanup,'_jobs',lambda _:jobs)
    calls=[]
    monkeypatch.setattr(queue,'slurm_job_record',lambda ident:dict(JobId=ident,JobState='PENDING',
        WorkDir=str(review/'a'),ExcNodeList='node16,bad' if not calls else 'bad'))
    def run(command,**kw):
        calls.append(command)
        return SimpleNamespace(stdout='node16\nbad\n' if command[-1]=='node16,bad' else 'bad\n',
                               returncode=0,stderr='')
    monkeypatch.setattr(queue.subprocess,'run',run)
    queue.expand_pending(review,['node16'],{})
    assert calls==[['scontrol','show','hostnames','node16,bad'],
                   ['scontrol','update','JobId=1','ExcNodeList=bad'],['scontrol','show','hostnames','bad']]


@pytest.mark.parametrize('replace_owner', [False, True])
def test_handoff_repeat_adopts_identical_job_ids_and_keeps_production_unchanged(tree,monkeypatch,replace_owner):
    production,review,destination,source,*_=tree
    before={p:p.read_bytes() for p in production.rglob('*') if p.is_file()}
    events=[]
    monkeypatch.setattr(queue,'preflight',lambda *a:events.append('preflight'))
    jobs=[dict(job_id='live-42',directory=str(review/'relaxations'),state='RUNNING')]
    def stop(prod,rev,marker,**kwargs):
        events.append('stop-review')
        atomic_json(marker,dict(retained_jobs=jobs))
        return jobs
    monkeypatch.setattr(queue,'stop_owners',stop)
    real=relaxation.prepare
    def original(prod,plan,out,**kwargs):
        if out==review/'relaxations':
            events.append('adopt-original');return {'graphs':[]}
        return real(prod,plan,out,**kwargs)
    monkeypatch.setattr(relaxation,'prepare',original)
    monkeypatch.setattr(queue,'expand_pending',lambda *a:None)
    submitted={}
    def supervise(rev,pin,allowed):
        events.append('all-phases')
        for system in queue.read(destination/'candidates/review_launch.json')['systems']:
            path=Path(system['manifest']);m=queue.read(path)
            for job in m['jobs']:
                key=str(path)+job['name']
                if key not in submitted:submitted[key]=str(len(submitted)+100)
                job.setdefault('job_id',submitted[key])
                assert job['job_id']==submitted[key]
            atomic_json(path,m)
        return {'retained':0}
    monkeypatch.setattr(queue,'supervise',supervise)
    for _ in range(2):queue.perform(production,review,Path('/pin'),['node16'],replace_owner=replace_owner)
    assert len(submitted)==24
    assert events==['preflight','stop-review','adopt-original','all-phases']*2
    assert queue.read(review/'queue_all_refresh.json')['retained_jobs']==jobs
    assert len(list((review/'.queue_all_history').glob('*.json')))==1
    assert all(p.read_bytes()==data for p,data in before.items())


def test_review_supervisor_identity_and_self_exclusion(tmp_path,monkeypatch):
    production=tmp_path/'production';review=tmp_path/'review'
    args=[production,'--allow-node','node16','--code',production/'.salvage_code.abc']
    previous=process(1001,'queue_progress_review.py',args,review)
    current=process(os.getpid(),'queue_progress_review.py',args+['--replace-owner'],review)
    other=process(1002,'queue_progress_review.py',[tmp_path/'other'],review)
    child=process(1003,'vrc_relaxation_recon.py',['run',review/'relaxations','--queue-all'],review)
    child['ppid']=1001
    rows=[previous,current,other,child]
    monkeypatch.setattr(queue.run_cleanup,'_processes',lambda _:rows)
    assert queue.inventory(production,review)==[previous,child]
    for extra in (['--unknown'], ['--code','/unrelated/pin'], ['--allow-node']):
        assert not queue.owner(process(1004,'queue_progress_review.py',[production,*extra],review),production,review)
    rows.append(dict(pid=1005,ppid=1001,starttime='10050',command='squeue',cwd=str(review)))
    with pytest.raises(run_refresh.RefreshBusyError):queue.inventory(production,review)


def test_explicit_replacement_releases_old_supervisor_lock_before_adoption(tree,monkeypatch):
    production,review,destination,*_=tree
    events=[]
    old=(review/'.queue_all.lock').open('a+')
    fcntl.flock(old,fcntl.LOCK_EX|fcntl.LOCK_NB)
    def preflight(*args):
        assert not old.closed
        events.append('preflight')
    def stop(prod,rev,marker,**kwargs):
        events.append('verified-stop')
        old.close()
        atomic_json(marker,dict(retained_jobs=[dict(job_id='42')]))
        return [dict(job_id='42')]
    monkeypatch.setattr(queue,'preflight',preflight)
    monkeypatch.setattr(queue,'stop_owners',stop)
    monkeypatch.setattr(relaxation,'prepare',lambda *args:events.append('adopt'))
    monkeypatch.setattr(followup,'stage',lambda *args:events.append('adopt'))
    monkeypatch.setattr(queue,'expand_pending',lambda *args:None)
    def supervise(*args):
        # The transition lock must not stay held for the new owner's lifetime:
        # another explicit replacement must be able to reach preflight later.
        with (review/'.queue_transition.lock').open('a+') as f:
            fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with (review/'.queue_all.lock').open('a+') as f:
            with pytest.raises(BlockingIOError):fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
        events.append('new-owner')
        return {'relaxations':0}
    monkeypatch.setattr(queue,'supervise',supervise)
    try:
        with pytest.raises(RuntimeError,match='lock still held'):
            queue.perform(production,review,Path('/pin'),[])
        assert events==[]
        queue.perform(production,review,Path('/pin'),[],replace_owner=True)
    finally:
        old.close()
    assert events==['preflight','verified-stop','adopt','adopt','new-owner']


def test_real_panel_retries_scratch_only_and_adopts_live_and_complete(tree,monkeypatch):
    from kinbot.anl import dispatch
    from tests.test_anl_dispatch import _mock_execution
    production,review,destination,*_=tree
    prepared=followup.stage(production,review,destination)
    output=destination/'relaxations'
    graphs=[Path(g) for g in prepared['graphs']]
    monkeypatch.setenv('KINBOT_QC_SCRATCH_RETRIES','2')
    _mock_execution(graphs[0],'relaxed_scan',geometry='final.xyz')
    unchanged={p:p.read_bytes() for base in (production,graphs[0]/'tasks',graphs[1]/'tasks')
               for p in base.rglob('*') if p.is_file()}
    science={g:(g/'workflow.json').read_bytes() for g in graphs}
    for i,g in enumerate(graphs[1:],1):
        s=dispatch._load(g)[2];entry=s['tasks']['relaxed_scan']
        entry.update(status='submitted' if i==1 else 'failed',job_id=str(100+i))
        atomic_json(g/'state.json',s)
        if i>1:
            atomic_json(g/'tasks/relaxed_scan/execution.json',dict(schema=1,task_id='relaxed_scan',
                status='failed',failure_phase='scratch_setup',retryable_scratch_node=True,
                geometry_sha256=entry['geometry_sha256'],slurm_job_id=str(100+i),
                execution_node='bad_scratch',error='No verified local QC scratch is usable: permission denied'))
    submitted=[]
    monkeypatch.setattr(dispatch,'_job_active',lambda j:j=='101' or j in submitted)
    monkeypatch.setattr(dispatch,'preflight_until_scheduler_ready',lambda g:None)
    def submit(command,*,cwd,**kwargs):
        assert 'bad_scratch' in run_policy.read(output)['excluded_nodes']
        assert '--exclude=' in (cwd/'job.slurm').read_text()
        assert 'bad_scratch' in (cwd/'job.slurm').read_text()
        ident=str(200+len(submitted));submitted.append(ident)
        return SimpleNamespace(returncode=0,stdout=ident,stderr='')
    monkeypatch.setattr(run_policy,'submit',submit)
    class PollFinished(Exception):pass
    def end_poll(_):raise PollFinished()
    monkeypatch.setattr(relaxation.time,'sleep',end_poll)
    for _ in range(2):
        with pytest.raises(PollFinished):
            relaxation.run(output,queue_all=True,allow_nodes=['bad_scratch','extra_node'])
    assert len(submitted)==8
    assert dispatch._load(graphs[0])[2]['tasks']['relaxed_scan']['status']=='complete'
    assert dispatch._load(graphs[1])[2]['tasks']['relaxed_scan']['job_id']=='101'
    for g in graphs[2:]:
        s=dispatch._load(g)[2];entry=s['tasks']['relaxed_scan']
        assert entry['status']=='submitted' and entry['attempt']==2
        assert len(s['scratch_node_retries']['relaxed_scan'])==1
        assert Path(entry['previous_attempt'],'execution.json').exists()
    assert all((g/'workflow.json').read_bytes()==data for g,data in science.items())
    assert all(p.read_bytes()==data for p,data in unchanged.items())


@pytest.fixture
def full_review(tree):
    production,review,destination,source,request,case=tree
    # The real two-channel plan, plus the real compact follow-up. Native
    # artifacts are synthetic; only orchestration/provenance is tested here.
    positions=[[1.2000100119,-.4655510032,.1013882575],[.9455572501,.7518250541,1.3204837046],
        [.2070132629,-.3355235642,2.8547670677],[-.134624981,1.3246843886,.7341794664],
        [.3345123258,-.9791964859,-.2861472967],[2.1494724455,-.9711253636,.1677215698]]
    for entry in queue.read(queue.PLAN/'relaxations.json')['panels']:
        recipe=queue.PLAN.parent/'r7_additional_20261021'/(Path(entry['source']).name+'.json')
        req=production/entry['request']
        atomic_json(req,dict(queue.read(recipe)['request_match'],schema=1,status='needs_review',positions=positions))
        completed_screen((req,case,production.parent/entry['source']),recipe)
    record=screens._saved(source)
    atomic_json(review/'screens/review_launch.json',dict(schema=1,files={},reviews=[str(source)],
        production_policy=record['production_policy'],systems=record['systems']))
    original=relaxation.prepare(production,queue.PLAN/'relaxations.json',review/'relaxations')
    compact=legacy_metadata_change(tree)
    # Complete the compact comparison once, retaining native result identities.
    for system in queue.read(destination/'candidates/review_launch.json')['systems']:
        path=Path(system['manifest']);m=queue.read(path)
        for i,j in enumerate(m['jobs']):j.update(status='complete',job_id=str(900+i))
        atomic_json(path,m)
    return tree,[Path(g) for g in [*original['graphs'],*compact['graphs']]]


def process_doubles(monkeypatch,production,review):
    rows=[process(1234,'vrc_relaxation_recon.py',['run',review/'relaxations','--queue-all'],review)]
    signals=[]
    monkeypatch.setattr(queue.run_cleanup,'_processes',lambda _:rows)
    monkeypatch.setattr(run_refresh,'_probe_process_handles',lambda:None)
    monkeypatch.setattr(run_refresh,'_pidfd_open',lambda pid:os.open(os.devnull,os.O_RDONLY))
    monkeypatch.setattr(run_refresh,'_starttime',lambda pid:str(pid*10))
    monkeypatch.setattr(run_refresh,'_stopped',lambda pid:True)
    def send(fd,sig):
        signals.append(sig)
        if sig==signal.SIGKILL:rows.clear()
    monkeypatch.setattr(run_refresh,'_pidfd_send_signal',send)
    return signals


def test_full_saved_review_handoff_recovers_legacy_metadata_and_fifteen_scratch_failures(full_review,monkeypatch):
    from kinbot.anl import dispatch
    from tests.test_anl_dispatch import _mock_execution
    tree,graphs=full_review
    production,review,destination,source,*_=tree
    _mock_execution(graphs[0],'relaxed_scan',geometry='final.xyz')
    for i,g in enumerate(graphs[1:],1):
        state=dispatch._load(g)[2];e=state['tasks']['relaxed_scan']
        e.update(status='submitted' if i<3 else 'failed',job_id=str(100+i))
        atomic_json(g/'state.json',state)
        if i>=3:
            atomic_json(g/'tasks/relaxed_scan/execution.json',dict(schema=1,task_id='relaxed_scan',
                status='failed',failure_phase='scratch_setup',retryable_scratch_node=True,
                geometry_sha256=e['geometry_sha256'],slurm_job_id=e['job_id'],execution_node='node34',
                error='No verified local QC scratch is usable: permission denied'))
    retained={p:p.read_bytes() for base in (production,destination/'candidates',graphs[0]/'tasks')
              for p in base.rglob('*') if p.is_file()}
    workflows={g:(g/'workflow.json').read_bytes() for g in graphs}
    native_before=(source/'native/reconnaissance.json').read_bytes()
    jobs=[dict(job_id=str(100+i),directory=str(graphs[i]),state='RUNNING') for i in (1,2)]
    monkeypatch.setattr(queue.run_cleanup,'_jobs',lambda _:jobs)
    signals=process_doubles(monkeypatch,production,review)
    checked=[]
    monkeypatch.setattr(dispatch,'preflight',lambda g:checked.append(g))
    monkeypatch.setattr(dispatch,'preflight_until_scheduler_ready',lambda g:None)
    monkeypatch.setenv('KINBOT_QC_SCRATCH_RETRIES','2')
    submitted=[]
    monkeypatch.setattr(dispatch,'_job_active',lambda j:j in ['101','102',*submitted])
    def submit(command,*,cwd,**kwargs):
        assert 'node34' in (cwd/'job.slurm').read_text()
        ident=str(500+len(submitted));submitted.append(ident)
        return SimpleNamespace(returncode=0,stdout=ident,stderr='')
    monkeypatch.setattr(run_policy,'submit',submit)
    class PollFinished(Exception):pass
    monkeypatch.setattr(relaxation.time,'sleep',lambda _:(_ for _ in ()).throw(PollFinished()))
    launches=[]
    class Child:
        def __init__(self,command,**kwargs):
            assert kwargs['env']['PYTHONPATH']=='/verified-new-pin'
            self.command=command;self.pid=len(launches)+1;self.returncode=None;launches.append(command)
        def poll(self):
            assert len(launches)%4==0
            if self.returncode is None and Path(self.command[3]).name=='vrc_relaxation_recon.py':
                with pytest.raises(PollFinished):
                    relaxation.run(Path(self.command[5]),queue_all=True,allow_nodes=['node16','node21','node24'])
            self.returncode=0
            return 0
    monkeypatch.setattr(queue.subprocess,'Popen',Child)
    monkeypatch.setattr(screens,'evidence',lambda p:'Synthetic orchestration test, not native acceptance.')
    for _ in range(2):
        assert not any(queue.perform(production,review,Path('/verified-new-pin'),
            ['node16','node21','node24'],replace_owner=True).values())
    assert signals.count(signal.SIGKILL)==1
    assert len(checked)==36 and len(submitted)==15 and len(launches)==8
    assert queue.read(review/'queue_all_refresh.json')['retained_jobs']==jobs
    assert all((g/'workflow.json').read_bytes()==data for g,data in workflows.items())
    assert all(p.read_bytes()==data for p,data in retained.items())
    assert (source/'native/reconnaissance.json').read_bytes()==native_before
    assert [dispatch._load(g)[2]['tasks']['relaxed_scan']['job_id'] for g in graphs[1:3]]==['101','102']
    for g in graphs[3:]:
        s=dispatch._load(g)[2];e=s['tasks']['relaxed_scan']
        assert e['attempt']==2 and len(s['scratch_node_retries']['relaxed_scan'])==1
        assert Path(e['previous_attempt'],'execution.json').exists()
    assert len(list(review.glob('**/probe_*/workflow.json')))==18


@pytest.mark.parametrize('during_freeze',[False,True])
def test_actual_panel_mismatch_aborts_before_kill_and_resumes_frozen_owner(full_review,monkeypatch,during_freeze):
    tree,graphs=full_review
    production,review,destination,source,request,case=tree
    monkeypatch.setattr(queue.run_cleanup,'_jobs',lambda _:[])
    monkeypatch.setattr(queue.dispatch,'preflight',lambda g:None)
    signals=process_doubles(monkeypatch,production,review)
    def change():
        data=queue.read(case);data['vrc_tst_scan_fmax']=.002;atomic_json(case,data)
    if during_freeze:
        original=run_refresh._freeze
        @contextmanager
        def frozen(*args,**kwargs):
            with original(*args,**kwargs) as handles:
                change()
                yield handles
        monkeypatch.setattr(run_refresh,'_freeze',frozen)
    else:change()
    monkeypatch.setattr(queue,'supervise',lambda *a:pytest.fail('Invalid panels must not launch'))
    with pytest.raises(ValueError,match='plan or evidence changed'):
        queue.perform(production,review,Path('/pin'),[],replace_owner=True)
    assert signal.SIGKILL not in signals
    assert signals==([0,signal.SIGSTOP,signal.SIGCONT] if during_freeze else [])


def test_staging_lock_race_never_freezes_writer(full_review,monkeypatch):
    tree,graphs=full_review
    production,review,*_=tree
    monkeypatch.setattr(queue.run_cleanup,'_jobs',lambda _:[])
    signals=process_doubles(monkeypatch,production,review)
    with (graphs[0]/'.manifest.lock').open('r') as writer:
        fcntl.flock(writer,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with pytest.raises(RuntimeError,match='stayed busy'):
            queue.stop_owners(production,review,review/'receipt.json',timeout=0,
                              validate_saved=lambda:pytest.fail('Writer still holds manifest'))
    assert signals==[]
