"""Live operations preserve allocation ownership and scientific provenance."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import live_control as control, pool_drain, qc_live, run_control, run_policy
from kinbot.anl import dispatch, live_control as owner
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec, molpro_task
from tests.test_anl_dispatch import _mock_execution
from tests.test_timeout_recovery import pool, ordinary


@pytest.fixture
def root(tmp_path):
    control.configure(tmp_path)
    return tmp_path


def graph(root, names=('change','live','done')):
    spec = dispatch_spec()
    spec['tasks'] = [molpro_task(n,'rhf\n',geometry_from='initial') for n in names]
    return dispatch.prepare(spec,root/'graph')


def allocation(directory, job='42', state='FAILED', priority='100'):
    return dict(JobId=job,WorkDir=str(directory),JobState=state,Priority=priority,
                TimeLimit='1-00:00:00',Partition='long')


def test_named_rules_scope_order_remove_history_and_legacy(tmp_path):
    assert control.overrides(tmp_path/'any') == {}
    assert not control.holding(tmp_path/'any')
    control.configure(tmp_path)
    control.resources(tmp_path,'all',['*'],{'cores':4,'memory_mb':20000})
    control.resources(tmp_path,'one',['graph/tasks/foo'],{'cores':8})
    assert control.overrides(tmp_path/'graph/tasks/foo') == {'cores':8,'memory_mb':20000}
    assert control.overrides(tmp_path/'other/foo')['cores'] == 4
    control.resources(tmp_path,'one',remove=True)
    assert control.overrides(tmp_path/'graph/tasks/foo')['cores'] == 4
    assert len(list((tmp_path/'.kinbot_controls/history').glob('*.json'))) == 4
    (tmp_path/'child').mkdir()
    with pytest.raises(ValueError,match='existing'):
        control.configure(tmp_path/'child')


@pytest.mark.parametrize('values',[{'cores':0},{'cores':True},{'memory_mb':-1},
                                   {'walltime':'00:00:00'},{'walltime':'forever'},{'method':'HF'}])
def test_invalid_resources_never_published(root,values):
    before=(root/control.NAME).read_bytes()
    with pytest.raises(ValueError): control.resources(root,'bad',['*'],values)
    assert (root/control.NAME).read_bytes()==before


def test_resource_restage_preserves_complete_and_running_native_files(root,monkeypatch):
    g=graph(root)
    _mock_execution(g,'done')
    dispatch.advance(g)
    state=dispatch._load(g)[2]
    state['tasks']['live'].update(status='submitted',job_id='42')
    atomic_json(g/'state.json',state)
    before={p:file_sha256(p) for n in ('live','done') for p in (g/'tasks'/n).iterdir() if p.is_file()}
    control.resources(root,'resources',['*'],{'cores':8,'memory_mb':32000,'walltime':'7-00:00:00'})
    monkeypatch.setattr(dispatch,'_job_active',lambda job:True)
    owner.process(g)
    state=dispatch.advance(g)
    assert state['tasks']['done']['status']=='complete'
    assert state['tasks']['live']['job_id']=='42'
    assert before=={p:file_sha256(p) for p in before}
    changed=json.loads((g/'tasks/change/task.json').read_text())['task']
    assert changed['resources']['cores']==8
    assert '--ntasks=8' in (g/'tasks/change/job.slurm').read_text()
    assert state['tasks']['change']['attempt']==1  # resource restage is not a QC attempt
    assert (Path(state['tasks']['change']['previous_attempt'])/'task.json').exists()
    snapshots=len(list((g/'attempts').iterdir()))
    dispatch.advance(g)
    assert len(list((g/'attempts').iterdir()))==snapshots
    control.resources(root,'resources',remove=True)
    owner.process(g)
    dispatch.advance(g)
    assert dispatch._load(g)[1]['tasks'][0]['resources']['cores']==4


def test_legacy_pending_worker_defers_manifest_change(root):
    g=graph(root)
    state=dispatch._load(g)[2]
    p=g/'tasks/live/task.json'; m=json.loads(p.read_text()); m.pop('live_control_schema')
    atomic_json(p,m)
    state['tasks']['live'].update(status='submitted',job_id='42',task_sha256=file_sha256(p))
    atomic_json(g/'state.json',state)
    before=(g/'workflow.json').read_bytes()
    control.resources(root,'new',['*/change'],{'cores':8,'memory_mb':32000})
    owner.process(g)
    assert (g/'workflow.json').read_bytes()==before
    assert json.loads((g/'.live_control_status.json').read_text())['status']=='deferred'


def test_invalid_native_resources_leave_valid_work_running(root):
    g=graph(root)
    before=(g/'workflow.json').read_bytes()
    control.resources(root,'new',['*'],{'cores':999})
    owner.process(g)
    assert (g/'workflow.json').read_bytes()==before
    assert json.loads((g/'.live_control_status.json').read_text())['status']=='blocked'


def test_transaction_interruption_allows_unchanged_worker_and_resumes_once(root,monkeypatch):
    from kinbot import run_recover
    g=graph(root)
    before=dispatch._load(g)
    control.resources(root,'one',['*/change'],{'cores':8,'memory_mb':32000})
    original=run_recover.atomic_json
    def fault(p,data):
        if p==g/'state.json': raise OSError('interrupted')
        return original(p,data)
    monkeypatch.setattr(run_recover,'atomic_json',fault)
    with pytest.raises(OSError): owner.process(g)
    loaded=owner.worker_load(g,'live')
    assert loaded[1:]==before[1:]
    with pytest.raises(RuntimeError,match='Retired'): owner.worker_load(g,'change')
    monkeypatch.setattr(run_recover,'atomic_json',original)
    owner.process(g)
    dispatch.advance(g)
    assert dispatch._load(g)[2]['tasks']['change']['attempt']==1
    assert len(list((g/'attempts').iterdir()))==1


def request_path(root,record):
    return root/'.kinbot_controls/requests'/(record['id']+'.json')


def test_manual_retry_exhausted_budget_is_one_attempt_and_stale_request_rejected(root,monkeypatch):
    g=graph(root,('change',)); t=g/'tasks/change'
    state=dispatch._load(g)[2]
    state['tasks']['change'].update(status='failed',job_id='42',attempt=4)
    state['timeout_retries']={'change':[1,2,3]}
    atomic_json(g/'state.json',state)
    monkeypatch.setattr(owner,'slurm_job_record',lambda job:allocation(t,job))
    r=control.request(root,t,reason='operator reviewed failure')
    owner.process(g)
    dispatch.advance(g)
    state=dispatch._load(g)[2]
    assert state['tasks']['change']['attempt']==5
    assert state['timeout_retries']['change']==[1,2,3]
    assert json.loads(request_path(root,r).read_text())['status']=='applied'
    owner.process(g)
    assert len(list((g/'attempts').iterdir()))==1
    r2=control.request(root,t,reason='stale click',expected=r['expected'])
    owner.process(g)
    assert json.loads(request_path(root,r2).read_text())['status']=='blocked'


@pytest.mark.parametrize('phase',['RUNNING','COMPLETED','UNKNOWN'])
def test_retry_never_duplicates_live_completed_or_unknown_allocation(root,monkeypatch,phase):
    g=graph(root,('change',)); t=g/'tasks/change'
    state=dispatch._load(g)[2]; state['tasks']['change'].update(status='failed',job_id='42')
    atomic_json(g/'state.json',state)
    monkeypatch.setattr(owner,'slurm_job_record',lambda job:{} if phase=='UNKNOWN' else allocation(t,job,phase))
    r=control.request(root,t,reason='check fencing')
    before=(g/'state.json').read_bytes()
    owner.process(g)
    assert (g/'state.json').read_bytes()==before
    assert not (g/'attempts').exists()
    assert json.loads(request_path(root,r).read_text())['status'] in ('blocked','retained','queued')


def test_release_failure_hold_is_attempt_specific(root):
    g=graph(root,('change',)); t=g/'tasks/change'
    control.configure(root,hold_failures=True)
    assert control.holding(t)
    r=control.request(root,t,action='release',reason='continue without this optional result')
    owner.process(g)
    assert not control.holding(t)
    state=dispatch._load(g)[2];state['tasks']['change']['attempt']=2
    atomic_json(g/'state.json',state)
    assert control.holding(t)


def test_pending_cancel_start_race_retains_running_job(root,monkeypatch):
    import kinbot.execution as execution
    from kinbot import run_refresh
    g=graph(root,('change',));t=g/'tasks/change'
    record=allocation(t,state='PENDING'); calls=[]
    monkeypatch.setattr(owner,'slurm_job_record',lambda job:dict(record))
    monkeypatch.setattr(execution,'slurm_job_record',lambda job:dict(record))
    def command(args):
        calls.append(args)
        if args[1]=='hold': record.update(JobState='RUNNING',Priority='0')
        if args[1]=='release': record.update(Priority='100')
        return ''
    monkeypatch.setattr(run_refresh,'_pending_request',command)
    r=control.request(root,t,pending_only=True,reason='change resource')
    assert not owner.stop_attempt(t,'42',request_path(root,r),r)
    assert all(c[0]!='scancel' for c in calls)
    assert record['Priority']=='100'
    assert r['status']=='retained'


def test_native_qc_resources_unicode_and_shell_are_consistent():
    source="title = 'Å'; kwargs = {'method':'B2PLYP','scf':'conver=8','nprocshared':4}\n"
    edited=qc_live.worker_resources(source,'gauss',{'cores':8,'memory_mb':32000})
    env={};exec(edited,env)
    assert env['title']=='Å' and env['kwargs']['scf']=='conver=8'
    assert env['kwargs']['nprocshared']==8 and env['kwargs']['mem']=='22400MB'
    script=qc_live.script_resources('#!/bin/bash\n#SBATCH --exclusive\n#SBATCH -c4\n#SBATCH --mem=0\n#SBATCH -t 01:00:00\nexport OMP_NUM_THREADS=4\npython job.py\n',
        {'cores':8,'memory_mb':32000,'walltime':'7-00:00:00'})
    assert '--mem=32000M' in script and 'OMP_NUM_THREADS=4' not in script
    assert script.index('--time=7')<script.index('export')
    assert '--exclusive' in script
    with pytest.raises(ValueError,match='one SBATCH'):
        qc_live.script_resources('#SBATCH --cpus-per-task=4 --exclusive\ntrue\n',{'cores':8})


def test_pending_time_can_decrease_and_preserves_user_holds(root,monkeypatch):
    records={'42':allocation(root,state='PENDING',priority='0')}
    monkeypatch.setattr(owner,'slurm_job_record',lambda job:dict(records[job]))
    plan=dict(job_id='42',directory=str(root),before='1-00:00:00',requested='12:00:00',status='reviewed')
    calls=[]
    def command(args):
        calls.append(args)
        if args[1]=='update':records['42']['TimeLimit']='12:00:00'
    run_control._pending_time_one(plan,root/'receipt.json',command)
    assert plan['status']=='applied'
    assert [c[1] for c in calls]==['update']
    run_control._pending_time_one(plan,root/'receipt.json',command)
    assert len(calls)==1


def test_time_update_denial_releases_our_hold_and_is_reported(root,monkeypatch):
    record=allocation(root,state='PENDING')
    monkeypatch.setattr(owner,'slurm_job_record',lambda job:dict(record))
    def command(args):
        if args[1]=='hold':record['Priority']='0'
        if args[1]=='release':record['Priority']='100'
        return 'access denied'
    plan=dict(job_id='42',directory=str(root),before='1-00:00:00',requested='2-00:00:00',status='reviewed')
    run_control._pending_time_one(plan,root/'receipt.json',command)
    assert plan['status']=='unconfirmed' and record['Priority']=='100'


def test_pending_time_uncertain_ownership_keeps_visible_resumable_receipt(root,monkeypatch):
    from kinbot import run_walltime
    monkeypatch.setattr(run_walltime,'jobs',lambda r:[dict(
        id='42',state='PENDING',limit='1-00:00:00',partition='long',directory=str(root))])
    monkeypatch.setattr(run_walltime,'partition_limits',lambda:{'long':31*86400})
    monkeypatch.setattr(owner,'slurm_job_record',lambda job:{})
    plans=run_control.pending_time(root,['42'],'+02:00:00',delta=True,apply=True)
    assert plans[0]['status']=='unconfirmed'
    assert Path(plans[0]['receipt']).is_file()
    shown=control.status(root)['scheduler_updates']
    assert len(shown)==1 and shown[0]['requested']==plans[0]['requested']
    record=allocation(root,state='PENDING')
    monkeypatch.setattr(owner,'slurm_job_record',lambda job:dict(record))
    calls=[]
    def command(args):
        calls.append(args)
        if args[1]=='hold': record['Priority']='0'
        elif args[1]=='release': record['Priority']='100'
        else: record['TimeLimit']=plans[0]['requested']
    run_control._pending_time_one(plans[0],Path(plans[0]['receipt']),command)
    assert plans[0]['status']=='applied'
    assert run_walltime.seconds(record['TimeLimit'])==26*3600
    run_control._pending_time_one(plans[0],Path(plans[0]['receipt']),command)
    assert sum(args[1]=='update' for args in calls)==1


def test_pool_lower_drain_raise_release_and_never_release_user_hold(root,monkeypatch):
    from kinbot import execution,run_refresh
    run_policy.update(root,max_jobs=1,manage_pending=True)
    records={str(i):allocation(root,str(i),'RUNNING' if i==1 else 'PENDING',
                               priority='0' if i==4 else '100') for i in range(1,5)}
    monkeypatch.setattr(execution,'slurm_job_record',lambda j:dict(records[j]))
    monkeypatch.setattr(owner,'slurm_job_record',lambda j:dict(records[j]))
    calls=[]
    def command(args):
        calls.append(args)
        records[args[-1]]['Priority']='0' if args[1]=='hold' else '100'
        return ''
    monkeypatch.setattr(run_refresh,'_pending_request',command)
    path=root/run_policy.NAME;jobs={j:1 for j in records}
    pool_drain.reconcile(path,run_policy.read(root),jobs,{})
    assert records['1']['Priority']=='100'
    assert records['2']['Priority']==records['3']['Priority']=='0'
    run_policy.update(root,max_jobs=2)
    pool_drain.reconcile(path,run_policy.read(root),jobs,{})
    assert records['2']['Priority']=='100' and records['3']['Priority']=='0'
    run_policy.update(root,max_jobs=3)
    pool_drain.reconcile(path,run_policy.read(root),jobs,{})
    assert records['3']['Priority']=='100' and records['4']['Priority']=='0'
    assert not any(c[-1] in ('1','4') for c in calls)


def test_input_defaults_do_not_overwrite_later_live_edits(tmp_path):
    options={'hold_failures':True,'resource_rules':[
        dict(name='all',selectors=['*'],resources={'walltime':'7-00:00:00'})]}
    control.initialize(tmp_path,options)
    control.resources(tmp_path,'all',['*'],{'walltime':'10-00:00:00'})
    control.initialize(tmp_path,options)
    assert control.overrides(tmp_path/'job')['walltime']=='10-00:00:00'
    with pytest.raises(ValueError):control.validate_options({'hold_failures':1})


def test_failed_search_gate_closes_before_handoff_and_races_queue(root):
    target=root/'reaction'
    with pytest.raises(ValueError,match='registered'):
        control.request(root,target,kind='reaction',reason='retry')
    control.search_scope(target,'reaction',-999)
    r=control.request(root,target,kind='reaction',reason='retry')
    assert not control.search_scope(target,'reaction',-999,close=True)
    control.acknowledge(request_path(root,r),r,'blocked',error='review')
    assert control.search_scope(target,'reaction',-999,close=True)
    with pytest.raises(ValueError,match='handed off'):
        control.request(root,target,kind='reaction',reason='too late')


def test_manual_qc_retry_has_own_budget_resources_and_idempotent_receipt(ordinary,monkeypatch):
    import hashlib
    from kinbot import qc_timeout
    qc,calls=ordinary;qc.qc='gauss'
    root=Path.cwd();control.configure(root)
    worker=root/'well.py';worker.write_text("kwargs = {'method':'B3LYP','nprocshared':4,'scf':'conver=8'}\n")
    token=hashlib.sha256(b'well').hexdigest()[:20]
    guard=root/'.kinbot_slurm_io_retries'/token/'42/execution.json'
    m=json.loads(guard.read_text());(guard.parent/'initial/well.py').write_bytes(worker.read_bytes())
    m['initial_artifacts']['well.py']=file_sha256(worker);atomic_json(guard,m)
    monkeypatch.setattr(owner,'slurm_job_record',lambda job:allocation(root,job))
    control.observe_qc(root/'well','42','error')
    control.resources(root,'edit',['well'],{'cores':8,'memory_mb':32000,'walltime':'3-00:00:00'})
    r=control.request(root,'well',kind='qc',reason='manual retry after exhausted automatic limit')
    assert qc_live.process(qc,'well','error')=='running'
    assert qc.job_ids['well']=='45' and len(calls)==1
    receipt=json.loads(request_path(root,r).read_text())
    assert receipt['status']=='applied' and receipt['job_id']=='45'
    record=json.loads(qc_timeout._marker('well').read_text())
    retry=Path(record['script']).read_text()
    assert '--cpus-per-task=8' in retry and '--time=3-00:00:00' in retry
    updated=(Path(receipt['archive'])/'retry_worker.py').read_text()
    assert "'nprocshared': 8" in updated and "'scf': 'conver=8'" in updated
    assert (root/'well.chk').read_bytes()==b'pre-allocation checkpoint'
    assert record['attempts'][-1]['source']=='operator'
    qc_live.process(qc,'well','running')
    assert len(calls)==1


def test_reaction_retry_keeps_shared_products_and_longer_reaction_names(root,monkeypatch):
    from kinbot import search_control
    from ase import Atoms
    from ase.db import connect
    monkeypatch.chdir(root)
    db=connect(str(root/'kinbot.db'))
    names=['rxn','rxn_long','product_well']
    for name in names:
        (root/(name+'.py')).write_text('kwargs = {}\n')
        (root/(name+'.log')).write_text('old native file')
        db.write(Atoms('H'),name=name,data={'status':'error'})
    qc=SimpleNamespace(db=db,job_ids={'rxn':'42'},job_generation_tagged={})
    from kinbot import qc as qc_module
    monkeypatch.setattr(qc_module,'_slurm_job_active',lambda j:False)
    control.search_scope(root/'rxn','reaction',-999)
    r=control.request(root,'rxn',kind='reaction',reason='try again')
    search_control._retire_reaction(qc,'rxn',{'rxn','rxn_long'},request_path(root,r),r)
    assert not (root/'rxn.py').exists()
    assert (root/'rxn_long.py').exists() and (root/'product_well.log').exists()
    assert list(db.select(name='rxn'))[-1].data['status']==0
    assert list(db.select(name='product_well'))[-1].data['status']=='error'
    count=len(list(db.select(name='rxn')))
    search_control._retire_reaction(qc,'rxn',{'rxn','rxn_long'},request_path(root,r),r)
    assert len(list(db.select(name='rxn')))==count


def test_rotd_resources_freeze_without_editing_input_pickle(root,monkeypatch):
    from kinbot.batch_live import callback
    script=root/'surf0_face0_samp0.sh'
    script.write_text('#!/bin/bash\n#SBATCH --ntasks=4\n#SBATCH --mem=32000M\n'
                      'export OMP_NUM_THREADS=1\npython -m kinbot.rotdpy_worker surf0_face0_samp0\n')
    payload=root/'surf0_face0_samp0.pkl';payload.write_bytes(b'immutable pending sample')
    control.resources(root,'sampling',['surf*'],{'cores':8,'memory_mb':64000})
    prepare=callback(script,kind='rotd');prepare()
    text=script.read_text()
    assert '#SBATCH --ntasks=8' in text and '#SBATCH --cpus-per-task=1' in text
    assert 'KINBOT_ROTD_RESOURCES=' in text and 'OMP_NUM_THREADS=1' in text
    assert payload.read_bytes()==b'immutable pending sample'
    control.resources(root,'sampling',remove=True);prepare()
    assert script.read_text().find('KINBOT_ROTD_RESOURCES')==-1


def test_admission_refresh_runs_only_after_capacity_is_available(root,monkeypatch):
    run_policy.update(root,max_jobs=1)
    monkeypatch.setattr(run_policy,'_queue',lambda r:{'42':1})
    calls=[]
    assert run_policy.submit(['sbatch','job.slurm'],cwd=root,before_submit=lambda:calls.append('render')) is None
    assert calls==[]
    monkeypatch.setattr(run_policy,'_queue',lambda r:{})
    monkeypatch.setattr(run_policy.subprocess,'run',lambda *a,**kw:SimpleNamespace(returncode=0,stdout='Submitted batch job 43',stderr=''))
    run_policy.submit(['sbatch','job.slurm'],cwd=root,before_submit=lambda:calls.append('render'))
    assert calls==['render']


def test_completed_during_explicit_stop_is_retained(root,monkeypatch):
    from kinbot import run_refresh
    g=graph(root,('change',));t=g/'tasks/change'
    records=iter([allocation(t,state='RUNNING'),allocation(t,state='COMPLETED')])
    monkeypatch.setattr(owner,'slurm_job_record',lambda j:next(records))
    monkeypatch.setattr(run_refresh,'_pending_request',lambda args:'')
    r=control.request(root,t,reason='restart',cancel_running=True)
    assert not owner.stop_attempt(t,'42',request_path(root,r),r)
    assert r['status']=='retained'


def test_readonly_reconciliation_does_not_apply_operator_requests(root,monkeypatch):
    g=graph(root,('change',)); t=g/'tasks/change'
    r=control.request(root,t,reason='cancel',cancel_running=True)
    monkeypatch.setattr(owner,'process',lambda *a:pytest.fail('read-only reconciliation consumed request'))
    dispatch.refresh_status(g)
    assert json.loads(request_path(root,r).read_text())['status']=='queued'


def test_resource_edit_does_not_shorten_automatic_timeout_escalation(root):
    g=graph(root,('change',))
    control.resources(root,'time',['*'],{'walltime':'7-00:00:00'})
    owner.process(g); dispatch.advance(g)
    state=dispatch._load(g)[2]
    state['timeout_retries']={'change':[{'walltime':'14-00:00:00','partition':'long'}]}
    atomic_json(g/'state.json',state)
    owner.process(g); dispatch.advance(g)
    assert dispatch._load(g)[1]['tasks'][0]['resources']['walltime']=='14-00:00:00'
    control.resources(root,'time',remove=True)
    owner.process(g)
    assert dispatch._load(g)[1]['tasks'][0]['resources']['walltime']=='14-00:00:00'


def test_failed_hold_prevents_fallback_before_operator_retry_or_release(root):
    control.configure(root,hold_failures=True)
    spec=dispatch_spec(); spec['tasks']=[
        molpro_task('primary','rhf\n',geometry_from='initial'),
        molpro_task('fallback','rhf\n',geometry_from='initial')]
    spec['tasks'][1]['fallback_on']=['primary']
    g=dispatch.prepare(spec,root/'graph')
    state=dispatch._load(g)[2]
    state['tasks']['primary'].update(status='unavailable',failure_phase='execution')
    atomic_json(g/'state.json',state)
    state=dispatch.advance(g)
    assert state['tasks']['fallback']['status']=='blocked'
    control.request(root,g/'tasks/primary',action='release',reason='allow fallback')
    owner.process(g)
    state=dispatch.advance(g)
    assert state['tasks']['fallback']['status']=='staged'


def test_qc_cancel_poll_accepts_new_exit_marker_for_same_allocation(ordinary,monkeypatch):
    qc,calls=ordinary; qc.qc='gauss'
    root=Path.cwd();control.configure(root)
    control.observe_qc(root/'well','42','running')
    r=control.request(root,'well',kind='qc',reason='explicit retry',cancel_running=True)
    (root/'well.exitcode').write_text('42 143')
    monkeypatch.setattr(owner,'slurm_job_record',lambda job:allocation(root,job,state='CANCELLED'))
    assert qc_live.process(qc,'well','error')=='running'
    assert json.loads(request_path(root,r).read_text())['status']=='applied'
    assert len(calls)==1


def test_second_manual_qc_retry_updates_the_actual_previous_worker(ordinary,monkeypatch):
    import hashlib
    import shutil
    from kinbot import qc_timeout
    qc,_=ordinary;qc.qc='gauss'
    root=Path.cwd();control.configure(root)
    (root/'well.py').write_text("kwargs = {'method':'B3LYP','nprocshared':4}\n")
    token=hashlib.sha256(b'well').hexdigest()[:20]
    seed=root/'.kinbot_slurm_io_retries'/token/'42'
    guard=json.loads((seed/'execution.json').read_text())
    (seed/'initial/well.py').write_bytes((root/'well.py').read_bytes())
    guard['initial_artifacts']['well.py']=file_sha256(root/'well.py')
    atomic_json(seed/'execution.json',guard)
    monkeypatch.setattr(owner,'slurm_job_record',lambda job:allocation(root,job))
    submitted=[]
    def submit(command,**kwargs):
        submitted.append(command)
        return SimpleNamespace(returncode=0,stdout=f'Submitted batch job {44+len(submitted)}\n',stderr='')
    monkeypatch.setattr(run_policy,'submit',submit)
    control.observe_qc(root/'well','42','error')
    control.resources(root,'edit',['well'],{'cores':8,'memory_mb':32000})
    control.request(root,'well',kind='qc',reason='first retry')
    qc_live.process(qc,'well','error')
    assert qc.job_ids['well']=='45'
    shutil.copytree(seed,seed.with_name('45'))
    guard['slurm_job_id']='45';atomic_json(seed.with_name('45')/'execution.json',guard)
    control.observe_qc(root/'well','45','error')
    control.resources(root,'edit',['well'],{'cores':6,'memory_mb':24000})
    r=control.request(root,'well',kind='qc',reason='second reviewed retry')
    qc_live.process(qc,'well','error')
    assert qc.job_ids['well']=='46' and len(submitted)==2
    state=json.loads(qc_timeout._marker('well').read_text())
    worker=Path(state['attempts'][-1]['archive'])/'retry_worker.py'
    assert "'nprocshared': 6" in worker.read_text()
    assert str(worker) in Path(state['script']).read_text()
    assert json.loads(request_path(root,r).read_text())['status']=='applied'
    qc_live.process(qc,'well','running')
    assert len(submitted)==2
