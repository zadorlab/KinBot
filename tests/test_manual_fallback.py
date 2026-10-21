from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess

import pytest

from kinbot import manual_fallback as fallback, live_control as control, qc_live, qc_timeout, run_policy
from kinbot.anl import dispatch, live_control as owner
from kinbot.execution import atomic_json, file_sha256
from tests.test_live_control import root, graph, allocation, request_path
from tests.test_timeout_recovery import pool, ordinary


def test_literal_worker_edits_preserve_identity_and_allow_removal():
    original="title='Å'; kwargs = {'method':'B3LYP','guess':'read','scf':'xqc'}\nsella_kwargs={'internal':True}\nsteps=100\n"
    options=dict(schema=1, worker_parameters={'kwargs':{'scf':'conver=8,CDIIS','guess':None},
                 'sella_kwargs':{'internal':False},'steps':300})
    changed=fallback.worker_options(original,options)
    namespace={};exec(changed,namespace)
    assert namespace['kwargs']=={'method':'B3LYP','scf':'conver=8,CDIIS'}
    assert namespace['steps']==300 and namespace['sella_kwargs']=={'internal':False}
    with pytest.raises(ValueError,match='no literal'):
        fallback.worker_options('pass\n',options)
    with pytest.raises(ValueError,match='Ambiguous'):
        fallback.worker_options(original+'\nsteps=20\n',options)


@pytest.mark.parametrize('options',[
    {'environment':{'USER':'someone'}}, {'environment':{'SLURM_JOB_ID':'4'}},
    {'environment':{'KINBOT_QC_ACCOUNT':'someone'}}, {'environment':{'PYTHONPATH':'other_code'}},
    {'calculator_kwargs':{'basis':'cc-pVDZ'}}, {'optimizer':{'order':0}},
    {'optimizer':{'sella_kwargs':{'order':0}}},
    {'resources':{'cores':0}}, {'policy':{'optimization_step_retries':-1}},
    {'input_edits':[{'old':'','new':'x'}]}, {'setup':['module load molpro\nscancel 4']},
])
def test_invalid_fallbacks_are_rejected(options):
    with pytest.raises(ValueError):fallback.validate(dict(schema=1,**options),'task')


def test_runtime_values_are_shell_quoted_and_setup_is_after_site(tmp_path):
    options=dict(schema=1,setup=['module load molpro'],environment={
        'TEST_VALUE':'literal $(touch unexpected); value', 'GAUSS_MEMDEF':None,
        'KINBOT_MOLPRO_EXECUTABLE':'/opt/molpro/bin/molpro'},policy={'optimization_step_retries':3})
    script=fallback.worker_script('#!/bin/bash\n#SBATCH --exclusive\nsource site.sh\npython job.py\n',
                                  'job.py',tmp_path/'retry worker.py',options)
    assert script.index('source site.sh')<script.index('module load molpro')<script.index('python ')
    assert "'literal $(touch unexpected); value'" in script
    assert 'unset GAUSS_MEMDEF' in script and 'KINBOT_OPTIMIZATION_STEP_RETRIES=3' in script
    p=tmp_path/'job.sh';p.write_text(script)
    subprocess.run(['bash','-n',str(p)],check=True)
    with pytest.raises(ValueError,match='uniquely'):
        fallback.worker_script('python job.py\npython job.py\n','job.py','x',options)


def test_task_retry_keeps_siblings_and_vrc_binding(root,monkeypatch):
    from kinbot.vrc_graph import scientific_spec
    g=graph(root);target=g/'tasks/change'
    original=dispatch._load(g)[1]
    state=dispatch._load(g)[2]
    state['tasks']['change'].update(status='failed',job_id='42',attempt=3)
    state['tasks']['live'].update(status='submitted',job_id='43')
    atomic_json(g/'state.json',state)
    preserved={p:file_sha256(p) for p in (g/'tasks/live').iterdir() if p.is_file()}
    monkeypatch.setattr(owner,'slurm_job_record',lambda job:allocation(target,job))
    monkeypatch.setattr(dispatch,'_job_active',lambda job:job=='43')
    options=dict(schema=1,input_edits=[dict(old='rhf\n',new='{rhf;maxit,300}\n')],
                 environment={'TMP_VALUE':'literal $HOME'},resources={'walltime':'7-00:00:00'})
    request=control.request(root,target,reason='Reviewed SCF fallback',fallback=options)
    owner.process(g);dispatch.advance(g)
    _,spec,state=dispatch._load(g)
    assert scientific_spec(spec)==scientific_spec(original)
    assert '{rhf;maxit,300}' in (target/'change.inp').read_text()
    assert state['tasks']['change']['attempt']==4 and state['tasks']['live']['job_id']=='43'
    assert preserved=={p:file_sha256(p) for p in preserved}
    assert json.loads(request_path(root,request).read_text())['status']=='applied'
    owner.process(g)
    assert len(list((g/'attempts').iterdir()))==1
    modified=deepcopy(spec);modified['tasks'][0]['input_template']+='! unapproved change\n'
    with pytest.raises(ValueError,match='receipt'):scientific_spec(modified)


def test_mismatched_edit_does_not_cancel_live_attempt(root,monkeypatch):
    g=graph(root,('change',));target=g/'tasks/change'
    state=dispatch._load(g)[2];state['tasks']['change'].update(status='submitted',job_id='42')
    atomic_json(g/'state.json',state)
    monkeypatch.setattr(owner,'stop_attempt',lambda *a,**k:pytest.fail('Cancelled before validation'))
    r=control.request(root,target,reason='mistyped keyword',cancel_running=True,
        fallback=dict(schema=1,input_edits=[dict(old='absent',new='changed')]))
    owner.process(g)
    assert json.loads(request_path(root,r).read_text())['status']=='blocked'
    assert not (g/'attempts').exists()


def test_qc_manual_retry_applies_runtime_and_keywords_once(ordinary,monkeypatch):
    qc,calls=ordinary;root=Path.cwd();qc.qc='gauss'
    control.configure(root)
    Path('well.py').write_text("kwargs={'scf':'xqc','method':'B3LYP','nprocshared':4}\nsteps=100\n")
    token=hashlib.sha256(b'well').hexdigest()[:20]
    guard=root/'.kinbot_slurm_io_retries'/token/'42/execution.json'
    saved=json.loads(guard.read_text())
    seed=guard.parent/'initial/well.py';seed.write_bytes(Path('well.py').read_bytes())
    saved['initial_artifacts']['well.py']=file_sha256(seed);atomic_json(guard,saved)
    control.observe_qc(root/'well','42','error')
    monkeypatch.setattr(owner,'slurm_job_record',lambda job:allocation(root,job))
    options=dict(schema=1,worker_parameters={'kwargs':{'scf':'conver=8,CDIIS'},'steps':400},
                 environment={'GAUSS_MEMDEF':None},resources={'cores':8,'walltime':'7-00:00:00'})
    r=control.request(root,'well',kind='qc',reason='manual SCF retry',fallback=options)
    assert qc_live.process(qc,'well','error')=='running'
    assert qc.job_ids['well']=='45' and len(calls)==1
    m=json.loads(qc_timeout._marker('well').read_text())
    worker=Path(m['attempts'][-1]['archive'])/'retry_worker.py'
    namespace={};exec(worker.read_text(),namespace)
    assert namespace['steps']==400 and namespace['kwargs']['nprocshared']==8
    assert namespace['kwargs']['scf']=='conver=8,CDIIS'
    assert 'unset GAUSS_MEMDEF' in Path(m['script']).read_text()
    assert Path('other.chk').read_text()=='sibling'
    assert json.loads(request_path(root,r).read_text())['status']=='applied'
    qc_timeout.resume(qc,'well');assert len(calls)==1


def test_task_only_failure_hold_does_not_freeze_reaction_searches(root):
    control.configure(root,hold_task_failures=True)
    assert control.holding(root/'graph/tasks/failed')
    assert not control.holding(root/'rejected',kind='qc')
    assert not control.holding(root/'rejected',kind='reaction')


def test_reaction_fallback_applies_only_to_its_new_workers(root,monkeypatch):
    monkeypatch.chdir(root)
    from types import SimpleNamespace
    value=dict(schema=1,worker_parameters={'kwargs':{'scf':'conver=8,CDIIS'}},
               environment={'GAUSS_MEMDEF':None},resources={'cores':2})
    fallback.install_reaction(root/'reaction',dict(id='reviewed',fallback=value))
    assert fallback.for_reaction_job(root/'reaction_high')==value
    assert fallback.for_reaction_job(root/'reaction2_high')=={}
    assert fallback.for_reaction_job(root/'shared_product_well')=={}
    Path('reaction_high.py').write_text("kwargs={'method':'B3LYP','nprocshared':4}\n")
    script=qc_live.prepare_resources(SimpleNamespace(qc='gauss'),'reaction_high',
        '#!/bin/bash\n#SBATCH --exclusive\npython reaction_high.py\n')
    namespace={};exec(Path('reaction_high.py').read_text(),namespace)
    assert namespace['kwargs']['scf']=='conver=8,CDIIS'
    assert namespace['kwargs']['nprocshared']==2
    assert 'unset GAUSS_MEMDEF' in script and '--exclusive' in script


def test_ase_retry_preserves_geometry_order_and_checks_command():
    from tests.anl_fixture import dispatch_spec
    task=dispatch_spec()['tasks'][0]
    changed=fallback.task_options(task,dict(schema=1,
        calculator_kwargs={'command':'/opt/molpro/bin/molpro','scf_convergence':1e-8},
        optimizer={'steps':300,'sella_kwargs':{'internal':False}}))
    assert changed['profile']['command']=='/opt/molpro/bin/molpro'
    assert changed['optimizer']['steps']==300
    assert changed['optimizer'].get('order')==task['optimizer'].get('order')
    assert fallback.logical_task(changed)==task


def test_successive_fallbacks_keep_resources_and_run_setup_only_in_allocation(root,monkeypatch):
    g=graph(root,('change',));task=dispatch._load(g)[1]['tasks'][0]
    first=fallback.task_options(task,dict(schema=1,resources={'walltime':'7-00:00:00'},
        setup=['echo allocation-only'],environment={'CUSTOM_VALUE':'quoted $value'}))
    second=fallback.task_options(first,dict(schema=1,input_edits=[dict(old='rhf\n',new='{rhf;maxit,300}\n')]))
    assert second['manual_fallback']['options']['resources']['walltime']=='7-00:00:00'
    assert 'echo allocation-only' not in second.get('setup',[])
    spec=dispatch._load(g)[1];spec['tasks']=[second]
    state=dispatch._load(g)[2];state['tasks']['change'].update(status='staged')
    owner.commit(g,spec,state,{'change'},{'change':'Fallback before submission'})
    dispatch.advance(g)
    script=(g/'tasks/change/job.slurm').read_text()
    assert script.index('echo allocation-only')>script.index('source ../../site_setup.sh')
    assert script.index('echo allocation-only')<script.index('-m kinbot.anl.dispatch run-task')
    import subprocess
    monkeypatch.setattr(dispatch.shutil,'which',lambda name:'/usr/bin/'+name)
    calls=[]
    def probe(args,**kwargs):
        calls.append(args)
        assert args[:2]==['sbatch','--test-only']
        return subprocess.CompletedProcess(args,0,'','')
    monkeypatch.setattr(dispatch.subprocess,'run',probe)
    report=dispatch.preflight(g)
    assert report['runtime_environment_deferred']==['change'] and len(calls)==1


def test_fallback_can_edit_unsent_task_without_scheduler_mutation(root,monkeypatch):
    g=graph(root,('change',));target=g/'tasks/change'
    monkeypatch.setattr(owner,'stop_attempt',lambda *a,**k:pytest.fail('Unsent job cancellation'))
    r=control.request(root,target,reason='Correct setup before first submission',
        fallback=dict(schema=1,environment={'NEW_VALUE':'data'}))
    owner.process(g);dispatch.advance(g)
    assert json.loads(request_path(root,r).read_text())['status']=='applied'
    assert 'export NEW_VALUE=data' in (target/'job.slurm').read_text()


def test_retry_preview_uses_effective_worker_in_nested_directory(root):
    directory=root/'well';directory.mkdir()
    target=directory/'conf/job';target.parent.mkdir()
    token=hashlib.sha256(b'conf/job').hexdigest()[:24]
    archive=directory/'.kinbot_timeout_retries'/token/'old'
    archive.mkdir(parents=True)
    effective=archive/'retry_worker.py';effective.write_text('steps=500\n')
    batch=archive/'retry.sbatch';batch.write_text('python '+str(effective)+'\n')
    atomic_json(archive.parent/'current.json',dict(replacement_ids=['43'],
        script=str(batch),attempts=[{'archive':str(archive)}]))
    assert qc_live.retry_sources(target,'43')==(effective,batch,str(effective))
