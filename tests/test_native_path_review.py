"""Native-excerpt classification plus simulated-site owner orchestration.

Synthetic workers exercise real scratch/dispatch/receipt adoption. They are
not evidence that Gaussian's DVV or Molpro's numerical audit converges.
"""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import qc_scratch, run_policy, run_recover
from kinbot.anl import dispatch, live_control, path_review, protocol
from kinbot.execution import atomic_json, file_sha256
from tests.test_anl_dispatch import _mock_execution
from tests.test_gaussian_checkpoint_inputs import graph, old_failure, resume_panel
from tests.test_r8_validation_panel import load, site


def corrector(kind='plateau'):
    lines = []
    for i in range(1,101):
        dx = .000166 if kind == 'plateau' else .000012*(1.014**(100-i))
        lines.extend([f' Pt 31 Step number {i} out of a maximum of 100',
                      f' New End-Old End Dist. = {dx:.8f}',
                      ' Recorrection delta-x convergence threshold: 0.000010'])
    return '\n'.join(lines)+'\n Maximum number of corrector steps exceded.\n Error termination via Lnk1e\n'


def one_point():
    return ''' CURRENT STRUCTURE
 Cartesian Coordinates (Ang):
 1 8 0.0 0.0 0.0
 2 1 0.9 0.0 0.0
 3 1 4.0 0.0 0.0
 ---------------------------------------------------------------------
 # OF POINTS ALONG THE PATH = 1
 NET REACTION COORDINATE UP TO THIS POINT = 0.09991
 PES minimum detected on this side of the pathway.
 Magnitude of the gradient = 0.0001286
 Reaction path calculation complete.
 Normal termination of Gaussian
 Standard orientation:
 1 8 0 0.0 0.0 0.0
 2 1 0 0.9 0.0 0.0
 3 1 0 0.0 0.9 0.0
 ---------------------------------------------------------------------
 Optimization completed.
 Normal termination of Gaussian
'''


def failure(graph, kind='plateau'):
    old_failure(graph)
    directory = graph/'tasks/forward_10'
    (directory/'forward_10.log').write_text(corrector(kind))
    p = directory/'execution.json'; record = json.loads(p.read_text())
    record.update(failure_phase='execution', error='Native corrector integration did not converge.')
    atomic_json(p, record)


def retained_attempt(graph, monkeypatch, status):
    """Run the source through scratch publication, including imported files."""
    _, spec, state = dispatch._load(graph)
    task = next(t for t in spec['tasks'] if t['id'] == 'forward_10')
    executable = Path(task['command'][0])
    output = corrector() if status == 'failed' else one_point()
    executable.write_text(executable.read_text().replace(
        'print("Normal termination of Gaussian")', f'print({output!r})\nsys.exit({int(status == "failed")})'))
    state['tasks']['forward_10'].update(status='submitted', job_id='42')
    state['tasks']['reverse_10'].update(status='submitted', job_id='77')
    atomic_json(graph/'state.json', state)
    monkeypatch.setenv('SLURM_JOB_ID', '42')
    directory = graph/'tasks/forward_10'
    if status == 'failed':
        with pytest.raises(RuntimeError, match='exited with status 1'):
            dispatch.run_task(directory/'task.json')
    else:
        assert dispatch.run_task(directory/'task.json')['status'] == 'executed'
    assert len(list(directory.glob('kinbot_parent_*.chk'))) == 1
    assert (directory/'forward_10.com.native').is_file()
    dispatch.advance(graph)


def make_plan(resume, root, monkeypatch, recipe='dvv'):
    module = resume.reconciliation
    monkeypatch.setattr(module, 'IRC_PILOTS', (('example','forward_10',recipe),))
    monkeypatch.setattr(module, 'CC_AUDITS', ())
    plan = module.prepare(root)
    # Fixture uses a neutral graph name; the example decision's hold list is
    # otherwise the exact set of original independent IRC graphs.
    plan['held_irc_graphs'] = ['example']
    plan['sha256'] = path_review.digest({k:v for k,v in plan.items() if k!='sha256'})
    atomic_json(root/module.DECISION, plan)
    return plan


def test_native_stop_and_post_optimization_are_separate():
    report = path_review.gaussian_path(one_point().splitlines())
    assert report['classification'] == 'early_gradient_stop'
    assert report['points'] == 1 and report['gradient'] == .0001286
    assert report['possible_fragment_recapture']
    assert report['irc_components'] == [[0,1],[2]]
    assert report['optimized_components'] == [[0,1,2]]
    assert report['scientific_approval'] is False
    report = path_review.gaussian_path(one_point().replace('PATH = 1','PATH = 22').splitlines())
    assert report['classification'] == 'gradient_endpoint_review'
    # A normal termination without a native reason never becomes a minimum.
    assert path_review.gaussian_path(['Normal termination of Gaussian'])['classification'] == 'unclassified'


def test_contraction_is_required_before_a_cycle_extension(graph):
    failure(graph)
    with pytest.raises(ValueError, match='contraction'):path_review.irc_probe(graph, 'forward_10', 'pc150')
    probe = path_review.irc_probe(graph, 'forward_10', 'dvv')
    text = probe['task']['input_template']
    assert 'GradientOnly,DVV,MassWeighted' in text and ',VeryTight,' in text
    assert '--Link1--' not in text and 'Opt=' not in text and 'MaxCycle=100' in text
    assert 'StepSize=-10,MaxPoints=600' in text
    assert probe['task']['depends_on'] == ['hessian']
    failure(graph, 'contracting')
    probe = path_review.irc_probe(graph, 'forward_10', 'pc150')
    assert 'EulerPC' in probe['task']['input_template'] and 'MaxCycle=150,Tight' in probe['task']['input_template']
    with pytest.raises(ValueError, match='applicable'):path_review.irc_probe(graph, 'forward_10', 'dvv')


@pytest.mark.parametrize('defect', ['active','wrong_job','modified_input','modified_checkpoint','unknown_failure'])
def test_probes_refuse_unverified_or_inapplicable_sources(graph, monkeypatch, defect):
    failure(graph)
    directory = graph/'tasks/forward_10'
    if defect == 'active':monkeypatch.setattr(dispatch,'_job_active',lambda _:True)
    elif defect == 'wrong_job':
        p = directory/'execution.json'; record=json.loads(p.read_text());record['slurm_job_id']='44';atomic_json(p,record)
    elif defect == 'modified_input':(directory/'forward_10.com').write_text('changed')
    elif defect == 'modified_checkpoint':(graph/'tasks/hessian/hessian.chk').write_text('changed')
    else:(directory/'forward_10.log').write_text('Some unrelated error\n')
    with pytest.raises((ValueError,RuntimeError)):
        path_review.irc_probe(graph, 'forward_10', 'dvv')


def test_completed_branch_is_retained_as_source_not_reset(graph):
    directory = graph/'tasks/forward_10'
    (directory/'forward_10.log').write_text(one_point())
    _mock_execution(graph,'forward_10')
    _,_,state=dispatch._load(graph);state['tasks']['forward_10'].update(status='submitted',job_id='42');atomic_json(graph/'state.json',state)
    p=directory/'execution.json';record=json.loads(p.read_text());record['slurm_job_id']='42';atomic_json(p,record)
    probe=path_review.irc_probe(graph,'forward_10','dvv')
    assert probe['report']['classification']=='early_gradient_stop'
    assert record['status']=='executed'


@pytest.mark.parametrize('source_status', ['synthetic', 'failed', 'complete'])
def test_owner_applies_decision_preserves_live_work_and_adopts_probe(resume_panel, monkeypatch, source_status):
    resume, root, graph = resume_panel
    if source_status == 'synthetic': failure(graph)
    else: retained_attempt(graph, monkeypatch, source_status)
    source_state = dispatch._load(graph)[2]['tasks']['forward_10']['status']
    preserved={p:file_sha256(p) for p in (graph/'tasks').glob('*/*') if p.is_file()}
    plan=make_plan(resume,root,monkeypatch)
    assert preserved=={p:file_sha256(p) for p in preserved}
    monkeypatch.setattr(protocol,'preflight_ready',lambda _:True)
    submitted=[]
    def submit(command, **kwargs):
        kwargs['before_submit']();submitted.append(str(kwargs['cwd']))
        return SimpleNamespace(returncode=0,stdout='88',stderr='')
    monkeypatch.setattr(run_policy,'submit',submit)
    class Tick(Exception):pass
    monkeypatch.setattr(resume.time,'sleep',lambda _:(_ for _ in ()).throw(Tick()))
    with pytest.raises(Tick):resume.run(root)
    ident='forward_10_dvv_v1'; state=dispatch._load(graph)[2]
    assert state['tasks'][ident]['job_id']=='88'
    assert state['tasks']['forward_10']['status']==source_state
    assert state['tasks']['reverse_10']['job_id']=='77'
    assert preserved=={p:file_sha256(p) for p in preserved}
    assert len(list((graph/'attempts').glob('live_*')))==1
    assert json.loads((root/'reconciliation_applied.json').read_text())['sha256']==plan['sha256']
    # Repeated owner startup adopts the allocation without another trial.
    monkeypatch.setattr(dispatch,'_job_active',lambda j:str(j) in ('77','88'))
    with pytest.raises(Tick):resume.run(root)
    assert len(submitted)==1
    # A native receipt arriving while allocation is live must wait for exit.
    _mock_execution(graph,ident)
    p=graph/'tasks'/ident/'execution.json';record=json.loads(p.read_text());record['slurm_job_id']='88';atomic_json(p,record)
    with pytest.raises(Tick):resume.run(root)
    assert dispatch._load(graph)[2]['tasks'][ident]['status']=='submitted'
    monkeypatch.setattr(dispatch,'_job_active',lambda j:str(j)=='77')
    with pytest.raises(Tick):resume.run(root)
    assert dispatch._load(graph)[2]['tasks'][ident]['status']=='complete'
    assert len(submitted)==1 and not json.loads((root/'progress.json').read_text())['scientific_approval']
    assert resume.reconciliation.prepare(root)['sha256']==plan['sha256']


@pytest.mark.parametrize('source_status', ['synthetic', 'failed', 'complete'])
def test_real_probe_worker_stages_read_only_parent_checkpoint(graph,monkeypatch,source_status):
    executable=Path(dispatch._load(graph)[1]['tasks'][1]['command'][0])
    program=executable.read_text()
    if source_status == 'synthetic': failure(graph)
    else: retained_attempt(graph, monkeypatch, source_status)
    retained={p:file_sha256(p) for p in (graph/'tasks/forward_10').iterdir() if p.is_file()}
    probe=path_review.irc_probe(graph,'forward_10','dvv')
    _,spec,state=dispatch._load(graph);spec['tasks'].append(probe['task'])
    live_control.commit(graph,spec,state,set(),{})
    dispatch.advance(graph)
    executable.write_text(program.replace('IRC=(RCFC,GradientOnly,EulerPC,MassWeighted,Forward','IRC=(RCFC,GradientOnly,DVV,MassWeighted,Forward'))
    monkeypatch.setenv('SLURM_JOB_ID','88')
    ident=probe['task']['id'];p=graph/'tasks'/ident/'task.json'
    before=file_sha256(graph/'tasks/hessian/hessian.chk')
    outcome=dispatch.run_task(p)
    assert outcome['status']=='executed'
    assert outcome['details']['checkpoint_inputs'][0]['checkpoints'][0]['sha256']==before
    assert file_sha256(graph/'tasks/hessian/hessian.chk')==before
    assert retained=={p:file_sha256(p) for p in retained}


def test_failed_probe_exhausts_budget_and_held_sampling_does_not_launch(resume_panel,monkeypatch):
    resume,root,graph=resume_panel;failure(graph);plan=make_plan(resume,root,monkeypatch)
    resume.reconciliation.apply(root,plan);dispatch.advance(graph)
    ident=plan['irc'][0]['task']['id']
    _,_,state=dispatch._load(graph)
    # End the live sibling as well, allowing the finite review wave to finish.
    for name,job in [('reverse_10','77'),(ident,'88')]:
        state['tasks'][name].update(status='submitted',job_id=job)
        atomic_json(graph/'tasks'/name/'execution.json',dict(schema=1,task_id=name,
            geometry_sha256=state['tasks'][name]['geometry_sha256'],slurm_job_id=job,
            status='failed',failure_phase='execution',error='No convergence'))
    atomic_json(graph/'state.json',state)
    monkeypatch.setattr(dispatch,'_job_active',lambda _:False)
    monkeypatch.setattr(run_policy,'submit',lambda *a,**k:pytest.fail('exhausted or held settings submitted'))
    resume.run(root);resume.run(root)
    assert dispatch._load(graph)[2]['tasks'][ident]['status']=='failed'
    assert len(list((graph/'attempts').glob('live_*')))==1
    assert not json.loads((root/'completion.json').read_text())['scientific_approval']


def test_interrupted_decision_commit_resumes_once(resume_panel,monkeypatch):
    resume,root,graph=resume_panel;failure(graph);plan=make_plan(resume,root,monkeypatch)
    write=run_recover.atomic_json
    def interrupt(path,value):
        if Path(path)==graph/'state.json':raise OSError('publication interrupted')
        return write(path,value)
    monkeypatch.setattr(run_recover,'atomic_json',interrupt)
    with pytest.raises(OSError):resume.reconciliation.apply(root,plan)
    monkeypatch.setattr(run_recover,'atomic_json',write)
    live_control.resume(graph)
    resume.reconciliation.apply(root,plan);resume.reconciliation.apply(root,plan)
    assert len(list((graph/'attempts').glob('live_*')))==1
    assert dispatch._load(graph)[2]['tasks']['reverse_10']['job_id']=='77'


def test_changed_proof_prevents_any_decision_application(resume_panel,monkeypatch):
    resume,root,graph=resume_panel;failure(graph);plan=make_plan(resume,root,monkeypatch)
    (graph/'tasks/forward_10/forward_10.log').write_text('changed')
    with pytest.raises(ValueError,match='evidence changed'):resume.reconciliation.apply(root,plan)
    assert not (graph/'attempts').exists()


@pytest.mark.parametrize('race', ['none','submitting','staging','archival','unowned_job','scheduler_outage'])
def test_decision_refresh_checks_before_retiring_controller(resume_panel,monkeypatch,race):
    from contextlib import contextmanager, nullcontext
    import subprocess
    resume,root,graph=resume_panel;retained_attempt(graph,monkeypatch,'failed');plan=make_plan(resume,root,monkeypatch)
    if race=='submitting':
        _,_,state=dispatch._load(graph);state['tasks']['reverse_10']['status']='submitting';atomic_json(graph/'state.json',state)
    elif race=='staging':(graph/'tasks/unfinished').mkdir()
    elif race=='archival':atomic_json(graph/'.live_control_transaction.json',{'status':'applying'})
    records=[dict(pid=123,starttime='456')];events=[]
    monkeypatch.setattr(resume,'inventory',lambda _:records)
    @contextmanager
    def freeze(*a,**k):
        events.append('freeze')
        try:yield [(records[0],12)]
        finally:events.append('thaw')
    monkeypatch.setattr(resume.run_refresh,'_freeze',freeze)
    monkeypatch.setattr(resume.run_refresh,'_submission_guard',lambda _:nullcontext())
    monkeypatch.setattr(resume.run_refresh,'_pidfd_send_signal',lambda *a:events.append('retire'))
    monkeypatch.setattr(resume,'launch',lambda _:events.append('launch') or dict(pid=456))
    def jobs(_):
        if race=='scheduler_outage':raise subprocess.CalledProcessError(1,'squeue')
        return [dict(job_id='99' if race=='unowned_job' else '77',directory=str(graph/'tasks/reverse_10'))]
    monkeypatch.setattr(resume.run_cleanup,'_jobs',jobs)
    before=(graph/'state.json').read_bytes()
    if race=='none':
        assert resume.refresh(root)['pid']==456
        assert events==['freeze','retire','thaw','launch']
        atomic_json(root/'native_owner.json',dict(pid=123,signature=resume.code_signature(),decision=plan['sha256']))
        assert resume.refresh(root)['status']=='already_running'
        assert len(events)==4
    else:
        with pytest.raises((RuntimeError,ValueError,subprocess.CalledProcessError)):resume.refresh(root)
        assert events==['freeze','thaw']
    assert (graph/'state.json').read_bytes()==before


def test_owner_holds_unsent_sampling_but_adopts_allocated_sampling(resume_panel,monkeypatch):
    resume,root,graph=resume_panel
    _,spec,state=dispatch._load(graph)
    extra=deepcopy(spec['tasks'][1]);extra['id']='forward_5'
    extra['input_template']=extra['input_template'].replace('StepSize=-10','StepSize=-5')
    spec['tasks'].append(extra)
    live_control.commit(graph,spec,state,set(),{})
    monkeypatch.setattr(resume.campaign,'graphs',lambda _:{'example':deepcopy(spec)})
    atomic_json(graph/'requested.json',spec);dispatch.advance(graph)
    failure(graph);make_plan(resume,root,monkeypatch)
    monkeypatch.setattr(protocol,'preflight_ready',lambda _:True)
    submitted=[]
    def submit(command,**kwargs):
        kwargs['before_submit']();submitted.append(Path(kwargs['cwd']).name)
        return SimpleNamespace(returncode=0,stdout='88',stderr='')
    monkeypatch.setattr(run_policy,'submit',submit)
    class Tick(Exception):pass
    monkeypatch.setattr(resume.time,'sleep',lambda _:(_ for _ in ()).throw(Tick()))
    with pytest.raises(Tick):resume.run(root)
    assert submitted==['forward_10_dvv_v1']
    assert dispatch._load(graph)[2]['tasks']['forward_5']['status']=='staged'
    # Model a pre-existing allocation, retaining its task/geometry hashes.
    _,_,state=dispatch._load(graph);state['tasks']['forward_5'].update(status='submitted',job_id='99');atomic_json(graph/'state.json',state)
    _mock_execution(graph,'forward_5')
    p=graph/'tasks/forward_5/execution.json';record=json.loads(p.read_text());record['slurm_job_id']='99';atomic_json(p,record)
    monkeypatch.setattr(dispatch,'_job_active',lambda j:str(j) in ('77','88','99'))
    with pytest.raises(Tick):resume.run(root)
    assert dispatch._load(graph)[2]['tasks']['forward_5']['status']=='submitted'
    monkeypatch.setattr(dispatch,'_job_active',lambda j:str(j) in ('77','88'))
    with pytest.raises(Tick):resume.run(root)
    assert dispatch._load(graph)[2]['tasks']['forward_5']['status']=='complete'
    assert submitted==['forward_10_dvv_v1']


def test_dead_owner_setup_recovers_journal_before_refresh(resume_panel,monkeypatch):
    resume,root,graph=resume_panel;failure(graph);plan=make_plan(resume,root,monkeypatch)
    write=run_recover.atomic_json
    def interrupt(path,value):
        if Path(path)==graph/'state.json':raise OSError('interrupted')
        return write(path,value)
    monkeypatch.setattr(run_recover,'atomic_json',interrupt)
    with pytest.raises(OSError):resume.reconciliation.apply(root,plan)
    monkeypatch.setattr(run_recover,'atomic_json',write)
    monkeypatch.setattr(resume,'inventory',lambda _:[])
    resume.recover_setup(root)
    assert json.loads((graph/'.live_control_transaction.json').read_text())['status']=='applied'
    resume.reconciliation.apply(root,plan)
    assert len(list((graph/'attempts').glob('live_*')))==1


def test_native_imaginary_modes_are_not_parser_artifacts():
    text='Normal Modes\nWavenumbers [cm-1] 942.71 1006.92\nNormal Modes of imaginary frequencies\nWavenumbers [cm-1] 428.54 408.26 214.84 152.08 5.16\nNormal Modes of low/zero frequencies\nWavenumbers [cm-1] 0.0 0.0\n'
    assert path_review.molpro_modes(text)==[942.71,1006.92,-428.54,-408.26,-214.84,-152.08,-5.16]


def test_curvature_audit_uses_failed_geometry_without_promoting_it(tmp_path,monkeypatch):
    campaign=load('campaign');point=dict(symbols=['O','H','H'],positions=[[0.,0,0],[1.,0,0],[0,1.,0]],charge=1,multiplicity=2)
    spec=campaign.cc_refinement('example',point,True,True);site(spec)
    monkeypatch.setattr(dispatch,'assign_partitions',lambda _:None)
    monkeypatch.setattr(dispatch,'render_site_setup',lambda _:'# simulated site')
    monkeypatch.setattr(dispatch,'_job_active',lambda _:False)
    graph=tmp_path/'example';dispatch.prepare(spec,graph)
    directory=graph/'tasks/refine'
    (directory/'optimized.xyz').write_text('3\nRejected candidate\nO 0 0 0\nH 0.9 0 0\nH 0 0.9 0\n')
    (directory/'refine.out').write_text('Normal Modes\nWavenumbers [cm-1] 1000\nNormal Modes of imaginary frequencies\nWavenumbers [cm-1] 191.15 38.59\n')
    _,_,state=dispatch._load(graph);state['tasks']['refine'].update(status='failed',job_id='22');atomic_json(graph/'state.json',state)
    atomic_json(directory/'execution.json',dict(schema=1,task_id='refine',slurm_job_id='22',status='failed',
        geometry_sha256=state['tasks']['refine']['geometry_sha256'],error='Molpro transition state must contain exactly one imaginary mode.'))
    probe=path_review.curvature_probe(graph,'refine');new=probe['spec']
    assert new['molecule']['positions'][1][0]==.9
    assert spec['molecule']['positions'][1][0]==1.
    dispatch.validate_spec(new)
    assert len(new['tasks'])==3 and not new['intent']['scientific_approval']
    for task in new['tasks']:
        assert 'optg' not in task['input_template'] and 'nocheck' not in task['input_template']
        assert 'result_parser' not in task and 'geometry_output' not in task
    module=load('reconcile');root=tmp_path/'new';root.mkdir()
    module.stage_extra(root,'audit',new);module.stage_extra(root,'audit',new)
    assert dispatch._load(root/'audit')[1]['molecule']==new['molecule']
    assert dispatch._load(graph)[2]['tasks']['refine']['status']=='failed'
    # Interrupted preparation may be recovered only with our exact marker.
    partial=root/'partial';partial.mkdir();(partial/'workflow.json').write_text('{}')
    with pytest.raises(ValueError,match='Unowned'):module.stage_extra(root,'partial',new)
    atomic_json(root/'partial.preparing.json',new)
    module.stage_extra(root,'partial',new)
    assert len(list((root/'interrupted_preparations').glob('partial.*')))==1
    module.stage_extra(root,'partial',new)
    # A foreign electronic/core setting must never silently disappear.
    _,source_spec,state=dispatch._load(graph)
    source_spec['tasks'][0]['input_template']=source_spec['tasks'][0]['input_template'].replace('\nrhf\n','\ncore,0\nrhf\n')
    # The production builder is strict even before a new attempt is run.
    monkeypatch.setattr(path_review,'source',lambda *a:(source_spec,source_spec['tasks'][0],
        dict(error='Molpro transition state must contain exactly one imaginary mode.'),probe['proof']))
    with pytest.raises(ValueError,match='Unknown electronic'):path_review.curvature_probe(graph,'refine')


def test_read_only_watch_and_export_include_supplemental_evidence(tmp_path,monkeypatch,capsys):
    watcher=load('watch');collector=load('collect')
    root=tmp_path/'run'/'validation'/'native';root.mkdir(parents=True)
    atomic_json(root/'stage.json',dict(workflows=['irc_example']))
    atomic_json(root/'reconciliation_v1.json',dict(irc=[dict(graph='irc_example',task=dict(id='probe'))],curvature={'audit':{}}))
    atomic_json(root/'native_owner.json',dict(pid=0,status='running'))
    atomic_json(root/'progress.json',dict(execution_finished=False,held_for_review={'irc_example':['old']}))
    for name,ident,job in [('irc_example','probe','88'),('audit','hessian_005','89')]:
        graph=root/name;(graph/'tasks'/ident).mkdir(parents=True)
        atomic_json(graph/'workflow.json',dict(tasks=[dict(id=ident)]))
        atomic_json(graph/'state.json',dict(tasks={ident:dict(status='submitted',job_id=job)}))
        (graph/'tasks'/ident/(ident+'.out')).write_text('native output')
    before={p:p.read_bytes() for p in root.rglob('*') if p.is_file()}
    monkeypatch.setattr(watcher.subprocess,'run',lambda *a,**k:SimpleNamespace(stdout=f'88|RUNNING|1:02|{root}/irc_example/tasks/probe\n77|RUNNING|2:01|{root.parent.parent}/ch3s/ordinary\n'))
    watcher.snapshot(root)
    output=capsys.readouterr().out
    assert 'process absent/unverified' in output and 'ALL R8 allocations 2' in output
    assert 'hessian_005' in output and 'HELD unsubmitted tasks 1' in output
    selected=collector.selected(root.parent)
    assert root/'audit/tasks/hessian_005/hessian_005.out' in selected
    assert before=={p:p.read_bytes() for p in before}
    assert set(before)=={p for p in root.rglob('*') if p.is_file()}
