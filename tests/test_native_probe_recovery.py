"""Evidence-specific recovery plus full owner/worker lifecycle simulations."""
from copy import deepcopy
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from kinbot import run_policy
from kinbot.anl import dispatch, live_control, path_review, probe_recovery, protocol
from kinbot.anl.tasks import molpro_task
from kinbot.execution import atomic_json, file_sha256
from tests.test_anl_dispatch import _mock_execution
from tests.test_gaussian_checkpoint_inputs import graph, resume_panel
from tests.test_native_path_review import failure, make_plan
from tests.test_r8_validation_panel import load, site


INVALID = '''Integration scheme = DVV
RxnCyc: Step size is -0.10000.
RxnCyc: Invalid step size found!
Error termination via Lnk1e in l123.exe
'''


def legacy_failure(resume, root, graph, monkeypatch):
    failure(graph)
    plan = make_plan(resume, root, monkeypatch)
    task = plan['irc'][0]['task']
    task['input_template'] = task['input_template'].replace('StepSize=10', 'StepSize=-10')
    task['native_path_probe'].pop('step_units')
    plan['sha256'] = path_review.digest({k:v for k,v in plan.items() if k!='sha256'})
    atomic_json(root/resume.reconciliation.DECISION, plan)
    resume.reconciliation.apply(root, plan); dispatch.advance(graph)
    ident = task['id']; directory = graph/'tasks'/ident
    executable = Path(task['command'][0])
    executable.write_text(executable.read_text().replace(
        'IRC=(RCFC,GradientOnly,EulerPC', 'IRC=(RCFC,GradientOnly,DVV').replace(
        'print("Normal termination of Gaussian")',
        f'assert "StepSize=-10" in text\nprint({INVALID!r})\nsys.exit(1)'))
    _,_,state = dispatch._load(graph)
    state['tasks'][ident].update(status='submitted', job_id='88')
    atomic_json(graph/'state.json', state); monkeypatch.setenv('SLURM_JOB_ID','88')
    with pytest.raises(RuntimeError): dispatch.run_task(directory/'task.json')
    dispatch.advance(graph)
    return plan, ident


def repaired_plan(resume, root, graph, monkeypatch):
    plan, ident = legacy_failure(resume, root, graph, monkeypatch)
    before = (root/resume.reconciliation.DECISION).read_bytes()
    merged = resume.reconciliation.prepare_repairs(root)
    assert (root/resume.reconciliation.DECISION).read_bytes() == before
    assert len(merged['irc'])==2 and merged['sha256'] != plan['sha256']
    assert resume.reconciliation.prepare_repairs(root) == merged
    return merged


def test_repaired_dvv_actual_worker_and_owner_adoption(resume_panel, monkeypatch):
    resume,root,graph = resume_panel
    plan = repaired_plan(resume, root, graph, monkeypatch)
    before = {p:file_sha256(p) for p in (graph/'tasks').glob('*/*') if p.is_file()}
    monkeypatch.setattr(protocol,'preflight_ready',lambda _:True)
    submitted=[]
    def submit(command, **kwargs):
        kwargs['before_submit'](); submitted.append(Path(kwargs['cwd']).name)
        return SimpleNamespace(returncode=0, stdout='99', stderr='')
    monkeypatch.setattr(run_policy,'submit',submit)
    class Tick(Exception): pass
    monkeypatch.setattr(resume.time,'sleep',lambda _:(_ for _ in ()).throw(Tick()))
    with pytest.raises(Tick): resume.run(root)
    ident = plan['irc'][-1]['task']['id']
    assert submitted == [ident]
    assert dispatch._load(graph)[2]['tasks']['reverse_10']['job_id']=='77'
    assert before == {p:file_sha256(p) for p in before}
    executable = Path(plan['irc'][-1]['task']['command'][0])
    executable.write_text(executable.read_text().replace('assert "StepSize=-10" in text',
        'assert "StepSize=10" in text').replace(f'print({INVALID!r})\nsys.exit(1)',
        'print("Normal termination of Gaussian")'))
    monkeypatch.setenv('SLURM_JOB_ID','99')
    outcome = dispatch.run_task(graph/'tasks'/ident/'task.json')
    assert outcome['status']=='executed'
    monkeypatch.setattr(dispatch,'_job_active',lambda j:str(j) in ('77','99'))
    with pytest.raises(Tick): resume.run(root)
    assert dispatch._load(graph)[2]['tasks'][ident]['status']=='submitted'
    monkeypatch.setattr(dispatch,'_job_active',lambda j:str(j)=='77')
    with pytest.raises(Tick): resume.run(root)
    assert dispatch._load(graph)[2]['tasks'][ident]['status']=='complete'
    assert submitted==[ident] and before=={p:file_sha256(p) for p in before}
    assert resume.reconciliation.prepare_repairs(root)==plan
    assert not json.loads((root/'progress.json').read_text())['scientific_approval']


@pytest.mark.parametrize('defect', ['wrong_job','active','different_failure','already_propagated','wrong_step','changed_parent'])
def test_dvv_repair_rejects_missing_or_stale_evidence(resume_panel, monkeypatch, defect):
    resume,root,graph=resume_panel; _,ident=legacy_failure(resume,root,graph,monkeypatch)
    directory=graph/'tasks'/ident
    if defect=='wrong_job':
        p=directory/'execution.json';r=json.loads(p.read_text());r['slurm_job_id']='89';atomic_json(p,r)
    elif defect=='active':monkeypatch.setattr(dispatch,'_job_active',lambda _:True)
    elif defect=='changed_parent':(graph/'tasks/hessian/hessian.chk').write_text('changed')
    else:
        native = {'different_failure':'SCF failed', 'already_propagated':INVALID+'Point Number: 1',
                  'wrong_step':INVALID.replace('-0.10000','-0.05000')}[defect]
        (directory/(ident+'.log')).write_text(native)
    with pytest.raises((RuntimeError,ValueError)):resume.reconciliation.prepare_repairs(root)
    assert not (root/resume.reconciliation.REPAIRS).exists()


def cc_native():
    # Native VAR(S), VAR(P), iterations 20..39 from the reported tight solve.
    values=[(2.2e-9,2.7e-10),(1.1e-9,1.6e-10),(5.1e-10,7.1e-11),(2.1e-10,3.8e-11),
            (8.8e-11,2.2e-11),(3.7e-11,1.3e-11),(2.2e-11,6.4e-12),(8.2e-12,1.2e-12),
            (1.4e-12,5.4e-13),(6.e-13,3.7e-13),(1.8e-13,1.3e-13),(1.e-13,4.2e-14),
            (8.3e-14,1.1e-14),(7.5e-14,6.2e-15),(6.8e-14,5.1e-15),(5.8e-14,4.7e-15),
            (4.4e-14,3.8e-15),(3.4e-14,2.7e-15),(2.7e-14,1.8e-15),(1.7e-14,1.7e-15)]
    text='ITER. SQ.NORM CORR.ENERGY TOTAL ENERGY ENERGY CHANGE DEN1 VAR(S) VAR(P) DIIS TIME\n'
    for i,(s,p) in enumerate(values,20):text+=f'{i} 1.23696451 -0.47022936 -835.12552056 0.00000000 -0.00000000 {s:.3e} {p:.3e} 6 5 114.85\n'
    return text+'?CONVERGENCE NOT REACHED AFTER MAX. ITERATIONS\n?The problem occurs in ccsd\n'


def scf_native():
    text=''
    for i in range(40,61):text+=f'{i} -834.56992354 -0.00001576 0.39D-04 0.82D-04 9 0 0.07 4.13 diag2\n'
    return text+'?No convergence in rhfpr\n?Error: RHF not converged.\n'


def test_electronic_failure_distinguishes_stagnation_from_contraction():
    assert probe_recovery.electronic_failure(cc_native())['recipe']=='contracting_cc100'
    assert probe_recovery.electronic_failure(scf_native())['recipe']=='rohf_level_shift'
    for text in [cc_native().replace('1.700e-14','9.900e-10'), cc_native()+'out of memory',
                 cc_native()+'PROGRAM * FREQUENCIES', scf_native().replace('0.39D-04','0.39D-10')]:
        with pytest.raises(ValueError):probe_recovery.electronic_failure(text)


@pytest.fixture(params=['cc','scf'])
def audit(tmp_path, monkeypatch, request):
    monkeypatch.setattr(dispatch,'assign_partitions',lambda _:None)
    monkeypatch.setattr(dispatch,'render_site_setup',lambda _:'# mock site')
    monkeypatch.setattr(dispatch,'_job_active',lambda _:False)
    tasks=[molpro_task(ident, 'basis=cc-pVDZ-F12\nrhf\nuccsd(t)-f12b,scale_trip=1\n'+end+'\n',
        geometry_from='initial',scf_convergence=1.e-12,cc_convergence=1.e-11)
        for ident,end in probe_recovery.DERIVATIVES.items()]
    spec=dict(schema=1,name='audit',molecule=dict(symbols=['O','H','H'],positions=[[0,0,0],[1,0,0],[0,1,0]],
              charge=1,multiplicity=2),limits=dict(max_nodes=3),tasks=tasks,
              intent=dict(integration='none',scientific_approval=False))
    site(spec);graph=tmp_path/'audit';dispatch.prepare(spec,graph)
    _,_,state=dispatch._load(graph)
    for i,task in enumerate(tasks,20):
        ident=task['id'];entry=state['tasks'][ident];entry.update(status='failed',job_id=str(i))
        directory=graph/'tasks'/ident
        (directory/(ident+'.out')).write_text(cc_native() if request.param=='cc' else scf_native())
        atomic_json(directory/'execution.json',dict(schema=1,task_id=ident,slurm_job_id=str(i),
            geometry_sha256=entry['geometry_sha256'],status='failed',failure_phase='execution'))
    atomic_json(graph/'state.json',state)
    return graph, request.param


def test_electronic_gate_blocks_derivatives_and_preserves_source(audit, tmp_path, monkeypatch):
    graph,kind=audit
    before={p:file_sha256(p) for p in graph.rglob('*') if p.is_file()}
    repair=probe_recovery.curvature_repair(graph)
    assert before=={p:file_sha256(p) for p in before}
    spec=repair['spec'];assert len(spec['tasks'])==4
    for task in spec['tasks']:
        text=task['input_template']
        assert 'thrden=1.d-11,thrvar=1.d-11' in text and 'energy=1.d-12' in text
        assert 'nocheck' not in text and 'ignore_error' not in text
        assert 'result_parser' not in task and 'geometry_output' not in task
        if kind=='cc':assert ';maxit,100}' in text and '\nrhf\n' in text
        else:assert '{rhf,maxit=60;shift,-1.0,-0.5}' in text and 'maxit,100' not in text
        if task['id']!='electronic':assert task['depends_on']==['electronic']
    module=load('reconcile'); root=tmp_path/'root';root.mkdir()
    module.stage_extra(root,'repaired',spec);module.stage_extra(root,'repaired',spec)
    target=root/'repaired'
    monkeypatch.setattr(protocol,'preflight_ready',lambda _:True)
    submitted=[]
    def submit(command,**kwargs):
        kwargs['before_submit']();submitted.append(Path(kwargs['cwd']).name)
        return SimpleNamespace(returncode=0,stdout=str(100+len(submitted)),stderr='')
    monkeypatch.setattr(run_policy,'submit',submit)
    state=dispatch.advance(target,submit=True);assert submitted==['electronic']
    # A failed gate cannot release any derivative calculation or auto-retry.
    entry=state['tasks']['electronic'];directory=target/'tasks/electronic'
    (directory/'electronic.out').write_text(cc_native())
    atomic_json(directory/'execution.json',dict(schema=1,task_id='electronic',slurm_job_id='101',
        geometry_sha256=entry['geometry_sha256'],status='failed',failure_phase='execution'))
    dispatch.advance(target,submit=True);dispatch.advance(target,submit=True)
    assert submitted==['electronic']
    with pytest.raises(ValueError):probe_recovery.curvature_repair(target)


@pytest.mark.parametrize('drift',['geometry','state'])
def test_curvature_source_cannot_silently_change_geometry_or_state(audit,drift):
    graph,_=audit
    p=graph/'workflow.json';spec=json.loads(p.read_text())
    if drift=='geometry':spec['molecule']['positions'][0][0]+=.1
    else:spec['molecule']['multiplicity']+=2
    atomic_json(p,spec)
    state=json.loads((graph/'state.json').read_text())
    state['spec_sha256']=file_sha256(p);atomic_json(graph/'state.json',state)
    with pytest.raises(ValueError,match='differs from the native attempt'):
        probe_recovery.curvature_repair(graph)


def test_combined_owner_adopts_electronic_gate_before_derivatives(audit, resume_panel, monkeypatch, tmp_path):
    resume,root,graph=resume_panel
    monkeypatch.setattr(dispatch,'_job_active',lambda j:str(j)=='77')
    base,old_ident=legacy_failure(resume,root,graph,monkeypatch)
    source,kind=audit;source.rename(root/'audit_old');source=root/'audit_old'
    spec=dispatch._load(source)[1];atomic_json(source/'requested.json',spec)
    proof=path_review.source(source,'gradient',['gradient.out'])[3]
    base['curvature']['audit_old']=dict(spec=spec,proof=proof)
    base['sha256']=path_review.digest({k:v for k,v in base.items() if k!='sha256'})
    atomic_json(root/resume.reconciliation.DECISION,base)
    plan=resume.reconciliation.prepare_repairs(root)
    retained={p:file_sha256(p) for p in root.glob('*/tasks/*/*') if p.is_file()}
    monkeypatch.setattr(protocol,'preflight_ready',lambda _:True)
    submitted=[]
    def submit(command,**kwargs):
        kwargs['before_submit']();submitted.append(Path(kwargs['cwd']))
        return SimpleNamespace(returncode=0,stdout=str(100+len(submitted)),stderr='')
    monkeypatch.setattr(run_policy,'submit',submit)
    class Tick(Exception):pass
    monkeypatch.setattr(resume.time,'sleep',lambda _:(_ for _ in ()).throw(Tick()))
    with pytest.raises(Tick):resume.run(root)
    assert [p.name for p in submitted]==[plan['irc'][-1]['task']['id'],'electronic']
    target=root/'audit_old_electronic_v2';directory=target/'tasks/electronic'
    # Use the actual worker/scratch/output validators with a stand-in executable.
    executable=tmp_path/'molpro_mock.py'
    executable.write_text(f'#!{sys.executable}\nfrom pathlib import Path\nimport sys\n'
        'p=Path(sys.argv[-1]);s=p.read_text()\n'
        'assert "thrden=1.d-11,thrvar=1.d-11" in s\n'
        'assert "forces," not in s and "frequencies," not in s\n'
        'p.with_suffix(".out").write_text("Molpro calculation terminated\\n")\n'
        'Path("electronic.molden").write_text("[Molden Format]\\n")\n')
    executable.chmod(0o755)
    original=dispatch._external_command
    monkeypatch.setattr(dispatch,'_external_command',lambda task,*a,**kw:
        [str(executable),task['input_name']] if task['backend']=='molpro' else original(task,*a,**kw))
    monkeypatch.setenv('SLURM_JOB_ID','102')
    dispatch.run_task(directory/'task.json')
    monkeypatch.setattr(dispatch,'_job_active',lambda j:str(j) in ('77','101','102'))
    with pytest.raises(Tick):resume.run(root)
    assert len(submitted)==2
    monkeypatch.setattr(dispatch,'_job_active',lambda j:str(j) in ('77','101','103','104','105'))
    with pytest.raises(Tick):resume.run(root)
    assert {p.name for p in submitted[2:]}==set(probe_recovery.DERIVATIVES)
    with pytest.raises(Tick):resume.run(root)
    assert len(submitted)==5 and retained=={p:file_sha256(p) for p in retained}
    assert resume.reconciliation.prepare_repairs(root)==plan
    assert not json.loads((root/'progress.json').read_text())['scientific_approval']


@pytest.mark.parametrize('race', ['none','submitting','staging','archival','unowned_job','scheduler_outage'])
def test_repair_refresh_reuses_complete_race_suite(resume_panel, monkeypatch, race):
    # Exercise the established refresh test with the second immutable decision
    # installed, including retention of native checkpoint imports on sources.
    from tests import test_native_path_review as original
    normal=original.make_plan
    def second_wave(resume,root,patcher):
        monkeypatch.setattr(original,'make_plan',normal)
        # The original test has already recorded the source attempt.
        # legacy_failure needs an executable with a normal terminal stub.
        executable=Path(dispatch._load(resume_panel[2])[1]['tasks'][1]['command'][0])
        program=executable.read_text()
        program=program[:program.index('print(')]+'print("Normal termination of Gaussian")\n'
        executable.write_text(program)
        return repaired_plan(resume,root,resume_panel[2],patcher)
    monkeypatch.setattr(original,'make_plan',second_wave)
    original.test_decision_refresh_checks_before_retiring_controller(resume_panel,monkeypatch,race)
