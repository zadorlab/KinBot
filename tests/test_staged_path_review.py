"""Reviewed half-step release through the real dispatch/owner and shared pool."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import run_capacity, run_policy
from kinbot.anl import dispatch, live_control, path_review, protocol, staged_path_review
from kinbot.execution import atomic_json, file_sha256
from tests.test_gaussian_checkpoint_inputs import graph, resume_panel
from tests.test_native_probe_recovery import audit
from tests.test_native_path_review import retained_attempt, one_point, corrector
from tests.test_dvv_departure_probe import departure_plan
from tests.test_r8_validation_panel import spectrum, load
from tests.test_anl_dispatch import _mock_execution


def add_sampling(resume, root, graph, monkeypatch):
    _, spec, state = dispatch._load(graph)
    fine = staged_path_review._rename(deepcopy(spec['tasks'][1]), 'forward_10', 'forward_5')
    fine['input_template'] = fine['input_template'].replace('StepSize=-10', 'StepSize=-5').replace('MaxPoints=600', 'MaxPoints=1200')
    spec['tasks'].append(fine)
    live_control.commit(graph, spec, state, set(), {})
    monkeypatch.setattr(resume.campaign, 'graphs', lambda _: {'example':deepcopy(spec)})
    atomic_json(graph/'requested.json', spec)
    dispatch.advance(graph)
    _, _, state = dispatch._load(graph)
    state['tasks']['hessian']['job_id'] = '55'
    atomic_json(graph/'state.json', state)
    native = graph/'tasks/hessian/hessian.log'
    native.write_text(spectrum([-350., 800., 1000.]))
    path = graph/'tasks/hessian/execution.json'; record = json.loads(path.read_text())
    record['slurm_job_id'] = '55'; record['artifacts']['hessian.log'] = file_sha256(native)
    atomic_json(path, record)
    return fine


@pytest.fixture
def sampling(resume_panel, monkeypatch):
    resume, root, graph = resume_panel
    fine = add_sampling(resume, root, graph, monkeypatch)
    retained_attempt(graph, monkeypatch, 'failed')
    return resume, root, graph, fine


def test_release_is_one_new_sampling_with_proofs_and_no_mutation(sampling):
    _, root, graph, fine = sampling
    before = {p:file_sha256(p) for p in graph.glob('tasks/*/*') if p.is_file()}
    result = staged_path_review.prepare(graph, 'forward_10', 'forward_5')
    assert result['reason'] == 'corrector_limit_untried_half_step'
    assert result['definition'] == fine and result['budget'] == 1
    assert not result['scientific_approval']
    assert before == {p:file_sha256(p) for p in before}
    assert staged_path_review.prepare(graph, 'forward_10', 'forward_5') == result


@pytest.mark.parametrize('defect', ['attempted', 'orphan_output', 'wrong_job', 'active', 'checkpoint',
    'input', 'model', 'step', 'budget', 'limit', 'unknown_failure', 'parent_modes', 'parent_job', 'alias'])
def test_release_rejects_changed_or_exhausted_work(sampling, monkeypatch, defect):
    _, root, graph, _ = sampling
    _, spec, state = dispatch._load(graph)
    path = graph/'tasks/forward_10/forward_10.log'
    if defect == 'attempted':
        state['tasks']['forward_5']['attempt'] = 2; atomic_json(graph/'state.json', state)
    elif defect == 'orphan_output': (graph/'tasks/forward_5/forward_5.log').write_text('partial native output')
    elif defect in ('wrong_job', 'parent_job'):
        name = 'hessian' if defect == 'parent_job' else 'forward_10'
        p = graph/'tasks'/name/'execution.json'; record = json.loads(p.read_text())
        record['slurm_job_id'] = '999'; atomic_json(p, record)
    elif defect == 'active': monkeypatch.setattr(dispatch, '_job_active', lambda j:str(j) in ('42','77'))
    elif defect == 'checkpoint': (graph/'tasks/hessian/hessian.chk').write_text('changed')
    elif defect == 'input': (graph/'tasks/forward_5/forward_5.com').write_text('changed')
    elif defect in ('model', 'step', 'budget', 'limit'):
        fine = spec['tasks'][-1]
        old, new = {'model':('cc-pVTZ','cc-pVDZ'), 'step':('StepSize=-5','StepSize=-2'),
                    'budget':('MaxPoints=1200','MaxPoints=2400'),
                    'limit':('MaxCycle=100','MaxCycle=150')}[defect]
        fine['input_template'] = fine['input_template'].replace(old, new)
        atomic_json(graph/'workflow.json', spec)
    elif defect == 'unknown_failure':
        path.write_text('SCF did not converge\nError termination of Gaussian\n')
    elif defect == 'parent_modes':
        native = graph/'tasks/hessian/hessian.log'; native.write_text(spectrum([-350., -50., 1000.]))
        p = native.parent/'execution.json'; record = json.loads(p.read_text())
        record['artifacts']['hessian.log'] = file_sha256(native); atomic_json(p, record)
    elif defect == 'alias':
        duplicate=staged_path_review._rename(deepcopy(spec['tasks'][-1]),'forward_5','earlier_half_step')
        spec['tasks'].append(duplicate);atomic_json(graph/'workflow.json',spec)
        state['tasks']['earlier_half_step']=dict(status='failed',job_id='999');atomic_json(graph/'state.json',state)
    with pytest.raises((ValueError, RuntimeError)):
        staged_path_review.prepare(graph, 'forward_10', 'forward_5')


@pytest.mark.parametrize('points,accepted', [(1,False), (22,True)])
def test_completed_path_only_supports_sampling_after_actual_departure(resume_panel, monkeypatch, points, accepted):
    resume, root, graph = resume_panel
    add_sampling(resume, root, graph, monkeypatch)
    executable = Path(dispatch._load(graph)[1]['tasks'][1]['command'][0])
    output = one_point().replace('PATH = 1', 'PATH = '+str(points))
    executable.write_text(executable.read_text().replace('print("Normal termination of Gaussian")', f'print({output!r})'))
    _, _, state = dispatch._load(graph); state['tasks']['forward_10'].update(status='submitted',job_id='42')
    atomic_json(graph/'state.json', state); monkeypatch.setenv('SLURM_JOB_ID','42')
    dispatch.run_task(graph/'tasks/forward_10/task.json'); dispatch.advance(graph)
    if accepted:
        assert staged_path_review.prepare(graph, 'forward_10', 'forward_5')['reason'] == 'raw_path_endpoint_sensitivity'
    else:
        with pytest.raises(staged_path_review.Inapplicable, match='early_gradient_stop'):
            staged_path_review.prepare(graph, 'forward_10', 'forward_5')


def sampling_plan(resume, root, graph, audit, monkeypatch):
    add_sampling(resume, root, graph, monkeypatch)
    departure_plan(resume, root, graph, audit, monkeypatch)
    before = {p:p.read_bytes() for p in root.glob('reconciliation*.json')}
    plan = resume.reconciliation.prepare_sampling(root)
    assert len(plan['released']) == 1
    assert before == {p:p.read_bytes() for p in before}
    assert resume.reconciliation.prepare_sampling(root) == plan
    return plan


@pytest.mark.parametrize('audit', ['scf'], indirect=True)
@pytest.mark.parametrize('race', ['none','submitting','staging','archival','unowned_job','scheduler_outage'])
def test_release_refresh_checks_full_snapshot_before_owner_retirement(audit, resume_panel, monkeypatch, race):
    from tests import test_native_path_review as original
    def fifth_wave(resume, root, patcher):
        executable = Path(dispatch._load(resume_panel[2])[1]['tasks'][1]['command'][0])
        program = executable.read_text()
        executable.write_text(program[:program.index('print(')]+'print("Normal termination of Gaussian")\n')
        return sampling_plan(resume, root, resume_panel[2], audit, patcher)
    monkeypatch.setattr(original, 'make_plan', fifth_wave)
    original.test_decision_refresh_checks_before_retiring_controller(resume_panel, monkeypatch, race)


@pytest.mark.parametrize('audit', ['scf'], indirect=True)
@pytest.mark.parametrize('native_success', [True, False])
def test_owner_release_capacity_backfill_adoption_and_repeat(audit, resume_panel, monkeypatch, native_success):
    resume, root, graph = resume_panel
    plan = sampling_plan(resume, root, graph, audit, monkeypatch)
    before = {p:file_sha256(p) for p in graph.glob('tasks/*/*') if p.is_file()}
    run_policy.update(root, max_jobs=2, exclude=['bad'])
    jobs = {'77':1, 'unrelated':1}; calls = []
    monkeypatch.setattr(run_policy, '_queue', lambda _:dict(jobs))
    monkeypatch.setattr(run_policy, 'slurm_job_record', lambda _: {})
    monkeypatch.setattr(dispatch, '_job_active', lambda job: str(job) in jobs)
    monkeypatch.setattr(protocol, 'preflight_ready', lambda _:True)
    native_run = run_policy.subprocess.run
    def run(command, **kwargs):
        if command[0] == 'sinfo':
            return SimpleNamespace(stdout='n1|allocated|96|126000|up\nn2|allocated|96|126000|up\nn3|idle|96|126000|up\n')
        if command[0] == 'sbatch':
            calls.append((command, Path(kwargs['cwd']))); ident = str(400+len(calls))
            jobs[ident] = 1
            return SimpleNamespace(returncode=0, stdout=ident+'\n', stderr='')
        return native_run(command, **kwargs)
    monkeypatch.setattr(run_policy.subprocess, 'run', run)
    class Tick(Exception): pass
    monkeypatch.setattr(resume.time, 'sleep', lambda _:(_ for _ in ()).throw(Tick()))
    with pytest.raises(Tick): resume.run(root)
    assert not calls  # Both slots are occupied by retained allocations.
    assert resume.reconciliation.sampling_capacity(root)['max_jobs'] == 3
    with pytest.raises(Tick): resume.run(root)
    assert len(calls) == 1 and calls[0][1].name == 'forward_5'
    assert '--exclude=bad' in calls[0][0] and '--exclusive' in calls[0][0]
    assert before == {p:file_sha256(p) for p in before}
    executable = Path(plan['released'][0]['definition']['command'][0])
    program = executable.read_text()
    program = program[:program.index('print(')].replace('GradientOnly,DVV', 'GradientOnly,EulerPC').replace('StepSize=10', 'StepSize=-5')
    output = 'Normal termination of Gaussian' if native_success else 'Error termination of Gaussian'
    executable.write_text(program+f'print({output!r})\nsys.exit({int(not native_success)})\n')
    monkeypatch.setenv('SLURM_JOB_ID','401')
    if native_success: dispatch.run_task(graph/'tasks/forward_5/task.json')
    else:
        with pytest.raises(RuntimeError): dispatch.run_task(graph/'tasks/forward_5/task.json')
    with pytest.raises(Tick): resume.run(root)
    assert dispatch._load(graph)[2]['tasks']['forward_5']['status'] == 'submitted'
    jobs.pop('401')
    # Clear the ledger only with a terminal scheduler record, as on the site.
    monkeypatch.setattr(run_policy, 'slurm_job_record', lambda j:{'JobState':'COMPLETED'} if j=='401' else {})
    with pytest.raises(Tick): resume.run(root)
    assert dispatch._load(graph)[2]['tasks']['forward_5']['status'] == ('complete' if native_success else 'failed')
    assert len(calls) == 2 and calls[-1][1].name == 'electronic'
    with pytest.raises(Tick): resume.run(root)
    assert len(calls) == 2 and set(jobs) == {'77', 'unrelated', '402'}
    jobs.pop('unrelated')
    with pytest.raises(Tick): resume.run(root)
    assert len(calls) == 3 and calls[-1][1].name.endswith('_departure_v3')
    assert set(jobs) == {'77', '402', '403'}
    assert dispatch._load(graph)[2]['tasks']['reverse_10']['job_id'] == '77'
    assert resume.reconciliation.prepare_sampling(root) == plan
    assert not resume.reconciliation.sampling_capacity(root)['applied']
    assert not json.loads((root/'progress.json').read_text())['scientific_approval']
    assert before == {p:file_sha256(p) for p in before}


@pytest.mark.parametrize('audit', ['scf'], indirect=True)
def test_release_target_drift_blocks_refresh_and_repeat(audit, resume_panel, monkeypatch):
    resume, root, graph = resume_panel
    sampling_plan(resume, root, graph, audit, monkeypatch)
    (graph/'tasks/forward_5/forward_5.com').write_text('changed')
    with pytest.raises(ValueError, match='sampling files changed'): resume.reconciliation.prepare_sampling(root)
    with pytest.raises(ValueError, match='sampling files changed'): resume.snapshot(root)


@pytest.mark.parametrize('audit', ['scf'], indirect=True)
def test_watcher_reports_release_and_capacity_without_changes(audit, resume_panel, monkeypatch, capsys):
    resume, root, graph = resume_panel
    sampling_plan(resume, root, graph, audit, monkeypatch)
    run_policy.update(root, max_jobs=21, exclude=['bad'])
    before={p:file_sha256(p) for p in root.rglob('*') if p.is_file()}
    watch=load('watch')
    monkeypatch.setattr(watch.subprocess,'run',lambda *a,**k:SimpleNamespace(stdout='77|RUNNING|2:00|'+str(graph/'tasks/reverse_10')+'\n'))
    watch.snapshot(root)
    output=capsys.readouterr().out
    assert 'example forward_5 state=staged' in output
    assert 'POOL max_jobs 21 excluded_nodes 1' in output
    assert before=={p:file_sha256(p) for p in before}


def test_full_sampling_wave_fills_21_slots_retains_seven_jobs_and_backfills(tmp_path, monkeypatch):
    resume = load('resume'); root = tmp_path/'native'; root.mkdir()
    monkeypatch.setattr(dispatch, 'assign_partitions', lambda _:None)
    monkeypatch.setattr(dispatch, 'render_site_setup', lambda _: '# simulated site\n')
    jobs = {str(i):1 for i in range(7)}; calls = []
    monkeypatch.setattr(dispatch, '_job_active', lambda j:str(j) in jobs)
    point = dict(symbols=['O','H','H'], positions=[[0.,0,0],[1.,0,0],[0,1.,0]], charge=0, multiplicity=1)
    specs = {}
    for index in range(12):
        name = 'irc_case_'+str(index); graph = root/name
        spec = resume.campaign.native_irc(name, point, graph)
        for task in spec['tasks']:
            task['resources'] = dict(cores=8, memory_mb=98000, walltime='1-00:00:00', partition='month-long-cpu')
        specs[name] = spec
        dispatch.prepare(spec, graph); atomic_json(graph/'requested.json', spec)
        _mock_execution(graph, 'hessian')
        native = graph/'tasks/hessian/hessian.log'; native.write_text(spectrum([-350.,800.,1000.]))
        p = native.parent/'execution.json'; record=json.loads(p.read_text())
        record['artifacts']['hessian.log']=file_sha256(native); record['slurm_job_id']=str(100+index)
        atomic_json(p, record)
        _, _, state = dispatch._load(graph)
        state['tasks']['hessian']['job_id']=str(100+index); atomic_json(graph/'state.json', state)
        dispatch.advance(graph)
        for offset, ident in enumerate(('forward_10', 'reverse_10')):
            job = str(200+2*index+offset)
            _mock_execution(graph, ident)
            p=graph/'tasks'/ident/'execution.json'; record=json.loads(p.read_text())
            native=p.parent/(ident+'.log')
            native.write_text(one_point() if index==0 else corrector())
            record.update(status='executed' if index==0 else 'failed', slurm_job_id=job)
            record['artifacts'][ident+'.log']=file_sha256(native); atomic_json(p, record)
            _, _, state=dispatch._load(graph)
            state['tasks'][ident].update(status='complete' if index==0 else 'failed',job_id=job)
            atomic_json(graph/'state.json',state)
    monkeypatch.setattr(resume.campaign, 'graphs', lambda _:deepcopy(specs))
    atomic_json(root/'stage.json', dict(sources=resume.campaign.fingerprint(), workflows=list(specs)))
    module=resume.reconciliation
    base=dict(schema=1,root=str(root),frozen_sources=resume.campaign.fingerprint(),held_irc_graphs=list(specs),irc=[],curvature={})
    base['sha256']=path_review.digest(base);atomic_json(root/module.DECISION,base)
    for filename in (module.REPAIRS,module.ELECTRONIC,module.DEPARTURE):
        extension=dict(schema=1,root=str(root),parent_sha256=module.load(root)['sha256'],irc=[],curvature={})
        extension['sha256']=path_review.digest(extension);atomic_json(root/filename,extension)
    plan=module.prepare_sampling(root)
    assert len(plan['released'])==22
    assert len([r for r in plan['sampling_review'].values() if r['status']=='held'])==2
    before={p:file_sha256(p) for p in root.glob('*/tasks/*/*') if p.is_file()}
    run_policy.update(root,max_jobs=20,exclude=['excluded'])
    monkeypatch.setattr(run_policy,'_queue',lambda _:dict(jobs))
    monkeypatch.setattr(run_policy,'slurm_job_record',lambda _:{})
    monkeypatch.setattr(protocol,'preflight_ready',lambda _:True)
    def run(command,**kwargs):
        if command[0]=='sinfo':
            return SimpleNamespace(stdout=''.join(f'n{i}|idle|96|126000|up\n' for i in range(21)))
        assert command[0]=='sbatch'
        calls.append((command,Path(kwargs['cwd'])));job=str(500+len(calls));jobs[job]=1
        return SimpleNamespace(returncode=0,stdout=job+'\n',stderr='')
    monkeypatch.setattr(run_policy.subprocess,'run',run)
    assert module.sampling_capacity(root)['max_jobs']==21
    class Tick(Exception):pass
    monkeypatch.setattr(resume.time,'sleep',lambda _:(_ for _ in ()).throw(Tick()))
    with pytest.raises(Tick):resume.run(root)
    assert len(calls)==14 and len(jobs)==21
    assert all('--exclude=excluded' in command for command,_ in calls)
    with pytest.raises(Tick):resume.run(root)
    assert len(calls)==14
    for _, directory in calls[:]:
        job=dispatch._load(directory.parents[1])[2]['tasks'][directory.name]['job_id']
        _mock_execution(directory.parents[1],directory.name)
        p=directory/'execution.json';record=json.loads(p.read_text());record['slurm_job_id']=job;atomic_json(p,record)
        jobs.pop(job)
    with pytest.raises(Tick):resume.run(root)
    # Fair admission may defer a prior graph until its next controller poll.
    with pytest.raises(Tick):resume.run(root)
    assert len(calls)==22 and len(jobs)==15
    assert all(str(i) in jobs for i in range(7))
    assert len({d for _,d in calls})==22
    progress=json.loads((root/'progress.json').read_text())
    assert sum(len(v) for v in progress['held_for_review'].values())==2
    assert not progress['scientific_approval']
    assert module.prepare_sampling(root)==plan
    selected=load('collect').selected(tmp_path)
    assert root/module.RELEASE in selected and root/'sampling_capacity.json' in selected
    assert before=={p:file_sha256(p) for p in before}
