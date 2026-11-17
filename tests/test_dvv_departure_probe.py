"""Evidence-bound DVV departure trial and combined electronic/path handoff."""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import run_policy
from kinbot.anl import dispatch, live_control, probe_recovery, protocol
from kinbot.execution import atomic_json, file_sha256
from tests.test_gaussian_checkpoint_inputs import graph, resume_panel
from tests.test_native_probe_recovery import audit, legacy_failure, INVALID
from tests.test_quadratic_probe_recovery import quadratic_plan
from tests.test_anl_dispatch import _mock_execution

pytestmark = pytest.mark.parametrize('audit', ['scf'], indirect=True)

EARLY = '''Integration scheme = DVV
Step size = 0.100 fs
Current Structure is TS -> form Hessian eigenvectors.
NET REACTION COORDINATE UP TO THIS POINT = 0.00400
# OF POINTS ALONG THE PATH = 1
PES minimum detected on this side of the pathway.
Magnitude of the gradient = 0.0000093
Normal termination of Gaussian
'''


def departure_plan(resume, root, graph, audit, monkeypatch):
    plan, _ = quadratic_plan(resume, root, graph, audit, monkeypatch,
        base_factory=lambda *a: legacy_failure(*a)[0])
    task = plan['irc'][-1]['task']; ident = task['id']
    executable = Path(task['command'][0])
    executable.write_text(executable.read_text().replace('assert "StepSize=-10" in text',
        'assert "StepSize=10" in text').replace(f'print({INVALID!r})\nsys.exit(1)', f'print({EARLY!r})'))
    dispatch.advance(graph)
    _, _, state = dispatch._load(graph)
    state['tasks'][ident].update(status='submitted', job_id='211')
    atomic_json(graph/'state.json', state); monkeypatch.setenv('SLURM_JOB_ID', '211')
    dispatch.run_task(graph/'tasks'/ident/'task.json'); dispatch.advance(graph)
    assert dispatch._load(graph)[2]['tasks'][ident]['status'] == 'complete'
    monkeypatch.setattr(resume.reconciliation, 'DEPARTURE_PILOT', (graph.name, ident))
    prior = {p:p.read_bytes() for p in root.glob('reconciliation*.json')}
    merged = resume.reconciliation.prepare_departure(root)
    assert prior == {p:p.read_bytes() for p in prior}
    assert resume.reconciliation.prepare_departure(root) == merged
    return merged, ident


@pytest.mark.parametrize('native_success', [True, False])
def test_combined_departure_and_so_wave_adopts_outputs_without_releasing_grid(audit, resume_panel, monkeypatch, native_success):
    resume, root, graph = resume_panel
    _, spec, state = dispatch._load(graph)
    extra = deepcopy(spec['tasks'][1]); extra['id'] = 'forward_5'
    spec['tasks'].append(extra)
    live_control.commit(graph, spec, state, set(), {})
    monkeypatch.setattr(resume.campaign, 'graphs', lambda _: {'example':deepcopy(spec)})
    atomic_json(graph/'requested.json', spec); dispatch.advance(graph)
    plan, source = departure_plan(resume, root, graph, audit, monkeypatch)
    task = plan['irc'][-1]['task']; ident = task['id']
    text = task['input_template']
    assert 'StepSize=10,MaxPoints=64,MaxCycle=100' in text
    assert 'IOp(1/108=-1)' in text and '--Link1--' not in text
    assert '1/109' not in text and '1/106' not in text
    assert 'result_parser' not in task and 'geometry_output' not in task
    before = {p:file_sha256(p) for p in root.glob('*/tasks/*/*') if p.is_file()}
    monkeypatch.setattr(protocol, 'preflight_ready', lambda _:True)
    submitted = []
    def submit(command, **kwargs):
        kwargs['before_submit'](); submitted.append(Path(kwargs['cwd']))
        return SimpleNamespace(returncode=0, stdout=str(400+len(submitted)), stderr='')
    monkeypatch.setattr(run_policy, 'submit', submit)
    class Tick(Exception): pass
    monkeypatch.setattr(resume.time, 'sleep', lambda _:(_ for _ in ()).throw(Tick()))
    with pytest.raises(Tick): resume.run(root)
    assert [p.name for p in submitted] == [ident, 'electronic']
    assert before == {p:file_sha256(p) for p in before}
    executable = Path(task['command'][0])
    native = ('# OF POINTS ALONG THE PATH = 64\nMaximum number of points reached.\n'
              'NET REACTION COORDINATE UP TO THIS POINT = 0.256\nNormal termination of Gaussian\n')
    if not native_success: native = 'Error termination via Lnk1e in l123.exe\n'
    executable.write_text(executable.read_text().replace(f'print({EARLY!r})',
        'assert "MaxPoints=64" in text and "IOp(1/108=-1)" in text\n'+f'print({native!r})'))
    monkeypatch.setenv('SLURM_JOB_ID', '401')
    if native_success: dispatch.run_task(graph/'tasks'/ident/'task.json')
    else:
        with pytest.raises(RuntimeError): dispatch.run_task(graph/'tasks'/ident/'task.json')
    monkeypatch.setattr(dispatch, '_job_active', lambda j:str(j) in ('77','401','402'))
    with pytest.raises(Tick): resume.run(root)
    assert dispatch._load(graph)[2]['tasks'][ident]['status'] == 'submitted'
    monkeypatch.setattr(dispatch, '_job_active', lambda j:str(j) in ('77','402'))
    with pytest.raises(Tick): resume.run(root)
    with pytest.raises(Tick): resume.run(root)
    state = dispatch._load(graph)[2]
    assert state['tasks'][ident]['status'] == ('complete' if native_success else 'failed')
    assert state['tasks']['reverse_10']['job_id'] == '77'
    assert state['tasks']['forward_5']['status'] == 'staged'
    # An existing sampling allocation is adopted even while unsent siblings
    # are held; it is neither relaunched nor mistaken for a new probe.
    state['tasks']['forward_5'].update(status='submitted', job_id='499')
    atomic_json(graph/'state.json', state); _mock_execution(graph, 'forward_5')
    receipt=graph/'tasks/forward_5/execution.json'
    entry=json.loads(receipt.read_text()); entry['slurm_job_id']='499'; atomic_json(receipt, entry)
    monkeypatch.setattr(dispatch, '_job_active', lambda j:str(j) in ('77','402','499'))
    with pytest.raises(Tick): resume.run(root)
    assert dispatch._load(graph)[2]['tasks']['forward_5']['status'] == 'submitted'
    monkeypatch.setattr(dispatch, '_job_active', lambda j:str(j) in ('77','402'))
    with pytest.raises(Tick): resume.run(root)
    assert dispatch._load(graph)[2]['tasks']['forward_5']['status'] == 'complete'
    assert len(submitted) == 2 and before == {p:file_sha256(p) for p in before}
    assert resume.reconciliation.prepare_departure(root) == plan
    assert resume.reconciliation.prepare_quadratic(root) == plan
    assert not json.loads((root/'progress.json').read_text())['scientific_approval']
    with pytest.raises(ValueError): probe_recovery.dvv_departure_probe(graph, ident)


@pytest.mark.parametrize('defect', ['active','wrong_job','long_path','late_stop','zero_gradient',
    'wrong_time','not_ts','changed_parent','changed_ancestor','input_drift'])
def test_departure_refuses_unproven_or_changed_sources(audit, resume_panel, monkeypatch, defect):
    resume, root, graph = resume_panel
    _, ident = departure_plan(resume, root, graph, audit, monkeypatch)
    directory = graph/'tasks'/ident
    if defect == 'active': monkeypatch.setattr(dispatch, '_job_active', lambda j:str(j) in ('77','211'))
    elif defect == 'wrong_job':
        p=directory/'execution.json';e=json.loads(p.read_text());e['slurm_job_id']='212';atomic_json(p,e)
    elif defect == 'changed_parent': (graph/'tasks/hessian/hessian.chk').write_text('changed')
    elif defect == 'changed_ancestor':
        (graph/'tasks/forward_10_dvv_v1/forward_10_dvv_v1.log').write_text('another failure')
    elif defect == 'input_drift':
        p=directory/(ident+'.com');p.write_text(p.read_text().replace('StepSize=10','StepSize=20'))
    else:
        native = EARLY
        replacements = {'long_path':('0.00400','0.40000'), 'late_stop':('PATH = 1','PATH = 2'),
                        'zero_gradient':('0.0000093','0.0000000'), 'wrong_time':('0.100 fs','0.200 fs'),
                        'not_ts':('Current Structure is TS','Current Structure is downhill')}
        native=native.replace(*replacements[defect]);(directory/(ident+'.log')).write_text(native)
    with pytest.raises((ValueError, RuntimeError)): probe_recovery.dvv_departure_probe(graph, ident)


@pytest.mark.parametrize('race', ['none','submitting','staging','archival','unowned_job','scheduler_outage'])
def test_fourth_wave_refresh_preserves_live_work_and_transient_races(audit, resume_panel, monkeypatch, race):
    from tests import test_native_path_review as original
    def fourth_wave(resume, root, patcher):
        executable=Path(dispatch._load(resume_panel[2])[1]['tasks'][1]['command'][0])
        program=executable.read_text()
        executable.write_text(program[:program.index('print(')]+'print("Normal termination of Gaussian")\n')
        return departure_plan(resume, root, resume_panel[2], audit, patcher)[0]
    monkeypatch.setattr(original, 'make_plan', fourth_wave)
    original.test_decision_refresh_checks_before_retiring_controller(resume_panel, monkeypatch, race)
