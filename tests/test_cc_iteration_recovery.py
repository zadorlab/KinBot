from copy import deepcopy
import json
from pathlib import Path

import pytest
from ase import Atoms
from kinbot.anl import dispatch, cc_iteration_recovery as recovery
from kinbot.ase_modules.calculators.molpro import render_input
from kinbot.execution import atomic_json
from tests.anl_fixture import dispatch_spec, molpro_task

NATIVE = ('?CONVERGENCE NOT REACHED AFTER MAX. ITERATIONS\n'
          '?ERROR: NO CONVERGENCE. TRIPLES WILL NOT BE DONE.\nThe problem occurs in ccsd\n')


def fail(graph):
    _, spec, state = dispatch._load(graph)
    task = spec['tasks'][0]; ident = task['id']; d = graph/'tasks'/ident
    entry = state['tasks'][ident]; entry.update(status='failed', job_id='42')
    atomic_json(d/'execution.json', dict(schema=1, task_id=ident, status='failed',
        geometry_sha256=entry['geometry_sha256'], slurm_job_id='42', failure_phase='execution',
        error='Molpro force evaluation l3_geometry_step_0001 exited with status 1; '
              'diagnostics: l3_geometry_step_0001.failure.json.'))
    atomic_json(d/'l3_geometry_step_0001.failure.json', dict(returncode=1))
    (d/'l3_geometry_step_0001.out').write_text(NATIVE)
    atomic_json(graph/'state.json', state)


def test_dispatch_bounded_cc_ladder_and_retained_contract(tmp_path, monkeypatch):
    spec = dispatch_spec(); task = deepcopy(spec['tasks'][1])
    task['geometry_from'] = 'initial'
    task['profile']['calculator_kwargs'] = dict(scf_convergence=1.e-11, correlation_convergence=1.e-11)
    spec['tasks'] = [task]
    graph = dispatch.prepare(spec, tmp_path/'graph')
    monkeypatch.setattr(dispatch, '_job_active', lambda job: False)
    for limit in (100, 200):
        fail(graph)
        state = dispatch.advance(graph)
        entry = state['tasks'][task['id']]
        assert entry['status'] == 'staged'
        assert (Path(entry['previous_attempt'])/'l3_geometry_step_0001.out').read_text() == NATIVE
        actual = dispatch._load(graph)[1]['tasks'][0]
        from kinbot.manual_fallback import logical_task
        assert logical_task(actual) == task
        kwargs = actual['profile']['calculator_kwargs']
        assert kwargs == dict(task['profile']['calculator_kwargs'], cc_max_iterations=limit)
        assert actual['optimizer'] == task['optimizer']
        text = render_input(Atoms('CH4', positions=spec['molecule']['positions']),
            basis='cc-pVTZ', geometry_name='geom.xyz', method='CCSD(T)', mult=1, **kwargs)
        assert f'maxit={limit}' in text and 'nocheck' not in text.lower()
    fail(graph)
    state = dispatch.advance(graph)
    assert state['tasks'][task['id']]['status'] == 'failed'
    assert [h['limit'] for h in state['cc_iteration_retries'][task['id']]] == [100, 200]


@pytest.mark.parametrize('method', ['CCSD(T)', 'CCSD(T)-F12b'])
def test_single_point_uses_native_output_not_launcher(tmp_path, method):
    task = molpro_task('sp', '{uccsd(t)-f12b,scale_trip=1,thrden=1.d-11,thrvar=1.d-11}\n',
        geometry_from='initial', result_parser=dict(kind='molpro_energy', method=method))
    (tmp_path/'launcher.stdout').write_text('MPI abort')
    (tmp_path/'sp.out').write_text(NATIVE)
    outcome = dict(status='failed', failure_phase='execution')
    change, proof = recovery.options(task, tmp_path, outcome)
    assert proof['limit'] == 100 and 'sp.out' in proof['native']
    new = change['input_edits'][0]['new']
    assert 'thrden=1.d-11,thrvar=1.d-11' in new and ';maxit,100}' in new
    (tmp_path/'sp.out').write_text('SCF did not converge')
    assert recovery.options(task, tmp_path, outcome) is None


def test_cfour_cc_only_and_mrcc_compiled_limit(tmp_path):
    task = dict(kind='external', backend='cfour', input_template='CC_PROG=VCC\nCC_MAXCYC=300\n',
        result_parser=dict(method='CCSDT(Q)',file='cfour.out'))
    (tmp_path/'cfour.out').write_text('CCSD equations failed to converge\n')
    outcome = dict(status='failed', failure_phase='execution')
    assert recovery.options(task,tmp_path,outcome)[1]['limit'] == 600
    task['input_template'] = task['input_template'].replace('300', '600')
    assert recovery.options(task,tmp_path,outcome)[1]['limit'] == 1000
    task['backend'] = 'mrcc'
    assert recovery.options(task,tmp_path,outcome) is None


@pytest.mark.parametrize('value', [True, 0, -1, 1001, 10.5])
def test_invalid_iteration_limits(value):
    with pytest.raises(ValueError):
        render_input(Atoms('H2'), basis='cc-pVDZ', geometry_name='x.xyz',
                     method='CCSD(T)', mult=1, cc_max_iterations=value)


INPUT_ERROR = ('Input error: Attempt to set Variable MAXIT in command block\n'
               '? Error\n? Input error\n? The problem occurs in check_input\n')


def fail_external(graph, text=NATIVE):
    _, spec, state = dispatch._load(graph)
    task = spec['tasks'][0]; ident = task['id']; directory = graph/'tasks'/ident
    state['tasks'][ident].update(status='failed', job_id='42')
    (directory/(ident+'.out')).write_text(text)
    atomic_json(directory/'execution.json', dict(schema=1, task_id=ident,
        status='failed', failure_phase='execution', slurm_job_id='42',
        geometry_sha256=state['tasks'][ident]['geometry_sha256'], error='Native failure'))
    atomic_json(graph/'state.json', state)


def external_spec(body):
    spec = dispatch_spec()
    spec['tasks'] = [molpro_task('sp', 'basis=cc-pVTZ-F12\nrhf\n'+body,
        geometry_from='initial', result_parser=dict(kind='molpro_energy',
            method='CCSD(T)-F12b', basis='cc-pVTZ-F12', reference='RHF'))]
    return spec


def test_external_dispatch_preserves_directive_on_second_retry(tmp_path, monkeypatch):
    task_body = '{uccsd(t)-f12b,scale_trip=1,thrden=1.d-11,thrvar=1.d-11}\n'
    spec = external_spec(task_body)
    graph = dispatch.prepare(spec, tmp_path/'graph')
    monkeypatch.setattr(dispatch, '_job_active', lambda job: False)
    for limit in (100, 200):
        fail_external(graph)
        state = dispatch.advance(graph)
        assert state['tasks']['sp']['status'] == 'staged'
        template = dispatch._load(graph)[1]['tasks'][0]['input_template']
        assert f';maxit,{limit}}}' in template
        assert ';maxit=' not in template
        assert 'thrden=1.d-11,thrvar=1.d-11' in template
        assert (Path(state['tasks']['sp']['previous_attempt'])/'sp.out').read_text() == NATIVE
        assert dispatch.advance(graph)['tasks']['sp']['status'] == 'staged'
    fail_external(graph)
    assert dispatch.advance(graph)['tasks']['sp']['status'] == 'failed'


def test_repair_rejected_second_retry_keeps_budget_and_live_work(tmp_path, monkeypatch):
    spec = external_spec('{uccsd(t)-f12b,scale_trip=1;maxit=200}\n')
    other = deepcopy(spec['tasks'][0]); other['id'] = 'live'
    other['input_name'] = 'live.inp'; other['required_outputs'] = ['live.out']
    other['result_parser']['file'] = 'live.out'; other['success_marker']['file'] = 'live.out'
    spec['tasks'].append(other)
    graph = dispatch.prepare(spec, tmp_path/'graph')
    fail_external(graph, INPUT_ERROR)
    _, _, state = dispatch._load(graph)
    history = [dict(limit=100), dict(limit=200)]
    state['cc_iteration_retries'] = {'sp':history}
    state['tasks']['live'].update(status='submitted', job_id='43')
    atomic_json(graph/'state.json', state)
    live = deepcopy(state['tasks']['live'])
    monkeypatch.setattr(dispatch, '_job_active', lambda job: job == '43')
    state = dispatch.advance(graph)
    assert state['tasks']['sp']['status'] == 'staged'
    assert state['tasks']['live'] == live
    assert state['cc_iteration_retries']['sp'] == history
    assert state['cc_iteration_syntax_repairs']['sp']['limit'] == 200
    assert (Path(state['tasks']['sp']['previous_attempt'])/'sp.out').read_text() == INPUT_ERROR
    assert ';maxit,200}' in dispatch._load(graph)[1]['tasks'][0]['input_template']
    assert dispatch.advance(graph)['tasks']['sp']['status'] == 'staged'
    # A subsequent genuine CC failure remains exhausted at the same limit.
    fail_external(graph)
    assert dispatch.advance(graph)['tasks']['sp']['status'] == 'failed'


@pytest.mark.parametrize('history,text', [([],INPUT_ERROR), ([dict(limit=100)],INPUT_ERROR),
                                         ([dict(limit=200)],NATIVE)])
def test_syntax_repair_requires_matching_retry_and_native_input_error(tmp_path, history, text):
    task = external_spec('{uccsd(t)-f12b;maxit=200}\n')['tasks'][0]
    (tmp_path/'sp.out').write_text(text)
    assert recovery.syntax_repair(task,tmp_path,
        dict(status='failed',failure_phase='execution'),history) is None


def test_targeted_drive_repairs_and_submits_once(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from kinbot import run_policy
    from kinbot.anl import protocol
    spec = external_spec('{uccsd(t)-f12b,scale_trip=1;maxit=200}\n')
    other = molpro_task('other', 'basis=cc-pVTZ-F12\nrhf\n{uccsd(t)-f12b,scale_trip=1}\n',
        geometry_from='initial', result_parser=dict(kind='molpro_energy',
            method='CCSD(T)-F12b', basis='cc-pVTZ-F12', reference='RHF'))
    spec['tasks'].append(other)
    graph = dispatch.prepare(spec, tmp_path/'graph')
    fail_external(graph, INPUT_ERROR)
    _, _, state = dispatch._load(graph)
    state['cc_iteration_retries'] = {'sp':[dict(limit=100),dict(limit=200)]}
    atomic_json(graph/'state.json',state)
    before=(graph/'tasks/other/other.inp').read_bytes()
    monkeypatch.setattr(dispatch,'_job_active',lambda job:str(job)=='91')
    monkeypatch.setattr(protocol,'preflight_ready',lambda _:True)
    submitted=[]
    def submit(command, **kwargs):
        kwargs['before_submit']()
        submitted.append(Path(kwargs['cwd']).name)
        return SimpleNamespace(returncode=0,stdout='91',stderr='')
    monkeypatch.setattr(run_policy,'submit',submit)
    for _ in range(2):
        dispatch.main(['drive',str(graph),'--once','--only','sp'])
    state=dispatch._load(graph)[2]
    assert submitted==['sp']
    assert state['tasks']['sp']['job_id']=='91'
    assert state['tasks']['other']['status']=='staged'
    assert (graph/'tasks/other/other.inp').read_bytes()==before
    assert [h['limit'] for h in state['cc_iteration_retries']['sp']]==[100,200]
