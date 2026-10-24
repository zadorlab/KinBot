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
