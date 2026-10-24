from copy import deepcopy
import json
from pathlib import Path

import pytest
from ase import Atoms

from kinbot.anl import dispatch, molpro_solver_recovery as recovery
from kinbot.ase_modules.calculators.molpro import render_input
from kinbot.execution import atomic_json
from tests.anl_fixture import dispatch_spec


def fail(graph, job='42'):
    graph, spec, state = dispatch._load(graph)
    entry = state['tasks']['l3_geometry']
    entry.update(status='failed', job_id=job)
    d = graph/'tasks/l3_geometry'
    atomic_json(d/'execution.json', dict(schema=1, task_id='l3_geometry', status='failed',
        geometry_sha256=entry['geometry_sha256'], slurm_job_id=job, failure_phase='execution',
        error='Molpro force evaluation l3_geometry_step_0006 exited with status 1; '
              'diagnostics: l3_geometry_step_0006.failure.json.'))
    atomic_json(d/'l3_geometry_step_0006.failure.json', dict(returncode=1, files={
        'work/native.out': {'error_context': 'Number of p-space configurations: 3\n'
            'Catastrophic failure in diagonalization(hsdel)\nThe expansion set has become singular',
            'tail': 'GLOBAL ERROR fehler on processor 0'}}))
    atomic_json(graph/'state.json', state)


@pytest.fixture
def failed(tmp_path, monkeypatch):
    spec = dispatch_spec()
    task = next(t for t in spec['tasks'] if t['id'] == 'l3_geometry')
    task['geometry_from'] = 'initial'
    task['profile']['method'] = 'caspt2-f12(4,4)'
    task['profile']['calculator_kwargs'] = dict(caspt2_shift=.1, caspt2_ipea=.25,
        state_model={'state_average': [{'multiplicity': 1, 'roots': 3, 'weights': [1/3]*3}],
                     'selected_roots': {'1': 2}, 'canonical_orbitals': True,
                     'correlation_options': {'reference_policy': 'state_average',
                                             'max_reference_vectors': 60}})
    spec['tasks'] = [task]
    graph = dispatch.prepare(spec, tmp_path/'graph')
    fail(graph)
    monkeypatch.setattr(dispatch, '_job_active', lambda job: False)
    return graph


@pytest.mark.parametrize('bounded', [False, True])
def test_two_bounded_retries_keep_model_and_archive(failed, bounded):
    if bounded:
        _, spec, state = dispatch._load(failed)
        spec['intent'] = {'vrc_fragment_retry': True}
        atomic_json(failed/'workflow.json', spec)
        from kinbot.execution import file_sha256
        state['spec_sha256'] = file_sha256(failed/'workflow.json')
        atomic_json(failed/'state.json', state)
    original = deepcopy(dispatch._load(failed)[1]['tasks'][0])
    for index in range(2):
        result = dispatch.advance(failed)
        entry = result['tasks']['l3_geometry']
        assert entry['status'] == 'staged'
        assert (Path(entry['previous_attempt'])/'execution.json').is_file()
        task = dispatch._load(failed)[1]['tasks'][0]
        before = original['profile']['calculator_kwargs']
        after = task['profile']['calculator_kwargs']
        assert after['caspt2_shift'] == before['caspt2_shift']
        assert after['caspt2_ipea'] == before['caspt2_ipea']
        model = deepcopy(after['state_model'])
        controls = model.pop('correlation_options')
        unchanged = deepcopy(before['state_model']); unchanged.pop('correlation_options')
        assert model == unchanged
        assert controls['reference_policy'] == 'state_average'
        assert controls['pspace'] == dict(energy_threshold=0., minimum_configurations=6)
        if index:
            assert controls['max_reference_vectors'] == 30
            assert controls['energy_convergence'] == 1.e-10
            assert controls['coefficient_convergence'] == 1.e-7
        assert task['optimizer'] == original['optimizer']
        assert task['resources'] == original['resources']
        from kinbot.manual_fallback import logical_task
        assert logical_task(task) == original
        fail(failed, str(43+index))
    result = dispatch.advance(failed)
    assert result['tasks']['l3_geometry']['status'] == 'failed'
    assert len(result['molpro_pspace_retries']['l3_geometry']) == 2


@pytest.mark.parametrize('fault', ['active', 'other', 'io', 'pspace', 'generation'])
def test_no_unproven_or_unrequested_retry(failed, monkeypatch, fault):
    graph, spec, state = dispatch._load(failed)
    d = graph/'tasks/l3_geometry'
    p = d/'l3_geometry_step_0006.failure.json'
    data = json.loads(p.read_text())
    if fault == 'active': monkeypatch.setattr(dispatch, '_job_active', lambda job: True)
    elif fault == 'other': data['files']['work/native.out']['error_context'] = 'MPI_Abort'
    elif fault == 'io': data['explicit_io_failure'] = True
    elif fault == 'pspace':
        spec['tasks'][0]['profile']['calculator_kwargs']['state_model']['correlation_options']['pspace'] = {
            'energy_threshold': 0., 'minimum_configurations': 10}
    else:
        e = json.loads((d/'execution.json').read_text()); e['geometry_sha256'] = 'wrong'
        atomic_json(d/'execution.json', e)
    atomic_json(p, data)
    if fault == 'generation':
        with pytest.raises(RuntimeError): recovery.recover(graph, spec, state)
    else:
        assert recovery.recover(graph, spec, state) == (spec, state)


def test_full_native_output_can_supply_signature_missing_from_short_excerpt(failed):
    d = failed/'tasks/l3_geometry'
    p = d/'l3_geometry_step_0006.failure.json'
    data = json.loads(p.read_text())
    full = data['files']['work/native.out']['error_context']
    data['files']['work/native.out']['error_context'] = 'Otherwise, try changing P space'
    atomic_json(p, data)
    (d/'l3_geometry_step_0006.out').write_text(full)
    assert dispatch.advance(failed)['tasks']['l3_geometry']['status'] == 'staged'


@pytest.mark.parametrize('method', ['caspt2-f12(4,4)', 'mrci-f12+q(4,4)'])
def test_native_pspace_card_does_not_change_root_or_nstati(failed, method):
    dispatch.advance(failed)
    task = dispatch._load(failed)[1]['tasks'][0]
    atoms = Atoms('CH4', positions=dispatch_spec()['molecule']['positions'])
    text = render_input(atoms, basis='cc-pVDZ-F12', geometry_name='geom.xyz',
                        method=method, mult=1,
                        **task['profile']['calculator_kwargs'])
    assert 'pspace,0,6;' in text
    assert 'state,1,2;refstate,3,1,2,3;' in text
    if method.startswith('caspt2'):
        assert 'nstati=' not in text
        assert 'shift=0.1,ipea=0.25' in text
    else:
        assert 'option,nstati=3;' in text
