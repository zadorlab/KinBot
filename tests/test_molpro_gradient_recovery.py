from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write

from kinbot.anl import dispatch, molpro_gradient_recovery as recovery
from kinbot.execution import atomic_json
from tests.test_molpro_solver_recovery import failed


def stagnate(graph, job='42'):
    graph, spec, state = dispatch._load(graph)
    ident = 'l3_geometry'
    directory = graph/'tasks'/ident
    task = spec['tasks'][0]
    entry = state['tasks'][ident]
    entry.update(status='failed', job_id=job)
    atomic_json(directory/'execution.json', dict(schema=1, task_id=ident, status='failed',
        geometry_sha256=entry['geometry_sha256'], slurm_job_id=job, failure_phase='execution',
        error='ASE geometry optimization did not converge.'))
    atomic_json(directory/'optimization_budget.json', dict(status='exhausted', chunk_steps=100,
        chunks_used=3, fmax=task['optimizer']['fmax']))
    (directory/'optimization.log').write_text(''.join(
        f'Sella {step} 12:00:00 -40.000000 0.3600 0.0000 0.0001 -20.0 1\n'
        for step in range(281, 301)))
    atoms = read(directory/'geometry.xyz')
    from kinbot.anl.optimization_recovery import load_seed
    load_seed(directory, json.loads((directory/'task.json').read_text()), atoms)
    atoms.calc = SinglePointCalculator(atoms, energy=-40., forces=np.full((len(atoms), 3), .36))
    write(directory/'optimization.traj', [atoms, atoms])
    atomic_json(graph/'state.json', state)


def test_owner_archives_and_retries_once_with_same_science(failed):
    from kinbot.manual_fallback import logical_task
    from kinbot.vrc_graph import identity
    original = deepcopy(dispatch._load(failed)[1])
    # This is a saved exhausted optimization, not a native eigensolver failure.
    stagnate(failed)
    result = dispatch.advance(failed)
    task = dispatch._load(failed)[1]['tasks'][0]
    entry = result['tasks']['l3_geometry']
    assert entry['status'] == 'staged'
    assert (Path(entry['previous_attempt'])/'execution.json').is_file()
    assert (Path(entry['previous_attempt'])/'optimization.traj').is_file()
    assert logical_task(task) == original['tasks'][0]
    assert identity(dispatch._load(failed)[1]) == identity(original)
    assert task['optimizer'] == original['tasks'][0]['optimizer']
    before = original['tasks'][0]['profile']['calculator_kwargs']
    after = task['profile']['calculator_kwargs']
    for key in ('caspt2_shift', 'caspt2_ipea'):
        assert before[key] == after[key]
    model = deepcopy(after['state_model'])
    for group, fields in (('casscf_options', ('energy_convergence', 'gradient_convergence')),
                          ('correlation_options', ('energy_convergence', 'coefficient_convergence'))):
        for key in fields:
            model[group].pop(key)
        if not model[group]:
            model.pop(group)
    assert model == before['state_model']
    assert after['numerical_gradient']['stencil'] == 'fourpoint'
    assert after['state_model']['casscf_options']['gradient_convergence'] == 1.e-7
    assert result['optimizer_continuations']['l3_geometry']['chunks_used'] == 0
    stagnate(failed, '43')
    second = dispatch.advance(failed)
    assert second['tasks']['l3_geometry']['status'] == 'failed'
    assert len(second['molpro_gradient_retries']['l3_geometry']) == 1


@pytest.mark.parametrize('fault', ['active', 'other_error', 'progressing', 'large_radius',
                                  'short_log', 'not_exhausted', 'generation', 'trajectory',
                                  'bounded_probe', 'already_tight'])
def test_requires_specific_proven_stagnation(failed, monkeypatch, fault):
    stagnate(failed)
    graph, spec, state = dispatch._load(failed)
    directory = graph/'tasks/l3_geometry'
    if fault == 'active':
        monkeypatch.setattr(dispatch, '_job_active', lambda job: True)
    elif fault in ('other_error', 'generation'):
        p = directory/'execution.json'; value = json.loads(p.read_text())
        value['error' if fault == 'other_error' else 'geometry_sha256'] = 'different'
        atomic_json(p, value)
    elif fault in ('progressing', 'large_radius', 'short_log'):
        p = directory/'optimization.log'; text = p.read_text()
        if fault == 'progressing':
            text = '\n'.join(s.replace('0.3600', '0.0600') if i > 14 else s
                             for i, s in enumerate(text.splitlines()))
        elif fault == 'large_radius': text = text.replace('0.0001', '0.1000')
        else: text = '\n'.join(text.splitlines()[-5:])
        p.write_text(text)
    elif fault == 'not_exhausted':
        p = directory/'optimization_budget.json'; value = json.loads(p.read_text())
        value['status'] = 'running'; atomic_json(p, value)
    elif fault == 'trajectory':
        atoms = read(directory/'optimization.traj', -1); atoms.positions += 1
        write(directory/'optimization.traj', [atoms, atoms])
    elif fault == 'bounded_probe': spec['intent'] = {'vrc_fragment_retry': True}
    else:
        spec['tasks'][0]['profile']['calculator_kwargs'].update(
            recovery.accuracy_options(spec['tasks'][0])['calculator_kwargs'])
    if fault in ('generation', 'trajectory'):
        with pytest.raises(RuntimeError): recovery.recover(graph, spec, state)
    else:
        assert recovery.recover(graph, spec, state) == (spec, state)
    assert not list((graph/'attempts').glob('live_*'))


@pytest.mark.parametrize('dependent', [False, True])
def test_active_other_work_is_preserved_or_blocks_source_retry(failed, tmp_path, dependent):
    from tests.anl_fixture import molpro_task
    spec = deepcopy(dispatch._load(failed)[1])
    spec['tasks'].append(molpro_task('other', 'rhf\n',
        geometry_from='l3_geometry' if dependent else 'initial'))
    graph = dispatch.prepare(spec, tmp_path/'with_other')
    stagnate(graph)
    graph, spec, state = dispatch._load(graph)
    state['tasks'].setdefault('other', {}).update(status='submitted', job_id='99')
    atomic_json(graph/'state.json', state)
    before = deepcopy(state['tasks']['other'])
    replacement, result = recovery.recover(graph, spec, state)
    assert result['tasks']['other'] == before
    assert replacement['tasks'][1] == spec['tasks'][1]
    if dependent:
        assert result['tasks']['l3_geometry']['status'] == 'failed'
        assert not result.get('molpro_gradient_retries')
    else:
        assert 'l3_geometry' not in result['tasks']  # Cleared for normal staging by advance.
        assert len(result['molpro_gradient_retries']['l3_geometry']) == 1
