from copy import deepcopy
import json
from pathlib import Path

import pytest

from kinbot.anl import dispatch, molpro_mpi_recovery as recovery
from kinbot.execution import atomic_json, file_sha256
from tests.anl_fixture import dispatch_spec


@pytest.fixture
def failed(tmp_path, monkeypatch):
    spec = dispatch_spec()
    task = next(t for t in spec['tasks'] if t['id'] == 'l3_geometry')
    task['geometry_from'] = 'initial'
    spec['tasks'] = [task]
    graph = dispatch.prepare(spec, tmp_path/'graph')
    state = dispatch._load(graph)[2]
    entry = state['tasks']['l3_geometry']
    entry.update(status='failed', job_id='42')
    d = graph/'tasks/l3_geometry'
    atomic_json(d/'execution.json', dict(schema=1, task_id=task['id'], status='failed',
        geometry_sha256=entry['geometry_sha256'], slurm_job_id='42', failure_phase='execution',
        error='Molpro force evaluation l3_geometry_step_0025 exited with status 255; '
              'diagnostics: l3_geometry_step_0025.failure.json.'))
    atomic_json(d/'l3_geometry_step_0025.failure.json', dict(returncode=255, host='node31',
        explicit_io_failure=False, files={'work/error.stderr': {'tail':
            '[5] Received an Error in Communication: (1)\napplication called MPI_Abort'}}))
    atomic_json(graph/'state.json', state)
    monkeypatch.setattr(dispatch, '_job_active', lambda job: False)
    return graph


def test_once_resource_only_transaction_archives_and_restages(failed):
    original = dispatch._load(failed)[1]['tasks'][0]
    geom = (failed/'tasks/l3_geometry/geometry.xyz').read_bytes()
    result = dispatch.advance(failed)
    entry = result['tasks']['l3_geometry']
    assert entry['status'] == 'staged' and entry['attempt'] == 2
    assert (Path(entry['previous_attempt'])/'execution.json').exists()
    task = dispatch._load(failed)[1]['tasks'][0]
    assert task['resources']['cores'] == 1
    assert task['profile'] == dict(original['profile'], calculator_kwargs=original['profile'].get('calculator_kwargs', {}))
    assert task['optimizer'] == original['optimizer']
    assert (failed/'tasks/l3_geometry/geometry.xyz').read_bytes() == geom
    from kinbot.manual_fallback import logical_task
    assert logical_task(task) == original
    assert dispatch.advance(failed)['molpro_serial_retries'] == result['molpro_serial_retries']


@pytest.mark.parametrize('reason', ['oom', 'io', 'active', 'budget', 'unproven', 'wrong_generation', 'singular'])
def test_does_not_retry_unproven_or_other_failures(failed, monkeypatch, reason):
    graph, spec, state = dispatch._load(failed)
    d = graph/'tasks/l3_geometry'
    p = d/'l3_geometry_step_0025.failure.json'
    data = json.loads(p.read_text())
    if reason == 'oom': data['files']['work/error.stderr']['tail'] += '\nout of memory'
    if reason == 'io': data['explicit_io_failure'] = True
    if reason == 'singular':
        data['files']['work/error.stderr']['error_context'] = 'The expansion set has become singular'
    if reason == 'unproven': data['files']['work/error.stderr']['tail'] = 'KILLED BY SIGNAL: 9'
    if reason == 'active': monkeypatch.setattr(dispatch, '_job_active', lambda job: True)
    if reason == 'budget': spec['intent'] = {'vrc_fragment_retry': True}
    if reason == 'wrong_generation':
        e = json.loads((d/'execution.json').read_text()); e['geometry_sha256'] = 'wrong'
        atomic_json(d/'execution.json', e)
    atomic_json(p, data)
    before = file_sha256(graph/'workflow.json')
    if reason == 'wrong_generation':
        with pytest.raises(RuntimeError): recovery.recover(graph, spec, state)
    else:
        assert recovery.recover(graph, spec, state) == (spec, state)
    assert file_sha256(graph/'workflow.json') == before
