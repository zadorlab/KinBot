"""Native-output components require complete, hash-checked task records."""

import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

from kinbot.anl.dispatch import prepare
from kinbot.anl.model import IncompleteRecipeError
from kinbot.anl.recipes import recipe
from kinbot.anl.results import parse_result
from kinbot.anl.workflow import cbs_task_component, task_component
from tests.anl_fixture import dispatch_spec


_F12 = """basis=cc-pVTZ-F12
ccsd(t)-f12,scale_trip=1
 !CCSD(T)-F12b total energy -40.454906199189
 CCSD(T)-F12/cc-pVTZ-F12 energy= -40.454906199189
 Molpro calculation terminated
"""

_F12_QZ = """basis=cc-pVQZ-F12
ccsd(t)-f12,scale_trip=1
 !CCSD(T)-F12b total energy -40.456608306474
 CCSD(T)-F12/cc-pVQZ-F12 energy= -40.456608306474
 Molpro calculation terminated
"""


def _complete_task(run_dir, task, output):
    task_id = task['id']
    directory = run_dir / 'tasks' / task_id
    (directory / f'{task_id}.out').write_text(output)
    (directory / 'launcher.stdout').write_text('')
    (directory / 'launcher.stderr').write_text('')
    state_path = run_dir / 'state.json'
    state = json.loads(state_path.read_text())
    state['tasks'][task_id]['status'] = 'complete'
    state_path.write_text(json.dumps(state))
    artifacts = {
        name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
        for name in ('geometry.xyz', 'task.json', f'{task_id}.inp',
                     f'{task_id}.out', 'launcher.stdout', 'launcher.stderr')
    }
    (directory / 'execution.json').write_text(json.dumps({
        'schema': 1, 'task_id': task_id, 'status': 'executed',
        'geometry_sha256': state['tasks'][task_id]['geometry_sha256'],
        'artifacts': artifacts,
        'details': {'parsed_result': parse_result(output, task['result_parser'])},
    }))
    return directory


def _completed_f12(root):
    spec = dispatch_spec()
    task = next(task for task in spec['tasks'] if task['id'] == 'f12_tz')
    task['geometry_from'] = 'initial'
    spec['tasks'] = [task]
    spec_file = root / 'spec.json'
    spec_file.write_text(json.dumps(spec))
    run_dir = prepare(spec_file, root / 'run')
    directory = _complete_task(run_dir, task, _F12)
    return run_dir, directory


def _completed_f12_pair(root):
    spec = dispatch_spec()
    tasks = [task for task in spec['tasks']
             if task['id'] in ('f12_tz', 'f12_qz')]
    for task in tasks:
        task['geometry_from'] = 'initial'
    spec['tasks'] = tasks
    spec_file = root / 'spec.json'
    spec_file.write_text(json.dumps(spec))
    run_dir = prepare(spec_file, root / 'run')
    for task, output in zip(tasks, (_F12, _F12_QZ)):
        _complete_task(run_dir, task, output)
    return run_dir


def test_read_completed_f12_component_reparses_and_checks_artifacts():
    with TemporaryDirectory() as temporary:
        run_dir, directory = _completed_f12(Path(temporary))
        result = task_component(run_dir, 'f12_tz', key='f12_tz',
                                state_id='state-A')
        assert result.value_hartree == pytest.approx(-40.454906199189)
        assert result.method == 'CCSD(T)-F12b'
        assert result.settings['scale_trip'] == 1
        assert result.source_sha256 == hashlib.sha256(
            (directory / 'f12_tz.out').read_bytes()).hexdigest()

        (directory / 'f12_tz.out').write_text(_F12.replace('-40.454906199189',
                                                          '-40.454000000000'))
        with pytest.raises(RuntimeError, match='artifact'):
            task_component(run_dir, 'f12_tz', key='f12_tz',
                           state_id='state-A')


def test_old_unparsed_dispatch_result_cannot_enter_recipe():
    with TemporaryDirectory() as temporary:
        run_dir, directory = _completed_f12(Path(temporary))
        execution_path = directory / 'execution.json'
        execution = json.loads(execution_path.read_text())
        del execution['details']['parsed_result']
        execution_path.write_text(json.dumps(execution))
        with pytest.raises(ValueError, match='differs from native'):
            task_component(run_dir, 'f12_tz', key='reference_cbs',
                           state_id='state-A')
        state_path = run_dir / 'state.json'
        state = json.loads(state_path.read_text())
        state['tasks']['f12_tz']['status'] = 'submitted'
        state_path.write_text(json.dumps(state))
        with pytest.raises(IncompleteRecipeError, match='not complete'):
            task_component(run_dir, 'f12_tz', key='reference_cbs',
                           state_id='state-A')


def test_verified_f12_task_pair_builds_recipe_cbs_component():
    with TemporaryDirectory() as temporary:
        run_dir = _completed_f12_pair(Path(temporary))
        requirement = next(item for item in recipe('ANL0-F12').requirements
                           if item.key == 'reference_cbs')
        kwargs = {'requirement': requirement, 'state_id': 'state-A',
                  'lower_basis': 'cc-pVTZ-F12',
                  'upper_basis': 'cc-pVQZ-F12'}
        result = cbs_task_component(run_dir, 'f12_tz', 'f12_qz', **kwargs)
        assert result.value_hartree == pytest.approx(-40.457504545051066)
        assert result.backend == 'composite'
        assert result.basis == 'CBS(cc-pVTZ-F12,cc-pVQZ-F12)'
        assert result.settings['scale_trip'] == 1
        assert result.settings['upper_cardinal'] == 4
        assert len(result.source_sha256) == 64

        with pytest.raises(ValueError, match='basis order'):
            cbs_task_component(run_dir, 'f12_qz', 'f12_tz', **kwargs)
        with pytest.raises(ValueError, match='distinct'):
            cbs_task_component(run_dir, 'f12_tz', 'f12_tz', **kwargs)
        with pytest.raises(ValueError, match='basis pair'):
            cbs_task_component(run_dir, 'f12_tz', 'f12_qz',
                               **{**kwargs, 'upper_basis': 'cc-pV5Z-F12'})
        harmonic = next(item for item in recipe('ANL1').requirements
                        if item.key == 'harmonic_zpe')
        with pytest.raises(ValueError, match='electronic CBS'):
            cbs_task_component(run_dir, 'f12_tz', 'f12_qz',
                               **{**kwargs, 'requirement': harmonic})


def test_cbs_pair_rejects_incomplete_or_changed_native_output():
    with TemporaryDirectory() as temporary:
        run_dir = _completed_f12_pair(Path(temporary))
        requirement = next(item for item in recipe('ANL0-F12').requirements
                           if item.key == 'reference_cbs')
        kwargs = {'requirement': requirement, 'state_id': 'state-A',
                  'lower_basis': 'cc-pVTZ-F12',
                  'upper_basis': 'cc-pVQZ-F12'}
        state_path = run_dir / 'state.json'
        state = json.loads(state_path.read_text())
        state['tasks']['f12_qz']['status'] = 'submitted'
        state_path.write_text(json.dumps(state))
        with pytest.raises(IncompleteRecipeError, match='not complete'):
            cbs_task_component(run_dir, 'f12_tz', 'f12_qz', **kwargs)
        state['tasks']['f12_qz']['status'] = 'complete'
        state_path.write_text(json.dumps(state))
        output = run_dir / 'tasks' / 'f12_qz' / 'f12_qz.out'
        output.write_text(output.read_text().replace('-40.456608306474',
                                                    '-40.456000000000'))
        with pytest.raises(RuntimeError, match='artifact'):
            cbs_task_component(run_dir, 'f12_tz', 'f12_qz', **kwargs)
