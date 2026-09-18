"""Native-output components require complete, hash-checked task records."""

import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

from kinbot.anl.dispatch import prepare
from kinbot.anl.model import IncompleteRecipeError
from kinbot.anl.results import parse_result
from kinbot.anl.workflow import task_component
from tests.anl_fixture import dispatch_spec


_F12 = """basis=cc-pVTZ-F12
ccsd(t)-f12,scale_trip=1
 !CCSD(T)-F12b total energy -40.454906199189
 CCSD(T)-F12/cc-pVTZ-F12 energy= -40.454906199189
 Molpro calculation terminated
"""


def _completed_f12(root):
    spec = dispatch_spec()
    task = next(task for task in spec['tasks'] if task['id'] == 'f12_tz')
    task['geometry_from'] = 'initial'
    spec['tasks'] = [task]
    spec_file = root / 'spec.json'
    spec_file.write_text(json.dumps(spec))
    run_dir = prepare(spec_file, root / 'run')
    directory = run_dir / 'tasks' / 'f12_tz'
    (directory / 'f12_tz.out').write_text(_F12)
    (directory / 'launcher.stdout').write_text('')
    (directory / 'launcher.stderr').write_text('')
    state_path = run_dir / 'state.json'
    state = json.loads(state_path.read_text())
    state['tasks']['f12_tz']['status'] = 'complete'
    state_path.write_text(json.dumps(state))
    artifacts = {
        name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
        for name in ('geometry.xyz', 'task.json', 'f12_tz.inp',
                     'f12_tz.out', 'launcher.stdout', 'launcher.stderr')
    }
    (directory / 'execution.json').write_text(json.dumps({
        'schema': 1, 'task_id': 'f12_tz', 'status': 'executed',
        'geometry_sha256': state['tasks']['f12_tz']['geometry_sha256'],
        'artifacts': artifacts,
        'details': {'parsed_result': parse_result(_F12, task['result_parser'])},
    }))
    return run_dir, directory


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
