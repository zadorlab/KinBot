"""The example transition must preserve results and refuse duplicate launches."""
import importlib.util
import json
from pathlib import Path

import pytest

PACKAGE = Path(__file__).parents[1]/'examples/anl/caspt2_production'


@pytest.fixture
def handoff(monkeypatch):
    monkeypatch.syspath_prepend(str(PACKAGE))
    spec = importlib.util.spec_from_file_location('caspt2_handoff', PACKAGE/'handoff.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_prepared_package_is_portable_immutable_and_records_resolved_settings(handoff, tmp_path):
    root = handoff.prepare(tmp_path/'production')
    assert handoff.prepare(root) == root
    record = handoff.verify(root)
    for case in handoff.CASES:
        resolved = json.loads((root/case/'resolved_parameters.json').read_text())
        assert resolved['epsilon_unit'] == 'cm-1'
        assert resolved['queue_job_limit'] == 10
        assert str(Path(resolved['vrc_tst_surface_model_file']).parent) == str(root/case)
        assert f'{case}/resolved_parameters.json' in record['sha256']
    model = root/'assets/ch3s_cs_state_model.json'
    model.write_text('{}')
    with pytest.raises(ValueError, match='Prepared artifact changed'):
        handoff.verify(root)
    with pytest.raises(ValueError, match='Existing prepared artifact differs'):
        handoff.prepare(root)


def test_retirement_inventories_first_and_preserves_unrelated_runs(handoff, monkeypatch, tmp_path):
    from kinbot import run_cleanup, orbital_campaign
    root = tmp_path/'production'
    root.mkdir()
    parent = tmp_path/'old'
    parent.mkdir()
    expected = [parent/name for name in handoff.OLD_NAMES[:2]]
    for old, marker in zip(expected, ('campaign.json', 'confirmation.json')):
        old.mkdir()
        (old/marker).write_text('{}')
        (old/'native.out').write_text('retained evidence')
    unrelated = parent/'unrelated_running_calculation'
    unrelated.mkdir()
    (unrelated/'native.out').write_text('unrelated evidence')
    events = []
    monkeypatch.setattr(run_cleanup, 'inspect', lambda path: events.append(('inspect', path)) or {})
    monkeypatch.setattr(orbital_campaign, 'stop_controller', lambda path: events.append(('controller', path)))

    def stop(path, *, archive, timeout):
        assert archive and timeout == 300
        events.append(('stop', path))
        destination = path.with_name(path.name+'.stopped.test')
        path.rename(destination)
        return {'archive': str(destination)}

    monkeypatch.setattr(run_cleanup, 'stop', stop)
    reports = handoff.retire(root, parent)
    assert events[:2] == [('inspect', path) for path in expected]
    assert [path for action, path in events if action == 'stop'] == expected
    assert (unrelated/'native.out').read_text() == 'unrelated evidence'
    for report in reports:
        assert (Path(report['archive'])/'native.out').read_text() == 'retained evidence'
    assert len(list(root.glob('retirement_inventory.*.json'))) == 1
    assert len(list(root.glob('retirement_results.*.json'))) == 1


def test_unrecognized_old_run_aborts_before_any_stop(handoff, monkeypatch, tmp_path):
    from kinbot import run_cleanup
    root = tmp_path/'production'
    root.mkdir()
    bad = tmp_path/handoff.OLD_NAMES[0]
    bad.mkdir()
    monkeypatch.setattr(run_cleanup, 'stop', lambda *a, **k: pytest.fail('must not stop'))
    with pytest.raises(ValueError, match='Unrecognized old campaign'):
        handoff.retire(root, tmp_path)


def test_previous_launch_refuses_duplicate_even_if_queue_empty(handoff, monkeypatch, tmp_path):
    from kinbot import run_cleanup
    case = tmp_path/'ch3s'
    case.mkdir()
    (case/'launch.json').write_text('{}')
    monkeypatch.setattr(handoff, 'verify', lambda root: {'inputs': [str(case/'input.json')]})
    monkeypatch.setattr(run_cleanup, 'inspect', lambda path: {'jobs': [], 'processes': []})
    with pytest.raises(RuntimeError, match='already started'):
        handoff.assert_not_started(tmp_path)
