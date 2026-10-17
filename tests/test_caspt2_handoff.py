"""The example transition must preserve results and refuse duplicate launches."""
import importlib.util
import json
import os
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


def test_rendered_preflight_uses_both_real_production_profiles(handoff, tmp_path, monkeypatch):
    monkeypatch.setenv('SBATCH_EXCLUDE', 'existing-exclusion')
    root = handoff.prepare(tmp_path/'production')
    result = handoff.render_gaussian_inputs(root, handoff.verify(root))
    assert os.environ['SBATCH_EXCLUDE'] == 'existing-exclusion'
    assert len(result['files']) == 36
    for item in result['files']:
        text = (root/item['path']).read_text()
        assert 'B2PLYP/cc-pVTZ' in text
        assert 'scf(xqc,conver=10)' in text
        assert 'scf(xqc,(' not in text
        assert 'EmpiricalDispersion(GD3BJ)' in text
        assert 'integral(UltraFine)' in text
        assert '\nopt(' not in text and '\nirc(' not in text
        if item['path'].endswith('_energy.com'):
            assert '\nfreq' not in text
    assert not (root/'ch3s/launch.json').exists()


@pytest.mark.parametrize('stop_fails', [False, True])
def test_clean_restart_stops_both_before_archiving_parent(handoff, monkeypatch, tmp_path, stop_fails):
    from kinbot import run_cleanup
    old = handoff.prepare(tmp_path/'previous')
    root = handoff.prepare(tmp_path/'replacement')
    (old/'ch3s/evidence.log').write_text('failed native input evidence')
    events = []
    def inspect(path):
        events.append(('inspect', path.name))
        return {'directory': str(path), 'jobs': [], 'processes': []}
    def stop(path, **kwargs):
        events.append(('stop', path.name))
        assert kwargs == {'archive': False, 'timeout': 300.}
        if stop_fails and path.name == 'ch3ss':
            raise RuntimeError('scheduler unavailable')
        assert old.exists()
        return {'status': 'stopped'}
    monkeypatch.setattr(run_cleanup, 'inspect', inspect)
    monkeypatch.setattr(run_cleanup, 'stop', stop)
    if stop_fails:
        with pytest.raises(RuntimeError, match='scheduler unavailable'):
            handoff.retire_production(root, old)
        assert old.exists()
        assert not (root/'ch3s/launch.json').exists()
    else:
        result = handoff.retire_production(root, old)
        assert not old.exists()
        assert (Path(result['archive'])/'ch3s/evidence.log').read_text() == 'failed native input evidence'
        assert events == [('inspect', 'ch3s'), ('inspect', 'ch3ss'),
                          ('stop', 'ch3s'), ('stop', 'ch3ss'),
                          ('inspect', 'ch3s'), ('inspect', 'ch3ss')]


def test_cleanup_refuses_unrecognized_or_overlapping_production(handoff, monkeypatch, tmp_path):
    from kinbot import run_cleanup
    old = handoff.prepare(tmp_path/'previous')
    root = handoff.prepare(tmp_path/'replacement')
    monkeypatch.setattr(run_cleanup, 'stop', lambda *a, **k: pytest.fail('must not stop'))
    with pytest.raises(ValueError, match='separate directories'):
        handoff.retire_production(old, old)
    data = json.loads((old/'prepared.json').read_text())
    data['inputs'][0] = '/unrelated/input.json'
    (old/'prepared.json').write_text(json.dumps(data))
    with pytest.raises(ValueError, match='Unrecognized previous'):
        handoff.retire_production(root, old)
