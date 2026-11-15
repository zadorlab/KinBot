"""New control roots require a complete, bounded ownership check."""
import json
from pathlib import Path

import pytest

from kinbot import live_control as control, run_control


@pytest.mark.parametrize('limit', ['directories', 'entries', 'seconds'])
def test_incomplete_discovery_never_enables_controls(tmp_path, monkeypatch, limit):
    (tmp_path/'checkpoint').mkdir()
    (tmp_path/'checkpoint/data').write_text('retained calculation')
    if limit == 'directories':
        monkeypatch.setattr(control, '_CONTROL_SCAN_DIRECTORIES', 1)
    elif limit == 'entries':
        monkeypatch.setattr(control, '_CONTROL_SCAN_ENTRIES', 1)
    else:
        ticks = iter([0., 0., 16.])
        monkeypatch.setattr(control.time, 'monotonic', lambda: next(ticks))
    with pytest.raises(ValueError, match='discovery limit.*'):
        control.configure(tmp_path)
    assert not (tmp_path/control.NAME).exists()
    assert not (tmp_path/'.kinbot_controls').exists()
    assert not (tmp_path/'.kinbot_live_control.lock').exists()
    assert (tmp_path/'checkpoint/data').read_text() == 'retained calculation'


def test_hidden_child_control_is_still_detected(tmp_path):
    child = tmp_path/'.hidden/owner'; child.mkdir(parents=True)
    control.configure(child, hold_failures=True)
    saved = (child/control.NAME).read_bytes()
    with pytest.raises(ValueError, match='Child live controls already exist'):
        control.configure(tmp_path)
    assert not (tmp_path/control.NAME).exists()
    assert (child/control.NAME).read_bytes() == saved


def test_unreadable_directory_fails_closed(tmp_path, monkeypatch):
    child = tmp_path/'child'; child.mkdir()
    scan = control.os.scandir
    def unreadable(path):
        if Path(path) == child:
            raise PermissionError('permission denied')
        return scan(path)
    monkeypatch.setattr(control.os, 'scandir', unreadable)
    with pytest.raises(ValueError, match='Cannot verify live-control ownership'):
        control.configure(tmp_path)
    assert not (tmp_path/control.NAME).exists()
    assert not (tmp_path/'.kinbot_controls').exists()


def test_directory_symlinks_do_not_expand_control_scope(tmp_path):
    other = tmp_path/'other'; other.mkdir(); control.configure(other)
    root = tmp_path/'root'; root.mkdir()
    (root/'linked').symlink_to(other, target_is_directory=True)
    control.configure(root)
    assert control.find(root) == root/control.NAME
    assert control.find(root/'linked') == other/control.NAME


def test_existing_root_preserves_settings_without_descendant_scan(tmp_path, monkeypatch):
    control.configure(tmp_path, hold_failures=True, hold_task_failures=True)
    control.resources(tmp_path, 'small', ['graph/*'], {'cores': 2})
    before = control.read(tmp_path)
    def unexpected(_):
        pytest.fail('Existing control root must not rescan checkpoint history.')
    monkeypatch.setattr(control, '_check_control_descendants', unexpected)
    after = control.configure(tmp_path)
    assert after['hold_failures'] and after['hold_task_failures']
    assert after['rules'] == before['rules']
    assert after['revision'] == before['revision'] + 1
    child = tmp_path/'child'; child.mkdir()
    with pytest.raises(ValueError, match='Use the existing live-control root'):
        control.configure(child)
    assert not (child/control.NAME).exists()


def test_ancestor_enabled_during_discovery_prevents_nested_root(tmp_path, monkeypatch):
    child = tmp_path/'child'; child.mkdir()
    scan = control._check_control_descendants
    def concurrent_ancestor(root):
        scan(root)
        monkeypatch.setattr(control, '_check_control_descendants', scan)
        control.configure(tmp_path)
    monkeypatch.setattr(control, '_check_control_descendants', concurrent_ancestor)
    with pytest.raises(ValueError, match='Use the existing live-control root'):
        control.configure(child)
    assert not (child/control.NAME).exists()
    assert control.read(tmp_path)['revision'] == 1


def test_enable_reports_scope_before_discovery_and_keeps_json_stdout(tmp_path, monkeypatch, capsys):
    scan = control._check_control_descendants
    def observed(root):
        output = capsys.readouterr()
        assert str(tmp_path) in output.err
        assert 'bounded' in output.err
        assert output.out == ''
        scan(root)
    monkeypatch.setattr(control, '_check_control_descendants', observed)
    run_control.main(['enable', str(tmp_path)])
    assert json.loads(capsys.readouterr().out)['revision'] == 1
