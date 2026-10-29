"""Ordinary launch registration is metadata-only and recoverable."""
import json
import socket

import pytest

from kinbot import run_adopt
from kinbot.execution import atomic_json


@pytest.fixture
def ordinary(tmp_path):
    root = tmp_path/'run'
    for name, pid, state in [('first', 123, 'complete'), ('second', 456, 'running')]:
        case = root/name
        case.mkdir(parents=True)
        atomic_json(case/'input.json', {'theory': 'unchanged'})
        atomic_json(case/'pes_status.json', dict(schema=1, pid=pid, status=state,
                    host=socket.gethostname(), started_utc='2026-01-01T00:00:00+00:00'))
        (case/'launch.pid').write_text(str(pid))
        (case/'native.log').write_text('native evidence')
    return root


def test_preview_writes_nothing_and_registration_preserves_native_files(ordinary):
    files = {p: p.read_bytes() for p in ordinary.rglob('*') if p.is_file()}
    assert run_adopt.register(ordinary)['phase'] == 'reviewed'
    assert {p: p.read_bytes() for p in ordinary.rglob('*') if p.is_file()} == files
    receipt = run_adopt.register(ordinary, apply=True)
    assert receipt['phase'] == 'registered'
    assert {p: p.read_bytes() for p in files} == files
    assert json.loads((ordinary/'salvage.json').read_text())['code'] is None
    assert run_adopt.register(ordinary, apply=True) == receipt


def test_partial_registration_can_resume(ordinary, monkeypatch):
    def interrupted(path, data):
        if path.name == 'salvage.json':
            raise OSError('interrupted metadata publication')
        atomic_json(path, data)
    monkeypatch.setattr(run_adopt, 'atomic_json', interrupted)
    with pytest.raises(OSError):
        run_adopt.register(ordinary, apply=True)
    assert (ordinary/'prepared.json').exists()
    monkeypatch.setattr(run_adopt, 'atomic_json', atomic_json)
    assert run_adopt.register(ordinary, apply=True)['phase'] == 'registered'


@pytest.mark.parametrize('bad', ['host', 'pid', 'metadata', 'escaping_input'])
def test_registration_refuses_inconsistent_launch(ordinary, bad, tmp_path):
    if bad == 'host':
        p = ordinary/'first/pes_status.json'
        atomic_json(p, dict(json.loads(p.read_text()), host='unrelated-host'))
    elif bad == 'pid':
        (ordinary/'first/launch.pid').write_text('999')
    elif bad == 'metadata':
        atomic_json(ordinary/'prepared.json', {})
    else:
        outside = tmp_path/'external.json'
        outside.write_text('{}')
        (ordinary/'first/input.json').unlink()
        (ordinary/'first/input.json').symlink_to(outside)
    with pytest.raises(ValueError):
        run_adopt.register(ordinary, apply=True)
    assert not (ordinary/'ordinary_adoption.json').exists()


def test_partial_registration_refuses_changed_input(ordinary, monkeypatch):
    def interrupted(path, data):
        if path.name == 'prepared.json':
            raise OSError('interrupted')
        atomic_json(path, data)
    monkeypatch.setattr(run_adopt, 'atomic_json', interrupted)
    with pytest.raises(OSError):
        run_adopt.register(ordinary, apply=True)
    (ordinary/'first/input.json').write_text('{"theory":"different"}')
    with pytest.raises(ValueError, match='Input changed'):
        run_adopt.register(ordinary, apply=True)
