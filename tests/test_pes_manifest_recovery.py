"""Storage faults must retain coordination and cannot authorize chemistry."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import pytest
from kinbot import pes_phase as phase
from kinbot.execution import atomic_json
from kinbot.input_refresh import owner_directory


def test_atomic_json_concurrent_publishers_never_share_temporary(tmp_path):
    path=tmp_path/'record.json'
    payloads=[dict(writer=i,data=str(i)*10000) for i in range(20)]
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(lambda data:atomic_json(path,data),payloads))
    assert json.loads(path.read_text()) in payloads
    assert not list(tmp_path.glob('*.tmp'))


def test_failed_flush_preserves_previous_json(tmp_path, monkeypatch):
    path=atomic_json(tmp_path/'record.json',dict(epoch=1))
    def fail(*a):raise OSError('simulated storage fault')
    monkeypatch.setattr('kinbot.execution.os.fsync',fail)
    with pytest.raises(OSError):atomic_json(path,dict(epoch=2))
    assert json.loads(path.read_text())==dict(epoch=1)
    assert not list(tmp_path.glob('*.tmp'))


def test_manifest_corruption_recovers_same_epoch_and_preserves_bytes(tmp_path):
    m=phase.initialize(tmp_path)
    phase.register(tmp_path,'well',123)
    (tmp_path/phase.NAME).write_text('')
    recovered=phase.read(tmp_path)
    assert recovered['token']==m['token'] and recovered['workers']['well']['pid']==123
    phase.tick(tmp_path,['well'],max_searches=1)
    assert phase.read(tmp_path)['workers']['well']['status']=='exploring'
    archives=list(tmp_path.glob('.pes_exploration.json.corrupt.*'))
    assert len(archives)==1 and archives[0].read_bytes()==b''


def test_legacy_empty_manifest_allows_locked_native_driver_startup(tmp_path,monkeypatch):
    source=atomic_json(tmp_path/'input.json',dict(smiles='CC'))
    (tmp_path/phase.NAME).write_text('')
    monkeypatch.chdir(tmp_path)
    assert owner_directory(source) is None
    from kinbot.kb import _driver_lock
    with _driver_lock():
        phase.initialize(tmp_path,source,'parent')
    assert phase.read(tmp_path)['initial_worker']=='parent'
    assert len(list(tmp_path.glob('.pes_exploration.json.corrupt.*')))==1


def test_corruption_without_backup_is_not_an_empty_exploration(tmp_path):
    (tmp_path/phase.NAME).write_text('')
    with pytest.raises(phase.ManifestUnavailable):phase.read(tmp_path)
