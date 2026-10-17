import json
import os
import sys
from types import SimpleNamespace

import pytest
from kinbot.driver_status import run_recorded
from kinbot import pes


@pytest.mark.parametrize('result,status', [(None, 'complete'), (0, 'complete'), (1, 'failed')])
def test_explicit_controller_return_is_recorded(tmp_path, result, status):
    path = tmp_path/'status.json'
    assert run_recorded(lambda: result, path) == result
    record = json.loads(path.read_text())
    assert record['status'] == status
    assert record['pid'] == os.getpid()
    assert 'finished_utc' in record


def test_exception_preserves_actual_failure(tmp_path):
    path = tmp_path/'status.json'
    def fail():
        raise ValueError('bad input')
    with pytest.raises(ValueError):
        run_recorded(fail, path)
    assert json.loads(path.read_text())['error'] == 'ValueError: bad input'


@pytest.mark.parametrize('record,code', [({'status': 'failed', 'returncode': 1}, 1),
    ({'status': 'running'}, -9), ({'status': 'complete', 'returncode': 0}, 1)])
def test_pes_refuses_failed_or_interrupted_workers(tmp_path, monkeypatch, record, code):
    (tmp_path/'kinbot_status.json').write_text(json.dumps(dict(record, pid=42)))
    monkeypatch.setattr(pes, '_worker_processes', {42: SimpleNamespace(poll=lambda: code)})
    with pytest.raises(RuntimeError, match='postprocessing was not started'):
        pes.require_worker_success(str(tmp_path), 42)


def test_pes_reaps_its_child_and_requires_matching_success(tmp_path, monkeypatch):
    monkeypatch.setattr(pes, '_worker_processes', {42: SimpleNamespace(poll=lambda: 0)})
    assert pes.check_status(str(tmp_path), 42) == 0
    with pytest.raises(RuntimeError, match='without a readable completion'):
        pes.require_worker_success(str(tmp_path), 42)
    (tmp_path/'kinbot_status.json').write_text(json.dumps({'status': 'complete', 'pid': 42, 'returncode': 0}))
    pes.require_worker_success(str(tmp_path), 42)


def test_pes_launch_uses_own_interpreter_without_shell_background_token(tmp_path, monkeypatch):
    from kinbot.run_format import ensure_current_run
    monkeypatch.chdir(tmp_path)
    ensure_current_run('123', create=True)
    par = {name: '' for name in ('queue_template', 'q_temp_am1', 'q_temp_mp2',
        'q_temp_hi', 'q_temp_l3', 'single_point_template',
        'barrierless_saddle_single_point_template', 'vrc_tst_scan_molpro_tpl')}
    def popen(command, **kwargs):
        assert command == [sys.executable, '-m', 'kinbot.kb', '123.json']
        assert kwargs['cwd'] == '123'
        return SimpleNamespace(pid=42, poll=lambda: None)
    monkeypatch.setattr(pes.subprocess, 'Popen', popen)
    monkeypatch.setattr(pes.time, 'sleep', lambda _: None)
    monkeypatch.setattr(pes, '_worker_processes', {})
    assert pes.submit_job('123', par) == 42
