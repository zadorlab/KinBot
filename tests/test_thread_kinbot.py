import subprocess
import sys
from types import SimpleNamespace

import pytest
from kinbot import thread_kinbot as launcher


@pytest.mark.parametrize('pes,module', [(False, 'kinbot.kb'), (True, 'kinbot.pes')])
def test_launcher_uses_current_python_and_files_not_undrained_pipes(tmp_path, monkeypatch, pes, module):
    def popen(command, **kwargs):
        assert command == [sys.executable, '-m', module, 'input file.json']
        assert kwargs['cwd'] == tmp_path
        assert kwargs['stdin'] == subprocess.DEVNULL
        kwargs['stdout'].write(b'worker progress\n')
        kwargs['stderr'].write(b'diagnostic\n')
        return SimpleNamespace(pid=42, poll=lambda: 0)
    monkeypatch.setattr(launcher.subprocess, 'Popen', popen)
    monkeypatch.setattr(launcher, '_processes', {})
    assert launcher.submit_job(str(tmp_path), 'input file.json', pes) == 42
    assert launcher.check_status(str(tmp_path), 42) == 0
    assert (tmp_path/'thread_kinbot.out').read_text() == 'worker progress\n'


def test_failure_not_counted_as_finished(tmp_path, monkeypatch):
    monkeypatch.setattr(launcher, '_processes', {42: SimpleNamespace(poll=lambda: 1)})
    with pytest.raises(RuntimeError, match='exited with status 1'):
        launcher.check_status(str(tmp_path), 42)
    with pytest.raises(ValueError, match='positive integer'):
        launcher.run_threads({}, 'invalid', max_running=0)
