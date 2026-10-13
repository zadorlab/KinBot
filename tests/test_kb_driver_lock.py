import fcntl
from pathlib import Path

import pytest

from kinbot import kb
from kinbot.run_format import MARKER, ensure_current_run


def test_second_driver_is_rejected_before_work_starts(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    ensure_current_run(create=True)
    with Path(MARKER).open('r') as owner:
        fcntl.flock(owner, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(RuntimeError, match='Another KinBot driver'):
            with kb._driver_lock():
                pass


def test_driver_lock_is_released_on_exit(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with kb._driver_lock():
        pass
    with kb._driver_lock():
        pass
