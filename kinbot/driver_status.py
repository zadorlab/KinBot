"""Explicit controller outcomes: a process exit alone is not success."""
from datetime import datetime, timezone
import os
from pathlib import Path
import socket

from kinbot.execution import atomic_json


def run_recorded(function, filename):
    path = Path(filename)
    record = {'schema': 1, 'status': 'running', 'pid': os.getpid(),
              'host': socket.gethostname(),
              'started_utc': datetime.now(timezone.utc).isoformat()}
    atomic_json(path, record)
    try:
        result = function()
    except BaseException as error:
        code = (error.code or 0) if isinstance(error, SystemExit) else 1
        record.update(status='complete' if code in (None, 0) else 'failed',
                      returncode=code, error=f'{type(error).__name__}: {error}')
        raise
    else:
        record.update(status='complete' if result in (None, 0) else 'failed',
                      returncode=result or 0)
        if record['status'] == 'failed':
            record['error'] = f'Controller returned {result}; inspect its kinbot.log.'
        return result
    finally:
        record['finished_utc'] = datetime.now(timezone.utc).isoformat()
        atomic_json(path, record)
