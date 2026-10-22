"""In-process VRC surface controllers; QC stays in the shared Slurm pool.

No working-directory changes or child controller processes are used. Workers
join before the PES owner can apply a live revision. A failed surface does not
stop independent surfaces; an epoch change quiesces every surface.
"""
from contextlib import contextmanager
from queue import Queue, Empty
import threading
from pathlib import Path


_LOCAL = threading.local()
_LOCKS = {}
_LOCKS_GUARD = threading.Lock()


class SurfaceInterrupted(RuntimeError):
    pass


def checkpoint():
    stop = getattr(_LOCAL, 'stop', None)
    if stop is not None and stop.is_set():
        raise SurfaceInterrupted('VRC surface controllers are quiescing; QC jobs retained.')


@contextmanager
def guard(key):
    """Serialize shared files/graphs among surfaces of the same KB owner."""
    key = str(Path(key).resolve())
    with _LOCKS_GUARD:
        lock = _LOCKS.setdefault(key, threading.RLock())
    from kinbot.pes_phase import checkpoint as check
    check()
    while not lock.acquire(timeout=.2):
        check()
    try:
        check()
        yield
    finally:
        lock.release()


def run(jobs, *, limit, completed, started=lambda name: None):
    """Run named callables, reporting results on the caller thread only.

    ``completed(name, value, error)`` sees each result once. Scientific/native
    failures are returned to that callback; interruption is raised only after
    all local controllers have released their graph ownership.
    """
    from kinbot.pes_phase import checkpoint as check, ExplorationReopened
    if type(limit) is not int or limit < 1:
        raise ValueError('Parallel VRC controller limit must be a positive integer.')
    jobs = list(jobs)
    if not jobs:
        return
    if limit == 1:
        for name, function in jobs:
            check()
            started(name)
            try:
                value = function()
            except ExplorationReopened:
                raise
            except Exception as error:
                completed(name, None, error)
            else:
                completed(name, value, None)
        return
    pending, results = Queue(), Queue()
    for job in jobs:
        pending.put(job)
    stop = threading.Event()

    def worker():
        _LOCAL.stop = stop
        try:
            while True:
                try:
                    name, function = pending.get_nowait()
                except Empty:
                    return
                try:
                    check()
                    results.put(('started', name, None, None))
                    value = function()
                except BaseException as error:
                    results.put(('completed', name, None, error))
                else:
                    results.put(('completed', name, value, None))
        finally:
            del _LOCAL.stop

    # SIGTERM during run_refresh exits the KB process without waiting for
    # Python executor shutdown. Submitted jobs keep their existing ledger.
    workers = []
    interrupt = None
    for i in range(min(limit, len(jobs))):
        thread = threading.Thread(target=worker, name=f'kinbot-vrc-{i}', daemon=True)
        try:
            thread.start()
        except BaseException as error:
            stop.set()
            if not workers:
                raise
            interrupt = error
            break
        workers.append(thread)
    remaining = len(jobs)
    while remaining:
        try:
            if interrupt is None:
                check()
            kind, name, value, error = results.get(timeout=.2)
        except Empty:
            continue
        except BaseException as error:
            interrupt = error
            stop.set()
            continue
        if kind == 'completed':
            remaining -= 1
        if isinstance(error, ExplorationReopened) or (error is not None and not isinstance(error, Exception)):
            interrupt = interrupt or error
            stop.set()
        try:
            if kind == 'started':
                started(name)
            else:
                completed(name, value, error)
        except BaseException as error:
            interrupt = interrupt or error
            stop.set()
    for thread in workers:
        thread.join()
    if interrupt is not None:
        raise interrupt
