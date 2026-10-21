"""PES-wide exploration barrier and cooperative reopening of live searches.

Parked well controllers retain their in-memory accepted chemistry. They do not
occupy exploration slots or submit downstream work until the PES owner opens
the barrier. QC allocations are never signalled by this coordinator.
"""
from contextlib import contextmanager
import fcntl
import json
import os
from pathlib import Path
import time
import uuid
import socket
import subprocess

from kinbot.execution import atomic_json

NAME = 'pes_exploration.json'
_EPOCH = None
_ROOT = None
_WORKER = None


class ExplorationReopened(RuntimeError):
    pass


@contextmanager
def locked(root):
    with (Path(root)/'.pes_exploration.lock').open('a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        yield


def process_identity(pid):
    try:
        stat = Path(f'/proc/{pid}/stat')
        if Path('/proc').is_dir():
            fields=stat.read_text().rsplit(')',1)[1].split()
            if fields[0]=='Z':
                return None
            start=fields[19]
        else:
            start=subprocess.check_output(['ps','-p',str(pid),'-o','lstart='],text=True,timeout=5).strip()
        boot=Path('/proc/sys/kernel/random/boot_id')
        return {'host':socket.gethostname(),'pid':pid,'start':start,
                'boot':boot.read_text().strip() if boot.exists() else None} if start else None
    except (OSError,subprocess.SubprocessError,IndexError):
        return None


def owner_alive(record):
    identity=record.get('owner_identity')
    return bool(identity and process_identity(record['owner_pid'])==identity)


def initialize(root, input_file=None, initial_worker=None):
    root = Path(root).resolve()
    record = dict(schema=1, token=uuid.uuid4().hex, epoch=0,
                  phase='exploring', owner_pid=os.getpid(),
                  owner_identity=process_identity(os.getpid()), workers={})
    if initial_worker is not None:
        record['initial_worker']=initial_worker
    if input_file is not None:
        from kinbot.input_refresh import capture
        record.update(capture(root,input_file))
    atomic_json(root/NAME, record)
    return record


def register(root, worker, pid):
    with locked(root):
        m = json.loads((Path(root)/NAME).read_text())
        m['workers'][worker] = dict(pid=pid, identity=process_identity(pid), status='exploring')
        atomic_json(Path(root)/NAME, m)


def attach():
    global _ROOT, _WORKER, _EPOCH
    root = os.environ.get('KINBOT_PES_PHASE_ROOT')
    token = os.environ.get('KINBOT_PES_PHASE_TOKEN')
    if not root or not token or Path.cwd().parent != Path(root):
        return False
    m = json.loads((Path(root)/NAME).read_text())
    if m['token'] != token:
        raise RuntimeError('PES exploration owner changed; refusing a stale controller.')
    _ROOT, _WORKER, _EPOCH = Path(root), Path.cwd().name, None
    return True


def _read():
    if _ROOT is None:
        return None
    m = json.loads((_ROOT/NAME).read_text())
    if m['token'] != os.environ.get('KINBOT_PES_PHASE_TOKEN'):
        raise RuntimeError('PES exploration generation changed underneath a worker.')
    return m


def active_worker(directory):
    """A closed search can reopen only under a live, revision-aware PES owner."""
    directory = Path(directory).resolve()
    path = directory.parent/NAME
    if not path.exists():
        return False
    m = json.loads(path.read_text())
    if m['phase'] not in ('exploring', 'postprocessing') or directory.name not in m['workers']:
        return False
    worker=m['workers'][directory.name]
    return (owner_alive(m) and bool(worker.get('identity'))
            and process_identity(worker['pid'])==worker['identity'])


def checkpoint():
    """Called at safe owner polling boundaries, outside native QC processes."""
    if _EPOCH is None:
        return
    m = _read()
    if m['epoch'] != _EPOCH or m['phase'] == 'exploring':
        raise ExplorationReopened('PES exploration reopened; retaining submitted QC allocations.')


def _pending(root):
    from kinbot import live_control as control
    return [(p, r) for p, r in control.requests(root, kinds={'reaction','network'})
            if (Path(r['target']).parent if r['kind']=='network' else Path(r['target']).parent.parent) == Path(root)]


def tick(root, jobs, *, max_searches, known_jobs=None):
    """PES owner only. Return workers occupying exploration slots and phase."""
    root = Path(root).resolve()
    with locked(root):
        m = json.loads((root/NAME).read_text())
        from kinbot.input_refresh import consume
        consume(root,m)
        pending = _pending(root)
        targets = {(Path(r['target']).name if r['kind']=='network' else Path(r['target']).parent.name) for _, r in pending}
        workers = m['workers']
        # A worker can publish another well between the owner's directory
        # read and this tick. Its ready acknowledgement uses the same lock.
        discovered = (root/'chemids').read_text().splitlines() if known_jobs is not None else []
        current_inventory = known_jobs is None or set(discovered) <= set(known_jobs)
        if m['phase'] == 'postprocessing' and targets:
            m.update(phase='exploring', epoch=m['epoch']+1, grants=[])
        busy = [n for n, w in workers.items() if w['status'] == 'exploring']
        # All branch threads must first surrender their graph locks. Waiting
        # for this acknowledgement prevents editing a plan while it is used.
        if m['phase'] == 'exploring' and not any(w['status']=='postprocessing' for w in workers.values()):
            grants = sorted(n for n in targets if n in workers
                            and workers[n]['status'] in ('ready','assembled'))[:max(0,max_searches-len(busy))]
            m['grants'] = grants
            busy += grants
            if (current_inventory and not targets and set(jobs) <= set(workers)
                    and all(workers[n]['status'] in ('ready','assembled') for n in jobs)):
                m.update(phase='postprocessing', grants=[])
        if (current_inventory and m['phase']=='postprocessing' and not targets and set(jobs)<=set(workers)
                and all(workers[n]['status']=='assembled' for n in jobs)):
            m['phase'] = 'closing'
        atomic_json(root/NAME, m)
    return busy, m['phase']


def _status(status):
    with locked(_ROOT):
        m = _read()
        record = m['workers'].get(_WORKER)
        if record is None:  # Parent has not published Popen's PID yet.
            return False
        if record['pid'] != os.getpid():
            raise RuntimeError('PES worker ownership changed.')
        record.update(status=status, epoch=m['epoch'])
        atomic_json(_ROOT/NAME, m)
        return True


def wait_for_exploration(generator=None):
    global _EPOCH
    from kinbot.pes_snapshot import publish
    publish(generator)
    if _ROOT is None:
        return
    _EPOCH = None
    while not _status('ready'):
        time.sleep(.1)
    while True:
        m = _read()
        if m['phase']=='postprocessing':
            _EPOCH = m['epoch']
            _status('postprocessing')
            checkpoint()
            return
        if _WORKER in m.get('grants', ()) and any(
                (Path(r['target']).name if r['kind']=='network' else Path(r['target']).parent.name)==_WORKER
                for _,r in _pending(_ROOT)):
            if generator is None:
                raise RuntimeError('Cannot revisit a worker without a reaction search.')
            _status('exploring')
            authorize_revision(Path.cwd(), m['token']+'_'+str(m['epoch']))
            from kinbot.search_control import reactions, revise_network
            revise_network(generator)
            reactions(generator)
            generator.generate()
            publish(generator)
            _status('ready')
        time.sleep(1)


def assembled():
    """Retain the worker until the whole PES has joined, permitting live edits."""
    if _ROOT is None:
        return
    _status('assembled')
    while True:
        checkpoint()
        if _read()['phase']=='closing':
            return
        time.sleep(1)


def authorize_revision(worker, revision):
    """Grant one append-only publication; called only after branches quiesce."""
    from kinbot.anl.ts_reuse import digest
    worker = Path(worker)
    # The input determines this directory, including non-default locations.
    inp = worker/(worker.name+'.json')
    parameters = json.loads(inp.read_text()) if inp.exists() else {}
    base = worker/parameters.get('anl_automatic_directory', 'anl_automatic')
    plan = base/'network_request.json'
    if plan.exists():
        from kinbot.anl.network_update import archive_assembly
        archive = worker/'.kinbot_network_revisions'/revision
        archive_assembly(base, archive)
        atomic_json(base/'network_revisit.json', dict(schema=1,
            previous_sha256=digest(json.loads(plan.read_text())),
            archive=str(archive), status='authorized'))
