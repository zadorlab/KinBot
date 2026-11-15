"""Resume this frozen panel with live controls, preserving every QC allocation.

The scientific generator is intentionally unchanged: existing stage and request
fingerprints still apply. Refresh replaces only its verified local Python owner.
"""
import argparse
from collections import Counter
from contextlib import ExitStack
from datetime import datetime, timezone
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from kinbot import run_cleanup, run_refresh
from kinbot.anl import dispatch, live_control
from kinbot.execution import atomic_json, file_sha256

HERE = Path(__file__).resolve().parent
MODULE = importlib.util.spec_from_file_location('frozen_panel', HERE/'campaign.py')
campaign = importlib.util.module_from_spec(MODULE); MODULE.loader.exec_module(campaign)


def inventory(root):
    known = []
    for p in Path('/proc').iterdir():
        if not p.name.isdecimal(): continue
        try:
            if p.stat().st_uid != os.getuid(): continue
            args = [s.decode() for s in (p/'cmdline').read_bytes().split(b'\0') if s]
            if not args or not Path(args[0]).name.startswith('python'): continue
            args = [args[0], *[a for a in args[1:] if a not in ('-B', '-u')]]
            if len(args) != 4 or args[2:] != ['run', str(root)]: continue
            script = Path(args[1]); script = script if script.is_absolute() else (p/'cwd').resolve()/script
            if script.resolve() not in (HERE/'campaign.py', Path(__file__).resolve()): continue
            stat = (p/'stat').read_text().rsplit(')',1)[1].split()
            if stat[0] == 'Z': continue
            known.append(dict(pid=int(p.name), starttime=stat[19]))
        except (OSError, ValueError, IndexError, UnicodeError): continue
    records = run_cleanup._processes(root, known=known)
    allowed = {(r['pid'],r['starttime']) for r in known}
    for r in records:
        if (r['pid'],r['starttime']) not in allowed:
            raise run_refresh.RefreshBusyError(r)
    if len(records) > 1:
        raise run_refresh.RefreshTransientError('More than one panel launcher is present; waiting for lock ownership.')
    return records


def snapshot(root):
    stage = json.loads((root/'stage.json').read_text())
    if stage['sources'] != campaign.fingerprint():
        raise ValueError('Frozen panel source fingerprints changed.')
    specs = campaign.graphs(root)
    if list(specs) != stage['workflows']: raise ValueError('Panel graph inventory changed.')
    owned, hashes = {}, {}
    for name in stage['workflows']:
        graph = root/name
        if not all((graph/f).is_file() for f in ('state.json','workflow.json','requested.json')):
            raise run_refresh.RefreshTransientError('Panel preparation is in flight.')
        transaction = graph/'.live_control_transaction.json'
        if transaction.is_file() and json.loads(transaction.read_text()).get('status') == 'applying':
            raise run_refresh.RefreshTransientError('Attempt archival is in flight.')
        _, spec, state = dispatch._load(graph)
        if json.loads((graph/'requested.json').read_text()) != specs[name]:
            raise ValueError('Panel scientific definition changed.')
        for directory in (graph/'tasks').iterdir():
            if directory.is_dir() and directory.name not in state['tasks']:
                raise run_refresh.RefreshTransientError('Task staging is in flight.')
        if any(e.get('status') == 'submitting' for e in state['tasks'].values()):
            raise run_refresh.RefreshTransientError('Submission bookkeeping is in flight.')
        for task in spec['tasks']:
            entry = state['tasks'].get(task['id'], {})
            if entry.get('task_sha256'):
                dispatch._verify_stage_files(graph, task, entry)
            if entry.get('status') == 'submitted':
                job = str(entry['job_id'])
                if job in owned: raise ValueError('Duplicate allocation ownership.')
                owned[job] = str(graph/'tasks'/task['id'])
        for filename in ('workflow.json', 'state.json', 'requested.json'):
            hashes[name+'/'+filename] = file_sha256(graph/filename)
    jobs = run_cleanup._jobs(root)
    if any(owned.get(j['job_id']) != str(Path(j['directory']).resolve()) for j in jobs):
        raise ValueError('A live allocation lacks matching panel ownership.')
    return dict(manifests=hashes, retained_jobs=jobs, scientific_approval=False)


def run(root, interval=20):
    root = Path(root).resolve()
    with ExitStack() as stack:
        lock = stack.enter_context((root/'controller.lock').open('a'))
        try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError: return
        specs = campaign._stage_unlocked(root)
        for name in specs:
            lock = stack.enter_context((root/name/'drive.lock').open('a'))
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        atomic_json(root/'native_owner.json', dict(pid=os.getpid(), code=file_sha256(__file__), status='running'))
        while True:
            counts, infrastructure, terminal = {}, {}, True
            for name in specs:
                try:
                    live_control.process(root/name)
                    state = dispatch.advance(root/name, submit=True)
                except dispatch.SchedulerUnavailableError as error:
                    infrastructure[name] = str(error); terminal=False; continue
                spec = dispatch._load(root/name)[1]
                current = dispatch._summary(spec, state)
                counts[name] = dict(Counter(current.values()))
                terminal &= all(s in ('complete','failed','blocked','skipped','unavailable') for s in current.values())
            atomic_json(root/'progress.json', dict(workflows=counts, infrastructure=infrastructure,
                execution_finished=terminal, scientific_approval=False))
            if terminal:
                atomic_json(root/'completion.json', dict(status='native_attempts_finished_review_required',
                    workflows=counts, scientific_approval=False, mess_updated=False))
                atomic_json(root/'native_owner.json', dict(pid=os.getpid(), code=file_sha256(__file__), status='finished'))
                return
            time.sleep(interval)


def launch(root):
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    log = root/('native.resume.'+stamp+'.log')
    with log.open('ab') as output:
        child = subprocess.Popen([sys.executable,'-B','-u',str(Path(__file__).resolve()),'run',str(root)],
            stdout=output, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True)
    return dict(pid=child.pid, log=str(log), jobs_modified=False)


def refresh(root):
    root = Path(root).resolve(strict=True)
    records = inventory(root)
    receipt = root/'native_owner.json'
    old = json.loads(receipt.read_text()) if receipt.is_file() else {}
    if records and old.get('pid') == records[0]['pid'] and old.get('code') == file_sha256(__file__):
        return dict(status='already_running', pid=old['pid'], jobs_modified=False)
    with run_refresh._submission_guard(root):
        with run_refresh._freeze(root, records, inventory=inventory) as handles:
            proof = snapshot(root)
            stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
            atomic_json(root/'controller_handoffs'/(stamp+'.json'), dict(previous=records, **proof))
            for _, fd in handles: run_refresh._pidfd_send_signal(fd, signal.SIGKILL)
        # Acquiring the original lock inode proves the retired owner released
        # it. A foreign owner or incomplete handoff leaves all QC work intact.
        with (root/'controller.lock').open('a') as lock:
            deadline = time.monotonic()+5
            while True:
                try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB); break
                except BlockingIOError:
                    if time.monotonic() >= deadline:
                        raise run_refresh.RefreshTransientError('Panel owner has not released its lock.')
                    time.sleep(.05)
            # Launch only after state was proven safe, never scancel anything.
        return launch(root)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('refresh','run'))
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    if args.action == 'run': run(args.root); return
    deadline = time.monotonic()+120
    while True:
        try: result = refresh(args.root); break
        except run_refresh.RefreshTransientError as error:
            if time.monotonic() >= deadline: raise
            print(str(error), flush=True); time.sleep(2)
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
