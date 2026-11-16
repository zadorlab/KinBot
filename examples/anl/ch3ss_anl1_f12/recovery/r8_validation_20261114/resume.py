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
MODULE = importlib.util.spec_from_file_location('panel_reconciliation', HERE/'reconcile.py')
reconciliation = importlib.util.module_from_spec(MODULE); MODULE.loader.exec_module(reconciliation)


def code_signature():
    from kinbot.anl import path_review, probe_recovery
    return reconciliation.path_review.digest({str(p):file_sha256(p) for p in
        (Path(__file__), HERE/'reconcile.py', Path(path_review.__file__), Path(probe_recovery.__file__))})


def decision_signature(root):
    plan = reconciliation.load(root)
    return plan['sha256'] if plan else None


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
    plan = reconciliation.load(root)
    if plan:
        for proof in reconciliation.proofs(plan): reconciliation.path_review.verify(proof)
    extra = reconciliation.extra_specs(plan)
    specs.update(extra)
    owned, hashes = {}, {}
    for name in specs:
        graph = root/name
        if name in extra and not graph.exists(): continue
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
    return dict(manifests=hashes, retained_jobs=jobs, decision=plan['sha256'] if plan else None,
                scientific_approval=False)


def run(root, interval=20):
    try:
        return _run(root, interval)
    except Exception as error:
        atomic_json(Path(root)/'native_owner.json', dict(pid=os.getpid(), code=file_sha256(__file__),
            status='error', error=str(error), updated_utc=datetime.now(timezone.utc).isoformat()))
        raise


def _run(root, interval=20):
    root = Path(root).resolve()
    with ExitStack() as stack:
        lock = stack.enter_context((root/'controller.lock').open('a'))
        try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError: return
        plan = reconciliation.load(root)
        specs = {**campaign.graphs(root), **reconciliation.extra_specs(plan)}
        locked = set()
        def lock_graphs():
            for name in specs:
                if name not in locked and (root/name).is_dir():
                    handle = stack.enter_context((root/name/'drive.lock').open('a'))
                    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB); locked.add(name)
                    live_control.resume(root/name)
        lock_graphs()
        campaign._stage_unlocked(root)
        lock_graphs()
        reconciliation.apply(root, plan)
        lock_graphs()
        owner = dict(pid=os.getpid(), code=file_sha256(__file__), signature=code_signature(),
                     decision=plan['sha256'] if plan else None)
        def owner_status(status, **fields):
            atomic_json(root/'native_owner.json', dict(**owner, status=status,
                updated_utc=datetime.now(timezone.utc).isoformat(), **fields))
        owner_status('running')
        while True:
            counts, infrastructure, held, terminal = {}, {}, {}, True
            for name in specs:
                try:
                    live_control.process(root/name)
                    allowed = reconciliation.submit_only(plan, name)
                    state = dispatch.advance(root/name, submit=True, submit_only=allowed)
                except dispatch.SchedulerUnavailableError as error:
                    infrastructure[name] = str(error); terminal=False; continue
                except Exception as error:
                    owner_status('error', graph=name, error=str(error))
                    raise
                spec = dispatch._load(root/name)[1]
                current = dispatch._summary(spec, state)
                counts[name] = dict(Counter(current.values()))
                held[name] = [n for n,s in current.items() if allowed is not None and n not in allowed
                              and s not in ('complete','failed','blocked','skipped','unavailable','submitted','submitting')]
                terminal &= all(s in ('complete','failed','blocked','skipped','unavailable')
                                or n in held[name] for n,s in current.items())
            atomic_json(root/'progress.json', dict(workflows=counts, infrastructure=infrastructure,
                held_for_review=held, decision=owner['decision'],
                execution_finished=terminal, scientific_approval=False))
            owner_status('running')
            if terminal:
                atomic_json(root/'completion.json', dict(status='native_attempts_finished_review_required',
                    workflows=counts, held_for_review=held, decision=owner['decision'],
                    scientific_approval=False, mess_updated=False))
                owner_status('finished')
                return
            time.sleep(interval)


def launch(root):
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    log = root/('native.resume.'+stamp+'.log')
    with log.open('ab') as output:
        child = subprocess.Popen([sys.executable,'-B','-u',str(Path(__file__).resolve()),'run',str(root)],
            stdout=output, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True)
    return dict(pid=child.pid, log=str(log), jobs_modified=False)


def recover_setup(root):
    """Repair our journalled setup only when there is no live panel owner.

    Worker allocations remain untouched; manifest journals fence their reads.
    Unowned task directories and ambiguous submissions still fail snapshot().
    """
    plan = reconciliation.load(root)
    if not plan: return
    specs = {**campaign.graphs(root), **reconciliation.extra_specs(plan)}
    with ExitStack() as stack:
        handle = stack.enter_context((root/'controller.lock').open('a'))
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if inventory(root):
            raise run_refresh.RefreshTransientError('A panel owner appeared during setup recovery.')
        for name in specs:
            graph = root/name
            if not graph.is_dir(): continue
            handle = stack.enter_context((graph/'drive.lock').open('a'))
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            live_control.resume(graph)
        for name, spec in reconciliation.extra_specs(plan).items():
            # Only preparation that was already started needs repair here.
            # Missing diagnostics are prepared by the replacement owner.
            if (root/name).exists(): reconciliation.stage_extra(root,name,spec)


def refresh(root):
    root = Path(root).resolve(strict=True)
    records = inventory(root)
    receipt = root/'native_owner.json'
    old = json.loads(receipt.read_text()) if receipt.is_file() else {}
    current_code = old.get('signature') == code_signature()
    # Old no-decision callers retain their idempotent refresh behavior.
    if old.get('signature') is None and decision_signature(root) is None:
        current_code = old.get('code') == file_sha256(__file__)
    if (records and old.get('pid') == records[0]['pid'] and current_code
            and old.get('decision') == decision_signature(root)):
        return dict(status='already_running', pid=old['pid'], jobs_modified=False)
    with run_refresh._submission_guard(root):
        if not records: recover_setup(root)
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
