"""Queue this finite diagnostic review concurrently, retaining native jobs.

Only the review owners are replaced. Production owners, inputs and allocations
are never stopped. All diagnostic seeds are validated before the handoff; no
model is published and orbital-predecessor dependencies remain enforced.
"""
import argparse
from contextlib import ExitStack
import fcntl
import json
import os
from pathlib import Path
import shlex
import signal
import subprocess
import sys
import tempfile
import time
import uuid

import continue_ch2_sh_s_review as followup
from kinbot import requested_orbital_recon as screens
from kinbot import vrc_relaxation_recon as relaxation
from kinbot import run_cleanup, run_refresh, run_salvage, run_policy
from kinbot.anl import dispatch
from kinbot.execution import atomic_json, file_sha256, slurm_job_record
from kinbot.rotdpy_worker_code import validate


TOOLS = Path(__file__).resolve().parent
PLAN = TOOLS.parent/'orbital_review_batches/r7_progress_20261024'


def read(path):
    return json.loads(Path(path).read_text())


def campaigns(review):
    return [('screens', 'requested_orbital_recon.py', review/'screens'),
            ('relaxations', 'vrc_relaxation_recon.py', review/'relaxations'),
            ('compact', 'requested_orbital_recon.py', review/'ch2_sh_s_followup/candidates'),
            ('compact_relaxations', 'vrc_relaxation_recon.py', review/'ch2_sh_s_followup/relaxations')]


def owner(record, production, review):
    """An exact known script/argument pair, never a PID or substring alone."""
    args = shlex.split(record['command'])
    if '-m' in args:
        index = args.index('-m')
        module = args[index+1] if len(args) > index+1 else ''
        from kinbot.review_campaign import MODULES
        if module in MODULES.values():
            script = module.rsplit('.', 1)[1]+'.py'
            return owner(dict(record, command=shlex.join(
                [*args[:index], str(TOOLS/script), *args[index+2:]])), production, review)
    scripts = {name for _, name, _ in campaigns(review)} | {
        'refresh_progress_20261024.sh', 'continue_ch2_sh_s_review.py',
        'queue_progress_review.py'}
    for i, arg in enumerate(args):
        name = Path(arg).name
        if name not in scripts:
            continue
        if not ((Path(arg).is_absolute() and Path(arg).resolve()==TOOLS/name)
                or arg=='examples/anl/tools/'+name):
            continue
        tail = args[i+1:]
        if name == 'queue_progress_review.py':
            if (not tail or not Path(tail[0]).is_absolute()
                    or Path(tail[0]).resolve() != production):
                return False
            remaining = tail[1:]
            while remaining:
                option = remaining.pop(0)
                if option == '--replace-owner':
                    continue
                if option == '--allow-node' and remaining:
                    run_policy._nodes([remaining.pop(0)])
                    continue
                if option == '--code' and remaining:
                    pin = Path(remaining.pop(0))
                    if (pin.is_absolute() and pin.parent.resolve() == production
                            and pin.name.startswith('.salvage_code.')):
                        continue
                return False
            return True
        if name in ('refresh_progress_20261024.sh', 'continue_ch2_sh_s_review.py'):
            return len(tail)==1 and Path(tail[0]).is_absolute() and Path(tail[0]).resolve()==production
        targets = {p for _, script, p in campaigns(review) if script == name}
        if (len(tail)<2 or tail[0]!='run' or not Path(tail[1]).is_absolute()
                or Path(tail[1]).resolve() not in targets):
            return False
        remaining = tail[2:]
        while remaining:
            option = remaining.pop(0)
            if option in ('--queue-all', '--report-failures'):
                continue
            if option == '--allow-node' and remaining:
                run_policy._nodes([remaining.pop(0)])
                continue
            return False
        return True
    return False


def inventory(production, review):
    records = run_cleanup._processes(Path('/'))
    selected = [p for p in records if p['pid'] != os.getpid()
                and owner(p, production, review)]
    ids = {p['pid'] for p in selected}
    # Local scheduler/preflight children must finish their bookkeeping. They
    # are not killed together with an owner, even if they run outside review.
    for p in records:
        if p['ppid'] in ids and p['pid'] not in ids:
            raise run_refresh.RefreshBusyError(p)
    return selected


def submission_check(review):
    for _, script, folder in campaigns(review):
        marker = folder/('review_launch.json' if script.startswith('requested') else 'panel.json')
        if not marker.exists():
            continue
        record = read(marker)
        if script.startswith('requested'):
            rows = [j for s in record['systems'] for j in read(s['manifest'])['jobs']]
        else:
            from kinbot.anl.submission_recovery import reconcile_pre_admission
            # Called with admission guards held and owners frozen. Repair
            # only a proven legacy write-before-admission, never a live RPC.
            rows = [t for g in record['graphs']
                    for t in reconcile_pre_admission(g, pool_locked=True)['tasks'].values()]
        if any(j.get('submission_status') == 'submitting' or j.get('status') == 'submitting' for j in rows):
            raise run_refresh.RefreshTransientError('A review submission is being journaled; retrying the idle boundary.')


def stop_owners(production, review, marker, *, timeout=180, validate_saved=None):
    run_refresh._probe_process_handles()
    deadline = time.monotonic()+timeout
    while True:
        try:
            records = inventory(production, review)
            with ExitStack() as guards:
                for _, script, folder in campaigns(review):
                    if script.startswith('vrc') and folder.exists():
                        guards.enter_context(run_refresh._submission_guard(folder, independent=True))
                        # Never freeze an owner halfway through a panel or
                        # workflow/state transaction, then wait for its lock.
                        paths=[folder/'.prepare.lock']
                        if (folder/'panel.json').exists():
                            paths += [Path(g)/'.manifest.lock' for g in read(folder/'panel.json')['graphs']]
                        for path in paths:
                            if path.exists():
                                stream=guards.enter_context(path.open('r'))
                                try:
                                    fcntl.flock(stream,fcntl.LOCK_SH|fcntl.LOCK_NB)
                                except BlockingIOError:
                                    raise run_refresh.RefreshTransientError('Review panel is being staged: '+str(path)) from None
                with run_refresh._freeze(review, records,
                        inventory=lambda _: inventory(production, review)) as handles:
                    submission_check(review)
                    if validate_saved is not None:
                        validate_saved()
                    jobs = run_cleanup._jobs(review)
                    atomic_json(marker, dict(phase='stopping_review_owners', owners=records,
                                             retained_jobs=jobs, production=str(production)))
                    for _, fd in handles:
                        run_refresh._pidfd_send_signal(fd, signal.SIGKILL)
            print('Review controllers replaced; retained Slurm IDs: '+
                  ', '.join(j['job_id'] for j in jobs), flush=True)
            return jobs
        except run_refresh.RefreshTransientError as error:
            if time.monotonic() >= deadline:
                raise RuntimeError('Review owner stayed busy; existing work was preserved. '+str(error)) from error
            print('Waiting for review submission bookkeeping: '+str(error), flush=True)
            time.sleep(2)


def acquire(stack, path, timeout=30):
    stream = stack.enter_context(path.open('a+'))
    deadline = time.monotonic()+timeout
    while True:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
            return
        except BlockingIOError:
            if time.monotonic() >= deadline:
                raise RuntimeError(f'Owner lock still held: {path}; do not delete it.')
            time.sleep(.1)


def preflight(production, review):
    validate_existing(production, review)
    screens._saved(review/'screens')
    # Temporary decks are never submitted. Preparation verifies the completed
    # native seed hashes, request identity, constraints and both candidate CASs.
    with tempfile.TemporaryDirectory(prefix='.queue_preflight.', dir=review) as temp:
        staged = followup.stage(production, review, Path(temp))
        original = relaxation.prepare(production, PLAN/'relaxations.json', Path(temp)/'original')
        for graph in [*staged['graphs'], *original['graphs']]:
            dispatch.preflight(graph)
    run_cleanup._jobs(review)  # An unavailable scheduler cannot mean no jobs.


def validate_existing(production, review):
    """Read the actual destinations, not just freshly reconstructed decks."""
    relaxation.prepare(production, PLAN/'relaxations.json', review/'relaxations', validate_only=True)
    followup.stage(production, review, review/'ch2_sh_s_followup', validate_only=True)


def expand_pending(review, allowed, journal):
    """Change only this review's still-pending jobs; retain every other node."""
    for job in run_cleanup._jobs(review):
        if not Path(job['directory']).resolve().is_relative_to(review):
            raise ValueError('Pending update escaped the review directory.')
        if job['state'] != 'PENDING':
            continue
        ident = job['job_id']; before = slurm_job_record(ident)
        if before.get('JobState') != 'PENDING':
            continue
        if (before.get('JobId') != ident or not before.get('WorkDir')
                or Path(before['WorkDir']).resolve() != Path(job['directory']).resolve()):
            raise RuntimeError('Review pending-job identity changed: '+ident)
        old = before.get('ExcNodeList')
        if old in ('', '(null)', 'None'):
            continue
        if old is None:
            raise RuntimeError('Scheduler omitted exclusions for '+ident)
        expanded = subprocess.run(['scontrol','show','hostnames',old],capture_output=True,
                                  text=True,check=True,timeout=20)
        nodes = set(expanded.stdout.split()); desired = nodes-set(allowed)
        if desired == nodes:
            continue
        journal.setdefault('pending_changes', {})[ident] = dict(before=before, excluded=sorted(desired))
        atomic_json(review/'queue_all_refresh.json', journal)
        result = subprocess.run(['scontrol','update','JobId='+ident,'ExcNodeList='+','.join(sorted(desired))],
                                capture_output=True,text=True,timeout=20)
        after = slurm_job_record(ident)
        if result.returncode:
            if after.get('JobState') in ('RUNNING','CONFIGURING','COMPLETING','COMPLETED'):
                continue
            raise RuntimeError('Pending exclusion update failed: '+result.stderr.strip())
        if after.get('JobState')=='PENDING':
            if after.get('JobId')!=ident or after.get('WorkDir')!=before['WorkDir']:
                raise RuntimeError('Pending identity changed after update: '+ident)
            actual=after.get('ExcNodeList')
            if actual is None:
                raise RuntimeError('Pending exclusions missing after update: '+ident)
            nodes=set() if actual in ('','(null)','None') else set(subprocess.run(
                ['scontrol','show','hostnames',actual],capture_output=True,text=True,check=True,timeout=20).stdout.split())
            if nodes!=desired:
                raise RuntimeError('Pending exclusion update was not confirmed: '+ident)
        journal['pending_changes'][ident]['after']=after
        atomic_json(review/'queue_all_refresh.json',journal)


def supervise(review, pin, allowed):
    from kinbot.review_campaign import supervise as run_campaigns
    records = [dict(name=name, kind='orbitals' if script.startswith('requested') else 'relaxations',
                    directory=str(folder)) for name, script, folder in campaigns(review)]
    return run_campaigns(review, records, code=pin, allow_nodes=allowed, queue_all=True)


def perform(production, review, pin, allowed, *, replace_owner=False):
    allowed = run_policy._nodes(allowed)
    destination = review/'ch2_sh_s_followup'; destination.mkdir(exist_ok=True)
    with ExitStack() as locks, ExitStack() as transition:
        # Serialize replacement handoffs separately from the supervisor lock.
        # Older supervisors hold only .queue_all.lock. A normal duplicate
        # invocation must still fail without stopping anything.
        acquire(transition, review/'.queue_transition.lock', timeout=0)
        if not replace_owner:
            acquire(locks, review/'.queue_all.lock', timeout=0)
        print('Checking the finite review plans and native seeds before replacing review owners.',flush=True)
        preflight(production, review)
        marker = review/'queue_all_refresh.json'
        if marker.exists():
            atomic_json(review/'.queue_all_history'/(uuid.uuid4().hex+'.json'),read(marker))
        jobs = stop_owners(production, review, marker,
                           validate_saved=lambda: validate_existing(production, review))
        if replace_owner:
            acquire(locks, review/'.queue_all.lock')
        acquire(locks, review/'.handoff.lock')
        acquire(locks, destination/'.owner.lock')
        # Preparation is idempotent and retains saved graph/task/job identities.
        relaxation.prepare(production, PLAN/'relaxations.json', review/'relaxations')
        followup.stage(production, review, destination)
        journal = read(marker)
        journal.update(phase='prepared',code=str(pin),allowed_nodes=allowed)
        expand_pending(review, allowed, journal)
        journal.update(phase='owners_starting',retained_jobs=jobs)
        atomic_json(marker,journal)
        transition.close()
        return supervise(review, pin, allowed)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('production')
    parser.add_argument('--allow-node',action='append',default=[])
    parser.add_argument('--replace-owner',action='store_true',
                        help='Replace this review supervisor and owners after preflight; retain Slurm jobs.')
    parser.add_argument('--code',help=argparse.SUPPRESS)
    args = parser.parse_args()
    production = Path(args.production).expanduser().resolve(strict=True)
    review = production.parent/'r7_progress_review_20261024'
    read(review/'refresh_complete.json')
    if not args.code:
        with ExitStack() as locks:
            acquire(locks,review/'.queue_pin.lock')
            pin = run_salvage._pin_code(production,refresh=True)
        os.execve(sys.executable,[sys.executable,'-B','-u',str(Path(__file__).resolve()),
            *sys.argv[1:],'--code',str(pin)],dict(os.environ,PYTHONPATH=str(pin),
            OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1'))
    pin = validate(args.code)
    import kinbot
    if Path(kinbot.__file__).resolve().parent.parent != pin:
        raise RuntimeError('Review controller must load the verified code snapshot.')
    result=perform(production,review,pin,args.allow_node,replace_owner=args.replace_owner)
    if any(result.values()):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
