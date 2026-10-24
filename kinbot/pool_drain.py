"""Conservative pending holds under the shared admission lock.

Running jobs are never stopped. Uncertain scheduler replies stop admission;
only holds recorded here are eligible for automatic release.
"""
import json
from pathlib import Path
import time

from kinbot.execution import atomic_json, SLURM_TERMINAL_STATES


def reconcile(path, policy, jobs, claims):
    from kinbot.anl.live_control import _ownership
    from kinbot.run_refresh import _pending_request
    from kinbot.execution import slurm_job_record
    root = path.parent
    marker = root/'.kinbot_policy_holds.json'
    journal = json.loads(marker.read_text()) if marker.exists() else {'jobs': {}}
    holds = journal['jobs']
    if (sum(jobs.values())+len(claims) <= policy['max_jobs']
            and not any(h['status'] not in ('released','terminal') for h in holds.values())):
        return journal  # The common case needs no per-job scontrol queries.
    records = {}
    for job in set(jobs) | {j for j,h in holds.items() if h['status'] not in ('released','terminal')}:
        current = slurm_job_record(job)
        if not current.get('WorkDir') or not Path(current['WorkDir']).resolve().is_relative_to(root):
            raise RuntimeError(f'{job}: pending-drain ownership is unavailable.')
        current = _ownership(Path(current['WorkDir']), job)
        if not str(current.get('Priority', '')).isdecimal():
            raise RuntimeError(f'{job}: pending-drain priority is unavailable.')
        records[job] = current
    # Running/configuring/completing jobs and ambiguous submissions reserve
    # capacity. Pending jobs count too; no speculative dependency discount.
    running = sum(jobs[j] for j,r in records.items() if j in jobs and r['JobState'] != 'PENDING')
    spare = max(0, policy['max_jobs']-running-len(claims))
    for job,r in sorted(records.items(), key=lambda pair:int(pair[0])):
        phase = r['JobState'].split('+')[0]
        h = holds.get(job)
        ours = h is not None and h['status'] not in ('released','terminal')
        if ours and Path(h['directory']).resolve() != Path(r['WorkDir']).resolve():
            raise RuntimeError(f'{job}: held allocation changed ownership.')
        if phase in SLURM_TERMINAL_STATES:
            if ours: h['status'] = 'terminal'
            continue
        if phase != 'PENDING':
            # A hold raced with start: retain the calculation, undo only ours.
            want_hold = False
        elif r['Priority'] == '0' and not ours:
            continue  # Never release an administrator/user's pre-existing hold.
        else:
            cost = jobs.get(job, 1)
            want_hold = spare < cost
            if not want_hold: spare -= cost
        if want_hold:
            if not ours:
                h = holds[job] = dict(directory=r['WorkDir'], status='holding', original=r)
                atomic_json(marker, journal)
            if r['Priority'] != '0':
                _pending_request(['scontrol', 'hold', job])
                r = _ownership(Path(h['directory']), job)
            if r['JobState'] == 'PENDING' and r.get('Priority') == '0':
                h['status'] = 'held'
            else:
                # Re-evaluate the running count on the next poll, not on an
                # optimistic snapshot. No additional job may be admitted.
                atomic_json(marker, journal)
                raise RuntimeError(f'{job}: pending hold raced with start or was not confirmed.')
        elif ours:
            h['status'] = 'releasing'; atomic_json(marker, journal)
            if r['Priority'] == '0':
                _pending_request(['scontrol', 'release', job])
                r = _ownership(Path(h['directory']), job)
                if r.get('Priority') == '0':
                    raise RuntimeError(f'{job}: pending hold release is unconfirmed.')
            h['status'] = 'released'
    journal.update(updated_at=time.time(), status='reconciled', error=None)
    atomic_json(marker, journal)
    return journal


def pulse(directory, *, force=False):
    """One pool-wide scheduler refresh per 20 seconds, independent of sbatch."""
    from kinbot import run_policy
    path = run_policy.find(directory)
    if path is None or not run_policy._read(path).get('manage_pending'):
        return None
    with run_policy._locked(path):
        policy = run_policy._read(path)
        stamp = path.parent/'.kinbot_policy_drain_poll.json'
        previous = json.loads(stamp.read_text()) if stamp.exists() else {}
        if not force and (previous.get('revision') == policy['revision']
                          and time.time()-previous.get('time',0) < 20):
            return previous
        result = dict(time=time.time(),revision=policy['revision'])
        try:
            policy = run_policy.effective(path, policy)
            ledger = path.parent/'.kinbot_policy_admissions.json'
            claims = json.loads(ledger.read_text()) if ledger.exists() else {}
            jobs = run_policy._queue(path.parent)
            run_policy._reconcile(path.parent,jobs,claims,archive=True)
            atomic_json(ledger,claims)
            result['holds'] = reconcile(path,policy,jobs,claims)
            result['status'] = 'reconciled'
        except (OSError,RuntimeError) as error:
            result.update(status='blocked',error=str(error))
        atomic_json(stamp,result)
        return result
