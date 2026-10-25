"""Refresh code, reconcile accepted VRC models, then run an explicit review plan.

This orchestrates existing owners and live controls. It never selects active
spaces, changes resource policy, cancels allocations or publishes new models.
"""
import argparse
import fcntl
import json
import math
from pathlib import Path
import time

from kinbot import input_refresh, pes_phase, requested_orbital_recon as review
from kinbot import run_policy, run_refresh


def _ready(inputs):
    waiting = []
    for source in inputs:
        root = input_refresh.owner_directory(source)
        if root is None:
            waiting.append(f'{source}: waiting for the PES owner')
            continue
        phase = json.loads((root/pes_phase.NAME).read_text())
        journal = root/input_refresh.JOURNAL
        if journal.exists():
            record = json.loads(journal.read_text())
            if record.get('token') == phase['token']:
                receipt = input_refresh.receipt(root)
                if receipt['status'] == 'blocked':
                    raise RuntimeError('Existing input revision is blocked: '+json.dumps(receipt))
                if receipt['status'] != 'applied':
                    waiting.append(f'{source}: input revision {receipt["status"]}')
                    continue
        if phase['phase'] not in ('exploring', 'postprocessing'):
            waiting.append(f'{source}: owner phase {phase["phase"]}')
            continue
        missing = [name for name, worker in phase['workers'].items()
                   if worker['status'] != 'exploring' and not pes_phase.active_worker(root/name)]
        if missing:
            waiting.append(f'{source}: waiting for workers {", ".join(sorted(missing))}')
    return waiting


def _wait_owners(inputs, seconds):
    deadline = time.monotonic() + seconds
    while waiting := _ready(inputs):
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise RuntimeError('Owner readiness wait expired; rerun this handoff to resume. '
                               + '; '.join(waiting))
        print('WAIT: '+'; '.join(waiting), flush=True)
        time.sleep(min(5., remaining))


def perform(production, plan, output, case_inputs, *, apply=False,
            wait_idle=600., wait_owner=120., expected_max_jobs=None):
    production = Path(production).expanduser().resolve(strict=True)
    plan = Path(plan).expanduser().resolve(strict=True)
    output = Path(output).expanduser().resolve()
    inputs = sorted({Path(p).expanduser().resolve(strict=True) for p in case_inputs})
    if not inputs or any(production not in p.parents for p in inputs):
        raise ValueError('Reconciliation inputs must be explicit case inputs within production.')
    if any(not math.isfinite(v) or v < 0 for v in (wait_idle, wait_owner)):
        raise ValueError('Wait durations must be finite and nonnegative.')
    policy = run_policy.read(production) or {}
    if (expected_max_jobs is not None
            and policy.get('configured_max_jobs', policy.get('max_jobs')) != expected_max_jobs):
        raise ValueError('Configured shared job limit differs from the requested handoff; '
                         'no resource settings were changed.')
    def preflight():
        result = review.batch(production, plan, output, shared_pool=True,
                              preflight=True, resume=True)
        print('REVIEW PREFLIGHT '+json.dumps(result), flush=True)
        return result
    if not apply:
        return dict(status='reviewed', panel=preflight(), case_inputs=list(map(str, inputs)),
                    production_models_accepted=False)
    # Keep duplicate invocations out before either can replace controllers.
    # The campaign also has its own existing submission/manifest lock.
    with (production/'.kinbot_orbital_handoff.lock').open('a+') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError('Another orbital handoff owns this run; inspect its log.') from error
        preflight()  # Every request/deck is checked before controller replacement.
        # Share run_refresh's lock with standalone refresh/salvage invocations.
        with production.with_suffix('.salvage.lock').open('a') as refresh_lock:
            fcntl.flock(refresh_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            refreshed = run_refresh.refresh_when_idle(production, apply=True, wait_idle=wait_idle)
        if refreshed['phase'] != 'launched':
            raise RuntimeError('Controller refresh has not reached launched; review was not started.')
        print('CONTROLLER REFRESH '+json.dumps({k: refreshed.get(k)
              for k in ('phase', 'code', 'archive')}), flush=True)
        print('Retained Slurm IDs: '+', '.join(j['job_id'] for j in refreshed.get('jobs', [])), flush=True)
        _wait_owners(inputs, wait_owner)
        # Preflight all revisions before applying any. submit repeats validation
        # under the owner's locks; a concurrent scientific edit remains a refusal.
        for source in inputs:
            input_refresh.submit(source, reconcile_models=True)
        receipts = {str(source): input_refresh.submit(source, reconcile_models=True, apply=True)
                    for source in inputs}
        print('MODEL RECONCILIATION '+json.dumps(receipts), flush=True)
        prepared = review.batch(production, plan, output, shared_pool=True, resume=True)
        print('REVIEW PANEL '+json.dumps(prepared), flush=True)
        result = review.run(output, report_failures=True)
        return dict(status='review_finished', review=result, reconciliation=receipts,
                    production_models_accepted=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('production', 'plan', 'output'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--case-input', action='append', required=True, type=Path)
    parser.add_argument('--apply', action='store_true')
    parser.add_argument('--wait-idle', type=float, default=600.)
    parser.add_argument('--wait-owner', type=float, default=120.)
    parser.add_argument('--expected-max-jobs', type=int)
    args = parser.parse_args()
    print(json.dumps(perform(args.production, args.plan, args.output, args.case_input,
        apply=args.apply, wait_idle=args.wait_idle, wait_owner=args.wait_owner,
        expected_max_jobs=args.expected_max_jobs), indent=2))


if __name__ == '__main__':
    main()
