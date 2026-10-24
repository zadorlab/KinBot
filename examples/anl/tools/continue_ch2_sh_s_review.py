"""Continue the finite CH2(SH)S review after the existing diagnostic owner.

This example never refreshes production controllers or publishes surface models.
The old handoff lock serializes access to its three reserved review slots. Its
campaigns are reconciled to termination before the new comparisons are run.
"""
import argparse
from contextlib import ExitStack
import fcntl
import json
import os
from pathlib import Path
import sys

import requested_orbital_recon as screens
import vrc_relaxation_recon as relaxation
from kinbot.execution import atomic_json, file_sha256


PLAN = Path(__file__).resolve().parents[1]/'orbital_review_batches/r7_progress_20261024'


def read(path):
    return json.loads(Path(path).read_text())


def stage(production, review, destination, plan=PLAN, *, validate_only=False):
    source = review/'screens/ch2_sh_s_cs'
    recipe = plan/'ch2_sh_s_cs_candidates.json'
    target = destination/'candidates'
    current = screens._saved(source)
    if (target/'review_launch.json').exists():
        saved = screens._saved(target)
        if (saved.get('source_review') != str(source)
                or read(target/'recipe.json') != read(recipe)):
            raise ValueError('Existing candidate provenance or recipe differs; preserve and inspect it.')
        old = read(target/'source_manifest.json')
        now = read(current['systems'][0]['manifest'])
        if relaxation.native_evidence(old) != relaxation.native_evidence(now):
            raise ValueError('Source native evidence changed since the candidates were prepared.')
    else:
        if target.exists():
            raise ValueError('Incomplete candidate preparation; inspect the saved directory before replacing owners.')
        if not validate_only:
            screens.candidates(source, recipe, target)
    # This repeats native output/WFU checks and rejects operator-rejected paths.
    panel = relaxation.prepare(production, plan/'ch2_sh_s_cs_relaxations.json',
                               destination/'relaxations', validate_only=validate_only,
                               evidence_snapshots=([target/'source_manifest.json']
                                   if (target/'review_launch.json').exists() else []))
    receipt = dict(schema=1, production=str(production), source=str(source),
                   candidates=str(target), graphs=panel['graphs'],
                   recipe_sha256=file_sha256(recipe), production_model_accepted=False)
    if not validate_only:
        atomic_json(destination/'staged.json', receipt)
    return receipt


def perform(production, review, destination, plan=PLAN):
    destination.mkdir(parents=True, exist_ok=True)
    plan_files = [plan/name for name in ('ch2_sh_s_cs_candidates.json',
        'ch2_sh_s_cs_relaxations.json', 'ch2_sh_s_cs_probe5_4.json', 'ch2_sh_s_cs_probe7_6.json')]
    expected = {p:file_sha256(p) for p in plan_files}
    with ExitStack() as locks:
        owner = locks.enter_context((destination/'.owner.lock').open('a+'))
        try:
            fcntl.flock(owner, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('This follow-up already has an owner; nothing was submitted.') from None
        gate = locks.enter_context((review/'.handoff.lock').open('a+'))
        print('Waiting for the existing review owner; production continues.', flush=True)
        fcntl.flock(gate, fcntl.LOCK_EX)
        if any(file_sha256(p)!=digest for p,digest in expected.items()):
            raise ValueError('Diagnostic recipe changed while waiting; nothing new was submitted.')
        # Locks alone do not prove old Slurm allocations stopped. Reconcile and
        # adopt any surviving jobs through their original owners before reusing
        # the capacity. These runs retain completed work on repeated invocation.
        print('Reconciling the prior screens and eight relaxed probes.', flush=True)
        screens.run(review/'screens', report_failures=True)
        relaxation.run(review/'relaxations')
        staged = stage(production, review, destination, plan)
        print('STAGED', json.dumps(staged), flush=True)
        target = destination/'candidates'
        result = screens.run(target, report_failures=True)
        (destination/'orbital_evidence.txt').write_text(screens.evidence(target)+'\n')
        record = screens._saved(target)
        jobs = [j for s in record['systems'] for j in read(s['manifest'])['jobs']]
        if result.get('status') == 'failed' or any(
                j.get('status') != 'complete' or j.get('validation_errors') for j in jobs):
            status = dict(status='candidate_review_required',
                          reason='Incomplete or failed fixed-geometry comparison; relaxed probes remain staged.',
                          production_model_accepted=False)
        else:
            report = relaxation.run(destination/'relaxations')
            complete = all(p['task']['status']=='complete' for p in report['panels'])
            status = dict(status='review_required' if complete else 'relaxation_review_required',
                          production_model_accepted=False)
        atomic_json(destination/'followup_status.json', status)
        print(json.dumps(status), flush=True)
        return status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('production')
    args = parser.parse_args()
    production = Path(args.production).expanduser().resolve(strict=True)
    review = production.parent/'r7_progress_review_20261024'
    pin = Path(read(review/'refresh_complete.json')['code']).resolve(strict=True)
    from kinbot.rotdpy_worker_code import validate
    validate(pin)
    # Keep the previously verified production/retry implementation for these
    # workers. Only the example tools/recipes changed for this follow-up.
    import kinbot
    if Path(kinbot.__file__).resolve().parent.parent != pin:
        env = dict(os.environ, PYTHONPATH=str(pin), PYTHONDONTWRITEBYTECODE='1',
                   OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
        os.execve(sys.executable, [sys.executable, '-B', '-u', str(Path(__file__).resolve()),
                                  str(production)], env)
    perform(production, review, review/'ch2_sh_s_followup')


if __name__ == '__main__':
    main()
