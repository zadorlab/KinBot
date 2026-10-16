"""Prepare the corrected confirmation, then stop only its superseded run.

All old evidence remains in place. A partially prepared replacement may be
reused only before any jobs have been submitted. No SSH or user-wide kills.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

from kinbot import run_cleanup
from kinbot.execution import atomic_json, file_sha256
from prepare_confirmation import build
from run_confirmation import validated_record


def prepare_and_stop(old, new, catalog):
    old, new, catalog = (Path(p).resolve() for p in (old, new, catalog))
    if old == new or old in new.parents or new in old.parents:
        raise ValueError('Old and replacement runs must be separate, non-nested directories.')
    validated_record(old)
    run_cleanup.inspect(old)  # Scheduler failure must not look like an empty queue.
    if not new.exists():
        build(catalog, new)  # Seed/input checks precede cancellation.
    replacement = validated_record(new)
    if replacement['catalog_sha256'] != file_sha256(catalog):
        raise ValueError('Replacement catalogue differs; use a fresh directory.')
    if replacement['status'] != 'prepared_review_required' or replacement.get('controller_pid'):
        raise ValueError('Replacement already started; use its run/resume command.')
    for system in replacement['systems']:
        for job in json.loads(Path(system['manifest']).read_text())['jobs']:
            if job.get('job_id') or job.get('submission_status') or job.get('status') == 'complete':
                raise ValueError('Replacement already contains submitted/completed jobs.')
            inp = Path(job['directory']) / (job['name'] + '.inp')
            if file_sha256(inp) != job['input_sha256']:
                raise ValueError('Replacement input changed: ' + str(inp))
    stopped = run_cleanup.stop(old, archive=False, timeout=300.)
    original = validated_record(old)
    manifests = [Path(s['manifest']) for s in original['systems']]
    jobs = [j for p in manifests for j in json.loads(p.read_text())['jobs']]
    report = {'schema': 1, 'time_utc': datetime.now(timezone.utc).isoformat(),
              'old_directory': str(old), 'replacement_directory': str(new),
              'status': stopped['status'], 'old_files_deleted': False,
              'completed_results_preserved': sum(j.get('status') == 'complete' for j in jobs),
              'previous_job_counts': dict(Counter(j.get('status', 'unknown') for j in jobs)),
              'previous_registry': original,
              'previous_manifests': {str(p): json.loads(p.read_text()) for p in manifests},
              'reuse_policy': 'Retain all old results as evidence; rerun the corrected protocol consistently.'}
    atomic_json(new / 'previous_confirmation.json', report)
    return {k: v for k, v in report.items() if k not in {'previous_registry', 'previous_manifests'}}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('old', type=Path)
    p.add_argument('new', type=Path)
    p.add_argument('--catalog', type=Path,
                   default=Path(__file__).with_name('geometry_and_seed_catalog.json'))
    args = p.parse_args()
    print(json.dumps(prepare_and_stop(args.old, args.new, args.catalog), indent=2))
