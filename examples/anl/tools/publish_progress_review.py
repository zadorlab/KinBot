#!/usr/bin/env python3
"""Retire this saved review, refresh retained drivers, publish five selections.

This example only orchestrates shared controls and the explicit scientific plan.
It never launches reconnaissance, accepts scan results, or edits native evidence.
Run from the updated checkout using its virtual environment and PYTHONPATH.
"""
import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess

from kinbot import requested_orbital_recon as recon
from kinbot.execution import atomic_json
import queue_progress_review as queue
import export_orbital_evidence as exporter

CHECKOUT = Path(__file__).resolve().parents[3]
PLAN = CHECKOUT/'examples/anl/orbital_review_batches/r7_progress_20261024'


def perform(production, export_directory):
    production = Path(production).expanduser().resolve(strict=True)
    review = production.parent/'r7_progress_review_20261024'
    export_directory = Path(export_directory).expanduser().resolve()
    if export_directory.exists():
        raise FileExistsError('Use a new evidence export directory.')
    if any(export_directory.is_relative_to(p) for p in (CHECKOUT,production.parent)):
        raise ValueError('Evidence destination must be outside checkout and run roots.')
    with ExitStack() as locks:
        queue.acquire(locks, production/'.orbital_publication.lock', timeout=0)
        print('Stopping standalone reviews and restoring production capacity to 23.', flush=True)
        stopped = queue.retire(production, review)
        print('Checking all five selections before replacing production owners.', flush=True)
        preview = recon.review_batch(production, production.parent,
                                    PLAN/'all_selections.json', require_all=True)
        if any(r['status'].startswith('waiting') for r in preview['cases'].values()):
            raise RuntimeError('Selections are not all ready; drivers were not replaced: '+json.dumps(preview))
        print('Refreshing production code; existing QC jobs and sampling owners are retained.', flush=True)
        subprocess.run(['bash', str(CHECKOUT/'examples/anl/tools/refresh_retained_run.sh'),
                        str(production), '23'], cwd=CHECKOUT, check=True)
        published = recon.review_batch(production, production.parent,
                                       PLAN/'all_selections.json', require_all=True, apply=True)
        atomic_json(review/'publication_handoff.json', dict(schema=1, retirement=stopped,
            publication=published, export_directory=str(export_directory)))
        print('PUBLICATION RECEIPT '+json.dumps({case:dict(status=r['status'],
            channels=len(r.get('channels',[])),refresh=r.get('refresh'))
            for case,r in published['cases'].items()}), flush=True)
        if any(r['status']=='blocked' or r['status'].startswith('waiting')
               for r in published['cases'].values()):
            raise RuntimeError('Publication needs attention; see the saved receipt.')
        export_directory.mkdir(parents=True)
        report = exporter.export(PLAN/'export.json',
            dict(repo=CHECKOUT, runs=production.parent, production=production),
            export_directory/'orbital-evidence.tar.gz')
        print('EVIDENCE DIRECTORY '+str(export_directory), flush=True)
        if report['issues']:
            raise RuntimeError('Archive retained with explicit collection issues: '+json.dumps(report['issues']))
        print('Selection publication submitted and evidence archived. A queued receipt still awaits '
              'the owning drivers; this is not scan or rate completion.', flush=True)
        return published


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('production')
    parser.add_argument('--export-directory')
    args = parser.parse_args()
    destination = args.export_directory or (Path.home()/'kinbot-evidence'/
        ('vrc-selection-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+f'-{os.getpid()}'))
    perform(args.production, destination)


if __name__ == '__main__':
    main()
