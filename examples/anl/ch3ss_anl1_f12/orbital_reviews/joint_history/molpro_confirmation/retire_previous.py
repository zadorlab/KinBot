"""Preserve and audit a superseded campaign before the fixed Molpro check.

Run on Blodgett, where the user handles SSH. Only the explicitly supplied old
campaign is stopped. All native outputs and checkpoints stay at their original
paths. A possible duplicate of a new point stops the handoff for reuse review.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shutil

from kinbot.execution import atomic_json, file_sha256
from kinbot.orbital_campaign import stop_controller
from kinbot import molpro_reconnaissance as recon
from kinbot import run_cleanup
from prepare_confirmation import configurations


def matching_points(job, configs):
    """Identify candidates for reuse review; counts alone do not certify MOs."""
    if job.get('status') != 'complete':
        return []
    matches = []
    space = job.get('probe_active_space') or {}
    for config in configs:
        target = config['probe_definition']
        if (job.get('symbols') != config['symbols']
                or job.get('basis') != config['bases'][0]['basis']
                or [space.get('active_electrons'), space.get('active_orbitals')] !=
                [target['active_electrons'], target['active_orbitals']]):
            continue
        old = job.get('geometry_angstrom', [])
        for point in config['anchors']:
            new = point['geometry_angstrom']
            if (len(old) == len(new) and old
                    and all(len(row) == 3 for row in [*old, *new]) and all(
                    abs(a-b) <= 1.e-8 for row_a, row_b in zip(old, new)
                    for a, b in zip(row_a, row_b))):
                matches.append({'chain': config['title'], 'point': point['name'],
                    'review_required': 'Check actual method/shift, roots, orbital subspace, '
                        'convergence and checkpoint provenance before reusing energies.'})
    return matches


def retire(root, destination, catalog_path, *, timeout=300.):
    root = Path(root).resolve(strict=True)
    destination = Path(destination).resolve()
    if destination == root or root in destination.parents:
        raise ValueError('Preservation report must be outside the campaign being stopped.')
    if destination.exists():
        raise FileExistsError(destination)
    if not (root / 'campaign.json').is_file() or not (root / '.kinbot_run.json').is_file():
        raise ValueError('Expected an explicitly selected KinBot campaign directory.')
    catalog = json.loads(Path(catalog_path).read_text())
    seed_evidence = []
    for case in catalog['cases']:
        seed = case['seed']
        if file_sha256(seed['path']) != seed['sha256']:
            raise ValueError('Required confirmation seed missing or changed: ' + seed['path'])
        seed_evidence.append({'case': case['id'], 'path': seed['path'], 'sha256': seed['sha256']})
    configs = configurations(catalog)
    campaign = json.loads((root / 'campaign.json').read_text())
    manifests = sorted({Path(s[k]).resolve() for s in campaign['systems']
                        for k in ('pilot_manifest', 'manifest') if s.get(k)})
    for path in manifests:
        path.relative_to(root)
        if not path.is_file():
            raise FileNotFoundError(path)
    # A failed scheduler query stops here, before signals or cancellation.
    inventory = run_cleanup.inspect(root)
    destination.mkdir(parents=True)
    atomic_json(destination / 'before_stop.json', inventory)
    atomic_json(destination / 'confirmation_catalog.json', catalog)
    for path in [root / 'campaign.json', root / 'campaign_settings.json', *manifests]:
        if path.exists():
            target = destination / 'metadata_before' / path.relative_to(root)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)
    atomic_json(destination / 'transition.json', {
        'status': 'stopping_superseded_campaign', 'source': str(root),
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'Only jobs/processes inside this old campaign; no data moved or deleted.',
        'reviewed_seeds_reused': seed_evidence})
    # Prevent new submissions, then use verified WorkDir/PID identity scoping.
    stop_controller(root)
    stopped = run_cleanup.stop(root, archive=False, timeout=timeout)
    atomic_json(destination / 'after_stop.json', stopped)
    observations = []
    candidates = []
    for path in manifests:
        manifest = recon.audit(path)
        atomic_json(destination / 'audited' / path.relative_to(root), manifest)
        for job in manifest['jobs']:
            directory = Path(job['directory'])
            row = {key: job.get(key) for key in ('name', 'directory', 'status', 'basis',
                'geometry_angstrom', 'symbols', 'probe_active_space', 'state_model',
                'root_states', 'correlation_checks', 'state_energies_hartree',
                'correlated_energies_hartree', 'natural_occupations', 'input_sha256',
                'output_sha256', 'wfu', 'wfu_sha256', 'orbital_xml',
                'validation_errors', 'validation_warnings')}
            row['manifest'] = str(path)
            execution = directory / 'execution.json'
            row['execution'] = json.loads(execution.read_text()) if execution.exists() else None
            output = directory / (job['name'] + '.out')
            text = output.read_text(errors='replace') if output.exists() else ''
            row['printed_partial_reconc_values'] = re.findall(
                r'(?im)^\s*RECONC\d+\s*=.*$', text)
            row['partial_values_are_validated_results'] = job.get('status') == 'complete'
            matches = matching_points(job, configs)
            if matches:
                candidates.append({'source': str(directory), 'possible_reuse': matches})
            if job.get('status') != 'complete' and text:
                lines = text.splitlines()
                chosen = set(range(max(0, len(lines)-45), len(lines)))
                for i, line in enumerate(lines):
                    if re.search(r'NO CONVERGENCE|NOT CONVERGED|\?\s*Error|GLOBAL ERROR', line, re.I):
                        chosen.update(range(max(0, i-5), min(len(lines), i+9)))
                detail = destination / 'diagnostics' / path.parent.relative_to(root) / (job['name'] + '.txt')
                detail.parent.mkdir(parents=True, exist_ok=True)
                detail.write_text('\n'.join(f'{i+1}: {lines[i]}' for i in sorted(chosen)) + '\n')
            observations.append(row)
    result = {'schema': 1, 'source': str(root),
        'status': 'reuse_review_required' if candidates else 'retired_ready_for_confirmation',
        'completed_utc': datetime.now(timezone.utc).isoformat(),
        'counts': dict(Counter(row['status'] for row in observations)),
        'reviewed_seeds_reused': seed_evidence, 'possible_duplicate_points': candidates,
        'native_files_preserved_in_place': True,
        'production_approved': False, 'observations': observations}
    atomic_json(destination / 'evidence.json', result)
    atomic_json(destination / 'transition.json', {k: v for k,v in result.items() if k != 'observations'})
    print(json.dumps({k: v for k,v in result.items() if k != 'observations'}, indent=2))
    print('Retained results and failure diagnostics:', destination)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('old_campaign', type=Path)
    parser.add_argument('preservation_directory', type=Path)
    parser.add_argument('--catalog', type=Path,
        default=Path(__file__).with_name('geometry_and_seed_catalog.json'))
    parser.add_argument('--timeout', type=float, default=300.)
    args = parser.parse_args()
    result = retire(args.old_campaign, args.preservation_directory, args.catalog, timeout=args.timeout)
    raise SystemExit(2 if result['possible_duplicate_points'] else 0)
