"""Durable, bounded fragment-seeded retries without publishing scan gaps."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re

import numpy as np

from kinbot.execution import atomic_json, file_sha256
from kinbot.vrc_multistate import SurfaceModelRequired


class ScanPointFailure(RuntimeError):
    def __init__(self, message, graph, *, connectivity=False):
        super().__init__(message)
        self.graph = str(graph)
        self.connectivity = connectivity


class ScanReviewRequired(SurfaceModelRequired):
    def __init__(self, channels):
        self.channels = channels
        super().__init__('VRC scan paths require review: ' + ', '.join(
            c['reaction'] for c in channels))


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False,
                                    separators=(',', ':')).encode()).hexdigest()


class ScanProgress:
    """One owner per surface; persist seeds BEFORE dispatch for replay safety.

    Seeds are immutable for this scientific context, even if an earlier gap
    is repaired later. Already accepted/allocated later points then reattach
    to their original graph instead of being recreated with a different seed.
    The existing graph and dispatcher locks own all native job submissions.
    """
    def __init__(self, name, context):
        self.name = name
        safe = re.sub(r'[^A-Za-z0-9_.-]', '_', name)
        self.key = digest(context)
        self.directory = Path('vrctst/scan_progress') / safe / self.key
        self.directory.mkdir(parents=True, exist_ok=True)
        self.request = Path('vrctst/scan_review') / safe / 'scan_review.json'
        self.points = {}
        previous = self.directory/'progress.json'
        saved = json.loads(previous.read_text()) if previous.exists() else {}
        self.history = saved.get('history', saved.get('points', {}))
        self.context = context
        self.publish()

    def seed(self, base, geometry, model):
        path = self.directory / (base + '.json')
        if path.exists():
            saved = json.loads(path.read_text())
            payload = saved['seed']
            if saved.get('sha256') != digest(payload) or payload['context'] != self.key:
                raise RuntimeError(f'VRC saved scan seed changed: {path}')
        else:
            payload = dict(context=self.key, geometry=np.asarray(geometry).tolist(),
                           state_model=deepcopy(model))
            atomic_json(path, dict(seed=payload, sha256=digest(payload)))
        return np.asarray(payload['geometry'], dtype=float), deepcopy(payload['state_model'])

    def record(self, base, distance, status, **details):
        attempts = deepcopy(self.points.get(base, self.history.get(base, {})).get('attempts', []))
        for attempt in details.pop('attempts', []):
            if attempt not in attempts:
                attempts.append(attempt)
        if attempts:
            details['attempts'] = attempts
        self.points[base] = dict(distance_angstrom=float(distance), status=status, **details)
        self.publish()

    def publish(self, *, finished=False):
        gaps = {k: v for k, v in self.points.items() if v['status'] == 'needs_review'}
        flagged = any(a.get('connectivity') for p in self.points.values()
                      for a in p.get('attempts', []))
        evidence = dict(context_sha256=self.key, points=self.points)
        review_key = digest(evidence)
        approval = self.request.with_name('approval.json')
        approved = (finished and approval.exists()
                    and json.loads(approval.read_text()).get('evidence_sha256') == review_key)
        held = bool(gaps or (flagged and not approved))
        record = dict(schema=1, kind='scan_connectivity', reaction=self.name,
                      context_sha256=self.key, context=self.context, points=self.points,
                      history={**self.history, **self.points},
                      evidence_sha256=review_key, complete_grid=not gaps and finished,
                      status=('needs_review' if held else 'resolved') if finished else 'running',
                      error=('Unresolved VRC scan points: ' + ', '.join(gaps)) if gaps else
                      'Fragment connectivity changed: review the entire scan before rates.' if held else None)
        atomic_json(self.request, record)
        atomic_json(self.directory/'progress.json', record)
        return dict(kind='scan_connectivity', reaction=self.name,
                    request=str(self.request), error=record['error']) if held else None


def reject_point(graph, reason):
    """Bind a scientific rejection to a native-complete attempt for live retry."""
    graph = Path(graph)
    task = graph/'tasks/relaxed_scan'
    atomic_json(task/'scan_rejection.json', dict(schema=1, reason=reason,
        hashes={str(p.relative_to(graph)): file_sha256(p) for p in (
            graph/'workflow.json', task/'execution.json', task/'final.xyz')}))


def bind_correction_review(correction):
    """Attach provenance in a sidecar, preserving existing correction hashes."""
    correction = Path(correction)
    name = correction.stem.removeprefix('corr_')
    review = correction.parent/'scan_review'/name/'scan_review.json'
    if not review.is_file():
        return
    record = json.loads(review.read_text())
    if not record.get('complete_grid'):
        raise ValueError(f'{name}: incomplete scan cannot be published as a correction.')
    atomic_json(correction.with_suffix('.scan.json'), dict(schema=1,
        scan_evidence_sha256=record['evidence_sha256'], correction_sha256=file_sha256(correction)))


def require_ready_for_rates(correction, *, expected_sha256=None):
    """An old correction/result cannot bypass a newer scan review hold."""
    correction = Path(correction)
    name = correction.stem.removeprefix('corr_')
    review = correction.parent/'scan_review'/name/'scan_review.json'
    if not review.is_file():
        return  # Legacy or non-Molpro scan without a recovery record.
    record = json.loads(review.read_text())
    evidence = dict(context_sha256=record.get('context_sha256'), points=record.get('points'))
    if (record.get('reaction') != name or record.get('status') != 'resolved'
            or not record.get('complete_grid') or digest(evidence) != record.get('evidence_sha256')):
        raise ScanReviewRequired([dict(kind='scan_connectivity', reaction=name,
            request=str(review), error='Scan must finish and any path review must be resolved before rates.')])
    binding = correction.with_suffix('.scan.json')
    payload = json.loads(binding.read_text()) if binding.is_file() else {}
    current_sha256 = file_sha256(correction)
    if (payload.get('scan_evidence_sha256') != record['evidence_sha256']
            or payload.get('correction_sha256') != current_sha256):
        raise ValueError(f'{name}: correction is not bound to the current completed scan.')
    if expected_sha256 is not None and current_sha256 != expected_sha256:
        raise ValueError(f'{name}: sampled correction differs from the current reviewed correction.')


def rejected_point(graph, ident):
    graph = Path(graph)
    path = graph/'tasks'/ident/'scan_rejection.json'
    if ident != 'relaxed_scan' or not path.is_file():
        return False
    record = json.loads(path.read_text())
    expected = {'workflow.json', 'tasks/relaxed_scan/execution.json', 'tasks/relaxed_scan/final.xyz'}
    return (record.get('schema') == 1 and set(record.get('hashes', {})) == expected
            and all((graph/p).is_file() and file_sha256(graph/p) == value
                    for p, value in record['hashes'].items()))


def review_changed(channels):
    """Wake a held owner for explicit approval or targeted live task retries."""
    from kinbot import live_control
    for channel in channels:
        if channel.get('kind') != 'scan_connectivity':
            continue
        path = Path(channel['request'])
        record = json.loads(path.read_text())
        approval = path.with_name('approval.json')
        if (record.get('complete_grid') and approval.exists()
                and json.loads(approval.read_text()).get('evidence_sha256') == record['evidence_sha256']):
            return True
        for point in record['points'].values():
            for attempt in point.get('attempts', []):
                if live_control.requests(Path(attempt['graph'])/'tasks', kinds={'task'}):
                    return True
    return False


def approve(path, reason):
    """Explicit scientific review, scoped to the exact completed scan evidence."""
    path = Path(path).resolve(strict=True)
    record = json.loads(path.read_text())
    if (record.get('kind') != 'scan_connectivity' or not record.get('complete_grid')
            or record.get('status') not in ('needs_review', 'resolved')
            or not record.get('points')
            or any(p['status'] != 'accepted' for p in record['points'].values())):
        raise ValueError('Only a complete accepted grid can be approved; scan gaps cannot be waived.')
    if not reason.strip():
        raise ValueError('Record the scientific review and its evidence in --reason.')
    evidence = dict(context_sha256=record['context_sha256'], points=record['points'])
    if digest(evidence) != record['evidence_sha256']:
        raise ValueError('Scan review evidence changed.')
    for point in record['points'].values():
        hashes = point.get('hashes', {})
        if not hashes or any(not Path(p).is_file() or file_sha256(p) != value
                             for p, value in hashes.items()):
            raise ValueError('Accepted scan artifacts changed or lack provenance; refresh the review first.')
    from kinbot.live_control import now
    atomic_json(path.with_name('approval.json'), dict(schema=1, reason=reason,
        evidence_sha256=record['evidence_sha256'], approved_at=now()))
    return dict(status='approved', request=str(path), drivers_started=0)


def retry(path, distance, reason, fallback=None):
    """Queue an explicit retry of the last failed attempt at one scan gap."""
    from kinbot import live_control
    path = Path(path).resolve(strict=True)
    record = json.loads(path.read_text())
    if record.get('status') != 'needs_review' or record.get('kind') != 'scan_connectivity':
        raise ValueError('Wait for the scan to finish and publish its review request.')
    points = [p for p in record['points'].values() if p['distance_angstrom'] == distance
              and p['status'] == 'needs_review']
    if len(points) != 1 or not points[0].get('attempts'):
        raise ValueError('Select exactly one unresolved point; accepted points are retained.')
    graph = Path(points[0]['attempts'][-1]['graph']).resolve(strict=True)
    if not graph.is_relative_to(path.parents[2]/'molpro_scan'):
        raise ValueError('Retry graph is outside this well\'s VRC scan directory.')
    policy = live_control.find(graph)
    if policy is None:
        raise ValueError('Live controls are not enabled for this run.')
    return live_control.request(policy.parent, graph/'tasks/relaxed_scan', reason=reason,
                                fallback=fallback)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['status', 'approve', 'retry'])
    parser.add_argument('request', type=Path)
    parser.add_argument('--reason')
    parser.add_argument('--distance', type=float)
    parser.add_argument('--fallback', type=Path, help='Optional schema-1 manual fallback JSON.')
    args = parser.parse_args()
    if args.action == 'status':
        record = json.loads(args.request.read_text())
        result = {k: record.get(k) for k in ('reaction', 'status', 'error', 'complete_grid', 'points')}
    else:
        if not args.reason:
            parser.error('--reason is required for a review/retry action')
        if args.action == 'retry':
            if args.distance is None:
                parser.error('--distance is required for retry')
            fallback = json.loads(args.fallback.read_text()) if args.fallback else None
            result = retry(args.request, args.distance, args.reason, fallback)
        else:
            result = approve(args.request, args.reason)
    print(json.dumps(result, indent=2))
