"""Persistent operator exclusions, applied by the live network owner.

Raw accepted chemistry stays intact. Inclusion removes an operator exclusion;
it never turns a failed calculation into accepted chemistry. Selection and VRC
both consume this policy, independently of search checkpoint/retry state.
"""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re

from kinbot.execution import atomic_json, file_sha256

NAME = '.kinbot_reaction_decisions.json'
MARKER = '.kinbot_reaction_rejection.json'


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def _validate_marker(record, spec):
    tasks = {t['id']:t for t in spec['tasks']}
    if (record['molecule_sha256'] != _digest(spec['molecule'])
            or any(n not in tasks or _digest(tasks[n]) != h for n,h in record['task_contracts'].items())):
        raise ValueError('Excluded task changed after its decision; reconcile before inclusion.')


def read(worker):
    path = Path(worker)/NAME
    return json.loads(path.read_text()) if path.exists() else {'schema': 1, 'reactions': {}}


def rejected(worker, name):
    return read(worker)['reactions'].get(name, {}).get('decision') == 'reject'


def validate(worker, decisions):
    if (not isinstance(decisions, dict) or not decisions or any(
            not isinstance(n, str) or not re.fullmatch(r'[A-Za-z0-9_.-]+', n)
            or d not in ('reject', 'include') for n, d in decisions.items())):
        raise ValueError('Reaction decisions must map exact names to reject or include.')
    evidence = json.loads((Path(worker)/'accepted_l2_network.json').read_text())
    names = {r['name'] for r in evidence['network']['reactions']}
    if set(decisions) - names:
        raise ValueError('Decisions require previously accepted reaction names; failed chemistry cannot be waived.')
    return decisions


def owns(name, channel):
    return name == channel or name.startswith(channel+'__') or name.startswith(channel+'_vts_')


def graphs(worker, channel):
    # Only dedicated VRC graphs; shared ANL/species/fragment wells are retained.
    for folder in ('molpro_scan', 'molpro', 'molpro_pivots'):
        base = Path(worker)/'vrctst'/folder
        for p in sorted(base.glob('*/workflow.json')):
            spec = json.loads(p.read_text())
            if owns(p.parent.name, channel) or any(owns(t['id'], channel) for t in spec['tasks']):
                yield p.parent


def campaign_inputs(worker, channel):
    for p in sorted(Path(worker).parent.glob('*/'+channel+'*.campaign.json')):
        record = json.loads(p.read_text())
        inp = Path(record['input']).resolve()
        from kinbot.rotdpy_campaign import receipt
        if receipt(inp) != p.resolve() or not owns(inp.stem, channel):
            continue
        yield inp


def apply(worker, decisions, path, request):
    import fcntl
    from contextlib import ExitStack
    with ExitStack() as stack:
        for graph in sorted({g for name in decisions for g in graphs(worker, name)}):
            handle = stack.enter_context((graph/'drive.lock').open('a+'))
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                return False
        return _apply(worker, decisions, path, request)


def _apply(worker, decisions, path, request):
    """Apply after the existing PES epoch barrier releases branch ownership.

    Returns False while cancellation is settling. Journals and exact process /
    scheduler ownership checks make re-entry safe after an interrupted refresh.
    Native artifacts remain in place; derived selections are rebuilt by PES.
    """
    from kinbot import live_control as control, rotdpy_campaign as campaigns
    from kinbot.anl import dispatch, live_control
    worker = Path(worker).resolve()
    validate(worker, decisions)
    archive = worker/'.kinbot_controls/reaction_decisions'/request['id']
    graph_list, campaign_list = [], []
    for name, decision in decisions.items():
        graph_list.extend((name, decision, g) for g in graphs(worker, name))
        campaign_list.extend((name, decision, p) for p in campaign_inputs(worker, name))
    # Preflight before stopping anything. Existing accepted chemistry is the
    # authority for namespace ownership, never an arbitrary path prefix.
    for name, decision, graph in graph_list:
        _, _, state = dispatch._load(graph)
        for ident, entry in state['tasks'].items():
            if not (owns(graph.name, name) or owns(ident, name)):
                continue
            if entry.get('status') == 'submitting':
                raise ValueError(f'{graph}: ambiguous submission must be reconciled first.')
            job = entry.get('job_id')
            if decision == 'reject' and job and dispatch._job_active(str(job)):
                live_control._ownership(graph/'tasks'/ident, job)
                if not request['cancel_running']:
                    raise ValueError('Channel has live work; specify --cancel-running to reject it.')
    for _, decision, inp in campaign_list:
        record = json.loads(campaigns.receipt(inp).read_text())
        campaigns.validate(inp, record)
        if decision == 'reject' and campaigns.active(inp, record) and not request['cancel_running']:
            raise ValueError('Channel has a live sampling owner; specify --cancel-running.')
    control.acknowledge(path, request, 'applying', archive=str(archive))
    archive.mkdir(parents=True, exist_ok=True)
    if not (archive/'previous.json').exists():
        atomic_json(archive/'previous.json', read(worker))
    for name, decision, inp in campaign_list:
        if decision == 'reject' and not campaigns.stop_for_rejection(inp, archive, request):
            return False
        if decision == 'include':
            with campaigns.locked(inp):
                receipt = campaigns.receipt(inp)
                record = json.loads(receipt.read_text())
                if record.pop('operator_rejected', None):
                    record['operator_reinclude'] = request['id']
                atomic_json(receipt, record)
    for name, decision, graph in graph_list:
        _, spec, state = dispatch._load(graph)
        saved = archive/graph.relative_to(worker)
        if not (saved/'state.json').exists():
            atomic_json(saved/'state.json', state)
            atomic_json(saved/'workflow.json', spec)
        markers = json.loads((graph/MARKER).read_text()) if (graph/MARKER).exists() else {'schema': 1, 'reactions': {}}
        if decision == 'reject' and name in markers['reactions']:
            prior = markers['reactions'][name]
            _validate_marker(prior, spec)
            # Finish an interrupted marker-before-state publication.
            for ident in prior['tasks']:
                state['tasks'].setdefault(ident, {}).update(status='skipped', manual_rejection=prior['request'])
            atomic_json(graph/'state.json', state)
            continue
        if decision == 'reject':
            for ident, entry in state['tasks'].items():
                if not (owns(graph.name, name) or owns(ident, name)):
                    continue
                job = entry.get('job_id')
                if not job or not dispatch._job_active(str(job)):
                    continue
                cancel_path = saved/(str(job)+'.cancel.json')
                cancel = (json.loads(cancel_path.read_text()) if cancel_path.exists()
                          else dict(request, kind='reaction_decision', cancel_requested=False))
                if not live_control.stop_attempt(graph/'tasks'/ident, job,
                        cancel_path, cancel, allow_completed=True):
                    return False
            targets = [t['id'] for t in spec['tasks'] if owns(graph.name, name) or owns(t['id'], name)]
            markers['reactions'][name] = dict(reaction=name, tasks=targets,
                request=request['id'], reason=request['reason'], archive=str(saved),
                workflow_sha256=file_sha256(graph/'workflow.json'), molecule_sha256=_digest(spec['molecule']),
                task_contracts={t['id']:_digest(t) for t in spec['tasks'] if t['id'] in targets})
            atomic_json(graph/MARKER, markers)
            for ident in targets:
                state['tasks'].setdefault(ident, {}).update(status='skipped', manual_rejection=request['id'])
            atomic_json(graph/'state.json', state)
        elif name in markers['reactions']:
            prior = markers['reactions'].pop(name)
            _validate_marker(prior, spec)
            atomic_json(saved/'rejection.json', prior)
            old = json.loads((Path(prior['archive'])/'state.json').read_text())
            for ident in prior['tasks']:
                if ident in old['tasks']:
                    state['tasks'][ident] = deepcopy(old['tasks'][ident])
                else:
                    state['tasks'].pop(ident, None)
            atomic_json(graph/'state.json', state)
            atomic_json(graph/MARKER, markers)
    updated = deepcopy(read(worker))
    for name, decision in decisions.items():
        updated['reactions'][name] = dict(decision=decision, reason=request['reason'],
            request=request['id'], created_at=request['created_at'])
    atomic_json(worker/NAME, updated)
    return True


def graph_rejected(graph):
    marker = Path(graph)/MARKER
    if not marker.exists():
        return False
    records = json.loads(marker.read_text())['reactions']
    from kinbot.anl.dispatch import _load
    _, spec, _ = _load(graph)
    targets = set()
    for record in records.values():
        _validate_marker(record, spec)
        targets.update(record['tasks'])
    return targets == {t['id'] for t in spec['tasks']}
