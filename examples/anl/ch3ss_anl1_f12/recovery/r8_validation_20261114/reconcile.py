"""Prepare/adopt the finite evidence-driven reconciliation panel.

prepare writes one immutable decision; only resume.py's locked owner applies
it. review is read-only. Original attempts, scientific gates and MESS stay put.
"""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import fcntl
import json
from pathlib import Path
import uuid

from kinbot.anl import dispatch, live_control, path_review
from kinbot.execution import atomic_json, file_sha256

DECISION = 'reconciliation_v1.json'
IRC_PILOTS = (
    ('irc_native_ch3ss_rxn_12_2', 'forward_10', 'dvv'),
    ('irc_native_ch3ss_rxn_12_2', 'reverse_10', 'dvv'),
    ('irc_paper_ts2', 'forward_5', 'dvv'),
    ('irc_paper_ts1', 'forward_5', 'pc150'),
    ('irc_paper_ts5', 'reverse_10', 'dvv'),
)
CC_AUDITS = ('cc_si_TS2', 'cc_si_TS3')


def load(root):
    p = Path(root)/DECISION
    if not p.exists(): return None
    result = json.loads(p.read_text()); body = {k:v for k,v in result.items() if k != 'sha256'}
    if result.get('schema') != 1 or result.get('sha256') != path_review.digest(body):
        raise ValueError('Reconciliation decision checksum is invalid.')
    if result['root'] != str(Path(root).resolve()): raise ValueError('Decision belongs to another panel.')
    return result


def proofs(plan):
    return [x['proof'] for x in plan['irc']] + [x['proof'] for x in plan['curvature'].values()]


def prepare(root):
    root = Path(root).resolve(strict=True)
    existing = load(root)
    if existing:
        for proof in proofs(existing): path_review.verify(proof)
        return existing  # same decision, even after probes finish/fail
    stage = json.loads((root/'stage.json').read_text())
    probes = [dict(graph=name, **path_review.irc_probe(root/name, task, recipe))
              for name,task,recipe in IRC_PILOTS]
    audits = {name+'_curvature_v1':path_review.curvature_probe(root/name, 'refine') for name in CC_AUDITS}
    # Explicitly stop launching the remaining original independent IRC grid.
    # Already allocated siblings are still polled and their output is adopted.
    held = [name for name in stage['workflows'] if name.startswith('irc_')]
    result = dict(schema=1, root=str(root), source_stage_sha256=file_sha256(root/'stage.json'),
        frozen_sources=stage['sources'], irc=probes, curvature=audits, held_irc_graphs=held,
        created_utc=datetime.now(timezone.utc).isoformat(), scientific_approval=False)
    result['sha256'] = path_review.digest(result)
    # Refuse simultaneous preparations, rather than overwriting a reviewed
    # decision that another owner could already be using.
    with (root/'.reconciliation_prepare.lock').open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        existing = load(root)
        if existing: return existing
        atomic_json(root/DECISION, result)
    return result


def extra_specs(plan):
    return {name:entry['spec'] for name,entry in plan['curvature'].items()} if plan else {}


def stage_extra(root, name, spec):
    """Recover only our own interrupted prepare; never reset native attempts."""
    graph = root/name; marker = root/(name+'.preparing.json')
    if graph.exists() and not (graph/'state.json').exists():
        if (not marker.is_file() or json.loads(marker.read_text()) != spec
                or any(graph.glob('tasks/*/execution.json')) or any(graph.glob('tasks/*/*.jobid'))):
            raise ValueError('Unowned or executed partial curvature preparation.')
        archive = root/'interrupted_preparations'; archive.mkdir(exist_ok=True)
        graph.rename(archive/(name+'.'+uuid.uuid4().hex))
    if not graph.exists():
        atomic_json(marker, spec)
        dispatch.prepare(deepcopy(spec), graph)
    _, saved, _ = dispatch._load(graph)
    normalized = deepcopy(saved)
    resources = {t['id']:t['resources'] for t in spec['tasks']}
    for task in normalized['tasks']: task['resources'] = resources[task['id']]
    if normalized != spec: raise ValueError('Curvature diagnostic definition changed.')
    receipt = graph/'requested.json'
    if receipt.exists() and json.loads(receipt.read_text()) != spec:
        raise ValueError('Curvature request changed.')
    atomic_json(receipt, spec); marker.unlink(missing_ok=True)


def apply(root, plan):
    """Called under controller and all affected graph drive locks."""
    if plan is None: return
    stage = json.loads((root/'stage.json').read_text())
    if stage['sources'] != plan['frozen_sources']:
        raise ValueError('Frozen source definitions differ from reconciliation decision.')
    for proof in proofs(plan): path_review.verify(proof)
    grouped = {}
    for probe in plan['irc']: grouped.setdefault(probe['graph'], []).append(probe['task'])
    # Preflight all graphs before committing any. Per-graph journals then make
    # interrupted multi-graph adoption restartable without repeating trials.
    changes = []
    for name, tasks in grouped.items():
        graph, spec, state = dispatch._load(root/name)
        replacement = deepcopy(spec); current = {t['id']:t for t in spec['tasks']}
        for task in tasks:
            if task['id'] in current:
                saved = deepcopy(current[task['id']]); saved['resources'] = task['resources']
                if saved != task: raise ValueError('Existing probe differs from the reviewed decision.')
            else: replacement['tasks'].append(deepcopy(task))
        if replacement != spec:
            live_control.check_commit(graph, replacement, state, set())
            changes.append((graph,replacement,state))
    for graph, replacement, state in changes:
        live_control.commit(graph, replacement, state, set(), {})
    for name, spec in extra_specs(plan).items(): stage_extra(root, name, spec)
    atomic_json(root/'reconciliation_applied.json', dict(sha256=plan['sha256'], scientific_approval=False))


def submit_only(plan, name):
    if plan and name in plan['held_irc_graphs']:
        return {x['task']['id'] for x in plan['irc'] if x['graph'] == name}
    return None


def review(root):
    """Read logs/receipts without dispatch.advance (which stages/reconciles)."""
    root = Path(root).resolve(); plan = load(root)
    stage = json.loads((root/'stage.json').read_text())
    result = dict(scientific_approval=False, decision=plan['sha256'] if plan else None, branches={}, curvature={})
    for name in stage['workflows']:
        if not name.startswith('irc_'): continue
        _, spec, state = dispatch._load(root/name)
        for task in spec['tasks']:
            if 'IRC=(' not in task.get('input_template',''): continue
            directory = root/name/'tasks'/task['id']; path = directory/task['stdout']
            entry = state['tasks'].get(task['id'], {})
            row = dict(state=entry.get('status', 'waiting'), job_id=entry.get('job_id'))
            if path.is_file():
                with path.open(errors='replace') as stream: report = path_review.gaussian_path(stream)
                rows = report.pop('last_corrector')
                row.update(report, last_corrector=rows[-1:] , log_sha256=file_sha256(path))
            result['branches'][name+'/'+task['id']] = row
    for name in extra_specs(plan):
        graph = root/name
        if not (graph/'state.json').is_file(): continue
        _, spec, state = dispatch._load(graph)
        result['curvature'][name] = dict(tasks=state['tasks'])
        for task in spec['tasks']:
            path = graph/'tasks'/task['id']/(task['id']+'.out')
            if path.is_file():
                result['curvature'][name][task['id']] = dict(modes=path_review.molpro_modes(path.read_text(errors='replace')),
                    output_sha256=file_sha256(path))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare','review'))
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    result = prepare(args.root) if args.action == 'prepare' else review(args.root)
    if args.action == 'prepare':
        result = dict(decision=result['sha256'], irc_trials=len(result['irc']),
            curvature_jobs=sum(len(x['spec']['tasks']) for x in result['curvature'].values()),
            status='prepared_for_owner', scientific_approval=False)
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
