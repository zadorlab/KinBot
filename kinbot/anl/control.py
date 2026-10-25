"""User-facing live ANL sharing and protocol selection (preview by default)."""
import argparse
import fnmatch
import json
from pathlib import Path

from kinbot.anl import protocol, sharing


def graphs(root, *, patterns=(), species=(), all_species=False):
    if not (patterns or species or all_species):
        raise ValueError('Select --graph, --species or --all explicitly.')
    result = []
    for p in sorted(root.rglob('workflow.json')):
        relative = p.parent.relative_to(root)
        if any(x == 'attempts' or x.startswith('.') for x in relative.parts):
            continue
        spec = json.loads(p.read_text())
        if spec.get('name') != 'anl-automatic-post-l2':
            continue
        mol = spec['molecule']
        if (all_species or any(fnmatch.fnmatch(str(relative), s) for s in patterns)
                or mol.get('smiles') in species):
            result.append(p.parent)
    if not result:
        raise ValueError('No automatic ANL graphs matched the selection.')
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    s = sub.add_parser('share', help='Configure an explicit source PES → receiving PES route.')
    s.add_argument('root', type=Path); s.add_argument('--source', type=Path, required=True)
    s.add_argument('--target', type=Path, required=True); s.add_argument('--apply', action='store_true')
    p = sub.add_parser('protocol', help='Preview or publish an owner-applied selection of graph changes.')
    p.add_argument('root', type=Path)
    p.add_argument('--graph', action='append', default=[], help='Relative graph glob; repeat for subsets.')
    p.add_argument('--species', action='append', default=[], help='Exact accepted SMILES; includes its declared states.')
    p.add_argument('--all', action='store_true')
    p.add_argument('--tier', choices=protocol.TIERS)
    p.add_argument('--geometry', choices=('default', 'tz', 'qz', 'l2'))
    p.add_argument('--harmonic', choices=('default', 'tz', 'qz', 'cbs', 'l2'))
    p.add_argument('--omit-task', action='append')
    p.add_argument('--omit-component', action='append', choices=sorted(protocol.FAMILIES))
    p.add_argument('--clear-omissions', action='store_true')
    p.add_argument('--anharmonic-qz', action=argparse.BooleanOptionalAction, default=None)
    p.add_argument('--cancel-running', action='store_true')
    p.add_argument('--reason', default='Explicit user ANL protocol selection')
    p.add_argument('--apply', action='store_true')
    r = sub.add_parser('reopen', help='Ask an existing PES owner to reassemble a handed-off worker.')
    r.add_argument('root', type=Path); r.add_argument('--worker', type=Path, action='append', required=True)
    r.add_argument('--reason', required=True); r.add_argument('--apply', action='store_true')
    t = sub.add_parser('status'); t.add_argument('root', type=Path)
    t = sub.add_parser('wake-owners', help='Retry owner wakeups for queued protocol requests.')
    t.add_argument('root', type=Path)
    args = parser.parse_args(argv); root = args.root.resolve(strict=True)
    if args.command == 'share':
        result = sharing.configure(root, root/args.source, root/args.target, apply=args.apply)
    elif args.command == 'protocol':
        options = {k:getattr(args, k) for k in ('tier', 'geometry', 'harmonic', 'anharmonic_qz')
                   if getattr(args, k) is not None}
        for key, value in (('omit_tasks', args.omit_task), ('omit_components', args.omit_component)):
            if value is not None or args.clear_omissions:
                options[key] = value or []
        selected = graphs(root, patterns=args.graph, species=args.species, all_species=args.all)
        plans = [protocol.plan(g, options, allow_closed=True) for g in selected]
        if args.apply:
            result = [dict(id=r['id'], graph=r['graph'], status=r['status']) for r in
                      protocol.publish(root, plans, reason=args.reason, cancel_running=args.cancel_running)]
        else:
            result = [{k:v for k,v in row.items() if k not in ('replacement', 'expected')}
                      | dict(sources=protocol.sources(row['replacement']),
                             omitted_tasks=sorted(protocol.omitted(row['replacement']))) for row in plans]
    elif args.command == 'reopen':
        from kinbot import live_control
        workers = [(root/w).resolve(strict=True) for w in args.worker]
        generations = [live_control.generation(w, 'network') for w in workers]
        result = ([live_control.request(root, w, kind='network', reason=args.reason, options={}, expected=e)
                   for w,e in zip(workers, generations)] if args.apply else
                  [dict(worker=str(w), status='preview', expected=e) for w,e in zip(workers, generations)])
    else:
        if args.command == 'wake-owners':
            protocol.wake_owners(root)
        result = dict(sharing=sharing.configuration(root), requests=[
            {k:v for k,v in json.loads(p.read_text()).items() if k != 'plan'}
            for p in sorted((root/protocol.QUEUE).glob('*.json'))])
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
