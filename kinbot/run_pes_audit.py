"""Read-only topology and duplicate evidence, including unfinished PES workers.

Writes reports only to an explicit directory outside the run. SQLite is opened
in mode=ro; no QC wrapper, dispatcher, Slurm command, or recovery is invoked.
Cache comparisons are diagnostic candidates, never authorization for reuse.
"""
import argparse
from collections import Counter, defaultdict
from copy import deepcopy
from datetime import datetime, timezone
from itertools import combinations
import json
import math
from pathlib import Path
import sqlite3

from kinbot.anl import ts_reuse


def endpoints(row):
    return tuple(sorted(((row['reactant'],), tuple(sorted(row['products'])))))


def numerical_differences(a, b):
    """Small-molecule diagnostic only: exhaustive proper-rotation alignment.

    Stop before combinatorial growth. These numbers do not change the reuse
    policy and cannot justify treating distinct saddles as the same reaction.
    """
    import networkx as nx
    import numpy as np
    result = dict(energy_gap_kcal_mol=abs(a['energy_hartree']-b['energy_hartree'])*627.509474,
        zero_k_gap_kcal_mol=abs(a['energy_hartree']+a['harmonic_zpe_hartree']-
                              b['energy_hartree']-b['harmonic_zpe_hartree'])*627.509474)
    mx,my = a['molecule'],b['molecule']
    counts = Counter(mx['symbols'])
    if counts!=Counter(my['symbols']):
        return result
    permutations = math.prod(math.factorial(n) for n in counts.values())
    if permutations>10000:
        return dict(result, alignment='not_enumerated; more than 10000 element permutations')
    px,py = (np.asarray(m['positions'],float) for m in (mx,my))
    if px.shape!=py.shape or px.shape!=(len(mx['symbols']),3) or not np.isfinite(px).all() or not np.isfinite(py).all():
        raise ValueError('Cannot align malformed or nonfinite coordinates.')
    px,py = px-px.mean(axis=0),py-py.mean(axis=0)
    def graph(symbols):
        g = nx.Graph()
        g.add_nodes_from((i,{'element':s}) for i,s in enumerate(symbols))
        return g
    matcher = nx.algorithms.isomorphism.GraphMatcher(graph(mx['symbols']),graph(my['symbols']),
        node_match=lambda x,y:x['element']==y['element'])
    best = float('inf')
    for mapping in matcher.isomorphisms_iter():
        target = py[[mapping[i] for i in range(len(px))]]
        u,_,vt = np.linalg.svd(px.T@target)
        correct = np.eye(3)
        correct[-1,-1] = np.linalg.det(u@vt)
        best = min(best,float(np.sqrt(np.mean(np.sum((px@u@correct@vt-target)**2,axis=1)))))
    return dict(result,minimum_proper_rotation_rmsd_angstrom=best,element_permutations=permutations)


def read(path, errors, default=None):
    try:
        if not path.is_file():
            return {} if default is None else default
        return json.loads(path.read_text())
    except (OSError, ValueError) as error:
        errors.append(dict(file=str(path), error=str(error)))
        return {} if default is None else default


def input_label(inp):
    if inp.get('smiles'):
        return inp['smiles'], 'input_smiles'
    if not inp.get('structure'):
        return '?', 'unknown'
    from kinbot.stationary_pt import StationaryPoint
    point = StationaryPoint('audit', inp.get('charge', 0), inp.get('mult', 1),
                            structure=inp['structure'])
    point.characterize()
    return point.smiles, 'input_geometry_connectivity_only'


def database_record(database, name, inp, label):
    """Latest row only; never resurrect an older success after a newer failure."""
    from ase.db.sqlite import SQLite3Database
    from kinbot import constants
    class ReadOnly(SQLite3Database):
        def _connect(self):
            return sqlite3.connect(Path(self.filename).resolve().as_uri()+'?mode=ro',
                                   uri=True, timeout=3)
    db = ReadOnly(str(database), use_lock_file=False)
    rows = list(db.select(name=name, sort='-id', limit=1))
    if not rows:
        return None
    row = rows[0]
    data = row.data
    if data.get('status') != 'normal' or any(k not in data for k in ('energy','zpe','frequencies')):
        return None
    atom = row.toatoms()
    parameters = dict(row.get('calculator_parameters') or {})
    for key in ('label','chk','oldchk','mem','nprocshared','command','directory'):
        parameters.pop(key,None)
    return dict(smiles=label, molecule=dict(symbols=atom.get_chemical_symbols(),
        positions=atom.positions.tolist(), charge=inp.get('charge',0),
        multiplicity=inp.get('mult',1)), energy_hartree=float(data['energy'])*constants.EVtoHARTREE,
        harmonic_zpe_hartree=float(data['zpe']),
        frequencies_cm_inverse=list(map(float,data['frequencies'])),
        l2_theory=dict(calculator=row.get('calculator'), parameters=parameters),
        source_row_id=row.id, source_job=name)


def failure_inventory(worker, errors):
    """Bounded reads of recorded outcomes; an exit marker is not a live status."""
    failures = []
    for p in sorted(worker.glob('*.exitcode')):
        try:
            if int(p.read_text().strip()) == 0:
                continue
            result = read(p.with_suffix('.result.json'),errors)
            reason = result.get('error','')
            log = p.with_name(p.stem+'_sella.log')
            if log.exists():
                with log.open('rb') as f:
                    f.seek(max(0,log.stat().st_size-8192))
                    lines = f.read().decode(errors='replace').splitlines()
                reasons = [s for s in lines if 'Worker failure:' in s]
                if reasons:
                    reason = reasons[-1].split('Worker failure:',1)[1].strip()
            failures.append(dict(job=p.stem,source=str(p),recorded_reason=reason or 'native diagnostic required'))
        except (OSError, ValueError) as error:
            errors.append(dict(file=str(p),error=str(error)))
    return failures


def audit_case(case, *, use_database=False):
    errors, labels, workers, reactions, records, failures, soft = [], {}, [], [], {}, [], []
    local_aliases = {}
    registered = read(case/'pes_anl_reuse.json', errors)
    for inp_path in sorted(case.glob('[0-9]*/*.json')):
        if inp_path.name != inp_path.parent.name+'.json':
            continue
        worker = inp_path.parent
        inp = read(inp_path, errors)
        base = worker/inp.get('anl_automatic_directory','anl_automatic')
        plan = read(base/'network_request.json', errors)
        saved = read(worker/'accepted_l2_network.json', errors)
        network = plan.get('network', {})
        label = network.get('route_to_smiles', {}).get(worker.name)
        origin = 'network_manifest'
        if not label:
            try:
                label, origin = input_label(inp)
            except Exception as error:
                label, origin = '?', 'unknown'
                errors.append(dict(file=str(inp_path), error=str(error)))
        labels.update(network.get('route_to_smiles', {}))
        labels.setdefault(worker.name,label)
        labels.update(saved.get('network',{}).get('route_to_smiles',{}))
        aliases = read(base/'ts_reuse.json',errors).get('aliases',{})
        local_aliases.update(aliases)
        failures.extend(dict(worker=worker.name,**f) for f in failure_inventory(worker,errors))
        rows = []
        monitor = worker/'kinbot_monitor.out'
        try:
            before = monitor.read_bytes() if monitor.exists() else b''
        except OSError as error:
            errors.append(dict(file=str(monitor),error=str(error)))
            before = b''
        for line in before.decode(errors='replace').splitlines():
            a = line.split()
            if len(a)>=4 and a[0]=='-1':
                rows.append(dict(name=a[2],reactant=worker.name,products=sorted(a[3:]),
                    transition_state=None if '_hom_sci_' in a[2] else a[2],
                    source='live_monitor', worker=worker.name))
        # Summary/plan alone does not overwrite the current acceptance status.
        by_name = {r['name']: r for r in network.get('reactions',[])}
        snap_rows = {r['name']:r for r in saved.get('network',{}).get('reactions',[])}
        counts = Counter()
        for row in rows:
            name = row['name']
            prior = by_name.get(name)
            snap = snap_rows.get(name)
            if snap and endpoints(row)==endpoints(snap):
                row.update(path_class=snap.get('path_class'))
                record = saved['network'].get('transition_states',{}).get(name)
                origin_record = 'accepted_l2_snapshot'
            elif prior and endpoints(row)==endpoints(prior):
                row.update(path_class=prior.get('path_class'))
                record = network.get('transition_states',{}).get(name)
                origin_record = 'network_manifest'
            else:
                record = None
                origin_record = 'unavailable'
            if record is None and row['transition_state'] and use_database and (worker/'kinbot.db').is_file():
                job = name+('_high' if inp.get('high_level',1) else '')
                try:
                    record = database_record(worker/'kinbot.db',job,inp,label)
                    origin_record = 'latest_successful_cache_row_unverified' if record else 'unavailable'
                except Exception as error:
                    errors.append(dict(file=str(worker/'kinbot.db'),name=name,error=str(error)))
            if record:
                records[name] = deepcopy(record)
                row['record_source'] = origin_record
                negative = [v for v in record.get('frequencies_cm_inverse',[]) if v<0]
                threshold = inp.get('imagfreq_threshold',50.)
                if negative and min(negative)>=-threshold:
                    soft.append(dict(name=name,imaginary_cm_inverse=negative,
                        threshold_cm_inverse=threshold,record_source=origin_record,
                        action='IRC/mode review; not automatically rejected'))
            counts[origin_record if row['transition_state'] else 'barrierless'] += 1
        try:
            if monitor.exists() and before != monitor.read_bytes():
                errors.append(dict(file=str(monitor),error='acceptance_changed_during_read; repeat audit'))
        except OSError as error:
            errors.append(dict(file=str(monitor),error=str(error)))
        workers.append(dict(name=worker.name,label=label,label_source=origin,
            accepted_search_routes=len(rows), endpoint_groups=len({endpoints(r) for r in rows}),
            network_materialized=bool(network), snapshot_available=bool(saved),
            record_sources=dict(counts), local_aliases=aliases))
        reactions.extend(rows)
    groups = defaultdict(list)
    for row in reactions:
        groups[endpoints(row)].append(row)
    comparisons, decisions = [], Counter()
    routes = {r['name'] for r in reactions}
    alias_map = dict(local_aliases)
    for a,b in registered.get('aliases',{}).items():
        a,b = a.rsplit('/',1)[-1],b.rsplit('/',1)[-1]
        if a in routes and b in routes:
            alias_map[a] = b
    def canonical(name):
        seen = set()
        while name in alias_map and name not in seen:
            seen.add(name)
            name = alias_map[name]
        return name
    for group in groups.values():
        for a,b in combinations(group,2):
            if not a['transition_state'] or not b['transition_state']:
                continue
            pair = dict(left=a['name'],right=b['name'], reverse=a['reactant']!=b['reactant'])
            if a['name'] not in records or b['name'] not in records:
                answer = dict(equivalent=False,reason='accepted_L2_record_unavailable')
            elif ('cache' in a.get('record_source','')) != ('cache' in b.get('record_source','')):
                answer = dict(equivalent=False,reason='different_evidence_sources; requires_validation')
            else:
                p = dict(network=dict(reactions=[a,b], transition_states=records))
                try:
                    answer = ts_reuse.comparison(p,a['name'],b['name'],ts_reuse.options({}))
                except (KeyError, ValueError, TypeError) as error:
                    answer = dict(equivalent=False,reason='incomplete_comparison_evidence',error=str(error))
            pair.update(answer)
            pair['already_mapped'] = canonical(a['name'])==canonical(b['name'])
            if a['name'] in records and b['name'] in records:
                try:
                    pair['diagnostic_differences'] = numerical_differences(records[a['name']],records[b['name']])
                except (KeyError,ValueError,TypeError) as error:
                    pair['diagnostic_error'] = str(error)
            pair['reusable'] = False  # audit never grants ownership or scientific acceptance
            pair['record_sources'] = [a.get('record_source'),b.get('record_source')]
            comparisons.append(pair)
            decisions[answer['reason']] += 1
    return dict(case=case.name,labels=labels,workers=workers,reactions=reactions,
        endpoint_pairs=[dict(endpoints=k,routes=[r['name'] for r in v]) for k,v in groups.items()],
        accepted_search_routes=len(reactions), endpoint_pair_count=len(groups),
        local_aliases=local_aliases,pes_aliases=registered.get('aliases',{}),
        owners=len(registered.get('owners',[])),comparison_counts=dict(decisions),
        comparisons=comparisons,soft_saddles=soft,recorded_worker_failures=failures,errors=errors)


def summary(report):
    out = ['Read-only PES audit | '+report['captured_utc']]
    for c in report['cases']:
        out += ['', c['case']+': '+str(c['accepted_search_routes'])+' accepted search routes; '+
                str(c['endpoint_pair_count'])+' unordered endpoint pairs (not unique TS count)',
                f"Reuse mappings: local={len(c['local_aliases'])}; across wells={len(c['pes_aliases'])}; owners={c['owners']}",
                'Pair comparison reasons: '+json.dumps(c['comparison_counts'],sort_keys=True)]
        for w in c['workers'][:12]:
            out.append(f"  {w['name']} {w['label']}: routes={w['accepted_search_routes']}; "
                       f"groups={w['endpoint_groups']}; evidence={w['record_sources']}")
        if len(c['workers'])>12:
            out.append(f"  {len(c['workers'])-12} additional workers in report.json")
        candidates = [p for p in c['comparisons'] if p['equivalent'] and not p['already_mapped']]
        for p in candidates[:6]:
            out.append('  MATCH CANDIDATE: '+p['left']+' <-> '+p['right'])
        near = sorted((p for p in c['comparisons'] if not p['equivalent'] and 'diagnostic_differences' in p),
                      key=lambda p:p['diagnostic_differences'].get('minimum_proper_rotation_rmsd_angstrom',1.e99))
        for p in near[:3]:
            out.append('  UNMERGED: '+p['left']+' <-> '+p['right']+' '+
                       json.dumps(p['diagnostic_differences'],sort_keys=True))
        out.append(f"Read/consistency errors={len(c['errors'])}; match candidates={len(candidates)}")
        causes = Counter(f['recorded_reason'] for f in c.get('recorded_worker_failures',[]))
        out.append(f"Recorded nonzero worker exits={sum(causes.values())} (may precede retries)")
        for reason,count in causes.most_common(3):
            out.append(f'  {count}x '+reason[:220])
        for row in c.get('soft_saddles',[])[:3]:
            out.append('  SOFT SADDLE REVIEW: '+row['name']+' '+str(row['imaginary_cm_inverse'])+' cm^-1')
    out += ['Endpoint grouping does not merge different saddles, stereoisomers, or complexes.',
            'Cache matches require accepted geometry/pathway verification before cancellation or reuse.',
            'No jobs, run files, or task states changed.']
    return '\n'.join(out)+'\n'


def dot(case):
    """Topology only; no invented barriers or relative-energy ordering."""
    q = json.dumps
    lines = ['graph PES {', 'graph [overlap=false];', 'node [shape=box];']
    for pair in case['endpoint_pairs']:
        for end in pair['endpoints']:
            key = '+'.join(end)
            label = ' + '.join(case['labels'].get(x,x) for x in end)
            lines.append(f'{q(key)} [label={q(label)}];')
        a,b = ('+'.join(e) for e in pair['endpoints'])
        lines.append(f'{q(a)} -- {q(b)} [label={q(str(len(pair["routes"]))+" searches")}];')
    lines.append('}')
    return '\n'.join(lines)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('--database',action='store_true',help='Read latest ASE cache rows for missing L2 snapshots (not reuse proof).')
    parser.add_argument('--output',type=Path,required=True,help='New report directory outside the run')
    args = parser.parse_args()
    root, output = args.root.resolve(strict=True), args.output.resolve()
    if output.is_relative_to(root):
        parser.error('Report output must be outside the run directory.')
    output.mkdir(parents=True,exist_ok=True)
    if any(output.iterdir()):
        parser.error('Report directory must be empty; existing reports are preserved.')
    cases = [root] if (root/'input.json').exists() else sorted(p for p in root.iterdir()
        if p.is_dir() and (p/'input.json').exists() and not p.name.startswith('.'))
    report = dict(schema=1,root=str(root),captured_utc=datetime.now(timezone.utc).isoformat(),
        cases=[audit_case(c,use_database=args.database) for c in cases])
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    (output/'summary.txt').write_text(summary(report))
    for c in report['cases']:
        (output/(c['case']+'.dot')).write_text(dot(c))
    print(summary(report),end='')
    print('Full report and topology graphs: '+str(output))


if __name__ == '__main__':
    main()
