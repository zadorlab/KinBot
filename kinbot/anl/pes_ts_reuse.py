"""Immutable composite graph ownership within one PES (never across PESs).

Registration is serialized independently of dispatcher locks. An alias finishes
its private work first, then waits for or adopts the canonical graph under its
existing drive lock. No driver holds private graph locks while waiting for a
borrowed graph, so reverse discoveries cannot create a lock cycle.
"""
from contextlib import contextmanager
from copy import deepcopy
import fcntl
import hashlib
import json
from pathlib import Path
import time

from kinbot.execution import atomic_json
from kinbot.anl import ts_reuse

NAME = 'pes_anl_reuse.json'


def scope(base):
    base = Path(base).resolve()
    for worker in base.parents:
        if ((worker/(worker.name+'.json')).is_file()
                and (worker.parent/'pes_exploration.json').is_file()):
            return worker.parent
    return None


@contextmanager
def locked(case):
    with (case/'.pes_ts_reuse.lock').open('a+') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def scientific(plan):
    resources = {k:v for k,v in plan['resources'].items()
                 if k not in ('partition', 'cores_per_job', 'global_max_jobs')}
    return dict(resources=resources, l2_theory=plan.get('l2_theory'),
                requested_method=plan.get('requested_method'))


def signature(record):
    value=deepcopy(record)
    value.pop('roles',None)  # A product can later also be an explored well.
    return ts_reuse.digest(value)


def catalog(case, local_maps=None):
    result = {}
    for inp in sorted(case.glob('[0-9]*/*.json')):
        if inp.name != inp.parent.name+'.json':
            continue
        parameters = json.loads(inp.read_text())
        base = (inp.parent/parameters.get('anl_automatic_directory','anl_automatic')).resolve()
        if not base.is_relative_to(case.resolve()):
            raise ValueError('TS reuse directory escapes its PES.')
        p = base/'network_request.json'
        if not p.is_file():
            continue
        plan = json.loads(p.read_text())
        mapping = (local_maps or {}).get(str(base))
        if mapping is None:
            mapping = json.loads((base/'ts_reuse.json').read_text()) if (base/'ts_reuse.json').exists() else {}
        for name, record in plan['network']['transition_states'].items():
            if name in mapping.get('aliases',{}):
                continue
            key = str(base.relative_to(case))+'/'+name
            result[key] = dict(base=str(base.relative_to(case)), name=name, kind='ts',
                              record=record, entry=record, plan=plan,
                              graph=str(ts_reuse.directory(base,name)))
        from kinbot.anl.automatic import _state_key, _safe_state_directory, _model_record
        for entry in plan['cbh']['calculated_states']:
            state=_state_key(entry['smiles'],entry['charge'],entry['multiplicity'])
            ident=_safe_state_directory(entry['smiles'],entry['charge'],entry['multiplicity'])
            record=entry.get('network_record')
            if record is None:
                model=base/'reference_l2'/ident/'thermochemistry_model.json'
                if not model.is_file():
                    continue
                record=_model_record(model,entry,plan['l2_theory'])
            name='stable:'+state
            revision=entry.get('graph_revision','')
            key=str(base.relative_to(case))+'/'+name+'/'+revision
            graph=base/'states'/(ident+('_'+revision if revision else ''))/'composite'
            result[key]=dict(base=str(base.relative_to(case)),name=name,kind='stable',
                             record=record,entry=entry,plan=plan,graph=str(graph))
        for route, entry in plan['network'].get('product_complexes', {}).items():
            name = 'complex:'+route
            graph = base/'states'/('complex_'+hashlib.sha256(route.encode()).hexdigest()[:16])/'composite'
            result[str(base.relative_to(case))+'/'+name] = dict(
                base=str(base.relative_to(case)), name=name, kind='complex',
                record=entry['record'], entry=entry['record'], products=entry['products'],
                plan=plan, graph=str(graph))
    return result


def compare(a, b, limits):
    if a['kind']!=b['kind'] or scientific(a['plan'])!=scientific(b['plan']):
        return None
    if a['kind'] in ('stable', 'complex'):
        if any(a['entry'].get(k)!=b['entry'].get(k) for k in ('smiles','charge','multiplicity','highest_tier','fallback_tiers')):
            return None
        if a['kind']=='complex' and a['products']!=b['products']:
            return None
        merged={'network':{'transition_states':{'a':a['record'],'b':b['record']},
            'reactions':[{'transition_state':n,'reactant':'same','products':['same']} for n in ('a','b')]}}
        # Monatomic minima have no vibrational modes; exact atom/state/theory
        # and energy identity suffice for their geometry operation.
        if len(a['record']['molecule']['symbols'])==1:
            merged=deepcopy(merged)
            for r in merged['network']['transition_states'].values():
                r['frequencies_cm_inverse']=[0.]
        return ts_reuse.equivalent(merged,'a','b',limits,stationary_point='minimum')
    # TS SMILES is a reactant connectivity anchor, not the saddle's identity.
    merged = deepcopy(a['plan'])
    merged['network']['transition_states'] = {'a':a['record'],'b':b['record']}
    merged['network']['reactions'] = []
    for label, item in [('a',a),('b',b)]:
        rows = [r for r in item['plan']['network']['reactions'] if r.get('transition_state')==item['name']]
        if len(rows)!=1:
            return None
        merged['network']['reactions'].append(dict(rows[0],transition_state=label))
    # A vdW endpoint has a discovery-specific route. Normalize it only after
    # proving the two accepted complex geometries/states are equivalent.
    left, right = (i['plan']['network'].get('product_complexes', {}) for i in (a,b))
    for route, entry in left.items():
        for other, candidate in right.items():
            x=dict(a,kind='complex',entry=entry['record'],record=entry['record'],products=entry['products'])
            y=dict(b,kind='complex',entry=candidate['record'],record=candidate['record'],products=candidate['products'])
            if compare(x,y,limits) is None:
                continue
            row=merged['network']['reactions'][1]
            if row['reactant']==other:
                row['reactant']=route
            row['products']=[route if p==other else p for p in row['products']]
    return ts_reuse.equivalent(merged,'a','b',limits)


def history(case, result, current):
    """Keep old graph revisions usable when a reference becomes a network well."""
    items = dict(current)
    for key, expected in result['records'].items():
        item = items.get(key, result.get('catalog', {}).get(key))
        if item is None or signature(item['record']) != expected:
            raise ValueError(f'{key}: PES ANL ownership provenance changed.')
        graph = Path(item['graph']).resolve()
        if not graph.is_relative_to(case.resolve()):
            raise ValueError('PES ANL graph escapes its owning PES.')
        items[key] = item
    return items


def discover(case, limits, local_maps=None, *, reconcile_owners=False):
    p = case/NAME
    result = json.loads(p.read_text()) if p.exists() else dict(schema=1,aliases={},records={},evidence={})
    items = history(case, result, catalog(case,local_maps))
    def rank(key):
        item=items[key]; graph=Path(item['graph'])
        state=json.loads((graph/'state.json').read_text()) if (graph/'state.json').exists() else {}
        entries=list(state.get('tasks',{}).values())
        return (-sum(e.get('status')=='complete' for e in entries),
                -sum(e.get('status')=='submitted' for e in entries),-int(graph.exists()),key)
    # Reserve even unique owners before releasing the registry lock. Otherwise
    # two simultaneously preparing wells can both see an uncreated graph and
    # assign ownership in opposite directions.
    canonical=sorted(set(result.get('owners',result['aliases'].values())),key=rank)
    if reconcile_owners:
        # Only the stopped-controller repair may reconsider reserved owners.
        # Normal registration must keep reservations immutable to avoid two
        # live planners assigning ownership in opposite directions.
        survivors = []
        for key in canonical:
            match = next(((other, proof) for other in survivors
                          if (proof := compare(items[key], items[other], limits)) is not None), None)
            if match is None:
                survivors.append(key)
                continue
            other, proof = match
            # Re-prove every redirected alias against the final owner. Close
            # numerical matches are not necessarily transitive.
            children = [alias for alias, owner in result['aliases'].items() if owner == key]
            proofs = {alias: compare(items[alias], items[other], limits) for alias in children}
            if any(value is None for value in proofs.values()):
                survivors.append(key)
                continue
            for alias, evidence in {key: proof, **proofs}.items():
                result['aliases'][alias] = other
                result['evidence'][alias] = dict(canonical=other, measurements=evidence,
                    limits=limits, scientific=scientific(items[alias]['plan']))
        canonical = survivors
    for key in sorted(items,key=rank):
        if key in result['aliases'] or key in canonical:
            continue
        match=next(((other,evidence) for other in canonical
                    if (evidence:=compare(items[key],items[other],limits)) is not None),None)
        if match:
            other,evidence=match
            result['aliases'][key]=other
            result['evidence'][key]=dict(canonical=other,measurements=evidence,limits=limits,
                                         scientific=scientific(items[key]['plan']))
        else:
            canonical.append(key)
    result['owners']=canonical
    result['records']={key:signature(items[key]['record']) for key in
                       set(result['aliases'])|set(canonical)}
    result['catalog']={key:deepcopy(items[key]) for key in result['records']}
    return result,items


def install(base, value):
    case=scope(base)
    if case is None:
        return
    with locked(case):
        p=case/NAME
        previous=json.loads(p.read_text()) if p.exists() else {'aliases':{}}
        result,items=discover(case,ts_reuse.options(value))
        for key in set(result['aliases'])-set(previous['aliases']):
            item=items[key]
            if Path(item['graph']).exists():
                raise ValueError('Equivalent ANL graph already exists; use --deduplicate-anl with stopped controllers.')
        atomic_json(p,result)


def load(base):
    case=scope(base)
    if case is None or not (case/NAME).exists():
        return {}
    result=json.loads((case/NAME).read_text()); current=catalog(case)
    items=history(case,result,current)
    aliases=result['aliases']
    if any(k==v or v in aliases for k,v in aliases.items()):
        raise ValueError('PES TS ownership must be flat.')
    answer={}
    for key,other in aliases.items():
        if key not in current:
            continue
        for k in (key,other):
            if k not in items or signature(items[k]['record'])!=result['records'][k]:
                raise ValueError('PES TS reuse provenance changed.')
        if str((case/items[key]['base']).resolve())!=str(Path(base).resolve()):
            continue
        evidence=result['evidence'][key]
        if (scientific(items[key]['plan'])!=evidence['scientific'] or
                compare(items[key],items[other],evidence['limits']) is None):
            raise ValueError('PES TS reuse scientific settings changed.')
        target=items[other]
        answer[items[key]['name']] = dict(target,
            graph=Path(target['graph']))
    return answer


def finish(target, interval=20.):
    """Reuse a whole canonical graph; if its owner stopped, safely drive it."""
    from kinbot.anl import automatic
    from kinbot.pes_phase import checkpoint
    graph=target['graph']
    while True:
        checkpoint()
        # Callers have released all their private graph locks. Ownership of
        # this graph's drive lock is sufficient; waiting for the peer's whole
        # network could wait forever after selection omitted this saddle or
        # an unrelated CBH reference blocked that peer.
        if target['kind']=='ts':
            automatic._prepare_ts_graph(target['plan'],target['record'],graph)
        else:
            automatic._prepare_stable_graph(target['plan'],target['entry'],target['record'],graph)
        with (graph/'drive.lock').open('a+') as lock:
            try:
                fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError:
                pass
            else:
                automatic.preflight_until_scheduler_ready(graph,skip_optional=True,interval=interval)
                automatic._drive_graphs_locked({'ts:'+target['name']:graph},
                    target['plan']['resources']['global_max_jobs'],interval)
                return graph
        time.sleep(interval)
