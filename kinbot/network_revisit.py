"""Stopped-controller network/TS migration, with exact ownership and archives."""
from copy import deepcopy
import json
import math
from pathlib import Path
import uuid

from kinbot.execution import atomic_json, file_sha256
from kinbot.anl import ts_reuse


def settings(vdw=None, rejected=False, deduplicate=False, barrier_l2=None):
    if vdw is not None and (isinstance(vdw, bool) or not math.isfinite(vdw) or vdw < 0):
        raise ValueError('vdW_detection must be finite and nonnegative.')
    if barrier_l2 is not None and (isinstance(barrier_l2,bool) or not math.isfinite(barrier_l2) or barrier_l2<=0):
        raise ValueError('L2 barrier ceiling must be positive and finite.')
    return dict(vdW_detection=vdw, rejected=bool(rejected), deduplicate_ts=bool(deduplicate), barrier_threshold_L2=barrier_l2)


def retirement(root, live, graph, canonical_graph, source):
    from kinbot.anl import dispatch
    cg, spec, state = dispatch._load(canonical_graph)
    if any(spec['molecule'].get(k)!=source.get(k) for k in ('symbols','positions','charge','multiplicity')):
        raise RuntimeError(f'{cg}: canonical workflow geometry differs from the accepted L2 record.')
    for task in spec['tasks']:
        e = state['tasks'].get(task['id'],{})
        if e.get('status')=='complete':
            dispatch._verify_stage_files(cg,task,e)
            execution = json.loads((cg/'tasks'/task['id']/'execution.json').read_text())
            if dispatch._verify_execution(cg,task,e,execution)!=e.get('final_geometry_sha256'):
                raise RuntimeError(f'{cg}: canonical geometry verification failed.')
    _, spec, state = dispatch._load(graph)
    selected=[]
    for task in spec['tasks']:
        e=state['tasks'].get(task['id'],{})
        if e.get('status')=='submitting':
            raise RuntimeError(f'{graph}: ambiguous submission; reconcile before retiring.')
        job=live.get(str(graph/'tasks'/task['id']))
        if job:
            if str(e.get('job_id'))!=str(job['job_id']):
                raise RuntimeError(f'{graph}: live allocation owner differs.')
            selected.append(job)
    return dict(graph=str(graph),workflow_sha256=file_sha256(graph/'workflow.json')),selected


def plan(root, jobs, options):
    from kinbot import run_salvage
    from kinbot.anl import dispatch
    prepared, inputs = run_salvage._cases(root)
    edits, workers, aliases, selected = [], [], [], []
    live = {str(Path(j['directory']).resolve()): j for j in jobs}
    for case in inputs:
        from kinbot.barrier_policy import for_case
        reference = for_case(case.parent,options['barrier_threshold_L2']) if options.get('barrier_threshold_L2') is not None else None
        paths = [case]
        if (case.parent/'resolved_parameters.json').exists():
            paths.append(case.parent/'resolved_parameters.json')
        paths += sorted(case.parent.glob('[0-9]*/*.json'))
        paths += sorted(case.parent.glob('[0-9]*/.kinbot_search_options.json'))
        paths = [p for p in paths if p in (case,case.parent/'resolved_parameters.json')
                 or p.name==p.parent.name+'.json' or p.name=='.kinbot_search_options.json']
        for p in paths:
            previous = json.loads(p.read_text())
            revised = deepcopy(previous)
            if p.name=='.kinbot_search_options.json':
                if options['vdW_detection'] is not None:
                    revised['vdW_detection']=options['vdW_detection']
                if revised!=previous:
                    edits.append(dict(path=str(p.relative_to(root)),before=file_sha256(p),after=revised))
                continue
            if reference is not None:
                revised.update(pes_barrier_reference=reference, barrier_threshold=None,
                               barrier_threshold_L2=reference['ceiling_kcal_mol'])
            revised['pes_phase_policy'] = 'exploration_first'
            if options['deduplicate_ts']:
                revised['anl_ts_reuse'] = {}
            if options['vdW_detection'] is not None:
                revised['vdW_detection'] = options['vdW_detection']
            # Preserve successful L1 work while permitting new reaction work.
            revised['profiled_l1_reuse_only'] = 0
            if revised != previous:
                edits.append(dict(path=str(p.relative_to(root)),before=file_sha256(p),after=revised))
            if p.name!=p.parent.name+'.json':
                continue
            worker = p.parent
            monitor = worker/'kinbot_monitor.out'
            rows = [line.split() for line in monitor.read_text().splitlines()] if monitor.exists() else []
            rejected = [r[2] for r in rows if len(r)>=3 and r[0]=='-999'] if options['rejected'] else []
            workers.append(dict(path=str(worker.relative_to(root)), rejected=rejected,
                                automatic_directory=previous.get('anl_automatic_directory','anl_automatic')))
            base = worker/previous.get('anl_automatic_directory','anl_automatic')
            request = base/'network_request.json'
            if not options['deduplicate_ts'] or not request.exists():
                continue
            network = json.loads(request.read_text())
            old = json.loads((base/'ts_reuse.json').read_text()) if (base/'ts_reuse.json').exists() else None
            mapping = ts_reuse.discover(network, base, ts_reuse.options({}), old, ts_reuse.protected(base))
            remove = []
            for alias, canonical in mapping['aliases'].items():
                graph = ts_reuse.directory(base,alias)
                if not graph.exists():
                    continue
                canonical_graph = ts_reuse.directory(base,canonical)
                item, owned = retirement(root,live,graph,canonical_graph,
                    network['network']['transition_states'][canonical]['molecule'])
                selected.extend(owned)
                remove.append(item)
            aliases.append(dict(base=str(base),network_sha256=file_sha256(request),mapping=mapping,remove=remove))
    external=[]
    if options['deduplicate_ts']:
        from kinbot.anl import pes_ts_reuse
        local_maps={e['base']:e['mapping'] for e in aliases}
        for case in inputs:
            scope=case.parent
            if not (scope/'pes_exploration.json').exists():
                continue
            mapping,items=pes_ts_reuse.discover(scope,ts_reuse.options({}),local_maps)
            remove=[]
            for key,other in mapping['aliases'].items():
                item,target=items[key],items[other]
                graph=Path(item['graph'])
                if not graph.exists():
                    continue
                canonical=Path(target['graph'])
                retired,owned=retirement(root,live,graph,canonical,target['record']['molecule'])
                remove.append(retired);selected.extend(owned)
            external.append(dict(case=str(scope),mapping=mapping,remove=remove))
    return dict(settings=options, inputs=edits, workers=workers, aliases=aliases,
                external_ts=external, selected_jobs=selected, prepared_before=prepared,
                prepared_sha256=file_sha256(root/'prepared.json'), revision=uuid.uuid4().hex)


def apply(root, record, archive, journal, marker):
    from kinbot import vrc_restart, run_refresh
    root, archive = Path(root), Path(archive)
    # Same fail-closed cancellation protocol, restricted to verified TS aliases.
    vrc_restart.cancel_scans(record,journal,marker)
    remaining = run_refresh._owned_jobs(root)
    for entry in record['aliases']:
        base = Path(entry['base'])
        if file_sha256(base/'network_request.json') != entry['network_sha256']:
            raise RuntimeError('Accepted network changed during TS deduplication.')
        for item in entry['remove']:
            graph = Path(item['graph'])
            saved = archive/'duplicate_ts'/graph.relative_to(root)
            if any(Path(j['directory']).is_relative_to(graph) for j in remaining):
                raise RuntimeError(f'{graph}: duplicate allocation is still present.')
            if graph.exists():
                if saved.exists() or file_sha256(graph/'workflow.json')!=item['workflow_sha256']:
                    raise RuntimeError(f'{graph}: archive/input collision.')
                saved.parent.mkdir(parents=True,exist_ok=True)
                graph.rename(saved)
            elif not saved.exists():
                raise RuntimeError(f'{graph}: duplicate and archive are both missing.')
        atomic_json(base/'ts_reuse.json',entry['mapping'])
    for entry in record.get('external_ts',[]):
        for item in entry['remove']:
            graph=Path(item['graph']);saved=archive/'duplicate_ts'/graph.relative_to(root)
            if any(Path(j['directory']).is_relative_to(graph) for j in remaining):
                raise RuntimeError(f'{graph}: duplicate allocation still present.')
            if graph.exists():
                if saved.exists() or file_sha256(graph/'workflow.json')!=item['workflow_sha256']:
                    raise RuntimeError(f'{graph}: archive/input collision.')
                saved.parent.mkdir(parents=True,exist_ok=True);graph.rename(saved)
            elif not saved.exists():
                raise RuntimeError(f'{graph}: duplicate and archive both missing.')
        atomic_json(Path(entry['case'])/'pes_anl_reuse.json',entry['mapping'])
    # Reuse the exact-hash, idempotent input publisher; it edits no QC graph.
    payload = {**record,'scans':[]}
    vrc_restart.apply_inputs(root,payload,archive/'inputs', record_name='network_restart.json')
    from kinbot.pes_phase import authorize_revision
    for entry in record['workers']:
        worker = root/entry['path']
        authorize_revision(worker,record['revision'])
        pending = worker/'.kinbot_revisit.json'
        if pending.exists() and json.loads(pending.read_text()).get('revision')==record['revision']:
            continue
        atomic_json(pending,dict(schema=1, revision=record['revision'],
            rejected=entry['rejected'], attempted=[], reassess=False))
    record['archive'] = str(archive)
    atomic_json(marker,journal)
