"""Live retry gates for the ordinary reaction and conformer search phases."""
import json
from pathlib import Path

from kinbot import live_control as control
from kinbot.execution import atomic_json


def _backend(qc, job):
    return qc._backend_for(job) if hasattr(qc, '_backend_for') else qc


def reactions(generator):
    """Rebuild a rejected reaction while its search still owns the network."""
    if control.find() is None:
        return
    species, qc = generator.species, generator.qc
    objects = {obj.instance_name:(i,obj) for i,obj in enumerate(species.reac_obj)}
    for name,(i,obj) in objects.items():
        control.search_scope(Path.cwd()/name,'reaction',species.reac_ts_done[i])
    for path,r in control.requests(Path.cwd(), kinds={'reaction'}):
        target = Path(r['target'])
        if target.parent != Path.cwd() or target.name not in objects:
            continue
        i,obj = objects[target.name]
        try:
            if r['action'] == 'release':
                control.acknowledge(path,r,'applied')
                continue
            if species.reac_ts_done[i] != -999:
                raise ValueError('Only a rejected reaction can be restarted; accepted or active reactions are retained.')
            if getattr(obj, 'stereochemical_rejection', None):
                raise ValueError('Unsupported stereochemistry cannot be waived by retrying.')
            if control.generation(target,'reaction') != r['expected'] and r['status'] != 'applying':
                raise ValueError('Reaction inputs changed after the request.')
            # Construct first so an unsupported reaction constructor cannot
            # leave the saved search partially retired.
            new = type(obj)(species, qc, generator.par, obj.instance, obj.instance_name)
            new.prod_done, new.valid_prod = 0, []
            _retire_reaction(qc, target.name, set(objects), path, r)
            species.reac_obj[i] = new
            species.reac_ts_done[i], species.reac_step[i] = 0, 0
            control.acknowledge(path,r,'applied', detail='Reaction search restarted from its original species geometry.')
        except (ValueError, RuntimeError, OSError) as error:
            control.acknowledge(path,r,'blocked' if isinstance(error,ValueError) else r['status'], error=str(error))


def _retire_reaction(qc, name, reaction_names, path, record):
    """Archive only this reaction's namespace; retain shared product wells."""
    from kinbot.qc import _slurm_job_active
    root = Path.cwd()
    folder = root/'.kinbot_controls/attempts'/record['id']
    journal = folder/'reaction.json'
    if journal.exists():
        plan = json.loads(journal.read_text())
    else:
        def owned(value):
            value = Path(value).name
            matches = [n for n in reaction_names if value==n or value.startswith(n+'_') or value.startswith(n+'.')]
            return bool(matches) and max(matches,key=len)==name
        db = qc.db
        rows = [row for row in db.select() if owned(str(row.get('name','')))]
        jobs = {row.name for row in rows}
        files = []
        for directory in (root,root/'perm',root/'hir',root/'conf',root/'perm/conf',root/'perm/hir'):
            if directory.is_dir():
                for p in directory.iterdir():
                    if owned(p.name):
                        if p.is_symlink() or p.is_dir():
                            raise ValueError(f'Reaction namespace contains a directory/symlink requiring review: {p}')
                        files.append(str(p.relative_to(root)))
                        if p.suffix=='.py':
                            jobs.add(str(p.relative_to(root).with_suffix('')))
        for job in jobs:
            backend = _backend(qc,job)
            ident = backend.job_ids.get(job)
            if ident and _slurm_job_active(ident) is not False:
                raise ValueError(f'{job}: allocation is live or scheduler status is unknown; retained.')
        plan = dict(files=files, rows=[row.id for row in rows], jobs=sorted(jobs))
        atomic_json(journal,plan)
    control.acknowledge(path,record,'applying',archive=str(folder))
    for relative in plan['files']:
        source, saved = root/relative, folder/'files'/relative
        if source.exists():
            if saved.exists():
                raise RuntimeError(f'Artifact reappeared during reaction retry: {source}')
            saved.parent.mkdir(parents=True,exist_ok=True)
            source.rename(saved)
    for row_id in plan['rows']:
        row = qc.db.get(id=row_id)
        if not any(x.data.get('manual_retry')==record['id'] for x in qc.db.select(name=row.name)):
            qc.db.write(row.toatoms(),name=row.name,data={'status':0,'manual_retry':record['id'],
                                                         'supersedes_row_id':row_id})
    for job in plan['jobs']:
        if hasattr(qc,'_clear_job_id'):
            qc._clear_job_id(job)
        backend = _backend(qc,job)
        backend.job_ids.pop(job,None)
        getattr(backend,'job_generation_tagged',{}).pop(job,None)


def close_reactions(generator):
    outcomes = [control.search_scope(Path.cwd()/o.instance_name,'reaction',status,close=True)
                for o,status in zip(generator.species.reac_obj,generator.species.reac_ts_done)]
    return all(outcomes)


def conformers(search, statuses, *, ring=False, close=False):
    if control.find() is None:
        return False
    from kinbot.qc_live import process
    held = False
    for i, outcome in enumerate(statuses):
        add = 'semi_emp_' if search.semi_emp else ''
        job = search.get_job_name(i,cyc=int(ring),add=add)
        target = Path.cwd()/job
        closed = control.search_scope(target,'conformer',outcome,close=close)
        if close:
            held |= not closed
            continue
        if outcome != 1:
            continue
        backend = _backend(search.qc,job)
        if process(backend,job,'error',kinds={'conformer'},allow_completed=True)=='running':
            statuses[i] = -1
        held |= control.holding(target,kind='conformer') and statuses[i]==1
    return held
