"""Live retry gates for the ordinary reaction and conformer search phases."""
import json
from pathlib import Path

from kinbot import live_control as control
from kinbot.execution import atomic_json


def _backend(qc, job):
    return qc._backend_for(job) if hasattr(qc, '_backend_for') else qc


def reactions(generator):
    """Rebuild a rejected reaction while its search still owns the network."""
    species, qc = generator.species, generator.qc
    objects = {obj.instance_name:(i,obj) for i,obj in enumerate(species.reac_obj)}
    _offline_revisit(generator, objects)
    if control.find() is None:
        return
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


def _offline_revisit(generator, objects):
    """A refresh grants one additional attempt, without resetting auto budgets."""
    marker = Path('.kinbot_revisit.json')
    if not marker.exists():
        return
    m = json.loads(marker.read_text())
    for name in m['rejected']:
        if name in m['attempted'] or name not in objects:
            continue
        i, obj = objects[name]
        if generator.species.reac_ts_done[i] != -999:
            continue  # Normal replay must first establish the rejected result.
        if getattr(obj, 'stereochemical_rejection', None):
            m.setdefault('blocked', {})[name] = obj.stereochemical_rejection
            m['attempted'].append(name)
            atomic_json(marker,m)
            continue
        new = type(obj)(generator.species, generator.qc, generator.par, obj.instance, obj.instance_name)
        new.prod_done, new.valid_prod = 0, []
        token = m['revision']+'_'+name
        path = Path('.kinbot_controls/revisits')/(token+'.json')
        r = json.loads(path.read_text()) if path.exists() else dict(id=token,status='queued')
        _retire_reaction(generator.qc,name,set(objects),path,r)
        generator.species.reac_obj[i] = new
        generator.species.reac_ts_done[i],generator.species.reac_step[i] = 0,0
        m['attempted'].append(name)
        atomic_json(marker,m)
        control.acknowledge(path,r,'applied')


def pending_revisit(generator):
    marker=Path('.kinbot_revisit.json')
    if not marker.exists():
        return False
    m=json.loads(marker.read_text())
    rejected={o.instance_name for o,s in zip(generator.species.reac_obj,generator.species.reac_ts_done) if s==-999}
    return bool(rejected.intersection(set(m['rejected'])-set(m['attempted'])))


def revise_network(generator):
    """Apply operator search-policy changes only after global branch quiescence."""
    from kinbot.execution import file_sha256
    for path, r in control.requests(Path.cwd(),kinds={'network'}):
        if Path(r['target'])!=Path.cwd():
            continue
        inp = Path(generator.inp).resolve()
        try:
            if r['status']!='applying' and file_sha256(inp)!=r['expected']['input_sha256']:
                raise ValueError('Worker input changed after network revisit was requested.')
            options=r['options']
            species=generator.species
            models=options.get('vrc_tst_surface_models')
            if models is not None:
                from kinbot.vrc_multistate import extend_surface_model, normalize_inline_model_paths, load_surface_model
                from kinbot.species_routing import routing_name
                normalize_inline_model_paths(options,inp.parent)
                existing=generator.par.get('vrc_tst_surface_models',{})
                if set(existing)-set(models):
                    raise ValueError('A live model extension cannot remove existing species models.')
                current_route=routing_name(species)
                for selector, model in models.items():
                    prior=existing.get(selector)
                    # Adding an override of the current default must preserve
                    # its accepted channels, just as updating an explicit map.
                    if (prior is None and selector in (current_route, species.smiles)
                            and generator.par.get('vrc_tst_surface_model_scope') in (None,'',current_route)):
                        loaded=load_surface_model(generator.par)
                        prior=loaded[0] if loaded else None
                    extend_surface_model(prior,model)
            replay=options.get('vdW_detection') is not None or options.get('rejected',False)
            rejected=[o.instance_name for o,status in zip(species.reac_obj,species.reac_ts_done)
                      if status==-999] if options.get('rejected') else []
            copies=([type(o)(species,generator.qc,generator.par,o.instance,o.instance_name)
                     for o in species.reac_obj] if replay else [])
            for obj in copies:
                obj.prod_done,obj.valid_prod=0,[]
            control.acknowledge(path,r,'applying')
            saved=json.loads(inp.read_text())
            override_path=Path('.kinbot_search_options.json')
            overrides=json.loads(override_path.read_text()) if override_path.exists() else {}
            if options.get('vdW_detection') is not None:
                saved['vdW_detection']=options['vdW_detection']
                overrides['vdW_detection']=options['vdW_detection']
            if models is not None:
                saved['vrc_tst_surface_models']=models
                overrides['vrc_tst_surface_models']=models
            atomic_json(inp,saved)
            atomic_json(override_path,overrides)
            generator.par.update(overrides)
            if replay:
                atomic_json(Path('.kinbot_revisit.json'),dict(schema=1,revision=r['id'],rejected=rejected,attempted=[]))
                species.reac_obj=copies
                species.reac_ts_done=[0]*len(copies)
                species.reac_step=[0]*len(copies)
            control.acknowledge(path,r,'applied',detail=(
                'Replaying cached accepted searches and retrying selected failures.' if replay
                else 'Reviewed VRC models accepted; existing searches retained.'))
        except OSError:
            # Keep the applying receipt for reconciliation. Do not continue a
            # partially persisted scientific edit as though it were rejected.
            raise
        except (RuntimeError,ValueError) as error:
            control.acknowledge(path,r,'blocked',error=str(error))
