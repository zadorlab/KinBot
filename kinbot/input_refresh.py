"""Explicit file resubmission to an existing PES owner; never starts a driver.

The edited input is a proposal. New wells read the last accepted snapshot until
the owner publishes a journalled revision at its polling boundary.
"""
from copy import deepcopy
import fcntl
import json
from pathlib import Path
import uuid

from kinbot.execution import atomic_json, file_sha256

SNAPSHOT = '.kinbot_accepted_input.json'
PARAMETERS_SNAPSHOT = '.kinbot_accepted_parameters.json'
JOURNAL = '.kinbot_input_refresh.json'
FIELDS = {'vdW_detection', 'vrc_tst_surface_models'}


def capture(root, source):
    root, source = Path(root), Path(source).resolve()
    snapshot(root, source, json.loads(source.read_text()))
    return dict(input_file=str(source), input_sha256=file_sha256(source))


def snapshot(root, source, values):
    """Keep raw proposals and preserve the origin of external relative paths."""
    root, source = Path(root).resolve(), Path(source).resolve()
    atomic_json(root/SNAPSHOT, values)
    if source.parent != root:
        parsed = deepcopy(values)
        parsed.setdefault('input_file', str(source))
        atomic_json(root/PARAMETERS_SNAPSHOT, parsed)


def owner_directory(source, *, directory=None, require_live=True):
    from kinbot import pes_phase
    source = Path(source).resolve()
    candidates = {Path(directory).resolve()} if directory else {source.parent, Path.cwd().resolve()}
    found = []
    for root in candidates:
        marker = root/pes_phase.NAME
        if not marker.exists():
            continue
        m = json.loads(marker.read_text())
        if m.get('input_file') == str(source) and (not require_live or pes_phase.owner_alive(m)):
            found.append(root)
    if len(found) > 1:
        raise ValueError('Several PES owners use this input; specify --run-directory.')
    return found[0] if found else None


def accepted_source(source, *, directory=None):
    from kinbot import pes_phase
    source = Path(source).resolve()
    root = Path(directory).resolve() if directory else source.parent
    marker = root/pes_phase.NAME
    if marker.exists():
        m = json.loads(marker.read_text())
        if m.get('input_file') == str(source) and pes_phase.owner_alive(m):
            return root/(SNAPSHOT if source.parent == root else PARAMETERS_SNAPSHOT)
    return source


def validate_models(previous, current, base, *, scope=None):
    from kinbot.vrc_multistate import normalize_inline_model_paths, extend_surface_model, load_surface_model
    a, b = deepcopy(previous), deepcopy(current)
    normalize_inline_model_paths(a, base)
    normalize_inline_model_paths(b, base)
    old, new = a.get('vrc_tst_surface_models', {}), b.get('vrc_tst_surface_models', {})
    if set(old)-set(new):
        raise ValueError('A live update cannot remove reviewed species models.')
    for key, model in new.items():
        prior = old.get(key)
        if prior is None and key == scope:
            loaded = load_surface_model(dict(a, input_file=str(Path(base)/'input.json')))
            prior = loaded[0] if loaded else None
        extend_surface_model(prior, model)
        from kinbot.vrc_reference_review import validate_additions
        validate_additions(prior, model, b, base)
    return new


def submit(source, *, apply=False, directory=None, proposal=None, expected=None):
    """Submit an on-disk input, or atomically install a prepared proposal.

    Programmatic editors provide both proposal and expected. Validation and
    the pending-revision check then precede any input write under the same
    control/phase locks used by publication. CLI file resubmission is unchanged.
    """
    from kinbot import pes_phase as phase, live_control as control
    source = Path(source).resolve(strict=True)
    root = owner_directory(source, directory=directory)
    if root is None:
        raise ValueError('No live PES owner for this input; no driver was launched.')
    control_path = control.find(root)
    if control_path is None:
        raise ValueError('Enable live controls before resubmitting an input.')
    with control.locked(control_path.parent), phase.locked(root):
        m = json.loads((root/phase.NAME).read_text())
        if m.get('input_file') != str(source) or not phase.owner_alive(m):
            raise ValueError('No live PES owner for this input; no driver was launched.')
        current = json.loads(source.read_text())
        if proposal is not None and (expected is None or current != expected):
            raise ValueError('Case input changed during review preparation; no input was overwritten.')
        proposed = deepcopy(proposal) if proposal is not None else current
        prior = json.loads((root/SNAPSHOT).read_text())
        old = json.loads((root/JOURNAL).read_text()) if (root/JOURNAL).exists() else None
        if (old and old['token']==m['token'] and proposed == current
                and old['source_sha256']==file_sha256(source)):
            return receipt(root)
        if old and old['token']==m['token'] and receipt(root)['status'] not in ('applied', 'blocked'):
            raise ValueError('An input revision is still pending; inspect its receipt first.')
        changed = {k for k in prior.keys() | proposed.keys() if prior.get(k)!=proposed.get(k)}
        if not changed:
            return dict(status='unchanged', input_file=str(source), drivers_started=0)
        if changed-FIELDS:
            raise ValueError('Unsupported live scientific changes: '+', '.join(sorted(changed-FIELDS))+
                             '. Use resource controls or an explicit scientific migration.')
        model_only = changed == {'vrc_tst_surface_models'}
        # Additive orbital reviews do not change L1/L2 exploration. Existing
        # branches quiesce; searching/new workers inherit at their ready gate.
        if m['phase'] not in ('exploring', 'postprocessing') or (not model_only and (
                m['phase'] != 'postprocessing' or any(w['status'] not in ('postprocessing','assembled')
                                                    for w in m['workers'].values()))):
            raise ValueError('Resubmit after PES exploration is complete; no files or jobs changed.')
        options = {k: proposed[k] for k in changed if k in proposed}
        if changed-set(options):
            raise ValueError('Live settings must be specified explicitly, not removed.')
        if 'vdW_detection' in options:
            from kinbot.network_revisit import settings
            settings(options['vdW_detection'])
        scope = prior.get('vrc_tst_surface_model_scope') or m.get('initial_worker')
        if not scope and (root/'chemids').exists():
            scope = (root/'chemids').read_text().splitlines()[0]
        if 'vrc_tst_surface_models' in options:
            options['vrc_tst_surface_models'] = validate_models(prior, proposed, source.parent, scope=scope)
        ident = uuid.uuid4().hex
        requests = []
        for name in sorted(m['workers']):
            if model_only and m['workers'][name]['status'] == 'exploring':
                continue
            worker = root/name
            if not phase.active_worker(worker):
                raise ValueError(f'{name}: no live worker to accept this revision.')
            if control.requests(worker, kinds={'reaction','network'}):
                raise ValueError(f'{name}: another search revision is pending.')
            inp = worker/(name+'.json')
            saved = json.loads(inp.read_text())
            if 'vrc_tst_surface_models' in options:
                validate_models(saved, dict(saved, **options), worker,
                                scope=saved.get('vrc_tst_surface_model_scope'))
            requests.append(dict(schema=1, id=ident+'_'+name, kind='network', target=str(worker),
                action='retry', expected=dict(input_sha256=file_sha256(inp)), options=options,
                cancel_running=False, pending_only=False, reason='Explicit PES input resubmission',
                input_revision=dict(path=str(root/JOURNAL), id=ident),
                status='queued', created_at=control.now()))
        record = dict(schema=1, id=ident, token=m['token'], status='queued',
            input_file=str(source), source_sha256=file_sha256(source), before=prior, after=proposed,
            previous_source_sha256=m['input_sha256'], options=options,
            requests=requests, control_root=str(control_path.parent), drivers_started=0,
            model_only=model_only)
        resolved=root/'resolved_parameters.json'
        if resolved.exists():
            values=json.loads(resolved.read_text())
            record['resolved']=dict(before=file_sha256(resolved), after=dict(values,**options))
        if apply:
            if proposal is not None and proposed != current:
                # Detect edits by tools which do not take these locks too.
                if json.loads(source.read_text()) != expected:
                    raise ValueError('Case input changed during review preparation; no input was overwritten.')
                atomic_json(source, proposed)
                record['source_sha256'] = file_sha256(source)
            if old:
                atomic_json(root/'.kinbot_input_history'/(old['id']+'.json'), old)
            atomic_json(root/JOURNAL, record)
        return {k: record[k] for k in ('id','status','input_file','drivers_started')} | {
            'status': 'queued' if apply else 'reviewed', 'changed': sorted(changed),
            'workers': len(requests), 'receipt': str(root/JOURNAL)}


def consume(root, m):
    """PES owner, under the phase lock. Replay partial publication idempotently."""
    root = Path(root)
    path = root/JOURNAL
    if not path.exists():
        return
    r = json.loads(path.read_text())
    if r['token']==m['token'] and r['status']=='published':
        m['input_sha256']=r['source_sha256']
        return
    if r['token'] != m['token'] or r['status'] not in ('queued','publishing'):
        return
    if r['status']=='queued' and file_sha256(r['input_file']) != r['source_sha256']:
        r.update(status='blocked', error='Input was edited again before the owner accepted it.')
        atomic_json(path, r)
        return
    if r['status'] == 'queued' and r.get('model_only'):
        # A search can become ready between submit and this locked owner tick.
        # Include it before publishing; later ready workers inherit below.
        from kinbot import live_control as control
        existing = {Path(q['target']).name for q in r['requests']}
        for name, worker in sorted(m['workers'].items()):
            if name in existing or worker['status'] == 'exploring':
                continue
            directory = root/name
            if control.requests(directory, kinds={'reaction','network'}):
                r.update(status='blocked',error=f'{name}: another search revision arrived before publication.')
                atomic_json(path,r)
                return
            inp = directory/(name+'.json')
            saved = json.loads(inp.read_text())
            validate_models(saved, dict(saved, **r['options']), directory,
                            scope=saved.get('vrc_tst_surface_model_scope'))
            r['requests'].append(dict(schema=1, id=r['id']+'_'+name, kind='network',
                target=str(directory), action='retry', expected=dict(input_sha256=file_sha256(inp)),
                options=r['options'], cancel_running=False, pending_only=False,
                reason='Explicit PES orbital-model resubmission',
                input_revision=dict(path=str(path), id=r['id']), status='queued',
                created_at=control.now()))
    if r['status'] == 'queued':
        from kinbot import live_control as control
        for request in r['requests']:
            directory=Path(request['target'])
            inp=directory/(directory.name+'.json')
            if (control.requests(directory,kinds={'reaction','network'})
                    or file_sha256(inp)!=request['expected']['input_sha256']):
                r.update(status='blocked',error=f'{directory.name}: worker revision changed before publication.')
                atomic_json(path,r)
                return
    r['status']='publishing'; atomic_json(path,r)
    for request in r['requests']:
        dest=Path(r['control_root'])/'.kinbot_controls/requests'/(request['id']+'.json')
        if not dest.exists():
            atomic_json(dest, request)
    snapshot(root,r['input_file'],r['after'])
    # Preserve the prepared-run verifier's exact-byte contract. Separate PES
    # owners serialize this update; unrelated hashes are never recomputed.
    manifests=[p/'prepared.json' for p in (root,*root.parents) if (p/'prepared.json').exists()]
    if manifests:
        manifest=manifests[0]
        with (manifest.parent/'.kinbot_prepared.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            prepared=json.loads(manifest.read_text())
            source=Path(r['input_file'])
            name=str(source.relative_to(manifest.parent)) if source.is_relative_to(manifest.parent) else None
            expected=prepared['sha256'].get(name)
            if expected not in (None,r['previous_source_sha256'],r['source_sha256']):
                raise ValueError('Prepared input hash changed outside this revision.')
            if expected is not None:
                prepared['sha256'][name]=r['source_sha256']
                if 'resolved' in r:
                    resolved=root/'resolved_parameters.json'
                    key=str(resolved.relative_to(manifest.parent))
                    current=file_sha256(resolved)
                    if current!=r['resolved']['before'] and json.loads(resolved.read_text())!=r['resolved']['after']:
                        raise ValueError('Resolved settings changed outside this revision.')
                    if prepared['sha256'].get(key) not in (None,r['resolved']['before'],current):
                        raise ValueError('Prepared resolved-settings hash changed outside this revision.')
                    atomic_json(resolved,r['resolved']['after'])
                    if key in prepared['sha256']:
                        prepared['sha256'][key]=file_sha256(resolved)
                atomic_json(manifest,prepared)
    m['input_sha256']=r['source_sha256']
    r['status']='published'; atomic_json(path,r)


def inherit_models(root, generator):
    """At a worker's locked ready boundary, inherit accepted additive reviews.

    The owner lock excludes partial publication. Never copy a proposal or
    change a running reaction/branch; only its future VRC model is extended.
    """
    root = Path(root)
    if generator is None or not (root/JOURNAL).exists():
        return
    r = json.loads((root/JOURNAL).read_text())
    from kinbot import pes_phase, live_control as control
    if (r['token'] != json.loads((root/pes_phase.NAME).read_text())['token']
            or r['status'] != 'published' or 'vrc_tst_surface_models' not in r['options']):
        return
    models = r['options']['vrc_tst_surface_models']  # normalized absolute paths
    prior = generator.par
    if prior.get('vrc_tst_surface_models', {}) == models:
        return
    inp = Path(generator.inp).resolve()
    if any(Path(q['target']) == inp.parent for _, q in control.requests(inp.parent, kinds={'network'})):
        return  # Its generation-fenced request owns this edit and acknowledgement.
    validate_models(prior, dict(prior, vrc_tst_surface_models=models), inp.parent,
                    scope=prior.get('vrc_tst_surface_model_scope'))
    saved = json.loads(inp.read_text())
    saved['vrc_tst_surface_models'] = deepcopy(models)
    atomic_json(inp, saved)
    path = inp.parent/'.kinbot_search_options.json'
    overrides = json.loads(path.read_text()) if path.exists() else {}
    overrides['vrc_tst_surface_models'] = deepcopy(models)
    atomic_json(path, overrides)
    generator.par['vrc_tst_surface_models'] = deepcopy(models)
    from kinbot.search_state import rebase
    rebase(inp, r['id'])
    if getattr(generator, 'search_checkpoint', None):
        generator.search_checkpoint.adopt_revision()
    from kinbot.vrc_review_status import acknowledge_installed
    acknowledge_installed(generator.par, generator.species, inp.parent)


def receipt(root):
    root=Path(root)
    r=json.loads((root/JOURNAL).read_text())
    statuses=[]
    for child in r['requests']:
        p=Path(r['control_root'])/'.kinbot_controls/requests'/(child['id']+'.json')
        statuses.append(json.loads(p.read_text()) if p.exists() else child)
    status=r['status']
    if status=='published':
        if any(v['status']=='blocked' for v in statuses):
            status='blocked'
        elif all(v['status']=='applied' for v in statuses):
            status='applied'
    return dict(id=r['id'], status=status, input_file=r['input_file'], drivers_started=0,
                workers=[{k:v.get(k) for k in ('target','status','error')} for v in statuses])
