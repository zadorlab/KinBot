"""Finite verification using existing recon and ANL owners. No model publication.

Stage is offline/read-only with respect to all source runs. Run is an explicit
separate operation; repeated runs adopt the same manifests and dependency graphs.
"""
import argparse
from contextlib import ExitStack
import fcntl
from copy import deepcopy
import json
from pathlib import Path
import time
import tempfile
import numpy as np
from ase.io import read as read_atoms
from kinbot import molpro_reconnaissance as recon, run_policy
from kinbot.anl import dispatch
from kinbot.execution import atomic_json, file_sha256
from kinbot.vrc_graph import select_graph

HERE=Path(__file__).resolve().parent
MAX_JOBS=23

def read(p):return json.loads(Path(p).read_text())

def scientific_plan(plan):
    value=deepcopy(plan)
    value['protocol'].pop('max_jobs',None)
    return value

def partition(electrons, source, target_e, indices):
    """Pack explicitly specified source-natural indices, preserving electron count."""
    se,so=source['active_electrons'],source['active_orbitals'];c=(electrons-se)//2
    indices=sorted(indices);active=list(range(c+1,c+so+1))
    opened=[i for i in indices if i<=c];external=[i for i in indices if i>c+so]
    keep=[i for i in indices if i in active]
    close=(se+2*len(opened)-target_e)//2
    if se+2*len(opened)-target_e!=2*close or close<0:raise ValueError('Invalid electron partition.')
    other=[i for i in active if i not in keep]
    if close>len(other) or (opened and close):raise ValueError('Incompatible closure/promotion.')
    if external!=list(range(c+so+1,c+so+1+len(external))):raise ValueError('External trial window must be explicit and contiguous.')
    desired=other[:close]+keep+other[close:];order=active[:];swaps=[]
    for n,want in enumerate(desired):
        k=order.index(want)
        if k!=n:swaps.append([n+1,k+1]);order[n],order[k]=order[k],order[n]
    closed_order=list(range(1,c+1));closed_swaps=[]
    for n,want in zip(range(c-len(opened),c),opened):
        k=closed_order.index(want)
        if n!=k:closed_swaps.append([n+1,k+1]);closed_order[n],closed_order[k]=closed_order[k],closed_order[n]
    p=dict(source_active_electrons=se,source_active_orbitals=so,close_lowest=close,
           discard_highest=len(other)-close,orbital_swaps=swaps)
    if opened:p.update(open_highest_closed=len(opened),closed_orbital_swaps=closed_swaps)
    if external:p['open_lowest_external']=len(external)
    if external and p['discard_highest']:raise ValueError('Cannot replace active space with unspecified external orbitals.')
    return p


def anchors(source):
    xyz=np.array(source['geometry']);coord=source['coordinate'];i,j=[k-1 for k in coord['bond']]
    near=float(np.linalg.norm(xyz[j]-xyz[i]));mid=max(3.5,near+.75)
    if mid>=6:raise ValueError('Unexpected near geometry; review the finite grid.')
    base=dict(geometry_angstrom=xyz.tolist())
    def displaced(r,rotation=None):
        a=dict(coord,distance_angstrom=r,name='point')
        if rotation:a['_rotation']=rotation
        return recon._anchor_geometry(base,a).tolist()
    result=[]
    for name,r in [('near',near),('mid',mid),('r6',6.),('r20',20.),('r30',30.),('back6',6.),('back_near',near)]:
        result.append(dict(name=name,orbital_chain='radial_return',geometry_angstrom=displaced(r)))
    moving=[k-1 for k in coord['moving_atoms']];fixed=sorted(set(range(len(xyz)))-set(moving))
    direction=(xyz[j]-xyz[i])/near;axis=np.cross(direction,np.eye(3)[np.argmin(abs(direction))]);axis/=np.linalg.norm(axis)
    axes=[axis,np.cross(direction,axis)];fragments=[f for f in (fixed,moving) if len(f)>1]
    # Two independent fragment orientations, with their own separated references.
    # If one fragment is atomic rotate its partner about two transverse axes.
    turns=[(fragments[n%len(fragments)],axes[n]) for n in range(2)]
    for n,(atoms,axis) in enumerate(turns):
        pivot=i if i in atoms else j
        rotation=dict(atoms=atoms,pivot=pivot,axis=axis.tolist(),angle_degrees=60.)
        for label,r in [('angle',mid),('angle_far',30.)]:
            result.append(dict(name=f'{label}_{n+1}',orbital_chain=f'angular_{n+1}',geometry_angstrom=displaced(r,rotation)))
    if len({tuple(np.round(np.array(a['geometry_angstrom']).ravel(),8)) for a in result})!=9:
        raise ValueError('Angular test is a duplicate geometry; select a different rotation.')
    return result


def configuration(case,variant,basis,resolve, *, geometries=None, seed=None):
    s=case['sources'][basis];e=variant['electrons'];inds=variant['indices'][basis];o=len(inds)
    points=anchors(s) if geometries is None else geometries
    model=deepcopy(case['state_model']);model.pop('basis_models',None);model.pop('orbital_selection',None)
    part=partition(s['electrons'],s['source_space'],e,inds)
    entry={'path':str(resolve(s['wfu'])),'sha256':file_sha256(resolve(s['wfu'])),'basis':basis} if seed is None else seed
    if seed is not None:
        part=dict(source_active_electrons=e,source_active_orbitals=o,close_lowest=0,discard_highest=0)
    checks=[dict(method=f'caspt2-f12({e},{o})',shift=.1,ipea=.25,reference_policy='state_average')]
    if case['h_exit'] and variant['name']=='candidate':
        checks.append(dict(method=f'mrci+q({e},{o})',reference_policy='state_average'))
    cfg=dict(schema=1,title=f'v{case["id"]:02d}_{variant["name"]}_{basis.replace("-", "_")}',
        symbols=s['symbols'],geometry_angstrom=s['geometry'],charge=0,multiplicity=2,
        bases=[dict(name='dz' if basis=='avdz-f12' else 'tz',basis=basis)],anchors=points,
        probe_active_space=True,probe_method='casscf',orbital_initialization='propagated',
        probe_definition=dict(active_electrons=e,active_orbitals=o,state_model=model,
            chemical_note=variant['purpose'],orbital_partition=part),
        correlation_checks=checks,initial_records={},ppn=8,memory_mw=1000,memory_mb=104000,
        partition='month-long-cpu',walltime='partition_max',scf_convergence=1.e-8,
        scf_options=dict(solver='so_sci',max_iterations=2000),orbital_print_virtuals=12,orbital_print_threshold=.1)
    seen=set()
    for a in points:
        chain=a.get('orbital_chain',a['name'])
        if chain not in seen:cfg['initial_records'][cfg['bases'][0]['name']+'__'+a['name']]=entry
        seen.add(chain)
    return cfg


def reuse_fixed(case,variant,basis,cfg,resolve):
    """Reuse only the explicitly reconciled historical trials, never failed branches."""
    matches={}
    for anchor in cfg['anchors']:
        if anchor['name'] not in ('near','mid','r6','r20','r30'):continue
        for old in case.get('reuse_fixed',[]):
            j=old['job']
            if old['variant']==variant['name'] and old['basis']==basis and np.allclose(j['geometry_angstrom'],anchor['geometry_angstrom'],atol=1.e-8,rtol=0):
                copied=deepcopy(j)
                copied['directory']=str(resolve(old['files']['.out']).parent)
                copied['wfu']=str(resolve(old['files']['.wfu']))
                copied['anchor']=anchor['name']
                matches[anchor['name']]=copied
                break
    if not matches:return cfg,{}
    if 'near' not in matches:raise ValueError('Baseline reuse needs an established near seed for all new chains.')
    target=cfg['probe_definition'];target['orbital_partition']=dict(source_active_electrons=target['active_electrons'],source_active_orbitals=target['active_orbitals'],close_lowest=0,discard_highest=0)
    prefix=cfg['bases'][0]['name'];cfg['initial_records']={};points=[];last=matches['near'];segment=0;seen=set()
    for anchor in cfg['anchors']:
        if anchor['name'] in matches:
            last=matches[anchor['name']];segment+=1;continue
        a=deepcopy(anchor)
        if a['orbital_chain']=='radial_return':a['orbital_chain']+=str(segment);seed=last
        else:seed=matches['near']
        chain=a['orbital_chain']
        if chain not in seen:
            cfg['initial_records'][prefix+'__'+a['name']]=dict(path=seed['wfu'],sha256=seed['wfu_sha256'],basis=basis)
        seen.add(chain);points.append(a)
    cfg['anchors']=points
    return cfg,matches


def source_resolver(plan, export_root=None):
    files=plan['source_files']
    def resolve(member):
        f=files[member]
        return (Path(export_root)/member) if export_root else Path(f['source'])
    return resolve


def validate_sources(plan,resolve):
    for member,f in plan['source_files'].items():
        p=resolve(member)
        if not p.is_file() or file_sha256(p)!=f['sha256']:raise ValueError('Missing/changed source: '+str(p))


def prepare_configuration(root,cfg,policy):
    path=root/'input.json'
    if path.exists() and read(path)!=cfg:raise ValueError('Existing scientific configuration differs: '+str(path))
    if not path.exists():atomic_json(path,cfg)
    manifest=recon.prepare(path,root/'native')
    m=read(manifest)
    if policy:
        if m.get('shared_policy_file',str(policy))!=str(policy):raise ValueError('Changed scheduling owner.')
        m.update(shared_policy_file=str(policy),startup_retries=2)
        atomic_json(manifest,m)
    return str(manifest)


def stage(output, *, export_root=None, pool=None):
    plan=read(HERE/'plan.json')
    validate_sources(plan,source_resolver(plan,export_root))
    with recon._manifest_lock(output):
        return _stage(output,export_root=export_root,pool=pool)


def _stage(output, *, export_root=None, pool=None):
    output=Path(output).resolve();plan=read(HERE/'plan.json');resolve=source_resolver(plan,export_root)
    validate_sources(plan,resolve) # all inputs before creating any calculation directory
    policy=None
    if pool is not None:
        pool=Path(pool).resolve(strict=True);policy=run_policy.find(pool)
        p=run_policy.read(pool)
        if policy!=pool/run_policy.NAME or not p or not 0<p['max_jobs']<=MAX_JOBS or not output.is_relative_to(pool) or output==pool:
            raise ValueError(f'Require a separate existing pool capped at {MAX_JOBS} and a child output directory.')
    elif export_root is None:raise ValueError('Native staging requires the explicit detached pool.')
    # Validate every deck before committing the first calculation directory.
    configs=[]
    for case in plan['cases']:
        for variant in case['variants']:
            for basis in plan['protocol']['bases']:
                cfg=configuration(case,variant,basis,resolve)
                cfg,reused=reuse_fixed(case,variant,basis,cfg,resolve)
                configs.append((case,variant,basis,cfg,reused))
    for case in plan['cases']:
        if 'far_sources' not in case:continue
        for basis in plan['protocol']['bases']:
            for variant in case['variants']:
                alternatives=[18,19] if variant['name']=='candidate' else [None]
                for donor in alternatives:
                    v=deepcopy(variant)
                    if donor is not None:v['indices'][basis]=[donor,20,21,22]
                    v['name']+='_'+'far'+(str(donor) if donor is not None else '')
                    c=deepcopy(case);c['sources'][basis]['wfu']=case['far_sources'][basis]['wfu']
                    original=anchors(case['sources'][basis]);byname={a['name']:a for a in original}
                    points=[dict(name='in_'+k,orbital_chain='independent_inward',geometry_angstrom=byname[k]['geometry_angstrom']) for k in ('r30','r6','near')]
                    cfg=configuration(c,v,basis,resolve,geometries=points)
                    configs.append((case,v,basis,cfg,{}))
    with tempfile.TemporaryDirectory(prefix='kinbot-verification-preflight-') as tmp:
        for n,(_,_,_,cfg,_) in enumerate(configs):
            p=Path(tmp)/(str(n)+'.json');atomic_json(p,cfg);loaded=recon._load(p)
            basis=loaded['bases'][0]
            for a in loaded['anchors']:
                key=basis['name']+'__'+a['name'];previous=cfg['initial_records'].get(key,'predecessor')
                recon.validate_input_blocks(recon._deck(loaded,basis['basis'],a,key,previous=previous),'verification')
    receipt=output/'panel.json'
    key=file_sha256(HERE/'plan.json')
    if receipt.exists():
        old=read(receipt)
        if old['policy']!=(str(policy) if policy else None):raise ValueError('Existing panel/owner differs.')
        if old['plan_sha256']!=key:
            saved=output/'plan.json'
            if scientific_plan(read(saved))!=scientific_plan(plan):
                raise ValueError('Existing scientific panel differs.')
            atomic_json(output/'plan_history'/(file_sha256(saved)+'.json'),read(saved))
    systems=[]
    for case,variant,basis,cfg,reused in configs:
        folder=output/'fixed'/cfg['title'];manifest=prepare_configuration(folder,cfg,policy)
        systems.append(dict(case_id=case['id'],variant=variant['name'],basis=basis,configuration=str(folder/'input.json'),manifest=manifest,reused=reused,diagnostic_branch=variant['name'] not in {v['name'] for v in case['variants']}))
    record=dict(schema=1,plan_sha256=key,source_export_root=str(Path(export_root).resolve()) if export_root else None,
        policy=str(policy) if policy else None,preview=policy is None,systems=systems,
        status='staged_not_scientifically_approved',production_model_accepted=False)
    atomic_json(receipt,record)
    atomic_json(output/'plan.json',plan)
    return dict(directory=str(output),configurations=len(systems),fixed_jobs=sum(len(read(s['manifest'])['jobs']) for s in systems),
                covered_channels=29,excluded_channels=5,production_model_accepted=False)


def verified_job(system,name):
    m=read(system['manifest']);j=system.get('reused',{}).get(name)
    if j is None:j=next(j for j in m['jobs'] if j['anchor']==name)
    if j.get('status')!='complete' or j.get('validation_errors'):raise ValueError('Native prerequisite incomplete: '+j['name'])
    for path,key in [(Path(j['directory'])/(j['name']+'.out'),'output_sha256'),(Path(j['directory'])/(j['name']+'.inp'),'input_sha256'),(Path(j['wfu']),'wfu_sha256')]:
        if file_sha256(path)!=j[key]:raise ValueError('Changed native prerequisite: '+str(path))
    return j


def relaxation_spec(case,job,distance_label):
    """Diagnostic trial, not a fabricated reviewed selection or production model."""
    variant=case['variants'][0];e=variant['electrons'];o=len(variant['indices']['avdz-f12']);c=(case['sources']['avdz-f12']['electrons']-e)//2
    model=deepcopy(case['state_model']);model['orbital_selection']=dict(mode='molpro_record',path=job['wfu'],sha256=job['wfu_sha256'],
        record='2140.2',orbital_type='natural',align_active_orbitals=True,active_indices=list(range(c+1,c+o+1)))
    pair=[i-1 for i in case['sources']['avdz-f12']['coordinate']['bond']]
    return dict(schema=1,name=f'verify_{case["id"]}_{distance_label}',molecule=dict(symbols=job['symbols'],positions=job['geometry_angstrom'],charge=0,multiplicity=2),
        limits=dict(max_nodes=1),intent=dict(production_model_accepted=False,verification_case=case['id'],anchor=distance_label,
            source_output_sha256=job['output_sha256'],vrc_fragment_retry=True,cap_walltime_to_explicit_partition=True),
        tasks=[dict(id='relaxed_scan',kind='ase_optimize',geometry_from='initial',geometry_output='final.xyz',
            profile=dict(calculator='molpro',method=f'caspt2-f12({e},{o})',basis='avdz-f12',command='molpro',optimizer='sella',
                calculator_kwargs=dict(state_model=model,caspt2_shift=.1,caspt2_ipea=.25,gradient_mode='numerical',scf_convergence=1.e-8,correlation_convergence=1.e-8)),
            optimizer=dict(fmax=.0005,steps=50,step_retries=0,order=0,sella_kwargs=dict(internal=True),constraints=dict(fixed_bonds=[pair])),
            resources=dict(cores=8,memory_mb=104000,min_stack_mw=1000,partition='month-long-cpu',walltime='31-00:00:00'))])


def stage_relaxations(output):
    output=Path(output);record=read(output/'panel.json');plan=read(output/'plan.json');result=[]
    resolve=source_resolver(plan,record['source_export_root']);validate_sources(plan,resolve)
    holds=[]
    for case in plan['cases']:
        system=next(s for s in record['systems'] if s['case_id']==case['id'] and s['variant']=='candidate' and s['basis']=='avdz-f12')
        for label in ('near','mid'):
            if label in case.get('reuse_relaxations',{}):
                result.append(dict(case_id=case['id'],anchor=label,reused=case['reuse_relaxations'][label]))
                continue
            try:job=verified_job(system,label)
            except ValueError as error:
                holds.append(dict(case_id=case['id'],anchor=label,reason=str(error)));continue
            spec=relaxation_spec(case,job,label);dispatch.validate_spec(spec)
            graph=select_graph(output/'relaxed',f'v{case["id"]:02d}_{label}',spec)
            result.append(dict(case_id=case['id'],anchor=label,graph=str(graph)))
    atomic_json(output/'relaxations.json',dict(graphs=result,holds=holds,production_model_accepted=False))
    return result


def stage_relaxed_points(output):
    output=Path(output);record=read(output/'panel.json');plan=read(output/'plan.json');cases={c['id']:c for c in plan['cases']}
    resolve=source_resolver(plan,record['source_export_root']);systems=[];holds=[]
    for r in read(output/'relaxations.json')['graphs']:
        if r.get('reused'):
            final=resolve(r['reused']['geometry'])
        else:
            try:final=verified_final(Path(r['graph']))
            except ValueError as error:
                holds.append(dict(case_id=r['case_id'],anchor=r['anchor'],reason=str(error)));continue
        case=cases[r['case_id']];a=read_atoms(final)
        if a.get_chemical_symbols()!=case['sources']['avdz-f12']['symbols']:raise ValueError('Relaxation changed atom identity.')
        for variant in case['variants']:
            for basis in plan['protocol']['bases']:
                system=next(s for s in record['systems'] if s['case_id']==case['id'] and s['variant']==variant['name'] and s['basis']==basis)
                try:job=verified_job(system,r['anchor'])
                except ValueError as error:
                    holds.append(dict(case_id=case['id'],variant=variant['name'],basis=basis,anchor=r['anchor'],reason=str(error)));continue
                seed=dict(path=job['wfu'],sha256=job['wfu_sha256'],basis=basis)
                points=[dict(name='relaxed',geometry_angstrom=a.positions.tolist())]
                cfg=configuration(case,variant,basis,resolve,geometries=points,seed=seed)
                cfg['title']+='_'+r['anchor']+'_relaxed'
                folder=output/'relaxed_points'/cfg['title'];manifest=prepare_configuration(folder,cfg,record['policy'])
                systems.append(dict(case_id=case['id'],variant=variant['name'],basis=basis,anchor=r['anchor'],manifest=manifest,configuration=str(folder/'input.json'),geometry_sha256=file_sha256(final)))
    atomic_json(output/'relaxed_points.json',dict(systems=systems,holds=holds,production_model_accepted=False))
    return systems


def verified_final(g):
    state=read(g/'state.json')
    if state['tasks']['relaxed_scan']['status']!='complete':raise ValueError('Relaxation incomplete: '+str(g))
    d=g/'tasks/relaxed_scan';ex=read(d/'execution.json');final=d/'final.xyz'
    if ex.get('status')!='executed' or ex.get('artifacts',{}).get('final.xyz')!=file_sha256(final):raise ValueError('Unverified final geometry: '+str(final))
    return final

def run(output):
    output=Path(output).resolve(strict=True);record=read(output/'panel.json')
    if record['preview'] or run_policy.find(output)!=Path(record['policy']):raise ValueError('Preview or changed pool: execution refused.')
    with recon._manifest_lock(output):
        # Recheck immutable inputs and all original sources before adoption/submission.
        _stage(output,pool=Path(record['policy']).parent)
        failures=[]
        try:recon.run_many([s['manifest'] for s in record['systems']],max_jobs=MAX_JOBS,keep_going=True,allocation_retries=2)
        except recon.ProbeFailure as error:failures.append(str(error))
        graphs=[g for g in stage_relaxations(output) if 'graph' in g]
        with ExitStack() as locks:
            for r in graphs:
                handle=locks.enter_context((Path(r['graph'])/'drive.lock').open('a'))
                fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
            for r in graphs:dispatch.preflight_until_scheduler_ready(r['graph'])
            while True:
                try:states=[dispatch.advance(r['graph'],submit=True) for r in graphs]
                except dispatch.SchedulerUnavailableError:
                    time.sleep(20);continue
                statuses=[e['status'] for s in states for e in s['tasks'].values()]
                if all(s in ('complete','failed','unavailable','blocked','skipped') for s in statuses):break
                time.sleep(20)
        systems=stage_relaxed_points(output)
        if systems:
            try:recon.run_many([s['manifest'] for s in systems],max_jobs=MAX_JOBS,keep_going=True,allocation_retries=2)
            except recon.ProbeFailure as error:failures.append(str(error))
        holds=read(output/'relaxations.json')['holds']+read(output/'relaxed_points.json')['holds']
        atomic_json(output/'completion.json',dict(status='incomplete_review_required' if failures or holds else 'native_complete_review_required',failures=failures,holds=holds,production_model_accepted=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['stage','run']);p.add_argument('output',type=Path)
    p.add_argument('--export-root',type=Path);p.add_argument('--pool',type=Path)
    a=p.parse_args()
    if a.action=='stage':print(json.dumps(stage(a.output,export_root=a.export_root,pool=a.pool),indent=2))
    else:run(a.output)
