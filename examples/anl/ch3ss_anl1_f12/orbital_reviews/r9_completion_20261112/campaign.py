"""One diagnostic campaign, using KinBot's native recon owner and shared pool.

The adaptive steps enumerate test spaces and choose further diagnostic trials.
They never approve a production CAS model. Original evidence is read-only.
"""
import argparse
from copy import deepcopy
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
from ase.data import atomic_numbers

from kinbot import molpro_reconnaissance as recon, run_policy
from kinbot.execution import atomic_json, file_sha256
from kinbot.molpro import automatic_vrc_state_model
from kinbot.radical_pair_evidence import read_orbitals, projection

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('prior_verification', HERE.parent/'verification_20261027/panel.py')
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)
KCAL = 627.5094740631


def read(path):
    return json.loads(Path(path).read_text())


def immutable(path, value):
    value=json.loads(json.dumps(value,allow_nan=False))
    if path.exists():
        if read(path) != value:
            raise ValueError('Existing scientific input differs: '+str(path))
    else:
        atomic_json(path, value)


def state(roots=1):
    return dict(state_average=[dict(multiplicity=2, roots=roots, weights=[1./roots]*roots)],
                selected_roots={'2': 1}, canonical_orbitals=True,
                casscf_options=dict(solver='so_sci', max_iterations=999,
                                    energy_convergence=1.e-8, gradient_convergence=1.e-5),
                correlation_options=dict(reference_policy='state_average', max_iterations=2000,
                                         max_reference_iterations=2000, max_reference_vectors=60,
                                         energy_convergence=1.e-8, coefficient_convergence=1.e-5))


def base_config(title, symbols, geometry, bases, mult=2):
    return dict(schema=1, title=title, symbols=symbols, geometry_angstrom=geometry,
                charge=0, multiplicity=mult, bases=[dict(name='dz' if b=='avdz-f12' else 'tz', basis=b) for b in bases],
                ppn=8, memory_mw=1000, memory_mb=104000, partition='month-long-cpu',
                walltime='partition_max', scf_convergence=1.e-8,
                scf_options=dict(solver='so_sci', max_iterations=2000),
                orbital_print_virtuals=30, orbital_print_threshold=.05)


def h_anchors(case):
    source=dict(geometry=case['geometry'], coordinate=case['coordinate'])
    source['geometry']=recon._anchor_geometry(
        {'geometry_angstrom': source['geometry']},
        dict(source['coordinate'], name='near', distance_angstrom=2.)).tolist()
    points=prior.anchors(source)
    byname={a['name']: a for a in points}
    # Begin with separated ground fragments. Every angular branch also starts far.
    result=[]
    for name in ('r30','r20','r6','mid','near','back6','back_near'):
        result.append(byname[name])
    for number in (1,2):
        result.extend([byname[f'angle_far_{number}'],byname[f'angle_{number}']])
    result.append(dict(name='parent',orbital_chain='bound_parent',geometry_angstrom=case['parent_geometry']))
    return result


def reference_config(case, bases, *, fragment_mult=None, quartet=False):
    space=case['reference_space']
    symbols=case['symbols']; geometry=case['geometry']; mult=4 if quartet else 2
    title='h'+str(case['id'])+'_reference'
    anchors=h_anchors(case)
    if quartet:
        title+='q';anchors=[a for a in anchors if a['name'] in ('r30','mid')]
    if fragment_mult is not None:
        space=case['fragment_space'];symbols=[symbols[i] for i in case['fragment_map']]
        geometry=case['fragment_geometry'];mult=fragment_mult
        title='h'+str(case['id'])+'_fragment'+str(mult)
        anchors=[dict(name='fragment',geometry_angstrom=geometry)]
    model=automatic_vrc_state_model(symbols,0,[mult],active_space=space)
    model['casscf_options']=state()['casscf_options']
    model['correlation_options']=state()['correlation_options']
    e,o=space['active_electrons'],space['active_orbitals']
    cfg=base_config(title,symbols,geometry,bases,mult)
    cfg.update(probe_active_space=True,probe_method='casscf',orbital_initialization='propagated',
               anchors=anchors,probe_definition=dict(active_electrons=e,active_orbitals=o,
               state_model=model,chemical_note='Heavy-atom valence reference plus the departing H 1s for the combined system; native character and spin checks required.'),
               correlation_checks=[dict(method=f'caspt2-f12({e},{o})',shift=.1,ipea=.25,reference_policy='state_average')])
    return cfg


def prepare(root, cfg, policy=None):
    directory=root/cfg['title'];directory.mkdir(parents=True,exist_ok=True)
    immutable(directory/'input.json',cfg)
    manifest=recon.prepare(directory/'input.json',directory/'native')
    data=read(manifest)
    if policy:
        if data.get('shared_policy_file',str(policy)) != str(policy):
            raise ValueError('Existing manifest uses a different shared pool.')
        data.update(shared_policy_file=str(policy),startup_retries=2)
        atomic_json(manifest,data)
    return str(manifest)


def resolve_seed(seed, source_map):
    result=deepcopy(seed)
    for key in ('path','output'):
        original=result[key];result[key]=source_map.get(original,original)
        expected=result['sha256' if key=='path' else 'output_sha256']
        if file_sha256(result[key])!=expected:
            raise ValueError('Changed original native evidence: '+original)
    return result


def existing_configs(plan, source_map):
    configs=[]
    for check in plan['existing_checks']:
        case=check['case'];cid=check['case_id']
        for basis in plan['bases']:
            seed=resolve_seed(check['seeds'][basis],source_map)
            seed={key: seed[key] for key in ('path','sha256','basis')}
            e,o=(7,6) if cid==7 else (5,4)
            source=deepcopy(case['sources'][basis])
            source['source_space']=dict(active_electrons=e,active_orbitals=o)
            source['wfu']=seed['path']
            variant=dict(name='candidate',electrons=e,indices={basis:list(range((source['electrons']-e)//2+1,(source['electrons']-e)//2+o+1))},purpose='Retain the established chemical subspace; verify corrected state tracking or exact endpoint binding.')
            current=deepcopy(case);current['sources'][basis]=source
            if cid==7:
                points={a['name']:a for a in prior.anchors(source)}
                for n in (1,2):
                    anchors=[]
                    for name,source_name in [('far',f'angle_far_{n}'),('in',f'angle_{n}'),('back_far',f'angle_far_{n}'),('back_in',f'angle_{n}')]:
                        a=deepcopy(points[source_name]);a.update(name=name,orbital_chain='angular_return');anchors.append(a)
                    cfg=prior.configuration(current,variant,basis,lambda path:Path(path),geometries=anchors,seed=seed)
                    cfg['title']=f'v07_farstart_a{n}_'+basis.replace('-','_');configs.append(cfg)
            else:
                source['geometry']=check['target_geometry'];source['coordinate']=check['target_coordinate']
                points={a['name']:a for a in prior.anchors(source)}
                anchors=[points[k] for k in ('near','mid','r30','back_near')]
                cfg=prior.configuration(current,variant,basis,lambda path:Path(path),geometries=anchors,seed=seed)
                cfg['title']='v13_r9_binding_'+basis.replace('-','_');configs.append(cfg)
    return configs


def stage(output, source_map=None, pool=None):
    with recon._manifest_lock(output):
        return _stage(output,source_map,pool)


def _stage(output, source_map=None, pool=None):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    plan=read(HERE/'plan.json');source_map=source_map or {};policy=None
    if pool is not None:
        pool=Path(pool).resolve(strict=True);policy=run_policy.find(pool)
        settings=run_policy.read(pool)
        if policy!=pool/run_policy.NAME or not settings or not 0<settings['max_jobs']<=plan['max_jobs'] or not output.is_relative_to(pool):
            raise ValueError('Require the explicit shared pool capped at 20 and a child campaign directory.')
    configs=existing_configs(plan,source_map)
    for case in plan['h_loss_cases']:
        configs.extend([reference_config(case,plan['bases']),reference_config(case,plan['bases'],quartet=True),
                        reference_config(case,plan['bases'],fragment_mult=1),reference_config(case,plan['bases'],fragment_mult=3)])
        far=next(a['geometry_angstrom'] for a in h_anchors(case) if a['name']=='r30')
        h=case['coordinate']['moving_atoms'][0]-1
        cfg=base_config('h'+str(case['id'])+'_hydrogen',['H'],[far[h]],plan['bases'])
        cfg.update(orbital_character_only=True,anchors=[dict(name='hydrogen',geometry_angstrom=[far[h]])]);configs.append(cfg)
    # Validate every independent configuration before preparing any job manifest.
    for cfg in configs:
        path=output/'preflight'/ (cfg['title']+'.json');path.parent.mkdir(exist_ok=True)
        immutable(path,cfg);loaded=recon._load(path)
        for basis in loaded['bases']:
            for anchor in loaded['anchors']:
                name=basis['name']+'__'+anchor['name']
                deck=recon._deck(loaded,basis['basis'],anchor,name,previous=cfg.get('initial_records',{}).get(name))
                recon.validate_input_blocks(deck,'R9 completion')
    manifests=[prepare(output/'initial',cfg,policy) for cfg in configs]
    record=dict(schema=1,plan_sha256=file_sha256(HERE/'plan.json'),preview=policy is None,
                policy=str(policy) if policy else None,manifests=manifests,production_model_accepted=False)
    immutable(output/'campaign.json',record);immutable(output/'plan.json',plan)
    return record


def completed_jobs(manifests):
    result={}
    for path in manifests:
        for job in read(path)['jobs']:
            if job.get('status')!='complete' or job.get('validation_errors'):
                continue
            directory=Path(job['directory'])
            for suffix,key in (('.inp','input_sha256'),('.out','output_sha256'),('.wfu','wfu_sha256')):
                file=Path(job['wfu']) if suffix=='.wfu' else directory/(job['name']+suffix)
                if file_sha256(file)!=job[key]:raise ValueError('Changed native result: '+str(file))
            result[(read(Path(path).parents[1]/'input.json')['title'],job['basis'],job['anchor'])]=job
    return result


def energy(job):
    matches=[x['energy_hartree'] for x in job.get('correlated_energies_hartree',[])
             if x['method'].lower().startswith('caspt2-f12') and x['multiplicity']==2 and x['root']==1]
    if len(matches)!=1:raise ValueError('Need one converged CASPT2-F12 doublet-ground energy.')
    return matches[0]


def hydrogen_partition(case, reference, atomic):
    d=read_orbitals(Path(reference['directory'])/(reference['name']+'.natural.orbitals.xml'))
    h=read_orbitals(Path(atomic['directory'])/(atomic['name']+'.scf.orbitals.xml'))
    atom=case['coordinate']['moving_atoms'][0]-1
    if h.symbols!=['H'] or h.atom_basis[0]!=d.atom_basis[atom] or not np.allclose(h.positions[0],d.positions[atom],rtol=0,atol=2.e-8):
        raise ValueError('Hydrogen and combined-system AO bases/coordinates differ.')
    singly=np.flatnonzero(abs(h.occupations-1.)<1.e-5)
    if len(singly)!=1 or abs(sum(h.occupations)-1.)>1.e-5:
        raise ValueError('Independent hydrogen reference must have one electron in one SOMO.')
    orbital=int(singly[0]);embedded=np.zeros(d.coefficients.shape[1])
    embedded[d.atom_columns[atom]]=h.coefficients[orbital,h.atom_columns[0]]
    amplitudes=embedded@np.linalg.pinv(d.coefficients,rcond=1.e-12)
    if abs(np.dot(amplitudes,amplitudes)-1.)>2.e-5:raise ValueError('Hydrogen orbital embedding failed orthonormality.')
    space=case['reference_space'];total=sum(atomic_numbers[x] for x in case['symbols'])
    closed=(total-space['active_electrons'])//2
    active=list(range(closed,closed+space['active_orbitals']))
    radical=max(active,key=lambda i:amplitudes[i]**2)
    # Degenerate singly occupied orbitals may rotate between fragments. Retain
    # their full span rather than demanding a localized individual natural MO.
    frontier=[i for i in active if abs(d.occupations[i]-1.)<.02]
    hydrogen_projection=float(sum(amplitudes[i]**2 for i in frontier))
    hydrogen_electrons=float(np.dot(amplitudes**2,d.occupations))
    if (hydrogen_projection<.99 or abs(hydrogen_electrons-1.)>.02
            or not frontier or len(frontier)%2!=1):
        raise ValueError('Separated reference lacks a verified neutral-H odd-electron frontier subspace.')
    remaining=sorted(set(active)-set(frontier),key=lambda i:-d.occupations[i])
    occupied_count=(space['active_electrons']-len(frontier))//2
    occupied,virtual=remaining[:occupied_count],remaining[occupied_count:]
    return dict(radical=radical+1,frontier=[i+1 for i in frontier],
                occupied=[i+1 for i in occupied],virtual=[i+1 for i in virtual],
                hydrogen_projection=hydrogen_projection,hydrogen_electrons=hydrogen_electrons,
                occupations=d.occupations[active].tolist(),
                source_sha256=file_sha256(Path(reference['directory'])/(reference['name']+'.natural.orbitals.xml')))


def candidate_config(case, bases, sources, indices, title, anchors, roots=1):
    # Each basis retains its own hash-bound native natural-orbital indices.
    e=len(sources[bases[0]]['partition']['frontier'])+2*len(set(indices[bases[0]]) & set(sources[bases[0]]['partition']['occupied']))
    o=len(indices[bases[0]])
    cfg=base_config(title,case['symbols'],case['geometry'],bases)
    cfg.update(probe_active_space=True,probe_method='casscf',orbital_initialization='propagated',anchors=anchors,
               correlation_checks=[dict(method=f'caspt2-f12({e},{o})',shift=.1,ipea=.25,reference_policy='state_average')])
    # One configuration per basis keeps the orbital partition explicit.
    if len(bases)!=1:raise ValueError('Candidate partitions require a basis-specific configuration.')
    basis=bases[0];s=sources[basis];space=case['reference_space']
    partition=prior.partition(sum(atomic_numbers[x] for x in case['symbols']),space,e,indices[basis])
    cfg['probe_definition']=dict(active_electrons=e,active_orbitals=o,state_model=state(roots),orbital_partition=partition,
        chemical_note='Hydrogen SOMO identified by an independent atomic orbital projection, plus explicitly enumerated fragment valence pair(s). Candidate only; final chemical/state/subspace review required.')
    cfg['initial_records']={}
    seen=set()
    for a in anchors:
        chain=a.get('orbital_chain',a['name'])
        if chain not in seen:
            cfg['initial_records'][cfg['bases'][0]['name']+'__'+a['name']]=dict(path=s['job']['wfu'],sha256=s['job']['wfu_sha256'],basis=basis)
        seen.add(chain)
    return cfg


def run_manifests(manifests, cap):
    try:
        recon.run_many(manifests,max_jobs=cap,keep_going=True,allocation_retries=2)
    except recon.ProbeFailure as error:
        return str(error)
    return None


def distinct_pair_trials(case, scores, jobs, bases):
    """Avoid full-path repetition of already identical optimized screen spaces."""
    ordered=sorted((r for r in scores if r['name'].startswith('pair_')),
                   key=lambda r:(r['max_potential_difference'],r['name']))
    representatives=[];equivalent={};cache={}
    total=sum(atomic_numbers[s] for s in case['symbols'])
    def orbitals(job):
        path=Path(job['directory'])/(job['name']+'.natural.orbitals.xml')
        if path not in cache:cache[path]=read_orbitals(path)
        return cache[path]
    for row in ordered:
        for prior_row in representatives:
            same=True
            for basis,point in itertools.product(bases,('r30','mid','near')):
                left=jobs[(f'h{case["id"]}_{prior_row["name"]}_{basis.replace("-","_")}',basis,point)]
                right=jobs[(f'h{case["id"]}_{row["name"]}_{basis.replace("-","_")}',basis,point)]
                if abs(energy(left)-energy(right))>1.e-7:
                    same=False;break
                try:
                    space=left['probe_active_space'];closed=(total-space['active_electrons'])//2
                    indices=list(range(closed,closed+space['active_orbitals']))
                    overlaps=projection(orbitals(left),orbitals(right),indices,indices)
                except ValueError:
                    same=False;break
                if min(overlaps)<1.-1.e-6:
                    same=False;break
            if same:
                equivalent[row['name']]=prior_row['name'];break
        else:
            representatives.append(row)
    return representatives[:3],equivalent


def run(output):
    output=Path(output).resolve(strict=True);record=read(output/'campaign.json');plan=read(output/'plan.json')
    if record['preview'] or run_policy.find(output)!=Path(record['policy']):raise ValueError('Execution requires the explicit native pool.')
    if file_sha256(HERE/'plan.json')!=record['plan_sha256']:raise ValueError('Campaign plan changed.')
    policy=Path(record['policy']);holds=[];all_manifests=list(record['manifests'])
    with recon._manifest_lock(output):
        print('Running hash-bound references and mandatory endpoint/state checks.',flush=True)
        failure=run_manifests(record['manifests'],plan['max_jobs'])
        if failure:holds.append(failure)
        initial=completed_jobs(record['manifests'])
        for case in plan['h_loss_cases']:
            cid=case['id'];sources={};variants={};screen_manifests=[]
            try:
                for basis in plan['bases']:
                    ref=initial[(f'h{cid}_reference',basis,'r30')];atom=initial[(f'h{cid}_hydrogen',basis,'hydrogen')]
                    sources[basis]=dict(job=ref,partition=hydrogen_partition(case,ref,atom))
                counts={(len(s['partition']['frontier']),len(s['partition']['occupied']),len(s['partition']['virtual'])) for s in sources.values()}
                if len(counts)!=1:raise ValueError('Reference partition dimensions differ between bases.')
                frontier_count,occupied,virtual=next(iter(counts))
                baseline='single' if frontier_count==1 else 'open_frontier'
                variants[baseline]={b:sources[b]['partition']['frontier'] for b in plan['bases']}
                for i,j in itertools.product(range(occupied),range(virtual)):
                    variants[f'pair_{i}_{j}']={b:sorted(sources[b]['partition']['frontier']+[sources[b]['partition']['occupied'][i],sources[b]['partition']['virtual'][j]]) for b in plan['bases']}
                points={a['name']:a for a in h_anchors(case)}
                screen=[points[k] for k in ('r30','mid','near')]
                for name,indices in variants.items():
                    for b in plan['bases']:
                        cfg=candidate_config(case,[b],sources,indices,f'h{cid}_{name}_{b.replace("-","_")}',screen)
                        screen_manifests.append(prepare(output/'screen',cfg,policy))
                immutable(output/f'h{cid}_partitions.json',dict(sources=sources,variants=variants))
                print(f'H-loss case {cid}: screening {len(variants)} explicit odd-electron spaces in both bases.',flush=True)
                failure=run_manifests(screen_manifests,plan['max_jobs'])
                if failure:holds.append(failure)
                all_manifests.extend(screen_manifests);screen_jobs=completed_jobs(screen_manifests)
                scores=[];screen_issues=[]
                for name in variants:
                    differences=[];far_identity={}
                    try:
                        for b in plan['bases']:
                            title=f'h{cid}_{name}_{b.replace("-","_")}'
                            far_job=screen_jobs[(title,b,'r30')]
                            candidate_case=dict(case,reference_space=far_job['probe_active_space'])
                            far_identity[b]=hydrogen_partition(candidate_case,far_job,initial[(f'h{cid}_hydrogen',b,'hydrogen')])
                            far=energy(far_job);ref_far=energy(initial[(f'h{cid}_reference',b,'r30')])
                            for p in ('near','mid'):
                                differences.append(KCAL*((energy(screen_jobs[(title,b,p)])-far)-(energy(initial[(f'h{cid}_reference',b,p)])-ref_far)))
                    except (KeyError,ValueError) as error:
                        screen_issues.append(dict(name=name,reason=str(error)));continue
                    scores.append(dict(name=name,max_potential_difference=max(map(abs,differences)),differences_kcal_mol=differences,far_identity=far_identity))
                ranked,equivalent=distinct_pair_trials(case,scores,screen_jobs,plan['bases'])
                if not ranked:raise ValueError('No complete pair trial; inspect native records.')
                complete_names={r['name'] for r in scores}
                names=([baseline] if baseline in complete_names else [])+[r['name'] for r in ranked]
                selected={name:variants[name] for name in names}
                for count in range(2,len(ranked)+1):
                    union={b:sorted(set().union(*(variants[r['name']][b] for r in ranked[:count]))) for b in plan['bases']}
                    if union not in selected.values():selected[f'union_{count}']=union
                immutable(output/f'h{cid}_followup_trials.json',dict(scores=scores,selected=selected,screen_issues=screen_issues,equivalent_screen_spaces=equivalent,
                    meaning='Diagnostic allocation based on matched own-asymptote potentials. This ranking does not approve an active space.',production_model_accepted=False))
                print(f'H-loss case {cid}: testing full paths for {list(selected)}.',flush=True)
                follow=[]
                for name,indices in selected.items():
                    for b in plan['bases']:
                        cfg=candidate_config(case,[b],sources,indices,f'h{cid}_{name}_full_{b.replace("-","_")}',h_anchors(case))
                        # Adopt exact completed screen points and seed independent new branches from their far record.
                        if name in variants:
                            title=f'h{cid}_{name}_{b.replace("-","_")}';reused={k:screen_jobs[(title,b,k)] for k in ('r30','mid','near') if (title,b,k) in screen_jobs}
                            # Store source identities for review. Only new geometries/return checks are submitted.
                            cfg['anchors']=[a for a in cfg['anchors'] if a['name'] not in reused]
                            cfg['initial_records']={}
                            seed=reused['r30'];seen=set()
                            for a in cfg['anchors']:
                                chain=a.get('orbital_chain',a['name'])
                                if chain not in seen:cfg['initial_records'][cfg['bases'][0]['name']+'__'+a['name']]=dict(path=seed['wfu'],sha256=seed['wfu_sha256'],basis=b)
                                seen.add(chain)
                            e=cfg['probe_definition']['active_electrons'];o=cfg['probe_definition']['active_orbitals']
                            cfg['probe_definition']['orbital_partition']=dict(source_active_electrons=e,source_active_orbitals=o,close_lowest=0,discard_highest=0)
                            immutable(output/f'h{cid}_{name}_{b}_reused.json',reused)
                        follow.append(prepare(output/'followup',cfg,policy))
                # Test state-average dependence for the leading pair and its full retained union.
                largest=max(selected,key=lambda name:len(selected[name][plan['bases'][0]]))
                for name in dict.fromkeys((ranked[0]['name'],largest)):
                    for b in plan['bases']:
                        cfg=candidate_config(case,[b],sources,selected[name],f'h{cid}_{name}_twoD_{b.replace("-","_")}',screen,roots=2)
                        follow.append(prepare(output/'state_support',cfg,policy))
                failure=run_manifests(follow,plan['max_jobs'])
                if failure:holds.append(failure)
                all_manifests.extend(follow)
            except (KeyError,ValueError) as error:
                holds.append(f'h{cid}: {error}')
        atomic_json(output/'completion.json',dict(schema=1,holds=holds,manifests=all_manifests,
                    status='native_complete_review_required' if not holds else 'native_partial_review_required',production_model_accepted=False))
        print(json.dumps(dict(completion=str(output/'completion.json'),holds=holds,production_model_accepted=False)),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['stage','run']);parser.add_argument('output',type=Path)
    parser.add_argument('--pool',type=Path);parser.add_argument('--source-map',type=Path)
    args=parser.parse_args()
    if args.action=='stage':print(json.dumps(stage(args.output,read(args.source_map) if args.source_map else None,args.pool),indent=2))
    else:run(args.output)
