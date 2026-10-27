"""Read-only numerical/native evidence report. Never publishes a CAS selection."""
import argparse
from copy import deepcopy
import gzip
import importlib.util
import json
from pathlib import Path
import numpy as np
from kinbot.execution import file_sha256
from kinbot.molpro_diagnostics import parse_native_output

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('panel',HERE/'panel.py');panel=importlib.util.module_from_spec(spec);spec.loader.exec_module(panel)
compare_path=HERE.parent/'batches/r7_compact_union_20261025/compare_orbitals.py'
if not compare_path.is_file():compare_path=HERE.parents[2]/'orbital_review_batches/r7_compact_union_20261025/compare_orbitals.py'
if not compare_path.is_file():compare_path=HERE/'compare_orbitals.py'
spec=importlib.util.spec_from_file_location('same_geometry_orbitals',compare_path);orbitals=importlib.util.module_from_spec(spec);spec.loader.exec_module(orbitals)
HARTREE_KCAL=627.5094740631


def native(job):
    d=Path(job['directory']);p=d/(job['name']+'.out');inp=d/(job['name']+'.inp')
    if job.get('status')!='complete':return dict(status=job.get('status'),reason='Native calculation incomplete.')
    if file_sha256(p)!=job['output_sha256'] or file_sha256(inp)!=job['input_sha256']:raise ValueError('Changed native file: '+str(p))
    parsed=parse_native_output(p.read_text(),inp.read_text())
    multi=[s for s in parsed['native_sections'] if s['program']=='MULTI']
    result=dict(status=job['status'],output=str(p),output_sha256=job['output_sha256'],input_sha256=job['input_sha256'],
        normal=parsed['normal_termination_seen'],errors=parsed['native_error_lines'],issues=parsed['issues'],
        geometry=job['geometry_angstrom'],root_states=job['root_states'],state_model=job['state_model'],
        active_space=job['probe_active_space'],state_energies=job.get('state_energies_hartree'),
        correlation=job.get('correlated_energies_hartree'),orbital_xml=job.get('orbital_xml'),
        native_multi=multi,production_model_accepted=False)
    return result


def energy(job,method='caspt2-f12',mult=2,root=1):
    values=[v['energy_hartree'] for v in job.get('correlated_energies_hartree',[]) if v['method'].lower().split('(')[0]==method and v['multiplicity']==mult and v['root']==root]
    return values[0] if len(values)==1 else None


def analyze(output):
    output=Path(output).resolve();record=panel.read(output/'panel.json');plan=panel.read(output/'plan.json')
    systems=record['systems']+panel.read(output/'relaxed_points.json').get('systems',[]) if (output/'relaxed_points.json').exists() else record['systems']
    jobs={};evidence={};comparisons=[];missing=[];same_geometry=[];state_checks=[]
    for system in systems:
        m=panel.read(system['manifest']);alljobs={j['anchor']:j for j in m['jobs']};alljobs.update(system.get('reused',{}))
        for anchor,j in alljobs.items():
            point=anchor if 'anchor' not in system else system['anchor']+'_relaxed'
            key=(system['case_id'],system['variant'],system['basis'],point)
            jobs[key]=j;label='|'.join(map(str,key));evidence[label]=native(j)
            if j.get('status')!='complete':missing.append(label)
    for key,j in jobs.items():
        cid,name,basis,point=key
        if j.get('status')!='complete':continue
        for v in j.get('correlated_energies_hartree',[]):
            method=v['method'].lower().split('(')[0];mult=v['multiplicity'];root=v['root']
            row=dict(case_id=cid,variant=name,basis=basis,point=point,method=method,multiplicity=mult,root=root,energy_hartree=v['energy_hartree'])
            refpoint='near' if point in ('back_near','in_near') else ('r6' if point in ('back6','in_r6') else ('r30' if point in ('r20','in_r30') else None))
            refname=name
            if '_far' in name:refname=name.split('_far',1)[0]
            if refpoint:
                ref=energy(jobs.get((cid,refname,basis,refpoint),{}),method,mult,root)
                if ref is not None:row['same_state_difference_kcal_mol']=(v['energy_hartree']-ref)*HARTREE_KCAL;row['compared_to']=refpoint
            d1=energy(j,method,2,1)
            if d1 is not None:row['gap_from_D1_kcal_mol']=(v['energy_hartree']-d1)*HARTREE_KCAL
            far=energy(jobs.get((cid,name,basis,'r30'),{}),method,mult,root)
            if far is not None:row['relative_to_own_asymptote_kcal_mol']=(v['energy_hartree']-far)*HARTREE_KCAL
            state_checks.append(row)
    for case in plan['cases']:
        cid=case['id'];bybasis={}
        for basis in plan['protocol']['bases']:
            points={key[3] for key in jobs if key[:3]==(cid,'candidate',basis)}
            for variant in case['variants']:
                name=variant['name'];far=jobs.get((cid,name,basis,'r30'),{})
                ef=energy(far)
                for point in sorted(points):
                    j=jobs.get((cid,name,basis,point),{});e=energy(j)
                    if ef is None or e is None:continue
                    # Every angular potential uses that orientation's own asymptote.
                    ref='angle_far_'+point.split('_')[-1] if point.startswith('angle_') and not point.startswith('angle_far_') else 'r30'
                    er=energy(jobs.get((cid,name,basis,ref),{}))
                    if er is None:continue
                    bybasis[(name,basis,point)]=(e-er)*HARTREE_KCAL
                    base=jobs.get((cid,'candidate',basis,point),{});basefar=jobs.get((cid,'candidate',basis,ref),{})
                    eb,ebf=energy(base),energy(basefar)
                    comparisons.append(dict(case_id=cid,variant=name,basis=basis,point=point,
                        relative_kcal_mol=(e-er)*HARTREE_KCAL,
                        difference_from_candidate_kcal_mol=None if eb is None or ebf is None else ((e-er)-(eb-ebf))*HARTREE_KCAL,
                        reference_point=ref))
            for variant in case['variants']:
                name=variant['name'];j=jobs.get((cid,name,basis,'near'))
                if not j or j.get('status')!='complete':continue
                source=case['sources'][basis];left=panel.source_resolver(plan,record['source_export_root'])(source['wfu']).with_suffix('.natural.orbitals.xml')
                right=Path(j['directory'])/(j['name']+'.natural.orbitals.xml')
                e=variant['electrons'];o=len(variant['indices'][basis]);c=(source['electrons']-e)//2
                try:
                    a,atoms_a,basis_a=orbitals.read(left);b,atoms_b,basis_b=orbitals.read(right)
                    if atoms_a!=atoms_b or basis_a!=basis_b:raise ValueError('Not exactly the same geometry/AO basis.')
                    t=b@np.linalg.pinv(a);res=float(np.linalg.norm(t@a-b)/np.linalg.norm(b));orth=float(np.max(abs(t@t.T-np.eye(len(t)))))
                    if res>1.e-9 or orth>1.e-8:raise ValueError('Invalid complete-MO transformation.')
                    selected=np.array(variant['indices'][basis])-1;target=np.arange(c,c+o)
                    block=t[np.ix_(target,selected)]
                    same_geometry.append(dict(case_id=cid,variant=name,basis=basis,left=str(left),right=str(right),
                        left_sha256=file_sha256(left),right_sha256=file_sha256(right),
                        singular_values=np.linalg.svd(block,compute_uv=False).tolist(),
                        source_projection_into_active={str(k+1):float(sum(t[target,k]**2)) for k in selected},
                        reconstruction_error=res,orthogonality_error=orth))
                except (ValueError,OSError) as error:same_geometry.append(dict(case_id=cid,variant=name,basis=basis,unavailable=str(error)))
        for (name,basis,point),v in list(bybasis.items()):
            if basis!='avtz-f12' or (name,'avdz-f12',point) not in bybasis:continue
            dz=bybasis[(name,'avdz-f12',point)];power=plan['protocol']['cbs_power'];ratio=(3/2)**power
            comparisons.append(dict(case_id=cid,variant=name,basis='CBS_D_T',point=point,relative_kcal_mol=(ratio*v-dz)/(ratio-1)))
    return dict(schema=1,archive_sha256=plan['archive_sha256'],evidence=evidence,matched_potentials=comparisons,
        same_geometry_subspaces=same_geometry,state_checks=state_checks,missing=missing,production_model_accepted=False,
        requirements=['Inspect native CI, all supported spin/spatial roots, and fragment identity.',
            'Compare outgoing/return energies and active subspaces at identical geometries.',
            'Inspect 20-to-30 A plateau and the separate angular reference energies.',
            'Compare each CAS with its own asymptote in both bases and CBS; total energy lowering is not evidence of necessity.',
            'H-exit external trial orbitals must be characterized. A diffuse substitution does not establish valence-pair convergence.',
            'MRCI+Q provides independent pair-correlation evidence; it is not a replacement production potential.',
            'No finite panel guarantees perfection or a unique minimal CAS. Any failed scientific criterion remains unresolved.'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--report',required=True,type=Path)
    a=p.parse_args();result=analyze(a.output)
    with gzip.open(a.report,'wt') as f:json.dump(result,f)
    print(json.dumps(dict(report=str(a.report),native_points=len(result['evidence']),incomplete=len(result['missing']),production_model_accepted=False),indent=2))
