"""Bind explicit unresolved hypotheses to the completed detached export.

This example writes recipes, never runs QC and never accepts a CAS selection.
"""
from copy import deepcopy
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
OPTS = dict(solver='so_sci', quasi_newton=True, exact_diagonalization=True,
            energy_convergence=1e-8, gradient_convergence=1e-5, max_iterations=999)


def read(p):return json.loads(p.read_text())
def write(p,x):p.write_text(json.dumps(x,indent=2)+'\n')


def model(blocks):
    return dict(canonical_orbitals=True,
                selected_roots={str(m):1 for m,n in blocks},
                state_average=[dict(multiplicity=m,roots=n,weights=[1/n]*n,
                                    block_weight=1/len(blocks)) for m,n in blocks])


def build(export, output):
    pool=export/'runs/ch3s_ch3ss_orbital_reviews_20261027'
    output.mkdir(parents=True,exist_ok=True)
    screens=[]; sources=[]

    def source(panel, name, title=None):
        root=pool/panel/name;launch=read(root/'review_launch.json')
        systems=[s for s in launch['systems'] if title is None or s['title']==title]
        if len(systems)!=1:raise ValueError('Explicit source model required')
        system=systems[0]
        manifest=root/Path(system['manifest']).parent.name/'reconnaissance.json'
        native=read(manifest)
        if any(j['status']!='complete' or j.get('validation_errors') for j in native['jobs']):
            raise ValueError('Unfinished source')
        jobs={j['name']:j for j in native['jobs']}
        for j in jobs.values():
            for path,key in ((Path(j['wfu']),'wfu_sha256'),
                             (Path(j['directory'])/(j['name']+'.out'),'output_sha256')):
                local=pool/path.relative_to('/home/hram/kinbot-runs/ch3s_ch3ss_orbital_reviews_20261027')
                if hashlib.sha256(local.read_bytes()).hexdigest()!=j[key]:raise ValueError(local)
        sources.append(dict(manifest=str(manifest.relative_to(export)),
                            sha256=hashlib.sha256(manifest.read_bytes()).hexdigest()))
        q=read(root/'request.json')
        match={k:q[k] for k in ('parent','reaction','symbols','positions','bond','charge','multiplicity','products')}
        recipe=dict(schema=1,status='ready',max_jobs=3,request_match=match,candidates=[])
        if title:recipe['source_model']=title
        return recipe,jobs

    def anchors(jobs, *, scf=None):
        names=[j.split('__',1)[1] for j in jobs if j.startswith('adz__')]
        far=next(n for n in names if '_30' in n)
        result=[]
        for name in names:
            seed=far if scf else name
            a=dict(source=name,seed=seed,seed_sha256={b:jobs[b+'__'+seed]['wfu_sha256'] for b in ('adz','atz')})
            if scf:a['active_indices']={b:scf for b in ('adz','atz')}
            result.append(a)
        return result

    def candidate(name,e,o,states,points,note,partition=None,method='casscf'):
        definition=dict(active_electrons=e,active_orbitals=o,state_model=states,chemical_note=note)
        if partition:definition['orbital_partition']=partition
        else:definition['orbital_source']='scf'
        return dict(name=name,probe_method=method,orbital_initialization='independent',
                    anchors=points,casscf_options=deepcopy(OPTS) if method=='casscf' else {},
                    correlation_checks=[],probe_definition=definition)

    def partition(e,o,opened=0,swap=None):
        p=dict(source_active_electrons=e,source_active_orbitals=o,close_lowest=0,
               discard_highest=0,open_highest_closed=opened)
        if swap:p['closed_orbital_swaps']=[swap]
        return p

    def save(name,panel,src,recipe):
        write(output/(name+'.json'),recipe)
        screens.append(dict(name=name,case='ch3ss',source_root='resource',
                            source=panel+'/'+src,recipe=name+'.json'))

    # Only these six SCF-only requests remain without a compact CAS execution.
    scf_specs=[
        ('new_csh_sh_a','state_support','new_csh_sh_a',[
            ('frontier',5,4,[(2,2)],[18,19,21,22],
             'Own HCSH pi/pi-star pair and complete SH pi pair. Use the earlier CH2S/SH exit selection as chemical guidance, never its ordinal checkpoint.')]),
        ('new_csh_sh_b','state_support','new_csh_sh_b',[
            ('frontier',5,4,[(2,2)],[18,19,21,22],
             'Own HCSH pi/pi-star pair and complete SH pi pair; independent stereochemical parent and WFU.')]),
        ('stereo_cs','additional_requests','stereo_cs',[
            ('frontier',3,3,[(4,1),(2,1)],[20,21,22],
             'CH2 two open shells plus HSS open shell. Control for the occupied S2 donor; not presumed sufficient for C-S bond formation.'),
            ('donor',5,4,[(4,1),(2,1)],[18,20,21,22],
             'Add the own-source occupied S2 donor 18 to CH2/HSS open-shell frontier 20-22. Earlier parent required this chemical direction; source ordinal differs.')]),
        ('stereo_ss','additional_requests','stereo_ss',[
            ('frontier',5,4,[(4,3),(2,3)],[19,20,21,22],
             'Complete atomic S 3p manifold 19-21 (four electrons) plus CH2SH radical 22. Three roots per spin preserve all atomic components.'),
            ('donor',7,5,[(4,3),(2,3)],[18,19,20,21,22],
             'Explicit occupied bound-sulfur 18 control for possible additional fragment component. This is not an extra atomic S component.')]),
        ('stereo_ch2_hss_exit','additional_requests','stereo_ch2_hss_exit',[
            ('frontier',3,3,[(4,1),(2,1)],[20,21,22],
             'Outer CH2/HSS separation from the verified reverse L2 product endpoint; CH2 two frontier directions plus HSS radical. Follow earlier outer-exit CAS(3,3) evidence.'),
            ('donor',5,4,[(4,1),(2,1)],[18,20,21,22],
             'Same outer exit with occupied S2 donor 18 control. Distinct from bound-parent C-S homolysis; do not assume its larger CAS is necessary here.')]),
        ('stereo_h_exit','additional_requests','stereo_h_exit',[
            ('somo',1,1,[(2,1)],[21],
             'Own departing H5 1s SOMO and closed-shell HCSSH fragment at the verified forward L2 endpoint. Confirm the formal radical minimum as for previously reviewed H exits; no one-hole donor comparison is repeated.')])]
    for name,panel,src,choices in scf_specs:
        recipe,jobs=source(panel,src)
        for suffix,e,o,blocks,indices,note in choices:
            recipe['candidates'].append(candidate(name+'_'+suffix,e,o,model(blocks),
                                                  anchors(jobs,scf=indices),note))
        save(name,panel,src,recipe)

    # Fixed orbital controls cannot silently exchange the named donor with a
    # different inactive orbital. They are diagnostics, not publishable seeds.
    for name,title,blocks in [('new_cs','new_cs_ch_frontier',[(2,2)]),
                               ('new_hcs_h2s_exit','new_hcs_h2s_pi',[(2,1)])]:
        recipe,jobs=source('state_support',name,title)
        points=anchors(jobs)
        for index in (18,19):
            note=('Fixed-orbital CASCI: add source inactive '+str(index)+
                  ' to the identical optimized frontier 20-22. Report the identity and CI participation of this direction at each endpoint; do not equate near/far ordinals with a conserved donor. '
                  'Separates CI-space response from the orbital substitution in the earlier optimized control.')
            recipe['candidates'].append(candidate(name+'_fixed_'+str(index),5,4,model(blocks),points,note,
                                                   partition(3,3,1,[18,19] if index==18 else None),'casci'))
        recipe['candidates'].append(candidate(name+'_occupied_union',7,5,model(blocks),points,
            'Union of both occupied alternatives 18/19 and the reviewed frontier 20-22. Native optimized character and CI must show whether a smaller chemically complete subset can be retained; total-energy lowering is not the acceptance rule.',partition(3,3,2)))
        save(name,'state_support',name,recipe)

    recipe,jobs=source('state_support','new_ss','new_ss_sh_pi')
    recipe['candidates'].append(candidate('new_ss_fixed_sigma',7,5,model([(2,2)]),anchors(jobs),
        'Fixed same-geometry CASCI adds occupied source 18 to the optimized SH plus HCSH frontier 19-22. Measures CI response without changing the original frontier virtual. The prior optimized comparison reorganized that virtual; do not mislabel its energy change an isolated sigma increment.',partition(5,4,1),'casci'))
    save('new_ss','state_support','new_ss',recipe)

    # All six states are retained. Branch identity is inferred from native CI
    # and character along the coordinate, never by demanding a fixed root rank.
    recipe,jobs=source('state_support','csh2_sh_exit','csh2_sh_six_state_frontier')
    states=deepcopy(read(pool/'state_support/csh2_sh_exit/csh2_sh_six_state_frontier.json')['probe_definition']['state_model'])
    points=[]
    for i,distance in enumerate((3.,3.5,4.,5.,6.,8.,12.,20.,30.)):
        a=dict(source='exit_l2',name='out_'+str(distance).replace('.','p'),distance_angstrom=distance,orbital_chain='out')
        if i==0:a.update(seed='exit_l2',seed_sha256={b:jobs[b+'__exit_l2']['wfu_sha256'] for b in ('adz','atz')})
        else:a['propagate']=True
        points.append(a)
    c=candidate('csh2_sh_branch_out',5,4,states,points,
        'SH retention is already established. Correlate the 4D/2Q states through separation: far D3/D4 are singlet-fragment plus SH whereas far D1/D2 coincide with quartet triplet-fragment branches. Record crossings/mixing, not a forced lowest-root assignment.',partition(5,4))
    c['orbital_initialization']='propagated';recipe['candidates'].append(c)
    # Return from the independently converged far seed checks initialization
    # dependence at the interacting anchors, without repeating old failed 2D averages.
    points=[]
    for i,distance in enumerate((8.,4.,3.,2.4857243464237877)):
        a=dict(source='exit_l2',name='back_'+str(i),distance_angstrom=distance,orbital_chain='back')
        if i==0:a.update(seed='exit_r_30p0',seed_sha256={b:jobs[b+'__exit_r_30p0']['wfu_sha256'] for b in ('adz','atz')})
        else:a['propagate']=True
        points.append(a)
    c=candidate('csh2_sh_branch_back',5,4,states,points,
        'Return branch check from independently optimized 30 A seed, retaining the identical six-state manifold. Test path dependence of the near optimized solution.',partition(5,4))
    c['orbital_initialization']='propagated';recipe['candidates'].append(c)
    save('csh2_sh_exit','state_support','csh2_sh_exit',recipe)
    write(output/'screens.json',dict(schema=1,status='ready',max_jobs=3,screens=screens))
    write(output/'source_evidence.json',dict(schema=1,
        archive_sha256='e4c5bbe1dc5346b95db594f41d9d59319173bbc3fc8a0790892adaba2bcfe2a8',sources=sources,
        note='Scientific hypotheses only. No CAS selection or production approval is implied.'))
    print(json.dumps(dict(screens=len(screens),tasks=sum(2*len(c['anchors']) for s in screens
                    for c in read(output/s['recipe'])['candidates']))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('export',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();build(a.export.resolve(strict=True),a.output)
