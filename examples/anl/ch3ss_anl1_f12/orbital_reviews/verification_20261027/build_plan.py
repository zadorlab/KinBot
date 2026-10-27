"""Build the finite, example-specific panel from the verified native export.

No scheduler access, orbital approval, or production input modification.
"""
import argparse
from copy import deepcopy
import json
import gzip
from pathlib import Path
import numpy as np
from kinbot import molpro_reconnaissance as recon
from kinbot.execution import file_sha256

HERE=Path(__file__).resolve().parent
ARCHIVE='15b6ef48db20e9aa97b901a3e15759003728bc3bea585865cffcef4c90f3a256'
# Source-natural indices, not a proposed automatic selection algorithm.
CONTROLS={
 1: [('sh_sigma_pair',7,[10,11,12,13,14,15],'Add the identified SH sigma/antibonding pair.')],
 2: [('h2s_acceptor',5,[8,12,13,14,15],'Add the identified H2S antibonding direction.'),
     ('h2s_pair',7,[8,11,12,13,14,15],'Add its occupied H2S partner as a genuine pair control.')],
 6: [('without_h2_pair',5,[19,20,21,22],'Determine whether the H2 sigma/antibonding pair is necessary.')],
 7: [('frontier_only',5,[19,20,21,22],'Test whether the occupied/antibonding pair beyond the SH/frontier manifold is necessary.')],
 8: [('sigma_pair',7,[18,19,20,21,22,23],'Test the identified additional valence sigma/antibonding pair.')],
 9: [('occupied_union',7,[18,19,20,21,22],'Test the next sulfur occupied direction while retaining the three radical frontier directions.')],
 10:[('occupied_union',7,[18,19,20,21,22],'Test the adjacent sulfur occupied direction; reject any frontier substitution.')],
 11:[('donor',5,[19,20,21,22],'Test sulfur donation into the three radical frontier directions at this outer endpoint.')],
 12:[('donor',5,[19,20,21,22],'Test sulfur donation at this independently mapped outer endpoint.')],
 13:[('occupied_union',7,[18,19,20,21,22],'Test the adjacent occupied direction without losing either SH pi component.')],
 14:[('pair_union',5,[19,20,21,22,23],'Two occupied valence directions, H SOMO and two external trial directions; inspect virtual character explicitly. This is not the previous one-hole CAS(3,2).')],
 15:[('pair_union',5,[19,20,21,22,23],'Own-endpoint two-donor/two-external pair control; virtuals are trial directions, not assumed valence orbitals.')],
 16:[('pair_union',5,[19,20,21,22,23],'Own-endpoint two-donor/two-external pair control; confirm resulting valence character.')],
 17:[('occupied_union',7,[18,19,20,21,22],'Test omitted sulfur donor; keep both SH pi components and the molecular frontier pair.')],
 18:[('occupied_union',7,[18,19,20,21,22],'Same chemical test at this distinct parent/atom mapping.')],
 19:[('sigma_pair',7,[18,19,20,21,22,23],'Identified molecular C-S occupied/antibonding pair beyond CH2/HSS frontier.')],
 20:[('donor',7,[18,19,20,21,22],'Test bound-sulfur donor while preserving all three departing-atom S p directions.')],
 21:[('occupied_union',7,[18,19,20,21,22],'Test second occupied sulfur direction beyond the already necessary bound-sulfur donor.')],
 22:[('occupied_union',9,[17,18,19,20,21,22],'Test the adjacent molecular valence sigma direction; preserve atomic p completeness.')],
 24:[('bound_s_donor',5,[18,20,21,22],'Open bound-S2 donor 18, not distal-S3 donor 19, at the reverse product endpoint.')],
 25:[('h2s_donor',5,[19,20,21,22],'Open the actual H2S donor 19. Its optimized retention must be checked, not inferred from fixed CI.')],
 27:[('frontier_only',3,[20,21,22],'Matched actual CASPT2-F12 reference for the donor effect.'),
     ('occupied_union',7,[18,19,20,21,22],'Bound and distal sulfur occupied union: test whether the proposed smaller space survives and is sufficient.')],
 28:[('occupied_union',7,[18,19,20,21,22],'Repeat the optimized union only with actual correlated potential and transport; prior substituted solution is not accepted evidence.')],
 29:[('pair_union',5,[19,20,21,22,23],'Two occupied valence directions plus H SOMO and two external trial directions; verify valence character and pair participation.')],
 33:[('h2s_donor',5,[19,20,21,22],'Test the actual H2S donor rather than the previously tested additional HCS occupied orbital.')],
}
H_EXITS={14,15,16,29}


def read(p):return json.loads(Path(p).read_text())

def build(root, manifest_path, seed_index):
    assert read(root/'VERIFIED.json')['archive_sha256']==ARCHIVE
    assert not read(root/'VERIFIED.json')['issues']
    manifest=read(manifest_path); bysource={f['source']:f for f in manifest['files']}
    bymember={f['member']:f for f in manifest['files']}
    seeds=read(seed_index);ledger=read(HERE/'prior_selection_ledger.json')
    proof={}
    def capture(member):
        f=bymember[member]
        if member not in proof:
            assert file_sha256(root/member)==f['sha256'],member
            proof[member]={k:f[k] for k in ('member','source','sha256','bytes') if k in f}
        return member
    native=read(HERE/'production_native_audit.json') if (HERE/'production_native_audit.json').exists() else json.load(gzip.open(HERE/'production_native_audit.json.gz','rt'))
    cases=[];coverage=[]
    for i,c in enumerate(ledger['channels']):
        status='excluded' if i in (3,23,30,31,32) else ('existing_native_review' if i in (0,4,5) else 'verification_staged')
        coverage.append(dict(reaction=c['reaction'],case=c['case'],selection=c['r8_selection'],action=status,
                             evidence=c.get('evidence'),ledger_row=i,verification_case=25 if i==26 else i))
        if i not in CONTROLS:continue
        model=deepcopy(c['models'][0]) if i!=27 else read(HERE/'unresolved_ch_hssh.json')
        sources={}
        for basis in ('avdz-f12','avtz-f12'):
            chosen=(model.get('basis_seeds') or {b:v['orbital_selection'] for b,v in model['state_model']['basis_models'].items()})[basis]
            h=chosen.get('wfu_sha256',chosen.get('sha256'));s=seeds[h];j=s['job'];folder=Path(s['manifest']).parent/j['name']
            capture(s['manifest']);cfgmember=bysource[s['source']]['member'];capture(cfgmember);cfg=read(root/cfgmember)
            if cfg.get('kinbot_input'):cfg['kinbot_input']=str(root/bysource[cfg['kinbot_input']]['member'])
            anchors=recon._coordinates(cfg,root/cfgmember) if 'coordinates' in cfg else cfg['anchors']
            anchor=next(a for a in anchors if a['name']==j['anchor'])
            assert 'bond' in anchor and 'moving_atoms' in anchor
            wf=capture(str(folder/(j['name']+'.wfu')));out=capture(str(folder/(j['name']+'.out')))
            assert proof[wf]['sha256']==h and proof[out]['sha256']==j['output_sha256']
            for ext in ('.inp','.natural.orbitals.xml'):
                member=str(folder/(j['name']+ext));capture(member)
            sources[basis]=dict(wfu=wf,output=out,configuration=cfgmember,manifest=s['manifest'],
                geometry=j['geometry_angstrom'],symbols=j['symbols'],electrons=read(root/s['manifest'])['electron_count'],
                source_space=j['probe_active_space'],coordinate={k:anchor[k] for k in ('bond','moving_atoms')},
                selected_indices=chosen['active_indices'],job=j['name'])
        sm=deepcopy(model['state_model'].get('basis_models',{}).get('avdz-f12',model['state_model']));sm.pop('orbital_selection',None)
        other=deepcopy(model['state_model'].get('basis_models',{}).get('avtz-f12',sm));other.pop('orbital_selection',None)
        assert sm==other, 'Review basis-dependent state definitions explicitly.'
        # Correlated state averages and root support are the existing choices.
        sm['correlation_options']={**sm.get('correlation_options',{}),'reference_policy':'state_average',
            'max_iterations':2000,'max_reference_iterations':2000,'max_reference_vectors':60}
        sm.setdefault('casscf_options',dict(solver='so_sci',max_iterations=999,exact_diagonalization=True,
            quasi_newton=True,energy_convergence=1.e-8,gradient_convergence=1.e-5))
        active=model['active_space'];variants=[dict(name='candidate',electrons=active['active_electrons'],
            indices={b:s['selected_indices'] for b,s in sources.items()},purpose=model.get('rationale',model.get('review_note','Existing selected active space.')))]
        for name,e,indices,note in CONTROLS[i]:variants.append(dict(name=name,electrons=e,indices={b:indices for b in sources},purpose=note))
        case=dict(id=i,reactions=[c['reaction']],case=c['case'],sources=sources,variants=variants,
            state_model=sm,reactive_states=deepcopy(model.get('reactive_states',model.get('states'))),h_exit=i in H_EXITS,
            status='unresolved_partition' if i==27 else 'candidate_requires_verification',production_model_accepted=False)
        assert case['reactive_states'], 'Retain the explicitly specified reactive-state definition.'
        if i==27:
            case['far_sources']={}
            for b,near_source in sources.items():
                far=next(v for v in seeds.values() if v['manifest']==near_source['manifest'] and v['job']['basis']==b and v['job'].get('distance_angstrom',0)>29)
                j=far['job'];folder=Path(far['manifest']).parent/j['name']
                files={}
                for ext in ('.wfu','.out','.inp','.natural.orbitals.xml'):
                    files[ext]=capture(str(folder/(j['name']+ext)))
                case['far_sources'][b]=dict(wfu=files['.wfu'],geometry=j['geometry_angstrom'],files=files,
                    note='Far occupied 18/19 are mixed HSSH sulfur combinations, not individually labeled bound/distal donors. Test both starts independently and follow inward.')
        if i in (1,2):
            case['reuse_relaxations']={}
            near=float(np.linalg.norm(np.array(sources['avdz-f12']['geometry'])[0]-np.array(sources['avdz-f12']['geometry'])[1]))
            for label,distance in [('near',near),('mid',max(3.5,near+.75))]:
                candidates=[g for g in native['channels'][i]['graphs'] if g['statuses'].get('relaxed_scan')=='complete' and abs(g['distance']-distance)<.151]
                g=min(candidates,key=lambda g:abs(g['distance']-distance))
                folder=Path(g['state']['path']).parent;task=folder/'tasks/relaxed_scan'
                workflow=read(root/folder/'workflow.json');execution=read(root/task/'execution.json')
                assert execution['status']=='executed' and execution['artifacts']['final.xyz']==file_sha256(root/task/'final.xyz')
                constraints=workflow['tasks'][0]['optimizer']['constraints']
                assert not constraints.get('minimum_bonds') and constraints['fixed_bonds']==[[0,1]]
                members=[str(folder/'workflow.json'),str(folder/'state.json'),str(task/'execution.json'),str(task/'final.xyz')]
                for member in members:capture(member)
                case['reuse_relaxations'][label]=dict(geometry=str(task/'final.xyz'),evidence=members,distance=g['distance'],
                    reason='Completed native CASPT2-F12(5,4)/avdz-f12 geometry, fmax=.0005, same sole separation constraint. Nearest completed mid-grid point is used explicitly; no new optimization.')
        if i==25:case['reactions'].append(ledger['channels'][26]['reaction'])
        if i in (1,7):
            case['reuse_fixed']=[]
            prefix=('runs/r7_ch2sh_cs_orbital_validation/' if i==1 else 'runs/r7_progress_review_20261024/ch2_sh_s_followup/candidates/')
            for saved in seeds.values():
                j=saved['job'];sp=j.get('probe_active_space') or {}
                if not saved['manifest'].startswith(prefix) or j.get('status')!='complete' or j.get('validation_errors'):continue
                variant=next((v for v in variants if (v['electrons'],len(v['indices'][j['basis']]))==(sp.get('active_electrons'),sp.get('active_orbitals'))),None)
                if variant is None:continue
                expected=[dict(method=f"caspt2-f12({sp['active_electrons']},{sp['active_orbitals']})",shift=.1,ipea=.25,reference_policy='state_average')]
                if j['correlation_checks']!=expected or j['state_model']['state_average']!=sm['state_average']:continue
                folder=Path(saved['manifest']).parent/j['name'];capture(saved['manifest'])
                files={}
                for suffix in ('.wfu','.out','.inp','.natural.orbitals.xml'):
                    member=str(folder/(j['name']+suffix));capture(member);files[suffix]=member
                assert proof[files['.out']]['sha256']==j['output_sha256'] and proof[files['.wfu']]['sha256']==j['wfu_sha256']
                case['reuse_fixed'].append(dict(variant=variant['name'],basis=j['basis'],job=j,files=files,
                    reason='Previously completed matched CASPT2-F12 protocol, spin/root average, basis and active-space trial from this channel. Retained as baseline; new backward/angular checks remain required.'))
        cases.append(case)
    return dict(schema=1,archive_sha256=ARCHIVE,source_files=proof,coverage=coverage,cases=cases,
        protocol=dict(bases=['avdz-f12','avtz-f12'],method='caspt2-f12',shift=.1,ipea=.25,
                      reference_policy='state_average',cbs_power=3.7,max_jobs=23),
        no_production_models_accepted=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('manifest',type=Path);p.add_argument('seed_index',type=Path)
    a=p.parse_args();plan=build(a.root,a.manifest,a.seed_index)
    (HERE/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    print('CASES',len(plan['cases']),'CHANNELS',len(plan['coverage']),'SOURCE_FILES',len(plan['source_files']))
