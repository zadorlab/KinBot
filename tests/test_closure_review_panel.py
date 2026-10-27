"""Simulated complete detached successor orchestration, not native QC approval."""
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.db import connect

from kinbot import requested_orbital_recon as review, run_policy, run_walltime
from kinbot.execution import atomic_json, file_sha256
from test_state_support_panel import (
    support_tree, support_plan, combined, spectrum, continuation_tree, resolution_tree,
    resolution_plan, tree, followup_plan, plan, original_tree, snapshot,
)

BASE=Path(__file__).resolve().parents[1]/'examples/anl/orbital_review_batches'
PLAN=BASE/'r7_closure_20261027'


@pytest.fixture
def closure(support_tree, support_plan, tmp_path):
    production,_=support_tree
    pool=tmp_path/'pool';pool.mkdir();run_policy.update(pool,max_jobs=3,exclude=['node9'],fair_submission=False)
    opts=dict(shared_pool=True,resource_root=pool)
    review.batch(production,support_plan/'screens.json',pool/'state_support',**opts)
    additional=BASE/'r7_additional_requests_20261027'
    evidence=review._read(additional/'captured_geometry.json')
    worker=production/'ch3ss'/evidence['parent'];network=[]
    for item in review._read(additional/'screens.json')['screens']:
        recipe=review._read(additional/item['recipe']);q=recipe['request_match']
        atomic_json(worker/'vrctst/orbital_review'/q['reaction']/'request.json',dict(q,schema=1,status='needs_review'))
        if 'geometry_source' in recipe:
            network.append(dict(name=q['reaction'],product_complex=recipe['geometry_source']['product_complex']))
    atomic_json(worker/'accepted_l2_network.json',dict(network=dict(reactions=network)))
    db=connect(worker/'kinbot.db')
    for row in evidence['endpoint_rows']:
        db.write(Atoms('CSSHHH',positions=row['positions']),name=row['name'],data=dict(status='normal'))
    review.batch(production,additional/'screens.json',pool/'additional_requests',**opts)
    audit=review._read(BASE/'r7_detached_review_20261027/audit.json')
    native=audit['jobs']
    for panel in ('state_support','additional_requests'):
        for system in review._saved(pool/panel)['systems']:
            path=Path(system['manifest']);state=review._read(path);screen=path.parents[1].name
            for job in state['jobs']:
                d=Path(job['directory']);out=d/(job['name']+'.out')
                if not job['probe_active_space']:
                    row=next(x for x in native if x['panel']==panel and x['screen']==screen and x['job']==job['name'])
                    orbitals=row['orbitals']
                    first=min(int(o['index'].split('.')[0]) for o in orbitals)
                    orbitals=[dict(index=f'{i}.1',occupation=2.,energy_Eh=-20.)
                              for i in range(1,first)]+orbitals
                    out.write_text('PROGRAM * RHF-SCF\n ELECTRON ORBITALS\n\n'+''.join(
                        f" {o['index']} {o['occupation']:.5f} {o['energy_Eh']:.5f} 1 1 px 0.50000\n"
                        for o in orbitals)+'\nEND\n')
                else:out.write_text('Synthetic natural-CAS completion for transport tests.\n')
                wfu=out.with_suffix('.wfu');wfu.write_bytes(('mock:'+screen+job['name']).encode())
                job.update(status='complete',wfu=str(wfu),wfu_sha256=file_sha256(wfu),output_sha256=file_sha256(out))
            state['status']='complete_review_required';atomic_json(path,state)
    dest=tmp_path/'recipes';dest.mkdir()
    for p in PLAN.glob('*.json'):(dest/p.name).write_bytes(p.read_bytes())
    for entry in review._read(dest/'screens.json')['screens']:
        path=dest/entry['recipe'];r=review._read(path)
        selected=review._selected_system(review._saved(pool/entry['source']),r)
        jobs={j['name']:j for j in review._read(selected['manifest'])['jobs']}
        for candidate in r['candidates']:
            for a in candidate['anchors']:
                if not a.get('propagate'):
                    a['seed_sha256']={b:jobs[b+'__'+a['seed']]['wfu_sha256'] for b in ('adz','atz')}
        atomic_json(path,r)
    return production,pool,dest,pool/'closure',opts


def test_all_90_preflight_resume_partitions_and_source_preservation(closure,monkeypatch):
    production,pool,plan,out,opts=closure
    monkeypatch.setattr(review.subprocess,'run',lambda *a,**k:pytest.fail('Staging cannot submit'))
    before=snapshot(production);sources={p:v for p,v in snapshot(pool).items() if p.name!=run_policy.NAME}
    assert review.batch(production,plan/'screens.json',out,preflight=True,**opts)['tasks']==90
    assert not out.exists() and snapshot(production)==before
    assert review.batch(production,plan/'screens.json',out,**opts)['tasks']==90
    casci=propagated=0
    for s in review._saved(out)['systems']:
        cfg=review.recon._load(s['configuration']);jobs=review._read(s['manifest'])['jobs']
        for j in jobs:
            text=(Path(j['directory'])/(j['name']+'.inp')).read_text()
            assert 'rs2' not in text and not j.get('job_id')
            if cfg['probe_method']=='casci':
                casci+=1;assert '{multi,cas-ci;' in text and 'noextra\n' in text
                assert not j['orbital_predecessor']
                if s['title'].endswith('fixed_18'):assert 'rotate,18.1,19.1,0' in text
                else:assert 'rotate,' not in text
                original=review._read(Path(s['manifest']).parents[1]/'source_manifest.json')
                old=next(x for x in original['jobs'] if x['name']==j['name'])
                assert np.array_equal(old['geometry_angstrom'],j['geometry_angstrom'])
            if j['orbital_predecessor']:
                propagated+=1;assert 'diab,2140.2' in text
            assert max(x['active_orbitals'] for x in [j['probe_active_space']])<=5
    assert (casci,propagated)==(20,22)
    staged=snapshot(out)
    assert review.batch(production,plan/'screens.json',out,resume=True,**opts)['tasks']==90
    assert snapshot(out)==staged and snapshot(production)==before
    assert all(p.read_bytes()==v for p,v in sources.items())


@pytest.mark.parametrize('change',['wrong_model','escape','changed_hash','closed_swap','no_resource'])
def test_bad_successor_rejected_before_any_staging(closure,change):
    production,pool,plan,out,opts=closure
    path=plan/'new_cs.json';r=review._read(path)
    if change=='wrong_model':r['source_model']='nonexistent'
    elif change=='changed_hash':r['candidates'][0]['anchors'][0]['seed_sha256']['adz']='0'*64
    elif change=='closed_swap':r['candidates'][0]['probe_definition']['orbital_partition']['closed_orbital_swaps']=[[19,20]]
    elif change=='escape':
        p=plan/'screens.json';v=review._read(p);v['screens'][-1]['source']='../production';atomic_json(p,v)
    else:opts=dict(shared_pool=True)
    atomic_json(path,r);before=snapshot(production)
    with pytest.raises((ValueError,FileNotFoundError)):
        review.batch(production,plan/'screens.json',out,**opts)
    assert not out.exists() and snapshot(production)==before


def test_successor_scheduler_preserves_other_work_and_repeated_invocation(closure,monkeypatch):
    production,pool,plan,out,opts=closure
    review.batch(production,plan/'screens.json',out,**opts);before=snapshot(production)
    launched=[];running={};done=set()
    monkeypatch.setattr(run_policy,'_queue',lambda root:running.copy())
    monkeypatch.setattr(run_policy,'slurm_job_record',lambda j:{'JobId':j,'JobState':'COMPLETED'} if j in done else {})
    monkeypatch.setattr(run_walltime,'partition_walltime',lambda p:'31-00:00:00')
    def scheduler(command,**kw):
        if command[0]=='bash':return SimpleNamespace(returncode=0,stdout='/mock/molpro\n',stderr='')
        if command[0]=='squeue':return SimpleNamespace(returncode=0,stdout='\n'.join([*running,'unrelated']),stderr='')
        assert command[0]=='sbatch' and len(running)<3
        assert Path(kw['cwd']).is_relative_to(out)
        job=str(2000+len(launched));launched.append(str(kw['cwd']));running[job]=1
        return SimpleNamespace(returncode=0,stdout=job+'\n',stderr='')
    monkeypatch.setattr(review.subprocess,'run',scheduler)
    def audit(path,**kw):
        state=review._read(path)
        for j in state['jobs']:
            j['status']='complete' if j.get('job_id') in done else 'incomplete'
            if j['status']=='complete':
                p=Path(j['directory'])/(j['name']+'.wfu');p.write_bytes(b'mock next checkpoint')
                j.update(wfu=str(p),wfu_sha256=file_sha256(p))
        state['status']='complete_review_required' if all(j['status']=='complete' for j in state['jobs']) else 'incomplete'
        atomic_json(path,state);return state
    monkeypatch.setattr(review.recon,'_audit',audit)
    monkeypatch.setattr(review.recon,'_recover_slurm_startup_failures',lambda *a:None)
    monkeypatch.setattr(review.recon,'_retry_allocation_failures',lambda p,m,*a:m)
    def finish(_):done.update(running);running.clear()
    monkeypatch.setattr(review.recon.time,'sleep',finish)
    assert review.run(out)['status']=='complete_review_required'
    assert len(launched)==len(set(launched))==90
    review.run(out);assert len(launched)==90 and snapshot(production)==before
    assert review._read(pool/'.kinbot_policy_admissions.json')=={}
