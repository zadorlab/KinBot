"""Finite native-verification staging/orchestration, simulated site only."""
import importlib.util
import json
from copy import deepcopy
from pathlib import Path
import numpy as np
import pytest
from ase import Atoms
from ase.io import write
from kinbot import run_policy
from kinbot.execution import atomic_json,file_sha256

BASE=Path(__file__).resolve().parents[1]/'examples/anl/ch3ss_anl1_f12/orbital_reviews/verification_20261027'
spec=importlib.util.spec_from_file_location('verification_panel',BASE/'panel.py')
panel=importlib.util.module_from_spec(spec);spec.loader.exec_module(panel)


@pytest.fixture
def synthetic(tmp_path,monkeypatch):
    plan=panel.read(BASE/'plan.json')
    # Exercise the unresolved donor's three alternatives and a genuine external pair control.
    plan['cases']=[c for c in plan['cases'] if c['id'] in (14,27)]
    needed={s['wfu'] for c in plan['cases'] for s in c['sources'].values()} | {s['wfu'] for c in plan['cases'] for s in c.get('far_sources',{}).values()}
    evidence=tmp_path/'immutable_source';evidence.mkdir();files={}
    for n,member in enumerate(sorted(needed)):
        p=evidence/(str(n)+'.wfu');p.write_text('Synthetic checkpoint; no native scientific evidence. '+str(n))
        files[member]=dict(source=str(p),member=member,sha256=file_sha256(p))
    plan['source_files']=files
    recipe=tmp_path/'recipe';recipe.mkdir();atomic_json(recipe/'plan.json',plan)
    monkeypatch.setattr(panel,'HERE',recipe)
    pool=tmp_path/'pool';pool.mkdir();run_policy.update(pool,max_jobs=23,exclude=['node9'],fair_submission=False)
    return plan,evidence,pool,pool/'verification'


def test_inventory_covers_every_channel_and_preserves_existing_models():
    p=panel.read(BASE/'plan.json');rows=p['coverage']
    assert len(rows)==34 and sum(r['action']=='excluded' for r in rows)==5
    assert sum(r['action']=='existing_native_review' for r in rows)==3
    assert len(p['cases'])==25 and sum(len(c['reactions']) for c in p['cases'])==26
    assert p['no_production_models_accepted']
    selected={r['reaction'] for r in rows if r['action']=='verification_staged'}
    assert selected=={r for c in p['cases'] for r in c['reactions']}
    for c in p['cases']:
        assert c['reactive_states']
        for state in c['reactive_states']:
            support=next(s for s in c['state_model']['state_average'] if s['multiplicity']==state['multiplicity'])
            assert 0<state['root']<=support['roots'] and state['weight']>0
        for b,s in c['sources'].items():
            a=panel.anchors(s)
            assert np.allclose(a[0]['geometry_angstrom'],a[-5]['geometry_angstrom']) # return near
            for v in c['variants']:
                d=panel.partition(s['electrons'],s['source_space'],v['electrons'],v['indices'][b])
                assert d['source_active_electrons']-2*d['close_lowest']+2*d.get('open_highest_closed',0)==v['electrons']
    h=next(c for c in p['cases'] if c['id']==14)
    assert h['variants'][1]['electrons']==5 and len(h['variants'][1]['indices']['avdz-f12'])==5


def test_stage_repeat_hash_failures_and_immutable_sources(synthetic,monkeypatch):
    plan,evidence,pool,out=synthetic
    monkeypatch.setattr(panel.recon.subprocess,'run',lambda *a,**k:pytest.fail('Offline staging must not query/submit.'))
    before={p:p.read_bytes() for p in evidence.iterdir()}
    first=panel.stage(out,pool=pool);record=panel.read(out/'panel.json')
    assert first['fixed_jobs']==134
    inputs={p:p.read_bytes() for p in out.rglob('*.inp')}
    assert panel.stage(out,pool=pool)==first
    assert all(p.read_bytes()==v for p,v in inputs.items())
    assert all(p.read_bytes()==v for p,v in before.items())
    entry=next(s for s in record['systems'] if s['case_id']==27 and s['variant']=='candidate')
    jobs=panel.read(entry['manifest'])['jobs'];deck=(Path(jobs[0]['directory'])/(jobs[0]['name']+'.inp')).read_text()
    assert 'rotate,19.1,18.1,0' in deck and 'occ,22\nclosed,18' in deck
    assert deck.count('{rs2-f12,')==2
    original=next(iter(before));original.write_text('changed')
    with pytest.raises(ValueError,match='Missing/changed source'):
        panel.stage(pool/'blocked',pool=pool)
    assert not (pool/'blocked').exists()


def test_complete_three_phase_owner_resume_and_failed_prerequisite(synthetic,monkeypatch):
    plan,evidence,pool,out=synthetic
    panel.stage(out,pool=pool)
    before={p:p.read_bytes() for p in evidence.iterdir()};new_jobs=[];graph_ids=[]
    def fake_many(paths,**options):
        assert options['max_jobs']==23 and options['allocation_retries']==2 and options['keep_going']
        for p in paths:
            m=panel.read(p)
            for j in m['jobs']:
                if j.get('status')=='complete':continue
                d=Path(j['directory']);o=d/(j['name']+'.out');w=d/(j['name']+'.wfu')
                o.write_text('SIMULATION: no native computation.');w.write_text('SIMULATED result '+j['name'])
                j.update(status='complete',wfu=str(w),wfu_sha256=file_sha256(w),output_sha256=file_sha256(o))
                new_jobs.append((str(p),j['name']))
            atomic_json(p,m)
        return [panel.read(p) for p in paths]
    def fake_advance(g,submit):
        assert submit and run_policy.find(g)==pool/run_policy.NAME
        g=Path(g);state=panel.read(g/'state.json');spec=panel.read(g/'workflow.json')
        if state['tasks']['relaxed_scan']['status']!='complete':
            graph_ids.append(str(g));d=g/'tasks/relaxed_scan';final=d/'final.xyz'
            write(final,Atoms(spec['molecule']['symbols'],positions=spec['molecule']['positions']))
            atomic_json(d/'execution.json',dict(status='executed',artifacts={'final.xyz':file_sha256(final)}))
            state['tasks']['relaxed_scan']['status']='complete';atomic_json(g/'state.json',state)
        return state
    monkeypatch.setattr(panel.recon,'run_many',fake_many)
    monkeypatch.setattr(panel.dispatch,'preflight_until_scheduler_ready',lambda g:None)
    monkeypatch.setattr(panel.dispatch,'advance',fake_advance)
    panel.run(out)
    assert len(graph_ids)==4 and len(new_jobs)==154
    assert panel.read(out/'completion.json')['production_model_accepted'] is False
    identity=list(graph_ids);jobs=list(new_jobs)
    panel.run(out)
    assert graph_ids==identity and new_jobs==jobs
    assert all(p.read_bytes()==v for p,v in before.items())
    system=panel.read(out/'panel.json')['systems'][0];m=panel.read(system['manifest']);m['jobs'][0]['status']='failed';atomic_json(system['manifest'],m)
    with pytest.raises(ValueError,match='prerequisite incomplete'):panel.verified_job(system,'near')


def test_submitting_journal_stops_before_followons_and_preserves_live_work(synthetic,monkeypatch):
    _,evidence,pool,out=synthetic
    panel.stage(out,pool=pool);s=panel.read(out/'panel.json')['systems'][0]
    m=panel.read(s['manifest']);m['jobs'][0].update(submission_status='submitting',job_id='1234')
    atomic_json(s['manifest'],m)
    def unknown(*a,**kw):raise RuntimeError('Interrupted Slurm submission requires scheduler reconciliation')
    monkeypatch.setattr(panel.recon,'run_many',unknown)
    monkeypatch.setattr(panel.dispatch,'advance',lambda *a,**kw:pytest.fail('Uncertain submission must not trigger follow-ons.'))
    with pytest.raises(RuntimeError,match='scheduler reconciliation'):panel.run(out)
    saved=panel.read(s['manifest'])['jobs'][0]
    assert saved['job_id']=='1234' and saved['submission_status']=='submitting'
    assert not (out/'relaxations.json').exists()


def test_failed_native_prerequisite_is_recorded_without_erasing_other_work(synthetic,monkeypatch):
    _,evidence,pool,out=synthetic
    panel.stage(out,pool=pool)
    for s in panel.read(out/'panel.json')['systems']:
        m=panel.read(s['manifest'])
        for j in m['jobs']:j['status']='failed'
        atomic_json(s['manifest'],m)
    def failures(*a,**k):raise panel.recon.ProbeFailure('Synthetic native failure')
    monkeypatch.setattr(panel.recon,'run_many',failures)
    monkeypatch.setattr(panel.dispatch,'advance',lambda *a,**kw:pytest.fail('Failed prerequisite cannot be optimized.'))
    panel.run(out)
    r=panel.read(out/'completion.json')
    assert r['status']=='incomplete_review_required' and len(r['holds'])==4
    assert not r['production_model_accepted']


def test_capacity_only_change_adopts_work_and_preserves_exclusions(synthetic):
    plan,_,pool,out=synthetic
    plan['protocol']['max_jobs']=3
    atomic_json(panel.HERE/'plan.json',plan)
    run_policy.update(pool,max_jobs=3)
    panel.stage(out,pool=pool)
    system=panel.read(out/'panel.json')['systems'][0]
    saved=panel.read(system['manifest'])
    saved['jobs'][0].update(status='submitted',job_id='98765')
    atomic_json(system['manifest'],saved)
    before={p:p.read_bytes() for p in out.rglob('*.inp')}
    run_policy.update(pool,max_jobs=23)
    plan['protocol']['max_jobs']=23
    atomic_json(panel.HERE/'plan.json',plan)
    panel.stage(out,pool=pool)
    assert run_policy.read(pool)['excluded_nodes']==['node9']
    after=panel.read(system['manifest'])['jobs'][0]
    assert after['status']=='submitted' and after['job_id']=='98765'
    assert list((out/'plan_history').glob('*.json'))
    assert all(p.read_bytes()==value for p,value in before.items())
    with panel.recon._manifest_lock(out):
        with pytest.raises(RuntimeError,match='Another reconnaissance controller'):
            panel.stage(out,pool=pool)
    plan['protocol']['shift']=.2
    atomic_json(panel.HERE/'plan.json',plan)
    with pytest.raises(ValueError,match='scientific panel differs'):
        panel.stage(out,pool=pool)


def test_pool_above_authorized_cap_refused(synthetic):
    _,_,pool,out=synthetic
    run_policy.update(pool,max_jobs=24)
    with pytest.raises(ValueError,match='capped at 23'):
        panel.stage(out,pool=pool)
