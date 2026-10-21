from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest
from kinbot.execution import atomic_json
from kinbot.anl import ts_reuse, pes_ts_reuse as reuse
from test_network_revision import saddle_plan


def test_permutations_and_reverse_endpoints_require_same_saddle():
    plan=saddle_plan()
    plan['network']['reactions'][1].update(reactant='p',products=['w'])
    record=plan['network']['transition_states']['b']
    record['smiles']='[CH2]S'  # reactant connectivity anchor changes in reverse
    mol=record['molecule']; permutation=[1,0,3,2]
    mol['symbols']=[mol['symbols'][i] for i in permutation]
    mol['positions']=[mol['positions'][i] for i in permutation]
    evidence=ts_reuse.equivalent(plan,'a','b',ts_reuse.options({}))
    assert evidence and evidence['rmsd_angstrom']<1e-10
    mol['positions'][0][0]+=.02
    assert ts_reuse.equivalent(plan,'a','b',ts_reuse.options({})) is None


def test_reflected_distinct_ligands_are_not_merged():
    plan=saddle_plan()
    for entry in plan['network']['transition_states'].values():
        entry['molecule']['symbols']=['C','S','H','F']
    plan['network']['transition_states']['b']['molecule']['positions'][3][2]=-1
    assert ts_reuse.equivalent(plan,'a','b',ts_reuse.options({})) is None


def fixture(case):
    atomic_json(case/'pes_exploration.json',{})
    bases=[]
    for worker,name,reactant,product in [('111','forward','w','p'),('222','reverse','p','w')]:
        base=case/worker/'anl_automatic'; bases.append(base)
        atomic_json(base.parent/(worker+'.json'),{'anl_ts_reuse':{}})
        plan=saddle_plan(); record=deepcopy(plan['network']['transition_states']['a'])
        plan['network']['transition_states']={name:record}
        plan['network']['reactions']=[dict(transition_state=name,reactant=reactant,products=[product])]
        plan.update(l2_theory='same',requested_method='ANL1')
        plan['resources']['global_max_jobs']=2
        atomic_json(base/'network_request.json',plan)
    return bases


def test_pes_registry_prefers_completed_owner_and_never_crosses_cases(tmp_path):
    a,b=fixture(tmp_path/'case_a'); other=fixture(tmp_path/'case_b')
    atomic_json(ts_reuse.directory(b,'reverse')/'state.json',{'tasks':{'sp':{'status':'complete'}}})
    reuse.install(a,{})
    target=reuse.load(a)['forward']
    assert target['graph']==ts_reuse.directory(b,'reverse')
    assert reuse.load(b)=={} and reuse.load(other[0])=={}
    reuse.install(b,{})
    assert reuse.load(a)['forward']['name']=='reverse'
    plan=json.loads((b/'network_request.json').read_text())
    plan['network']['transition_states']['reverse']['energy_hartree']+=.01
    atomic_json(b/'network_request.json',plan)
    with pytest.raises(ValueError,match='provenance'):
        reuse.load(a)


def test_existing_reverse_graph_needs_verified_retirement(tmp_path):
    a,b=fixture(tmp_path/'case')
    for base,name in [(a,'forward'),(b,'reverse')]:
        atomic_json(ts_reuse.directory(base,name)/'state.json',{'tasks':{}})
    with pytest.raises(ValueError,match='already exists'):
        reuse.install(a,{})
    assert not (a.parent.parent/reuse.NAME).exists()


def test_changed_theory_prevents_cross_well_sharing(tmp_path):
    a,b=fixture(tmp_path/'case')
    plan=json.loads((b/'network_request.json').read_text());plan['resources']['cc_convergence']=1e-6
    atomic_json(b/'network_request.json',plan)
    reuse.install(a,{})
    assert reuse.load(a)=={} and reuse.load(b)=={}


def test_shared_graph_waits_for_live_owner_private_work_then_adopts(tmp_path,monkeypatch):
    from kinbot import pes_phase
    from kinbot.anl import automatic
    a,b=fixture(tmp_path/'case');reuse.install(a,{})
    target=reuse.load(b)['reverse']; graph=target['graph']
    calls=[]
    monkeypatch.setattr(pes_phase,'checkpoint',lambda:None)
    monkeypatch.setattr(pes_phase,'active_worker',lambda w:True)
    def sleep(interval):
        calls.append('wait')
        atomic_json(graph.parents[2]/'ts_private_graphs_complete.json',
                    {'plan_sha256':automatic._sha256_json(target['plan'])})
    monkeypatch.setattr(reuse.time,'sleep',sleep)
    def prepare(*args):
        graph.mkdir(parents=True,exist_ok=True);calls.append('prepare')
    monkeypatch.setattr(automatic,'_prepare_ts_graph',prepare)
    monkeypatch.setattr(automatic,'preflight_until_scheduler_ready',lambda *a,**k:calls.append('verify'))
    monkeypatch.setattr(automatic,'_drive_graphs_locked',lambda *a:calls.append('drive'))
    assert reuse.finish(target)==graph
    assert calls==['wait','prepare','verify','drive']


def add_minimum(base, *, complex_route=None, revision=''):
    plan=json.loads((base/'network_request.json').read_text())
    record=deepcopy(next(iter(plan['network']['transition_states'].values())))
    record['frequencies_cm_inverse']=[100.,300.,400.]
    if complex_route:
        record['smiles']='[CH3].[SH]'
        plan['network'].setdefault('product_complexes',{})[complex_route]=dict(record=record,products=['CH3','SH'])
        name='complex:'+complex_route
    else:
        entry=dict(smiles=record['smiles'],charge=0,multiplicity=2,highest_tier='ANL1',fallback_tiers=['L3'],network_record=record)
        if revision: entry['graph_revision']=revision
        plan['cbh']['calculated_states']=[entry]
        name='stable:C[S]|0|2'
    atomic_json(base/'network_request.json',plan)
    return name


def test_stable_and_complex_reuse_but_not_different_conformer(tmp_path):
    a,b=fixture(tmp_path/'case')
    stable=add_minimum(a);add_minimum(b)
    add_minimum(a,complex_route='vdw_a');add_minimum(b,complex_route='vdw_b')
    reuse.install(a,{})
    shared=reuse.load(b)
    assert stable in shared and 'complex:vdw_b' in shared
    assert shared[stable]['kind']=='stable'
    assert shared['complex:vdw_b']['name']=='complex:vdw_a'
    # Distinct complex geometry is not a duplicate, even with identical fragments.
    x,y=shared['complex:vdw_b'],deepcopy(shared['complex:vdw_b'])
    y['record']['molecule']['positions'][1][0]+=.1
    assert reuse.compare(x,y,ts_reuse.options({})) is None


def test_immutable_unique_owner_prevents_simultaneous_registration_race(tmp_path):
    a,b=fixture(tmp_path/'case')
    saved=(b/'network_request.json').read_bytes();(b/'network_request.json').unlink()
    reuse.install(a,{})
    # A new peer with a lexically earlier task cannot demote a reserved owner,
    # even if the first driver has not yet created its graph.
    (b/'network_request.json').write_bytes(saved)
    reuse.install(b,{})
    assert reuse.load(a)=={}
    assert reuse.load(b)['reverse']['graph']==ts_reuse.directory(a,'forward')


def test_old_reference_graph_survives_promotion_and_role_growth(tmp_path):
    a,b=fixture(tmp_path/'case');name=add_minimum(a);add_minimum(b)
    reuse.install(a,{})
    old=reuse.load(b)[name]['graph']
    add_minimum(a,revision='promoted')
    plan=json.loads((a/'network_request.json').read_text())
    plan['cbh']['calculated_states'][0]['network_record']['roles']=['product','well']
    atomic_json(a/'network_request.json',plan)
    reuse.install(a,{})
    assert reuse.load(b)[name]['graph']==old
    assert reuse.load(a)[name]['graph']==old
    assert old.name=='composite'


def test_different_vdw_names_on_reverse_ts_are_matched_only_with_complex_proof(tmp_path):
    a,b=fixture(tmp_path/'case')
    add_minimum(a,complex_route='vdw_a');add_minimum(b,complex_route='vdw_b')
    for base,start,end in [(a,'w','vdw_a'),(b,'vdw_b','w')]:
        p=json.loads((base/'network_request.json').read_text())
        p['network']['reactions'][0].update(reactant=start,products=[end])
        atomic_json(base/'network_request.json',p)
    reuse.install(a,{})
    assert reuse.load(b)['reverse']['name']=='forward'


def test_shared_stable_graph_reaches_vpt2_and_final_handoff_bookkeeping(tmp_path,monkeypatch):
    from kinbot.anl import automatic
    a,b=fixture(tmp_path/'case');stable=add_minimum(a);add_minimum(b)
    reuse.install(a,{})
    monkeypatch.setattr(reuse,'finish',lambda target,interval:target['graph'])
    calls=[]
    def assemble(graph,record,destination,**options):
        calls.append((graph,options))
        return dict(thermochemistry={'ANL1':'profile.json'},barrier={})
    monkeypatch.setattr(automatic,'assemble_available_profiles',assemble)
    graphs={};external,profiles=automatic._shared_profiles(b,graphs,{'electronic_state_models':{}},.01)
    assert graphs[stable]==external[stable]['graph']
    assert profiles[stable]['thermochemistry']['ANL1']=='profile.json'
    stable_options=next(kw for graph,kw in calls if kw['stationary_point']=='minimum')
    assert stable_options['state_id']=='C[S]|0|2'
    assert json.loads((b/'pes_reuse_consumed.json').read_text())['graphs'][stable]['canonical_name']==stable


def test_migration_cancels_only_proven_duplicate_and_archives_without_touching_owner(tmp_path,monkeypatch):
    from kinbot import network_revisit as nr,run_salvage,run_refresh,vrc_restart
    from kinbot.anl import dispatch
    from tests.anl_fixture import dispatch_spec,molpro_task
    from tests.test_anl_dispatch import _mock_execution
    from kinbot.execution import file_sha256
    a,b=fixture(tmp_path/'case');case=a.parent.parent
    atomic_json(case/'input.json',{})
    prepared=dict(inputs=[str(case/'input.json')],sha256={'case/input.json':file_sha256(case/'input.json')})
    atomic_json(tmp_path/'prepared.json',prepared)
    monkeypatch.setattr(run_salvage,'_cases',lambda root:(prepared,[case/'input.json']))
    for base,name in [(a,'forward'),(b,'reverse')]:
        record=json.loads((base/'network_request.json').read_text())['network']['transition_states'][name]
        record['molecule']['multiplicity']=1
        plan=json.loads((base/'network_request.json').read_text());plan['network']['transition_states'][name]=record
        atomic_json(base/'network_request.json',plan)
        spec=dispatch_spec();spec['molecule']=record['molecule'];spec['tasks']=[molpro_task('sp','rhf\n',geometry_from='initial')]
        dispatch.prepare(spec,ts_reuse.directory(base,name))
    canonical=ts_reuse.directory(b,'reverse');duplicate=ts_reuse.directory(a,'forward')
    _mock_execution(canonical,'sp');dispatch.advance(canonical)
    state=dispatch._load(duplicate)[2];state['tasks']['sp'].update(status='submitted',job_id='41');atomic_json(duplicate/'state.json',state)
    old={p:file_sha256(p) for p in canonical.rglob('*') if p.is_file()}
    job=dict(job_id='41',directory=str(duplicate/'tasks/sp'),state='RUNNING')
    plan=nr.plan(tmp_path,[job],nr.settings(deduplicate=True))
    assert plan['selected_jobs']==[job]
    cancelled=[]
    monkeypatch.setattr(vrc_restart,'cancel_scans',lambda record,*a:cancelled.extend(record['selected_jobs']))
    monkeypatch.setattr(run_refresh,'_owned_jobs',lambda root:[])
    journal={'network_revision':plan};marker=tmp_path/'journal.json';atomic_json(marker,journal)
    nr.apply(tmp_path,plan,tmp_path/'archive',journal,marker)
    assert cancelled==[job] and not duplicate.exists()
    assert (tmp_path/'archive/duplicate_ts'/duplicate.relative_to(tmp_path)/'workflow.json').exists()
    assert old=={p:file_sha256(p) for p in old}
    assert reuse.load(a)['forward']['graph']==canonical
