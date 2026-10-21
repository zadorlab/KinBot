from copy import deepcopy
from dataclasses import asdict
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from kinbot.execution import atomic_json
from kinbot.anl import network_update, ts_reuse


def network():
    return dict(schema=1, resources={'cores':8}, network=dict(initial='w',
        route_to_smiles={'w':'C'}, stable_species=[dict(routing_name='w', roles=['initial_well'],smiles='C')],
        transition_states={}, reactions=[]), cbh=dict(targets={}, terminal_atct_species={},
        calculated_states=[], thermochemistry_states=[]))


def test_append_only_network_expansion_and_changes_refused():
    old=network(); new=deepcopy(old)
    new['network']['reactions'].append(dict(name='new',reactant='w',products=['new']))
    new['network']['route_to_smiles']['new']='[CH3]'
    result=network_update.extend(old,new)
    assert result==new and old!=new
    altered=deepcopy(new);altered['resources']['cores']=16
    with pytest.raises(ValueError,match='resources'):
        network_update.extend(old,altered)
    altered=deepcopy(new);altered['network']['stable_species'][0]['smiles']='N'
    with pytest.raises(ValueError,match='saved scientific'):
        network_update.extend(old,altered)
    altered=deepcopy(new);altered['network']['stable_species']=[]
    with pytest.raises(ValueError,match='saved scientific'):
        network_update.extend(old,altered)


def test_reopening_removes_derived_energies_without_changing_l2(monkeypatch):
    from kinbot.anl import handoff
    point=SimpleNamespace(energy=-100.,zpe=.01,geom=[[0,0,0]],freq=[100.],
        final_zero_k_energy='composite',barrier_zero_k_energy='barrier',
        formation_enthalpy_0k='formation',anl_thermochemistry_frequencies=[99.])
    monkeypatch.setattr(handoff,'_network_species',lambda root:({'w':point},{}))
    handoff.detach_handoff(object())
    assert vars(point)==dict(energy=-100.,zpe=.01,geom=[[0,0,0]],freq=[100.])


def test_reference_promoted_to_network_gets_new_geometry_graph():
    old=network(); new=deepcopy(old)
    ref=dict(smiles='O',charge=0,multiplicity=1,source='automatic_reference_l2',network_record=None)
    old['cbh']['calculated_states']=[ref]
    new['cbh']['calculated_states']=[dict(ref,source='network_l2',network_record={'positions':'new'})]
    result=network_update.extend(old,new)
    record=result['cbh']['calculated_states'][0]
    assert record['graph_revision']==ts_reuse.digest({'positions':'new'})[:16]
    assert network_update.extend(result,new)==result


def test_new_role_retains_geometry_graph_revision():
    old=network()
    ref=dict(smiles='C',charge=0,multiplicity=1,source='network_l2',
             network_record=deepcopy(old['network']['stable_species'][0]),graph_revision='retained')
    old['cbh']['calculated_states']=[ref]
    new=deepcopy(old)
    new['network']['stable_species'][0]['roles'].append('product')
    new['cbh']['calculated_states'][0]['network_record']['roles'].append('product')
    new['cbh']['calculated_states'][0].pop('graph_revision')
    assert network_update.extend(old,new)['cbh']['calculated_states'][0]['graph_revision']=='retained'


def test_live_network_revisit_replays_cache_and_grants_one_failed_attempt(tmp_path,monkeypatch):
    import os
    from kinbot import pes_phase,live_control,search_control
    live_control.configure(tmp_path)
    pes_phase.initialize(tmp_path)
    worker=tmp_path/'123';worker.mkdir()
    inp=worker/'123.json';atomic_json(inp,dict(vdW_detection=.5))
    pes_phase.register(tmp_path,'123',os.getpid())
    class Reaction:
        def __init__(self,species,qc,par,instance,name):
            self.instance,self.instance_name=instance,name
    species=SimpleNamespace(reac_ts_done=[-1,-999],reac_step=[0,0])
    par={'vdW_detection':.5}
    species.reac_obj=[Reaction(species,None,par,[],n) for n in ('accepted','rejected')]
    generator=SimpleNamespace(species=species,par=par,qc=None,inp=str(inp))
    r=live_control.request(tmp_path,worker,kind='network',reason='reassess',options={'vdW_detection':0.,'rejected':True})
    monkeypatch.chdir(worker)
    search_control.revise_network(generator)
    assert species.reac_ts_done==[0,0] and par['vdW_detection']==0
    m=json.loads(Path('.kinbot_revisit.json').read_text())
    assert m['rejected']==['rejected'] and m['attempted']==[]
    archived=[]
    monkeypatch.setattr(search_control,'_retire_reaction',lambda *a:archived.append(a[1]))
    species.reac_ts_done=[-1,-999]
    assert search_control.pending_revisit(generator)
    search_control.reactions(generator)
    assert archived==['rejected'] and species.reac_ts_done==[-1,0]
    species.reac_ts_done=[-1,-999]
    search_control.reactions(generator)
    assert archived==['rejected'] and not search_control.pending_revisit(generator)


def test_network_publication_requires_exact_authorization_and_archives_only_assembly(tmp_path):
    old=network(); new=deepcopy(old);new['network']['reactions']=[dict(name='new')]
    atomic_json(tmp_path/'network_request.json',old)
    assert not network_update.publish(old,new,tmp_path)
    atomic_json(tmp_path/'network_revisit.json',dict(previous_sha256='bad',archive=str(tmp_path/'archive')))
    with pytest.raises(ValueError,match='authorization'):
        network_update.publish(old,new,tmp_path)
    atomic_json(tmp_path/'network_revisit.json',dict(previous_sha256=ts_reuse.digest(old),archive=str(tmp_path/'archive')))
    atomic_json(tmp_path/'execution.json',{'status':'complete'})
    (tmp_path/'native.out').write_text('retained')
    assert network_update.publish(old,new,tmp_path)==new
    assert (tmp_path/'native.out').read_text()=='retained'
    assert (tmp_path/'archive/execution.json').exists()
    assert json.loads((tmp_path/'network_request.json').read_text())==new


def saddle_plan():
    m=network()
    r=dict(smiles='C[S]', l2_theory='theory',highest_tier='ANL1',fallback_tiers=['L3'],
        molecule=dict(symbols=['C','S','H','H'],charge=0,multiplicity=2,
                      positions=[[0,0,0],[1,0,0],[0,1,0],[0,0,1]]),
        energy_hartree=-20.,harmonic_zpe_hartree=.01,frequencies_cm_inverse=[-200,300,400])
    m['network']['transition_states']={'a':r,'b':deepcopy(r),'different':deepcopy(r)}
    m['network']['transition_states']['different']['molecule']['positions'][1][0]=1.1
    m['network']['reactions']=[dict(transition_state=n,reactant='w',products=['p']) for n in ('a','b','different')]
    return m


def test_ts_identity_not_shared_product_only_and_no_reflections(tmp_path):
    m=saddle_plan(); limits=ts_reuse.options({})
    assert ts_reuse.equivalent(m,'a','b',limits)
    assert ts_reuse.equivalent(m,'a','different',limits) is None
    m['network']['transition_states']['b']['molecule']['positions'][3][2]=-1
    assert ts_reuse.equivalent(m,'a','b',limits) is None


def test_ts_reuse_prefers_completed_then_active_and_freezes_canonical(tmp_path):
    m=saddle_plan()
    atomic_json(ts_reuse.directory(tmp_path,'b')/'state.json',dict(tasks={'geometry':{'status':'complete'}}))
    mapping=ts_reuse.discover(m,tmp_path,ts_reuse.options({}))
    assert mapping['aliases']=={'a':'b'}
    atomic_json(tmp_path/'ts_reuse.json',mapping)
    assert ts_reuse.load(m,tmp_path)=={'a':'b'}
    m['network']['transition_states']['a']['energy_hartree']+=.001
    with pytest.raises(ValueError,match='provenance'):
        ts_reuse.load(m,tmp_path)


def test_new_alias_refuses_to_strand_existing_graph(tmp_path):
    m=saddle_plan()
    for name in ('a','b'):
        atomic_json(ts_reuse.directory(tmp_path,name)/'state.json',dict(tasks={}))
    with pytest.raises(ValueError,match='already exists'):
        ts_reuse.install_new(m,tmp_path,{})
    assert not (tmp_path/'ts_reuse.json').exists()


def test_binding_formation_has_correct_sign_state_and_recursive_provenance(tmp_path):
    from test_anl_cbh_atct import _table
    from kinbot.anl.cbh import CBHReaction,ZeroKEnergy,solve_formation_enthalpy,formation_from_dict,HARTREE_TO_KJ_MOL
    from kinbot.anl.automatic import _terminal_formation
    _, table=_table()
    methane=ZeroKEnergy('C',-40.,'ANL0','methane.json')
    terminal=_terminal_formation('C',(0,1),{'C':1,'H':4},table,methane)
    reaction=CBHReaction(0,'C.C',{'C.C':-1,'C':2},{'C.C':{'C':2,'H':8},'C':{'C':1,'H':4}})
    energies={'C.C':ZeroKEnergy('C.C',-80.01,'ANL0','complex.json'),'C':methane}
    f=solve_formation_enthalpy(reaction,energies,binding_energy=True,reference_formations={'C':terminal})
    assert f.formation_0k_kj_mol==pytest.approx(2*terminal.formation_0k_kj_mol-.01*HARTREE_TO_KJ_MOL)
    assert f.construction=='binding_energy'
    assert formation_from_dict(asdict(f))==f
    bad=asdict(f);bad['construction']='cbh'
    with pytest.raises(ValueError,match='lower-rung'):
        formation_from_dict(bad)
    energies['C']=ZeroKEnergy('C',-40.,'other method','methane.json')
    with pytest.raises(ValueError,match='same'):
        solve_formation_enthalpy(reaction,energies,binding_energy=True,reference_formations={'C':terminal})


def test_complex_is_added_to_plan_and_gets_its_own_graph(tmp_path, monkeypatch):
    from test_anl_automatic import _species
    from kinbot.anl import automatic as a
    root=_species('root','CO',['C','O','H','H','H','H'])
    methyl=_species('methyl','[CH3]',['C','H','H','H'],2)
    oh=_species('oh','[OH]',['O','H'],2)
    complex_=_species('complex','[CH3].[OH]',root.atom)
    ts=_species('ts','CO',root.atom)
    root.reac_obj=[SimpleNamespace(instance_name='ts',family_name='addition',products=[methyl,oh],
                                   ts=ts,do_vdW=True,irc_prod_opt=SimpleNamespace(species=complex_))]
    root.reac_ts_done=[-1];root.reac_type=['addition']
    plan=a.build_network_request(root,{'composite_method':'ANL0-F12','queue_job_limit':2})
    assert plan['network']['product_complexes']['complex']['products']==['methyl','oh']
    # Planning must remain noncovalent; it is not a fake CBH covalent target.
    assert 'complex' not in plan['cbh']['targets']
    made=[]
    monkeypatch.setattr(a,'prepare',lambda spec,graph:made.append((spec,graph)))
    monkeypatch.setattr(a,'preflight_until_scheduler_ready',lambda *a,**k:{})
    sources={a._state_key(e['smiles'],e['charge'],e['multiplicity']):
             (e['network_record'] or plan['network']['stable_species'][0]) for e in plan['cbh']['calculated_states']}
    graphs=a._prepare_graphs(plan,sources,tmp_path/'states')
    assert 'complex:complex' in graphs
    spec=next(spec for spec,graph in made if graph==graphs['complex:complex'])
    assert set(spec['molecule']['smiles'].split('.'))=={'[CH3]','[OH]'}


def test_refresh_network_edit_is_idempotent_and_keeps_qc(tmp_path, monkeypatch):
    from kinbot import network_revisit as nr, run_salvage, run_refresh
    worker=tmp_path/'case'/'123';worker.mkdir(parents=True)
    case=worker.parent/'input.json'
    original=dict(vdW_detection=.5,profiled_l1_reuse_only=1)
    atomic_json(case,original);atomic_json(worker/'123.json',original)
    (worker/'kinbot_monitor.out').write_text('-1 0 accepted\n-999 0 rejected\n')
    (worker/'kinbot.db').write_text('untouched')
    prepared={'sha256':{'case/input.json':__import__('kinbot.execution',fromlist=['file_sha256']).file_sha256(case)}}
    atomic_json(tmp_path/'prepared.json',prepared)
    monkeypatch.setattr(run_salvage,'_cases',lambda root:(prepared,[case]))
    monkeypatch.setattr(run_refresh,'_owned_jobs',lambda root:[])
    plan=nr.plan(tmp_path,[],nr.settings(0,True,False))
    assert plan['workers'][0]['rejected']==['rejected']
    journal={'network_revision':plan}
    marker=tmp_path/'journal.json';atomic_json(marker,journal)
    archive=tmp_path/'archive'
    nr.apply(tmp_path,plan,archive,journal,marker)
    assert json.loads(case.read_text())['vdW_detection']==0
    assert json.loads(case.read_text())['pes_phase_policy']=='exploration_first'
    pending=worker/'.kinbot_revisit.json';m=json.loads(pending.read_text());m['attempted']=['rejected'];atomic_json(pending,m)
    nr.apply(tmp_path,plan,archive,journal,marker)
    assert json.loads(pending.read_text())['attempted']==['rejected']
    assert (worker/'kinbot.db').read_text()=='untouched'
