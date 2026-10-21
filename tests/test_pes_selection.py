"""The outer PES selects L2 routes before ANL; accepted QC is never erased."""
from copy import deepcopy
import itertools
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import pes_selection as selection, pes_phase as phase
from kinbot.execution import atomic_json, file_sha256


def evidence(worker, rows):
    stable = {worker, *(p for _, products, _, _ in rows for p in products)}
    def record(energy):
        return dict(energy_hartree=energy, harmonic_zpe_hartree=.01,
                    l2_theory={'method':'B2PLYP','basis':'cc-pVTZ'})
    return dict(errors=[], network=dict(initial=worker,
        stable_species=[dict(record(-100.), routing_name=n) for n in sorted(stable)],
        transition_states={name: record(energy) for name, _, energy, _ in rows
                           if 'hom_sci' not in name},
        reactions=[dict(name=name, reactant=worker, products=products,
                        path_class=path, barrierless='hom_sci' in name,
                        transition_state=None if 'hom_sci' in name else name)
                   for name, products, _, path in rows]))


def test_global_selection_matches_final_pes_including_reverse_and_stereo():
    a = evidence('a', [('forward_high',['b'],-99.,'ordinary'),
                       ('distinct_stereo',['b'],-98.,'stereopath:x'),
                       ('hom_sci_a',['b'],0.,None)])
    b = evidence('b', [('reverse_low',['a'],-99.5,'ordinary')])
    for order in itertools.permutations([('a',a),('b',b)]):
        result = selection.choose(dict(order))
        assert result['workers']['a']['selected'] == ['distinct_stereo']
        assert result['workers']['b']['selected'] == ['reverse_low']
        assert result['workers']['a']['omitted']['forward_high']['retained'] == ['reverse_low']
        # No claim that geometrically different high/low saddles are equal.
        assert result['workers']['a']['omitted']['forward_high']['reason'] == 'outer_pes_lowest_l2_pathway'


def test_same_products_from_different_wells_are_not_same_pathway():
    result = selection.choose({n:evidence(n, [(n+'_ts',['p','h'],-99.,'ordinary')])
                               for n in ('a','b')})
    assert result['workers']['a']['selected'] == ['a_ts']
    assert result['workers']['b']['selected'] == ['b_ts']


@pytest.mark.parametrize('damage', ['theory','energy','class','error'])
def test_incomplete_or_incomparable_evidence_fails_before_dispatch(damage):
    a=evidence('a',[('ts',['b'],-99.,'ordinary')])
    b=evidence('b',[('reverse',['a'],-99.,'ordinary')])
    if damage=='theory': b['network']['transition_states']['reverse']['l2_theory']['basis']='other'
    if damage=='energy': b['network']['transition_states']['reverse']['energy_hartree']=float('nan')
    if damage=='class': b['network']['reactions'][0].pop('path_class')
    if damage=='error': b['errors']=[dict(error='missing accepted energy')]
    with pytest.raises(ValueError): selection.choose({'a':a,'b':b})


def test_outer_owner_does_not_select_until_all_wells_ready(tmp_path):
    m=phase.initialize(tmp_path)
    for name in ('a','b'):
        atomic_json(tmp_path/name/(name+'.json'),{'composite_method':'ANL1'})
        row=evidence(name,[(name+'_ts',['b' if name=='a' else 'a'],-99. if name=='a' else -99.5,'ordinary')])
        row['generation']=dict(token=m['token'],epoch=0)
        atomic_json(tmp_path/name/'accepted_l2_network.json',row)
        phase.register(tmp_path,name,os.getpid())
    def ready(name):
        data=json.loads((tmp_path/phase.NAME).read_text())
        data['workers'][name]['status']='ready';atomic_json(tmp_path/phase.NAME,data)
    ready('a')
    assert phase.tick(tmp_path,['a','b'],max_searches=1)[1]=='exploring'
    assert not (tmp_path/selection.NAME).exists()
    ready('b')
    assert phase.tick(tmp_path,['a','b'],max_searches=1)[1]=='postprocessing'
    result=json.loads((tmp_path/selection.NAME).read_text())
    assert result['workers']['a']['selected']==[]
    assert result['workers']['b']['selected']==['b_ts']


def test_incomplete_selection_does_not_block_vrc_phase(tmp_path):
    m=phase.initialize(tmp_path)
    atomic_json(tmp_path/'a/a.json',{'composite_method':'ANL1'})
    phase.register(tmp_path,'a',os.getpid())
    m=json.loads((tmp_path/phase.NAME).read_text());m['workers']['a']['status']='ready'
    atomic_json(tmp_path/phase.NAME,m)
    assert phase.tick(tmp_path,['a'],max_searches=1)[1]=='postprocessing'
    assert json.loads((tmp_path/selection.NAME).read_text())['status']=='blocked'


@pytest.mark.parametrize('missing_input', [True, False])
def test_selection_does_not_silently_ignore_a_worker(tmp_path, missing_input):
    atomic_json(tmp_path/'a/a.json', {'composite_method':'ANL1'})
    if not missing_input:
        atomic_json(tmp_path/'b/b.json', {})
    selection.publish(tmp_path, dict(token='x', epoch=0), ['a','b'])
    result=json.loads((tmp_path/selection.NAME).read_text())
    assert result['status']=='blocked' and 'worker inventory' in result['error']


def test_selection_is_nonmutating_and_old_generations_are_rejected(tmp_path,monkeypatch):
    a=evidence('a',[('kept',['p'],-100.,'ordinary'),('omitted',['p'],-99.,'ordinary')])
    result=selection.choose({'a':a});result.update(token='generation',epoch=3)
    atomic_json(tmp_path/selection.NAME,result)
    monkeypatch.setattr(phase,'_ROOT',tmp_path);monkeypatch.setattr(phase,'_WORKER','a')
    monkeypatch.setattr(phase,'_read',lambda:dict(token='generation',epoch=3))
    monkeypatch.setattr('kinbot.pes_snapshot.snapshot',lambda *args:deepcopy(a))
    species=SimpleNamespace(reac_obj=[SimpleNamespace(instance_name=n) for n in ('kept','omitted')],
                            reac_ts_done=[-1,-1])
    view=selection.selected_view(species,{'composite_method':'ANL1'})
    assert view.reac_ts_done==[-1,-2] and species.reac_ts_done==[-1,-1]
    assert view.anl_accepted_source is species and view.reac_obj is species.reac_obj
    result['epoch']=2;atomic_json(tmp_path/selection.NAME,result)
    with pytest.raises(ValueError,match='another PES exploration generation'):
        selection.selected_view(species,{'composite_method':'ANL1'})


def test_anl_preflight_skips_unselected_saved_graph_and_reuses_selected_owner(tmp_path,monkeypatch):
    from kinbot.anl import automatic, ts_reuse
    from test_network_revision import saddle_plan
    plan=saddle_plan()
    plan['cbh']['calculated_states']=[]
    mapping=ts_reuse.discover(plan,tmp_path,ts_reuse.options({}))
    assert mapping['aliases']=={'b':'a'}
    atomic_json(tmp_path/'ts_reuse.json',mapping)
    for name in plan['network']['transition_states']:
        graph=ts_reuse.directory(tmp_path,name)
        atomic_json(graph/'workflow.json',{'native':name})
        atomic_json(graph/'state.json',{'tasks':{'sp':{'status':'submitted','job_id':name}}})
    files={p:file_sha256(p) for p in (tmp_path/'states').rglob('*.json')}
    checked=[]
    monkeypatch.setattr(automatic,'preflight_until_scheduler_ready',lambda graph,**kw:checked.append(graph) or {})
    graphs=automatic._prepare_graphs(plan,{},tmp_path/'states',selection={'selected':['b']})
    assert graphs=={'ts:a':ts_reuse.directory(tmp_path,'a')}
    assert checked==list(graphs.values())
    assert {p:file_sha256(p) for p in files}==files
    # Repeating preparation adopts precisely the same owner, with no new graph.
    assert automatic._prepare_graphs(plan,{},tmp_path/'states',selection={'selected':['b']})==graphs


def test_reselection_archives_assembly_without_touching_qc(tmp_path):
    plan={'network':{'transition_states':{'a':{},'b':{}}}}
    atomic_json(tmp_path/'states/saved/state.json',{'tasks':{'energy':{'job_id':'123','status':'submitted'}}})
    before=(tmp_path/'states/saved/state.json').read_bytes()
    atomic_json(tmp_path/'handoff.json',{'old':'all pathways'})
    selection.bind(tmp_path,plan,{'selected':['a']})
    assert not (tmp_path/'handoff.json').exists()
    assert list((tmp_path/'selection_history').glob('*/handoff.json'))
    atomic_json(tmp_path/'handoff.json',{'new':'selected a'})
    selection.bind(tmp_path,plan,{'selected':['a']})
    assert (tmp_path/'handoff.json').exists()
    selection.bind(tmp_path,plan,{'selected':['a','b']})
    assert not (tmp_path/'handoff.json').exists()
    assert (tmp_path/'states/saved/state.json').read_bytes()==before


@pytest.mark.parametrize('crash_at', ['marker', 'journal'])
def test_selection_replays_after_archive_or_marker_publication(tmp_path, monkeypatch, crash_at):
    plan={'network':{'transition_states':{'a':{},'b':{}}}}
    selection.bind(tmp_path, plan, {'selected':['a']})
    atomic_json(tmp_path/'handoff.json', {'selected':['a']})
    atomic_json(tmp_path/'states/native/state.json', {'job_id':'live'})
    original=(tmp_path/'states/native/state.json').read_bytes()
    def interrupted(path, data):
        if ((crash_at=='marker' and Path(path).name=='kinetic_selection.json') or
            (crash_at=='journal' and Path(path).name=='selection_update.json' and data['status']=='complete')):
            raise OSError('interrupted selection')
        atomic_json(path,data)
    with monkeypatch.context() as m:
        m.setattr(selection,'atomic_json',interrupted)
        with pytest.raises(OSError,match='interrupted selection'):
            selection.bind(tmp_path,plan,{'selected':['b']})
    assert not (tmp_path/'handoff.json').exists()
    selection.bind(tmp_path,plan,{'selected':['b']})
    assert json.loads((tmp_path/'selection_update.json').read_text())['status']=='complete'
    atomic_json(tmp_path/'handoff.json',{'selected':['b']})
    selection.bind(tmp_path,plan,{'selected':['a']})
    assert len(list((tmp_path/'selection_history').glob('*/handoff.json')))==2
    assert (tmp_path/'states/native/state.json').read_bytes()==original


def test_unmanaged_assembly_has_no_new_selection_state(tmp_path):
    atomic_json(tmp_path/'handoff.json', {'legacy':'complete'})
    selection.bind(tmp_path, {'network':{'transition_states':{}}}, None)
    assert (tmp_path/'handoff.json').exists()
    assert not (tmp_path/'kinetic_selection.json').exists()


def test_selected_cbh_projection_has_only_required_targets_and_keeps_qc_inputs():
    from kinbot.anl.automatic import build_network_request
    from test_anl_automatic import _species
    root=_species('well','CC=O',['C','C','O','H','H','H','H'],energy=-153.)
    methyl=_species('methyl','[CH3]',['C','H','H','H'],2)
    formyl=_species('formyl','[CH]=O',['C','O','H'],2)
    ketene=_species('ketene','C=C=O',['C','C','O','H','H'])
    hydrogen=_species('hydrogen','[H][H]',['H','H'])
    root.reac_obj=[SimpleNamespace(instance_name='hom_sci_1',products=[methyl,formyl]),
                   SimpleNamespace(instance_name='elimination',products=[ketene,hydrogen],
                                   ts=_species('ts','',root.atom))]
    root.reac_type=['hom_sci','h2_elim'];root.reac_ts_done=[-1,-1]
    plan=build_network_request(root,{'composite_method':'ANL0-F12'})
    original=deepcopy(plan)
    projected=selection.project(plan,{'selected':['hom_sci_1'],'complexes':[]})
    assert set(projected['cbh']['targets'])=={'well','methyl','formyl'}
    assert {r['routing_name'] for r in projected['network']['stable_species']}=={'well','methyl','formyl'}
    assert 'C=C=O' not in {r['smiles'] for r in projected['cbh']['calculated_states']}
    for entry in projected['cbh']['calculated_states']:
        assert entry in plan['cbh']['calculated_states']
    assert plan==original
    # A later revision reactivates the same stored scientific records.
    expanded=selection.project(plan,{'selected':['hom_sci_1','elimination'],'complexes':[]})
    assert expanded==plan


def test_vrc_receives_full_network_while_anl_receives_outer_selection(tmp_path,monkeypatch):
    from kinbot import kb
    from kinbot.anl import automatic
    from test_kb_vrc_selection import _branch_parameters
    monkeypatch.chdir(tmp_path)
    original=SimpleNamespace(reac_ts_done=[-1,-1])
    selected=SimpleNamespace(reac_ts_done=[-1,-2])
    monkeypatch.setattr(selection,'selected_view',lambda *a:selected)
    received=[]
    def prepare(well,*a,**kw):
        assert well is selected
        return {'schema':1,'status':'planned'}
    monkeypatch.setattr(automatic,'prepare_automatic_cbh_anl',prepare)
    def vrc(well,*a):
        assert well is original
        received.append('vrc')
    def anl(well,*a):
        assert well is selected
        received.append('anl')
        return {'status':'complete'}
    monkeypatch.setattr(kb,'_execute_vrc_branch',vrc)
    monkeypatch.setattr(kb,'_execute_automatic_branch',anl)
    kb._post_reaction_branches(original,_branch_parameters(tmp_path),object(),'input.json',True)
    assert sorted(received)==['anl','vrc']
    assert original.reac_ts_done==[-1,-1]


def test_selected_handoff_requires_only_selected_ts_and_does_not_mark_failure(tmp_path,monkeypatch):
    from kinbot.anl import handoff
    from test_anl_handoff import _composite
    point=lambda n:SimpleNamespace(name=n,chemid=n,smiles='C',charge=0,mult=1)
    root=point('a');product=point('b');good=point('good');omitted=point('omitted')
    reactions=[SimpleNamespace(instance_name=n,ts=ts,do_vdW=False,
                  prod_opt=[SimpleNamespace(species=product)]) for n,ts in [('good',good),('omitted',omitted)]]
    root.reac_obj=reactions;root.reac_ts_done=[-1,-1]
    from copy import copy
    view=copy(root);view.reac_ts_done=[-1,-2]
    native=tmp_path/'composite.json';_composite(native)
    manifest=dict(schema=1,mode='anl',strict=True,
                  species={n:{'composite':str(native)} for n in ('a','b')},
                  transition_states={'good':{'composite':str(native),'identity_smiles':'C'}})
    atomic_json(tmp_path/'handoff.json',manifest)
    monkeypatch.setattr(handoff,'routing_name',lambda p:p.name)
    result=handoff.apply_handoff(view,tmp_path/'handoff.json',report_file=tmp_path/'applied.json')
    assert set(result['transition_states'])=={'good'}
    assert not hasattr(omitted,'final_zero_k_energy')
    assert root.reac_ts_done==[-1,-1]


@pytest.mark.parametrize('managed', [False, True])
def test_full_cbh_assembly_preserves_pathway_choice_through_ladder_selection(tmp_path, monkeypatch, managed):
    """Exercise real plan/profiles/CBH/handoff together; replace native execution only."""
    from kinbot.anl import automatic
    from kinbot.anl.validation import molecule_from_smiles
    from test_anl_automatic import _species
    from test_anl_cbh_atct import _table
    root=_species('ethanol','CCO',['C','C','O']+['H']*6,energy=-150.)
    product=_species('ether','COC',root.atom,energy=-149.9)
    root.reac_obj=[SimpleNamespace(instance_name=name,products=[product],
        ts=_species(name,'',root.atom,energy=energy)) for name,energy in [('low',-149.8),('high',-149.7)]]
    root.reac_ts_done=[-1,-1];root.reac_type=['intra_H_migration']*2
    par={'composite_method':'ANL0-F12','anl_ts_reuse':False}
    if managed:
        root.anl_pathway_selection={'selected':['low'],'complexes':[]}
    plan=automatic.prepare_automatic_cbh_anl(root,par,directory=tmp_path)
    monkeypatch.setattr(automatic,'load_atct_references',lambda *a:_table()[1])
    def references(entries,*args):
        result={}
        for entry in entries:
            mol=molecule_from_smiles(entry['smiles'],charge=entry['charge'],multiplicity=entry['multiplicity'])
            obj=_species(entry['smiles'],entry['smiles'],mol['symbols'],entry['multiplicity'])
            obj.charge=entry['charge'];obj.geom=mol['positions']
            record=automatic._l2_species_record(obj,plan['l2_theory'])
            record['thermochemistry_model']={}
            result[automatic._state_key(entry['smiles'],entry['charge'],entry['multiplicity'])]=record
        return result
    monkeypatch.setattr(automatic,'_run_reference_l2',references)
    def native(record, graph):
        # An empty fixture graph lets real assembly read the accepted L2
        # fallback; no licensed calculation or scheduler is used in this test.
        atomic_json(graph/'workflow.json', {'name':'anl-automatic-post-l2','molecule':record['molecule'],'tasks':[]})
        atomic_json(graph/'state.json', {'tasks':{}})
    monkeypatch.setattr(automatic,'_prepare_stable_graph',lambda p,e,r,g:native(r,g))
    monkeypatch.setattr(automatic,'_prepare_ts_graph',lambda p,r,g:native(r,g))
    monkeypatch.setattr(automatic,'preflight_until_scheduler_ready',lambda *a,**kw:{})
    monkeypatch.setattr(automatic,'_drive_graphs',lambda *a:None)
    monkeypatch.setattr(automatic,'_load',lambda p:(p,json.loads((p/'workflow.json').read_text()),json.loads((p/'state.json').read_text())))
    result=automatic.run_automatic_cbh_anl(root,par,'unused',directory=tmp_path,plan=plan)
    expected={'low'} if managed else {'low','high'}
    handoff=json.loads(Path(result['handoff_manifest']).read_text())
    assert set(handoff['transition_states'])==expected
    published=json.loads(Path(result['plan']).read_text())
    assert set(published['network']['transition_states'])==expected
    assert all(Path(p).is_file() for p in result['formations'].values())
    assert root.reac_ts_done==[-1,-1]
    before={p:p.read_bytes() for p in (tmp_path/'states').rglob('*.json')}
    again=automatic.run_automatic_cbh_anl(root,par,'unused',directory=tmp_path,plan=plan)
    assert again==result
    assert before=={p:p.read_bytes() for p in before}


def test_real_workers_wait_for_outer_selection_and_restart_retains_qc(tmp_path):
    import subprocess
    import sys
    import time
    record=phase.initialize(tmp_path)
    for name in ('a','b'):
        atomic_json(tmp_path/name/(name+'.json'),{'composite_method':'ANL1'})
        row=evidence(name,[(name+'_ts',['b' if name=='a' else 'a'],-99. if name=='a' else -99.5,'ordinary')])
        row['generation']=dict(token=record['token'],epoch=0)
        atomic_json(tmp_path/name/'accepted_l2_network.json',row)
    native=tmp_path/'a/anl_automatic/states/retained/state.json'
    atomic_json(native,{'tasks':{'sp':{'status':'submitted','job_id':'retained'}}})
    original=native.read_bytes()
    child='''import json,time
from pathlib import Path
from kinbot import pes_phase as p
assert p.attach()
Path('explored').write_text(str(time.time_ns()))
p.wait_for_exploration()
m=json.loads((Path('..')/'pes_l2_selection.json').read_text())
assert m['status']=='selected'
Path('selected').write_text(json.dumps(m['workers'][Path.cwd().name]['selected']))
p.assembled()
'''
    env={**os.environ,'PYTHONPATH':str(Path(__file__).parents[1]),
         'KINBOT_PES_PHASE_ROOT':str(tmp_path),'KINBOT_PES_PHASE_TOKEN':record['token']}
    processes={}
    try:
        deadline=time.monotonic()+15
        while time.monotonic()<deadline:
            busy,status=phase.tick(tmp_path,['a','b'],max_searches=1)
            if not busy and len(processes)<2:
                name=('a','b')[len(processes)]
                proc=subprocess.Popen([sys.executable,'-c',child],cwd=tmp_path/name,env=env)
                processes[name]=proc;phase.register(tmp_path,name,proc.pid)
            if len(processes)==2 and all(p.poll() is not None for p in processes.values()):
                break
            time.sleep(.03)
        assert all(p.poll()==0 for p in processes.values()) and len(processes)==2
        assert json.loads((tmp_path/'a/selected').read_text())==[]
        assert json.loads((tmp_path/'b/selected').read_text())==['b_ts']
        assert native.read_bytes()==original
        # A new controller generation cannot silently consume an old choice.
        fresh=phase.initialize(tmp_path)
        selection.publish(tmp_path,fresh,['a','b'])
        assert json.loads((tmp_path/selection.NAME).read_text())['status']=='blocked'
        assert native.read_bytes()==original
    finally:
        for p in processes.values():
            if p.poll() is None:p.terminate()
            p.wait(timeout=5)
