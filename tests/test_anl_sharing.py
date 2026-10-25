"""Cross-PES adoption through real catalogs, graph manifests and owner joins."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from kinbot import live_control
from kinbot.anl import sharing, dispatch, automatic, pes_ts_reuse as reuse
from kinbot.execution import atomic_json, file_sha256
from tests.test_anl_live_protocol import specification, fixed_resources
from tests.test_pes_ts_reuse import fixture, add_minimum


@pytest.fixture
def pair(tmp_path):
    live_control.configure(tmp_path)
    cases = [tmp_path/'source', tmp_path/'receiver']
    graphs = []
    for case in cases:
        base = fixture(case)[0]
        name = add_minimum(base)
        spec = specification()
        plan = json.loads((base/'network_request.json').read_text())
        entry = plan['cbh']['calculated_states'][0]
        entry.update(smiles='C', multiplicity=1)
        entry['network_record'].update(smiles='C', molecule=deepcopy(spec['molecule']))
        atomic_json(base/'network_request.json', plan)
        name = 'stable:C|0|1'
        item = next(v for v in reuse.catalog(case).values() if v['name'] == name)
        graph = Path(item['graph'])
        dispatch.prepare(spec, graph)
        graphs.append(graph)
    state = dispatch._load(graphs[0])[2]
    state['tasks']['l3_geometry_tz'].update(status='submitted', job_id='42')
    atomic_json(graphs[0]/'state.json', state)
    return tmp_path, cases, graphs


def test_opt_in_direction_cycles_and_repeat(pair):
    root, (a,b), (source,target) = pair
    assert sharing.discover(target) is None
    sharing.configure(root,a,b)
    assert not (root/sharing.CONFIG).exists()
    first = sharing.configure(root,a,b,apply=True)
    assert first == sharing.configure(root,a,b,apply=True)
    with pytest.raises(ValueError,match='acyclic'):
        sharing.configure(root,b,a,apply=True)
    saved = {p:file_sha256(p) for p in source.rglob('*') if p.is_file()}
    alias = sharing.discover(target)
    assert alias['graph'] == source
    assert sharing.discover(target) == alias
    assert sharing.discover(source) is None
    assert sharing.consumers(source) == [str(target)]
    assert saved == {p:file_sha256(p) for p in saved}


@pytest.mark.parametrize('change', ['charge', 'geometry', 'theory', 'protocol', 'started'])
def test_mismatch_or_existing_receiver_execution_is_not_adopted(pair,change):
    root,(a,b),(source,target)=pair
    sharing.configure(root,a,b,apply=True)
    if change in ('charge','geometry','theory'):
        path=target.parents[2]/'network_request.json'
        plan=json.loads(path.read_text()); entry=plan['cbh']['calculated_states'][0]
        if change=='charge': entry['charge']=1
        elif change=='geometry': entry['network_record']['molecule']['positions'][0][0]+=.01
        else: plan['l2_theory']='different'
        atomic_json(path,plan)
    elif change=='protocol':
        from kinbot.anl import live_control as owner
        _,spec,state=dispatch._load(target)
        spec['intent']['live_protocol']=dict(omitted_tasks=['cfour_dboc'])
        owner.commit(target,spec,state,set(),{})
    else:
        state=dispatch._load(target)[2]
        state['tasks']['l3_geometry_tz'].update(status='failed',job_id='99')
        atomic_json(target/'state.json',state)
    assert sharing.discover(target) is None


def test_changed_evidence_blocks_use_and_detachment_preserves_owner(pair):
    root,(a,b),(source,target)=pair
    sharing.configure(root,a,b,apply=True); sharing.discover(target)
    path=source.parents[2]/'network_request.json'
    old=path.read_bytes(); data=json.loads(old)
    data['cbh']['calculated_states'][0]['network_record']['energy_hartree']+=.01
    atomic_json(path,data)
    with pytest.raises(ValueError,match='provenance'):
        sharing.load(target.parents[2])
    path.write_bytes(old)
    sharing.detach(target,reason='receiver protocol changed')
    assert not sharing.consumers(source)
    assert sharing.load(target.parents[2]) == {}
    assert sharing.discover(target) is None
    assert dispatch._load(source)[2]['tasks']['l3_geometry_tz']['job_id']=='42'


def test_owner_join_keeps_canonical_graph_in_downstream_handoff(pair,monkeypatch):
    root,(a,b),(source,target)=pair
    sharing.configure(root,a,b,apply=True); alias=sharing.discover(target)
    calls=[]
    monkeypatch.setattr(reuse,'finish',lambda row,interval: calls.append(row['graph']) or row['graph'])
    monkeypatch.setattr(automatic,'assemble_available_profiles',lambda *a,**k: {'thermochemistry':{'ANL1':'native.json'},'barrier':{}})
    graphs={}
    external,profiles=automatic._shared_profiles(target.parents[2],graphs,{},.01)
    assert calls==[source] and graphs[alias['name']]==source
    assert profiles[alias['name']]['thermochemistry']['ANL1']=='native.json'
    receipt=json.loads((target.parents[2]/'pes_reuse_consumed.json').read_text())
    assert receipt['graphs'][alias['name']]['graph']==str(source)


def test_live_receiver_revision_detaches_before_native_dispatch(pair,monkeypatch):
    from kinbot.anl import protocol
    root,(a,b),(source,target)=pair
    monkeypatch.setattr(dispatch, 'assign_partitions', fixed_resources)
    sharing.configure(root,a,b,apply=True); sharing.discover(target)
    source_hash=file_sha256(source/'workflow.json')
    requests=protocol.publish(root,[protocol.plan(source,dict(omit_components=['dboc']))],reason='source edit')
    protocol.process(source)
    receipt=json.loads((root/protocol.QUEUE/(requests[0]['id']+'.json')).read_text())
    assert receipt['status']=='queued' and 'receivers' in receipt['error']
    protocol.publish(root,[protocol.plan(target,dict(omit_components=['dboc']))],reason='receiver edit')
    protocol.process(target)
    assert not sharing.consumers(source) and not sharing.load(target.parents[2])
    assert file_sha256(source/'workflow.json')==source_hash
    assert dispatch._load(source)[2]['tasks']['l3_geometry_tz']['job_id']=='42'


def test_catalog_does_not_authorize_a_different_actual_graph(pair):
    root,(a,b),(source,target)=pair
    from kinbot.anl.live_control import commit
    sharing.configure(root,a,b,apply=True)
    _,spec,state=dispatch._load(target)
    spec['molecule']['positions'][0][0]+=.01
    from ase import Atoms
    mol=spec['molecule']
    mol['source']['geometry_sha256']=dispatch._geometry_hash(Atoms(mol['symbols'],positions=mol['positions']))
    commit(target,spec,state,set(),{})
    assert sharing.discover(target) is None
