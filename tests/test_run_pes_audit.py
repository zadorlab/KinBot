from copy import deepcopy
import json

from ase import Atoms
from ase.db import connect
import pytest

from kinbot import run_pes_audit as audit
from kinbot.anl import ts_reuse
from kinbot.execution import atomic_json
from test_network_revision import saddle_plan


def case_fixture(tmp_path):
    case = tmp_path/'case'
    for worker in ('111','222'):
        d = case/worker
        atomic_json(d/(worker+'.json'),dict(smiles='C[S]',charge=0,mult=2))
        plan = saddle_plan()
        name = worker+'_reaction'
        record = deepcopy(plan['network']['transition_states']['a'])
        plan['network']['transition_states'] = {name:record}
        plan['network']['reactions'] = [dict(name=name,transition_state=name,
            reactant=worker,products=['222' if worker=='111' else '111'])]
        atomic_json(d/'anl_automatic/network_request.json',plan)
        (d/'kinbot_monitor.out').write_text('-1 5 '+name+' '+('222' if worker=='111' else '111')+'\n')
    return case


def test_full_topology_includes_reverse_and_unmaterialized_workers(tmp_path):
    case = case_fixture(tmp_path)
    p = case/'222/anl_automatic/network_request.json'
    plan = json.loads(p.read_text())
    # A snapshot is available before its CBH graph materializes.
    atomic_json(case/'222/accepted_l2_network.json',dict(network=plan['network']))
    p.unlink()
    before = {p:p.read_bytes() for p in case.rglob('*') if p.is_file()}
    result = audit.audit_case(case)
    assert result['accepted_search_routes']==2 and result['endpoint_pair_count']==1
    assert result['comparisons'][0]['equivalent']
    assert not result['comparisons'][0]['reusable']
    assert result['comparisons'][0]['reverse']
    assert '2 searches' in audit.dot(result)
    assert before=={p:p.read_bytes() for p in case.rglob('*') if p.is_file()}


def test_old_plan_does_not_resurrect_rejected_current_reaction(tmp_path):
    case = case_fixture(tmp_path)
    (case/'222/kinbot_monitor.out').write_text('-999 5 222_reaction 111\n')
    result = audit.audit_case(case)
    assert result['accepted_search_routes']==1 and not result['comparisons']


def test_latest_database_failure_does_not_reuse_old_success(tmp_path):
    database = tmp_path/'kinbot.db'
    db = connect(database)
    mol = Atoms('SH',positions=[[0,0,0],[0,0,1.3]])
    db.write(mol,name='rx_high',data=dict(status='normal',energy=-100.,zpe=.01,frequencies=[-100.]))
    before = database.read_bytes()
    inp = dict(charge=0,mult=2)
    found = audit.database_record(database,'rx_high',inp,'[SH]')
    assert found['source_row_id']==1
    assert before==database.read_bytes()
    db.write(mol,name='rx_high',data=dict(status='error'))
    assert audit.database_record(database,'rx_high',inp,'[SH]') is None


def test_comparison_explains_but_does_not_loosen_thresholds():
    p = saddle_plan()
    p['network']['transition_states']['b']['energy_hartree']+=1.e-4
    result = ts_reuse.comparison(p,'a','b',ts_reuse.options({}))
    assert result['reason']=='numerical_tolerance_exceeded'
    assert result['exceeded']==['energy_hartree']
    assert ts_reuse.equivalent(p,'a','b',ts_reuse.options({})) is None
    p['network']['reactions'][1]['products']=['different']
    assert ts_reuse.comparison(p,'a','b',ts_reuse.options({}))['reason']=='different_endpoints'


def test_database_open_is_read_only(tmp_path):
    database = tmp_path/'missing.db'
    with pytest.raises(Exception):
        audit.database_record(database,'rx',{},'C')
    assert not database.exists()


def test_pre_downstream_snapshot_does_not_require_successful_cbh(monkeypatch):
    from types import SimpleNamespace
    from kinbot import pes_snapshot
    from test_anl_automatic import _species
    from kinbot.anl import automatic
    point = _species('111','C[S]',['C','S','H','H','H'],2)
    product = _species('222','C=[SH]',['C','S','H','H','H'],2)
    ts = _species('111_rx','',['C','S','H','H','H'],2)
    point.reac_obj = [SimpleNamespace(instance_name='111_rx',products=[product],ts=ts)]
    point.reac_ts_done = [-1]
    point.reac_type = ['test']
    monkeypatch.setattr(automatic,'_resolved_l2_theory',lambda p:{'method':'test'})
    monkeypatch.setattr('kinbot.reaction_path.reaction_path_id',lambda r:'ordinary')
    monkeypatch.setattr(automatic,'build_network_request',lambda *a:pytest.fail('Must not construct CBH'))
    result = pes_snapshot.snapshot(point,{})
    assert result['network']['reactions'][0]['path_class']=='ordinary'
    assert set(result['network']['transition_states'])=={'111_rx'}
    assert result['errors']==[]


def test_metadata_only_upgrade_retains_existing_qc_files(tmp_path,monkeypatch):
    from kinbot.anl import automatic
    old = saddle_plan()
    for row in old['network']['reactions']:
        row['name'] = row['transition_state']
    current = deepcopy(old)
    for row in current['network']['reactions']:
        row['path_class']='ordinary'
    atomic_json(tmp_path/'network_request.json',old)
    qc = tmp_path/'states/ts_a/composite/tasks/expensive/output'
    qc.parent.mkdir(parents=True)
    qc.write_text('expensive result')
    monkeypatch.setattr(automatic,'build_network_request',lambda *a:current)
    assert automatic.prepare_automatic_cbh_anl(None,{},directory=tmp_path)==current
    assert qc.read_text()=='expensive result'
    assert len(list((tmp_path/'network_metadata_history').glob('*.json')))==1
    assert json.loads((tmp_path/'network_request.json').read_text())==current


def test_failure_inventory_is_recorded_evidence_not_a_live_status(tmp_path):
    (tmp_path/'job.exitcode').write_text('1\n')
    (tmp_path/'job_sella.log').write_text('Worker failure: Optimization did not converge within 250 steps.\n')
    found = audit.failure_inventory(tmp_path,[])
    assert found[0]['recorded_reason'].startswith('Optimization did not converge')
    assert 'status' not in found[0]
