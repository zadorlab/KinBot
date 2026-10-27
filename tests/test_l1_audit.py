from copy import deepcopy
import json

from ase import Atoms
from ase.db import connect
import pytest

from kinbot import constants, l1_audit as audit, l1_screening as l1
from kinbot.execution import atomic_json
from test_l1_screening import case


def point(energy, symbols=('H','H')):
    return dict(symbols=list(symbols), energy_hartree=energy, zpe_hartree=0.)


def test_vdw_is_strict_depth_not_absolute_energy():
    fragments = [point(-.5, ['H'])]*2
    for depth, expected in [(0.,'exclude'), (.49,'exclude'), (.51,'retain')]:
        result = audit.depth_decision(point(-1-depth/constants.AUtoKCAL), fragments, .5)
        assert result['decision'] == expected
        assert result['depth_kcal_mol'] == pytest.approx(depth)
    assert audit.depth_decision(point(-1, ['C']), fragments, .5)['decision'] == 'unknown'


def worker(name, targets=(), height='retain'):
    return dict(worker=name, well=dict(decision=height), work=dict(graphs=[]), channels=[
        dict(reaction=name+'_'+target, accepted_l2=True, operator_rejected=False,
             product_ids=[target], direct_decision=decision, vdw=dict(decision='not_applicable'),
             vrc=dict(kind='homolytic', existing=True, review=dict(request_status='needs_review')))
        for target, decision in targets])


def test_alternate_reversible_paths_unknown_edges_and_fragment_stoichiometry():
    rows = [worker('a',[('b','exclude'),('c','retain'),('d','exclude')]),
            worker('b',[('c','retain')]), worker('c',[('e','unknown')]), worker('d'), worker('e'),
            worker('f',height='exclude'), worker('g')]
    rows[0]['channels'].append(dict(deepcopy(rows[0]['channels'][1]), product_ids=['g','H']))
    audit.reachability(rows,'a')
    statuses={w['worker']:w['reachability'] for w in rows}
    assert statuses['b']=='reachable_with_retained_edges'  # c -> b uses reverse accepted path
    assert statuses['d']=='disconnected_by_screen_in_saved_graph'
    assert statuses['e']=='reachable_with_unknown_edges'
    assert statuses['f']=='above_limit'
    assert statuses['g']=='not_connected_in_saved_graph'  # no single-fragment hop


def test_high_parent_removes_all_local_vrc_but_shallow_complex_keeps_tight_reaction():
    rows=[worker('a',[('b','retain')]),worker('b',[('c','retain')],height='exclude')]
    channel=rows[0]['channels'][0]
    channel['vrc']['kind']='vdw_exit'
    channel['vdw']=dict(decision='exclude',depth_kcal_mol=.2)
    audit.reachability(rows,'a')
    assert channel['decision']=='retain'
    assert channel['vrc']['decision']=='exclude'
    assert rows[1]['channels'][0]['vrc']['decision']=='exclude'


def test_report_readonly_depth_review_plan_and_work_counts(tmp_path,monkeypatch):
    w,par=case(tmp_path)
    routes=json.loads((w/'.kinbot_theory_jobs.json').read_text())
    routes['jobs'].update({'100_exit':dict(level='l1'),'100_IRC_F':dict(level='l1')})
    atomic_json(w/'.kinbot_theory_jobs.json',routes)
    db=connect(w/'kinbot.db')
    db.write(Atoms('H2',positions=[[0,0,0],[0,0,.7]]),name='100_well',
             data=dict(status='normal',energy=-1/constants.EVtoHARTREE,zpe=0.))
    db.write(Atoms('H2'),name='100_exit',data=dict(status='normal',energy=-.95/constants.EVtoHARTREE,zpe=0.))
    db.write(Atoms('H'),name='H_well',data=dict(status='normal',energy=-.48/constants.EVtoHARTREE,zpe=0.))
    endpoint=dict(point(-.96-.2/constants.AUtoKCAL),name='100_exit_IRC_F_prod',row_id=9)
    monkeypatch.setattr(audit,'product_side',lambda *a:(endpoint,['H','H']))
    atomic_json(w/'accepted_l2_network.json',dict(network=dict(reactions=[
        dict(name='100_exit',family='test',products=['H','H'],product_complex='complex') ])))
    atomic_json(w/'vrctst/orbital_review/100_exit/request.json',dict(status='needs_review',bond=[1,2]))
    atomic_json(w/'vrctst/molpro_scan/100_exit__one/state.json',dict(tasks={'x':dict(status='complete')}))
    plan=tmp_path/'reviews/screens.json'
    atomic_json(plan,dict(screens=[dict(name='probe',case=tmp_path.name,recipe='probe.json')]))
    atomic_json(plan.with_name('probe.json'),dict(request_match=dict(parent='100',reaction='100_exit')))
    before={str(p):p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    result=audit.report(tmp_path,150,.5,[plan])
    c=result['workers'][0]['channels'][0]
    assert c['decision']=='retain'
    assert c['vrc']['decision']=='exclude'
    assert c['vdw']['depth_kcal_mol']==pytest.approx(.2)
    assert result['review_plans'][0]['vrc_decision']=='exclude'
    assert result['counts']['recorded_task_statuses']=={'complete':1}
    assert result['workers'][0]['work']['graphs'][0]['dedicated_vrc_reactions']==['100_exit']
    assert before=={str(p):p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    audit.summarize(result)


def test_incomplete_endpoint_preserves_low_complex_when_fragments_high(tmp_path,monkeypatch):
    w,par=case(tmp_path)
    atomic_json(w/'accepted_l2_network.json',dict(network=dict(reactions=[
        dict(name='exit',family='test',products=['H','H'],product_complex='complex') ])))
    result=audit.report(tmp_path)
    c=result['workers'][0]['channels'][0]
    assert c['isolated_products']['decision']=='exclude'
    assert c['products']['decision']=='unknown'
    assert c['vrc']['decision']=='unknown'


def test_missing_initial_reference_is_reported_not_used_as_zero(tmp_path):
    w,par=case(tmp_path)
    (w/'.kinbot_theory_jobs.json').unlink()
    result=audit.report(tmp_path)
    assert result['reference'] is None
    assert result['workers'][0]['well']['decision']=='unknown'
    assert result['issues']


def test_irc_child_uses_recorded_parent_route_only(tmp_path):
    w,par=case(tmp_path)
    routes=json.loads((w/'.kinbot_theory_jobs.json').read_text())
    routes['jobs']['100_IRC_F']=dict(level='l1')
    atomic_json(w/'.kinbot_theory_jobs.json',routes)
    db=connect(w/'kinbot.db')
    db.write(Atoms('H2'),name='100_IRC_F_prod',data=dict(status='normal',energy=-1,zpe=.1))
    assert l1.read_points(w,['100_IRC_F_prod'],l1.theory(par))
    routes['jobs']['100_IRC_F']['level']='l2'
    atomic_json(w/'.kinbot_theory_jobs.json',routes)
    with pytest.raises(ValueError,match='unverified'):
        l1.read_points(w,['100_IRC_F_prod'],l1.theory(par))


def test_real_product_side_identity_selects_opposite_endpoint(tmp_path):
    w,par=case(tmp_path)
    routes=json.loads((w/'.kinbot_theory_jobs.json').read_text())
    routes['jobs'].update({n:dict(level='l1') for n in ('rxn_IRC_F','rxn_IRC_R')})
    atomic_json(w/'.kinbot_theory_jobs.json',routes)
    db=connect(w/'kinbot.db')
    for side,distance in [('F',4.),('R',.7)]:
        db.write(Atoms('H2',positions=[[0,0,0],[0,0,distance]]),name='rxn_IRC_'+side+'_prod',
            data=dict(status='normal',energy=-1,zpe=0.))
    parent=audit.species(dict(point(-1),name='parent',positions=[[0,0,0],[0,0,.7]]),par)
    chosen,fragments=audit.product_side(w,'rxn',par,l1.theory(par),parent)
    assert chosen['name']=='rxn_IRC_F_prod'
    assert len(fragments)==2
    db.write(Atoms('H2',positions=[[0,0,0],[0,0,.7]]),name='rxn_IRC_F_prod',
        data=dict(status='normal',energy=-1,zpe=0.))
    with pytest.raises(ValueError,match='exactly one'):
        audit.product_side(w,'rxn',par,l1.theory(par),parent)
