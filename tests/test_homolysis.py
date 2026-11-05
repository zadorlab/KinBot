"""Discovery, durable scan admission, and rate-export coverage contracts."""
import copy
import json
from pathlib import Path
from types import SimpleNamespace as NS
from unittest.mock import Mock

import numpy as np
import pytest

from kinbot import homolysis
from kinbot.parameters import Parameters
from kinbot.reaction_finder import ReactionFinder
from kinbot.reaction_generator import ReactionGenerator
from kinbot.reactions.reac_dissociation_saddle import DissociationSaddle
from kinbot.search_state import SearchState
from kinbot.execution import atomic_json


@pytest.fixture
def scan(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atomic_json('input.json', {'barrier_threshold':150})
    par = dict(Parameters('input.json').par, homolytic_saddle_scan_points=5, scan_step=30,
               barrier_threshold=150, delete_intermediate_files=0)
    species = NS(geom=np.array([[0., 0., 0.], [1.1, 0., 0.]]), natom=2,
                 linear=[], reac_type=[homolysis.FAMILY], reac_name=['x_dissociation_saddle_1_2'],
                 reac_inst=[[0, 1]], reac_step=[0], reac_ts_done=[0],
                 reac_scan_energy=[[]])
    qc = NS(job_ids={}, check_qc=Mock(return_value=0))
    obj = DissociationSaddle(species, qc, par, [0, 1], species.reac_name[0])
    obj.homolysis_parent='x_hom_sci_1_2'
    species.reac_obj=[obj]
    atomic_json('input.json', {})
    g=NS(species=species, qc=qc, par=par, inp=str(tmp_path/'input.json'))
    g.search_checkpoint=SearchState(g)
    admissions=[]
    def submit(obj, step, command, **kwargs):
        admissions.append(step)
        obj.search_checkpoint.prepare_submission(obj, step)
        if getattr(g, 'pool_full', False):
            obj.search_checkpoint.finish_submission(obj, step)
            return step
        qc.job_ids[obj.instance_name]=str(100+step)
        obj.search_checkpoint.submitted(obj.instance_name, str(100+step))
        obj.search_checkpoint.finish_submission(obj, step+1)
        qc.check_qc.return_value='running'
        return step+1
    monkeypatch.setattr('kinbot.reac_family.carry_out_reaction', submit)
    def complete(energy):
        step=species.reac_step[0]-1
        geom=np.array([[0.,0.,0.],[1.1+step*.1,0.,0.]])
        qc.check_qc.return_value='normal'
        qc.get_qc_energy=Mock(return_value=(0,energy))
        qc.get_qc_geom=Mock(return_value=(0,geom))
    g.admissions, g.complete = admissions, complete
    return g


def test_scan_starts_at_equilibrium_and_uses_absolute_targets(scan):
    obj=scan.species.reac_obj[0]
    for step in (0, 1, 3):
        _,fix,change,_=obj.get_constraints(step,np.zeros((2,3)))
        assert change==[[1,2,pytest.approx(1.1+step*.1)]] and fix==[]
    assert obj.max_step==5


def test_scan_adopts_live_and_staged_work_and_peak_survives_pool_full(scan):
    g=scan
    assert homolysis.advance(g,0)
    assert g.admissions==[0]
    for _ in range(2):homolysis.advance(g,0)  # live: no new rendering/admission
    assert g.admissions==[0]
    for energy in (-10.,-9.,-9.5):
        g.complete(energy)
        if energy==-9.5:g.pool_full=True
        homolysis.advance(g,0)
    assert g.admissions==[0,1,2,5]
    assert g.species.reac_step==[5]
    g.search_checkpoint.save(0)
    g.search_checkpoint=SearchState(g)  # controller replacement during admission
    g.species.reac_obj[0].search_checkpoint=g.search_checkpoint
    homolysis.advance(g,0)
    assert len(homolysis.record(g.species.reac_obj[0])['points'])==3
    assert g.species.reac_step==[5]
    assert np.array(g.species.reac_obj[0].scan_saddle_geometry)[1,0]==pytest.approx(1.2)
    g.pool_full=False
    homolysis.advance(g,0)
    assert g.species.reac_step==[6]
    assert homolysis.advance(g,0) is False  # native TS/frequency/IRC machinery owns this


def test_new_point_not_visible_does_not_advance(scan):
    homolysis.advance(scan,0)
    scan.qc.check_qc.return_value=0
    homolysis.advance(scan,0)
    assert scan.admissions==[0] and scan.species.reac_step==[1]


@pytest.mark.parametrize('energies', [[-10,-9,-8,-7,-6],[-10]*5])
def test_complete_finite_scan_is_not_a_barrierless_certificate(scan,energies):
    homolysis.advance(scan,0)
    for e in energies:
        scan.complete(e);homolysis.advance(scan,0)
    assert scan.species.reac_ts_done==[-999]
    r=homolysis.record(scan.species.reac_obj[0])
    assert r['outcome']=='no_maximum_on_sampled_path'
    assert len(r['points'])==5 and not r['proves_global_absence_of_saddle']
    assert scan.admissions==[0,1,2,3,4]  # no fake final TS/dummy H database row


def test_failed_or_wrong_coordinate_point_remains_inconclusive(scan):
    homolysis.advance(scan,0)
    scan.complete(-10)
    scan.qc.get_qc_geom.return_value=(0,np.array([[0,0,0],[8,0,0]]))
    homolysis.advance(scan,0)
    assert homolysis.record(scan.species.reac_obj[0])['outcome']=='scan_constraint_not_met'
    assert scan.species.reac_ts_done==[-999]


def legacy_wrong_target(g):
    homolysis.advance(g, 0)
    g.complete(-10.)
    homolysis.advance(g, 0)
    g.complete(-9.)
    guess = [[0.,0.,0.],[1.19887,0.,0.]]  # requested point is 1.20000 A
    g.qc.get_qc_geom.return_value = (0,np.array(guess))
    obj = g.species.reac_obj[0]
    Path(obj.instance_name+'.py').write_text(
        f'input_geom = np.asarray({guess!r}, dtype=float)\n'
        'base_0_fix = [[idx - 1 for idx in fix] for fix in [[1,2]]]\n'
        'strict_scan = True\n')
    homolysis.advance(g, 0)
    g.search_checkpoint.entries[obj.homolysis_parent] = dict(status=-1)
    g.search_checkpoint.save(0)
    return obj


@pytest.mark.parametrize('interrupt_before_checkpoint', [False, True])
def test_legacy_target_repair_resumes_only_bad_point_once(scan, interrupt_before_checkpoint):
    g=scan; obj=legacy_wrong_target(g)
    previous=homolysis.record(obj)['points']
    assert g.species.reac_ts_done == [-999]
    if interrupt_before_checkpoint:
        assert homolysis.repair_constraint_target(g.qc,obj,g.search_checkpoint.entries[obj.instance_name]) == 1
    g.search_checkpoint=SearchState(g)
    assert g.species.reac_ts_done == [0] and g.species.reac_step == [1]
    assert homolysis.record(obj)['points'] == previous
    repair=homolysis.record(obj)['constraint_target_repair']
    assert all(Path(f['archive']).is_file() for f in repair['evidence'].values())
    homolysis.advance(g,0)
    assert g.admissions == [0,1,1] and g.species.reac_step == [2]
    for _ in range(2):
        g.search_checkpoint=SearchState(g)
        homolysis.advance(g,0)
    assert g.admissions == [0,1,1]  # adopt live repaired point
    g.complete(-9.)
    homolysis.advance(g,0)
    assert g.admissions == [0,1,1,2]
    assert homolysis.record(obj)['points'][:1] == previous
    assert len(homolysis.record(obj)['points']) == 2
    # A later failure is not another entitlement to repeat exhausted settings.
    g.complete(-8.)
    g.qc.get_qc_geom.return_value = (0,np.array([[0.,0.,0.],[8.,0.,0.]]))
    homolysis.advance(g,0); g.search_checkpoint.save(0)
    g.search_checkpoint=SearchState(g)
    assert g.species.reac_ts_done == [-999]


@pytest.mark.parametrize('problem', ['live', 'qc_error', 'explicit_target', 'different_frozen', 'unknown', 'parent_rejected'])
def test_target_repair_requires_specific_native_proof(scan, problem):
    g=scan; obj=legacy_wrong_target(g)
    if problem in ('live','qc_error','unknown'):
        g.qc.check_qc.return_value = {'live':'running','qc_error':'error','unknown':0}[problem]
    elif problem == 'explicit_target':
        with Path(obj.instance_name+'.py').open('a') as stream:
            stream.write('constraint_targets = [[1,2,1.2]]\n')
    elif problem == 'parent_rejected':
        g.search_checkpoint.entries[obj.homolysis_parent]['status'] = -999
        g.search_checkpoint.write()
    else:
        g.qc.get_qc_geom.return_value = (0,np.array([[0.,0.,0.],[1.17,0.,0.]]))
    g.search_checkpoint=SearchState(g)
    assert g.species.reac_ts_done == [-999]
    assert not Path('.kinbot_scan_repairs').exists()


def test_discovery_creates_companions_only_after_user_selection(tmp_path):
    from tests.test_reaction_paths import peroxy
    p=peroxy()
    atomic_json(tmp_path/'input.json', {'barrier_threshold':150})
    par=dict(Parameters(tmp_path/'input.json').par, families=['hom_sci'],
             homolytic_bonds={str(p.chemid):[[0,1]]})
    ReactionFinder(p,par,None).find_reactions()
    assert p.reac_type==['hom_sci',homolysis.FAMILY]
    assert p.reac_obj[1].homolysis_parent==p.reac_name[0]
    assert p.reac_inst==[[0,1],[0,1]]
    assert 'hom_sci' not in p.reac_name[1]  # Gaussian must retain its TS optimization


def test_ordinary_discovery_precedes_check_and_matching_ts_is_reused(scan,monkeypatch):
    g=scan; s=g.species; obj=s.reac_obj[0]
    parent=NS(instance_name=obj.homolysis_parent,instance=[0,1],products=[NS(key='p'),NS(key='q')])
    s.bond=np.array([[0,1],[1,0]])
    ts=NS(instance_name='ordinary',products=list(reversed(parent.products)),irc_prod=NS(bond=np.zeros((2,2))))
    s.reac_obj.extend([parent,ts]);s.reac_type.extend(['hom_sci','ordinary'])
    s.reac_name.extend([parent.instance_name,ts.instance_name])
    s.reac_ts_done.extend([-1,1]);s.reac_step.extend([0,13])
    monkeypatch.setattr(homolysis,'routing_key',lambda p:p.key)
    monkeypatch.setattr('kinbot.reaction_path.reaction_path_id',lambda o:'validated-path')
    assert homolysis.deferred(g,0)  # general search still active
    assert not g.admissions
    s.reac_ts_done[2]=-1
    assert homolysis.deferred(g,0)
    assert s.reac_ts_done[0]==-999 and not g.admissions
    assert homolysis.record(obj)['saddles']==['ordinary']
    g.search_checkpoint.save(0)  # Normal generator checkpoint after dispatch.
    homolysis.finish(s)
    assert json.loads(Path('homolysis_checks.json').read_text())['checks'][0]['outcome']=='superseded_by_validated_saddle'
    s.reac_ts_done[2]=-999
    homolysis.finish(s)
    assert homolysis.unresolved_checks()[0]['outcome'] == 'no_validated_saddle'
    homolysis.ensure_export_ready()  # A terminal nonfinding is not an export hold.
    assert 'not proof of barrierlessness' in homolysis.export_notes()[0]


def test_real_generator_scan_exhaustion_retains_outcome_without_false_acceptance(scan,monkeypatch):
    g=scan;s=g.species
    # The fixture isolates the scan owner; defer and final review are covered
    # above. All actual stage-0 dispatch, checkpoints and loop logic run here.
    monkeypatch.setattr(homolysis,'deferred',lambda *a:False)
    monkeypatch.setattr(homolysis,'finish',lambda *a:None)
    monkeypatch.setattr('kinbot.reaction_generator.time.sleep',lambda *a:None)
    original=homolysis.advance
    def advance(owner,i):
        if s.reac_step[0]>0:g.complete(-10+s.reac_step[0])
        return original(owner,i)
    monkeypatch.setattr(homolysis,'advance',advance)
    generator=ReactionGenerator(s,g.par,g.qc,g.inp)
    generator.generate()
    assert s.reac_ts_done==[-999]
    assert homolysis.record(s.reac_obj[0])['outcome']=='no_maximum_on_sampled_path'
    assert g.admissions==[0,1,2,3,4]


def test_pending_revision_blocks_stale_pes_mess_before_reading_qc(tmp_path,monkeypatch):
    from kinbot import pes
    from kinbot.mess import MESS
    monkeypatch.chdir(tmp_path)
    atomic_json('homolysis_checks.json',dict(checks=[dict(parent='a',outcome='pending')]))
    with pytest.raises(RuntimeError,match='old summaries'):
        MESS.write_input(NS(),None)
    with pytest.raises(RuntimeError,match='old summaries'):
        pes.postprocess({},['.'],'all',[],1)


def test_l1_scan_and_l2_saddle_keep_requested_theories(tmp_path,monkeypatch):
    from kinbot.qc import QuantumChemistry
    monkeypatch.chdir(tmp_path)
    atomic_json('input.json',dict(barrier_threshold=150,method='B3LYP',basis='cc-pVDZ',
        high_level_method='B2PLYP',high_level_basis='cc-pVTZ',scan_method='AM1',scan_basis='',
        queuing='local',qc='gauss',mult=2))
    qc=QuantumChemistry(Parameters('input.json').par)
    name='well_dissociation_saddle_1_2'
    for step in (0,1,40):
        args=qc.get_qc_arguments(name,2,0,18,ts=True,scan=True,step=step,max_step=40)
        assert (args['method'],args['basis'])==('B3LYP','cc-pVDZ')
        assert 'opt' in args
    args=qc.get_qc_arguments(name+'_high',2,0,18,ts=True,high_level=True)
    assert (args['method'],args['basis'])==('B2PLYP','cc-pVTZ')


def test_accepted_tight_path_supersedes_placeholder_using_exporter_endpoint_rule(monkeypatch):
    parent=NS(instance=[0,1],products=[NS(key='p'),NS(key='q')])
    bond=np.array([[0,1,0],[1,0,1],[0,1,0]])
    other=NS(instance_name='rearrangement',products=parent.products,
             irc_prod=NS(bond=np.array([[0,0,1],[0,0,0],[1,0,0]])))
    species=NS(bond=bond,reac_obj=[other],reac_type=['migration'],reac_ts_done=[-1])
    monkeypatch.setattr(homolysis,'routing_key',lambda p:p.key)
    monkeypatch.setattr('kinbot.reaction_path.reaction_path_id',lambda o:'ordinary')
    assert homolysis.matching_saddles(species,parent)==['rearrangement']
    # A different configured fragment or an unaccepted TS does not replace it.
    other.products = [NS(key='p'), NS(key='q_other_state')]
    assert homolysis.matching_saddles(species,parent)==[]
    other.products = parent.products
    species.reac_ts_done[0] = -999
    assert homolysis.matching_saddles(species,parent)==[]


@pytest.mark.parametrize('outcome', ['qc_failed', 'scan_constraint_not_met',
    'no_maximum_on_sampled_path', 'saddle_submitting', 'geometry_or_submission_failed'])
def test_terminal_attempts_export_with_notes_without_changing_native_history(tmp_path, outcome):
    path = tmp_path/'homolysis_checks.json'
    atomic_json(path, dict(checks=[dict(parent='x_hom_sci_1_2', search='scan',
        outcome='unresolved', scan_outcome=outcome)]))
    state = tmp_path/'.kinbot_search_state.json'
    atomic_json(state, dict(reactions={'scan':dict(status=-999, step=41)}))
    before = [p.read_bytes() for p in (path, state)]
    for _ in range(2):
        homolysis.ensure_export_ready(tmp_path)
        assert outcome in homolysis.export_notes(tmp_path)[0]
    assert [p.read_bytes() for p in (path, state)] == before
    # A live/reopened search, unlike a terminal scan's last intent, must finish.
    atomic_json(state, dict(reactions={'scan':dict(status=1, step=41)}))
    with pytest.raises(RuntimeError, match='unfinished homolysis'):
        homolysis.ensure_export_ready(tmp_path)


@pytest.mark.parametrize('tail,expected', [([-2.,-2.,-2.],'supports_provisional_pst'),
    ([-2.1,-2.1,-2.1],'unresolved'),([-2.1,-2.05,-2.],'unresolved')])
def test_asymptote_requires_matching_state_flatness_and_energy(tail,expected,tmp_path,monkeypatch):
    from kinbot.stationary_pt import StationaryPoint
    monkeypatch.chdir(tmp_path)
    a=StationaryPoint('hydrogen',0,2,atom=['H'],geom=np.zeros((1,3)));a.characterize()
    pair=StationaryPoint('pair',0,1,atom=['H','H'],geom=np.array([[0.,0.,0.],[5.,0.,0.]]))
    parent=NS(products=[a,a])
    pair.reac_name=['hom'];pair.reac_obj=[parent]
    obj=NS(species=pair,homolysis_parent='hom')
    par={'homolytic_saddle_asymptote_tolerance':.5,'homolytic_saddle_flatness':.1}
    monkeypatch.setattr('kinbot.electronic_states.apply_electronic_state_model',lambda *a:None)
    monkeypatch.setattr('kinbot.l1_screening.theory',lambda p:{})
    monkeypatch.setattr('kinbot.l1_screening.read_points',lambda *a,**k:[{'energy_hartree':-1.}]*2)
    data={'points':[dict(energy_hartree=e,positions=pair.geom.tolist()) for e in tail]}
    result=homolysis.sampled_asymptote(NS(species=pair,par=par),obj,data)
    assert result['status']==expected
    # Same energies cannot license a different set of fragments.
    oxygen=StationaryPoint('oxygen',0,3,atom=['O'],geom=np.zeros((1,3)));oxygen.characterize()
    parent.products=[oxygen]
    assert homolysis.sampled_asymptote(NS(species=pair,par=par),obj,data)['status']=='unresolved'


@pytest.mark.parametrize('converged',[True,False])
def test_generated_strict_scan_worker_requires_optimizer_convergence(tmp_path,monkeypatch,converged):
    from ase.calculators.emt import EMT
    from ase.db import connect
    from kinbot import kb_path
    import kinbot.ase_modules.calculators.gaussian as gaussian
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(gaussian,'Gaussian',lambda **kw:EMT())
    calls=[]
    def optimize(opt,**kw):
        calls.append(kw)
        # The requested target must survive an approximate initial geometry.
        assert opt.pes.cons._targets['bonds'] == [pytest.approx(.7412)]
        return converged
    monkeypatch.setattr('kinbot.optimization_budget.run_optimizer',optimize)
    template=Path(kb_path,'tpl/ase_sella_ts_search.tpl.py').read_text()
    script=template.format(working_dir=str(tmp_path),label='point',code='gaussian',Code='Gaussian',
        atom=['H','H'],geom=[[0.,0.,0.],[.74,0.,0.]],fix=[[1,2]],kwargs={},sella_kwargs={'internal':False},
        strict_scan=True,fmax=.0005,steps=80,constraint_targets=[[1,2,.7412]])
    if converged:
        exec(compile(script,'generated_scan','exec'),{})
        row=connect('kinbot.db').get(name='point')
        assert row.data['status']=='normal' and np.isfinite(row.data['energy'])
        assert Path('point_sella.log').read_text().endswith('done\n')
    else:
        with pytest.raises(RuntimeError,match='did not converge'):
            exec(compile(script,'generated_scan','exec'),{})
        assert not list(connect('kinbot.db').select(name='point'))
    assert calls==[dict(fmax=.0005,steps=80)]


def test_real_sella_scan_converges_to_requested_not_modified_distance(tmp_path, monkeypatch):
    from ase.calculators.emt import EMT
    from ase.db import connect
    from kinbot import kb_path
    import kinbot.ase_modules.calculators.gaussian as gaussian
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(gaussian, 'Gaussian', lambda **kw: EMT())
    script=Path(kb_path,'tpl/ase_sella_ts_search.tpl.py').read_text().format(
        working_dir=str(tmp_path),label='point',code='gaussian',Code='Gaussian',
        atom=['H','H'],geom=[[0.,0.,0.],[.739,0.,0.]],fix=[[1,2]],kwargs={},
        sella_kwargs={'internal':False},strict_scan=True,fmax=.0005,steps=80,
        constraint_targets=[[1,2,.7412]])
    exec(compile(script,'generated_scan','exec'),{})  # actual Sella, local analytic QC substitute
    row=connect('kinbot.db').get(name='point')
    assert row.data['status']=='normal'
    assert row.toatoms().get_distance(0,1) == pytest.approx(.7412, abs=1.e-5)
