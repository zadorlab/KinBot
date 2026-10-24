"""Finite follow-up decks and owner sequencing; no native QC is simulated as science."""
import fcntl
import importlib.util
import json
import shutil
from pathlib import Path
import sys
import threading

import numpy as np
import pytest

from kinbot.execution import atomic_json
from tests.test_requested_orbital_recon import completed_screen, helper as screens
from tests.test_vrc_relaxation_recon import helper as relaxation

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT/'examples/anl/orbital_review_batches/r7_progress_20261024'
spec = importlib.util.spec_from_file_location('followup', ROOT/'examples/anl/tools/continue_ch2_sh_s_review.py')
followup = importlib.util.module_from_spec(spec)
sys.modules.setdefault('requested_orbital_recon', screens)
sys.modules.setdefault('vrc_relaxation_recon', relaxation)
spec.loader.exec_module(followup)


@pytest.fixture
def tree(tmp_path):
    production = tmp_path/'production'
    review = tmp_path/'r7_progress_review_20261024'
    recipe = PLAN/'ch2_sh_s_cs.json'
    match = json.loads(recipe.read_text())['request_match']
    case = production/'ch3ss/input.json'
    request = case.parent/match['parent']/'vrctst/orbital_review'/match['reaction']/'request.json'
    xyz = [[1.8122818351,.0029772478,-.0010645061], [2.2322433887,1.7849627863,-.0036765215],
           [.0180913054,.0086000664,-.0026139346], [1.7933644276,2.0266677332,-1.2480033488],
           [2.2475352842,-.4907718979,-.8623412662], [2.2307468144,-.3955252319,.9183373719]]
    atomic_json(request, dict(match,schema=1,status='needs_review',positions=xyz))
    atomic_json(case, dict(queue_name='test-partition'))
    atomic_json(production/'.kinbot_resource_policy.json',dict(schema=1,max_jobs=22,
        excluded_nodes=['badnode'],manage_pending=True,orbital_review_reserve=3))
    source = review/'screens/ch2_sh_s_cs'
    completed_screen((request,case,source),recipe)
    atomic_json(case,dict(queue_name='test-partition',smiles='[H][S][C]([H])([H])[S]',mult=2,
        vrc_tst_scan_basis='avdz-f12',vrc_tst_scan_shift=.1,vrc_tst_scan_ipea=.25,
        vrc_tst_scan_gradient='numerical',vrc_tst_scf_convergence=1e-8,
        vrc_tst_correlation_convergence=1e-7,vrc_tst_scan_fmax=.0005,barrier_threshold=400.))
    destination = review/'ch2_sh_s_followup'
    destination.mkdir()
    return production,review,destination,source,request,case


def test_exact_candidate_and_relaxation_panels_preserve_source_and_reuse_staging(tree,monkeypatch):
    production,review,destination,source,request,case = tree
    before={p:p.read_bytes() for base in (source,production) for p in base.rglob('*') if p.is_file()}
    monkeypatch.setattr(screens.subprocess,'run',lambda *a,**k:pytest.fail('Preparation must not submit'))
    first=followup.stage(production,review,destination)
    assert followup.stage(production,review,destination)==first
    assert len(first['graphs'])==10
    record=followup.read(destination/'candidates/review_launch.json')
    assert record['max_jobs']==3 and record['candidate_comparison']
    for system in record['systems']:
        m=followup.read(system['manifest'])
        assert len(m['jobs'])==12 and m['slurm_excluded_nodes']==['badnode']
        for job in m['jobs']:
            compact=job['probe_active_space']['active_electrons']==5
            deck=(Path(job['directory'])/(job['name']+'.inp')).read_text()
            assert ('closed,18' if compact else 'closed,17') in deck
            assert ('occ,22' if compact else 'occ,23') in deck
            assert 'wf,41,1,1\nstate,2' in deck
            assert deck.count('{rs2-f12,')==2 and 'uccsd' not in deck
            assert job['initial_record']['basis']==job['basis']
            assert not job.get('job_id')
    by_space={4:[],6:[]}
    for graph in first['graphs']:
        s=followup.read(Path(graph)/'workflow.json');t=s['tasks'][0]
        proof=s['intent']['relaxation_recon'];xyz=np.array(s['molecule']['positions'])
        distance=proof['distance_angstrom']
        assert np.linalg.norm(xyz[0]-xyz[1])==pytest.approx(distance)
        assert t['optimizer']['constraints']=={'fixed_bonds':[[0,1]]}
        indices=t['profile']['calculator_kwargs']['state_model']['orbital_selection']['active_indices']
        by_space[len(indices)].append(distance)
        assert indices==list(range(19,23)) if len(indices)==4 else indices==list(range(18,24))
        assert s['intent']['production_model_accepted'] is False
    assert all(sorted(ds)==[2.,2.75,3.5,6.,30.] for ds in by_space.values())
    for count in ('5_4','7_6'):
        with pytest.raises(ValueError,match='diagnostic-only'):
            screens._selection(source,case,PLAN/f'ch2_sh_s_cs_probe{count}.json',followup.read(case))
    assert all(p.read_bytes()==data for p,data in before.items())


@pytest.mark.parametrize('distances',[[],[True],[0],[float('nan')],[2.,2.]])
def test_invalid_relaxation_grid_is_refused_before_writes(tree,distances):
    production,review,destination,source,request,case=tree
    with pytest.raises(ValueError,match='distinct positive finite'):
        relaxation.specifications(source,request,PLAN/'ch2_sh_s_cs_probe5_4.json',case,
                                  distances_angstrom=distances)
    assert not list(destination.iterdir())


def test_changed_native_source_cannot_be_reused(tree):
    production,review,destination,source,request,case=tree
    followup.stage(production,review,destination)
    manifest=followup.read(source/'native/reconnaissance.json')
    job=manifest['jobs'][0]
    (Path(job['directory'])/(job['name']+'.out')).write_text('changed native evidence')
    with pytest.raises(ValueError,match='incomplete or changed'):
        followup.stage(production,review,destination)


def legacy_metadata_change(tree):
    """Reproduce panels prepared before queue-all rewrote source exclusions."""
    production,review,destination,source,*_=tree
    staged=followup.stage(production,review,destination)
    shutil.rmtree(destination/'relaxations/.native_evidence')
    from kinbot.execution import file_sha256
    old=file_sha256(destination/'candidates/source_manifest.json')
    path=source/'native/reconnaissance.json'
    assert old==file_sha256(path)
    for graph in staged['graphs']:
        proof=followup.read(Path(graph)/'workflow.json')['intent']['relaxation_recon']
        assert proof['evidence'][str(path)]==old
    m=followup.read(path)
    m.update(slurm_excluded_nodes=['node34'],slurm_inherited_exclusion_allow=['node16','node21','node24'])
    atomic_json(path,m)
    return staged


def test_legacy_manifest_bookkeeping_adopts_exact_graphs_and_native_records(tree):
    production,review,destination,source,*_=tree
    first=legacy_metadata_change(tree)
    retained={p:p.read_bytes() for g in first['graphs'] for p in Path(g).rglob('*') if p.is_file()}
    before={p:p.read_bytes() for p in destination.rglob('*') if p.is_file()}
    assert followup.stage(production,review,destination,validate_only=True)==first
    assert {p:p.read_bytes() for p in destination.rglob('*') if p.is_file()}==before
    for _ in range(2):assert followup.stage(production,review,destination)==first
    assert all(p.read_bytes()==data for p,data in retained.items())
    assert len(list((destination/'relaxations').glob('probe_*')))==10
    assert len(list((destination/'relaxations/.native_evidence').glob('*.json')))==2
    # Later resumes use the panel's own immutable snapshots, without depending
    # on a candidate directory or resetting the original graph bindings.
    path=source/'native/reconnaissance.json';m=followup.read(path)
    m['slurm_excluded_nodes'].append('another_bad_node');atomic_json(path,m)
    result=relaxation.prepare(production,PLAN/'ch2_sh_s_cs_relaxations.json',destination/'relaxations')
    assert result['graphs']==first['graphs']
    assert all(p.read_bytes()==data for p,data in retained.items())


def test_legacy_panel_requires_the_matching_historical_manifest(tree):
    production,review,destination,*_=tree
    legacy_metadata_change(tree)
    with pytest.raises(ValueError,match='Missing hash-matched'):
        relaxation.prepare(production,PLAN/'ch2_sh_s_cs_relaxations.json',
                           destination/'relaxations',validate_only=True)


@pytest.mark.parametrize('change',['snapshot','threshold','native','selection','request'])
def test_manifest_hash_exception_never_masks_scientific_changes(tree,change):
    production,review,destination,source,request,case=tree
    first=legacy_metadata_change(tree)
    if change=='snapshot':
        # Hash-bound legacy evidence cannot be replaced by a current manifest.
        snap=destination/'candidates/source_manifest.json'
        snap.write_bytes((source/'native/reconnaissance.json').read_bytes())
    elif change=='threshold':
        par=followup.read(case);par['vrc_tst_scan_fmax']=.001;atomic_json(case,par)
    elif change=='native':
        m=followup.read(source/'native/reconnaissance.json');j=m['jobs'][-1]
        out=Path(j['directory'])/(j['name']+'.out');out.write_text('changed output')
        from kinbot.execution import file_sha256
        j['output_sha256']=file_sha256(out);atomic_json(source/'native/reconnaissance.json',m)
    elif change=='request':
        data=followup.read(request);data['positions'][0][0]+=.001;atomic_json(request,data)
    else:
        g=Path(first['graphs'][0]);spec=followup.read(g/'workflow.json')
        spec['tasks'][0]['profile']['calculator_kwargs']['state_model']['orbital_selection']['active_indices']=[18,20,21,22]
        atomic_json(g/'workflow.json',spec)
        from kinbot.execution import file_sha256
        state=followup.read(g/'state.json');state['spec_sha256']=file_sha256(g/'workflow.json');atomic_json(g/'state.json',state)
    with pytest.raises(ValueError):followup.stage(production,review,destination,validate_only=True)


def test_recipe_change_during_wait_prevents_any_campaign_action(tree,tmp_path,monkeypatch):
    production,review,destination,*_=tree
    plan=tmp_path/'plan';plan.mkdir()
    for name in ('ch2_sh_s_cs_candidates.json','ch2_sh_s_cs_relaxations.json',
                 'ch2_sh_s_cs_probe5_4.json','ch2_sh_s_cs_probe7_6.json'):
        shutil.copyfile(PLAN/name,plan/name)
    real=fcntl.flock
    def lock(file,operation):
        if operation==fcntl.LOCK_EX and getattr(file,'name','')==str(review/'.handoff.lock'):
            path=plan/'ch2_sh_s_cs_candidates.json'
            path.write_text(path.read_text()+'\n')
        return real(file,operation)
    monkeypatch.setattr(followup.fcntl,'flock',lock)
    monkeypatch.setattr(screens,'run',lambda *a,**k:pytest.fail('Must stop before running campaigns'))
    with pytest.raises(ValueError,match='recipe changed while waiting'):
        followup.perform(production,review,destination,plan)
    assert not (destination/'candidates').exists()


def test_cli_uses_existing_verified_core_pin_without_controller_refresh(tree,monkeypatch):
    import kinbot
    from kinbot import rotdpy_worker_code
    production,review,destination,*_=tree
    pin=production/'.salvage_code.test';pin.mkdir()
    atomic_json(review/'refresh_complete.json',{'code':str(pin)})
    monkeypatch.setattr(sys,'argv',['continue_ch2_sh_s_review.py',str(production)])
    validated=[]
    monkeypatch.setattr(rotdpy_worker_code,'validate',lambda p:validated.append(p))
    class Reexec(Exception):pass
    def execute(p,args,env):
        assert p==sys.executable and args[:3]==[p,'-B','-u']
        assert args[-1]==str(production) and env['PYTHONPATH']==str(pin)
        assert env['OPENBLAS_NUM_THREADS']=='1'
        raise Reexec
    monkeypatch.setattr(followup.os,'execve',execute)
    with pytest.raises(Reexec):followup.main()
    monkeypatch.setattr(kinbot,'__file__',str(pin/'kinbot/__init__.py'))
    actions=[]
    monkeypatch.setattr(followup,'perform',lambda *args:actions.append(args))
    followup.main()
    assert actions==[(production,review,destination)]
    assert validated==[pin,pin]


@pytest.mark.parametrize('fail_candidates',[False,True])
def test_followup_waits_for_owner_drains_old_work_and_retains_repeat(tree,monkeypatch,fail_candidates):
    production,review,destination,source,request,case=tree
    events=[];started=threading.Event();attempts=[]
    def fixed(path,**options):
        events.append(('fixed',path))
        if path==destination/'candidates':
            for system in followup.read(path/'review_launch.json')['systems']:
                m=followup.read(system['manifest'])
                for job in m['jobs']:
                    if job.get('status')!='complete' and not job.get('job_id'):
                        attempts.append(job['name'])
                        job.update(status='failed' if fail_candidates else 'complete',job_id='saved-job')
                atomic_json(system['manifest'],m)
            return {'status':'failed' if fail_candidates else 'complete_review_required'}
        return {'status':'complete_review_required'}
    def relaxed(path):
        events.append(('relaxed',path))
        return {'panels':[{'task':{'status':'complete'}}]}
    monkeypatch.setattr(screens,'run',fixed)
    monkeypatch.setattr(screens,'evidence',lambda _: 'Native evidence placeholder for orchestration test')
    monkeypatch.setattr(relaxation,'run',relaxed)
    real_flock=fcntl.flock
    def lock(file,operation):
        if operation==fcntl.LOCK_EX and getattr(file,'name','')==str(review/'.handoff.lock'):
            started.set()
        return real_flock(file,operation)
    monkeypatch.setattr(followup.fcntl,'flock',lock)
    outcome=[]
    def work():
        try:outcome.append(followup.perform(production,review,destination))
        except Exception as e:outcome.append(e)
    with (review/'.handoff.lock').open('a+') as prior:
        real_flock(prior,fcntl.LOCK_EX|fcntl.LOCK_NB)
        thread=threading.Thread(target=work);thread.start()
        assert started.wait(5)
        assert events==[]
        with pytest.raises(RuntimeError,match='already has an owner'):
            followup.perform(production,review,destination)
        real_flock(prior,fcntl.LOCK_UN)
    thread.join(20)
    assert not thread.is_alive()
    expected_status = {'status':'candidate_review_required','reason':
        'Incomplete or failed fixed-geometry comparison; relaxed probes remain staged.',
        'production_model_accepted':False} if fail_candidates else {
            'status':'review_required','production_model_accepted':False}
    assert outcome==[expected_status]
    expected=[('fixed',review/'screens'),('relaxed',review/'relaxations'),('fixed',destination/'candidates')]
    if not fail_candidates:expected.append(('relaxed',destination/'relaxations'))
    assert events==expected and len(attempts)==24
    again=followup.perform(production,review,destination)
    assert again==outcome[0] and len(attempts)==24
    assert events==expected*2
