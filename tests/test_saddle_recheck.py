"""Explicit-theory checks through native owners, staging and controller restart."""
import ast
from copy import deepcopy
import json
from pathlib import Path

from ase import Atoms
import pytest

from kinbot import saddle_recheck as checks
from kinbot.execution import atomic_json
from kinbot.qc import QuantumChemistry
from kinbot.theory import THEORY_PRESETS


def manifest():
    profile=deepcopy(THEORY_PRESETS['b3lyp-b2plyp-anl']['l2'])
    return dict(schema=1,parameters=dict(l1_profile=profile,l2_profile=profile,
        method='B2PLYP',basis='cc-pVTZ',scan_method='B2PLYP',scan_basis='cc-pVTZ',
        barrier_threshold=150,queuing='local',sella_irc_steps=300,sella_irc_dx=.05,high_level=0),
        cases=[dict(name='seed',mode='optimize',seed=dict(symbols=['O','H','H'],
            positions=[[0.,0,0],[1,0,0],[0,1,0]],charge=0,multiplicity=1))])


def stage(tmp_path,data=None):
    data=data or manifest()
    source=tmp_path/'source.json'
    atomic_json(source,data)
    root=tmp_path/'run'
    checks.stage(source,root,tmp_path)
    return root,source


def test_stage_is_idempotent_and_does_not_rebase_inputs(tmp_path):
    root,source=stage(tmp_path)
    sentinel=root/'seed'/'accepted.log'
    sentinel.write_text('native calculation')
    before={p:p.read_bytes() for p in root.rglob('*') if p.is_file()}
    checks.stage(source,root,tmp_path)
    assert all(p.read_bytes()==data for p,data in before.items())
    changed=manifest();changed['cases'][0]['seed']['positions'][1][0]=1.1
    atomic_json(source,changed)
    with pytest.raises(ValueError,match='different scientific'):
        checks.stage(source,root,tmp_path)
    assert all(p.read_bytes()==data for p,data in before.items())


def test_native_source_hash_failure_prevents_copy(tmp_path):
    data=manifest();case=data['cases'][0]
    case.update(mode='reuse',native={'files':{'.log':dict(path='native.log',sha256='bad')}})
    (tmp_path/'native.log').write_text('source')
    with pytest.raises(ValueError,match='source hash differs'):
        stage(tmp_path,data)
    assert not (tmp_path/'run/seed/source.log').exists()


def test_complete_controller_restart_adopts_live_qc_and_all_irc_directions(tmp_path,monkeypatch):
    root,_=stage(tmp_path)
    allocations={};submitted=[];polls={}

    def submit(qc,job,nproc,**kwargs):
        assert job not in allocations
        contents=Path(job+'.py').read_text()
        ast.parse(contents)
        assert 'B2PLYP' in contents and 'cc-pVTZ' in contents and 'B3LYP' not in contents
        if '_IRC_' in job and '_freq_recovery_' not in job:
            assert 'steps=300' in contents and 'dx=0.05' in contents
        allocations[job]=str(len(allocations)+1)
        qc.job_ids[job]=allocations[job]
        submitted.append(job)
        return 1

    def status(qc,job,**kwargs):
        if job not in allocations:return 0
        polls[job]=polls.get(job,0)+1
        if polls[job]<3:return 'running'
        if not qc.db.count(name=job):
            atoms=Atoms('OHH',positions=[[0,0,0],[1,0,0],[0,1,0]])
            frequencies=[100.,200.,300.] if '_freq_recovery_' in job else [-100.,100.,200.]
            qc.db.write(atoms,name=job,data=dict(status='normal',energy=-1.,zpe=.01,
                                               frequencies=frequencies))
            if '_IRC_' in job and '_freq_recovery_' not in job:
                qc.db.write(atoms,name=job+'_prod',data=dict(status='normal',energy=-1.1))
                atomic_json(job+'.irc_path.json',dict(converged=True,moved=True))
        return 'normal'

    monkeypatch.setattr(QuantumChemistry,'submit_qc',submit)
    monkeypatch.setattr(QuantumChemistry,'check_qc',status)
    class Interrupted(Exception):pass
    monkeypatch.setattr(checks.time,'sleep',lambda _:(_ for _ in ()).throw(Interrupted()))
    with pytest.raises(Interrupted):checks.run(root)
    worker=root/'seed'/'saddle.py';before=worker.read_bytes()
    assert submitted==['saddle']
    # Recreate all controller objects from saved routing and case records.
    monkeypatch.setattr(checks.time,'sleep',lambda _:None)
    checks.run(root)
    assert submitted[:3]==['saddle','saddle_IRC_F','saddle_IRC_R']
    assert len(submitted)==5 and all('_freq_recovery_' in name for name in submitted[3:])
    assert worker.read_bytes()==before
    results=json.loads((root/'results.json').read_text())
    assert results['seed']['outcome']=='irc_endpoints_available'
    assert len(results['seed']['endpoint_frequencies'])==2
    checks.run(root)
    assert len(submitted)==5
    assert json.loads((root/'results.json').read_text())==results


def test_terminal_search_outcome_is_not_automatically_retried(tmp_path,monkeypatch):
    root,_=stage(tmp_path)
    monkeypatch.chdir(root/'seed')
    par=checks.Parameters('input.json').par
    qc=QuantumChemistry(par)
    monkeypatch.setattr(qc,'check_qc',lambda job:'error')
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(qc,'qc_opt_ts',lambda *a,**kw:pytest.fail('Unexpected retry'))
        first=checks.tick(manifest()['cases'][0],qc,par)
        assert first['outcome']=='saddle_not_validated'
        assert checks.tick(manifest()['cases'][0],qc,par)==first


def test_changed_staged_input_stops_before_submission(tmp_path,monkeypatch):
    root,_=stage(tmp_path)
    p=root/'seed/input.json';data=json.loads(p.read_text());data['method']='other'
    atomic_json(p,data)
    monkeypatch.setattr(QuantumChemistry,'submit_qc',lambda *a,**kw:pytest.fail('Unexpected submission'))
    with pytest.raises(ValueError,match='input changed'):
        checks.run(root)


def test_family_restart_uses_durable_submission_checkpoint(tmp_path,monkeypatch):
    from kinbot.stationary_pt import StationaryPoint
    species=StationaryPoint('parent',0,2,smiles='[CH2]S')
    species.characterize()
    carbon=list(species.atom).index('C');sulfur=list(species.atom).index('S')
    hydrogen=next(i for i,s in enumerate(species.atom) if s=='H' and species.bond[sulfur][i])
    data=manifest();case=data['cases'][0]
    case.update(mode='family',source_reaction='parent_R_Addition_MultipleBond_1_2_3',
                motif=[carbon,sulfur,hydrogen],guess=dict(bond_scale=1.6,angle='parent'),
                seed=dict(symbols=list(species.atom),positions=species.geom.tolist(),charge=0,multiplicity=2))
    root,_=stage(tmp_path,data);monkeypatch.chdir(root/'seed')
    parameters=checks.Parameters('input.json').par
    qc=QuantumChemistry(parameters)
    monkeypatch.setattr(QuantumChemistry,'check_qc',lambda *a:'running')
    class Interrupted(Exception):pass
    def interrupt(reaction,step,command):
        checkpoint=reaction.search_checkpoint
        checkpoint.prepare_submission(reaction,step)
        checkpoint.submitted(reaction.instance_name,'41')
        qc.l1.job_ids[reaction.instance_name]='41'
        qc._record(reaction.instance_name,'l1','41')
        raise Interrupted
    monkeypatch.setattr(checks,'carry_out_reaction',interrupt)
    with pytest.raises(Interrupted):checks.tick(case,qc,parameters)
    resumed=QuantumChemistry(parameters)
    observed=[]
    monkeypatch.setattr(checks,'carry_out_reaction',lambda reaction,step,command:observed.append(step) or step)
    checks.tick(case,resumed,parameters)
    assert observed==[1]
    saved=json.loads(Path('.kinbot_search_state.json').read_text())
    entry=next(iter(saved['reactions'].values()))
    assert entry['step']==1 and 'pending_submission' not in entry


def test_second_controller_preserves_live_owner(tmp_path,monkeypatch):
    import fcntl
    root,_=stage(tmp_path)
    atomic_json(root/'controller.json',dict(status='running',pid=12345))
    before=(root/'controller.json').read_bytes()
    monkeypatch.setattr(QuantumChemistry,'submit_qc',lambda *a,**kw:pytest.fail('Duplicate submission'))
    with (root/'controller.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert checks.run(root)==0
    assert (root/'controller.json').read_bytes()==before


def test_selected_path_run_adopts_forward_and_preserves_other_cases(tmp_path,monkeypatch):
    data=manifest()
    other=deepcopy(data['cases'][0]);other['name']='other'
    data['cases'].append(other)
    root,_=stage(tmp_path,data)
    atomic_json(root/'results.json',{'other':{'stage':'irc','irc_statuses':{'F':'running'}}})
    atomic_json(root/'controller.json',{'status':'running','pid':123})
    atoms=Atoms('OHH',positions=[[0,0,0],[1,0,0],[0,1,0]])
    with checks.in_directory(root/'seed'):
        db=QuantumChemistry(checks.Parameters('input.json').par).db
    db.write(atoms,name='saddle',data=dict(status='normal',energy=-1.,frequencies=[-100.,100.,200.]))
    db.write(atoms,name='saddle_IRC_F',data=dict(status='normal',energy=-1.1))
    db.write(atoms,name='saddle_IRC_F_prod',data=dict(status='normal',energy=-1.1))
    atomic_json(root/'seed/saddle_IRC_F.irc_path.json',dict(converged=True,moved=True))
    before={p:p.read_bytes() for p in (root/'other').rglob('*') if p.is_file()}
    full_owner=(root/'controller.json').read_bytes()
    allocations={};submitted=[];polls={}
    def submit(qc,job,nproc,**kwargs):
        assert Path.cwd()==root/'seed'
        assert job=='saddle_IRC_R' and job not in allocations
        allocations[job]='42';submitted.append(job);qc.job_ids[job]='42'
        return 1
    def status(qc,job,**kwargs):
        if qc.db.count(name=job):return 'normal'
        if job not in allocations:return 0
        polls[job]=polls.get(job,0)+1
        if polls[job]<3:return 'running'
        qc.db.write(atoms,name=job,data=dict(status='normal',energy=-1.1))
        qc.db.write(atoms,name=job+'_prod',data=dict(status='normal',energy=-1.1))
        atomic_json(job+'.irc_path.json',dict(converged=True,moved=True))
        return 'normal'
    monkeypatch.setattr(QuantumChemistry,'submit_qc',submit)
    monkeypatch.setattr(QuantumChemistry,'check_qc',status)
    monkeypatch.setattr(QuantumChemistry,'qc_freq',lambda *a,**kw:pytest.fail('Unrequested frequency'))
    class Interrupted(Exception):pass
    monkeypatch.setattr(checks.time,'sleep',lambda _:(_ for _ in ()).throw(Interrupted()))
    with pytest.raises(Interrupted):
        checks.main(['run',str(root),'--only','seed','--paths-only'])
    assert submitted==['saddle_IRC_R']
    assert json.loads((root/'results.json').read_text())['seed']['irc_statuses']['R']=='running'
    worker=(root/'seed/saddle_IRC_R.py').read_bytes()
    monkeypatch.setattr(checks.time,'sleep',lambda _:None)
    checks.main(['run',str(root),'--only','seed','--paths-only'])
    assert submitted==['saddle_IRC_R']
    assert (root/'seed/saddle_IRC_R.py').read_bytes()==worker
    result=json.loads((root/'results.json').read_text())
    assert result['other']=={'stage':'irc','irc_statuses':{'F':'running'}}
    assert result['seed']['outcome']=='irc_paths_available'
    assert result['seed']['complete'] is False
    assert (root/'controller.json').read_bytes()==full_owner
    assert all(p.read_bytes()==content for p,content in before.items())
    assert json.loads((root/'controller.selected.json').read_text())['status']=='complete'
    assert json.loads((root/'selected_execution.json').read_text())==dict(schema=1,
        cases=['seed'],paths_only=True,campaign_complete=False,scientific_approval=False)
    checks.run(root,only=['seed'],paths_only=True)
    assert submitted==['saddle_IRC_R']


@pytest.mark.parametrize('only',[[],['missing'],['seed','missing']])
def test_unknown_selection_fails_before_writing_or_submitting(tmp_path,monkeypatch,only):
    root,_=stage(tmp_path)
    before={p:p.read_bytes() for p in root.rglob('*') if p.is_file()}
    monkeypatch.setattr(checks,'QuantumChemistry',lambda *a:pytest.fail('QC initialized'))
    with pytest.raises(ValueError,match='selection'):
        checks.run(root,only=only,paths_only=True)
    assert {p:p.read_bytes() for p in root.rglob('*') if p.is_file()}==before


def test_selected_run_respects_existing_owner_lock(tmp_path,monkeypatch):
    import fcntl
    root,_=stage(tmp_path)
    monkeypatch.setattr(checks,'QuantumChemistry',lambda *a:pytest.fail('QC initialized'))
    with (root/'controller.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        assert checks.run(root,only=['seed'],paths_only=True)==0
    assert not (root/'selected_execution.json').exists()


def test_paths_only_receipt_can_later_continue_frequency_validation(tmp_path,monkeypatch):
    root,_=stage(tmp_path)
    monkeypatch.chdir(root/'seed')
    atomic_json('recheck.json',dict(schema=1,complete=False,outcome='irc_paths_available'))
    par=checks.Parameters('input.json').par
    qc=QuantumChemistry(par)
    seen=[]
    monkeypatch.setattr(qc,'check_qc',lambda job:'normal')
    monkeypatch.setattr(checks.IRC,'do_irc_calculations',lambda self:None)
    monkeypatch.setattr(qc,'qc_freq',lambda species,job,**kw:seen.append(job) or job+'_freq')
    def observation(qc,job,seed):
        return dict(job=job,symbols=['O','H','H'],positions=[[0,0,0],[1,0,0],[0,1,0]],
            data=dict(status='normal',frequencies=[-100.,100.,200.] if job=='saddle' else [100.,200.,300.]))
    monkeypatch.setattr(checks,'observation',observation)
    for direction in ('F','R'):
        atomic_json(f'saddle_IRC_{direction}.irc_path.json',dict(converged=True,moved=True))
    result=checks.tick(manifest()['cases'][0],qc,par)
    assert seen==['saddle_IRC_F_prod','saddle_IRC_R_prod']
    assert result['complete'] and result['outcome']=='irc_endpoints_available'
