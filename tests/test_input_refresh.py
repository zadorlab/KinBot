from copy import deepcopy
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import input_refresh as refresh, pes_phase as phase, live_control as control
from kinbot.execution import atomic_json, file_sha256


def surface(bond=(1,2)):
    return dict(schema=1,reviewed=True,review_note='Reviewed fixture orbitals',channels=[dict(
        bond=list(bond),combination='single',bottlenecks=[dict(name='single',coordinate='reactive',
            distances_angstrom=[3.,4.,5.],states=[dict(multiplicity=2,root=1,weight=1.,
                energy_offset_cm_inverse=0.,energy_reference_note='Ground asymptote')],
            parameters={'vrc_tst_scan_method':'caspt2-f12','vrc_tst_scan_basis':'avdz-f12'})])])


@pytest.fixture
def run(tmp_path):
    control.configure(tmp_path)
    root=tmp_path/'case';root.mkdir()
    inp=root/'input.json'
    original=dict(vdW_detection=.5,vrc_tst_surface_models={'123':surface()})
    atomic_json(inp,original)
    phase.initialize(root,inp)
    for name in ('123','456'):
        w=root/name;w.mkdir();atomic_json(w/(name+'.json'),original)
        phase.register(root,name,os.getpid())
    m=json.loads((root/phase.NAME).read_text())
    m['phase']='postprocessing'
    for w in m['workers'].values(): w['status']='postprocessing'
    atomic_json(root/phase.NAME,m)
    atomic_json(tmp_path/'prepared.json',dict(inputs=[str(inp)],sha256={'case/input.json':file_sha256(inp)}))
    return root,inp,original


def change(inp,original):
    updated=deepcopy(original)
    updated['vrc_tst_surface_models']['456']=surface((2,3))
    atomic_json(inp,updated)
    return updated


def test_external_input_retains_origin_and_resubmits_from_run_directory(tmp_path,monkeypatch):
    root=tmp_path/'run';root.mkdir()
    inputs=tmp_path/'inputs';inputs.mkdir()
    inp=inputs/'external.json'
    original=dict(vdW_detection=.5,vrc_tst_surface_models={'123':surface()})
    atomic_json(inp,original)
    control.configure(root)
    phase.initialize(root,inp)
    worker=root/'123';worker.mkdir();atomic_json(worker/'123.json',original)
    phase.register(root,'123',os.getpid())
    m=json.loads((root/phase.NAME).read_text())
    m['phase']='postprocessing';m['workers']['123']['status']='postprocessing'
    atomic_json(root/phase.NAME,m)
    monkeypatch.chdir(root)
    assert refresh.owner_directory(inp)==root
    accepted=refresh.accepted_source(inp,directory=root)
    assert accepted==root/refresh.PARAMETERS_SNAPSHOT
    assert json.loads(accepted.read_text())==dict(original,input_file=str(inp))
    # This is a proposal until the owner consumes it; relative orbital paths
    # continue to resolve beside external.json, never beside the run snapshot.
    updated=change(inp,original)
    assert json.loads(accepted.read_text())['vdW_detection']==.5
    assert refresh.submit(inp,apply=True)['status']=='queued'
    phase.tick(root,['123'],max_searches=1)
    assert json.loads(accepted.read_text())==dict(updated,input_file=str(inp))
    assert json.loads((root/refresh.SNAPSHOT).read_text())==updated
    monkeypatch.chdir(tmp_path)
    assert refresh.owner_directory(inp,directory=root)==root
    assert refresh.submit(inp,apply=True,directory=root)['id']==refresh.receipt(root)['id']


def test_preview_and_publication_never_launch_or_touch_qc(run,monkeypatch):
    root,inp,original=run
    native=root/'123/native.out';native.write_bytes(b'expensive result')
    updated=change(inp,original)
    monkeypatch.setattr(__import__('subprocess'),'Popen',lambda *a,**k: pytest.fail('Launched duplicate driver'))
    # Identity is tested separately; avoid its portable ps subprocess here.
    monkeypatch.setattr(phase,'owner_alive',lambda m:True)
    monkeypatch.setattr(phase,'active_worker',lambda p:True)
    assert refresh.submit(inp)['status']=='reviewed'
    assert not (root/refresh.JOURNAL).exists()
    assert json.loads(refresh.accepted_source(inp).read_text())==original
    request=refresh.submit(inp,apply=True)
    assert refresh.submit(inp,apply=True)['id']==request['id']
    assert json.loads((root/refresh.SNAPSHOT).read_text())==original
    assert phase.tick(root,['123','456'],max_searches=1)[1]=='exploring'
    assert json.loads((root/refresh.SNAPSHOT).read_text())==updated
    assert len(control.requests(root,kinds={'network'}))==2
    assert json.loads((root.parent/'prepared.json').read_text())['sha256']['case/input.json']==file_sha256(inp)
    assert native.read_bytes()==b'expensive result'
    assert refresh.receipt(root)['status']=='published'


def test_owner_applies_models_without_reconstructing_reactions(run,monkeypatch):
    from kinbot.search_control import revise_network
    root,inp,original=run
    change(inp,original);refresh.submit(inp,apply=True)
    phase.tick(root,['123','456'],max_searches=1)
    for name in ('123','456'):
        w=root/name;monkeypatch.chdir(w)
        # No reaction constructor is available: model-only edits must retain it.
        reaction=object()
        species=SimpleNamespace(reac_obj=[reaction],reac_ts_done=[-1],reac_step=[5],
                                name=name,chemid=int(name),smiles='fixture',
                                optical_reference={'status':'assigned','canonical_graphs':[]})
        generator=SimpleNamespace(inp=str(w/(name+'.json')),par=deepcopy(original),species=species)
        revise_network(generator)
        assert species.reac_obj==[reaction] and species.reac_ts_done==[-1]
        assert '456' in generator.par['vrc_tst_surface_models']
    assert refresh.receipt(root)['status']=='applied'
    # A second distinct proposal can follow a fully acknowledged revision.
    next_input=json.loads(inp.read_text());next_input['vrc_tst_surface_models']['789']=surface()
    atomic_json(inp,next_input)
    m=json.loads((root/phase.NAME).read_text());m['phase']='postprocessing';atomic_json(root/phase.NAME,m)
    assert refresh.submit(inp,apply=True)['status']=='queued'


@pytest.mark.parametrize('field,value', [('method','different'),('mult',4),('vdW_detection',-1.)])
def test_unsupported_input_changes_are_refused_without_publication(run,field,value):
    root,inp,original=run
    atomic_json(inp,dict(original,**{field:value}))
    with pytest.raises(ValueError): refresh.submit(inp,apply=True)
    assert not (root/refresh.JOURNAL).exists()
    assert not control.requests(root)


def test_existing_surface_cannot_be_replaced(run):
    root,inp,original=run
    changed=deepcopy(original)
    changed['vrc_tst_surface_models']['123']['channels'][0]['bottlenecks'][0]['parameters']['vrc_tst_scan_basis']='avtz-f12'
    atomic_json(inp,changed)
    with pytest.raises(ValueError,match='scientific migration'): refresh.submit(inp,apply=True)
    assert not control.requests(root)


def test_reedit_before_consumption_is_blocked(run):
    root,inp,original=run
    change(inp,original);refresh.submit(inp,apply=True)
    atomic_json(inp,dict(original,vdW_detection=0.))
    phase.tick(root,['123','456'],max_searches=1)
    assert refresh.receipt(root)['status']=='blocked'
    assert not control.requests(root)
    assert json.loads((root/refresh.SNAPSHOT).read_text())==original


def test_interrupted_batch_publication_resumes_without_duplicate_requests(run,monkeypatch):
    root,inp,original=run
    change(inp,original);refresh.submit(inp,apply=True)
    real=refresh.atomic_json
    def fail(path,value):
        if Path(path).name==refresh.SNAPSHOT: raise OSError('Injected interruption')
        return real(path,value)
    monkeypatch.setattr(refresh,'atomic_json',fail)
    with pytest.raises(OSError): phase.tick(root,['123','456'],max_searches=1)
    assert not control.requests(root)  # Incomplete batch is not visible to workers.
    ids=sorted(p.name for p in (root.parent/'.kinbot_controls/requests').glob('*.json'))
    monkeypatch.setattr(refresh,'atomic_json',real)
    phase.tick(root,['123','456'],max_searches=1)
    assert sorted(r['id']+'.json' for _,r in control.requests(root))==ids
    assert refresh.receipt(root)['status']=='published'


def test_stale_owner_and_exploration_science_edits_are_not_live_submission_targets(run):
    root,inp,original=run
    atomic_json(inp,dict(original,vdW_detection=0.))
    m=json.loads((root/phase.NAME).read_text());m['phase']='exploring';atomic_json(root/phase.NAME,m)
    with pytest.raises(ValueError,match='exploration'): refresh.submit(inp,apply=True)
    m['phase']='postprocessing';m['owner_identity']['start']='obsolete';atomic_json(root/phase.NAME,m)
    with pytest.raises(ValueError,match='No live PES'): refresh.submit(inp,apply=True)
    assert not (root/refresh.JOURNAL).exists()


def test_orbital_review_during_exploration_inherits_at_safe_boundary(run,monkeypatch):
    root, inp, original = run
    m=json.loads((root/phase.NAME).read_text())
    m['phase']='exploring'
    m['workers']['123']['status']='vrc'
    m['workers']['456']['status']='exploring'
    atomic_json(root/phase.NAME,m)
    change(inp, original)
    assert refresh.submit(inp,apply=True)['workers']==1
    phase.tick(root,['123','456'],max_searches=2)
    assert [Path(r['target']).name for _,r in control.requests(root)]==['123']
    w=root/'456'
    generator=SimpleNamespace(inp=str(w/'456.json'),par=deepcopy(original))
    refresh.inherit_models(root,generator)
    assert '456' in generator.par['vrc_tst_surface_models']
    assert json.loads((w/'456.json').read_text())['vrc_tst_surface_models']==generator.par['vrc_tst_surface_models']
    # A branch with a pending generation-fenced request must not preempt it.
    w=root/'123'
    generator=SimpleNamespace(inp=str(w/'123.json'),par=deepcopy(original))
    before=file_sha256(w/'123.json')
    refresh.inherit_models(root,generator)
    assert file_sha256(w/'123.json')==before
    # A newly discovered well whose input was written earlier inherits too.
    w=root/'789';w.mkdir();atomic_json(w/'789.json',original)
    generator=SimpleNamespace(inp=str(w/'789.json'),par=deepcopy(original))
    refresh.inherit_models(root,generator)
    assert '456' in generator.par['vrc_tst_surface_models']


def test_review_includes_worker_that_finishes_between_submit_and_consume(run):
    root, inp, original=run
    m=json.loads((root/phase.NAME).read_text());m['phase']='exploring'
    for w in m['workers'].values():w['status']='exploring'
    atomic_json(root/phase.NAME,m)
    change(inp,original)
    assert refresh.submit(inp,apply=True)['workers']==0
    m['workers']['456']['status']='vrc';atomic_json(root/phase.NAME,m)
    phase.tick(root,['123','456'],max_searches=2)
    assert [Path(r['target']).name for _,r in control.requests(root)]==['456']
    assert refresh.receipt(root)['status']=='published'


def test_review_conflict_before_publication_preserves_accepted_snapshot(run):
    root,inp,original=run
    m=json.loads((root/phase.NAME).read_text());m['phase']='exploring'
    for w in m['workers'].values():w['status']='exploring'
    atomic_json(root/phase.NAME,m)
    change(inp,original);refresh.submit(inp,apply=True)
    m['workers']['456']['status']='vrc';atomic_json(root/phase.NAME,m)
    other=control.request(root.parent,root/'456',kind='network',reason='Concurrent review',
                          options={'vdW_detection':0.,'rejected':False})
    phase.tick(root,['123','456'],max_searches=2)
    assert refresh.receipt(root)['status']=='blocked'
    assert json.loads((root/refresh.SNAPSHOT).read_text())==original
    assert [r['id'] for _,r in control.requests(root)]==[other['id']]


def test_repeated_pes_command_targets_owner_not_main(run,monkeypatch,capsys):
    from kinbot import pes
    root,inp,original=run
    change(inp,original)
    monkeypatch.setattr(pes.sys,'argv',['kinbot.pes',str(inp)])
    monkeypatch.setattr(pes,'_main',lambda:pytest.fail('Started second driver'))
    pes.main()
    assert json.loads(capsys.readouterr().out)['drivers_started']==0


def test_mixed_reviewed_and_new_vrc_channels_continue_independently(tmp_path,monkeypatch):
    import numpy as np
    from kinbot.vrc_tst_scan import VTS
    from kinbot import vrc_multistate as multistate, vrc_orbitals, pes
    monkeypatch.chdir(tmp_path)
    well=SimpleNamespace(name='123',chemid=123,smiles='C[S]',atom=['C','S'],
                         geom=np.zeros((2,3)),charge=0,mult=2)
    par=dict(vrc_tst_scan={'123':['reviewed','new']},vrc_tst_noscan={},pes=1)
    vts=VTS(well,par,None)
    vts.scan_reac={n:SimpleNamespace(scan_coo=[0,1],products=[]) for n in ('reviewed','new')}
    monkeypatch.setattr(vrc_orbitals,'asymptotic_fragments',lambda p:False)
    monkeypatch.setattr(vrc_orbitals,'pivot_backend',lambda p:'none')
    monkeypatch.setattr(vts,'configured_reactions',lambda p:['reviewed','new'] if p else [])
    calls=[]
    monkeypatch.setattr(vts,'opt_products',lambda names,submit=True:calls.append(('products',names,submit)))
    monkeypatch.setattr(vts,'save_products',lambda names:calls.append(('save',names)))
    monkeypatch.setattr(vts,'find_scan_coos',lambda names:None)
    monkeypatch.setattr(vts,'find_equiv',lambda names:None)
    def selection(par,reaction):
        if reaction is vts.scan_reac['new']: raise multistate.SurfaceModelRequired('new needs review')
        return ('reviewed',{})
    monkeypatch.setattr(multistate,'channel_for_reaction',selection)
    monkeypatch.setattr(multistate,'calculate_channels',lambda v,s:calls.extend(('calculated',n) for n in s))
    monkeypatch.setattr(pes,'create_rotdpy_inputs',lambda p,rows,*a,**k:calls.append(('rotd',[r[1] for r in rows])))
    with pytest.raises(multistate.SurfaceModelRequired,match='new'):
        vts.calculate_correction_potentials()
    assert ('calculated','reviewed') in calls
    assert ('products',['reviewed'],True) in calls
    assert ('save',['reviewed']) in calls
    assert ('rotd',['reviewed']) in calls
    assert json.loads((tmp_path/'vrctst/orbital_review/new/request.json').read_text())['status']=='needs_review'


def test_review_wait_releases_on_reopen_and_restores_qc_limits(tmp_path,monkeypatch):
    from kinbot import kb
    from kinbot.vrc_multistate import SurfaceModelRequired
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(kb.VTS,'calculate_correction_potentials',lambda s:(_ for _ in ()).throw(SurfaceModelRequired('review')))
    monkeypatch.setattr(control,'find',lambda *a:Path('control'))
    monkeypatch.setattr(phase,'active_worker',lambda *a:True)
    checks=[]
    def checkpoint():
        checks.append(1)
        if len(checks)>1: raise phase.ExplorationReopened()
    monkeypatch.setattr(phase,'checkpoint',checkpoint)
    monkeypatch.setattr(kb.time,'sleep',lambda *a:None)
    qc=SimpleNamespace(queue_job_limit=10)
    with pytest.raises(phase.ExplorationReopened):
        kb._execute_vrc_branch(SimpleNamespace(reac_obj=[]),{'vrc_tst_branch_slots':3},qc)
    assert qc.queue_job_limit==10
    assert json.loads(Path('vrc_model_review.json').read_text())['status']=='needs_review'


def test_prepared_proposal_checks_pending_before_writing_input(run):
    root,inp,original=run
    proposal=deepcopy(original)
    proposal['vrc_tst_surface_models']['456']=surface((2,3))
    assert refresh.submit(inp,proposal=proposal,expected=original)['status']=='reviewed'
    assert json.loads(inp.read_text())==original
    first=refresh.submit(inp,apply=True,proposal=proposal,expected=original)
    assert json.loads(inp.read_text())==proposal
    before=inp.read_bytes(); journal=(root/refresh.JOURNAL).read_bytes()
    second=deepcopy(proposal);second['vrc_tst_surface_models']['789']=surface()
    with pytest.raises(ValueError,match='still pending'):
        refresh.submit(inp,apply=True,proposal=second,expected=proposal)
    assert inp.read_bytes()==before and (root/refresh.JOURNAL).read_bytes()==journal
    assert refresh.submit(inp,apply=True,proposal=proposal,expected=proposal)['id']==first['id']
    phase.tick(root,['123','456'],max_searches=1)
    assert json.loads((root/refresh.SNAPSHOT).read_text())==proposal


def test_prepared_proposal_rejects_stale_editor_and_bad_science_without_writes(run):
    root,inp,original=run
    changed=change(inp,original)
    before=inp.read_bytes()
    with pytest.raises(ValueError,match='changed during'):
        refresh.submit(inp,apply=True,proposal=original,expected=original)
    bad=dict(changed,method='different')
    with pytest.raises(ValueError,match='Unsupported'):
        refresh.submit(inp,apply=True,proposal=bad,expected=changed)
    assert inp.read_bytes()==before and not (root/refresh.JOURNAL).exists()
