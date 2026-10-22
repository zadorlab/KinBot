"""Exercise the exact new-well screen recipe without submitting native work."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from kinbot.execution import atomic_json

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('requested_orbital_recon',
    ROOT/'examples/anl/tools/requested_orbital_recon.py')
helper = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(helper)
RECIPE = ROOT/'examples/anl/ch3s_anl1_f12/orbital_reviews/ch2sh_cs_screen.json'
CANDIDATES = RECIPE.with_name('ch2sh_cs_candidates.json')
SELECTION = RECIPE.with_name('ch2sh_cs_selected.json')


@pytest.fixture
def request_tree(tmp_path):
    pool = tmp_path/'production'
    case = pool/'case'
    directory = case/'471030810040000000002/vrctst/orbital_review/reaction'
    directory.mkdir(parents=True)
    request = directory/'request.json'
    recipe = json.loads(RECIPE.read_text())
    data = dict(recipe['request_match'], schema=1, status='needs_review',
                smiles='[H][S]=[C]([H])[H]', positions=[
        [1.03440581452718, -.3719293385761376, .03988588054708078],
        [.3381038973372055, .7267766509223542, 1.1679876585701665],
        [.7607169349022845, 1.8652195195547316, .6015953344491213],
        [.8725603978962576, -1.4191794072599508, .22000490140932294],
        [1.6514884746631548, -.030599788522718163, -.7700430150262825]])
    atomic_json(request, data)
    atomic_json(case/'input.json', {'queue_name': 'test-partition'})
    atomic_json(pool/'.kinbot_resource_policy.json', dict(schema=1, max_jobs=23,
        excluded_nodes=['badnode'], manage_pending=True, orbital_review_reserve=3))
    return request, case/'input.json', tmp_path/'screen'


def test_exact_request_prepares_eight_independent_probes_and_no_production_writes(request_tree, monkeypatch):
    request, case_input, output = request_tree
    original = {p: p.read_bytes() for p in case_input.parent.parent.rglob('*') if p.is_file()}
    monkeypatch.setattr(helper.subprocess, 'run', lambda *a, **k: pytest.fail('No launch during preparation'))
    result = helper.prepare(request, case_input, RECIPE, output)
    assert result['tasks'] == 8 and result['max_jobs'] == 3
    assert result['active_space']['active_electrons'] == 13
    assert result['active_space']['active_orbitals'] == 11
    manifest = json.loads((output/'native/reconnaissance.json').read_text())
    assert manifest['slurm_excluded_nodes'] == ['badnode']
    for job in manifest['jobs']:
        assert job['orbital_predecessor'] is None
        assert [(s['multiplicity'], s['root']) for s in job['root_states']] == [
            (4, 1), (4, 2), (2, 1), (2, 2), (2, 3), (2, 4)]
        d = Path(job['directory'])
        deck = (d/(job['name']+'.inp')).read_text()
        assert '{rhf,so-sci,maxit=256;wf,25,1,3;occ,14;closed,11;' in deck
        assert '{multi,cas-ci;' in deck and 'avas,nela=13' in deck
        assert 'center,5,1s' in deck and 'start,2100.2' in deck
        assert 'wf,25,1,1\nstate,4' in deck
        assert 'rs2' not in deck.lower()
        script = (d/(job['name']+'.sbatch')).read_text()
        assert '#SBATCH --exclusive' in script and '#SBATCH -p test-partition' in script
        assert '# KINBOT_WALLTIME=partition_max' in script
    assert all(p.read_bytes() == value for p, value in original.items())
    with pytest.raises(FileExistsError):
        helper.prepare(request, case_input, RECIPE, output)


@pytest.mark.parametrize('field,value', [('parent', 'wrong'), ('multiplicity', 4),
                                        ('status', 'model_accepted')])
def test_wrong_or_accepted_request_cannot_be_prepared(request_tree, field, value):
    request, case_input, output = request_tree
    data = json.loads(request.read_text())
    data[field] = value
    atomic_json(request, data)
    with pytest.raises(ValueError):
        helper.prepare(request, case_input, RECIPE, output)
    assert not output.exists()


def test_refuse_to_write_into_production_tree(request_tree):
    request, case_input, _ = request_tree
    with pytest.raises(ValueError, match='outside the production'):
        helper.prepare(request, case_input, RECIPE, case_input.parent/'screen')


def test_launch_uses_current_exclusions_reserved_cap_and_original_manifest(request_tree, monkeypatch):
    request, case_input, output = request_tree
    helper.prepare(request, case_input, RECIPE, output)
    policy_file = case_input.parent.parent/'.kinbot_resource_policy.json'
    policy = json.loads(policy_file.read_text())
    policy.update(orbital_review_reserve=2, excluded_nodes=['badnode', 'newbadnode'])
    atomic_json(policy_file, policy)
    monkeypatch.setattr(helper.subprocess, 'run', lambda *a, **k: None)
    from kinbot import run_walltime
    monkeypatch.setattr(run_walltime, 'partition_walltime', lambda _: '31-00:00:00')
    calls = []
    def launch(paths, **kwargs):
        calls.append((paths, kwargs))
        state = json.loads(paths[0].read_text())
        assert state['slurm_excluded_nodes'] == ['badnode', 'newbadnode']
        return [dict(status='complete_review_required', jobs=state['jobs'])]
    monkeypatch.setattr(helper.recon, 'run_many', launch)
    monkeypatch.chdir(output)  # Restore cwd automatically after this test.
    result = helper.run(output)
    assert len(calls) == 1 and calls[0][1] == dict(max_jobs=2, allocation_retries=2)
    assert not result['production_model_accepted']


def test_screen_monitor_reads_registry_without_pilot_or_manifest_writes(request_tree, monkeypatch):
    request, case_input, output = request_tree
    helper.prepare(request, case_input, RECIPE, output)
    from kinbot import orbital_monitor
    monkeypatch.setattr(orbital_monitor, '_queue', lambda: ({}, None))
    before = {p: p.read_bytes() for p in output.rglob('*') if p.is_file()}
    data = orbital_monitor.snapshot(output, registry=output/'review_launch.json', include_all=True)
    assert data['counts']['manifest']['WAITING'] == 8
    assert '8/8' not in orbital_monitor.render(data, details=True)
    assert all(p.read_bytes() == value for p, value in before.items())


def completed_screen(request_tree):
    """Synthetic native artifacts exercise provenance, never native chemistry."""
    from kinbot.execution import file_sha256
    request, case_input, output = request_tree
    helper.prepare(request, case_input, RECIPE, output)
    manifest = output/'native/reconnaissance.json'
    state = json.loads(manifest.read_text())
    for j in state['jobs']:
        d = Path(j['directory'])
        wfu, out = d/(j['name']+'.wfu'), d/(j['name']+'.out')
        wfu.write_bytes(b'synthetic checkpoint '+j['name'].encode())
        out.write_text('Synthetic output, not electronic structure evidence.\n')
        j.update(status='complete', wfu=str(wfu), wfu_sha256=file_sha256(wfu),
                 output_sha256=file_sha256(out))
    state['status'] = 'coverage_refinement_required'
    atomic_json(manifest, state)
    return output


def test_candidates_use_exact_native_partitions_matched_seeds_and_one_shared_cap(request_tree, monkeypatch):
    source = completed_screen(request_tree)
    output = source.parent/'candidate-checks'
    before = {p:p.read_bytes() for p in source.rglob('*') if p.is_file()}
    result = helper.candidates(source, CANDIDATES, output)
    assert result['tasks'] == 22 and result['models'] == 2 and result['max_jobs'] == 3
    assert all(p.read_bytes()==v for p,v in before.items())
    record = json.loads((output/'review_launch.json').read_text())
    assert len(record['systems']) == 2
    for system in record['systems']:
        config = helper.recon._load(system['configuration'])
        manifest = json.loads(Path(system['manifest']).read_text())
        compact = config['_probe_space']['active_electrons']==5
        assert len(manifest['jobs']) == (14 if compact else 8)
        for j in manifest['jobs']:
            assert j['orbital_predecessor'] is None
            assert j['initial_record']['basis'] == j['basis']
            assert j.get('status') in (None, 'prepared') and not j.get('job_id')
            deck = (Path(j['directory'])/(j['name']+'.inp')).read_text()
            assert 'start,2140.2,type=natural' in deck and 'diab,2140.2' not in deck
            assert ('closed,10' if compact else 'closed,9') in deck
            assert ('occ,14' if compact else 'occ,15') in deck
            assert deck.count('{rs2-f12,') == 6
            assert deck.count('{uccsd(t)-f12}') == int(compact)
            assert 'maxit,999;' in deck and 'refstate,4,1,2,3,4;' in deck
            if j['anchor']=='cs_r_20':
                assert 'cs_r_30' in j['initial_record']['path']
                a=next(a for a in config['anchors'] if a['name']==j['anchor'])
                xyz=helper.recon._anchor_geometry(config,a)
                assert np.linalg.norm(xyz[0]-xyz[1]) == pytest.approx(20.)
            if j['anchor']=='cs_r_3p5_sh_rot30':
                anchor=next(a for a in config['anchors'] if a['name']==j['anchor'])
                unrotated=next(a for a in config['anchors'] if a['name']=='cs_r_3p5')
                xyz=helper.recon._anchor_geometry(config,anchor)
                reference=helper.recon._anchor_geometry(config,unrotated)
                np.testing.assert_allclose(xyz[[0,1,3,4]],reference[[0,1,3,4]])
                assert not np.allclose(xyz[2],reference[2])
                assert np.linalg.norm(xyz[0]-xyz[1]) == pytest.approx(3.5)
                assert np.linalg.norm(xyz[2]-xyz[1]) == pytest.approx(
                    np.linalg.norm(reference[2]-reference[1]))
                assert 'cs_r_3p5' in j['initial_record']['path']
    with pytest.raises(FileExistsError):
        helper.candidates(source, CANDIDATES, output)
    monkeypatch.setattr(helper.subprocess,'run',lambda *a,**k:None)
    from kinbot import run_walltime
    monkeypatch.setattr(run_walltime,'partition_walltime',lambda _: '31-00:00:00')
    calls=[]
    def run_many(paths,**kwargs):
        calls.append((paths,kwargs))
        return [dict(status='complete_review_required',jobs=helper._read(p)['jobs']) for p in paths]
    monkeypatch.setattr(helper.recon,'run_many',run_many)
    assert helper.run(output)['tasks'] == 22
    assert len(calls)==1 and len(calls[0][0])==2
    assert calls[0][1]['max_jobs']==3


@pytest.mark.parametrize('corrupt', ['checkpoint','output','incomplete','counts','request'])
def test_candidate_preflight_preserves_source_and_refuses_changed_evidence(request_tree, corrupt):
    source=completed_screen(request_tree)
    manifest=source/'native/reconnaissance.json'
    state=helper._read(manifest)
    job=state['jobs'][0]
    if corrupt=='checkpoint':
        Path(job['wfu']).write_text('changed')
    elif corrupt=='output':
        (Path(job['directory'])/(job['name']+'.out')).write_text('changed')
    elif corrupt=='incomplete':
        job['status']='submitted'; atomic_json(manifest,state)
    elif corrupt=='counts':
        job['probe_active_space']['active_electrons']=11; atomic_json(manifest,state)
    else:
        p=Path(helper._read(source/'review_launch.json')['request'])
        m=helper._read(p); m['multiplicity']=4; atomic_json(p,m)
    destination=source.parent/'candidate-checks'
    with pytest.raises(ValueError):
        helper.candidates(source,CANDIDATES,destination)
    assert not destination.exists()


@pytest.mark.parametrize('rotation', [
    {'atoms':[2,3],'pivot_atom':2,'axis_bond':[1,2],'angle_degrees':0.},
    {'atoms':[2,3],'pivot_atom':2,'axis_bond':[1,1],'angle_degrees':30.},
    {'atoms':[2,3],'pivot_atom':2,'axis_bond':[1,6],'angle_degrees':30.},
    {'atoms':[2,3],'pivot_atom':1,'axis_bond':[1,2],'angle_degrees':30.},
])
def test_candidate_rotation_refuses_invalid_or_duplicate_geometry(request_tree, rotation):
    source=completed_screen(request_tree)
    recipe=helper._read(CANDIDATES)
    recipe['candidates'][0]['anchors'][-1]['rotation']=rotation
    recipe_path=source.parent/'invalid-rotation.json'
    atomic_json(recipe_path,recipe)
    output=source.parent/'invalid-candidates'
    with pytest.raises(ValueError):
        helper.candidates(source,recipe_path,output)
    assert not output.exists()


def test_saved_evidence_is_read_only_and_labels_theory(request_tree, monkeypatch):
    source=completed_screen(request_tree)
    manifest=source/'native/reconnaissance.json'
    state=helper._read(manifest)
    state['jobs'][0].update(
        state_energies_hartree=[dict(multiplicity=2,root=1,energy_hartree=-436.99)],
        correlated_energies_hartree=[dict(method='caspt2-f12(5,4)',multiplicity=2,
                                        root=1,energy_hartree=-437.1)])
    atomic_json(manifest,state)
    before={p:p.read_bytes() for p in source.rglob('*') if p.is_file()}
    monkeypatch.setattr(helper.subprocess,'run',lambda *a,**k:pytest.fail('read-only report'))
    report=helper.evidence(source)
    assert 'CASCI Eh :m2r1=-436.9900000000' in report
    assert 'CORRELATED Eh caspt2-f12(5,4):m2r1=-437.1000000000' in report
    assert 'no production model has been accepted' in report
    assert all(p.read_bytes()==v for p,v in before.items())


def test_review_publishes_only_selected_parent_and_reuses_basis_specific_seeds(request_tree, monkeypatch):
    from kinbot import input_refresh
    from kinbot.molpro import _render_vrc_state_model
    source=completed_screen(request_tree)
    _,case_input,_=request_tree
    selection=helper._read(SELECTION)
    existing=dict(schema=1,reviewed=True,review_note='Existing approved channel',channels=[dict(
        bond=[1,2],combination='single',bottlenecks=[dict(name='single',coordinate='reactive',
        distances_angstrom=[2.,3.,4.],parameters={},states=selection['states'])])])
    previous=helper._read(case_input)
    previous.update(rotdpy_dist=[2.,3.,4.],vrc_tst_surface_models={'initial_parent':existing},
                    vrc_tst_scan_basis='avdz-f12',rotdpy_flux_parameters={'pot_smp_max':1000})
    atomic_json(case_input,previous)
    source_before={p:p.read_bytes() for p in source.rglob('*') if p.is_file()}
    monkeypatch.setattr(input_refresh,'owner_directory',lambda _:case_input.parent)
    monkeypatch.setattr(helper.subprocess,'run',lambda *a,**k:pytest.fail('No scheduler or QC launch'))
    calls=[]
    def submit(path,**kw):
        calls.append((path,kw))
        assert path==case_input and kw['expected']==helper._read(path)
        if kw.get('apply'):
            atomic_json(path,kw['proposal'])
        return dict(status='queued',drivers_started=0)
    monkeypatch.setattr(input_refresh,'submit',submit)
    assert helper.review(source,case_input,SELECTION)['status']=='reviewed'
    assert helper._read(case_input)==previous and not calls
    result=helper.review(source,case_input,SELECTION,apply=True)
    assert result['reconnaissance_jobs_submitted']==0 and result['drivers_started']==0
    proposed=helper._read(case_input)
    assert {k for k in previous.keys()|proposed.keys() if previous.get(k)!=proposed.get(k)}=={'vrc_tst_surface_models'}
    assert proposed['vrc_tst_surface_models']['initial_parent']==existing
    new=proposed['vrc_tst_surface_models'][selection['request_match']['parent']]
    stage=new['channels'][0]['bottlenecks'][0]
    assert stage['distances_angstrom']==previous['rotdpy_dist']
    assert len(stage['states'])==1 and stage['states'][0]['root']==1
    assert stage['parameters']['vrc_tst_state_model_file']==''
    for basis,model in stage['parameters']['vrc_tst_state_model']['basis_models'].items():
        seed=selection['basis_seeds'][basis]
        assert seed['job'] in model['orbital_selection']['path']
        assert model['orbital_selection']['active_indices']==[11,12,13,14]
        _,deck,_=_render_vrc_state_model(model,nelectron=25,active_electrons=5,active_orbitals=4,spin=1)
        assert 'closed,10' in deck and 'occ,14' in deck
    assert all(p.read_bytes()==v for p,v in source_before.items())
    backups=list((case_input.parent/'orbital_reviews/live_refresh').glob('*/input_before.json'))
    assert len(backups)==1 and helper._read(backups[0])==previous
    # Reapplication reconciles the journal, never recreates jobs or overwrites
    # another parent. JSON int/string key normalization must be idempotent.
    helper.review(source,case_input,SELECTION,apply=True)
    assert len(calls)==4 and helper._read(case_input)==proposed
    live_request=Path(helper._read(source/'review_launch.json')['request'])
    atomic_json(live_request,dict(helper._read(live_request),status='model_accepted',surface_model={}))
    helper.review(source,case_input,SELECTION,apply=True)
    assert len(calls)==6 and helper._read(case_input)==proposed


@pytest.mark.parametrize('failure',['unreviewed','wrong_request','changed_wfu','changed_output','wrong_basis','no_owner'])
def test_review_rejects_invalid_selection_before_input_write(request_tree,monkeypatch,failure):
    from kinbot import input_refresh
    source=completed_screen(request_tree)
    _,case_input,_=request_tree
    atomic_json(case_input,dict(helper._read(case_input),rotdpy_dist=[2.,3.,4.]))
    before=case_input.read_bytes()
    recipe=helper._read(SELECTION)
    job=helper._read(source/'native/reconnaissance.json')['jobs'][0]
    if failure=='unreviewed': recipe['reviewed']=False
    elif failure=='wrong_request': recipe['request_match']['parent']='different'
    elif failure=='changed_wfu': Path(job['wfu']).write_text('changed')
    elif failure=='changed_output': (Path(job['directory'])/(job['name']+'.out')).write_text('changed')
    elif failure=='wrong_basis': recipe['basis_seeds']['avtz-f12']['job']=job['name']
    recipe_path=source.parent/'selection.json'
    atomic_json(recipe_path,recipe)
    monkeypatch.setattr(input_refresh,'owner_directory',lambda _:None if failure=='no_owner' else case_input.parent)
    monkeypatch.setattr(input_refresh,'submit',lambda *a,**kw:pytest.fail('Invalid selection'))
    with pytest.raises(ValueError): helper.review(source,case_input,recipe_path,apply=True)
    assert case_input.read_bytes()==before


def test_review_serializes_concurrent_input_edits(request_tree):
    import fcntl
    _,case_input,source=request_tree
    before=case_input.read_bytes()
    with (case_input.parent/'.kinbot_orbital_review.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(ValueError,match='Another orbital review'):
            helper.review(source,case_input,SELECTION,apply=True)
    assert case_input.read_bytes()==before


@pytest.mark.parametrize('worker_status',['exploring','ready'])
def test_review_reaches_existing_owner_and_preserves_native_work(request_tree,worker_status):
    import os
    from types import SimpleNamespace
    from kinbot import input_refresh, pes_phase, live_control
    source=completed_screen(request_tree)
    _,case_input,_=request_tree
    root=case_input.parent
    previous=dict(helper._read(case_input),rotdpy_dist=[2.,3.,4.])
    atomic_json(case_input,previous)
    live_control.configure(root.parent)
    pes_phase.initialize(root,case_input)
    parent=helper._read(SELECTION)['request_match']['parent']
    worker=root/parent
    atomic_json(worker/(parent+'.json'),previous)
    native=worker/'expensive.out';native.write_bytes(b'existing expensive work')
    pes_phase.register(root,parent,os.getpid())
    record=helper._read(root/pes_phase.NAME)
    record['workers'][parent]['status']=worker_status
    atomic_json(root/pes_phase.NAME,record)
    result=helper.review(source,case_input,SELECTION,apply=True)
    assert result['status']=='queued' and result['drivers_started']==0
    assert helper._read(root/input_refresh.SNAPSHOT)==previous
    pes_phase.tick(root,[parent],max_searches=1)
    proposed=helper._read(case_input)
    assert helper._read(root/input_refresh.SNAPSHOT)==proposed
    if worker_status=='exploring':
        generator=SimpleNamespace(inp=str(worker/(parent+'.json')),par=previous.copy())
        input_refresh.inherit_models(root,generator)
        assert generator.par['vrc_tst_surface_models']==proposed['vrc_tst_surface_models']
        assert helper._read(worker/(parent+'.json'))==proposed
    else:
        requests=live_control.requests(worker,kinds={'network'})
        assert len(requests)==1
        request=requests[0][1]
        assert not request['cancel_running']
        assert request['options']['vrc_tst_surface_models']==proposed['vrc_tst_surface_models']
    helper.review(source,case_input,SELECTION,apply=True)
    assert len(live_control.requests(worker,kinds={'network'}))==int(worker_status=='ready')
    assert native.read_bytes()==b'existing expensive work'
