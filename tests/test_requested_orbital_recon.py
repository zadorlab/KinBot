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


def test_ch_h2s_screen_uses_requested_geometry_and_only_doublet_casci(request_tree):
    request, case_input, output = request_tree
    recipe = RECIPE.with_name('ch_h2s_cs_screen.json')
    positions = [[.7467252173, -.7423275741, -.5474233962],
        [.1678425237, .4711769182, .5418506859],
        [1.0595729495, 1.3596853597, 1.1290564624],
        [.0422254396, -.2114253717, 1.6944421867],
        [1.1828471196, -.1760172776, -1.3830416653]]
    data = dict(json.loads(recipe.read_text())['request_match'],
                schema=1, status='needs_review', positions=positions,
                smiles='[H][C]=[S]([H])[H]')
    atomic_json(request, data)
    result = helper.prepare(request, case_input, recipe, output)
    assert result['tasks'] == 8
    assert json.loads((output/'request.json').read_text())['positions'] == positions
    manifest = json.loads((output/'native/reconnaissance.json').read_text())
    for job in manifest['jobs']:
        assert [(s['multiplicity'], s['root']) for s in job['root_states']] == [
            (2, 1), (2, 2), (2, 3), (2, 4)]
        deck = (Path(job['directory'])/(job['name']+'.inp')).read_text().lower()
        assert '{rhf,so-sci,maxit=256;wf,25,1,1;' in deck
        assert '{multi,cas-ci;' in deck
        assert 'rs2' not in deck and 'ccsd' not in deck


def test_ch_h2s_reference_panel_uses_six_matched_jobs_and_production_seed(request_tree, monkeypatch):
    request, _, source = request_tree
    recipe = RECIPE.with_name('ch_h2s_cs_screen.json')
    saved = helper._read(request)
    saved.update(helper._read(recipe)['request_match'])
    # CH has H5, H2S has H3/H4; retain the exact request coordinates.
    saved['positions'] = [[.7467252173,-.7423275741,-.5474233962],
        [.1678425237,.4711769182,.5418506859], [1.0595729495,1.3596853597,1.1290564624],
        [.0422254396,-.2114253717,1.6944421867], [1.1828471196,-.1760172776,-1.3830416653]]
    atomic_json(request, saved)
    source = completed_screen(request_tree, recipe)
    before = {p:p.read_bytes() for p in source.rglob('*') if p.is_file()}
    monkeypatch.setattr(helper.subprocess, 'run', lambda *a, **k: pytest.fail('No launch'))
    output = source.parent/'reference-check'
    result = helper.candidates(source, recipe.with_name('ch_h2s_cs_reference_check.json'), output)
    assert result['tasks'] == 6 and result['max_jobs'] == 3
    system = helper._read(output/'review_launch.json')['systems'][0]
    data = helper.recon._load(system['configuration'])
    selection = helper._read(recipe.with_name('ch_h2s_cs_selected.json'))
    assert data['probe_definition']['state_model'] == selection['state_model']
    manifest = helper._read(system['manifest'])
    assert {j['distance_angstrom'] for j in manifest['jobs']} == {2.,3.5,30.}
    for job in manifest['jobs']:
        assert job['initial_record']['path'].endswith('__cs_r_2.wfu')
        assert job['orbital_predecessor'] is None
        deck = (Path(job['directory'])/(job['name']+'.inp')).read_text()
        assert 'rotate,8.1,11.1,0' in deck and 'closed,10' in deck and 'occ,14' in deck
        assert 'wf,25,1,1\nstate,4' in deck and 'wf,25,1,3\nstate,1' in deck
        assert deck.count('{rs2-f12,') == 5 and deck.count('{uccsd(t)-f12b,') == 1
        assert 'orbital,3000.2\ndm\n' in deck
        assert 'load,dr1,den,3000.2,state=1.1,ms2=3' in deck
        assert deck.index('orbital,3000.2') < deck.index('canonical,2141.2')
        assert deck.index('load,dhf,den,2150.2') < deck.index('{uccsd(t)-f12b,')
        assert 'wf,25,1,3;occ,14;closed,11;orbital,2150.2' in deck
        assert '.reference_c6_m4.molden' in deck
        anchor = next(a for a in data['anchors'] if a['name'] == job['anchor'])
        xyz = helper.recon._anchor_geometry(data, anchor)
        assert np.linalg.norm(xyz[1]-xyz[0]) == pytest.approx(job['distance_angstrom'])
    assert all(p.read_bytes() == v for p,v in before.items())

    # A report is observational, including finite nonzero asymptotic gaps.
    from kinbot.execution import file_sha256
    for j in manifest['jobs']:
        distance, basis = j['distance_angstrom'], j['basis']
        offset = .01 if basis == 'avdz-f12' else 0.
        values = [-437.+offset-.2/distance, -437.1+offset-.3/distance,
                  -437.02+offset-.1/distance]
        j.update(status='complete', correlated_energies_hartree=[
            dict(method=method,multiplicity=mult,root=1,energy_hartree=e)
            for method,mult,e in zip(('uccsd(t)-f12b','caspt2-f12(5,4)','caspt2-f12(5,4)'),
                                     (4,2,4),values)])
        out = Path(j['directory'])/(j['name']+'.out')
        out.write_text('Synthetic native output; not chemical validation.\n')
        j['output_sha256'] = file_sha256(out)
    atomic_json(Path(system['manifest']), manifest)
    before = {p:p.read_bytes() for p in output.rglob('*') if p.is_file()}
    report = helper.reference_report(output)
    assert '-117.135102' in report  # (-.4/2+.4/30)*627.509474
    assert 'No production input changed' in report
    assert all(p.read_bytes() == v for p,v in before.items())


def test_production_evidence_uses_accepted_native_points_and_undisplaced_section(tmp_path, monkeypatch):
    import gzip
    from kinbot.execution import file_sha256
    root = tmp_path/'production'
    well = root/'ch3s/123'
    request = well/'vrctst/scan_review/123_hom_sci_1_2__single_m2_r1/scan_review.json'
    model = {'state_average':[dict(multiplicity=2,roots=3)], 'selected_roots':{'2':2}}
    native = gzip.open(ROOT/'tests/fixtures/molpro_recon/caspt2_native_legacy_selector.out.gz','rt').read()
    xyz = '2\nsynthetic geometry\nH 0 0 0\nH 0 0 2\n'
    points = {}
    for distance in (2.,4.,30.):
        graph = well/'vrctst/molpro_scan'/f'point_{distance}'
        d = graph/'tasks/relaxed_scan'
        atomic_json(d/'task.json', {'task':{'profile':dict(method='caspt2-f12(5,4)',
            basis='avdz-f12', calculator_kwargs=dict(state_model=model))}})
        (d/'final.xyz').write_text(xyz)
        (d/'relaxed_scan_step_0001.xyz').write_text(xyz)
        (d/'relaxed_scan_step_0001.inp').write_text('synthetic input')
        (d/'relaxed_scan_step_0001.out').write_text(native+'\n'+native.replace('PROGRAM * MULTI','PROGRAM * UNUSED'))
        atomic_json(graph/'workflow.json', {})
        atomic_json(d/'execution.json', {'status':'executed','details':{
            'generated_files':['relaxed_scan_step_0001.inp']}})
        points[str(distance)] = dict(status='accepted', graph=str(graph),distance_angstrom=distance,
            hashes={str(p):file_sha256(p) for p in (graph/'workflow.json',d/'execution.json',d/'final.xyz')})
    points['rejected'] = dict(status='needs_review',graph='must_not_read',distance_angstrom=5.)
    atomic_json(request,dict(status='running',reaction=request.parent.name,
        context={'separation_angstrom':30.},points=points))
    before = {p:p.read_bytes() for p in root.rglob('*') if p.is_file()}
    monkeypatch.setattr(helper.subprocess, 'run', lambda *a, **k: pytest.fail('Read only'))
    report = helper.production_evidence(root)
    assert 'POINT 4.0' in report and 'POINT 30.0' in report and 'POINT 2.0' not in report
    assert '"selected_roots": {"2": 2}' in report
    assert 'NATIVE STATES' in report and 'ACTIVE ' in report
    assert 'state agreement is NOT verified' in report
    assert all(p.read_bytes() == v for p,v in before.items())


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


def completed_screen(request_tree, recipe=RECIPE):
    """Synthetic native artifacts exercise provenance, never native chemistry."""
    from kinbot.execution import file_sha256
    request, case_input, output = request_tree
    helper.prepare(request, case_input, recipe, output)
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


def test_ch_h2s_selection_repacks_donor_and_uses_auxiliary_quartet_cbs(request_tree, monkeypatch):
    from kinbot import input_refresh
    from kinbot.molpro import _render_vrc_state_model
    from kinbot.vrc_correction import correction_components, expected_levels, spin_splitting_settings
    from kinbot.vrc_multistate import state_parameters
    _, case_input, output = request_tree
    recipe = RECIPE.with_name('ch_h2s_cs_screen.json')
    selected = RECIPE.with_name('ch_h2s_cs_selected.json')
    choice = helper._read(selected)
    match = choice['request_match']
    request = case_input.parent/match['parent']/'vrctst/orbital_review'/match['reaction']/'request.json'
    atomic_json(request, dict(match, schema=1, status='needs_review',
        smiles='[H][C]=[S]([H])[H]', positions=[
        [.7467252173, -.7423275741, -.5474233962],
        [.1678425237, .4711769182, .5418506859],
        [1.0595729495, 1.3596853597, 1.1290564624],
        [.0422254396, -.2114253717, 1.6944421867],
        [1.1828471196, -.1760172776, -1.3830416653]]))
    source = completed_screen((request,case_input,output), recipe)
    prior_channel = dict(schema=1, reviewed=True, review_note='Existing channel',
        channels=[dict(bond=[1,2], combination='single', bottlenecks=[dict(
            name='single', coordinate='reactive', distances_angstrom=[2.,3.,4.],
            parameters={}, states=choice['states'])])])
    previous = dict(helper._read(case_input), rotdpy_dist=[2.,3.,4.],
        vrc_tst_surface_models={'existing_parent':prior_channel},
        vrc_tst_scan_basis='avdz-f12', vrc_tst_sample_basis='avdz-f12',
        vrc_tst_high_basis=['avdz-f12','avtz-f12'], vrc_tst_high_extrapolation_power=3.7,
        vrc_tst_correction_model='spin_splitting',
        vrc_tst_spin_splitting={'reference_multiplicity':4},
        vrc_tst_state_average_roots={'2':3,'4':3}, vrc_tst_selected_roots={'2':1,'4':1})
    atomic_json(case_input, previous)
    before = {p:p.read_bytes() for p in source.rglob('*') if p.is_file()}
    monkeypatch.setattr(input_refresh, 'owner_directory', lambda _:case_input.parent)
    monkeypatch.setattr(helper.subprocess, 'run', lambda *a,**k:pytest.fail('No new QC or driver'))
    publications = []
    def submit(path, **kw):
        if kw.get('apply'):
            publications.append(kw['proposal'])
            atomic_json(path, kw['proposal'])
        return dict(status='queued', drivers_started=0)
    monkeypatch.setattr(input_refresh, 'submit', submit)
    result = helper.review(source, case_input, selected, apply=True)
    assert result['drivers_started'] == result['reconnaissance_jobs_submitted'] == 0
    proposed = helper._read(case_input)
    assert {k for k in previous.keys()|proposed.keys() if previous.get(k)!=proposed.get(k)} == {'vrc_tst_surface_models'}
    assert proposed['vrc_tst_surface_models']['existing_parent'] == prior_channel
    stage = proposed['vrc_tst_surface_models'][match['parent']]['channels'][0]['bottlenecks'][0]
    assert len(stage['states']) == 1 and stage['states'][0]['root'] == 1
    leaf = state_parameters(proposed, stage, stage['states'][0], [2,1])
    assert leaf['vrc_tst_correction_model'] == 'spin_splitting'
    assert leaf['vrc_tst_selected_roots'] == {'2':1,'4':1}
    spin = spin_splitting_settings(leaf, parent_multiplicity=2, fragment_multiplicities=[2,1])
    assert spin['reference_multiplicity'] == 4
    components = correction_components(leaf, symbols=match['symbols'], settings=spin,
                                       state_roots={2:4,4:1})
    assert [c['role'] for c in components] == ['sampling','reference:0','reference:1',
        'gap:caspt2_f12:target:0','gap:caspt2_f12:target:1',
        'gap:caspt2_f12:reference:0','gap:caspt2_f12:reference:1']
    for component in components:
        assert component['method'] == ('uccsd(t)-f12b' if component['role'].startswith('reference:')
                                       else 'caspt2-f12(5,4)')
        model = component['state_model']
        assert model['state_average'] == [dict(multiplicity=2,roots=4,weights=[.25]*4),
                                         dict(multiplicity=4,roots=1,weights=[1.])]
        assert model['orbital_selection']['active_indices'] == [8,12,13,14]
        assert '__cs_r_2.wfu' in model['orbital_selection']['path']
        _, deck, _ = _render_vrc_state_model(model, nelectron=25, active_electrons=5, active_orbitals=4, spin=1)
        assert 'closed,10' in deck and 'occ,14' in deck and 'wf,25,1,3' in deck
        assert 'rotate,8.1,11.1,0' in deck or 'rotate,11.1,8.1,0' in deck
    levels = expected_levels(leaf, symbols=match['symbols'], settings=spin, state_roots={2:4,4:1})
    assert levels['trusted_correction']['reference']['basis'] == previous['vrc_tst_high_basis']
    assert levels['trusted_correction']['reference']['extrapolation_power'] == 3.7
    assert stage['distances_angstrom'] == previous['rotdpy_dist']
    assert all(p.read_bytes()==value for p,value in before.items())
    helper.review(source, case_input, selected, apply=True)
    assert publications == [proposed,proposed]


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
def test_review_reaches_existing_owner_and_preserves_native_work(request_tree,worker_status,monkeypatch):
    import os
    from types import SimpleNamespace
    from kinbot import input_refresh, pes_phase, live_control
    # Exercise the real owner/worker identity comparisons without depending
    # on the host's ps latency or permission policy (Linux /proc is separate).
    monkeypatch.setattr(pes_phase,'process_identity',lambda pid:
        dict(host='test-host',pid=pid,start='fixed-start',boot='test-boot'))
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


@pytest.mark.parametrize('far_in_points', [True, False])
def test_production_checks_run_actual_saved_models_under_one_cap(request_tree, monkeypatch, far_in_points):
    from copy import deepcopy
    from ase import Atoms
    from ase.io import write
    from kinbot.execution import file_sha256
    from kinbot.molpro import validate_vrc_state_model
    request, case_input, source = request_tree
    source = completed_screen(request_tree)
    output = source.parent/'combined'
    helper.candidates(source, CANDIDATES, output)
    well = request.parents[3]
    # This fixture has an artificial accepted D2 leaf with two Q roots and
    # minimum-Q energy selection: staging must not silently turn it into D1/Q1.
    selection = helper._read(SELECTION)
    model = deepcopy(selection['state_model'])
    model['selected_roots']['2'] = 2
    model['energy_selection'] = {'4':'lowest'}
    model = validate_vrc_state_model(model,active_electrons=5,active_orbitals=4,nelectron=25)
    model = json.loads(json.dumps(model))
    state = dict(selection['states'][0], root=2)
    surface = dict(schema=1, reviewed=True, review_note='Synthetic test', channels=[dict(
        bond=[1,2],combination='single',bottlenecks=[dict(name='single',coordinate='reactive',
        distances_angstrom=[2.,4.,20.],parameters={},states=[state])])])
    par = dict(helper._read(case_input), vrc_tst_surface_model=surface,
        vrc_tst_active_space=selection['active_space'], vrc_tst_state_model=model,
        vrc_tst_selected_roots={'2':2,'4':1}, vrc_tst_correction_model='spin_splitting',
        vrc_tst_spin_splitting=dict(reference=dict(method='uccsd(t)-f12b',basis=['avdz-f12','avtz-f12'],
            extrapolation_power=3.7), gap_levels=[dict(label='caspt2',method='caspt2-f12(5,4)',
                basis=['avdz-f12','avtz-f12'],extrapolation_power=3.7,shift=.1,ipea=.25)],
            weighting={'model':'arithmetic_mean'}),
        vrc_tst_sample_method='caspt2-f12(5,4)',vrc_tst_sample_basis='avdz-f12')
    atomic_json(well/(well.name+'.json'),par)
    leaf = well.name+'_hom_sci_1_2__single_m2_r2'
    live = well/'vrctst/orbital_review'/(well.name+'_hom_sci_1_2')/'request.json'
    atomic_json(live,helper._read(request))
    points = {}
    for distance in (2.,4.,30.):
        graph = well/'vrctst/molpro_scan'/(f'{leaf}_vts_fragments_asymptote_test'
                                        if distance == 30. else f'point_{distance}')
        d = graph/'tasks/relaxed_scan'; d.mkdir(parents=True)
        wfu = d/'native.wfu';wfu.write_text('immutable synthetic checkpoint')
        checkpoint = dict(model, orbital_selection=dict(mode='molpro_record',path=str(wfu),
            record='2140.2',sha256=file_sha256(wfu),orbital_type='natural',active_indices=[11,12,13,14],
            align_active_orbitals=False))
        profile = dict(method='caspt2-f12(5,4)',basis='avdz-f12',calculator_kwargs=dict(
            state_model=checkpoint,caspt2_shift=.1,caspt2_ipea=.25,scf_convergence=1.e-8,
            correlation_convergence=1.e-7))
        atoms = Atoms('CSHHH',positions=helper._read(request)['positions'])
        for name in ('final.xyz','relaxed_scan_step_0001.xyz'): write(d/name,atoms,format='xyz')
        (d/'relaxed_scan_step_0001.inp').write_text('synthetic native force evaluation')
        (d/'relaxed_scan_step_0001.out').write_text('synthetic native force output')
        atomic_json(d/'task.json',dict(task=dict(profile=profile)))
        atomic_json(graph/'workflow.json',dict(molecule=dict(symbols=atoms.get_chemical_symbols(),
                                                           positions=atoms.positions.tolist(),charge=0,multiplicity=2)))
        artifacts = {p.name:file_sha256(p) for p in d.iterdir() if p.is_file()}
        atomic_json(d/'execution.json',dict(status='executed',artifacts=artifacts,
            details=dict(generated_files=['relaxed_scan_step_0001.inp','native.wfu'],orbital_checkpoint_model=checkpoint)))
        points[str(distance)] = dict(status='accepted',graph=str(graph),distance_angstrom=distance,
            hashes={str(p):file_sha256(p) for p in (graph/'workflow.json',d/'final.xyz',d/'execution.json')})
    points['bad'] = dict(status='needs_review',distance_angstrom=6.,graph='/never/use')
    if not far_in_points:
        points.pop('30.0')
    atomic_json(well/'vrctst/scan_review'/leaf/'scan_review.json',
        dict(reaction=leaf,status='running',context=dict(separation_angstrom=30.,
             state_model=checkpoint,fragment_geometry=atoms.positions.tolist()),points=points))
    production = case_input.parent.parent
    before = {p:p.read_bytes() for p in production.rglob('*') if p.is_file()}
    monkeypatch.setattr(helper.subprocess,'run',lambda *a,**k:pytest.fail('Staging must not launch anything'))
    result = helper.production_checks(production,output)
    assert result['added_tasks'] == 2 and result['skipped'] == [] and result['max_jobs'] == 3
    record = helper._saved(output)
    for system in record['systems'][-2:]:
        m = helper._read(system['manifest']);j = m['jobs'][0]
        assert len(m['jobs']) == 1
        assert j['initial_record']['path'].endswith('/native.wfu')
        assert j['state_model']['selected_roots']['2'] == 2
        assert j['state_model']['energy_selection']['4'] == 'lowest'
        assert j['state_model']['state_average'] == model['state_average']
        assert system['production_probe']['distance_angstrom'] in (4.,30.)
        deck = (Path(j['directory'])/(j['name']+'.inp')).read_text()
        assert 'rotate,' not in deck and 'diab,' not in deck and 'force' not in deck
        assert 'orbital,3000.2\ndm\n' in deck
        assert 'load,dr2,den,3000.2,state=2.1,ms2=3' in deck
        assert '{uccsd(t)-f12b,thrden=1.d-7,thrvar=1.d-10}' in deck
    assert helper.production_checks(production,output)['status'] == 'already_prepared'
    assert all(p.read_bytes()==value for p,value in before.items())
    monkeypatch.setattr(helper.subprocess,'run',lambda *a,**k:None)
    from kinbot import run_walltime
    monkeypatch.setattr(run_walltime,'partition_walltime',lambda _: '31-00:00:00')
    calls = []
    def launch(paths, **kwargs):
        calls.append((paths,kwargs))
        return [dict(status='complete_review_required',jobs=helper._read(p)['jobs']) for p in paths]
    monkeypatch.setattr(helper.recon,'run_many',launch)
    monkeypatch.chdir(output)
    monkeypatch.setattr(helper, '_reference_preflight', lambda *a: None)
    helper.run(output)
    assert len(calls) == 1 and len(calls[0][0]) == len(record['systems'])
    assert calls[0][1] == dict(keep_going=True,max_jobs=3,allocation_retries=2)

    # Extending an idle panel fills missing roles but must not duplicate or
    # rewrite the two already prepared production geometries/checkpoints.
    from types import SimpleNamespace
    monkeypatch.setattr(helper.subprocess,'run',lambda *a,**k:SimpleNamespace(stdout=''))
    original = [helper._read(s['manifest']) for s in record['systems']]
    extended = helper.production_checks(production,output,extend=True)
    assert extended['added_tasks'] == 2 and extended['total_tasks'] == result['total_tasks']
    assert [helper._read(s['manifest']) for s in record['systems']] == original
    path = Path(record['systems'][0]['manifest']); m=helper._read(path)
    m['jobs'][0]['job_id']='42'; atomic_json(path,m)
    monkeypatch.setattr(helper.subprocess,'run',lambda *a,**k:SimpleNamespace(stdout='42\n'))
    with pytest.raises(RuntimeError,match='allocations'):
        helper.production_checks(production,output,extend=True)


@pytest.mark.parametrize('maps', [[[0,3,4,5],[1,2]], [[0,1,3,4,5],[2]]])
def test_ch3ss_fragment_states_are_available_before_anl_network(maps):
    from kinbot.vrc_reference_probes import context_fragment_multiplicities
    from kinbot.vrc_scan_recovery import digest
    parameters = helper._read(ROOT/'examples/anl/ch3ss_anl1_f12/ch3ss.json')
    xyz = np.array([[1.0521,.0814,-.0820],[.4164,-.6196,-1.6231],
                    [1.1718,.6443,-3.0156],[.7032,-.5284,.7560],
                    [.6877,1.1027,.0633],[2.1461,.0805,-.0733]])
    xyz[maps[1]] += [30.,0.,0.]
    context = dict(symbols=['C','S','S','H','H','H'],charge=0,multiplicity=2,
                   maps=maps,fragment_geometry=xyz.tolist())
    saved = dict(context=context,context_sha256=digest(context))
    assert context_fragment_multiplicities(parameters,saved) == [2,3]
    with pytest.raises(ValueError,match='explicit electronic-state'):
        context_fragment_multiplicities({},saved)
    context['maps'] = [[0],[1,2,3,4,5]]
    saved['context_sha256'] = digest(context)
    with pytest.raises(ValueError,match='connectivity'):
        context_fragment_multiplicities(parameters,saved)


@pytest.mark.parametrize('selector', ['spin', 'ms2'])
def test_reference_export_retry_preserves_input_and_failed_evidence(request_tree, monkeypatch, selector):
    from types import SimpleNamespace
    from kinbot.execution import file_sha256
    source = completed_screen(request_tree)
    output = source.parent/'retry-panel'
    helper.candidates(source,CANDIDATES,output)
    record=helper._saved(output); path=Path(record['systems'][0]['manifest'])
    m=helper._read(path); j=m['jobs'][0]; d=Path(j['directory'])
    inp=d/(j['name']+'.inp')
    inp.write_text(inp.read_text()+'\n{multi;\nnatorb,3001.2,state=1.1,spin=3\n}\n'
                   f'{{matrop\nload,dr1,den,3001.2,state=1.1,{selector}=3\n}}\n')
    j.update(input_sha256=file_sha256(inp),job_id='42',status='incomplete')
    atomic_json(path,m)
    failed = '?REQUESTED DENSITY NOT FOUND IN RECORD    3001.2 FOR STATE= SYM= 1 MS2= 6\n'
    (d/(j['name']+'.out')).write_text(failed)
    atomic_json(d/'execution.json',dict(status='failed',error=failed,input_sha256=j['input_sha256']))
    original = inp.read_bytes()
    monkeypatch.setattr(helper.subprocess,'run',lambda *a,**k:SimpleNamespace(stdout='',stderr='',returncode=0))
    result=helper.retry_reference_exports(output)
    assert result['manifests'][0]['retry_ready'] == [j['name']]
    after=helper._read(path)['jobs'][0]
    assert after['attempt'] == 2 and 'job_id' not in after
    assert inp.read_bytes() == original
    assert (d/'attempts/1'/(j['name']+'.out')).read_text() == failed
    assert helper.retry_reference_exports(output)['manifests'] == []


@pytest.mark.parametrize('failure', [None, 'native', 'missing_root', 'bad_density'])
def test_reference_preflight_uses_one_job_and_checks_native_matrices(tmp_path, monkeypatch, failure):
    from kinbot.execution import file_sha256
    from kinbot import molpro_reconnaissance as recon
    path = tmp_path/'one/reconnaissance.json'
    other = tmp_path/'other/reconnaissance.json'
    d = path.parent/'probe'; d.mkdir(parents=True)
    job = dict(name='probe', status='complete', directory=str(d),
        correlation_checks=[dict(method='uccsd(t)-f12b',multiplicity=4)],
        state_model={'state_average':[dict(multiplicity=4,roots=2)]},
        probe_active_space=dict(active_electrons=3,active_orbitals=3))
    manifest = dict(jobs=[job],electron_count=5)
    atomic_json(path,manifest)
    record = {'systems':[dict(title='one',manifest=str(path),production_probe={'channel':'test'})]}
    matrices = {k:np.eye(4) for k in ('s','cas','hf')}
    matrices.update({k:np.diag([2.,1.,1.,1.]) for k in ('dhf','dr1','dr2')})
    if failure == 'missing_root': matrices.pop('dr2')
    if failure == 'bad_density': matrices['dr1'][0,0] = 0.
    files = {}
    for key, value in matrices.items():
        p=d/('probe.identity_'+key+'.out')
        p.write_text('BEGIN_DATA\n'+'\n'.join(' '.join(map(str,row)) for row in value)+'\nEND_DATA\n')
        files[p.name]=file_sha256(p)
    atomic_json(d/'execution.json',dict(reference_identity_matrices=files))
    calls=[]
    def launch(paths, **kwargs):
        calls.append((paths,kwargs))
        if failure == 'native': raise recon.ProbeFailure('native failure')
        return [manifest]
    monkeypatch.setattr(recon,'run_many',launch)
    if failure:
        with pytest.raises((RuntimeError, ValueError)):
            helper._reference_preflight(record,[path,other])
    else:
        helper._reference_preflight(record,[path,other])
    assert calls == [([path],dict(max_jobs=1,allocation_retries=2,capacity_manifests=[other]))]
