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
