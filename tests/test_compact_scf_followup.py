"""Supplied chemistry fixtures; simulated checkpoints/scheduler, never native QC."""
from pathlib import Path
import re
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.db import connect

from kinbot import requested_orbital_recon as review, recon_scf_seed
from kinbot.execution import atomic_json, file_sha256
from test_remaining_orbital_screens import tree as original_tree
from test_remaining_character_panel import PLAN as CHARACTER_PLAN
from test_remaining_orbital_screens import test_campaign_launch_and_resume_share_retained_production_capacity as simulate_campaign

PLAN = Path(__file__).resolve().parents[1]/'examples/anl/orbital_review_batches/r7_compact_followup_20261025'


@pytest.fixture
def plan():
    return CHARACTER_PLAN


@pytest.fixture
def tree(original_tree):
    production, _ = original_tree
    source = production/'orbital_recon/remaining_character_20261025'
    review.batch(production, CHARACTER_PLAN/'screens.json', source, shared_pool=True)
    supplied = (PLAN/'scf_character_evidence.txt').read_text()
    # Reconstruct the printed occupation lines solely to test source transport.
    # Checkpoints are mock bytes; no native CASSCF result is asserted here.
    for system in review._saved(source)['systems']:
        manifest_path = Path(system['manifest']); manifest = review._read(manifest_path)
        child = manifest_path.parents[1].name
        for job in manifest['jobs']:
            b = next(
                b for b in supplied.split('\nNATIVE ')[1:]
                if f'/{child}/native/{job["name"]}/' in b.splitlines()[0])
            rows = re.findall(r'(?m)^SCF (\d+\.1) occ=(\S+) energy_Eh=(\S+)', b)
            out = Path(job['directory'])/(job['name']+'.out')
            out.write_text('PROGRAM * RHF-SCF\n ELECTRON ORBITALS\n\n'+''.join(
                f' {i} {float(o):.5f} {float(e):.5f} 1 1 px 0.50000\n' for i,o,e in rows)+'\nEND\n')
            wfu = out.with_suffix('.wfu'); wfu.write_bytes(('mock:'+child+job['name']).encode())
            job.update(status='complete', wfu=str(wfu), wfu_sha256=file_sha256(wfu),
                       output_sha256=file_sha256(out))
        manifest['status'] = 'complete_review_required'; atomic_json(manifest_path, manifest)
    evidence = review._read(PLAN/'additional_endpoint_evidence.json')
    for entry in review._read(PLAN/'screens.json')['screens']:
        if 'source' in entry:
            continue
        recipe = review._read(PLAN/entry['recipe']); match = recipe['request_match']
        worker = production/'ch3ss'/match['parent']
        atomic_json(worker/'vrctst/orbital_review'/match['reaction']/'request.json',
                    dict(match, schema=1, status='needs_review', positions=evidence['positions']))
        if 'geometry_source' in recipe:
            source = recipe['geometry_source']; endpoint = evidence['endpoints'][match['reaction']]
            connect(worker/'kinbot.db').write(Atoms(match['symbols'], positions=endpoint['positions']),
                name=source['record_name'], data=dict(status=endpoint['status']))
    worker = production/'ch3ss/791671450790040000002'
    rows = []
    for path in PLAN.glob('well_*.json'):
        recipe = review._read(path)
        if 'geometry_source' in recipe:
            rows.append(dict(name=recipe['request_match']['reaction'],
                             product_complex=recipe['geometry_source']['product_complex']))
    atomic_json(worker/'accepted_l2_network.json', dict(network=dict(reactions=rows)))
    return production, production/'orbital_recon/compact_followup_20261025'


@pytest.fixture
def followup_plan(tree,tmp_path):
    # Bind test recipes to synthetic transport checkpoints, preserving all
    # scientific choices and the immutable recorded real hashes in examples.
    production,_=tree
    destination=tmp_path/'followup_plan';destination.mkdir()
    for source in PLAN.glob('*.json'):
        (destination/source.name).write_bytes(source.read_bytes())
    for entry in review._read(destination/'screens.json')['screens']:
        if 'source' not in entry:
            continue
        source=production/entry['source'];record=review._saved(source)
        jobs={j['name']:j for j in review._read(record['systems'][0]['manifest'])['jobs']}
        path=destination/entry['recipe'];recipe=review._read(path)
        for candidate in recipe['candidates']:
            for anchor in candidate['anchors']:
                anchor['seed_sha256']={basis:jobs[basis+'__'+anchor.get('seed',anchor['source'])]['wfu_sha256']
                                       for basis in ('adz','atz')}
        atomic_json(path,recipe)
    return destination


def snapshot(root):
    return {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}


def test_mixed_panel_all_56_jobs_preflight_stage_resume_and_atom_maps(tree,followup_plan, monkeypatch):
    plan=followup_plan
    production, output = tree
    monkeypatch.setattr(review.subprocess, 'run', lambda *a, **k: pytest.fail('No scheduler during staging'))
    before = snapshot(production)
    result = review.batch(production, plan/'screens.json', output, shared_pool=True, preflight=True)
    assert (result['screens'], result['tasks']) == (12, 56)
    assert not output.exists() and snapshot(production) == before
    result = review.batch(production, plan/'screens.json', output, shared_pool=True)
    assert (result['screens'], result['tasks'], result['max_jobs']) == (12, 56, 3)
    assert all(p.read_bytes() == v for p,v in before.items())
    n_cas = n_scf = 0
    evidence = review._read(PLAN/'additional_endpoint_evidence.json')
    for system in review._saved(output)['systems']:
        config = review.recon._load(system['configuration'])
        manifest = review._read(system['manifest'])
        assert manifest['shared_policy_file'] == str(production/'.kinbot_resource_policy.json')
        for job in manifest['jobs']:
            assert not job.get('job_id') and not job.get('orbital_predecessor')
            deck = (Path(job['directory'])/(job['name']+'.inp')).read_text()
            if config.get('_probe_space'):
                n_cas += 1
                space = config['_probe_space']; seed = job['initial_record']
                assert space['active_orbitals'] <= 4 and space['active_electrons'] <= 5
                assert '{multi,so-sci' in deck and 'start,2100.2,type=canonical,ovl' in deck
                assert 'rs2' not in deck.lower() and 'diab,' not in deck.lower()
                recon_scf_seed.validate(seed, 41, space['active_electrons'], space['active_orbitals'])
            else:
                n_scf += 1
                assert '{multi' not in deck.lower() and 'rs2' not in deck.lower()
        child = Path(system['configuration']).parent
        recipe = review._read(child/'recipe.json')
        if 'geometry_source' in recipe and 'candidates' not in recipe:
            endpoint = evidence['endpoints'][recipe['request_match']['reaction']]['positions']
            assert np.allclose(config['geometry_angstrom'], endpoint)
            assert np.allclose(review.recon._anchor_geometry(config, config['anchors'][0]), endpoint)
            for anchor in config['anchors']:
                xyz = review.recon._anchor_geometry(config, anchor)
                for f in recipe['geometry_source']['fragments']:
                    m = np.array(f['atoms'])-1
                    distances = lambda a: np.linalg.norm(a[:,None]-a[None,:], axis=2)
                    assert np.allclose(distances(xyz[m]), distances(np.array(endpoint)[m]))
    assert (n_cas,n_scf) == (24,32)
    staged = snapshot(production)
    assert review.batch(production,plan/'screens.json',output,shared_pool=True,resume=True)['tasks'] == 56
    assert snapshot(production) == staged


def test_bad_last_endpoint_leaves_everything_untouched(tree,followup_plan):
    plan=followup_plan
    production, output = tree
    recipe = review._read(PLAN/'well_r12_insertion_r_6_1_5.json')
    worker = production/'ch3ss'/recipe['request_match']['parent']
    xyz = review._read(PLAN/'additional_endpoint_evidence.json')['endpoints'][recipe['request_match']['reaction']]['positions']
    connect(worker/'kinbot.db').write(Atoms('CSSHHH', positions=xyz),
        name=recipe['geometry_source']['record_name'], data=dict(status='error'))
    before = snapshot(production)
    with pytest.raises(ValueError, match='normal completion'):
        review.batch(production,plan/'screens.json',output,shared_pool=True)
    assert not output.exists() and snapshot(production) == before


def test_changed_seed_blocks_preflight_and_resume(tree,followup_plan):
    plan=followup_plan
    production,output = tree
    review.batch(production,plan/'screens.json',output,shared_pool=True)
    system = review._saved(output)['systems'][0]
    config = review._read(system['configuration'])
    seed = next(iter(config['initial_records'].values()))
    Path(seed['source_output']).write_text('Changed source evidence')
    before = snapshot(production)
    with pytest.raises(ValueError, match='evidence changed'):
        review.batch(production,plan/'screens.json',output,shared_pool=True,resume=True)
    assert snapshot(production) == before


def test_changed_seed_blocks_run_before_any_scheduler_call(tree,followup_plan,monkeypatch):
    production,output=tree
    review.batch(production,followup_plan/'screens.json',output,shared_pool=True)
    system=review._saved(output)['systems'][0]
    seed=next(iter(review._read(system['configuration'])['initial_records'].values()))
    Path(seed['path']).write_bytes(b'changed checkpoint')
    monkeypatch.setattr(review.subprocess,'run',lambda *a,**k:pytest.fail('No allocation for a changed seed'))
    with pytest.raises(ValueError,match='evidence changed'):
        review.run(output)


def test_native_scf_seed_is_not_overwritten_by_new_rhf(tree,followup_plan, tmp_path, monkeypatch):
    plan=followup_plan
    production,output = tree
    review.batch(production,plan/'screens.json',output,shared_pool=True)
    manifest_path = Path(review._saved(output)['systems'][0]['manifest'])
    job = review._read(manifest_path)['jobs'][0]
    original = Path(job['initial_record']['path']).read_bytes()
    monkeypatch.setattr(review.recon, 'qc_runtime_environment', lambda *a,**k:
        ({}, {'scratch_directory': str(tmp_path/'scratch')}))
    monkeypatch.setattr(review.recon, 'cleanup_qc_runtime', lambda *a: None)
    (tmp_path/'scratch').mkdir()
    calls = []
    def execute(command, **kwargs):
        calls.append(command)
        local = kwargs['cwd']; deck = (local/command[-1]).read_text()
        source = re.search(r'(?m)^file,3,(\S+)$',deck)[1]
        target = re.search(r'(?m)^file,2,([^,\n]+),new$',deck)[1]
        assert (local/source).read_bytes() == original
        assert ',ovl' not in deck
        assert deck.count('start,2100.3,type=canonical') == 2
        assert 'rotate,19.1,17.1,0' in deck
        (local/target).write_bytes(b'new RHF and MULTI records')
        assert (local/source).read_bytes() == original
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(review.recon.subprocess,'run',execute)
    review.recon.execute(manifest_path,job['name'])
    assert len(calls)==1 and Path(job['initial_record']['path']).read_bytes()==original
    execution = review._read(Path(job['directory'])/'execution.json')
    assert execution['scf_seed']['source_record']=='2100.3'
    assert execution['scf_seed']['transport']=='same_basis_coefficients'
    assert execution['status']=='complete'


@pytest.mark.parametrize('indices',[[17,18,20,22],[17,18,21,21],[17,18,21,999]])
def test_scf_count_somo_and_missing_orbital_guards(tree,followup_plan,indices):
    plan=followup_plan
    production,output=tree
    recipe=review._read(plan/'ch_h2ss_cs_candidates.json')
    recipe['candidates'][0]['anchors'][0]['active_indices']['adz']=indices
    path=production/'bad_recipe.json';atomic_json(path,recipe)
    with pytest.raises(ValueError,match='electron count|distinct|absent'):
        review.candidates(production/'orbital_recon/remaining_character_20261025/ch_h2ss_cs',path,output)
    assert not output.exists()


def test_mixed_campaign_launch_and_resume_preserve_capacity(tree,followup_plan,monkeypatch):
    plan=followup_plan
    # Reuse the scheduler contract with this 56-job mixed campaign.
    simulate_campaign(tree,plan,monkeypatch)


def test_late_staging_failure_rolls_back_only_new_campaign_then_retry_works(tree,followup_plan,monkeypatch):
    plan=followup_plan
    production,output=tree
    before=snapshot(production)
    original=review.recon.prepare
    def failed(source,destination):
        destination=Path(destination)
        if output in destination.parents and 'well_r12_insertion_r_6_1_5' in str(destination):
            raise OSError('simulated staging interruption after earlier children prepared')
        return original(source,destination)
    monkeypatch.setattr(review.recon,'prepare',failed)
    with pytest.raises(OSError,match='staging interruption'):
        review.batch(production,plan/'screens.json',output,shared_pool=True)
    assert not output.exists() and snapshot(production)==before
    monkeypatch.setattr(review.recon,'prepare',original)
    assert review.batch(production,plan/'screens.json',output,shared_pool=True,resume=True)['tasks']==56


def test_competing_preparer_cannot_remove_an_owned_campaign(tmp_path):
    output=tmp_path/'campaign'
    with review._new_campaign(output):
        marker=output/'owned';marker.write_text('first owner')
        with pytest.raises(FileExistsError):
            with review._new_campaign(output):
                pytest.fail('Must not enter another preparer staging area')
        assert marker.read_text()=='first owner'


def test_recomputed_scf_cannot_silently_reuse_the_reviewed_indices(tree,followup_plan):
    production,output=tree
    source=production/'orbital_recon/remaining_character_20261025/ch_h2ss_cs'
    path=Path(review._saved(source)['systems'][0]['manifest']);manifest=review._read(path)
    job=next(j for j in manifest['jobs'] if j['name']=='adz__cs_r_30')
    Path(job['wfu']).write_bytes(b'different converged SCF solution')
    job['wfu_sha256']=file_sha256(job['wfu']);atomic_json(path,manifest)
    with pytest.raises(ValueError,match='explicitly reviewed seed hash'):
        review.candidates(source,followup_plan/'ch_h2ss_cs_candidates.json',output,preflight=True)
    assert not output.exists()


def test_example_hashes_are_the_user_supplied_native_records():
    supplied=(PLAN/'scf_character_evidence.txt').read_text()
    hashes={}
    for block in supplied.split('\nNATIVE ')[1:]:
        path=Path(block.splitlines()[0])
        hashes[(path.parts[-4],path.stem)]=re.search(r'(?m)^WFU \S+ SHA256 (\w+)$',block)[1]
    for path in PLAN.glob('*_candidates.json'):
        name=path.name.removesuffix('_candidates.json')
        for candidate in review._read(path)['candidates']:
            for anchor in candidate['anchors']:
                for basis,digest in anchor['seed_sha256'].items():
                    assert digest==hashes[(name,basis+'__'+anchor.get('seed',anchor['source']))]
