"""Compact native hypotheses and simulated orchestration; no native QC claims."""
from pathlib import Path
import importlib.util
import json
import re

import numpy as np
import pytest

from kinbot import requested_orbital_recon as review, recon_scf_seed
from kinbot.execution import atomic_json, file_sha256
from test_seeded_casci_panel import (
    combined, spectrum, continuation_tree, resolution_tree, resolution_plan, tree,
    followup_plan, plan, original_tree, snapshot,
)
from test_compact_scf_followup import simulate_campaign

PLAN = Path(__file__).resolve().parents[1]/'examples/anl/orbital_review_batches/r7_state_support_20261026'


def native_blocks():
    return (PLAN/'native_character_evidence.txt').read_text().split('\nNATIVE ')[1:]


@pytest.fixture
def support_tree(combined, tmp_path):
    production, old_output, old_plan = combined
    review.batch(production, old_plan/'screens.json', old_output, shared_pool=True)
    for system in review._saved(old_output)['systems']:
        path = Path(system['manifest']); state = review._read(path)
        child = path.parents[1].name
        for job in state['jobs']:
            directory = Path(job['directory']); out = directory/(job['name']+'.out')
            if not job['probe_active_space']:
                block = next(b for b in native_blocks() if f'/{child}/native/{job["name"]}/' in b.splitlines()[0])
                rows = re.findall(r'(?m)^SCF (\d+\.1) occ=(\S+) energy_Eh=(\S+)', block)
                out.write_text('PROGRAM * RHF-SCF\n ELECTRON ORBITALS\n\n'+''.join(
                    f' {i} {float(o):.5f} {float(e):.5f} 1 1 px 0.50000\n' for i,o,e in rows)+'\nEND\n')
            else:
                out.write_text('Synthetic fixed-orbital CASCI completion, not scientific evidence.\n')
            wfu = out.with_suffix('.wfu'); wfu.write_bytes(('mock:'+child+job['name']).encode())
            job.update(status='complete', wfu=str(wfu),wfu_sha256=file_sha256(wfu),output_sha256=file_sha256(out))
        state['status']='complete_review_required';atomic_json(path,state)
    for item in review._read(PLAN/'new_requests_evidence.json')['requests']:
        q = item['request']; w = production/'ch3ss'/q['parent']
        atomic_json(w/'vrctst/orbital_review'/q['reaction']/'request.json',q)
        atomic_json(w/'accepted_l2_network.json',{'network':{'reactions':item['network']}})
    return production,production/'orbital_recon/state_support_20261026'


@pytest.fixture
def support_plan(support_tree,tmp_path):
    production,_=support_tree;dest=tmp_path/'support_plan';dest.mkdir()
    for p in PLAN.glob('*.json'): (dest/p.name).write_bytes(p.read_bytes())
    for entry in review._read(dest/'screens.json')['screens']:
        if not entry.get('source'):continue
        jobs={j['name']:j for j in review._read(review._saved(production/entry['source'])['systems'][0]['manifest'])['jobs']}
        path=dest/entry['recipe'];recipe=review._read(path)
        for candidate in recipe['candidates']:
            for anchor in candidate['anchors']:
                anchor['seed_sha256']={b:jobs[b+'__'+anchor['seed']]['wfu_sha256'] for b in ('adz','atz')}
        atomic_json(path,recipe)
    return dest


def test_recorded_choices_have_native_hashes_somos_and_basis_specific_virtuals():
    assert file_sha256(PLAN/'native_character_evidence.txt')==review._read(PLAN/'evidence_provenance.json')['character_sha256']
    count=0
    for entry in review._read(PLAN/'screens.json')['screens']:
        recipe=review._read(PLAN/entry['recipe'])
        for c in recipe.get('candidates',[]):
            p=c['probe_definition'];count+=4
            for anchor in c['anchors']:
                for b in ('adz','atz'):
                    block=next(x for x in native_blocks() if f'/{entry["name"]}/' in x.splitlines()[0] and f'/{b}__{anchor["seed"]}/' in x.splitlines()[0])
                    assert re.search(r'(?m)^WFU .* SHA256 (\w+)',block)[1]==anchor['seed_sha256'][b]
                    if p.get('orbital_source')=='scf':
                        occ={int(i):float(o) for i,o in re.findall(r'(?m)^SCF (\d+)\.1 occ=(\S+)',block)}
                        indices=anchor['active_indices'][b]
                        assert sum(occ[i] for i in indices)==p['active_electrons']
                        assert all(i in indices for i,o in occ.items() if o==1)
    assert count==44
    ss=review._read(PLAN/'new_ss.json')['candidates'][0]['anchors'][0]['active_indices']
    assert ss=={'adz':[18,19,21,22],'atz':[18,19,21,23]}


def test_all_52_preflight_resume_geometry_state_support_and_immutable_sources(support_tree,support_plan,monkeypatch):
    production,output=support_tree;before=snapshot(production)
    monkeypatch.setattr(review.subprocess,'run',lambda *a,**k:pytest.fail('No scheduler in preflight/staging'))
    assert review.batch(production,support_plan/'screens.json',output,shared_pool=True,preflight=True)['tasks']==52
    assert snapshot(production)==before and not output.exists()
    assert review.batch(production,support_plan/'screens.json',output,shared_pool=True)['tasks']==52
    natural=0;scf=0;characters=0
    for system in review._saved(output)['systems']:
        cfg=review.recon._load(system['configuration']);space=cfg.get('_probe_space')
        manifest=review._read(system['manifest'])
        assert not cfg.get('correlation_checks')
        for job in manifest['jobs']:
            deck=(Path(job['directory'])/(job['name']+'.inp')).read_text()
            assert not job.get('job_id') and not job.get('orbital_predecessor')
            assert 'rs2' not in deck.lower()
            if not space:
                characters+=1;assert '{multi' not in deck
                assert job['distance_angstrom'] in (2.,30.)
                anchor=next(a for a in cfg['anchors'] if job['name'].endswith('__'+a['name']))
                assert sorted(anchor['moving_atoms']) == [2,4]
                xyz=np.array(job['geometry_angstrom']);original=np.array(cfg['geometry_angstrom'])
                assert np.allclose(xyz[[0,2,4,5]],original[[0,2,4,5]])
                assert np.allclose(xyz[[1,3]]-original[[1,3]],(xyz-original)[1])
                assert np.linalg.norm(xyz[0]-xyz[1])==pytest.approx(job['distance_angstrom'])
            elif space.get('orbital_partition'):
                natural+=1
                assert 'start,2140.2,type=natural' in deck and 'rotate,' not in deck
                assert ',ovl' not in deck.lower()
                assert (space['active_electrons'],space['active_orbitals']) in ((5,4),(9,6))
                assert 'closed,'+str((41-space['active_electrons'])//2) in deck
                assert job['root_states']==([{'multiplicity':2,'root':i} for i in range(1,5)]+[{'multiplicity':4,'root':i} for i in (1,2)])
                # Source and target anchors are identical, no molecular relaxation.
                old=review._read(Path(system['manifest']).parents[1]/'source_manifest.json')
                original=next(j for j in old['jobs'] if j['name']==job['name'])
                assert np.array_equal(job['geometry_angstrom'],original['geometry_angstrom'])
            else:
                scf+=1;assert space['active_orbitals']<=5
                recon_scf_seed.validate(job['initial_record'],41,space['active_electrons'],space['active_orbitals'])
                native=recon_scf_seed.separate(deck,'immutable_source.wfu')
                assert 'start,2100.3,type=canonical' in native and ',ovl' not in native
    assert (natural,scf,characters)==(8,36,8)
    after=snapshot(production)
    review.batch(production,support_plan/'screens.json',output,shared_pool=True,resume=True)
    assert snapshot(production)==after
    assert all(p.read_bytes()==v for p,v in before.items())


def test_bad_final_candidate_hash_prevents_any_staging(support_tree,support_plan):
    production,output=support_tree;p=support_plan/'new_hcssh_h_exit.json';recipe=review._read(p)
    recipe['candidates'][-1]['anchors'][-1]['seed_sha256']['atz']='0'*64;atomic_json(p,recipe)
    before=snapshot(production)
    with pytest.raises(ValueError,match='explicitly reviewed seed hash'):
        review.batch(production,support_plan/'screens.json',output,shared_pool=True)
    assert snapshot(production)==before and not output.exists()


def test_admission_retains_live_production_jobs_and_repeated_invocation(support_tree,support_plan,monkeypatch):
    simulate_campaign(support_tree,support_plan,monkeypatch)


def test_all_seven_live_requests_match_and_new_network_changes_are_refused(support_tree,support_plan):
    production,_=support_tree
    spec=importlib.util.spec_from_file_location('support_check',PLAN/'check_new_requests.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    before=snapshot(production)
    result=module.check(production,support_plan)
    assert len(result['requests'])==7 and snapshot(production)==before
    assert sum('orbital_screen_shared_with' in q for q in result['requests'].values())==1
    item=review._read(support_plan/'new_requests_evidence.json')['requests'][-1]
    path=production/'ch3ss'/item['request']['parent']/'accepted_l2_network.json'
    value=review._read(path);value['network']['reactions'][0]['products']=['changed']
    atomic_json(path,value)
    with pytest.raises(ValueError,match='Accepted reaction differs'):
        module.check(production,support_plan)
