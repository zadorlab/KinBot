"""Recorded orbital choices and simulated transport/admission, not native QC."""
from pathlib import Path
import re

import numpy as np
import pytest

from kinbot import requested_orbital_recon as review, recon_scf_seed
from kinbot.execution import atomic_json, file_sha256
from test_compact_scf_followup import (
    PLAN as PREVIOUS, tree, followup_plan, plan, original_tree, snapshot,
    simulate_campaign,
)

PLAN = PREVIOUS.parent/'r7_compact_resolution_20261025'


def native_blocks():
    return (PLAN/'native_character_evidence.txt').read_text().split('\nNATIVE ')[1:]


@pytest.fixture
def resolution_tree(tree, followup_plan):
    production, old_output = tree
    review.batch(production, followup_plan/'screens.json', old_output, shared_pool=True)
    # Complete only the source SCF screens. These checkpoints are mock bytes;
    # occupation records come from the user's native report for transport tests.
    for system in review._saved(old_output)['systems']:
        loaded = review.recon._load(system['configuration'])
        if loaded.get('_probe_space'):
            continue
        path = Path(system['manifest']); m = review._read(path)
        child = path.parents[1].name
        for job in m['jobs']:
            block = next(b for b in native_blocks()
                         if f'/{child}/native/{job["name"]}/' in b.splitlines()[0])
            rows = re.findall(r'(?m)^SCF (\d+\.1) occ=(\S+) energy_Eh=(\S+)', block)
            out = Path(job['directory'])/(job['name']+'.out')
            out.write_text('PROGRAM * RHF-SCF\n ELECTRON ORBITALS\n\n'+''.join(
                f' {i} {float(o):.5f} {float(e):.5f} 1 1 px 0.50000\n'
                for i,o,e in rows)+'\nEND\n')
            wfu = out.with_suffix('.wfu'); wfu.write_bytes(('mock:'+child+job['name']).encode())
            job.update(status='complete', wfu=str(wfu), wfu_sha256=file_sha256(wfu),
                       output_sha256=file_sha256(out))
        m['status'] = 'complete_review_required'; atomic_json(path, m)
    return production, production/'orbital_recon/compact_resolution_20261025'


@pytest.fixture
def resolution_plan(resolution_tree, tmp_path):
    production, _ = resolution_tree
    dest = tmp_path/'resolution_plan'; dest.mkdir()
    for p in PLAN.glob('*.json'):
        (dest/p.name).write_bytes(p.read_bytes())
    for entry in review._read(dest/'screens.json')['screens']:
        source = production/entry['source']; rec = review._saved(source)
        jobs = {j['name']:j for j in review._read(rec['systems'][0]['manifest'])['jobs']}
        path = dest/entry['recipe']; recipe = review._read(path)
        for c in recipe['candidates']:
            for a in c['anchors']:
                a['seed_sha256'] = {basis:jobs[basis+'__'+a['seed']]['wfu_sha256']
                                    for basis in ('adz','atz')}
        atomic_json(path, recipe)
    return dest


def test_new_seed_choices_bound_to_native_report_and_somos():
    assert file_sha256(PLAN/'native_character_evidence.txt') == review._read(PLAN/'evidence_provenance.json')['sha256']
    count = 0
    for entry in review._read(PLAN/'screens.json')['screens']:
        recipe = review._read(PLAN/entry['recipe'])
        if not entry['name'].startswith('well_'):
            continue
        for c in recipe['candidates']:
            p = c['probe_definition']; count += 4
            for anchor in c['anchors']:
                for basis in ('adz','atz'):
                    stem = basis+'__'+anchor['seed']
                    block = next(b for b in native_blocks() if f'/{entry["name"]}/native/{stem}/' in b.splitlines()[0])
                    assert re.search(r'(?m)^WFU .* SHA256 (\w+)', block)[1] == anchor['seed_sha256'][basis]
                    occ = {int(i):float(o) for i,o in re.findall(r'(?m)^SCF (\d+)\.1 occ=(\S+)', block)}
                    indices = anchor['active_indices'][basis]
                    assert sum(occ[i] for i in indices) == p['active_electrons']
                    assert all(i in indices for i,o in occ.items() if o == 1)
    assert count == 52


def test_all_60_preflight_stage_resume_preserve_geometry_states_and_existing_work(
        resolution_tree, resolution_plan, monkeypatch):
    production, output = resolution_tree
    monkeypatch.setattr(review.subprocess, 'run', lambda *a,**k: pytest.fail('Staging must not call scheduler'))
    before = snapshot(production)
    result = review.batch(production, resolution_plan/'screens.json', output, shared_pool=True, preflight=True)
    assert (result['screens'], result['tasks']) == (10,60)
    assert snapshot(production) == before and not output.exists()
    result = review.batch(production, resolution_plan/'screens.json', output, shared_pool=True)
    assert result['max_jobs'] == 3 and result['tasks'] == 60
    assert all(p.read_bytes() == data for p,data in before.items())
    source_entries = {e['name']:e for e in review._read(resolution_plan/'screens.json')['screens']}
    total = 0
    for child in review._saved(output)['reviews']:
        child = Path(child); entry = source_entries[child.name]
        source = production/entry['source']; source_record = review._saved(source)
        source_config = review.recon._load(source_record['systems'][0]['configuration'])
        if (source/'geometry_source.json').exists():
            assert (child/'geometry_source.json').read_bytes() == (source/'geometry_source.json').read_bytes()
        for system in review._saved(child)['systems']:
            config = review.recon._load(system['configuration'])
            manifest = review._read(system['manifest']); space = config['_probe_space']
            assert space['active_orbitals'] <= 5 and space['active_electrons'] <= 7
            assert manifest['shared_policy_file'] == str(production/'.kinbot_resource_policy.json')
            for anchor, original in zip(config['anchors'], source_config['anchors']):
                assert anchor['name'] == original['name']
                assert np.allclose(review.recon._anchor_geometry(config,anchor),
                                   review.recon._anchor_geometry(source_config,original))
            for job in manifest['jobs']:
                total += 1
                assert not job.get('job_id') and not job.get('orbital_predecessor')
                seed = job['initial_record']
                recon_scf_seed.validate(seed,41,space['active_electrons'],space['active_orbitals'])
                deck = (Path(job['directory'])/(job['name']+'.inp')).read_text()
                assert 'rs2' not in deck.lower() and 'diab,' not in deck.lower()
                # Use the actual executor transform, rather than assert only the template.
                transformed = recon_scf_seed.separate(deck, 'immutable_source.wfu')
                assert ',ovl' not in transformed and 'start,2100.3,type=canonical' in transformed
                assert f'closed,{(41-space["active_electrons"])//2}' in transformed
    assert total == 60
    after = snapshot(production)
    assert review.batch(production,resolution_plan/'screens.json',output,shared_pool=True,resume=True)['tasks'] == 60
    assert snapshot(production) == after


def test_bad_last_seed_prevents_partial_panel(resolution_tree, resolution_plan):
    production, output = resolution_tree
    entry = review._read(resolution_plan/'screens.json')['screens'][-1]
    recipe_path = resolution_plan/entry['recipe']; recipe = review._read(recipe_path)
    recipe['candidates'][0]['anchors'][-1]['seed_sha256']['atz'] = '0'*64
    atomic_json(recipe_path,recipe)
    before = snapshot(production)
    with pytest.raises(ValueError,match='explicitly reviewed seed hash'):
        review.batch(production,resolution_plan/'screens.json',output,shared_pool=True)
    assert not output.exists() and snapshot(production) == before


def test_all_60_admissions_retain_22_production_jobs_and_repeat_without_duplicates(
        resolution_tree, resolution_plan, monkeypatch):
    simulate_campaign(resolution_tree,resolution_plan,monkeypatch)
