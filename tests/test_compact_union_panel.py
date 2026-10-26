"""Ten explicit diagnostic jobs; simulated transport/site checks, not native QC."""
from pathlib import Path

import numpy as np
import pytest

from kinbot import requested_orbital_recon as review, recon_scf_seed
from kinbot.execution import atomic_json
from test_compact_resolution_panel import (
    resolution_tree, tree, followup_plan, plan, original_tree, snapshot, simulate_campaign,
)

PLAN = Path(__file__).resolve().parents[1]/'examples/anl/orbital_review_batches/r7_compact_union_20261025'


@pytest.fixture
def union_tree(resolution_tree):
    production, _ = resolution_tree
    return production, production/'orbital_recon/compact_union_20261025'


@pytest.fixture
def union_plan(union_tree, tmp_path):
    production, _ = union_tree
    dest = tmp_path/'union_plan'; dest.mkdir()
    for p in PLAN.glob('*.json'):
        (dest/p.name).write_bytes(p.read_bytes())
    for entry in review._read(dest/'screens.json')['screens']:
        saved = review._saved(production/entry['source'])
        jobs = {j['name']: j for j in review._read(saved['systems'][0]['manifest'])['jobs']}
        path = dest/entry['recipe']; recipe = review._read(path)
        for anchor in recipe['candidates'][0]['anchors']:
            anchor['seed_sha256'] = {b: jobs[b+'__'+anchor['seed']]['wfu_sha256']
                                     for b in ('adz', 'atz')}
        atomic_json(path, recipe)
    return dest


def test_ten_job_preflight_geometry_state_weights_and_resume(union_tree, union_plan, monkeypatch):
    production, output = union_tree
    before = snapshot(production)
    monkeypatch.setattr(review.subprocess, 'run', lambda *a, **k: pytest.fail('No scheduler during staging'))
    assert review.batch(production, union_plan/'screens.json', output,
                        shared_pool=True, preflight=True)['tasks'] == 10
    assert not output.exists() and snapshot(production) == before
    assert review.batch(production, union_plan/'screens.json', output, shared_pool=True)['tasks'] == 10
    for system in review._saved(output)['systems']:
        config = review.recon._load(system['configuration'])
        space = config['_probe_space']
        manifest = review._read(system['manifest'])
        cs = system['title'].startswith('well_')
        assert (space['active_electrons'], space['active_orbitals']) == ((5,4) if cs else (9,6))
        assert len(manifest['jobs']) == (6 if cs else 4)
        assert not config['correlation_checks']
        source = Path(review._saved(Path(system['manifest']).parents[1])['source_review'])
        old = review.recon._load(review._saved(source)['systems'][0]['configuration'])
        for a in (config['anchors'][0], config['anchors'][-1]):
            original = next(x for x in old['anchors'] if x['name'] == a['name'])
            assert np.array_equal(review.recon._anchor_geometry(config,a),
                                  review.recon._anchor_geometry(old,original))
        for job in manifest['jobs']:
            recon_scf_seed.validate(job['initial_record'], 41, space['active_electrons'], space['active_orbitals'])
            deck = (Path(job['directory'])/(job['name']+'.inp')).read_text()
            native = recon_scf_seed.separate(deck, 'unchanged_source.wfu')
            assert ',ovl' not in native and 'start,2100.3,type=canonical' in native
            assert f'closed,{18 if cs else 16}' in native
            assert f'occ,{22}' in native
            assert ('wf,41,1,3' in native) == cs
            assert ('state,2' in native) != cs
            assert not job.get('orbital_predecessor') and not job.get('job_id')
    after = snapshot(production)
    review.batch(production, union_plan/'screens.json', output, shared_pool=True, resume=True)
    assert snapshot(production) == after
    assert all(p.read_bytes() == v for p,v in before.items())


def test_bad_last_native_hash_blocks_entire_plan(union_tree, union_plan):
    production, output = union_tree
    path = union_plan/'csh2_sh_exit_resolution.json'; recipe = review._read(path)
    recipe['candidates'][0]['anchors'][-1]['seed_sha256']['atz'] = '0'*64
    atomic_json(path, recipe)
    before = snapshot(production)
    with pytest.raises(ValueError, match='explicitly reviewed seed hash'):
        review.batch(production, union_plan/'screens.json', output, shared_pool=True)
    assert snapshot(production) == before and not output.exists()


def test_ten_admissions_preserve_22_live_jobs_policy_and_repeat(union_tree, union_plan, monkeypatch):
    simulate_campaign(union_tree, union_plan, monkeypatch)
