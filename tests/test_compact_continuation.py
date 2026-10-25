"""Compact continuation provenance and orchestration; no native QC is simulated as science."""
from pathlib import Path
from types import SimpleNamespace
import json

import numpy as np
import pytest

from kinbot import requested_orbital_recon as review, run_policy, run_walltime
from kinbot.execution import atomic_json, file_sha256
from test_compact_resolution_panel import (
    resolution_tree, resolution_plan, tree, followup_plan, plan, original_tree, snapshot,
)

PLAN = Path(__file__).resolve().parents[1]/'examples/anl/orbital_review_batches/r7_compact_continuation_20261025'


@pytest.fixture
def continuation_tree(resolution_tree, resolution_plan, tmp_path):
    production, prior = resolution_tree
    review.batch(production, resolution_plan/'screens.json', prior, shared_pool=True)
    for system in review._saved(prior)['systems']:
        path = Path(system['manifest']); state = review._read(path)
        for job in state['jobs']:
            directory = Path(job['directory'])
            out = directory/(job['name']+'.out'); out.write_text('Synthetic CASSCF output.\n')
            wfu = out.with_suffix('.wfu'); wfu.write_bytes(str(directory).encode())
            job.update(status='complete', output_sha256=file_sha256(out),
                       wfu=str(wfu), wfu_sha256=file_sha256(wfu))
        atomic_json(path, state)
    prepared_plan = tmp_path/'continuation_plan'; prepared_plan.mkdir()
    for p in PLAN.glob('*.json'):
        (prepared_plan/p.name).write_bytes(p.read_bytes())
    for entry in review._read(prepared_plan/'screens.json')['screens']:
        record = review._saved(production/entry['source'])
        jobs = {j['name']:j for j in review._read(record['systems'][0]['manifest'])['jobs']}
        recipe_path = prepared_plan/entry['recipe']; recipe = review._read(recipe_path)
        for candidate in recipe['candidates']:
            for anchor in candidate['anchors']:
                if 'seed' in anchor:
                    anchor['seed_sha256'] = {
                        basis:jobs[basis+'__'+anchor['seed']]['wfu_sha256'] for basis in ('adz','atz')}
        atomic_json(recipe_path, recipe)
    return production, production/'orbital_recon/compact_continuation_20261025', prepared_plan


def test_recorded_seed_hashes_and_unchanged_compact_models():
    old = PLAN.parent/'r7_compact_resolution_20261025'
    evidence = (old/'resolution_character_evidence.txt').read_text()
    assert file_sha256(old/'resolution_character_evidence.txt') == review._read(
        old/'resolution_evidence_provenance.json')['sha256']
    for entry in review._read(PLAN/'screens.json')['screens']:
        recipe = review._read(PLAN/entry['recipe'])
        candidate = recipe['candidates'][0]
        previous = review._read(old/entry['recipe'])['candidates'][0]
        for key in ('active_electrons', 'active_orbitals', 'state_model'):
            assert candidate['probe_definition'][key] == previous['probe_definition'][key]
        assert candidate['casscf_options'] == previous['casscf_options']
        assert candidate['correlation_checks'] == []
        for basis, digest in candidate['anchors'][0]['seed_sha256'].items():
            block = next(b for b in evidence.split('\nNATIVE ')[1:]
                         if f'/{entry["name"]}/' in b.splitlines()[0]
                         and f'/{basis}__{candidate["anchors"][0]["seed"]}/' in b.splitlines()[0])
            assert 'SHA256 '+digest in block


def test_all_panel_preflight_order_hashes_geometry_and_repeat(continuation_tree, monkeypatch):
    production, output, plan = continuation_tree
    before = snapshot(production)
    monkeypatch.setattr(review.subprocess, 'run', lambda *a, **k: pytest.fail('No scheduler during staging'))
    assert review.batch(production, plan/'screens.json', output, shared_pool=True, preflight=True)['tasks'] == 14
    assert snapshot(production) == before and not output.exists()
    assert review.batch(production, plan/'screens.json', output, shared_pool=True)['tasks'] == 14
    for system in review._saved(output)['systems']:
        config = review.recon._load(system['configuration'])
        manifest = review._read(system['manifest'])
        assert config['orbital_initialization'] == 'propagated'
        assert len(config['initial_records']) == 2
        for basis in ('avdz-f12', 'avtz-f12'):
            jobs = [j for j in manifest['jobs'] if j['basis'] == basis]
            for index, job in enumerate(jobs):
                deck = (Path(job['directory'])/(job['name']+'.inp')).read_text()
                assert 'start,2140.2,type=natural' in deck
                assert 'rotate,' not in deck and ',ovl' not in deck and '{avas' not in deck
                if index == 0:
                    assert job['orbital_predecessor'] is None
                    assert job['initial_record']['sha256'] in deck
                    assert 'diab,' not in deck
                else:
                    assert job['orbital_predecessor'] == jobs[index-1]['name']
                    assert job['initial_record'] is None
                    assert 'diab,2140.2' in deck
            assert [j['distance_angstrom'] for j in jobs[:-1]] == (
                [6., 3.5, 2.5] if len(jobs) == 4 else [6., 3.5])
        # Final near geometry is the original independent target, not a newly relaxed/forced exit.
        child = Path(system['manifest']).parents[1]
        source = Path(review._saved(child)['source_review'])
        prior_config = review.recon._load(review._saved(source)['systems'][0]['configuration'])
        original_near = prior_config['anchors'][0]
        assert np.allclose(review.recon._anchor_geometry(config, config['anchors'][-1]),
                           review.recon._anchor_geometry(prior_config, original_near), atol=1e-12)
    after = snapshot(output)
    assert review.batch(production, plan/'screens.json', output, shared_pool=True, resume=True)['tasks'] == 14
    assert snapshot(output) == after
    assert all(p.read_bytes() == value for p, value in before.items())


@pytest.mark.parametrize('change,match', [
    ('first_propagates', 'preceding seed'), ('seed_on_successor', 'consume its predecessor'),
    ('bad_hash', 'explicitly reviewed seed hash'), ('new_chain', 'preceding seed'),
    ('reseed', 'one seed'), ('changed_coordinate', 'same coordinate'),
])
def test_invalid_chains_fail_before_any_staging(continuation_tree, change, match):
    production, output, plan = continuation_tree
    recipe_path = plan/'well_hom_sci_1_2.json'; recipe = review._read(recipe_path)
    anchors = recipe['candidates'][0]['anchors']
    if change == 'first_propagates':
        anchors[0].pop('seed'); anchors[0].pop('seed_sha256'); anchors[0]['propagate'] = True
    elif change == 'seed_on_successor':
        anchors[1]['seed'] = anchors[0]['seed']
    elif change == 'bad_hash':
        anchors[0]['seed_sha256']['atz'] = '0'*64
    elif change == 'new_chain':
        anchors[1]['orbital_chain'] = 'different'
    elif change == 'reseed':
        anchors[1].pop('propagate'); anchors[1]['seed'] = anchors[0]['seed']
    else:
        # A second coordinate in the source configuration must not silently join this chain.
        anchors[1]['source'] = 'cs_r_2'
        original = review.recon._load
        def changed(path):
            result = original(path)
            if str(path).endswith('well_hom_sci_1_2_radicals.json'):
                result['anchors'][0]['moving_atoms'] = [2]
            return result
        from unittest.mock import patch
        with patch.object(review.recon, '_load', changed):
            atomic_json(recipe_path, recipe)
            with pytest.raises(ValueError, match=match):
                review.batch(production, plan/'screens.json', output, shared_pool=True)
        assert not output.exists()
        return
    atomic_json(recipe_path, recipe)
    before = snapshot(production)
    with pytest.raises(ValueError, match=match):
        review.batch(production, plan/'screens.json', output, shared_pool=True)
    assert not output.exists() and snapshot(production) == before


def test_scheduler_waits_for_verified_predecessors_retains_live_work_and_repeats(
        continuation_tree, monkeypatch):
    production, output, plan = continuation_tree
    policy_path = production/'.kinbot_resource_policy.json'
    policy = review._read(policy_path); policy['fair_submission'] = False; atomic_json(policy_path, policy)
    review.batch(production, plan/'screens.json', output, shared_pool=True)
    policy_bytes = policy_path.read_bytes()
    running = {str(100+i):1 for i in range(23)}
    done, launched = set(), []
    monkeypatch.setattr(run_policy, '_queue', lambda root: running.copy())
    monkeypatch.setattr(run_policy, 'slurm_job_record', lambda j:
                        {'JobId':j, 'JobState':'COMPLETED'} if j in done else {})
    monkeypatch.setattr(run_walltime, 'partition_walltime', lambda p: '31-00:00:00')
    def scheduler(command, **kwargs):
        if command[0] == 'bash':
            return SimpleNamespace(returncode=0, stdout='/mock/molpro\n', stderr='')
        if command[0] == 'squeue':
            return SimpleNamespace(returncode=0, stdout='\n'.join(running), stderr='')
        assert command[0] == 'sbatch' and len(running) < 23
        directory = Path(kwargs['cwd']); manifest = review._read(directory.parent/'reconnaissance.json')
        job = next(j for j in manifest['jobs'] if j['directory'] == str(directory))
        if job['orbital_predecessor']:
            prior = next(j for j in manifest['jobs'] if j['name'] == job['orbital_predecessor'])
            assert prior['job_id'] in done and prior['status'] == 'complete'
            assert job['initial_wfu_sha256'] == file_sha256(prior['wfu'])
            assert (directory/(job['name']+'.wfu')).read_bytes() == Path(prior['wfu']).read_bytes()
            assert '{{PREVIOUS_WFU_SHA256}}' not in (directory/(job['name']+'.inp')).read_text()
        ident = str(1000+len(launched)); launched.append((ident, directory)); running[ident] = 1
        return SimpleNamespace(returncode=0, stdout=ident+'\n', stderr='')
    monkeypatch.setattr(review.subprocess, 'run', scheduler)
    def audit(path, **kwargs):
        state = review._read(path)
        for job in state['jobs']:
            job['status'] = 'complete' if job.get('job_id') in done else 'incomplete'
            if job['status'] == 'complete':
                wfu = Path(job['directory'])/(job['name']+'.wfu')
                wfu.write_bytes(('completed:'+job['name']).encode())
                job.update(wfu=str(wfu), wfu_sha256=file_sha256(wfu))
        state['status'] = ('complete_review_required' if all(j['status'] == 'complete' for j in state['jobs'])
                           else 'incomplete')
        atomic_json(path, state); return state
    monkeypatch.setattr(review.recon, '_audit', audit)
    monkeypatch.setattr(review.recon, '_recover_slurm_startup_failures', lambda *a: None)
    monkeypatch.setattr(review.recon, '_retry_allocation_failures', lambda p,m,*a: m)
    def finish(_):
        for ident in ['100', *[i for i,d in launched]]:
            if ident in running:
                running.pop(ident); done.add(ident)
    monkeypatch.setattr(review.recon.time, 'sleep', finish)
    assert review.run(output, report_failures=True)['status'] == 'complete_review_required'
    assert len(launched) == len({str(d) for i,d in launched}) == 14
    assert running == {str(100+i):1 for i in range(1,23)}
    assert policy_path.read_bytes() == policy_bytes
    assert review._read(production/'.kinbot_policy_admissions.json') == {}
    review.run(output, report_failures=True)
    assert len(launched) == 14
