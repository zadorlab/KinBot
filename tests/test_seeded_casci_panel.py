"""Fixed-orbital root diagnostics preserve seeds and cannot publish a model."""
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from kinbot import requested_orbital_recon as review
from kinbot.execution import atomic_json, file_sha256
from test_compact_continuation import (
    continuation_tree, resolution_tree, resolution_plan, tree, followup_plan,
    plan, original_tree, snapshot,
)
from test_orbital_candidate_selection import candidates, request_tree, completed_screen, helper, SELECTION

PLAN = Path(__file__).resolve().parents[1]/'examples/anl/orbital_review_batches/r7_root_spectrum_20261026'


@pytest.fixture
def spectrum(continuation_tree, tmp_path):
    production, _, _ = continuation_tree
    dest = tmp_path/'spectrum_plan'; dest.mkdir()
    for p in PLAN.glob('*.json'):
        (dest/p.name).write_bytes(p.read_bytes())
    entry = review._read(dest/'screens.json')['screens'][0]
    saved = review._saved(production/entry['source'])
    jobs = review._read(saved['systems'][0]['manifest'])['jobs']
    path = dest/entry['recipe']; recipe = review._read(path)
    for anchor in recipe['candidates'][0]['anchors']:
        anchor['seed_sha256'] = {b: next(j['wfu_sha256'] for j in jobs
            if j['name'] == b+'__'+anchor['seed']) for b in ('adz', 'atz')}
    atomic_json(path, recipe)
    return production, production/'orbital_recon/root_spectrum_20261026', dest


def test_four_fixed_subspace_diagnostics_preserve_geometry_roots_and_resume(spectrum, monkeypatch):
    production, output, plan = spectrum
    monkeypatch.setattr(review.subprocess, 'run', lambda *a, **k: pytest.fail('No scheduler in staging'))
    before = snapshot(production)
    assert review.batch(production, plan/'screens.json', output, shared_pool=True, preflight=True)['tasks'] == 4
    assert snapshot(production) == before and not output.exists()
    assert review.batch(production, plan/'screens.json', output, shared_pool=True)['tasks'] == 4
    system = review._saved(output)['systems'][0]
    config = review.recon._load(system['configuration'])
    jobs = review._read(system['manifest'])['jobs']
    old = review._read(Path(system['manifest']).parents[1]/'source_manifest.json')
    for job in jobs:
        deck = (Path(job['directory'])/(job['name']+'.inp')).read_text()
        assert '{multi,cas-ci;' in deck and 'start,2140.2,type=natural' in deck
        assert 'noextra\n' in deck and 'closed,16' in deck and 'occ,22' in deck
        assert all(s not in deck for s in ('diab,', '{avas', 'freeze,', 'frozen,', 'accuracy,'))
        assert job['root_states'] == ([{'multiplicity':2,'root':r} for r in range(1,5)]
                                     + [{'multiplicity':4,'root':r} for r in (1,2)])
        assert job['orbital_predecessor'] is None and job['initial_record']['sha256'] in deck
        source = next(j for j in old['jobs'] if j['name'] == job['name'])
        assert np.array_equal(job['geometry_angstrom'], source['geometry_angstrom'])
    assert config['correlation_checks'] == [] and not config['casscf_options']
    after = snapshot(production)
    review.batch(production, plan/'screens.json', output, shared_pool=True, resume=True)
    assert snapshot(production) == after
    assert all(p.read_bytes() == data for p, data in before.items())


@pytest.mark.parametrize('change', ['geometry', 'seed', 'propagate', 'optimizer', 'correlation', 'method'])
def test_inapplicable_fixed_diagnostic_leaves_no_staged_jobs(spectrum, change):
    production, output, plan = spectrum
    path = plan/'csh2_sh_exit.json'; recipe = review._read(path); c = recipe['candidates'][0]
    if change == 'geometry': c['anchors'][-1]['distance_angstrom'] = 29.
    elif change == 'seed': c['anchors'][-1]['seed_sha256']['atz'] = '0'*64
    elif change == 'propagate': c['orbital_initialization'] = 'propagated'
    elif change == 'optimizer': c['casscf_options'] = {'quasi_newton': True}
    elif change == 'correlation': c['correlation_checks'] = [{'method':'caspt2-f12(9,6)'}]
    else: c['probe_method'] = 'unknown'
    atomic_json(path, recipe); before = snapshot(production)
    with pytest.raises(ValueError):
        review.batch(production, plan/'screens.json', output, shared_pool=True)
    assert not output.exists() and snapshot(production) == before


def test_completed_seeded_casci_cannot_be_published_by_another_recipe(candidates):
    inp, source, recipe, plan, system, calls = candidates
    manifest = review._read(system['manifest'])
    for job in manifest['jobs']:
        job.update(probe_method='casci', probe_active_space={'orbital_partition': {'close_lowest':0}})
    atomic_json(system['manifest'], manifest); before = inp.read_bytes()
    with pytest.raises(ValueError, match='diagnostic only'):
        review.review_batch(inp.parent.parent, source.parent, plan, require_all=True, apply=True)
    assert inp.read_bytes() == before and not calls


def test_native_diagnostic_reads_checkpoint_copy_and_keeps_source(spectrum, monkeypatch, tmp_path):
    production, output, plan = spectrum
    review.batch(production, plan/'screens.json', output, shared_pool=True)
    path = Path(review._saved(output)['systems'][0]['manifest'])
    job = review._read(path)['jobs'][0]
    seed = Path(job['initial_record']['path']); before = seed.read_bytes()
    scratch = tmp_path/'scratch'; scratch.mkdir()
    monkeypatch.setattr(review.recon, 'qc_runtime_environment', lambda *a, **k:
                        ({}, {'scratch_directory':str(scratch)}))
    monkeypatch.setattr(review.recon, 'cleanup_qc_runtime', lambda *a: None)
    def execute(command, **kwargs):
        import re
        local = kwargs['cwd']; deck = (local/command[-1]).read_text()
        filename = re.search(r'(?m)^file,2,([^,\n]+)$', deck)[1]
        assert (local/filename).read_bytes() == before
        assert '{multi,cas-ci;' in deck and 'noextra\n' in deck
        assert 'start,2140.2,type=natural' in deck and 'state,4\n' in deck
        (local/filename).write_bytes(b'synthetic CASCI natural export')
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(review.recon.subprocess, 'run', execute)
    review.recon.execute(path, job['name'])
    assert seed.read_bytes() == before
    x = review._read(Path(job['directory'])/'execution.json')
    assert x['native_input_sha256'] == file_sha256(x['native_input'])


def test_diagnostic_scheduling_keeps_live_work_and_does_not_resubmit(spectrum, monkeypatch):
    # Exercise the same shared-cap runner with independent diagnostic jobs.
    from test_natural_union_panel import test_four_job_schedule_checks_predecessors_and_preserves_22_live_jobs
    test_four_job_schedule_checks_predecessors_and_preserves_22_live_jobs(spectrum, monkeypatch)


@pytest.fixture
def combined(spectrum):
    from ase import Atoms
    from ase.db import connect
    production,output,plan=spectrum
    rows=[]
    for item in review._read(plan/'new_requests_evidence.json')['requests']:
        q=item['request'];w=production/'ch3ss'/q['parent']
        atomic_json(w/'vrctst/orbital_review'/q['reaction']/'request.json',q)
        rows.append(item['network'])
        if item['geometry']:
            g=item['geometry']
            connect(w/'kinbot.db').write(Atoms(g['symbols'],positions=g['positions']),
                                        name=g['name'],data={'status':g['status']})
    atomic_json(w/'accepted_l2_network.json',{'network':{'reactions':rows}})
    atomic_json(plan/'screens.json',review._read(plan/'combined_screens.json'))
    return production,output,plan


def alias_checker():
    import importlib.util
    spec=importlib.util.spec_from_file_location('check_new_requests',PLAN/'check_new_requests.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def test_twenty_job_panel_checks_five_requests_and_preserves_all_fragment_geometries(combined,monkeypatch):
    production,output,plan=combined
    monkeypatch.setattr(review.subprocess,'run',lambda *a,**k:pytest.fail('No scheduler in staging'))
    before=snapshot(production)
    checks=alias_checker().check(production,plan)
    assert len(checks['requests'])==5
    assert sum('orbital_screen_shared_with' in r for r in checks['requests'].values())==1
    assert review.batch(production,plan/'screens.json',output,shared_pool=True,preflight=True)['tasks']==20
    assert before==snapshot(production) and not output.exists()
    review.batch(production,plan/'screens.json',output,shared_pool=True)
    characters=0
    for system in review._saved(output)['systems']:
        config=review.recon._load(system['configuration'])
        if not config.get('orbital_character_only'):continue
        manifest=review._read(system['manifest'])
        assert len(manifest['jobs'])==4
        for j in manifest['jobs']:
            characters+=1
            deck=(Path(j['directory'])/(j['name']+'.inp')).read_text()
            assert '{rhf' in deck and '{multi' not in deck and '{avas' not in deck and 'rs2' not in deck
            assert j['probe_active_space'] is None
        for anchor in config['anchors']:
            xyz=np.array(review.recon._anchor_geometry(config,anchor))
            original=np.array(config['geometry_angstrom'])
            i,j=np.array(anchor['bond'])-1
            assert np.linalg.norm(xyz[i]-xyz[j])==pytest.approx(anchor['distance_angstrom'])
            moving=np.array(anchor['moving_atoms'])-1
            fixed=np.array([i for i in range(6) if i not in moving])
            assert np.allclose(xyz[fixed],original[fixed])
            assert np.allclose(xyz[moving]-original[moving],(xyz[moving]-original[moving])[0])
    assert characters==16
    after=snapshot(production)
    review.batch(production,plan/'screens.json',output,shared_pool=True,resume=True)
    assert snapshot(production)==after
    assert all(p.read_bytes()==v for p,v in before.items())


def test_shared_shift_endpoint_refuses_changed_live_geometry(combined):
    from ase import Atoms
    from ase.db import connect
    production,output,plan=combined
    item=review._read(plan/'new_requests_evidence.json')['requests'][-1]
    q=item['request'];g=item['geometry']
    connect(production/'ch3ss'/q['parent']/'kinbot.db').write(
        Atoms(g['symbols'],positions=np.array(g['positions'])+.01),
        name=g['name'],data={'status':'normal'})
    before=snapshot(production)
    with pytest.raises(ValueError,match='endpoint differs'):
        alias_checker().check(production,plan)
    assert snapshot(production)==before and not output.exists()


def test_combined_panel_preserves_22_live_allocations_and_repeated_run(combined,monkeypatch):
    from test_remaining_orbital_screens import test_campaign_launch_and_resume_share_retained_production_capacity
    production,output,plan=combined
    test_campaign_launch_and_resume_share_retained_production_capacity((production,output),plan,monkeypatch)
