"""Eleven prepared model recipes; synthetic artifacts, never a chemistry verdict."""
from copy import deepcopy
from pathlib import Path

import pytest

from kinbot import input_refresh, requested_orbital_recon as review
from kinbot.execution import atomic_json, file_sha256
from kinbot.molpro import _render_vrc_state_model

PLAN = Path(__file__).resolve().parents[1]/'examples/anl/orbital_review_batches/r7_compact_union_20261025'


def test_eleven_explicit_recipes_and_combined_hold(tmp_path, monkeypatch):
    controls = review._read(PLAN/'publication_plan.json')
    assert controls['status'] == 'held' and len(controls['selections']) == 11
    root = tmp_path/'production'; case = root/'ch3ss'; case.mkdir(parents=True)
    inp = case/'input.json'
    previous = dict(rotdpy_dist=[3., 4., 7.], vrc_tst_scan_points=[2.3, 3.7, 6.1],
                    vrc_tst_fragment_geometry='asymptotic_scan', unrelated='retained')
    atomic_json(inp, previous)
    campaign = root/'orbital_recon'; campaign.mkdir()
    plan_dir = tmp_path/'plans'; plan_dir.mkdir()
    monkeypatch.setattr(input_refresh, 'owner_directory', lambda _: case)
    monkeypatch.setattr(review.subprocess, 'run', lambda *a, **k: pytest.fail('No scheduler'))
    published = []
    def submit(path, **kw):
        assert kw['expected'] == review._read(path)
        if kw.get('apply'):
            published.append(kw['proposal']); atomic_json(path, kw['proposal'])
        return dict(status='queued', drivers_started=0)
    monkeypatch.setattr(input_refresh, 'submit', submit)
    for item in controls['selections']:
        choice = review._read(PLAN/item['recipe']); match = choice['request_match']
        source = campaign/item['source']; source.mkdir(parents=True)
        xyz = [[float(i),0.,0.] for i in range(6)]
        request = case/match['parent']/'vrctst/orbital_review'/match['reaction']/'request.json'
        atomic_json(request, dict(match, positions=xyz, schema=1, status='needs_review'))
        atomic_json(source/'request.json', review._read(request))
        config = source/'configuration.json'; atomic_json(config, dict(geometry_angstrom=xyz))
        jobs = []
        for name in choice['required_complete_jobs']:
            directory = source/name; directory.mkdir()
            wfu, out = directory/(name+'.wfu'), directory/(name+'.out')
            wfu.write_text('synthetic '+name); out.write_text('synthetic native result')
            jobs.append(dict(name=name, directory=str(directory), status='complete',
                natural_orbital_probe=True, basis='avdz-f12' if name.startswith('adz') else 'avtz-f12',
                wfu=str(wfu), wfu_sha256=file_sha256(wfu), output_sha256=file_sha256(out)))
        for seed in choice['basis_seeds'].values():
            job = next(j for j in jobs if j['name'] == seed['job'])
            seed.update({k: job[k] for k in ('wfu_sha256', 'output_sha256')})
        manifest = source/'manifest.json'; atomic_json(manifest, dict(jobs=jobs, electron_count=41))
        atomic_json(source/'review_launch.json', dict(files={}, request=str(request), systems=[
            dict(title=choice['source_model'], manifest=str(manifest), configuration=str(config))]))
        recipe = plan_dir/item['recipe']; recipe.parent.mkdir(exist_ok=True)
        atomic_json(recipe, choice)
        item.pop('supporting_reviews', None)  # Supporting evidence has separate hash-failure tests.
    plan = plan_dir/'plan.json'; atomic_json(plan, controls)
    with pytest.raises(ValueError, match='not ready'):
        review.review_batch(root, campaign, plan, apply=True)
    assert not published and review._read(inp) == previous
    # Deliberately exercise only the prepared eleven in a synthetic test. The
    # real twelve-channel release remains held until the final chemistry review.
    controls['status'] = 'ready'; atomic_json(plan, controls)
    result = review.review_batch(root, campaign, plan, require_all=True, apply=True)
    assert len(result['cases']['ch3ss']['channels']) == 11 and len(published) == 1
    installed = review._read(inp)
    assert {k: installed[k] for k in previous} == previous
    channels = [c for m in installed['vrc_tst_surface_models'].values() for c in m['channels']]
    assert len(channels) == 11
    for channel in channels:
        stage = channel['bottlenecks'][0]; p = stage['parameters']; active = p['vrc_tst_active_space']
        ne, no = active['active_electrons'], active['active_orbitals']
        assert stage['coordinate'] == ('reactive' if '_hom_sci_' in channel['reaction'] else 'com')
        assert stage['distances_angstrom'] == previous['rotdpy_dist']
        assert len(stage['states']) == 1 and stage['states'][0]['root'] == 1
        for model in p['vrc_tst_state_model']['basis_models'].values():
            _, deck, _ = _render_vrc_state_model(model, nelectron=41,
                active_electrons=ne, active_orbitals=no, spin=1)
            assert 'start,2140.2,type=natural' in deck
            assert f'closed,{(41-ne)//2}' in deck
    review.review_batch(root, campaign, plan, require_all=True, apply=True)
    assert published[0] == published[1]
