"""Candidate publication uses an explicit model and reviewed artifact hashes."""
from copy import deepcopy
from pathlib import Path

import pytest

from kinbot import input_refresh
from kinbot.execution import atomic_json, file_sha256
from test_requested_orbital_recon import request_tree, completed_screen, helper, SELECTION


@pytest.fixture
def candidates(request_tree, monkeypatch):
    request, inp, source = request_tree
    completed_screen(request_tree)
    atomic_json(inp, dict(helper._read(inp), rotdpy_dist=[3., 4., 5.]))
    saved = helper._read(source/'review_launch.json')
    first = saved['systems'][0]
    manifest = helper._read(first['manifest'])
    second = dict(first, title='explicit_candidate', manifest=str(source/'selected.json'))
    for job in manifest['jobs']:
        d = source/'selected'/job['name']; d.mkdir(parents=True)
        wfu, out = d/(job['name']+'.wfu'), d/(job['name']+'.out')
        wfu.write_text('synthetic selected checkpoint '+job['name'])
        out.write_text('synthetic selected output')
        job.update(directory=str(d), wfu=str(wfu), wfu_sha256=file_sha256(wfu),
                   output_sha256=file_sha256(out))
    atomic_json(second['manifest'], manifest)
    old = helper._read(first['manifest']); old['jobs'][0]['status'] = 'failed'
    atomic_json(first['manifest'], old)
    saved['systems'].append(second); atomic_json(source/'review_launch.json', saved)
    choice = helper._read(SELECTION); choice['source_model'] = second['title']
    jobs = {j['name']: j for j in manifest['jobs']}
    for seed in choice['basis_seeds'].values():
        seed.update({k: jobs[seed['job']][k] for k in ('wfu_sha256', 'output_sha256')})
    recipe = source.parent/'selected_recipe.json'; atomic_json(recipe, choice)
    monkeypatch.setattr(input_refresh, 'owner_directory', lambda _: inp.parent)
    monkeypatch.setattr(helper.subprocess, 'run', lambda *a, **k: pytest.fail('No launch'))
    calls = []
    def submit(path, **kwargs):
        assert helper._read(path) == kwargs['expected']
        if kwargs.get('apply'):
            calls.append(kwargs['proposal'])
            atomic_json(path, kwargs['proposal'])
        return dict(status='queued', drivers_started=0)
    monkeypatch.setattr(input_refresh, 'submit', submit)
    plan = recipe.with_name('plan.json')
    atomic_json(plan, dict(schema=1, selections=[dict(name='selected', case=inp.parent.name,
        source=source.name, recipe=recipe.name)]))
    return inp, source, recipe, plan, second, calls


def test_explicit_second_candidate_preview_proposal_batch_and_repeat(candidates):
    inp, source, recipe, plan, system, calls = candidates
    before = {p: p.read_bytes() for p in source.rglob('*') if p.is_file()}
    result = helper.review(source, inp, recipe)
    assert result['source_model'] == 'explicit_candidate'
    assert all('/selected/' in b['path'] for b in result['basis_seeds'].values())
    out = source.parent/'proposal.json'
    helper.proposal(source, inp, recipe, out)
    assert system['manifest'] in helper._read(out)['evidence']
    assert str(source/'native/reconnaissance.json') not in helper._read(out)['evidence']
    result = helper.review_batch(inp.parent.parent, source.parent, plan, require_all=True, apply=True)
    assert result['cases'][inp.parent.name]['status'] == 'queued'
    assert len(calls) == 1
    helper.review_batch(inp.parent.parent, source.parent, plan, require_all=True, apply=True)
    assert calls[0] == calls[1]
    assert all(p.read_bytes() == v for p, v in before.items())


@pytest.mark.parametrize('failure', ['ambiguous', 'unknown', 'duplicate', 'changed_hash', 'validation', 'held'])
def test_invalid_candidate_never_publishes(candidates, failure):
    inp, source, recipe, plan, system, calls = candidates
    choice = helper._read(recipe)
    if failure == 'ambiguous': choice.pop('source_model')
    elif failure == 'unknown': choice['source_model'] = 'missing'
    elif failure == 'duplicate':
        saved = helper._read(source/'review_launch.json')
        saved['systems'].append(deepcopy(system)); atomic_json(source/'review_launch.json', saved)
    elif failure in ('changed_hash', 'validation'):
        manifest = helper._read(system['manifest'])
        job = next(j for j in manifest['jobs'] if j['name'] == choice['basis_seeds']['avdz-f12']['job'])
        if failure == 'changed_hash':
            Path(job['wfu']).write_text('new attempt with updated manifest')
            job['wfu_sha256'] = file_sha256(job['wfu'])
        else:
            choice['required_complete_jobs'] = []
            job['validation_errors'] = ['native evidence invalid']
        atomic_json(system['manifest'], manifest)
    else: choice['status'] = 'held'
    atomic_json(recipe, choice)
    before = inp.read_bytes()
    with pytest.raises(ValueError):
        helper.review_batch(inp.parent.parent, source.parent, plan, require_all=True, apply=True)
    assert inp.read_bytes() == before and not calls


def test_held_combined_plan_cannot_partially_publish(candidates):
    inp, source, _, plan, _, calls = candidates
    atomic_json(plan, dict(helper._read(plan), status='held', status_reason='one unresolved channel'))
    before = inp.read_bytes()
    with pytest.raises(ValueError, match='one unresolved channel'):
        helper.review_batch(inp.parent.parent, source.parent, plan, apply=True)
    assert inp.read_bytes() == before and not calls
