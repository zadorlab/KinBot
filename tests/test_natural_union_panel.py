"""Explicit natural-space promotion and a four-job campaign, not native QC."""
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from kinbot import requested_orbital_recon as review, run_policy, run_walltime
from kinbot.execution import atomic_json, file_sha256
from test_compact_continuation import (
    continuation_tree, resolution_tree, resolution_plan, tree, followup_plan,
    plan, original_tree, snapshot,
)
from test_molpro_reconnaissance import _partition_chain_configuration

PLAN = Path(__file__).resolve().parents[1]/'examples/anl/orbital_review_batches/r7_natural_union_20261025'


@pytest.fixture
def natural_union(continuation_tree, tmp_path):
    production, _, _ = continuation_tree
    dest = tmp_path/'natural_plan'; dest.mkdir()
    for p in PLAN.glob('*.json'):
        (dest/p.name).write_bytes(p.read_bytes())
    entry = review._read(dest/'screens.json')['screens'][0]
    saved = review._saved(production/entry['source'])
    jobs = review._read(saved['systems'][0]['manifest'])['jobs']
    path = dest/entry['recipe']; recipe = review._read(path)
    recipe['candidates'][0]['anchors'][0]['seed_sha256'] = {
        b:next(j['wfu_sha256'] for j in jobs if j['name'] == b+'__exit_l2')
        for b in ('adz','atz')}
    atomic_json(path, recipe)
    return production, production/'orbital_recon/natural_union_20261025', dest


def test_four_jobs_select_actual_natural_union_once_and_preserve_geometry(natural_union, monkeypatch):
    production, output, plan = natural_union
    monkeypatch.setattr(review.subprocess, 'run', lambda *a, **k: pytest.fail('No scheduler in staging'))
    before = snapshot(production)
    assert review.batch(production, plan/'screens.json', output, shared_pool=True, preflight=True)['tasks'] == 4
    assert not output.exists() and snapshot(production) == before
    assert review.batch(production, plan/'screens.json', output, shared_pool=True)['tasks'] == 4
    system = review._saved(output)['systems'][0]
    config = review.recon._load(system['configuration'])
    manifest = review._read(system['manifest'])
    child = Path(system['manifest']).parents[1]
    source = review._saved(Path(review._saved(child)['source_review']))
    prior = review.recon._load(source['systems'][0]['configuration'])
    assert config['_probe_model'] == prior['_probe_model']
    assert config['casscf_options'] == prior['casscf_options']
    assert config['correlation_checks'] == []
    for anchor, old in zip(config['anchors'], prior['anchors']):
        assert np.array_equal(review.recon._anchor_geometry(config, anchor),
                              review.recon._anchor_geometry(prior, old))
    for b in ('adz','atz'):
        jobs = [j for j in manifest['jobs'] if j['name'].startswith(b+'__')]
        for i, j in enumerate(jobs):
            deck = (Path(j['directory'])/(j['name']+'.inp')).read_text()
            assert 'start,2140.2,type=natural' in deck
            assert 'closed,16' in deck and 'occ,22' in deck
            assert 'state,2\nweight,0.5,0.5' in deck
            assert not any(s in deck for s in ('rotate,', ',ovl', '{avas', 'freeze,', 'frozen,'))
            assert 'orbital_source' not in j['probe_active_space']
            assert ('diab,2140.2' in deck) == bool(i)
            if i:
                assert j['initial_record'] is None
                assert j['orbital_predecessor'] == jobs[0]['name']
            else:
                assert j['initial_record']['path'].endswith(b+'__exit_l2.wfu')
                assert j['initial_record']['sha256'] in deck
                assert j['orbital_predecessor'] is None
    after = snapshot(production)
    review.batch(production, plan/'screens.json', output, shared_pool=True, resume=True)
    assert snapshot(production) == after
    assert all(p.read_bytes() == v for p, v in before.items())


def test_bad_last_basis_seed_blocks_all_staging(natural_union):
    production, output, plan = natural_union
    path = plan/'csh2_sh_exit.json'; recipe = review._read(path)
    recipe['candidates'][0]['anchors'][0]['seed_sha256']['atz'] = '0'*64
    atomic_json(path, recipe); before = snapshot(production)
    with pytest.raises(ValueError, match='explicitly reviewed seed hash'):
        review.batch(production, plan/'screens.json', output, shared_pool=True)
    assert not output.exists() and snapshot(production) == before


def test_native_execution_consumes_natural_checkpoint_without_replacing_source(natural_union, monkeypatch, tmp_path):
    production, output, plan = natural_union
    review.batch(production, plan/'screens.json', output, shared_pool=True)
    path = Path(review._saved(output)['systems'][0]['manifest'])
    job = review._read(path)['jobs'][0]
    source = Path(job['initial_record']['path']); before = source.read_bytes()
    scratch = tmp_path/'scratch'; scratch.mkdir()
    monkeypatch.setattr(review.recon, 'qc_runtime_environment', lambda *a, **k:
        ({}, {'scratch_directory':str(scratch)}))
    monkeypatch.setattr(review.recon, 'cleanup_qc_runtime', lambda *a: None)
    calls = []
    def execute(command, **kwargs):
        import re
        calls.append(command)
        local = kwargs['cwd']; deck = (local/command[-1]).read_text()
        filename = re.search(r'(?m)^file,2,([^,\n]+)$', deck)[1]
        assert (local/filename).read_bytes() == before
        assert 'start,2140.2,type=natural' in deck
        assert 'start,2100.3' not in deck and ',ovl' not in deck
        assert 'closed,16' in deck and 'occ,22' in deck
        assert 'state,2\nweight,0.5,0.5' in deck
        (local/filename).write_bytes(b'new synthetic optimized target checkpoint')
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(review.recon.subprocess, 'run', execute)
    review.recon.execute(path, job['name'])
    assert len(calls) == 1 and source.read_bytes() == before
    execution = review._read(Path(job['directory'])/'execution.json')
    assert execution['native_input_sha256'] == file_sha256(execution['native_input'])
    assert execution['status'] == 'complete'


@pytest.mark.parametrize('change', ['negative', 'bool', 'fraction', 'close_and_open', 'count', 'too_many'])
def test_invalid_natural_promotions_rejected(tmp_path, change):
    data = _partition_chain_configuration(tmp_path)
    definition = data['probe_definition']
    definition.update(active_electrons=8, active_orbitals=7)
    part = definition['orbital_partition']
    part.update(source_active_electrons=6, source_active_orbitals=6,
                close_lowest=0, discard_highest=0, open_highest_closed=1)
    if change in ('negative', 'bool', 'fraction'):
        part['open_highest_closed'] = {'negative':-1, 'bool':True, 'fraction':1.5}[change]
    elif change == 'close_and_open':
        part['close_lowest'] = 1
    elif change == 'count':
        definition['active_electrons'] = 6
    else:
        part.update(source_active_electrons=14, source_active_orbitals=8)
        definition.update(active_electrons=16, active_orbitals=9)
    path = tmp_path/'bad.json'; atomic_json(path, data)
    with pytest.raises(ValueError):
        review.recon.prepare(path, tmp_path/'out')
    assert not list((tmp_path/'out').rglob('*.inp'))
    assert not (tmp_path/'out/reconnaissance.json').exists()


def test_four_job_schedule_checks_predecessors_and_preserves_22_live_jobs(natural_union, monkeypatch):
    production, output, plan = natural_union
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
        directory = Path(kwargs['cwd']); state = review._read(directory.parent/'reconnaissance.json')
        job = next(j for j in state['jobs'] if j['directory'] == str(directory))
        if job['orbital_predecessor']:
            prior = next(j for j in state['jobs'] if j['name'] == job['orbital_predecessor'])
            assert prior['job_id'] in done and prior['status'] == 'complete'
            assert job['initial_wfu_sha256'] == file_sha256(prior['wfu'])
            assert (directory/(job['name']+'.wfu')).read_bytes() == Path(prior['wfu']).read_bytes()
        ident = str(1000+len(launched)); launched.append((ident,directory)); running[ident] = 1
        return SimpleNamespace(returncode=0, stdout=ident+'\n', stderr='')
    monkeypatch.setattr(review.subprocess, 'run', scheduler)
    def audit(path, **kwargs):
        state = review._read(path)
        for j in state['jobs']:
            j['status'] = 'complete' if j.get('job_id') in done else 'incomplete'
            if j['status'] == 'complete':
                wfu = Path(j['directory'])/(j['name']+'.wfu'); wfu.write_bytes(('target:'+j['name']).encode())
                j.update(wfu=str(wfu), wfu_sha256=file_sha256(wfu))
        state['status'] = 'complete_review_required' if all(j['status']=='complete' for j in state['jobs']) else 'incomplete'
        atomic_json(path,state); return state
    monkeypatch.setattr(review.recon, '_audit', audit)
    monkeypatch.setattr(review.recon, '_recover_slurm_startup_failures', lambda *a: None)
    monkeypatch.setattr(review.recon, '_retry_allocation_failures', lambda p,m,*a: m)
    def finish(_):
        for ident in ['100', *[i for i,d in launched]]:
            if ident in running:
                running.pop(ident); done.add(ident)
    monkeypatch.setattr(review.recon.time, 'sleep', finish)
    assert review.run(output, report_failures=True)['status'] == 'complete_review_required'
    assert len(launched) == len({str(d) for i,d in launched}) == 4
    assert running == {str(100+i):1 for i in range(1,23)}
    assert policy_path.read_bytes() == policy_bytes
    assert review._read(production/'.kinbot_policy_admissions.json') == {}
    review.run(output, report_failures=True)
    assert len(launched) == 4
