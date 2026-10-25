"""Failure-class recovery: synthetic records, no claim of native QC success."""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import molpro_reconnaissance as recon, requested_orbital_recon as review
from kinbot import recon_scf_seed
from kinbot.execution import atomic_json, file_sha256
from test_compact_scf_followup import tree, original_tree, plan, followup_plan, snapshot


def failed(job, identifier='51'):
    directory = Path(job['directory']); name = job['name']
    job.update(job_id=identifier, status='incomplete', submission_status='submitted')
    native = directory/(name+'.native.inp')
    native.write_text('file,2,new.wfu,new\nfile,3,source.wfu\n'
                      '{rhf;start,2100.3,type=canonical,ovl;wf,41,1,1}\n')
    (directory/(name+'.out')).write_text(
        ' PROGRAM * Restricted Hartree-Fock\n'
        ' ?ERROR IN VIRTORB: INCORRECT NUMBER OF ORBITALS: SYM=1 N=159 M=21 N-M=138 NV=155\n'
        ' GLOBAL ERROR fehler on processor 0\n')
    atomic_json(directory/'execution.json', dict(status='failed', slurm_job_id=identifier,
        input_sha256=job['input_sha256'], native_input=str(native), native_input_sha256=file_sha256(native),
        input_transform='native_xml_scf_seed_file_v7', scf_seed=job['initial_record'],
        error='Molpro exited with status 1.'))


def completed(job, identifier, *, invalid=False):
    d = Path(job['directory']); name = job['name']; space = job['probe_active_space']
    closed = (41-space['active_electrons'])//2
    occupations = [2.]*closed + [1.]*space['active_orbitals']
    occupations[-1] += space['active_electrons']-space['active_orbitals']
    if invalid: occupations[0] = 1.
    native = d/(name+'.native.inp'); native.write_text('file,2,final.wfu,new\n')
    (d/'final.wfu').write_bytes(b'completed native orbitals')
    header = 'state 1.1' if len(job['root_states']) == 1 else 'state averaged'
    (d/(name+'.out')).write_text(
        'PROGRAM * MULTI\nCONVERGENCE REACHED!\n'
        f'Number of active orbitals: {space["active_orbitals"]}\n'
        f'Number of active electrons: {space["active_electrons"]}\n'
        f'NATURAL ORBITALS ({header})\n'
        + ''.join(f'{i}.1 {v:.5f} -0.1\n' for i,v in enumerate(occupations, 1))
        + f'Natural orbital dump ({header})\n'
        'SETTING RECONENERGY = -834.64153362 HARTREE\n'
        + ''.join(f'SETTING RECONR{i} = -834.64153362 HARTREE\n'
                  for i in range(1,len(job['root_states'])+1))
        + 'SETTING RECONDONE = 1\nMolpro calculation terminated\n')
    exports = {}
    for kind in ('scf', 'natural'):
        p = d/(kind+'.xml'); p.write_text('<molpro/>'); exports[kind] = str(p)
    atomic_json(d/'execution.json', dict(status='complete', returncode=0, slurm_job_id=identifier,
        input_sha256=job['input_sha256'], native_input=str(native), native_input_sha256=file_sha256(native),
        input_transform='native_xml_scf_seed_file_v7', orbital_exports=exports))
    job.update(job_id=identifier, status='incomplete', terminal_failure='Invalid natural occupations: count=0')
    return occupations


def legacy_completed(job, identifier):
    completed(job, identifier)
    d = Path(job['directory']); native = d/(job['name']+'.native.inp')
    native.write_text('file,2,final.wfu,new\nfile,3,source.wfu\n'
        '{rhf;start,2100.3,type=canonical,ovl;wf,41,1,1}\n'
        '{multi,so-sci;\nstart,2100.3,type=canonical,ovl\n}\n')
    output = d/(job['name']+'.out')
    output.write_text('Program * Restricted Hartree-Fock\n'+output.read_text())
    record = d/'execution.json'; execution = review._read(record)
    execution.update(native_input_sha256=file_sha256(native), scf_seed=job['initial_record'])
    atomic_json(record, execution)


@pytest.fixture
def campaign(tree, followup_plan):
    production, output = tree
    review.batch(production, followup_plan/'screens.json', output, shared_pool=True)
    manifests = [Path(s['manifest']) for s in review._saved(output)['systems']]
    path = manifests[0]; manifest = review._read(path)
    failed(manifest['jobs'][0]); atomic_json(path, manifest)
    return production, output, manifests


def test_retry_archives_only_proven_failure_preserving_live_peers_and_repeat(campaign, monkeypatch):
    production, output, paths = campaign
    path = paths[0]; manifest = review._read(path)
    peer = manifest['jobs'][1]; peer.update(job_id='52', submission_status='submitted')
    atomic_json(path, manifest)
    protected = {p: p.read_bytes() for p in production.rglob('*') if p.is_file()
                 and Path(manifest['jobs'][0]['directory']) not in p.parents and p != path}
    monkeypatch.setattr(recon, 'slurm_job_record', lambda ident: dict(JobId=ident,
        WorkDir=manifest['jobs'][0]['directory'], JobState='FAILED'))
    recovered = recon._retry_seed_projection_failures(path, manifest, {'52'})
    job, live = recovered['jobs'][:2]
    assert job['attempt'] == 2 and job['submission_status'] == 'retry_ready'
    assert 'job_id' not in job
    assert all(live[k] == v for k,v in peer.items())
    history = job['attempt_history'][0]
    assert history['job_id'] == '51' and history['recovery']['kind'] == 'same_basis_seed_restart'
    assert history['input_transform'] == 'native_xml_scf_seed_file_v7'
    assert (Path(history['archive'])/'execution.json').is_file()
    assert all(p.read_bytes() == b for p,b in protected.items())
    before = snapshot(production)
    recon._retry_seed_projection_failures(path, recovered, {'52'})
    assert snapshot(production) == before
    # Even a stale executor producing the old failure cannot repeat this setting.
    failed(job, '53'); atomic_json(path, recovered)
    assert recon_scf_seed.projection_failure(job) is None


@pytest.mark.parametrize('change', ['live', 'new_transport', 'wrong_job', 'input_hash',
    'native_hash', 'source_provenance', 'different_basis', 'multi_failure', 'other_error', 'submitting'])
def test_recovery_never_resets_inapplicable_work(campaign, monkeypatch, change):
    production, output, paths = campaign
    path = paths[0]; manifest = review._read(path); job = manifest['jobs'][0]
    d = Path(job['directory']); record = d/'execution.json'; execution = review._read(record)
    live = set()
    if change == 'live': live = {'51'}
    elif change == 'new_transport': execution['input_transform'] = 'native_xml_scf_seed_same_basis_v8'
    elif change == 'wrong_job': execution['slurm_job_id'] = '50'
    elif change == 'input_hash': execution['input_sha256'] = 'wrong'
    elif change == 'native_hash': (d/(job['name']+'.native.inp')).write_text('changed input')
    elif change == 'source_provenance': execution['scf_seed'] = {}
    elif change == 'different_basis': job['basis'] = 'avtz-f12'
    elif change == 'multi_failure':
        p = d/(job['name']+'.out'); p.write_text(p.read_text()+'PROGRAM * MULTI\n')
    elif change == 'other_error': (d/(job['name']+'.out')).write_text('PROGRAM * RHF-SCF\nNo convergence')
    elif change == 'submitting': job['submission_status'] = 'submitting'
    atomic_json(record, execution); atomic_json(path, manifest)
    monkeypatch.setattr(recon, '_retry_failed', lambda *a,**k: pytest.fail('No retry'))
    before = snapshot(production)
    recon._retry_seed_projection_failures(path, manifest, live)
    assert snapshot(production) == before


def test_changed_seed_refuses_recovery_before_archiving(campaign):
    production, _, paths = campaign
    path = paths[0]; manifest = review._read(path); job = manifest['jobs'][0]
    Path(job['initial_record']['path']).write_bytes(b'new orbital solution')
    before = snapshot(production)
    with pytest.raises(ValueError, match='evidence changed'):
        recon._retry_seed_projection_failures(path, manifest, set())
    assert snapshot(production) == before


@pytest.fixture
def completed_pair(campaign, monkeypatch):
    production, output, paths = campaign
    paths = paths[4:6]
    for p, ident in zip(paths, ('61', '62')):
        m = review._read(p); legacy_completed(m['jobs'][0], ident); atomic_json(p, m)
    monkeypatch.setattr(review.subprocess, 'run', lambda *a, **k: SimpleNamespace(stdout='', returncode=0))
    monkeypatch.setattr(recon, 'slurm_job_record', lambda i: {})
    return production, output, paths


def test_explicit_seed_recheck_preserves_science_and_archives_completed_results_once(completed_pair):
    production, output, paths = completed_pair
    jobs = [review._read(p)['jobs'][0] for p in paths]
    original = [snapshot(Path(j['directory'])) for j in jobs]
    policy = production/'.kinbot_resource_policy.json'; policy_before = policy.read_bytes()
    review._recheck_seed_jobs(output, ['61', '62'])
    for path, before, files in zip(paths, jobs, original):
        after = review._read(path)['jobs'][0]
        assert after['attempt'] == 2 and 'job_id' not in after
        assert after['submission_status'] == 'retry_ready'
        for key in ('input_sha256', 'initial_record', 'root_states', 'probe_active_space', 'basis', 'anchor'):
            assert after[key] == before[key]
        archive = Path(after['attempt_history'][0]['archive'])
        assert after['attempt_history'][0]['recovery']['kind'] == 'reviewed_seed_recheck'
        assert all((archive/p.name).read_bytes() == b for p,b in files.items())
    review._recheck_seed_jobs(output, ['61', '62'])
    assert all(len(review._read(p)['jobs'][0]['attempt_history']) == 1 for p in paths)
    assert policy.read_bytes() == policy_before


@pytest.mark.parametrize('problem', ['unknown', 'live', 'hash', 'dependency', 'accepted', 'scheduler'])
def test_recheck_preflights_every_target_before_archiving_any(completed_pair, monkeypatch, problem):
    production, output, paths = completed_pair
    ids = ['61', '62']
    if problem == 'unknown': ids.append('999')
    if problem == 'live':
        monkeypatch.setattr(review.subprocess, 'run', lambda *a, **k: SimpleNamespace(stdout='62\n', returncode=0))
    if problem == 'hash':
        job = review._read(paths[1])['jobs'][0]
        record = Path(job['directory'])/'execution.json'; execution = review._read(record)
        execution['scf_seed'] = {}; atomic_json(record, execution)
    if problem == 'dependency':
        m = review._read(paths[1]); m['jobs'][1]['orbital_predecessor'] = m['jobs'][0]['name']; atomic_json(paths[1], m)
    if problem == 'accepted':
        child = next(Path(c) for c in review._read(output/'review_launch.json')['reviews']
                     if 'hcs_h2s_exit_compact' in c)
        request = Path(review._read(child/'review_launch.json')['request'])
        q = review._read(request); q['status'] = 'model_accepted'; atomic_json(request, q)
    if problem == 'scheduler':
        monkeypatch.setattr(recon, 'slurm_job_record', lambda i: dict(JobId=i, WorkDir='/wrong', JobState='COMPLETED'))
    with pytest.raises((ValueError, RuntimeError)):
        review._recheck_seed_jobs(output, ids)
    assert all(not review._read(p)['jobs'][0].get('attempt_history') for p in paths)
    assert all((Path(review._read(p)['jobs'][0]['directory'])/'execution.json').is_file() for p in paths)


def test_recheck_resumes_partial_archival_without_repeating_first_target(completed_pair, monkeypatch):
    _, output, paths = completed_pair
    real_retry = recon._retry_failed
    def partial(path, **kwargs):
        if path == paths[1]: raise OSError('simulated interruption')
        return real_retry(path, **kwargs)
    monkeypatch.setattr(recon, '_retry_failed', partial)
    with pytest.raises(OSError, match='simulated interruption'):
        review._recheck_seed_jobs(output, ['61', '62'])
    monkeypatch.setattr(recon, '_retry_failed', real_retry)
    review._recheck_seed_jobs(output, ['61', '62'])
    assert all(len(review._read(p)['jobs'][0]['attempt_history']) == 1 for p in paths)


@pytest.mark.parametrize('invalid', [False, True])
def test_single_state_native_completion_is_reaudited_without_reexecution(campaign, monkeypatch, invalid):
    production, _, paths = campaign
    path = paths[4]; manifest = review._read(path); job = manifest['jobs'][0]
    assert len(job['root_states']) == 1
    d = Path(job['directory'])
    occupations = completed(job, '61', invalid=invalid)
    atomic_json(path, manifest)
    before = snapshot(d)
    monkeypatch.setattr(recon.subprocess, 'run', lambda *a,**k: pytest.fail('No QC or submission needed'))
    after = recon._audit(path, live_job_ids=set())
    after = recon._retry_seed_projection_failures(path, after, set())
    assert after['jobs'][0]['status'] == ('incomplete' if invalid else 'complete')
    assert after['jobs'][0]['job_id'] == '61'
    assert 'attempt_history' not in after['jobs'][0]
    assert snapshot(d) == before
    if not invalid:
        assert after['jobs'][0]['natural_occupations'] == occupations
        assert 'terminal_failure' not in after['jobs'][0]


@pytest.mark.parametrize('explicit_recheck', [False, True])
def test_campaign_recovery_retains_live_work_releases_claims_and_resumes_once(campaign, monkeypatch, explicit_recheck):
    """Actual audit/archive/admission loop; only scheduler and QC are simulated."""
    from kinbot import run_policy, run_walltime
    production, output, paths = campaign
    paths = [paths[0], paths[4], paths[5]] if explicit_recheck else [paths[0], paths[4]]
    manifest = review._read(paths[0]); peer = manifest['jobs'][1]
    peer.update(job_id='52', submission_status='submitted'); atomic_json(paths[0], manifest)
    for p, identifier in zip(paths[1:], ('61', '62')):
        manifest = review._read(p)
        (legacy_completed if explicit_recheck else completed)(manifest['jobs'][0], identifier)
        atomic_json(p, manifest)
    protected = snapshot(Path(manifest['jobs'][0]['directory']))
    policy_path = production/'.kinbot_resource_policy.json'; policy = review._read(policy_path)
    policy['fair_submission'] = False; atomic_json(policy_path, policy); policy_before = policy_path.read_bytes()
    directories = {str(j['job_id']):j['directory'] for p in paths for j in review._read(p)['jobs'] if j.get('job_id')}
    claim = dict(job_id='51', status='submitted', directory=directories['51'],
                 command=['sbatch','--parsable','job.sbatch'])
    atomic_json(production/'.kinbot_policy_admissions.json', dict(old=claim))
    atomic_json(production/'.kinbot_policy_submissions/51.json', claim)
    running = {str(i):1 for i in range(100,122)} | {'52':1}
    done = {'51','61','62'}; launched = []
    def record(ident):
        if ident in directories:
            return dict(JobId=ident, WorkDir=directories[ident],
                        JobState='FAILED' if ident == '51' else 'COMPLETED' if ident in done else 'RUNNING')
        return {}
    monkeypatch.setattr(recon, 'slurm_job_record', record)
    monkeypatch.setattr(run_policy, 'slurm_job_record', record)
    monkeypatch.setattr(run_policy, '_queue', lambda root: running.copy())
    monkeypatch.setattr(run_walltime, 'partition_walltime', lambda p: '31-00:00:00')
    def scheduler(command, **kwargs):
        if command[0] == 'squeue':
            return SimpleNamespace(stdout='\n'.join(running), returncode=0)
        assert command[0] == 'sbatch' and len(running) < 23
        assert '--exclude=node34' in command
        ident = str(1000+len(launched)); directory = str(kwargs['cwd'])
        launched.append((ident,directory)); directories[ident] = directory; running[ident] = 1
        return SimpleNamespace(stdout=ident+'\n', stderr='', returncode=0)
    monkeypatch.setattr(recon.subprocess, 'run', scheduler)
    def finish(_):
        for p in paths:
            for j in review._read(p)['jobs']:
                ident = str(j.get('job_id'))
                if ident in running:
                    completed(j, ident); running.pop(ident); done.add(ident)
    monkeypatch.setattr(recon.time, 'sleep', finish)
    def runner(_):
        return recon.run_many(paths, max_jobs=3, keep_going=True, allocation_retries=2)
    if explicit_recheck:
        monkeypatch.setattr(review, '_run', runner)
        result = review.run(output, report_failures=True, recheck_seed_jobs=['61', '62'])
    else:
        result = runner(output)
    assert all(m['status'] in ('complete_review_required','coverage_refinement_required') for m in result)
    expected = 11 if explicit_recheck else 6
    assert len(launched) == len({d for _,d in launched}) == expected
    assert all(d != directories['52'] for _,d in launched)
    assert running == {str(i):1 for i in range(100,122)}
    assert policy_path.read_bytes() == policy_before
    if explicit_recheck:
        archive = Path(review._read(paths[-1])['jobs'][0]['attempt_history'][0]['archive'])
        assert all((archive/p.name).read_bytes() == b for p,b in protected.items())
    else:
        assert all(d != directories['61'] for _,d in launched)
        assert all(p.read_bytes() == b for p,b in protected.items())
    run_policy.reconcile_allocations(production)
    assert review._read(production/'.kinbot_policy_admissions.json') == {}
    if explicit_recheck:
        review.run(output, report_failures=True, recheck_seed_jobs=['61', '62'])
    else:
        runner(output)
    assert len(launched) == expected and len(review._read(paths[0])['jobs'][0]['attempt_history']) == 1


@pytest.mark.parametrize('heading', ['NATURAL ORBITALS', 'NATURAL ORBITALS (state 1.1)',
                                    'NATURAL ORBITALS (state averaged)'])
def test_single_state_occupations_and_table_boundaries(heading):
    text = (heading+'\n================\n Orbital Occupation Energy\n'
            ' 1.1 2.00000 -1.0 1 1 s 0.50000\n'
            '                 1 2 s 0.30000\n'
            ' 2.1 1.00000 -0.1\n Natural orbital dump\n 3.1 2.0 -0.5\n')
    assert recon._native_occupations(text, [dict(root=1, multiplicity=2)]) == [2., 1.]
    if 'averaged' not in heading:
        assert recon._native_occupations(text, [dict(root=i, multiplicity=2) for i in (1,2)]) == []
    if 'state 1.1' in heading:
        assert recon._native_occupations(text, [dict(root=2, multiplicity=2)]) == []


def test_supplied_single_state_native_excerpt():
    text = (Path(__file__).parent/'fixtures/molpro_single_state_natural_excerpt.out').read_text()
    occupations = recon._native_occupations(text, [dict(root=1, multiplicity=2)])
    assert occupations == [2.]*19 + [1.99965, 1., .00035, 0.]
    assert sum(occupations) == pytest.approx(41.)
    assert recon._native_occupations(text, [dict(root=i, multiplicity=2) for i in (1, 2)]) == []


def test_owner_wait_releases_no_jobs_and_acquires_once(tmp_path, monkeypatch, capsys):
    held = (tmp_path/'.reconnaissance.lock').open('a+')
    recon.fcntl.flock(held, recon.fcntl.LOCK_EX | recon.fcntl.LOCK_NB)
    sleeps = []
    def release(seconds):
        sleeps.append(seconds); held.close()
    monkeypatch.setattr(recon.time, 'sleep', release)
    with recon._manifest_lock(tmp_path, wait=60):
        with pytest.raises(RuntimeError, match='Another reconnaissance controller'):
            with recon._manifest_lock(tmp_path):
                pytest.fail('Second owner')
    assert sleeps == [5.] and 'its jobs continue' in capsys.readouterr().out


@pytest.mark.parametrize('wait', [-1., float('nan'), float('inf')])
def test_owner_wait_rejects_invalid_duration(tmp_path, wait):
    with pytest.raises(ValueError, match='finite and nonnegative'):
        with recon._manifest_lock(tmp_path, wait=wait): pass
