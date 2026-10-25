from copy import deepcopy
import json

import pytest
from ase import Atoms

from kinbot import live_control as control
from kinbot.anl import dispatch, protocol, anharmonic, live_control as owner
from kinbot.anl.validation import automatic_composite_spec
from kinbot.execution import atomic_json, file_sha256
from tests.test_anl_dispatch import _mock_execution


def specification(qz=True, tier='ANL1-F12', stationary='minimum'):
    a = Atoms('CH4', positions=[[0,0,0],[.63,.63,.63],[-.63,-.63,.63],[-.63,.63,-.63],[.63,-.63,-.63]])
    molecule = dict(symbols=a.get_chemical_symbols(), positions=a.positions.tolist(), charge=0, multiplicity=1,
        smiles='C', source=dict(geometry_sha256=dispatch._geometry_hash(a), record_sha256='2'*64,
        energy_hartree=-40, harmonic_zpe_hartree=.04,
        l2_theory=dict(calculator='gaussian', method='B2PLYP', basis='cc-pVTZ', dispersion='GD3BJ')))
    s = automatic_composite_spec(molecule, highest_tier=tier, f12_augmented=True,
                                anharmonic_qz=qz, stationary_point=stationary)
    fixed_resources(s)
    return s


def fixed_resources(spec):
    for t in spec['tasks']:
        if t['resources'].get('cores') == 'auto':
            t['resources'] = dict(cores=1, memory_mb=32000, walltime='02:00:00', partition='test')


@pytest.fixture
def graph(tmp_path, monkeypatch):
    control.configure(tmp_path)
    monkeypatch.setattr(dispatch, 'assign_partitions', fixed_resources)
    g = dispatch.prepare(specification(), tmp_path/'graph')
    _mock_execution(g, 'l3_geometry_tz', geometry='final.xyz')
    dispatch.advance(g)
    state = dispatch._load(g)[2]
    state['tasks']['l3_geometry'].update(status='submitted', job_id='42')
    state['tasks']['gaussian_vpt2'].update(status='submitted', job_id='99')
    atomic_json(g/'state.json', state)
    monkeypatch.setattr(dispatch, '_job_active', lambda j: j in ('42', '99'))
    return g


def queue(g, options, cancel=False):
    p = protocol.plan(g, options)
    r = protocol.publish(g.parent, [p], reason='reviewed protocol change', cancel_running=cancel)[0]
    return g.parent/protocol.QUEUE/(r['id']+'.json')


@pytest.mark.parametrize('tier', ['ANL1', 'ANL1-F12'])
def test_default_qz_branch_is_independent_l3_owned_and_idempotent(tier):
    s = specification(tier=tier); tasks = {t['id']:t for t in s['tasks']}
    assert tasks[anharmonic.GEOMETRY]['geometry_from'] == 'initial'
    assert tasks[anharmonic.GEOMETRY]['profile']['basis'] == 'cc-pVQZ'
    assert tasks[anharmonic.GEOMETRY]['optimizer']['fmax'] == 1e-5
    assert tasks[anharmonic.VPT2]['geometry_from'] == anharmonic.GEOMETRY
    assert tasks[anharmonic.VPT2]['depends_on'] == [anharmonic.HARMONIC]
    assert s['intent']['anharmonic_basis_extension']['extrapolation_power'] == 3.7
    assert anharmonic.extend(s) == s
    assert anharmonic.VPT2 not in {t['id'] for t in specification(stationary='transition_state')['tasks']}


def test_preview_preserves_tz_and_fences_qz_descendants(graph):
    p = protocol.plan(graph, dict(geometry='tz', harmonic='tz'))
    assert 'l3_geometry_tz' not in p['reset'] and 'harmonic_tz' not in p['reset']
    assert 'l3_geometry' in p['reset']
    revised = {t['id']:t for t in p['replacement']['tasks']}
    assert revised['f12_qz']['geometry_from'] == 'l3_geometry_tz'
    assert {'l3_geometry', 'harmonic_qz'} <= protocol.omitted(p['replacement'])
    q = queue(graph, dict(geometry='tz', harmonic='tz'))
    owner.process(graph)
    assert json.loads(q.read_text())['status'] == 'applying'
    state = dispatch.advance(graph)
    assert state['tasks']['l3_geometry']['job_id'] == '42'
    assert not (graph/'tasks/f12_qz').exists()


def test_explicit_cancel_archives_affected_only_and_preserves_allocations(graph, monkeypatch):
    keep = {p:file_sha256(p) for n in ('l3_geometry_tz', 'gaussian_vpt2')
            for p in (graph/'tasks'/n).iterdir() if p.is_file()}
    q = queue(graph, dict(geometry='tz', harmonic='tz'), cancel=True)
    calls = []
    monkeypatch.setattr(owner, 'stop_attempt', lambda d,j,*a,**k: calls.append(j) or True)
    owner.process(graph)
    assert calls == ['42']
    assert json.loads(q.read_text())['status'] == 'applied'
    state = dispatch.advance(graph)
    assert state['tasks']['l3_geometry']['status'] == 'skipped'
    assert state['tasks']['f12_qz']['status'] == 'staged'
    assert state['tasks']['gaussian_vpt2']['job_id'] == '99'
    assert keep == {p:file_sha256(p) for p in keep}
    owner.process(graph)
    assert calls == ['42']


def test_single_point_omission_leaves_shared_reference_allocations_untouched(graph, monkeypatch):
    monkeypatch.setattr(dispatch, '_job_active', lambda j: j == '99')
    _mock_execution(graph, 'l3_geometry', geometry='final.xyz')
    dispatch.advance(graph)
    state = dispatch._load(graph)[2]
    state['tasks']['ccsdtqp_dz'].update(status='submitted', job_id='43')
    state['tasks']['ccsdtq_tz'].update(status='submitted', job_id='44')
    atomic_json(graph/'state.json', state)
    reference = dispatch.prepare(specification(), graph.parent/'shared_reference')
    reference_state = dispatch._load(reference)[2]
    reference_state['tasks']['l3_geometry_tz'].update(status='submitted', job_id='73')
    atomic_json(reference/'state.json', reference_state)
    before = {p:file_sha256(p) for p in reference.rglob('*') if p.is_file()}
    monkeypatch.setattr(dispatch, '_job_active', lambda j: j in ('43', '44', '73', '99'))
    stopped = []
    monkeypatch.setattr(owner, 'stop_attempt', lambda d,j,*a,**kw: stopped.append(j) or True)
    assert protocol.plan(graph, dict(omit_tasks=['ccsdtqp_dz']))['reset'] == ['ccsdtqp_dz']
    q = queue(graph, dict(omit_tasks=['ccsdtqp_dz']), cancel=True)
    owner.process(graph)
    assert json.loads(q.read_text())['status'] == 'applied'
    assert stopped == ['43']
    assert dispatch._load(graph)[2]['tasks']['ccsdtq_tz']['job_id'] == '44'
    assert before == {p:file_sha256(p) for p in reference.rglob('*') if p.is_file()}
    assert dispatch._load(reference)[2]['tasks']['l3_geometry_tz']['job_id'] == '73'


def test_changed_generation_and_ambiguous_submission_are_not_cancelled(graph, monkeypatch):
    q = queue(graph, dict(geometry='tz'), cancel=True)
    s = dispatch._load(graph)[2]; s['tasks']['l3_geometry']['job_id'] = '43'
    atomic_json(graph/'state.json', s)
    monkeypatch.setattr(owner, 'stop_attempt', lambda *a, **k: pytest.fail('No cancellation'))
    owner.process(graph)
    assert json.loads(q.read_text())['status'] == 'blocked'
    s['tasks']['l3_geometry']['status'] = 'submitting'; atomic_json(graph/'state.json', s)
    with pytest.raises(RuntimeError, match='ambiguous'):
        protocol.plan(graph, dict(geometry='tz'))


def test_tier_and_component_options_are_composable(monkeypatch):
    monkeypatch.setattr(dispatch, 'assign_partitions', fixed_resources)
    source = specification()
    low = protocol.replacement(source, dict(tier='ANL0-F12'))
    assert protocol.sources(low) == ('l3_geometry_tz', ['harmonic_tz'])
    assert 'l3_geometry_tz' not in protocol.omitted(low)
    assert next(t for t in low['tasks'] if t['id']=='harmonic_tz') == next(t for t in source['tasks'] if t['id']=='harmonic_tz')
    reduced = protocol.replacement(source, dict(omit_components=['dboc'], omit_tasks=['ccsdtqp_dz']))
    assert {'cfour_dboc', 'ccsdtqp_dz'} <= protocol.omitted(reduced)
    assert 'cv_ae_qz' not in protocol.omitted(reduced)
    restored = protocol.replacement(reduced, dict(omit_components=[], omit_tasks=[]))
    assert 'cfour_dboc' not in protocol.omitted(restored)
    for tier in protocol.TIERS:
        changed = protocol.replacement(source, dict(tier=tier))
        assert changed['intent']['highest_tier'] == tier


def test_old_graph_migration_preserves_running_jobs_and_repeats_safely(tmp_path, monkeypatch):
    control.configure(tmp_path)
    g = dispatch.prepare(specification(qz=False), tmp_path/'graph')
    state = dispatch._load(g)[2]
    state['tasks']['l3_geometry_tz'].update(status='submitted', job_id='42')
    atomic_json(g/'state.json', state)
    files = {p:file_sha256(p) for p in (g/'tasks/l3_geometry_tz').iterdir() if p.is_file()}
    owner.process(g); owner.process(g)
    assert dispatch._load(g)[2]['tasks']['l3_geometry_tz']['job_id'] == '42'
    assert files == {p:file_sha256(p) for p in files}
    assert len(list((g/'attempts').glob('live_*'))) == 1
    monkeypatch.setattr(dispatch, '_job_active', lambda j: True)
    state = dispatch.advance(g)
    assert state['tasks'][anharmonic.GEOMETRY]['status'] == 'staged'


def test_legacy_worker_guard_precedes_any_cancellation(graph, monkeypatch):
    path = graph/'tasks/gaussian_vpt2/task.json'
    record = json.loads(path.read_text()); record.pop('live_control_schema')
    atomic_json(path, record)
    q = queue(graph, dict(geometry='tz'), cancel=True)
    monkeypatch.setattr(owner, 'stop_attempt', lambda *a, **k: pytest.fail('No cancellation'))
    owner.process(graph)
    receipt = json.loads(q.read_text())
    assert receipt['status'] == 'queued' and 'legacy workers' in receipt['error']


def test_interrupted_protocol_commit_retains_jobs_and_resumes_once(graph, monkeypatch):
    from kinbot import run_recover
    q = queue(graph, dict(geometry='tz'), cancel=True)
    monkeypatch.setattr(owner, 'stop_attempt', lambda *a, **k: True)
    save = run_recover.atomic_json
    def interrupt(path, data):
        if path == graph/'state.json':
            raise OSError('interrupted graph write')
        save(path, data)
    monkeypatch.setattr(run_recover, 'atomic_json', interrupt)
    with pytest.raises(OSError, match='interrupted'):
        owner.process(graph)
    monkeypatch.setattr(run_recover, 'atomic_json', save)
    owner.process(graph); owner.process(graph)
    assert json.loads(q.read_text())['status'] == 'applied'
    assert dispatch._load(graph)[2]['tasks']['gaussian_vpt2']['job_id'] == '99'
    assert len(list((graph/'attempts').glob('live_*'))) == 1


def test_closed_graph_queues_owner_revision_before_modification(tmp_path, monkeypatch):
    from kinbot import pes_phase
    from kinbot.anl import refinement
    control.configure(tmp_path)
    worker = tmp_path/'case/123'; atomic_json(worker/'123.json', {})
    g = dispatch.prepare(specification(qz=False), worker/'anl_automatic/states/a/composite')
    refinement.window(g, close=True)
    monkeypatch.setattr(pes_phase, 'active_worker', lambda w: w == worker)
    monkeypatch.setattr(dispatch, 'assign_partitions', fixed_resources)
    p = protocol.plan(g, dict(anharmonic_qz=True), allow_closed=True)
    before = file_sha256(g/'workflow.json')
    r = protocol.publish(tmp_path, [p], reason='enable new default branch')[0]
    path = tmp_path/protocol.QUEUE/(r['id']+'.json')
    assert json.loads(path.read_text())['owner_request']
    assert len(control.requests(worker, kinds={'network'})) == 1
    protocol.process(g)
    assert before == file_sha256(g/'workflow.json')
    # Model the normal owner epoch boundary, after derived assembly retires.
    refinement.window(g)
    owner.process(g)
    assert json.loads(path.read_text())['status'] == 'applied'
    assert anharmonic.VPT2 in {t['id'] for t in dispatch._load(g)[1]['tasks']}
    assert refinement.window(g, close=True)


def test_new_branch_site_preflight_is_cached_and_resumable(graph, monkeypatch):
    atomic_json(graph/'protocol_preflight.json', dict(status='pending'))
    seen = []
    def failed(*a, **k):
        seen.append('failed'); raise ValueError('site command unavailable')
    monkeypatch.setattr(dispatch, 'preflight', failed)
    assert not protocol.preflight_ready(graph)
    assert not protocol.preflight_ready(graph) and seen == ['failed']
    atomic_json(graph/'protocol_preflight.json', dict(status='pending'))
    monkeypatch.setattr(dispatch, 'preflight', lambda *a, **k: seen.append('verified'))
    assert protocol.preflight_ready(graph) and protocol.preflight_ready(graph)
    assert seen == ['failed', 'verified']
    assert dispatch._load(graph)[2]['tasks']['l3_geometry']['job_id'] == '42'


@pytest.mark.parametrize('version', [1, 2])
def test_every_tier_roundtrip_and_site_commands_are_preserved(monkeypatch, version):
    monkeypatch.setattr(dispatch, 'assign_partitions', fixed_resources)
    mol = specification()['molecule']
    for tier in protocol.TIERS:
        source = automatic_composite_spec(mol, highest_tier=tier, component_version=version)
        for t in source['tasks']:
            if t.get('result_parser', {}).get('kind', '').startswith('molpro'):
                t['command'][0] = '/site/molpro'
        revised = protocol.replacement(source, dict(tier='ANL1-F12'))
        dispatch.validate_spec(revised)
        if any(t.get('result_parser', {}).get('kind', '').startswith('molpro') for t in source['tasks']):
            for t in revised['tasks']:
                if t.get('result_parser', {}).get('kind', '').startswith('molpro'):
                    assert t['command'][0] == '/site/molpro'


def test_retained_mess_pair_and_explicit_l2_geometry(graph, monkeypatch):
    from kinbot.anl import thermal
    monkeypatch.setattr(dispatch, 'assign_partitions', fixed_resources)
    _, spec, state = dispatch._load(graph)
    for n in ('l3_geometry', 'l3_geometry_tz', 'harmonic_tz', 'harmonic_qz'):
        state['tasks'].setdefault(n, {})['status'] = 'complete'
    atomic_json(graph/'state.json', state)
    # This tests source selection, not the native-result validator exercised
    # separately by the thermal native fixtures.
    assert thermal.handoff_options(graph)['harmonic_task'] == 'harmonic_qz'
    revised = protocol.replacement(spec, dict(geometry='tz', harmonic='tz'))
    owner.commit(graph, revised, state, set(), {})
    assert thermal.handoff_options(graph) == dict(run=str(graph),
        geometry_task='l3_geometry_tz', harmonic_task='harmonic_tz')
    revised = protocol.replacement(revised, dict(geometry='l2', harmonic='qz'))
    assert protocol.sources(revised) == ('initial', ['harmonic_qz'])


def test_default_extension_does_not_relabel_custom_vpt2_basis():
    s = specification(qz=False)
    t = next(t for t in s['tasks'] if t['id'] == 'gaussian_vpt2')
    t['result_parser']['basis'] = 'custom-basis'
    assert anharmonic.extend(s) == s


def test_native_scheduler_ownership_is_checked_for_explicit_cancellation(graph, monkeypatch):
    from kinbot import run_refresh
    cancelled = []
    def record(job):
        return dict(JobId=job, WorkDir=str(graph/'tasks/l3_geometry'),
                    JobState='CANCELLED' if cancelled else 'RUNNING')
    monkeypatch.setattr(owner, 'slurm_job_record', record)
    monkeypatch.setattr(run_refresh, '_pending_request', lambda argv: cancelled.append(argv) or '')
    q = queue(graph, dict(geometry='tz'), cancel=True)
    owner.process(graph); owner.process(graph)
    assert cancelled == [['scancel', '42']]
    assert json.loads(q.read_text())['status'] == 'applied'
    assert dispatch._load(graph)[2]['tasks']['gaussian_vpt2']['job_id'] == '99'


@pytest.mark.parametrize('initial_phase', ['PENDING', 'RUNNING'])
@pytest.mark.parametrize('terminal_phase', ['CANCELLED', 'COMPLETED'])
def test_protocol_cancellation_preserves_full_pool_and_waits_for_terminal(
        graph, monkeypatch, initial_phase, terminal_phase):
    from types import SimpleNamespace
    from kinbot import run_policy as policy, run_refresh, run_walltime

    root = graph.parent
    policy.update(root, max_jobs=23, exclude=['bad1', 'bad2'],
                  minimum_walltime='7-00:00:00')
    policy_before = (root/policy.NAME).read_bytes()
    workflow_before = (graph/'workflow.json').read_bytes()
    retained = {p:file_sha256(p) for n in ('l3_geometry_tz', 'gaussian_vpt2')
                for p in (graph/'tasks'/n).iterdir() if p.is_file()}
    jobs = {'42': 1, '99': 1, **{str(100+i): 1 for i in range(21)}}
    phase = initial_phase
    cancelled, launched = [], []

    def record(job):
        assert job == '42'
        return dict(JobId=job, WorkDir=str(graph/'tasks/l3_geometry'), JobState=phase)

    def cancel(command):
        nonlocal phase
        cancelled.append(command)
        phase = 'COMPLETING'
        return ''

    def launch(command, **kwargs):
        assert command[0] == 'sbatch'
        launched.append(command)
        jobs['200'] = 1
        return SimpleNamespace(returncode=0, stdout='200\n', stderr='')

    monkeypatch.setattr(owner, 'slurm_job_record', record)
    monkeypatch.setattr(policy, 'slurm_job_record', record)
    monkeypatch.setattr(policy, '_queue', lambda root: dict(jobs))
    monkeypatch.setattr(dispatch, '_job_active', lambda job: job in jobs)
    monkeypatch.setattr(run_refresh, '_pending_request', cancel)
    monkeypatch.setattr(run_walltime, 'partition_limits', lambda **kw: {'test': 604800})
    q = queue(graph, dict(geometry='tz', harmonic='tz'), cancel=True)
    owner.process(graph)
    owner.process(graph)
    assert cancelled == [['scancel', '42']]
    assert json.loads(q.read_text())['status'] == 'applying'
    assert (graph/'workflow.json').read_bytes() == workflow_before
    assert not list((graph/'attempts').glob('live_*'))

    # A finishing allocation still occupies its slot. An old admission claim
    # is retired once visible in Slurm, without double-counting that job.
    atomic_json(root/'.kinbot_policy_admissions.json', {'old': dict(
        job_id='42', directory=str(graph/'tasks/l3_geometry'),
        command=['sbatch', '--parsable', 'job.slurm'], status='submitted')})
    next_job = root/'other_branch'; next_job.mkdir()
    (next_job/'job.slurm').write_text('#!/bin/bash\n#SBATCH --partition=test\n#SBATCH --time=02:00:00\n')
    with monkeypatch.context() as submission:
        submission.setattr(policy.subprocess, 'run', launch)
        assert policy.submit(['sbatch', '--parsable', 'job.slurm'], cwd=next_job,
                             capture_output=True, text=True) is None
    assert not launched and len(jobs) == 23
    assert json.loads((root/'.kinbot_policy_admissions.json').read_text()) == {}

    # Scheduler completion may win the cancellation race; both terminal
    # outcomes permit archival, but neither changes unrelated allocations.
    phase = terminal_phase
    owner.process(graph)
    assert json.loads(q.read_text())['status'] == 'applied'
    jobs.pop('42')
    owner.process(graph)
    state = dispatch._load(graph)[2]
    assert state['tasks']['gaussian_vpt2']['job_id'] == '99'
    assert retained == {p:file_sha256(p) for p in retained}
    assert len(list((graph/'attempts').glob('live_*'))) == 1
    assert (root/policy.NAME).read_bytes() == policy_before

    with monkeypatch.context() as submission:
        submission.setattr(policy.subprocess, 'run', launch)
        assert policy.submit(['sbatch', '--parsable', 'job.slurm'], cwd=next_job,
                             capture_output=True, text=True).returncode == 0
        assert policy.submit(['sbatch', '--parsable', 'job.slurm'], cwd=next_job,
                             capture_output=True, text=True) is None
    assert len(jobs) == 23 and len(launched) == 1
    assert '--exclude=bad1,bad2' in launched[0]
    assert '--time=7-00:00:00' in launched[0]
    assert '--nodes=1' in launched[0] and '--exclusive' in launched[0]
    assert cancelled == [['scancel', '42']]
    assert (root/policy.NAME).read_bytes() == policy_before


def test_cli_preview_and_batch_revalidation_are_write_free(graph, capsys):
    from kinbot.anl import control as cli
    capsys.readouterr()
    before = {p:file_sha256(p) for p in graph.rglob('*') if p.is_file()}
    cli.main(['protocol', str(graph.parent), '--species', 'C', '--geometry', 'tz'])
    preview = json.loads(capsys.readouterr().out)
    assert preview[0]['sources'] == ['l3_geometry_tz', ['harmonic_tz']]
    assert before == {p:file_sha256(p) for p in before}
    p = protocol.plan(graph, dict(geometry='tz'))
    state = dispatch._load(graph)[2]; state['tasks']['l3_geometry']['job_id'] = 'changed'
    atomic_json(graph/'state.json', state)
    with pytest.raises(ValueError, match='preview changed'):
        protocol.publish(graph.parent, [p], reason='stale preview')
    assert not list((graph.parent/protocol.QUEUE).glob('*.json'))


def test_graph_selectors_cover_one_state_subset_and_one_pes(tmp_path):
    from kinbot.anl.control import graphs
    names = ('pes_a/well/anl_automatic/states/target/composite',
             'pes_a/well/anl_automatic/states/reference/composite',
             'pes_b/well/anl_automatic/states/target/composite')
    paths = [tmp_path/n for n in names]
    for path, smiles in zip(paths, ('CCC', 'CC', 'CCC')):
        atomic_json(path/'workflow.json', dict(name='anl-automatic-post-l2',
                                               molecule=dict(smiles=smiles)))
    archived = paths[0]/'attempts/old/graphs'
    atomic_json(archived/'workflow.json', dict(name='anl-automatic-post-l2',
                                             molecule=dict(smiles='CCC')))
    assert graphs(tmp_path, patterns=[names[0]]) == paths[:1]
    assert graphs(tmp_path, patterns=[names[0], names[2]]) == [paths[0], paths[2]]
    assert graphs(tmp_path, patterns=['pes_a/*/anl_automatic/states/*/composite']) == [paths[1], paths[0]]
    assert graphs(tmp_path, species=['CCC']) == [paths[0], paths[2]]
    assert graphs(tmp_path, all_species=True) == [paths[1], paths[0], paths[2]]


def test_omitted_prerequisites_cannot_leave_unrunnable_consumers(monkeypatch):
    monkeypatch.setattr(dispatch, 'assign_partitions', fixed_resources)
    revised = protocol.replacement(specification(), dict(omit_tasks=[anharmonic.HARMONIC]))
    assert {anharmonic.HARMONIC, anharmonic.VPT2} <= protocol.omitted(revised)


def test_disabling_upper_anharmonic_branch_also_disables_recovery_descendants(monkeypatch):
    monkeypatch.setattr(dispatch, 'assign_partitions', fixed_resources)
    spec = specification()
    child = deepcopy(next(t for t in spec['tasks'] if t['id'] == anharmonic.GEOMETRY))
    child.update(id='upper_recovery_geometry', geometry_from=anharmonic.GEOMETRY)
    spec['tasks'].append(child)
    revised = protocol.replacement(spec, dict(anharmonic_qz=False))
    assert {anharmonic.GEOMETRY, anharmonic.HARMONIC, anharmonic.VPT2,
            'upper_recovery_geometry'} <= protocol.omitted(revised)
