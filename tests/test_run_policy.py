"""Cooperative Slurm scheduling; no native QC or actual Slurm is required."""
from copy import deepcopy
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from kinbot import run_policy as policy
from kinbot.execution import atomic_json, file_sha256

ORIGINAL_RUN = subprocess.run


@pytest.fixture
def pool(tmp_path, monkeypatch):
    policy.update(tmp_path, max_jobs=2, exclude=['bad1'])
    jobs, commands = {}, []
    monkeypatch.setattr(policy, '_queue', lambda root: dict(jobs))
    monkeypatch.setattr(policy, 'slurm_job_record', lambda job: {})
    def launch(command, **kwargs):
        if Path(command[0]).name != 'sbatch':
            return ORIGINAL_RUN(command, **kwargs)
        commands.append((command, kwargs))
        ident = str(100 + len(commands))
        jobs[ident] = 1
        text = ident + '\n' if '--parsable' in command else 'Submitted batch job ' + ident + '\n'
        return SimpleNamespace(returncode=0, stdout=text if kwargs.get('text') else text.encode(),
                               stderr='' if kwargs.get('text') else b'')
    monkeypatch.setattr(policy.subprocess, 'run', launch)
    return SimpleNamespace(root=tmp_path, jobs=jobs, commands=commands)


def submit(pool, name):
    cwd = pool.root/name
    cwd.mkdir(exist_ok=True)
    return policy.submit(['sbatch', '--parsable', 'job.slurm'], cwd=cwd,
                         capture_output=True, text=True, check=False, timeout=30)


def test_live_capacity_counts_retained_jobs_across_branches_and_can_pause(pool):
    pool.jobs['retained'] = 1
    assert submit(pool, 'anl').returncode == 0
    assert submit(pool, 'vrc') is None
    policy.update(pool.root, max_jobs=3)
    assert submit(pool, 'vrc').returncode == 0
    policy.update(pool.root, max_jobs=0)
    assert submit(pool, 'qc') is None
    assert len(pool.jobs) == 3  # Shrinking/pausing never cancels work.
    policy.update(pool.root, max_jobs=4)
    assert submit(pool, 'qc').returncode == 0
    assert len(pool.commands) == 3


def test_live_exclusions_override_stale_header_and_environment_without_editing(pool, monkeypatch):
    directory = pool.root/'case'; directory.mkdir()
    script = directory/'job.slurm'
    script.write_text('#!/bin/sh\n#SBATCH --exclude=old\n#SBATCH --cpus-per-task=8\n'
                      '#SBATCH --mem=104000M\npython worker.py\n')
    before = file_sha256(script)
    monkeypatch.setenv('SBATCH_EXCLUDE', 'stale')
    policy.update(pool.root, exclude=['new'], allow=['bad1'])
    submit(pool, 'case')
    command = pool.commands[-1][0]
    assert '--exclude=new' in command and '--nodes=1' in command and '--exclusive' in command
    assert '--cpus-per-task' not in ' '.join(command) and '--mem=' not in ' '.join(command)
    assert file_sha256(script) == before
    policy.exclude_failed_node(directory, 'newer')
    assert policy.read(directory)['excluded_nodes'] == ['new', 'newer']
    policy.update(pool.root, allow=['new', 'newer'])
    assert policy.read(directory)['excluded_nodes'] == []


def test_live_walltime_updates_existing_staged_scripts_without_changing_capacity_or_hashes(pool, monkeypatch):
    from kinbot import run_walltime
    monkeypatch.setattr(run_walltime, 'partition_limits', lambda **kw: {'long': 604800})
    for name in ('anl', 'vrc'):
        directory = pool.root/name; directory.mkdir()
        (directory/'job.slurm').write_text('#!/bin/bash\n#SBATCH --time=12:00:00\n'
                                         '#SBATCH --partition=long\nrun_qc\n')
    before = file_sha256(pool.root/'anl/job.slurm')
    policy.update(pool.root, minimum_walltime='7-00:00:00')
    assert policy.read(pool.root)['max_jobs'] == 2
    assert policy.read(pool.root)['excluded_nodes'] == ['bad1']
    assert submit(pool, 'anl').returncode == 0
    assert '--time=7-00:00:00' in pool.commands[-1][0]
    assert file_sha256(pool.root/'anl/job.slurm') == before
    claim = json.loads((pool.root/'.kinbot_policy_submissions/101.json').read_text())
    assert claim['walltime_override']['effective'] == '7-00:00:00'
    policy.update(pool.root, minimum_walltime='')
    assert submit(pool, 'vrc').returncode == 0
    assert not any(x.startswith('--time') for x in pool.commands[-1][0])


def test_walltime_lookup_outage_does_not_create_an_uncertain_submission(pool, monkeypatch):
    from kinbot import run_walltime
    from kinbot.anl.dispatch import SchedulerUnavailableError
    policy.update(pool.root, minimum_walltime='7-00:00:00')
    def offline(*args, **kwargs):
        raise subprocess.TimeoutExpired('sinfo', 20)
    monkeypatch.setattr(run_walltime, 'submission_override', offline)
    with pytest.raises(SchedulerUnavailableError, match='walltime'):
        submit(pool, 'anl')
    assert not pool.commands
    assert json.loads((pool.root/'.kinbot_policy_admissions.json').read_text()) == {}


def test_single_input_walltime_default_can_be_edited_live(tmp_path):
    policy.initialize(tmp_path, {'max_jobs': 3, 'minimum_walltime': '7-00:00:00'})
    policy.update(tmp_path, minimum_walltime='10-00:00:00')
    policy.initialize(tmp_path, {'max_jobs': 1, 'minimum_walltime': '1-00:00:00'})
    assert policy.read(tmp_path)['minimum_walltime'] == '10-00:00:00'
    with pytest.raises(ValueError, match='walltime'):
        policy.validate_options({'max_jobs': 3, 'minimum_walltime': 'invalid'})


def test_timeout_reserves_capacity_and_never_resubmits(pool, monkeypatch):
    policy.update(pool.root, max_jobs=1)
    def ambiguous(*args, **kwargs):
        raise subprocess.TimeoutExpired('sbatch', 30)
    monkeypatch.setattr(policy.subprocess, 'run', ambiguous)
    with pytest.raises(subprocess.TimeoutExpired):
        submit(pool, 'one')
    assert submit(pool, 'two') is None
    with pytest.raises(RuntimeError, match='unconfirmed'):
        submit(pool, 'one')
    assert len(json.loads((pool.root/'.kinbot_policy_admissions.json').read_text())) == 1


def test_scheduler_query_outage_leaves_no_submission_claim(pool, monkeypatch):
    from kinbot.anl.dispatch import SchedulerUnavailableError
    def down(root):
        raise SchedulerUnavailableError('offline')
    monkeypatch.setattr(policy, '_queue', down)
    with pytest.raises(SchedulerUnavailableError):
        submit(pool, 'one')
    assert not pool.commands
    assert not (pool.root/'.kinbot_policy_admissions.json').exists()


def test_waiting_caller_can_read_live_increase_and_scheduler_recovery(pool, monkeypatch):
    from kinbot.anl.dispatch import SchedulerUnavailableError
    policy.update(pool.root, max_jobs=0)
    calls = []
    def query(root):
        calls.append('query')
        if len(calls) == 1:
            raise SchedulerUnavailableError('offline')
        return {}
    def sleep(seconds):
        if len(calls) == 2:
            policy.update(pool.root, max_jobs=1)
    monkeypatch.setattr(policy, '_queue', query)
    monkeypatch.setattr(policy.time, 'sleep', sleep)
    result = policy.submit(['sbatch', 'qc.slurm'], cwd=pool.root, wait=True,
                           capture_output=True, text=True, check=False)
    assert result.returncode == 0 and len(pool.commands) == 1 and len(calls) == 3


@pytest.mark.parametrize('result', ['garbage', 'reject', 'exec_error'])
def test_submission_outcomes_are_not_confused_with_capacity(pool, monkeypatch, result):
    def launch(*args, **kwargs):
        if result == 'exec_error':
            raise FileNotFoundError('no sbatch')
        return SimpleNamespace(returncode=1 if result == 'reject' else 0, stdout='?', stderr='denied')
    monkeypatch.setattr(policy.subprocess, 'run', launch)
    if result == 'reject':
        assert submit(pool, 'one').returncode == 1
    else:
        with pytest.raises((FileNotFoundError, RuntimeError)):
            submit(pool, 'one')
    claims = json.loads((pool.root/'.kinbot_policy_admissions.json').read_text())
    assert len(claims) == (1 if result == 'garbage' else 0)


def test_fast_finished_job_releases_claim_only_with_terminal_evidence(pool, monkeypatch):
    submit(pool, 'one')
    pool.jobs.clear()  # Job disappeared before any subsequent squeue poll.
    policy.update(pool.root, max_jobs=1)
    assert submit(pool, 'two') is None
    monkeypatch.setattr(policy, 'slurm_job_record', lambda job: dict(JobId=job, JobState='COMPLETED'))
    assert submit(pool, 'two').returncode == 0


def test_initialize_and_update_preserve_input_and_learned_exclusions(tmp_path):
    case = tmp_path/'case'; case.mkdir()
    atomic_json(case/'input.json', {'slurm_exclude_nodes': ['old']})
    atomic_json(tmp_path/'prepared.json', {'inputs': [str(case/'input.json')]})
    atomic_json(case/'123/.kinbot_slurm_nodes.json', {'excluded_nodes': ['learned']})
    policy.update(tmp_path, max_jobs=4)
    assert policy.read(case)['excluded_nodes'] == ['learned', 'old']
    policy.initialize(case, {'max_jobs': 1, 'excluded_nodes': ['ignored_defaults']})
    assert policy.limit(1, case) == 4
    with pytest.raises(ValueError, match='ancestor'):
        policy.update(case, max_jobs=2)
    policy.update(tmp_path, allow=['old'])
    assert policy.read(case)['excluded_nodes'] == ['learned']
    assert len(list((tmp_path/'.kinbot_policy_history').glob('*.json'))) == 2


@pytest.mark.parametrize('options', [None, [], {'max_jobs': True}, {'max_jobs': -1},
                                    {'max_jobs': 2, 'unknown': 1},
                                    {'max_jobs': 2, 'excluded_nodes': ['node[1-5]']}])
def test_invalid_single_input_options_rejected(options):
    with pytest.raises(ValueError):
        policy.validate_options(options)


def test_independent_processes_share_one_capacity_and_file_lock(tmp_path):
    """Real subprocesses exercise flock across independent Python runtimes."""
    root = tmp_path/'campaign'; root.mkdir()
    policy.update(root, max_jobs=2)
    bin_dir = tmp_path/'bin'; bin_dir.mkdir()
    db = tmp_path/'scheduler.json'; db.write_text('{}')
    shim = bin_dir/'slurm_fake'
    shim.write_text(f'#!{sys.executable}\n' + '''import json,os,sys,time
from pathlib import Path
db=Path(os.environ['TEST_QUEUE'])
jobs=json.loads(db.read_text())
if Path(sys.argv[0]).name=='squeue':
    for ident,cwd in jobs.items(): print(ident+'|1|'+cwd)
elif Path(sys.argv[0]).name=='sbatch':
    time.sleep(.1)
    ident=str(700+len(jobs)); jobs[ident]=os.getcwd()
    db.write_text(json.dumps(jobs))
    print(ident)
else:
    sys.exit(1)
''')
    shim.chmod(0o755)
    for name in ('squeue', 'sbatch', 'scontrol'):
        (bin_dir/name).symlink_to(shim)
    env = dict(os.environ, PATH=str(bin_dir)+os.pathsep+os.environ['PATH'], TEST_QUEUE=str(db),
               PYTHONPATH=str(Path(__file__).resolve().parents[1]))
    children = []
    for n in range(4):
        directory = root/str(n); directory.mkdir()
        children.append(subprocess.Popen([sys.executable, '-c',
            "from kinbot.run_policy import submit; r=submit(['sbatch','--parsable','job.slurm'], "
            "capture_output=True,text=True,check=False,timeout=5); print('full' if r is None else r.stdout.strip())"],
            cwd=directory, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True))
    results = [p.communicate(timeout=30) for p in children]
    assert all(p.returncode == 0 for p in children), results
    assert sum(out.strip() == 'full' for out, err in results) == 2
    assert len(json.loads(db.read_text())) == 2


def test_dispatch_uses_shared_cap_instead_of_frozen_graph_cap(pool, monkeypatch):
    from kinbot.anl import dispatch
    from tests.anl_fixture import dispatch_spec, molpro_task
    spec = dispatch_spec(); spec['limits']['max_nodes'] = 1
    spec['tasks'] = [molpro_task('sp'+str(i), 'rhf\n', geometry_from='initial') for i in range(3)]
    graph = pool.root/'anl'; dispatch.prepare(spec, graph)
    before = file_sha256(graph/'workflow.json')
    state = dispatch.advance(graph, submit=True)
    assert [v['status'] for v in state['tasks'].values()] == ['submitted', 'submitted', 'staged']
    assert file_sha256(graph/'workflow.json') == before
    assert len(pool.commands) == 2


def test_dispatch_counting_outage_restores_staged_not_uncertain(pool, monkeypatch):
    from kinbot.anl import dispatch
    from tests.anl_fixture import dispatch_spec, molpro_task
    spec = dispatch_spec(); spec['tasks'] = [molpro_task('sp', 'rhf\n', geometry_from='initial')]
    graph = pool.root/'anl'; dispatch.prepare(spec, graph)
    monkeypatch.setattr(policy, '_queue', lambda root: (_ for _ in ()).throw(dispatch.SchedulerUnavailableError('off')))
    with pytest.raises(dispatch.SchedulerUnavailableError):
        dispatch.advance(graph, submit=True)
    assert dispatch._load(graph)[2]['tasks']['sp']['status'] == 'staged'


def test_automatic_branch_borrows_pool_slots_beyond_legacy_reservation(pool, monkeypatch):
    from kinbot.anl import automatic
    policy.update(pool.root, max_jobs=3)
    graphs = {str(i): pool.root/str(i) for i in range(3)}
    states = {n: {'tasks': {'sp': {'status': 'staged'}}} for n in graphs}
    calls = []
    def advance(graph, submit=False, submit_only=None):
        name = graph.name
        if submit:
            calls.append(name)
            states[name]['tasks']['sp']['status'] = 'complete'
        return deepcopy(states[name])
    monkeypatch.setattr(automatic, '_load', lambda graph: (None, {}, None))
    monkeypatch.setattr(automatic, '_required_geometries', lambda spec: set())
    monkeypatch.setattr(automatic, 'advance', advance)
    monkeypatch.setattr(automatic, '_retry_l3_execution_failures', lambda graph, spec, state: state)
    monkeypatch.setattr(automatic, '_retry_vpt2_framework_failures', lambda graph, state: state)
    sleeps = []
    monkeypatch.setattr(automatic.time, 'sleep', lambda seconds: sleeps.append(seconds))
    automatic._drive_graphs_locked(graphs, 1, 20)
    assert len(calls) == 3 and len(sleeps) == 1


def test_ordinary_qc_uses_same_gate_without_legacy_block(pool, monkeypatch):
    from kinbot.qc import QuantumChemistry
    monkeypatch.chdir(pool.root)
    qc = SimpleNamespace(check_qc=lambda job: 0, queue_job_limit=1,
        limit_jobs=lambda: pytest.fail('static cap must be overridden'),
        par={'queue_template': '', 'queue_walltime': '01:00:00'}, queuing='slurm', qc='fc',
        queue_name='short', slurm_feature='', job_ids={})
    assert QuantumChemistry.submit_qc(qc, 'job', 8) == 1
    assert qc.job_ids['job'] == '101'
    assert '--exclusive' in pool.commands[0][0]


def test_rotd_adapter_updates_admission_only(pool, monkeypatch):
    from kinbot.rotdpy_policy import cooperative_multi
    monkeypatch.chdir(pool.root)
    class Base:
        def __init__(self):
            self.sub_cmd = 'sbatch'
            self.calculator = {'max_jobs': 1, 'mem': 700, 'method': 'CASPT2'}
        def check_running_jobs(self):
            return 'checked'
        def submit_work(self, job):
            return job
    obj = cooperative_multi(Base)()
    assert 'kinbot.run_policy submit' in obj.sub_cmd
    policy.update(pool.root, max_jobs=9)
    assert obj.check_running_jobs() == 'checked'
    assert obj.calculator == {'max_jobs': 9, 'mem': 700, 'method': 'CASPT2'}
    policy.update(pool.root, max_jobs=0)
    assert obj.submit_work('sample') == 'sample'
    assert obj.calculator['max_jobs'] == 1  # Shared admission gates the actual pause.


def test_mess_joins_same_pool(pool, monkeypatch):
    from kinbot.mess_execution import _submit
    monkeypatch.chdir(pool.root)
    assert _submit('slurm', 'me/run.slurm') == '101'
    assert len(pool.commands) == 1 and '--exclusive' in pool.commands[0][0]


def test_nested_policy_fails_closed_even_if_created_outside_cli(pool):
    directory = pool.root/'child'; directory.mkdir()
    atomic_json(directory/policy.NAME, policy.read(pool.root))
    with pytest.raises(ValueError, match='Nested shared policies'):
        policy.submit(['sbatch', 'job.slurm'], cwd=directory)
    assert not pool.commands


def test_cli_live_edits_and_initial_single_input_defaults(tmp_path, capsys):
    policy.initialize(tmp_path, {'max_jobs': 3, 'excluded_nodes': ['bad']})
    policy.main(['set', str(tmp_path), '--max-jobs', '5', '--exclude', 'spare', '--allow', 'bad'])
    assert policy.read(tmp_path)['max_jobs'] == 5
    assert policy.read(tmp_path)['excluded_nodes'] == ['spare']
    assert 'Running and pending jobs were not modified' in capsys.readouterr().out
