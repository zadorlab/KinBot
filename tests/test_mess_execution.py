from pathlib import Path
from types import SimpleNamespace
import subprocess

import pytest

from kinbot import mess_execution as execution
from kinbot import kb_path


def test_pes_writes_all_uq_inputs_before_starting_mess(tmp_path, monkeypatch):
    import json
    import logging
    from unittest.mock import patch
    from kinbot import pes
    from kinbot.mess import MESS
    from kinbot.parameters import Parameters
    from tests.test_mess_conformers import point

    monkeypatch.chdir(tmp_path)
    Path('123').mkdir()
    Path('input.json').write_text(json.dumps(dict(barrier_threshold=100., smiles='O',
        multi_conf_tst=1, conformer_search=1, high_level=0, me=0, epsilon=100., sigma=3.)))
    par = Parameters('input.json', show_warnings=False).par
    par.update(pes=1, me=1, uq_n=2)
    species = point('water', 123)
    renderer = MESS(par, species)
    for index in range(2):
        Path(f'123/123_{index:04d}.mess').write_text(
            renderer.write_well(species, 0., 1., index))
    calls = []
    def run(writer):
        calls.append(writer)
        assert all(Path(f'me/mess_{index:04d}.inp').exists() for index in range(2))
    with patch.object(MESS, 'run', run), patch.object(pes, 'logger', logging.getLogger('KinBot'), create=True):
        pes.create_mess_input(par, ['123'], [], [], [], [], {'123': 0.}, {}, {'123': '123'}, 18., False)
    assert len(calls) == 1


def writer(tmp_path, monkeypatch, queue='slurm', count=1, limit=2, run=True):
    monkeypatch.chdir(tmp_path)
    Path('me').mkdir()
    return SimpleNamespace(par=dict(queuing=queue, uq_n=count, uq_max_runs=limit,
                                    run_me=run),
                           write_submitscript=lambda path, index: Path(path).write_text('script'))


def complete(index, status=0, output=True):
    Path(f'me/mess_{index:04d}.exitcode').write_text(str(status))
    if output:
        Path(f'me/mess_{index:04d}.out').write_text('calculated rates\n')


def test_queue_waits_throttles_and_drains(tmp_path, monkeypatch):
    obj = writer(tmp_path, monkeypatch, count=3, limit=2)
    submitted, states, sleeps = [], {}, []
    def submit(queue, path):
        index = len(submitted)
        # The first two can overlap; the third must wait for their completion.
        if index == 2:
            assert states['0'] == states['1'] == 2
        submitted.append(index)
        return str(index)
    def active(queue, pid):
        states[pid] = states.get(pid, 0) + 1
        if states[pid] == 1:
            return True
        complete(int(pid))
        return False
    monkeypatch.setattr(execution, '_submit', submit)
    monkeypatch.setattr(execution, '_active', active)
    monkeypatch.setattr(execution.time, 'sleep', lambda delay: sleeps.append(delay))
    assert execution.run_mess(obj) == 0
    assert states == {'0': 2, '1': 2, '2': 2}
    assert len(sleeps) == 4


@pytest.mark.parametrize('status,output,message', [(134, False, 'exit code 134'),
                                                   (0, False, 'no new, nonempty')])
def test_solver_failure_is_not_kinbot_success(tmp_path, monkeypatch, status, output, message):
    obj = writer(tmp_path, monkeypatch)
    monkeypatch.setattr(execution, '_submit', lambda *args: '12')
    def active(*args):
        complete(0, status, output)
        return False
    monkeypatch.setattr(execution, '_active', active)
    monkeypatch.setattr(execution.time, 'sleep', lambda delay: None)
    with pytest.raises(execution.MESSExecutionError, match=message):
        execution.run_mess(obj)


def test_cancelled_job_cannot_reuse_stale_success(tmp_path, monkeypatch):
    obj = writer(tmp_path, monkeypatch)
    complete(0)
    monkeypatch.setattr(execution, '_submit', lambda *args: '12')
    monkeypatch.setattr(execution, '_active', lambda *args: False)
    monkeypatch.setattr(execution.time, 'sleep', lambda delay: None)
    with pytest.raises(execution.MESSExecutionError, match='without a solver exit result'):
        execution.run_mess(obj)


def test_stale_rate_file_is_not_success(tmp_path, monkeypatch):
    writer(tmp_path, monkeypatch)
    complete(0)
    stamp = execution._stamp(Path('me/mess_0000.out'))
    with pytest.raises(execution.MESSExecutionError, match='no new, nonempty'):
        execution._verify(0, stamp)


@pytest.mark.parametrize('queue,command', [('pbs', 'qsub'), ('slurm', 'sbatch'), ('local', 'sh')])
def test_input_only_writes_correct_manual_commands(tmp_path, monkeypatch, queue, command):
    obj = writer(tmp_path, monkeypatch, queue=queue, run=False)
    monkeypatch.setattr(execution, '_submit', lambda *args: pytest.fail('must not submit'))
    assert execution.run_mess(obj) == 0
    assert Path('batch_me.sub').read_text().startswith(command + ' me/run_mess_0000')


@pytest.mark.parametrize('queue,output,active', [('slurm', 'RUNNING\n', True),
    ('slurm', 'COMPLETED\n', False), ('slurm', '', False),
    ('pbs', 'job_state = R\n', True), ('pbs', 'job_state = F\n', False)])
def test_scheduler_states(monkeypatch, queue, output, active):
    monkeypatch.setattr(execution, '_run', lambda command:
                        subprocess.CompletedProcess(command, 0, output, ''))
    assert execution._active(queue, '123') is active


def test_scheduler_outage_and_submission_failure_raise(monkeypatch):
    monkeypatch.setattr(execution, '_run', lambda command:
                        subprocess.CompletedProcess(command, 1, '', 'controller unavailable'))
    with pytest.raises(execution.MESSExecutionError, match='Cannot query'):
        execution._active('slurm', '123')
    with pytest.raises(execution.MESSExecutionError, match='submission failed'):
        execution._submit('slurm', 'script')


def test_submission_failure_still_drains_existing_jobs(tmp_path, monkeypatch):
    obj = writer(tmp_path, monkeypatch, count=3, limit=2)
    submitted, polled = [], []
    def submit(queue, path):
        submitted.append(path)
        if len(submitted) == 2:
            raise execution.MESSExecutionError('submission failed')
        return '123'
    def active(queue, pid):
        polled.append(pid)
        complete(0)
        return False
    monkeypatch.setattr(execution, '_submit', submit)
    monkeypatch.setattr(execution, '_active', active)
    monkeypatch.setattr(execution.time, 'sleep', lambda delay: None)
    with pytest.raises(execution.MESSExecutionError, match='submission failed'):
        execution.run_mess(obj)
    assert len(submitted) == 2
    assert polled == ['123']


@pytest.mark.parametrize('queue', ['pbs', 'slurm', 'local'])
@pytest.mark.parametrize('status', [0, 7])
def test_actual_shell_exit_and_output(tmp_path, monkeypatch, queue, status):
    obj = writer(tmp_path, monkeypatch, queue=queue)
    binary = tmp_path / 'mess'
    binary.write_text('#!/bin/sh\nout="${1%.inp}.out"\n'
                      'printf "calculated rates\\n" > "$out"\n'
                      f'exit {status}\n')
    binary.chmod(0o700)
    monkeypatch.setenv('PATH', f'{tmp_path}:/usr/bin:/bin')
    if queue == 'local':
        if status:
            with pytest.raises(execution.MESSExecutionError, match=f'exit code {status}'):
                execution.run_mess(obj)
        else:
            assert execution.run_mess(obj) == 0
    else:
        template = Path(kb_path, 'tpl', f'{queue}_mess_uq.tpl').read_text().format(n='0000')
        monkeypatch.setenv('PBS_O_WORKDIR', str(tmp_path))
        result = subprocess.run(['sh', '-c', template], check=False)
        assert result.returncode == status
        assert Path('me/mess_0000.exitcode').read_text().strip() == str(status)
        assert Path('me/mess_0000.out').stat().st_size


@pytest.mark.parametrize('queue', ['local', 'slurm', 'pbs'])
@pytest.mark.parametrize('custom_header', [False, True])
@pytest.mark.parametrize('fail_extra', [False, True])
def test_grouped_jobs_use_their_own_files_and_detect_extra_failure(
        tmp_path, monkeypatch, queue, custom_header, fail_extra):
    from kinbot.mess import MESS
    obj = writer(tmp_path, monkeypatch, queue=queue)
    obj.par.update(queue_template='', ppn=1, queue_name='test', slurm_feature='')
    if custom_header:
        Path('custom.tpl').write_text('#!/bin/sh\n# {name}\n')
        obj.par['queue_template'] = str(tmp_path / 'custom.tpl')
    obj.write_submitscript = lambda path, index: MESS.write_submitscript(obj, path, index)
    obj.mess_jobs = [dict(stem='mess_0000', status='ready'),
                     dict(stem='mess_0000_group_0002', status='ready')]
    binary = tmp_path / 'mess'
    binary.write_text('#!/bin/sh\nout="${1%.inp}.out"\n'
        'printf "rates\\n" > "$out"\npwd > "${1%.inp}.micro"\n'
        + ('case "$1" in *group*) exit 7;; esac\n' if fail_extra else '') + 'exit 0\n')
    binary.chmod(0o700)
    monkeypatch.setenv('PATH', f'{tmp_path}:/usr/bin:/bin')
    monkeypatch.setenv('PBS_O_WORKDIR', str(tmp_path))
    monkeypatch.setenv('SLURM_SUBMIT_DIR', str(tmp_path))
    submitted = []
    def submit(queue, script):
        subprocess.run(['sh', str(script)], check=False)
        submitted.append(script)
        return str(len(submitted))
    monkeypatch.setattr(execution, '_submit', submit)
    monkeypatch.setattr(execution, '_active', lambda *args: False)
    monkeypatch.setattr(execution.time, 'sleep', lambda *args: None)
    if fail_extra:
        with pytest.raises(execution.MESSExecutionError, match='0000_group_0002 failed with exit code 7'):
            execution.run_mess(obj)
    else:
        assert execution.run_mess(obj) == 0
    for job in obj.mess_jobs:
        assert Path('me', job['stem'] + '.out').exists()
        assert Path('me', job['stem'] + '.micro').read_text().strip() == str(tmp_path / 'me')
    assert Path('batch_me.sub').read_text().count('run_mess_') == 2


def test_primary_without_reactions_does_not_hide_ready_secondary(tmp_path, monkeypatch):
    obj = writer(tmp_path, monkeypatch, run=False)
    obj.mess_jobs = [dict(stem='mess_0000', status='no_reactions'),
                     dict(stem='mess_0000_group_0002', status='ready')]
    execution.run_mess(obj)
    batch = Path('batch_me.sub').read_text()
    assert batch.count('run_mess_') == 1
    assert 'run_mess_0000_group_0002' in batch
    # An explicitly empty current list must not reuse any existing inputs.
    obj.mess_jobs = []
    execution.run_mess(obj)
    assert Path('batch_me.sub').read_text() == ''
