from pathlib import Path
from types import SimpleNamespace
import shlex
import subprocess

import pytest

from kinbot import mess_execution as execution


@pytest.fixture(autouse=True)
def clock(monkeypatch):
    clock = SimpleNamespace(now=0., sleeps=[])
    def tick(delay):
        clock.now += delay
        clock.sleeps.append(delay)
    clock.tick = tick
    monkeypatch.setattr(execution.time, 'sleep', tick)
    monkeypatch.setattr(execution.time, 'monotonic', lambda: clock.now)
    return clock


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


def test_queue_waits_throttles_and_drains(tmp_path, monkeypatch, clock):
    obj = writer(tmp_path, monkeypatch, count=3, limit=2)
    submitted, states = [], {}
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
    assert execution.run_mess(obj) == 0
    assert states == {'0': 2, '1': 2, '2': 2}
    assert len(clock.sleeps) == 4


@pytest.mark.parametrize('status,output,message', [(134, False, 'exit code 134'),
                                                   (0, False, 'no new, nonempty')])
def test_solver_failure_is_not_kinbot_success(tmp_path, monkeypatch, clock, status, output, message):
    obj = writer(tmp_path, monkeypatch)
    monkeypatch.setattr(execution, '_submit', lambda *args: '12')
    def active(*args):
        complete(0, status, output)
        return False
    monkeypatch.setattr(execution, '_active', active)
    with pytest.raises(execution.MESSExecutionError, match=message):
        execution.run_mess(obj)
    assert clock.now == (5. if status else 65.)


def test_cancelled_job_cannot_reuse_stale_success(tmp_path, monkeypatch, clock):
    obj = writer(tmp_path, monkeypatch)
    complete(0)
    monkeypatch.setattr(execution, '_submit', lambda *args: '12')
    monkeypatch.setattr(execution, '_active', lambda *args: False)
    with pytest.raises(execution.MESSExecutionError, match='without a solver exit result'):
        execution.run_mess(obj)
    assert clock.now == 65.


def test_stale_rate_file_is_not_success(tmp_path, monkeypatch):
    writer(tmp_path, monkeypatch)
    complete(0)
    stamp = execution._stamp(Path('me/mess_0000.out'))
    with pytest.raises(execution.MESSExecutionError, match='no new, nonempty'):
        execution._verify(0, stamp)


@pytest.mark.parametrize('queue', ['slurm', 'pbs'])
def test_queued_result_grace_waits_for_marker_and_fresh_output(tmp_path, monkeypatch, clock, queue):
    obj = writer(tmp_path, monkeypatch, queue=queue, count=3, limit=2)
    complete(0)  # Neither this exit marker nor this rate file can prove success.
    submitted, queried = [], []
    def submit(*args):
        pid = str(len(submitted))
        submitted.append(clock.now)
        return pid
    def active(queue, pid):
        queried.append(pid)
        if pid != '0':
            complete(int(pid))
        return False
    def tick(delay):
        clock.tick(delay)
        if clock.now == 10.:
            Path('me/mess_0000.exitcode').write_text('0')
        if clock.now == 20.:
            Path('me/mess_0000.out').write_text('new calculated rate coefficients\n')
    monkeypatch.setattr(execution, '_submit', submit)
    monkeypatch.setattr(execution, '_active', active)
    monkeypatch.setattr(execution.time, 'sleep', tick)
    assert execution.run_mess(obj) == 0
    assert clock.now == 20.
    assert submitted == [0., 0., 5.]
    assert queried == ['0', '1', '2']


@pytest.mark.parametrize('persistent', [False, True])
def test_scheduler_query_retries_are_bounded_and_do_not_resubmit(
        tmp_path, monkeypatch, clock, persistent, caplog):
    obj = writer(tmp_path, monkeypatch, count=2, limit=1)
    submitted, queried = [], []
    states = iter([None, True, None, None, False])
    def submit(*args):
        pid = str(len(submitted))
        submitted.append(pid)
        return pid
    def active(queue, pid):
        queried.append(pid)
        state = None if persistent else next(states) if pid == '0' else False
        if state is None:
            raise execution.MESSExecutionError('scheduler controller unavailable')
        if not state:
            complete(int(pid))
        return state
    monkeypatch.setattr(execution, '_submit', submit)
    monkeypatch.setattr(execution, '_active', active)
    if persistent:
        with pytest.raises(execution.MESSExecutionError, match='3 consecutive.*still be active: 0'):
            execution.run_mess(obj)
        assert queried == ['0'] * 3
        assert submitted == ['0']
        assert clock.now == 15.
    else:
        assert execution.run_mess(obj) == 0
        # A successful status query resets the consecutive-failure count.
        assert queried == ['0'] * 5 + ['1']
        assert submitted == ['0', '1']
    assert 'Retrying scheduler query' in caplog.text


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
    monkeypatch.setattr(execution, '_run', lambda command, **kwargs:
                        subprocess.CompletedProcess(command, 0, output, ''))
    assert execution._active(queue, '123') is active


def test_scheduler_outage_and_submission_failure_raise(monkeypatch):
    monkeypatch.setattr(execution, '_run', lambda command, **kwargs:
                        subprocess.CompletedProcess(command, 1, '', 'controller unavailable'))
    with pytest.raises(execution.MESSExecutionError, match='Cannot query'):
        execution._active('slurm', '123')
    with pytest.raises(execution.MESSExecutionError, match='submission failed'):
        execution._submit('slurm', 'script')
    def timeout(command, **kwargs):
        assert kwargs['timeout'] == 30.
        raise subprocess.TimeoutExpired(command, kwargs['timeout'])
    monkeypatch.setattr(execution, '_run', timeout)
    for queue in ('slurm', 'pbs'):
        with pytest.raises(execution.MESSExecutionError, match='query timed out'):
            execution._active(queue, '123')
    def launch_error(command, **kwargs):
        raise OSError('temporarily unavailable')
    monkeypatch.setattr(execution, '_run', launch_error)
    for queue in ('slurm', 'pbs'):
        with pytest.raises(execution.MESSExecutionError, match='Cannot query.*temporarily unavailable'):
            execution._active(queue, '123')


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
    with pytest.raises(execution.MESSExecutionError, match='submission failed'):
        execution.run_mess(obj)
    assert len(submitted) == 2
    assert polled == ['123']


@pytest.mark.parametrize('queue', ['pbs', 'slurm', 'local'])
@pytest.mark.parametrize('status', [0, 7])
def test_actual_shell_exit_and_output(tmp_path, monkeypatch, queue, status):
    from kinbot.mess import MESS

    obj = writer(tmp_path, monkeypatch, queue=queue)
    binary = tmp_path / 'custom mess'
    argument = 'literal $HOME $(touch expanded); "quoted" argument'
    binary.write_text('#!/bin/sh\nprintf "%s" "$1" > received_argument\n'
                      'shift\nout="${1%.inp}.out"\n'
                      'printf "calculated rates\\n" > "$out"\n'
                      f'exit {status}\n')
    binary.chmod(0o700)
    obj.par.update(mess_command=shlex.join([str(binary), argument]), queue_template='',
                   ppn=1, queue_name='test', slurm_feature='')
    obj.write_submitscript = lambda path, index: MESS.write_submitscript(obj, path, index)
    monkeypatch.setenv('PBS_O_WORKDIR', str(tmp_path))
    monkeypatch.setenv('SLURM_SUBMIT_DIR', str(tmp_path))
    def submit(queue, script):
        result = subprocess.run(['sh', str(script)], check=False)
        assert result.returncode == status
        return '123'
    monkeypatch.setattr(execution, '_submit', submit)
    monkeypatch.setattr(execution, '_active', lambda *args: False)
    if status:
        with pytest.raises(execution.MESSExecutionError, match=f'exit code {status}'):
            execution.run_mess(obj)
    else:
        assert execution.run_mess(obj) == 0
    assert Path('me/received_argument').read_text() == argument
    assert not Path('me/expanded').exists()
    assert Path('me/mess_0000.out').stat().st_size
    if queue == 'local':
        # The generated manual script must run the same command as the API.
        assert subprocess.run(['sh', 'me/run_mess_0000.sh'], check=False).returncode == status
    else:
        assert Path('me/mess_0000.exitcode').read_text().strip() == str(status)
    obj.par.update(mess_command='/unavailable/compute-node-only/mess', run_me=False)
    assert execution.run_mess(obj) == 0  # Input writing needs no local MESS binary.
    if queue == 'local':
        obj.par['run_me'] = True
        with pytest.raises(execution.MESSExecutionError, match='Cannot start MESS.*mess_command'):
            execution.run_mess(obj)


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


@pytest.mark.parametrize('failure', ['routing', 'solver', 'unexpected'])
def test_cli_preserves_search_outputs_and_reports_failure(
        tmp_path, monkeypatch, caplog, capsys, failure):
    import json
    import logging
    import sys
    from ase.build import molecule
    from kinbot import kb
    from kinbot.stereo_routing import StereoRoutingError

    monkeypatch.chdir(tmp_path)
    atoms = molecule('H2O')
    structure = [value for atom, xyz in zip(atoms.get_chemical_symbols(), atoms.positions)
                 for value in (atom, *xyz)]
    Path('input.json').write_text(json.dumps(dict(structure=structure, barrier_threshold=100.,
        reaction_search=0, high_level=1, conformer_search=1, rotor_scan=1,
        me=1, run_me=1, uq=0, epsilon=100., sigma=3., queuing='local', do_clean=0)))
    monkeypatch.setattr(sys, 'argv', ['kinbot', 'input.json'])
    monkeypatch.setattr(kb, 'config_log', lambda *a, **k: logging.getLogger('KinBot'))
    events = []
    def optimize(*args, **kwargs):
        if failure == 'routing':
            raise StereoRoutingError('contradictory charge; saved in unsupported_observations/result.json')
        if failure == 'unexpected':
            raise ValueError('unrelated error')
    qc = SimpleNamespace(qc='nn_pes', qc_opt=optimize,
        get_qc_geom=lambda *a, **k: (0, atoms.positions.copy()),
        get_qc_freq=lambda *a, **k: (0, [1000., 2000., 3000.]),
        get_qc_energy=lambda *a: (0, -76.), get_qc_zpe=lambda *a: (0, .01))
    monkeypatch.setattr(kb, 'QuantumChemistry', lambda *a: qc)
    monkeypatch.setattr(kb, 'Optimize', lambda *a, **k:
                        SimpleNamespace(do_optimization=lambda: None, shigh=1))
    def save(name):
        def write(*args):
            events.append(name)
            Path(name).write_text('completed reaction-search output\n')
        return write
    def run():
        events.append('solver')
        raise execution.MESSExecutionError('MESS failed with exit code 134')
    monkeypatch.setattr(kb, 'MESS', lambda *a: SimpleNamespace(write_input=save('input'), run=run))
    for function, name in (('create_summary_file', 'summary'),
                           ('createPESViewerInput', 'pesviewer'), ('creatMLInput', 'ml')):
        monkeypatch.setattr(kb.postprocess, function, save(name))
    monkeypatch.setattr(kb, 'clean_files', lambda **k: events.append('cleanup'))
    expected = {'routing': SystemExit, 'solver': execution.MESSExecutionError,
                'unexpected': ValueError}[failure]
    with pytest.raises(expected) as error:
        kb.main()
    if failure == 'routing':
        assert error.value.code == 1
        assert 'contradictory charge; saved in unsupported_observations/result.json' in caplog.text
    if failure == 'solver':
        assert events == ['input', 'summary', 'pesviewer', 'ml', 'solver']
        assert all(Path(name).exists() for name in ('summary', 'pesviewer', 'ml'))
    else:
        assert events == []
    assert 'KinBot finished.' not in caplog.text
    assert 'Done!' not in capsys.readouterr().out
