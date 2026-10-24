"""Run the actual shell handoff across subprocesses with simulated site owners.

Real control publication, decision application, provenance checks, policy files
and shell ordering are exercised. Scheduler/process/native QC boundaries are
simulated; controller quiescence itself is covered by test_run_refresh.
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from kinbot import live_control, rotdpy, rotdpy_campaign, run_policy
from kinbot.execution import atomic_json, file_sha256
from tests.test_pes_selection import evidence


RUNNER = r'''
import importlib, json, os, sys, time
from pathlib import Path
sys.path.insert(0, os.environ['KB_TEST_SOURCE'])
from kinbot import live_control, reaction_decisions, pes_phase
from kinbot import run_refresh, run_control, run_monitor, rotdpy, rotdpy_campaign
from kinbot.execution import atomic_json
root = Path(os.environ['KB_RUN'])
def event(name):
    with (root/'events.jsonl').open('a') as f:
        f.write(json.dumps(name)+'\n')
original_iterdir = Path.iterdir
Path.iterdir = lambda p: iter(()) if str(p) == '/proc' else original_iterdir(p)
time.sleep = lambda seconds: None
pes_phase.active_worker = lambda p: (root/'controller_refresh.json').exists()
pes_phase.owner_alive = lambda p: p.get('fixture_alive', False)
def refresh(args):
    assert args == [str(root), '--apply', '--wait-idle', '1800']
    event('refresh')
    atomic_json(root/'controller_refresh.json', dict(phase='launched', code=str(root/'.salvage_code.test')))
    for case in ('first', 'second'):
        atomic_json(root/case/'pes_status.json', dict(status='running'))
        atomic_json(root/case/'pes_exploration.json', dict(phase='postprocessing', fixture_alive=True))
run_refresh.main = refresh
original_request = live_control.request
def request(*args, **kwargs):
    assert os.environ['PYTHONPATH'] == str(root/'.salvage_code.test')
    assert (root/'controller_refresh.json').exists()
    record = original_request(*args, **kwargs)
    path = root/'.kinbot_controls/requests'/(record['id']+'.json')
    # Simulate the replacement worker reaching its network barrier.
    assert reaction_decisions.apply(Path(record['target']), record['options']['reaction_decisions'], path, record)
    live_control.acknowledge(path, record, 'applied')
    event('decision_applied')
    return record
live_control.request = request
rotdpy_campaign.active = lambda inp, record: record.get('fixture_active', False)
original_spawn = rotdpy_campaign.subprocess.Popen
def spawn(args, **kwargs):
    if args[1:5] != ['-u', '-m', 'kinbot.rotdpy', 'run']:
        return original_spawn(args, **kwargs)
    assert reaction_decisions.rejected(root/'first/123', '123_channel')
    inp = Path(args[-1])
    record = json.loads(rotdpy_campaign.receipt(inp).read_text())
    assert not record.get('operator_rejected') and not record.get('fixture_active')
    # The real launcher obtains this identity and records it after spawning.
    from types import SimpleNamespace
    event('resume:'+inp.stem)
    return SimpleNamespace(pid=1234)
rotdpy_campaign.subprocess.Popen = spawn
pes_phase.process_identity = lambda pid: {'pid':1234} if pid == 1234 else None
original_active = rotdpy_campaign.active
rotdpy_campaign.active = lambda inp, record: (original_active(inp, record)
    or record.get('owner_identity') == {'pid':1234})
rotdpy.read_result = lambda inp: event('completed:'+inp.stem)
def reconcile(args):
    assert args == ['reconcile', str(root), '--apply']
    assert reaction_decisions.rejected(root/'first/123', '123_channel')
    event('reconcile')
    print(json.dumps({'graphs':[{'status':'reconciled','changes':['fixture result adopted']}]}))
run_control.main = reconcile
run_monitor.main = lambda args: event('monitor')
args = sys.argv[1:]
while args and args[0] in ('-B', '-u'):
    args.pop(0)
if args[0] == '-m':
    importlib.import_module(args[1]).main(args[2:])
else:
    source = sys.stdin.read() if args[0] == '-' else args[1]
    sys.argv = args if args[0] == '-' else ['-c', *args[2:]]
    exec(compile(source, '<handoff>', 'exec'), {'__name__':'__main__'})
'''


@pytest.fixture
def handoff(tmp_path):
    source = Path(__file__).resolve().parents[1]
    checkout = tmp_path/'checkout'
    script = checkout/'examples/anl/tools/refresh_retained_run.sh'
    script.parent.mkdir(parents=True)
    shutil.copyfile(source/'examples/anl/tools/refresh_retained_run.sh', script)
    python = checkout/'.venv/bin/python'
    python.parent.mkdir(parents=True)
    python.write_text('#!'+sys.executable+'\n'+RUNNER)
    python.chmod(0o755)
    bin_dir = tmp_path/'bin'
    bin_dir.mkdir()
    # flock is absent on some local development hosts. The inherited descriptor
    # stays open in bash, matching the real Linux utility's lock lifetime.
    flock = bin_dir/'flock'
    flock.write_text('#!'+sys.executable+'\nimport fcntl,sys\nfcntl.flock(int(sys.argv[-1]), fcntl.LOCK_EX|fcntl.LOCK_NB)\n')
    flock.chmod(0o755)
    scontrol = bin_dir/'scontrol'
    scontrol.write_text('#!/bin/sh\nprintf "%s\\n" "NodeName=node1 State=DOWN" "NodeName=node2 State=IDLE" "NodeName=node3 State=DRAIN" "NodeName=node4 State=UNKNOWN"\n')
    scontrol.chmod(0o755)
    root = tmp_path/'run'
    inputs = []
    for case in ('first', 'second'):
        worker = root/case/'123'
        worker.mkdir(parents=True)
        atomic_json(worker/'123.json', {})
        atomic_json(worker.parent/'input.json', {})
        inputs.append(str(worker.parent/'input.json'))
    pin = root/'.salvage_code.test'
    (pin/'kinbot').mkdir(parents=True)
    (pin/'kinbot/rotdpy_worker.py').write_text('# verified worker snapshot')
    atomic_json(pin/'manifest.json', {'kinbot/rotdpy_worker.py':file_sha256(pin/'kinbot/rotdpy_worker.py')})
    atomic_json(root/'prepared.json', dict(inputs=inputs))
    atomic_json(root/'first/123/accepted_l2_network.json', evidence('123', [('123_channel',['p','q'],0.,None)]))
    run_policy.update(root, max_jobs=23, exclude=['node9'])
    live_control.configure(root)
    decisions = root/'decisions.json'
    atomic_json(decisions, [dict(worker='first/123', reaction='123_channel', decision='reject',
        reason='Reviewed relaxed geometry', cancel_running=True)])
    sampling = root/'first/rotdPy'
    sampling.mkdir()
    (sampling/'site_setup.sh').write_text('# Generated from the environment visible during prepare.\n')
    (sampling/'qu.tpl').write_text('source '+str(sampling/'site_setup.sh')+'\npython -m kinbot.rotdpy_worker\n')
    retained = []
    for state in ('active', 'absent', 'complete', 'rejected'):
        inp = sampling/(state+'.py')
        inp.write_text('validation_metadata = {}\n')
        record = dict(input=str(inp), **rotdpy_campaign.contract(inp), fixture_active=state=='active')
        if state == 'rejected':
            record['operator_rejected'] = {'request':'previous'}
        atomic_json(rotdpy_campaign.receipt(inp), record)
        db = sampling/('kb_'+state)/'samples.db'
        db.parent.mkdir()
        db.write_bytes(b'completed samples retained')
        retained.extend((inp, db))
        if state == 'complete':
            atomic_json(rotdpy.execution_path(inp), dict(status='complete'))
    native = root/'first/123/kinbot.db'
    native.write_bytes(b'accepted chemistry retained')
    retained.append(native)
    env = dict(os.environ, PATH=str(bin_dir)+os.pathsep+os.environ['PATH'],
        KB_DECISIONS_FILE=str(decisions), KB_TEST_SOURCE=str(source))
    env.pop('KB_WAIT_PID', None)
    env.pop('KB_WAIT_LOG', None)
    def run():
        return subprocess.run(['bash', str(script), str(root), '22'], env=env,
            capture_output=True, text=True, timeout=30)
    run.env = env
    return root, decisions, {p:p.read_bytes() for p in retained}, run


def test_complete_shell_handoff_applies_decisions_then_adopts_and_reconciles(handoff):
    root, _, retained, run = handoff
    result = run()
    assert result.returncode == 0, result.stdout+result.stderr
    events = [json.loads(s) for s in (root/'events.jsonl').read_text().splitlines()]
    assert events == ['refresh', 'decision_applied', 'resume:absent', 'completed:complete', 'reconcile', 'monitor']
    policy = json.loads((root/run_policy.NAME).read_text())
    assert policy['max_jobs'] == 22
    assert policy['excluded_nodes'] == ['node1','node3','node4','node9']
    assert 'ROTD OWNER RETAINED active.py' in result.stdout
    assert 'ROTD CHANNEL EXCLUDED rejected.py' in result.stdout
    assert result.stdout.count('owner_alive True') == 4
    assert 'Handoff finished.' in result.stdout
    assert all(p.read_bytes() == old for p,old in retained.items())
    result = run()
    assert result.returncode == 0, result.stdout+result.stderr
    events = [json.loads(s) for s in (root/'events.jsonl').read_text().splitlines()]
    assert events.count('resume:absent') == events.count('decision_applied') == 1
    assert len(list((root/'.kinbot_controls/requests').glob('*.json'))) == 1


@pytest.mark.parametrize('fault', ['unknown_reaction', 'duplicate', 'invalid_reason'])
def test_invalid_decision_batch_stops_before_policy_or_controller_changes(handoff, fault):
    root, path, retained, run = handoff
    items = json.loads(path.read_text())
    if fault == 'unknown_reaction':
        items[0]['reaction'] = 'not_accepted'
    elif fault == 'duplicate':
        items.append(dict(items[0], decision='include'))
    else:
        items[0]['reason'] = ''
    atomic_json(path, items)
    before = (root/run_policy.NAME).read_bytes()
    result = run()
    assert result.returncode != 0
    assert (root/run_policy.NAME).read_bytes() == before
    assert not (root/'events.jsonl').exists()
    assert all(p.read_bytes() == old for p,old in retained.items())


def test_invalid_sampling_setup_stops_before_controller_replacement(handoff):
    root, _, retained, run = handoff
    (root/'first/rotdPy/qu.tpl').write_text('python -m custom_worker\n')
    before = (root/run_policy.NAME).read_bytes()
    result = run()
    assert result.returncode != 0 and 'generated setup' in result.stderr
    assert not (root/'events.jsonl').exists()
    assert (root/run_policy.NAME).read_bytes() == before
    assert all(p.read_bytes() == old for p,old in retained.items())


@pytest.mark.parametrize('state', ['running', 'needs_review'])
def test_scan_hold_preserves_campaigns_and_allows_full_handoff_and_later_resume(handoff, state):
    from kinbot import vrc_scan_recovery as scans
    root, _, retained, run = handoff
    native = root/'first/123/native.xyz'
    native.write_text('unchanged accepted geometry')
    evidence = dict(context_sha256='context', points={'pt':dict(status='accepted',
        hashes={str(native):file_sha256(native)})})
    reviews = []
    for name in ('active', 'absent', 'complete', 'rejected'):
        inp = root/'first/rotdPy'/(name+'.py')
        correction = root/'first/123/vrctst'/('corr_'+name+'.json')
        atomic_json(correction, {})
        inp.write_text('validation_metadata = '+repr(dict(name='production', correction={
            'source':str(correction), 'source_sha256':file_sha256(correction)}))+'\n')
        path = rotdpy_campaign.receipt(inp)
        record = json.loads(path.read_text())
        record.update(rotdpy_campaign.contract(inp))
        atomic_json(path, record)
        review = correction.parent/'scan_review'/name/'scan_review.json'
        atomic_json(review, dict(schema=1, kind='scan_connectivity', reaction=name,
            **evidence, evidence_sha256=scans.digest(evidence), complete_grid=state!='running',
            status=state))
        atomic_json(correction.with_suffix('.scan.json'), dict(schema=1,
            scan_evidence_sha256=scans.digest(evidence), correction_sha256=file_sha256(correction)))
        reviews.append(review)
        retained[inp] = inp.read_bytes()
    result = run()
    assert result.returncode == 0, result.stdout+result.stderr
    events = [json.loads(s) for s in (root/'events.jsonl').read_text().splitlines()]
    assert events == ['refresh', 'decision_applied', 'reconcile', 'monitor']
    assert 'ROTD OWNER RETAINED active.py' in result.stdout
    assert 'ROTD CHANNEL EXCLUDED rejected.py' in result.stdout
    assert 'ROTD HELD: SCAN NOT READY absent.py' in result.stdout
    assert 'ROTD HELD: SCAN NOT READY complete.py' in result.stdout
    assert 'Handoff finished.' in result.stdout
    receipt = rotdpy_campaign.receipt(root/'first/rotdPy/absent.py')
    pending = json.loads(receipt.read_text())
    assert pending['resume_requested'] and pending['resume_hold']
    result = run()
    assert result.returncode == 0, result.stdout+result.stderr
    assert json.loads(receipt.read_text()) == pending
    for review in reviews:
        record = json.loads(review.read_text())
        record.update(status='resolved', complete_grid=True)
        atomic_json(review, record)
    result = run()
    assert result.returncode == 0, result.stdout+result.stderr
    events = [json.loads(s) for s in (root/'events.jsonl').read_text().splitlines()]
    assert events.count('resume:absent') == events.count('completed:complete') == 1
    assert not json.loads(receipt.read_text()).get('resume_requested')
    assert all(p.read_bytes() == old for p,old in retained.items())
