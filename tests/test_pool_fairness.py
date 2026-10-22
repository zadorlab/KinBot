import json

from kinbot import pool_fairness as fair
from kinbot.execution import atomic_json


def test_busy_anl_cannot_take_waiting_search_or_vrc_turns(tmp_path, monkeypatch):
    clock = [1000.]
    monkeypatch.setattr(fair.time, 'time', lambda:clock[0])
    for case in ('a','b'):
        atomic_json(tmp_path/case/'pes_exploration.json', {})
    a=tmp_path/'a/111/anl_automatic/states/x/tasks/energy'
    b=tmp_path/'b/222'
    v=tmp_path/'a/111/vrctst/molpro_scan/x/tasks/relaxed_scan'
    for path in (a,b,v): path.mkdir(parents=True)
    command=['sbatch','job.slurm']
    assert fair.take_turn(tmp_path,a,command,1)
    for path in (a,b,v): assert not fair.take_turn(tmp_path,path,command,0)
    assert not fair.take_turn(tmp_path,a,command,1)
    # An unserved group gets the next opportunity, irrespective of how
    # frequently the already served ANL group polls.
    assert fair.take_turn(tmp_path,v,command,1)
    assert not fair.take_turn(tmp_path,a,command,1)
    assert fair.take_turn(tmp_path,b,command,1)
    assert fair.take_turn(tmp_path,a,command,1)
    assert fair.take_turn(tmp_path,a,command,4)  # No artificial per-case cap.


def test_abandoned_ticket_expires_and_full_pool_still_registers_waiters(tmp_path, monkeypatch):
    clock=[0.]; monkeypatch.setattr(fair.time,'time',lambda:clock[0])
    a=tmp_path/'a'; b=tmp_path/'b'; a.mkdir(); b.mkdir()
    command=['sbatch','job.slurm']
    assert not fair.take_turn(tmp_path,a,command,0)
    clock[0]=1.
    assert not fair.take_turn(tmp_path,b,command,1)
    clock[0]=fair.LEASE_SECONDS+1
    assert fair.take_turn(tmp_path,b,command,1)
    state=json.loads((tmp_path/'.kinbot_policy_waiters.json').read_text())
    assert state['waiting']=={}


def test_multiple_available_slots_can_be_used_by_different_waiters(tmp_path):
    command=['sbatch','job.slurm']
    paths=[tmp_path/str(i) for i in range(3)]
    for path in paths:
        path.mkdir()
        assert not fair.take_turn(tmp_path,path,command,0)
    assert fair.take_turn(tmp_path,paths[1],command,2)
    assert fair.take_turn(tmp_path,paths[0],command,1)
    assert fair.take_turn(tmp_path,paths[2],command,1)
