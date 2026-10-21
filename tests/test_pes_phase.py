import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot import pes_phase as phase, live_control as control
from kinbot.execution import atomic_json


@pytest.fixture
def root(tmp_path, monkeypatch):
    monkeypatch.setattr(phase, '_ROOT', None)
    monkeypatch.setattr(phase, '_WORKER', None)
    monkeypatch.setattr(phase, '_EPOCH', None)
    return tmp_path


def state(root, **values):
    m = json.loads((root/phase.NAME).read_text())
    m.update(values)
    atomic_json(root/phase.NAME, m)
    return m


def worker_status(root, name, status):
    m = state(root)
    m['workers'][name]['status'] = status
    state(root, workers=m['workers'])


def test_serial_exploration_does_not_serialize_downstream(root):
    phase.initialize(root)
    phase.register(root, 'first', os.getpid())
    busy, status = phase.tick(root, ['first','second'], max_searches=1)
    assert busy == ['first'] and status == 'exploring'
    worker_status(root,'first','ready')
    busy, status = phase.tick(root, ['first','second'], max_searches=1)
    assert busy == [] and status == 'exploring'  # Free slot for second well.
    phase.register(root,'second', os.getpid())
    worker_status(root,'second','ready')
    assert phase.tick(root,['first','second'],max_searches=1) == ([], 'postprocessing')
    # Both may postprocess even though exploration concurrency is one.
    worker_status(root,'first','postprocessing')
    worker_status(root,'second','postprocessing')
    assert phase.tick(root,['first','second'],max_searches=1)[1] == 'postprocessing'
    worker_status(root,'first','assembled')
    assert phase.tick(root,['first','second'],max_searches=1)[1] != 'closing'
    worker_status(root,'second','assembled')
    assert phase.tick(root,['first','second'],max_searches=1)[1] == 'closing'


def test_discovered_well_extends_barrier(root):
    phase.initialize(root)
    phase.register(root,'a',os.getpid())
    worker_status(root,'a','ready')
    assert phase.tick(root,['a','new'],max_searches=1)[1] == 'exploring'


def test_live_reopen_quiesces_both_branches_before_retry(root, monkeypatch):
    control.configure(root)
    record = phase.initialize(root)
    for name in ('a','b'):
        (root/name).mkdir()
        phase.register(root,name,os.getpid())
        worker_status(root,name,'postprocessing')
    state(root,phase='postprocessing')
    target = root/'a'/'failed_reaction'
    control.search_scope(target,'reaction',-999,close=True)
    request = control.request(root,target,kind='reaction',reason='retry rejected channel')
    assert phase.tick(root,['a','b'],max_searches=1)[1] == 'exploring'
    assert state(root)['epoch'] == 1 and state(root).get('grants') == []
    monkeypatch.setattr(phase,'_ROOT',root)
    monkeypatch.setattr(phase,'_WORKER','a')
    monkeypatch.setattr(phase,'_EPOCH',0)
    monkeypatch.setenv('KINBOT_PES_PHASE_TOKEN',record['token'])
    with pytest.raises(phase.ExplorationReopened):
        phase.checkpoint()
    worker_status(root,'a','ready')
    phase.tick(root,['a','b'],max_searches=1)
    assert state(root).get('grants') == []
    worker_status(root,'b','ready')
    assert phase.tick(root,['a','b'],max_searches=1)[0] == ['a']
    control.acknowledge(root/'.kinbot_controls/requests'/(request['id']+'.json'),request,'applied')
    assert phase.tick(root,['a','b','new'],max_searches=1)[1] == 'exploring'


def test_closed_search_refuses_dead_owner_or_finalized_pes(root):
    control.configure(root)
    phase.initialize(root)
    (root/'a').mkdir()
    phase.register(root,'a',os.getpid())
    target = root/'a'/'reaction'
    control.search_scope(target,'reaction',-999,close=True)
    state(root,phase='closing')
    with pytest.raises(ValueError,match='handed off'):
        control.request(root,target,kind='reaction',reason='late')


def test_checkpoint_noop_outside_managed_pes(root):
    phase.checkpoint()
    phase.wait_for_exploration()
    phase.assembled()


def test_concurrent_discovery_publishes_one_input_without_overwriting(root,monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    import threading
    from kinbot import pes
    (root/'chemids').write_text('old\n')
    monkeypatch.setattr(pes,'routing_name',lambda s:'new')
    writes=[]
    monkeypatch.setattr(pes,'write_input',lambda *a:writes.append(a))
    barrier=threading.Barrier(2)
    def publish():
        barrier.wait()
        return pes.publish_well('input.json',object(),10,20,root,2)
    with ThreadPoolExecutor(2) as pool:
        futures=[pool.submit(publish) for _ in range(2)]
        assert sorted(f.result(timeout=5) for f in futures)==[False,True]
    assert len(writes)==1 and (root/'chemids').read_text()=='old\nnew\n'


def test_new_well_published_between_inventory_and_tick_keeps_barrier_closed(root):
    phase.initialize(root)
    phase.register(root,'a',os.getpid())
    worker_status(root,'a','ready')
    (root/'chemids').write_text('a\nnew\n')
    assert phase.tick(root,['a'],max_searches=1,known_jobs=['a'])[1]=='exploring'


def test_real_serial_workers_all_explore_before_any_postprocessing(root):
    import subprocess
    import sys
    import time
    record=phase.initialize(root)
    child='''from pathlib import Path
import time
from kinbot import pes_phase as p
assert p.attach()
Path('explored').write_text(str(time.time_ns()))
p.wait_for_exploration()
Path('postprocessing').write_text(str(time.time_ns()))
p.assembled()
'''
    processes={}
    env={**os.environ,'PYTHONPATH':str(Path(__file__).parents[1]),
         'KINBOT_PES_PHASE_ROOT':str(root),'KINBOT_PES_PHASE_TOKEN':record['token']}
    try:
        deadline=time.monotonic()+15
        while time.monotonic()<deadline:
            busy,status=phase.tick(root,['a','b'],max_searches=1)
            if not busy and len(processes)<2:
                name=('a','b')[len(processes)]
                (root/name).mkdir()
                proc=subprocess.Popen([sys.executable,'-c',child],cwd=root/name,env=env)
                processes[name]=proc
                phase.register(root,name,proc.pid)
            if len(processes)==2 and all(p.poll() is not None for p in processes.values()):
                break
            time.sleep(.03)
        assert len(processes)==2
        assert all(p.poll()==0 for p in processes.values())
        assert max(int((root/n/'explored').read_text()) for n in processes) < min(
            int((root/n/'postprocessing').read_text()) for n in processes)
    finally:
        for p in processes.values():
            if p.poll() is None:
                p.terminate()
            p.wait(timeout=5)
